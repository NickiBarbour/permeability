#' Format data for permeability package
#' 
#' @description Function can run on 1+ individuals and will output two 
#' dataframes, one for the barrier and one for the movement steps,
#' annotated with any desired covariates. 
#' Covariate(s) can be categorical (in character format) or numeric (in numeric
#' format).
#' Users can choose whether to annotate steps with covariate values using the: 
#' first location of the step, end location of the step, or if numeric,
#' the mean value across the step (mean of start/end location values).
#' If the covariate values are specific to the barrier, with a unique covariate 
#' value for each barrier segment, these can be provided as a vector equal in 
#' length to the number of barrier segments.
#' 
#' @param move.df a dataframe with columns for: `ID`, `Time` 
#' (POSIXct format, date and time), spatial coordinates (`X`/`Y` - units in 
#' meters), 
#' and any covariate(s)
#' @param barrier.covar.values a dataframe with the covariate values 
#' for each unique segment of the barrier - default is NULL. 
#' @param move.covar.names a vector with the name(s) of the column(s) with  
#' covariate values for each location of the movement data - default is NULL
#' @param move.covar.values if covariate values are provided, a vector of the 
#' calculations to perform to get covariate values for each step - options 
#' include the `start` , `end`, or if numeric, optionally the
#' `mean` or `difference` of the start/end location values for each step. 
#' Default is `start`.  To use different options for different covariates, 
#' please provide desired values in the same order as in `move.covar.names` 
#' (e.g., if you provide covars "A" and "B", and want the `start` for A and
#'  `mean` for B, `move.covar.values = c("start", "mean)`).
#' @param barrier either an 'sf' LINESTRING or a matrix containing the spatial 
#' coordinates (X, Y) and `line_id`, with a unique numeric value for each line 
#' segment. 
#' @param plot.track default TRUE - whether to plot each track with the given 
#' barrier or not
#' @param buffer.m If NULL (default) will be set to the max distance from the
#' barrier observed for steps that cross. Can optionally also be set to a
#' user-provided value - just ensure that this value accurately reflects 
#' the distance from the barrier within which steps have a viable possibility
#' of crossing.  
#' @param sample.rate regular desired sample rate (in hours) of input data. 
#' This is an important assumption of the package. Steps with an irregular 
#' sample rate will not be used to create the crossing table. It's the user's 
#' responsibility to ensure a (mostly) regular sample rate (some missing steps
#' are tolerable as they simply not be used but you need enough data to 
#' robustly estimate permeability) for their data prior 
#' to uploading it to the package functions. 
#' @param sample.rate.tolerance tolerance level for sample rate (in hours) - 
#'  defaults to 1 hour if `sample.rate.tolerance = NULL`.
#' @param parallel whether to run function in parallel or not (recommended for 
#' very large datasets)
#' @param n.cores numbers of cores to use if running in parallel.  Defaults to 
#' 80% of the number of available cores. 
#' @param verbose whether to print progress or not - default TRUE
#' 
#' @return `permdata` object, which is a list containing two "stepwise" data 
#' frames. The output `track` object has rows corresponding to each observed 
#' step in the data and columns for: `ID`, `Time`, start/end location of each 
#' step (complex, `Z.start`, `Z.end`), the movement step (Step, complex, start 
#' and locations of step), the time step (in hours, `D.time`), whether the step 
#' was in the final buffer or not (`In.buffer`, TRUE/FALSE - will be FALSE also 
#' if not at the desired sample rate), movement metrics 
#'  (step length, m - `L`, absolute angle, radians - `Phi`, relative turning
#'  angle, radians - `Theta`), whether the step crossed the barrier or not 
#'  (`Crossed` - TRUE/FALSE), the distance of the first location in each step
#'  from the barrier (`Dist_toBarrier`, in meters), and any covariate(s) 
#'  (same name as provided in original data). The output `barrier` has rows 
#'  corresponding to each barrier segment and columns for the complex locations
#'  of the start and end of the barrier segments (`Z1` and `Z2`), the
#'  id of each barrier segment (`barrier.id`), the line id
#'  (`line_id`), and any covariate values (same names as input).
#' 
#' @example examples/examples_prepPermeability.R
#' @export


prepPermeability <- function(move.df, barrier, 
                             sample.rate, 
                             sample.rate.tolerance = NULL,
                             buffer.m = NULL, 
                             barrier.covar.values = NULL,
                             move.covar.names = NULL, 
                             move.covar.values = "start", 
                             plot.track = FALSE, verbose = FALSE,
                             parallel = FALSE, n.cores = NULL){
  # check if barrier is in sf format
  if("sf" %in% class(barrier)){
    barrier_xy <- st_coordinates(barrier)
    if(ncol(barrier_xy)==3){
      colnames(barrier_xy) <- c("X","Y","line_id")
      
      barrier_df <- data.frame(barrier_xy)
    }else{
      stop("Barrier is an sf object. It should be a LINESTRING.")
    }
  }else{
    barrier_df <- barrier |> data.frame()
    
    if(ncol(barrier_df)!=3){
      stop("Barrier should have 3 columns: X, Y, and line_id.")
    }
  }
  # convert barrier to "stepwise" matrix
  convert_xy_2_step <- function(x){
    z <- x$X + 1i*x$Y
    step <- cbind(z[-length(z)],z[-1])
    step <- data.frame(step)
    step$line_id <- unique(x$line_id)
    colnames(step) <- c("Z1","Z2","line_id")
    return(step)
  }
  barrier_step_df <- plyr::ddply(barrier_df, "line_id", convert_xy_2_step)
  
  # add covariate values
  if(!is.null(barrier.covar.values)){
    if(nrow(barrier.covar.values)!=nrow(barrier_step_df)){
      stop(paste("Covariate values for barrier do not match dimensions for the number of barrier segments. Barrier has", nrow(barrier_step_df), "unique segments and provided covariate data has", nrow(barrier.covar.values), "rows."))
    }
    barrier_step_df <- cbind(barrier_step_df, barrier.covar.values)
  }
  
  # Check if max step length < barrier dist
  ## add columns for step length and time step
  df <- move.df
  df.step <- df |> dplyr::mutate(z.move = X + 1i*Y) |> 
    dplyr::group_by(ID) |>
    mutate(L = c(NA,Mod(diff(X + 1i*Y))),
           dT = c(NA, round(difftime(Time[-1], 
                                   Time[-length(Time)], 
                                   units = "hours"))))
  
  if(is.null(sample.rate.tolerance)){
    sample.rate <- c(sample.rate-1, sample.rate, sample.rate+1)
  }else{
    sample.rate <- unique(c(sample.rate-sample.rate.tolerance, 
                     sample.rate, 
                     sample.rate+sample.rate.tolerance))
  }
  
  # retain og buffer value
  buffer_og <- buffer.m
  
  if(is.null(buffer.m)){
    # set temporarily to max steplength
    longstep <- max(df.step$L[df.step$dT %in% sample.rate], na.rm = TRUE)
    
    buffer.m <- longstep
  }
  
  # change ID to character
  df$ID <- as.character(df$ID)
  
  if(length(unique(df$ID)) > 1){
    id_list <- split(df, df$ID)
    
    if(parallel){
      if(is.null(n.cores)){
        n.cores <- round(parallel::detectCores()*0.80) # use 80% of cores
      } 
      if(verbose) print(paste(n.cores, "cores are being used: running in parallel."))
      
      plan(multisession, workers = n.cores)
      
      data_list <- future_lapply(id_list, 
                                 prepPermeability_individual,
                                 future.seed = TRUE,
                                 sample.rate = sample.rate,
                                 buffer.m = buffer.m,
                                 covar.names = move.covar.names, 
                                 covar.values = move.covar.values, 
                                 barrier = barrier_step_df, 
                                 plot.track = plot.track, 
                                 verbose = verbose)
      plan(sequential)
      
    }else{
      
      data_list <- lapply(id_list, 
                          prepPermeability_individual,
                          sample.rate = sample.rate,
                          buffer.m = buffer.m,
                          covar.names = move.covar.names, 
                          covar.values = move.covar.values, 
                          barrier = barrier_step_df, 
                          plot.track = plot.track, 
                          verbose = verbose)
    }
    # drop any NULL elements
    data_list2 <- data_list[!sapply(data_list,is.null)]
    data <- do.call("rbind", data_list2)
    
    if(verbose) print(paste(paste("Done.", sum(data$Crossed, na.rm = TRUE)), "steps crossed"))
  }else{
    data <- prepPermeability_individual(df = df,
                                        sample.rate = sample.rate,
                                        buffer.m = buffer.m,
                                        covar.names = move.covar.names, 
                                        covar.values = move.covar.values, 
                                        barrier = barrier_step_df, 
                                        plot.track = plot.track, 
                                        verbose = verbose)
    if(verbose) print(paste(paste("Done.", sum(data$Crossed, na.rm = TRUE)), "steps crossed"))
  }
  
  # change ID to factor
  data$ID <- as.factor(data$ID)
  
  if(is.null(buffer_og) & sum(data$Crossed, na.rm = TRUE) > 0){
    # set buffer to max crossing distance from barrier
    buffer.m <- max(subset(data, In.buffer & Crossed & D.time %in% sample.rate)$Dist_toBarrier)
    
    data$In.buffer <- FALSE
    
    data[which(data$Dist_toBarrier <= buffer.m),]$In.buffer <- TRUE
  }else if(sum(data$Crossed, na.rm = TRUE) == 0){
    warning("No steps crossed the barrier. The buffer distance is set to the maximum observed steplength in the input data.")
  }
  
  # add barrier id to barrier
  barrier_step_df$barrier.id <- 1:nrow(barrier_step_df)
  
  permdata <- list(track = data, 
                   barrier = barrier_step_df)
  
  class(permdata) <- "permdata"  
  
  attr(permdata, "buffer.size") <- buffer.m
  
  attr(permdata, "sample.rate") <- sample.rate
  
  return(permdata)
}

prepPermeability_individual <- function(df, covar.names = NULL, 
                                        covar.values = "start", 
                                        barrier, 
                                        buffer.m, sample.rate,
                                        plot.track = TRUE, verbose = TRUE){
  if(verbose) print(paste("Formatting Data for ID:", df$ID[1]))
  
  # make sure it's a data frame
  df <- df |>
    data.frame()
  
  # check that there are > 4 observations
  if(nrow(df) < 4){
    if(verbose) warning(paste(df$ID[1], "has less than 4 observations. It will be assigned NA's."))
    
    if(!is.null(covar.names)){
      
      covar_cols <- df[NA, covar.names]
      
      if(!is.data.frame(covar_cols)){
        covar_cols <- data.frame(covar = covar_cols)
        colnames(covar_cols) <- covar.names
      }
      
      df.step <- df |> 
        reframe(Z.start = NA,
                Z.end = NA,
                D.time = NA,
                ID = unique(ID),
                Time = unique(Time),
                Step = NA,
                L = NA,
                Phi = NA,
                Theta = NA,
                Step.ID = NA)
      
      # add covariate columns
      df.step <- cbind(df.step, covar_cols)
      
      # add additional columns
      df.step$In.buffer <- NA
      df.step$Dist_toBarrier <- NA
      df.step$Crossed <- NA
      
      return(df.step)
      
    }else{
      
      df.step <- df |> 
        reframe(Z.start = NA,
                Z.end = NA,
                D.time = NA,
                ID = unique(ID),
                Time = unique(Time),
                Step = NA,
                L = NA,
                Phi = NA,
                Theta = NA,
                Step.ID = NA,
                In.buffer = NA,
                Dist_toBarrier = NA,
                Crossed = NA)
      
      return(df.step)
      
    }
  }else{
    # extract Z coords of id and barrier to plot
    xy.move <- cbind(df$X, df$Y)
    z.move <- xy.move[,1] + 1i*xy.move[,2]
    
    if (plot.track){
      barrier_step <- cbind(barrier[,1], barrier[,2])
      
      plot(z.move, cex=0.7, asp=1, type="l", col="blue")
      points(z.move, cex=0.7, asp=1)
      apply(barrier_step, 1, function(x) lines(x,col = "red"))
      title(paste("ID:",unique(df$ID)))
    }
    
    if(!is.null(covar.names)){
      # create step-wise df
      df.step <- df |> mutate(Z = X + 1i*Y) |> 
        plyr::summarize(Z.start = Z[-length(Z)],
                        Z.end = Z[-1],
                        D.time = round(difftime(Time[-1], 
                                                Time[-length(Time)], 
                                                units = "hours")) |> 
                          as.numeric(),
                        ID = ID[-length(ID)],
                        Time = Time[-length(Time)]) |> 
        plyr::mutate(Step = Z.end - Z.start,
                     L = Mod(Step), 
                     Phi = Arg(Step), 
                     Theta = c(NA,diff(Phi))) |>
        plyr::mutate(Step.ID = cumsum(rep(1, length(Step))))
      
      # create dataframe of covar names and values
      covar.df <- data.frame(names = covar.names,
                             values = covar.values)
      df_covarcols <- df[, covar.names]
      # if only one column, will be vector
      if(is.data.frame(df_covarcols)){
        col_names <- colnames(df_covarcols)
        df.step.covars <- df_covarcols[-1,]
      }else{
        df.step.covars <- data.frame(covar = df_covarcols[-1])
        colnames(df.step.covars) <- covar.names
      }
      
      # use for-loop to determine values for each covar column
      for( i in c(1:ncol(df.step.covars))){
        
        if(is.data.frame(df_covarcols)){
          column <- df_covarcols[,i]
          
          colname <- col_names[i]
          
          colvalue <- subset(covar.df, names == colname)$values
        }else{
          column <- df_covarcols
          
          colname <- covar.names
          
          colvalue <- subset(covar.df, names == colname)$values
        }
        
        if(colvalue=="start"){
          
          new_col <- column[-length(column)]
        }else if(colvalue=="end"){
          
          new_col <- column[-1]
        }else if(colvalue=="mean"){
          
          if(!is.numeric(column)) stop("Can't take mean of non-numeric covariate column. Please specify a different value.")
          
          column_step <- cbind(column[-1], 
                               column[-length(column)])
          new_col <- apply(column_step, 1, function(x) mean(x))
        }else if(colvalue=="difference"){
          
          if(!is.numeric(column)) stop("Can't take difference of non-numeric covariate column. Please specify a different value.")
          
          column_step <- cbind(column[-1], 
                               column[-length(column)])
          new_col <- apply(column_step, 1, function(x) abs(x[1]-x[2]))
        }
        df.step.covars[,i] <- new_col
      }
      # add covariate columns to step-wise df
      df.step <- cbind(df.step, df.step.covars)
      
    }else{
      df.step <- df |> mutate(Z = X + 1i*Y) |> 
        plyr::summarize(Z.start = Z[-length(Z)],
                        Z.end = Z[-1],
                        D.time = round(difftime(Time[-1], 
                                                Time[-length(Time)], 
                                                units = "hours")) |> 
                          as.numeric(),
                        ID = ID[-length(ID)],
                        Time = Time[-length(Time)]) |> 
        plyr::mutate(Step = Z.end - Z.start,
                     L = Mod(Step), 
                     Phi = Arg(Step), 
                     Theta = c(NA,diff(Phi))) |>
        plyr::mutate(Step.ID = cumsum(rep(1, length(Step))))
    }
    # determine if steps are within buffer radius of barrier
    z.barrier <- barrier[,"Z1"] # use first loc of each barrier seg
    dist_starttobarrier <- getClosest(df.step$Z.start, z.barrier)
    dist_endtobarrier <- getClosest(df.step$Z.end, z.barrier)
    df.step$In.buffer <- dist_starttobarrier < buffer.m & 
      dist_endtobarrier < buffer.m
    
    # new column in step df: distance of first loc in step to barrier
    df.step$Dist_toBarrier <- dist_starttobarrier
    
    # check if any steps in buffer and at desired sample rate
    if(nrow(subset(df.step, In.buffer & 
                   D.time %in% c(sample.rate)))>0){
      max.step <- max(subset(df.step, In.buffer & 
                               D.time %in% c(sample.rate))$L)
      
      step.barrier <- cbind(barrier[,1], barrier[,2])
      step.move <- cbind(start = df.step$Z.start, end = df.step$Z.end)
      
      # set empty crossing column
      df.step$Crossed <- NA
      
      # loop through steps/rows and determine for each if they cross the barrier
      for(i in c(1:nrow(df.step))){
        my.step <- step.move[i,]
        
        # check if step is in buffer and has regular sample rate
        if(df.step[i,]$In.buffer & 
           df.step[i,]$D.time %in% c(sample.rate)){
          # crop barrier to only the segments near the observed movement step
          ## speeds up processing time
          my.barrier <- step.barrier[getClosest(step.barrier[,1], 
                                                my.step[2]) < max.step & 
                                       getClosest(step.barrier[,2], 
                                                  my.step[2]) < max.step,] |>
            matrix(ncol = 2)
          
          if(nrow(my.barrier)==0){
            df.step$Crossed[i] <- FALSE
          } else {
            crossing <- any(apply(my.barrier, 1, isIntersect, 
                                  z2 = c(my.step[1], my.step[2])))
            if(crossing){
              df.step$Crossed[i] <- TRUE
            }else{
              df.step$Crossed[i] <- FALSE
            }
          }
        }else{
          df.step$Crossed[i] <- NA
        }
      }
      
    }else{
      df.step$Crossed <- NA
    }
  }
  return(df.step)
}



#' @export
getClosest <- function(z1, z2){
  D <- outer(z1, z2, function(z1,z2) Mod(z1-z2))
  apply(D, 1, min)
}

#' @export
summary.permdata <- function(x){
  
  data <- subset(x$track, In.buffer)
  
  N_ID <- length(unique(data$ID))
  
  Sample_rate <- paste(attr(x, "sample.rate"), collapse = " ")
  
  data$Date <- as.Date(data$Time)
  
  Date_range <- paste(min(data$Date), max(data$Date), sep = "-")
  
  Steps_inbuffer <- nrow(subset(data, In.buffer))
  
  N_crossed <- sum(data$Crossed, na.rm = TRUE)
  
  if(N_crossed > 0){
    Mean_SL <- round(mean(subset(data, Crossed & In.buffer)$L, na.rm = TRUE), 
                     digits = 2)
  }else{
    Mean_SL <- NA
  }
  
  Buffer_d <- round(attr(x, "buffer.size"), digits = 2)
  
  data_summary <- data.frame("N.IDs in buffer" = N_ID,
                             "Sample rate (hrs)" = Sample_rate,
                             "Date range" = Date_range,
                             "N.steps in buffer" = Steps_inbuffer,
                             "N.steps to cross" = N_crossed,
                             "Mean crossing steplength (m)" = Mean_SL,
                             "Buffer dist (m)" = Buffer_d)
  
  data_summary
}

#'@export
print.permdata <- function(x){
  print(summary(x))
}
