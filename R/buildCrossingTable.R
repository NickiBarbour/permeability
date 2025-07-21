#' Create Crossing Table of Animal Crossings with Respect to Barrier
#' @name buildCrossingTable
#' @description Function runs on one individual at a time and outputs one large 
#' crossing table for all individuals. Simulates user-defined number of null 
#' steps around each observed, actual step and determines for 
#' all (actual and null) steps, whether each crossed the barrier or not.
#' Any individuals with less than 2 observations or without tracks within the
#' defined buffer radius of the barrier are dropped. Note that the user will 
#' need to annotate the output crossing table from this function with their 
#' desired covariates.
#' @param permdata object outputted from either \link{prepPermeability} or in 
#' case of simulated data, \link{simulatePermeability}, consisting of 
#' dataframes for the movement track(s) and barrier (`track` and `barrier`).
#' @param n.null number of null steps to simulate around each observed 
#' step - default 60
#' @param step.lengths vector of user-supplied step lengths (in meters) to 
#' sample from to create null steps. Step lengths can be obtained from the 
#' `permdata` and `track` object from \link{prepPermeability}, under the `L`
#' column. Ensure step lengths accurately represent steps that could cross the 
#' barrier and have a large enough sample size ( >= 30). If NULL 
#' (default), will use step lengths that cross the barrier (`Crossed = TRUE` 
#' from the `permdata` `track` object). If there are < 30 steplengths that cross 
#' the barrier or in the provided `step.lengths` vector, will use all
#' steplengths that cross the barrier (`step.lengths = NULL`) or the provided
#' vector of steplengths to fit and sample from a Weibull distribution, using 
#' the `fitdistr` function from the `MASS` package
#' with a starting scale value equal to the mean steplength for steps that 
#' cross the barrier.
#' @param turning.angles vector of user-supplied turning angles (in radians).
#'  Turning angles can be obtained from the `permdata` and `track` object from 
#' \link{prepPermeability} under the `Theta` column. If NULL
#' (default), will sample from turning angles within the buffer 
#' (`In.buffer = TRUE` from the `permdata` `track` object, a.k.a the max 
#' crossing steplength distance from the barrier.
#' @param parallel whether to run function in parallel or not (recommended for 
#' very large datasets)
#' @param n.cores numbers of cores to use if running in parallel.  Defaults to 
#' 80% of the number of available cores. 
#' @param verbose whether to print progress/messages or not - default TRUE
#' @return Data frame with rows corresponding to each observed step in the data 
#' and columns for: `ID`, `Step.ID` (id of the step, as output from
#'  \link{prepPermeability}), `barrier.id` (id of the barrier segment that the 
#'  majority of null steps crossed), start/end location of each step 
#'  (complex, `Z.start`, `Z.end`), `Time` (character format), number of 
#' null steps to cross the barrier (`null.crossed`), number of null 
#' steps to NOT cross the barrier (`null.stayed`), whether the observed step 
#' crossed the barrier or not (`crossed`), and if so, the barrier segment it 
#' crossed at (`barrier.crossed`).
#' @example examples/examples_buildCrossingTable.R
#' @export
buildCrossingTable <- function(permdata, parallel = FALSE, n.cores = NULL,
                               n.null = 60, 
                               step.lengths = NULL,
                               turning.angles = NULL, 
                               verbose = FALSE){
  if(!inherits(permdata, "permdata")){
    stop("permdata should be of class 'permdata`. Did you forget to run your data through the `prepPermeability` function?")
  }
  df <- permdata$track
  barrier <- permdata$barrier
  # drop ids not in buffer
  # and keep only those steps where:
  ## dtime == sample.rate ( default +/- 1 hr)
  sample.rate <- attr(permdata, "sample.rate")
  ## subset for steps within buffer
  ## and ensure there are no duplicate locations (step lengths > 0)
  df <- df |> 
    subset(D.time %in% c(sample.rate)) |> 
    subset(In.buffer==TRUE) |>
    subset(L>0) |>
    dplyr::mutate(Time = as.character(Time))
  if(nrow(df)==0){
    warning("No steps within buffer at desired sample rate (returning NULL).")
    cT <- NULL
    return(cT)
  }
  # if no steplengths provided:
  ## use steplengths that cross
  if(is.null(step.lengths)){
    if(nrow(subset(df, Crossed))==0) stop("No steps crossed the barrier. Please provide a vector of steplengths with magnitudes that would reasonably represent crossing steps to create null steps from.")
    if(length(subset(df, Crossed)$L) < 30){
      if(verbose) message("Less than 30 steps crossed. A Weibull distribution will be estimated for steplengths, using start parameters of shape 1 and scale equal to the mean crossing steplength, and fit to crossing steps to create the steplength null distribution.")
      step.lengths <- na.omit(subset(df, Crossed)$L)
      # find the mean crossing steplength
      mean_cross_SL <- mean(subset(df, Crossed)$L, na.rm = TRUE)
    }else{
      step.lengths <- na.omit(subset(df, Crossed)$L)
      # find the mean crossing steplength
      mean_cross_SL <- NA
    }
  }else{
    mean_cross_SL <- NA
  }
  # if no turning angles provided, use turning angles within buffer
  if(is.null(turning.angles)){
    turning.angles <- na.omit(df$Theta)
    
    if(length(turning.angles) < 30){
      warning("Less than 30 steps within the buffer, resulting in < 30 turning angles to sample from to create null steps. Consider simulating turning angles and providing them to the `turning.angles` argument.")
    }
  }
  
  if(length(unique(df$ID)) > 1){
    df$ID <- as.character(df$ID)
    id_list <- split(df, df$ID)
    
    if(parallel){
      if(is.null(n.cores)){
        n.cores <- round(parallel::detectCores()*0.80) # use 80% of cores
      } 
      plan(multisession, workers = n.cores)
      if(verbose) message(paste(n.cores, "cores are being used: running in parallel."))
      cT_list <- future_lapply(id_list, 
                               FUN = buildCrossingTable_individual,
                               future.seed = TRUE,
                               barrier = barrier,
                               n.null = n.null,  
                               step.lengths = step.lengths,
                               turning.angles = turning.angles,
                               verbose = verbose,
                               mean_cross_SL = mean_cross_SL)
      plan(sequential)
    }else{
      cT_list <- lapply(id_list, 
                        buildCrossingTable_individual,
                        barrier = barrier,
                        n.null = n.null,  
                        step.lengths = step.lengths,
                        turning.angles = turning.angles,
                        verbose = verbose,
                        mean_cross_SL = mean_cross_SL)
    }
    # drop any NULL elements
    cT_list2 <- cT_list[!sapply(cT_list,is.null)]
    cT <- do.call("rbind", cT_list)
  }else{
    df$ID <- as.character(df$ID)
    
    cT <- buildCrossingTable_individual(df = df,
                                        barrier = barrier,
                                        n.null = n.null,  
                                        step.lengths = step.lengths,
                                        turning.angles = turning.angles,
                                        verbose = verbose,
                                        mean_cross_SL = mean_cross_SL)
  }
  if(is.null(cT)) stop("No crossing table to return. Check your data.")
  # check for rows where actual crossed and possible did not:
  which.bad <- subset(cT, null.crossed == 0 & crossed == 1)
  
  if(verbose & nrow(which.bad) > 0) {
    # remove rows where crossed but no possible crossed
    message("You may need to either a) increase the null step sample size (n.null) or b) check that the distribution of step lengths sampled from accurately reflects the distribution of actual steps that crossed the barrier; ", 
            "there were ", nrow(which.bad), " crossings where no null crossings occurred.")
  }
  
  # drop any rows where there are no null crossings
  ## don't want any rows where crossed = 0 & null crossed = 0
  ## or rows where crossed = 1 and null crossed = 0
  cT <- subset(cT, null.crossed > 0)
  
  if(nrow(cT)==0){
    stop("No null steps crossed barrier, resulting in an empty crossing table. Check that step length distributions are correct and ensure that you have data near the barrier.")
  }
  
  cT$ID <- as.factor(cT$ID)
  # convert time to POSIXct, accounting for times at 00:00:00
  Time <- as.POSIXct(cT$Time, format = "%Y-%m-%d %H:%M:%S")
  Time2 <- as.POSIXct(cT$Time)
  cT$Time <- as.POSIXct(ifelse(is.na(Time), Time2, Time))
  
  if(length(unique(cT$ID)) > 1 & verbose) message("Done. ", length(unique(cT$ID)), " track(s) had steps with null crossings of the barrier and are retained in the crossing table.")
  
  return(cT)
}


buildCrossingTable_individual <- function(df, barrier, n.null, 
                                          step.lengths, turning.angles,
                                          verbose = TRUE, 
                                          mean_cross_SL){
  
  if(verbose) print(paste("Building crossing table for ID:", df$ID[1]))
    # extract max step and use 2x this to only retain barrier segments
    ## near each movement step
    max_step <- max(df$L, na.rm=TRUE)*2
    
    step.move <- cbind(df$Z.start, df$Z.end)
    # create blank crossing table
    crossingTable_df <- data.frame(matrix(NA, nrow = nrow(df), 
                                          ncol = 10))
    colnames(crossingTable_df) <- c("ID","Step.ID", "barrier.id","Z.start",
                                    "Z.end","Time",
                                    "null.crossed", "null.stayed", "crossed", 
                                    "barrier.crossed")

    # loop through each step in the dataframe to populate crossing table
    for(i in 1:nrow(df)){
      my.step <- df[i,]
      # extract just complex locs
      step.barrier <- cbind(barrier$Z1, barrier$Z2)
      # assign barrier id name using row id
      row.names(step.barrier) <- 1:nrow(step.barrier)
      # find barrier rows that are close to current
      close_rows <- unique(which(getClosest(step.barrier[,1], 
                                            my.step$Z.end) < max_step & 
                                   getClosest(step.barrier[,2], 
                                              my.step$Z.end) < max_step))
      # if no segments are within max step distance, set Nc to 0
      if(length(close_rows)==0){
        Nc <- 0
        Ns <- n.null
        switch <- FALSE
        barrier.id <- NA
        barrier.crossed <- NA
      }else{
        # check for step.barrier = 1 segment
        if(length(close_rows)==1){
          crossings <- isIntersect(step.barrier,
                                   z2 = c(my.step$Z.start, my.step$Z.end))
          switch_value <- any(crossings)
          # check if there crossings, if so find barrier segment where it
          ## crossed - if not, set to NA
          barrier.crossed <- ifelse(switch_value, 
                                    close_rows,
                                    NA)
          
          crossing.list <- getCandidateCrosses(my.step, df, step.barrier, 
                                               n.null, 
                                               step.lengths, turning.angles,
                                               close_rows, 
                                               mean_cross_SL = mean_cross_SL)
        }else{
          step.barrier <- step.barrier[getClosest(step.barrier[,1], 
                                                  my.step$Z.end) < max_step & 
                                         getClosest(step.barrier[,2], 
                                                    my.step$Z.end) < max_step,]
          
          crossings <- apply(step.barrier, 1, isIntersect, 
                             z2 = c(my.step$Z.start, my.step$Z.end))
          switch_value <- any(crossings)
          # check if there crossings, if so find barrier segment where it
          ## crossed - if not, set to NA
          barrier.crossed <- 
            ifelse(switch_value,
                   as.numeric(names(crossings[crossings==TRUE])),NA)
          
          crossing.list <- getCandidateCrosses(my.step, df, step.barrier, 
                                               n.null, step.lengths,
                                               turning.angles,
                                               close_rows = NULL,
                                               mean_cross_SL = mean_cross_SL)
        }
        Nc <- crossing.list$N.crossed
        Ns <- crossing.list$N.stayed
        switch <- as.numeric(switch_value)
        barrier.id <- crossing.list$barrier.id
        
        }
      crossingTable_df[i,1] <- my.step$ID
      crossingTable_df[i,2] <- my.step$Step.ID
      crossingTable_df[i,3] <- barrier.id
      crossingTable_df[i,4] <- my.step$Z.start
      crossingTable_df[i,5] <- my.step$Z.end
      crossingTable_df[i,6] <- my.step$Time
      crossingTable_df[i,7] <- Nc
      crossingTable_df[i,8] <- Ns
      crossingTable_df[i,9] <- switch
      crossingTable_df[i,10] <- barrier.crossed
      }
      
      crossingTable <- subset(crossingTable_df, is.na(crossed)==FALSE)
      
      return(crossingTable)
}

# define function to simulate null steps
## determine if each null step crosses barrier or not
# null steps are simulated using provided step length distributions
## if null steps cross barrier, barrier segment that crossed the most
## is also extracted as "barrier.id"
#' @export
getCandidateCrosses <- function(my.step, df, step.barrier, n.null, 
                                step.lengths, turning.angles,
                                close_rows = NULL,
                                plot_candidates = FALSE,
                                mean_cross_SL){
  
  if(is.na(mean_cross_SL)==FALSE){
    weib_try <- try(fitdistr(step.lengths, "weibull", 
                             start = list(shape = 1, scale = mean_cross_SL)), 
                    silent = TRUE)
    
    if(!inherits(weib_try, "try-error")){
      weib.fit <- fitdistr(step.lengths, "weibull", 
                           start = list(shape = 1, scale = mean_cross_SL))
      scale_fit <- weib.fit$estimate["scale"]
      shape_fit <- weib.fit$estimate["shape"]
      
      Ls <- rweibull(n = n.null, shape = shape_fit, scale = scale_fit)
    } else{
      message("Not enough crossing steps to estimate Weibull distribution parameters. Weibull distribution will be fit to crossing steps using shape 2 and scale equal to the mean crossing steplength.")
      Ls <- rweibull(n = n.null, shape = 2, 
                     scale = mean_cross_SL/(sqrt(pi)/2))
    }
    
    thetas <- sample(turning.angles, n.null, replace = TRUE)
    
    # generate null steps
    step.candidate <- cbind(start = my.step$Z.start, 
                            end = my.step$Z.start + 
                              complex(mod = Ls, arg = my.step$Phi + thetas))
  }else{
    Ls <- sample(step.lengths, n.null, replace = TRUE)
    
    thetas <- sample(turning.angles, n.null, replace = TRUE)
    
    # generate null steps
    step.candidate <- cbind(start = my.step$Z.start, 
                            end = my.step$Z.start + 
                              complex(mod = Ls, arg = my.step$Phi + thetas))
  }
  if(plot_candidates){
    plot(step.candidate, col="white", type = "l")
    apply(step.barrier, 1, function(x) lines(x, lwd = 2))
    segments(Re(step.candidate[,1]), Im(step.candidate[,1]), 
             Re(step.candidate[,2]), Im(step.candidate[,2]), lwd = 2, 
             col = "orange")
    segments(Re(my.step$Z.start), Im(my.step$Z.start), 
             Re(my.step$Z.end), Im(my.step$Z.end), lwd = 2, 
             col = "green")
  }
  # check for step.barrier = 1 segment only
  if(!is.null(close_rows)){
    switch.candidate <- apply(step.candidate, 1, isIntersect, z2 = step.barrier)
    
    N.crossed <- sum(switch.candidate)
    N.stayed <- sum(switch.candidate==0)
    
    if(N.crossed==0 | is.na(N.crossed)){
      barrier.id <- NA
    }else{
      barrier.id <- close_rows
    }
    
  }else{
    switch.candidate <- apply(step.barrier, 1, function(z2){
      apply(step.candidate, 1, isIntersect, z2 = z2)
    })
    
    N.crossed <- sum(rowSums(switch.candidate)!=0)
    N.stayed <- sum(rowSums(switch.candidate)==0)
    
    if(N.crossed==0 | is.na(N.crossed)){
      barrier.id <- NA
    }else{
      barrier.id <-
        rownames(step.barrier)[which.max(colSums(switch.candidate))]|> 
        as.numeric()
    }
  }
  
  return(candidate_crosses = list(N.crossed = N.crossed, 
                                  N.stayed = N.stayed, 
                                  barrier.id = barrier.id))
}

#' @export
getClosest <- function(z1, z2){
  D <- outer(z1, z2, function(z1,z2) Mod(z1-z2))
  apply(D, 1, min)
}

#' @export
isIntersect <- function(z1, z2){
  theta <- Arg(z1[2]-z1[1])
  z1.t <- (z1-z1[1]) * exp(-1i*theta)
  z2.t <- (z2-z1[1]) * exp(-1i*theta)
  slope.r <- (Re(z2.t[2]) - Re(z2.t[1]))/(Im(z2.t[2]) - Im(z2.t[1]))
  x.intercept <-  Re(z2.t[1]) - slope.r * Im(z2.t[1])
  intersect.x <- x.intercept > 0 & x.intercept < Re(z1.t[2])
  intersect.x & (Im(z2.t[2]) * Im(z2.t[1])) < 0
}
