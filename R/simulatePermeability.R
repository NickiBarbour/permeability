#' Simulate track(s) with a barrier
#' 
#' @description Simulate track(s), where each track is either a correlated 
#' random walk (default) or biased random walk, with a set permeability with
#' respect to a linear barrier and the option to have permeability vary
#' as a function of covariate values associated with barrier segments.
#' @param xmax maximum x-value of bounding box
#' @param ymax maximum y-value of bounding box
#' @param N.steps number of steps to simulate (length of simulated track)
#' @param theta.rho turning angle concentration parameter, to simulate track
#' turning angles using a wrapped Cauchy distribution
#' @param step.shape shape parameter to simulate track step lengths, using a
#' Weibull distribution
#' @param step.scale scale parameter to simulate track step lengths, using a
#' Weibull distribution
#' @param barrier three column matrix or data frame of barrier vertices 
#' (named `X` and `Y`) and line ID (named `line_id` - should be equal to 1 if 
#' there is only one line segment in the barrier)
#' @param n.segments if `barrier = NULL`, number of segments for simulated 
#' barriers  - for a barrier with multiple lines, assumes an equal 
#' number of segments per line.
#' @param n.null number of null steps to simulate around each actual step
#' to estimate the null probability of crossing. Default 30.
#' @param kappa desired permeability value of the barrier -
#' values closer to 0 represent a less permeable barrier and values
#' closer to 1 represent a more permeable barrier.  Should be `NULL` if using
#' a covariate and `beta`.. 
#' @param beta a vector of two coefficients (intercept and main effect) for the 
#' impact of the covariate on permeability - default `NULL` (no covariate). 
#' @param covar a numeric vector of covariate values for each barrier segment, 
#' equal in length to `n.segments` - default `NULL` (if `NULL`, function will 
#' simulate covar values on a small-to-big numeric scale)
#' @param n.tracks number of tracks to simulate - default 1
#' @param type whether to simulate a correlated random walk (`crw`) or biased 
#' correlated random walk (`bcrw`) - if the latter, will generate a random
#' target location, unless `target.location` is provided. Default type is `crw`.
#' @param target.location target location (two column matrix of `X`, `Y` 
#' locations) for biased random walk, if `type = "bcrw"`. 
#' @param plot.track whether to plot resulting track - default `FALSE`.
#' @param parallel if generating multiple tracks, whether to run in parallel or
#' not (recommended). Default `FALSE`.
#' 
#' @return `permdata` object, which is a list containing two "stepwise" data 
#' frames object. The `track` object has rows corresponding to each observed 
#' step in the data and columns for: `ID`, `Time`, start/end location of each 
#' step (complex, `Z.start`, `Z.end`), the movement step (`Step`, complex, start 
#' and locations of step), the time step (in hours, `D.time`), whether the step 
#' was in the buffer (`In.buffer`, TRUE/FALSE), movement metrics 
#'  (step length, m - `L`, absolute angle, radians - `Phi`, relative turning
#'  angle, radians - `Theta`), and whether the step crossed the barrier or not 
#'  (`Crossed` - TRUE/FALSE). The `barrier` object has rows 
#'  corresponding to each barrier segment and columns for the complex locations
#'  of the start and end of the barrier segments (`Z1` and `Z2`), the line id
#'  (`line_id`), covariate values (`covar`), the
#'  id of each barrier segment (`barrier.id`), and "true" kappa values 
#'  (`kappa`).
#' @example examples/examples_simulatePermeability.R
#' @export

simulatePermeability <- function(xmax = 50, ymax = 50, 
                              N.steps = 500, theta.rho = 0.5, 
                              step.shape = 5, step.scale = 10, 
                              barrier = NULL,
                              n.segments = NULL,
                              n.null = 30,
                              kappa = NULL, beta = NULL,
                              covar = NULL,
                              n.tracks = 1,
                              type = "crw",
                              target.location = NULL,
                              plot.track = FALSE,
                              parallel = FALSE){
  
  # extract (or simulate) barrier coordinates
  if(is.null(barrier)){
    if(is.null(n.segments)){
      stop("Barrier needs to be simulated but n.segments has not been defined.")
    }
    barrier <- cbind(X = xmax * rnorm(n.segments+1, sd = xmax / 5e3),
                     Y = ymax * seq(-1,1,length = n.segments+1),
                     line_id = 1)
    barrier <- data.frame(barrier)
  } else{ 
    
    barrier <- data.frame(barrier)
    if("X" %in% colnames(barrier) & "Y" %in% colnames(barrier) &
       "line_id" %in% colnames(barrier)) {
      
      } else stop("Barrier must have at least three columns (X, Y, and line_id).") 
    segments_df <- barrier |> dplyr::group_by(line_id) |> 
      dplyr::reframe(n = length(X)-1)
    n.segments <- unique(segments_df$n)
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
  barrier_step_df <- plyr::ddply(barrier, "line_id", convert_xy_2_step)
  
  # simulate covariate if no "covar" supplied
  if(!is.null(beta)){
    if(is.null(covar)){
      covar <- sort(runif(n.segments, 1, 20)) # small to big (south to north)
      # annotate each segment for each barrier line with covar values
      barrier_step_df <- barrier_step_df |> 
        dplyr::group_by(line_id) |>
        dplyr::mutate(covar = covar) |>
        data.frame()
    }else{
      # annotate each segment for each barrier line with covar values
      barrier_step_df <- barrier_step_df |> 
        dplyr::group_by(line_id) |>
        dplyr::mutate(covar = covar) |>
        data.frame()
    }
  }
  # if beta is provided, use inverse logistic function to find kappa from
  ## provided (and scaled) covar values
  if(!is.null(beta)){
    if(length(beta) < 2) 
      stop("You need two coefficients (intercept and slope) in beta to simulate 
         covariate effect.")
    beta <- matrix(beta, nrow = 1) |> as.vector()

    # design matrix - extract and scale
    formula <- ~covar
    X.design <- model.matrix(formula, data = 
                               data.frame(covar = barrier_step_df[,"covar"]))
    X.scaled <- scaleX(X.design) # scale X to account for large values
    
    # get untransformed kappa values
    kappa <- exp(X.scaled %*% beta) |> as.vector()
    
    barrier_step_df$kappa <- kappa
    
  } else{ # if not, just set kappa at one value (no covars)
    barrier_step_df$kappa <- kappa
    
    kappa <- barrier_step_df$kappa
  }
  
  if(type == "bcrw"){
    if(is.null(target.location)){
      # make target negative
      boundary_amount_y <- ymax - 0.5*(ymax)
      boundary_amount_x <- xmax - 0.1*(xmax)
      
      target.location <- -(boundary_amount_x) +
        1i*runif(1, -(boundary_amount_y), boundary_amount_y)
    }else{
      target.location <- target.location[,"X"] + 1i*target.location[,"Y"]
    }
  }
  
  if(n.tracks > 1){
    # make list of ids to loop through
    list <- c(1:n.tracks)
    if(parallel){
      n_cores <- round(parallel::detectCores()*0.80) # use 80% of cores
      plan(multisession, workers = n_cores)
      data_list <- future.apply::future_lapply(list, 
                               FUN = simulatePermeability_individual,
                               future.seed = TRUE,
                               barrier_df = barrier_step_df,
                               xmax = xmax, ymax = ymax,
                               N.steps = N.steps, 
                               theta.rho = theta.rho,
                               step.shape = step.shape,
                               step.scale = step.scale,
                               n.null = n.null,
                               type = type,
                               target.location = target.location,
                               plot.track = plot.track, 
                               kappa = kappa)
      plan(sequential)
      
    }else{
      data_list <- lapply(list,
                          function(x)
                            simulatePermeability_individual(
                              barrier_df = barrier_step_df,
                              xmax = xmax, ymax = ymax,
                              N.steps = N.steps, 
                              theta.rho = theta.rho,
                              step.shape = step.shape,
                              step.scale = step.scale,
                              n.null = n.null,
                              type = type,
                              target.location = target.location,
                              plot.track = plot.track, 
                              kappa = kappa))
    }
    
    for(i in c(1:length(data_list))){
      data_list[[i]]$ID <- as.character(paste("SimID",i, sep="_"))
    }
    
    data <- do.call("rbind", data_list)
    
    message(paste(paste("Simulation done.", sum(data$Crossed, na.rm = TRUE)), "steps crossed"))
  }else{
    data <- simulatePermeability_individual(barrier_df = barrier_step_df,
                                             xmax = xmax, ymax = ymax,
                                             N.steps = N.steps, 
                                             theta.rho = theta.rho,
                                             step.shape = step.shape,
                                             step.scale = step.scale,
                                             n.null = n.null,
                                             type = type,
                                             target.location = target.location,
                                             plot.track = plot.track, 
                                             kappa = kappa)
    data$ID <- as.character(paste("SimID",1, sep="_"))
    message(paste(paste("Simulation done.", sum(data$Crossed, na.rm = TRUE)), "steps crossed"))
  }
  
  # change ID to factor
  data$ID <- as.factor(data$ID)
  
  # add barrier id to barrier
  barrier_step_df$barrier.id <- 1:nrow(barrier_step_df)
  
  # add buffer.m as max crossing distance
  buffer.m <- max(subset(data, In.buffer & Crossed)$Dist_toBarrier)
  
  data$In.buffer <- FALSE
  
  data[which(data$Dist_toBarrier <= buffer.m),]$In.buffer <- TRUE
  
  permdata <- list(track = data, 
                   barrier = barrier_step_df)
  
  class(permdata) <- "permdata"  
  
  attr(permdata, "buffer.size") <- buffer.m
  
  attr(permdata, "sample.rate") <- 1
  
  return(permdata)
}

#' @export
simulatePermeability_individual <- function(track.id = 1,
                                            barrier_df,
                                            N.steps, 
                                            theta.rho, step.shape,
                                            step.scale, xmax, ymax,
                                            kappa, 
                                            n.null, 
                                            type = "crw",
                                            target.location = NULL,
                                            plot.track){
  # initialize track simulation:
  ## sampling set of turning angles/steps and start point, Z0
  N <- N.steps
  angles <- suppressWarnings(circular::rwrappedcauchy(N, theta.rho))
  steps <- rweibull(N, step.shape, step.scale)
  
  # if bcrw, keep starting point positive: on x-axis:
  if(type=="bcrw"){
  
    boundary_amount_x <- xmax - 0.1*(xmax)
    Z0 <- (boundary_amount_x) + 1i*runif(1, -ymax, ymax)
  }else{
    Z0 <- runif(1, -xmax, xmax) + 1i*runif(1, -ymax, ymax)
  }
  
  # check for covariate-version of kappa:
  if(is.null(kappa)){
    kappa <- barrier_step_df$kappa
  }
  
  # set empty vector of locs for track, Z
  Z <- rep(NA, N) 
  Z[1] <- Z0 
  Phi <- runif(1, 0, 2*pi) # starting direction
  # extract matrix of X/Y coords
  step.barrier <- cbind(barrier_df[,1], barrier_df[,2])
  # generate track, determining for each sim step if it "bounces"
  ## off barrier or crosses, based on probability kappa
  crossing.index <- c()
  for(i in 1:(N-1)){
    if(i > 1) Phi <- Arg(Z[i] - Z[i-1]) 
    if(type == "bcrw"){
      Z.new.options <- Z[i] + complex(arg = Arg(target.location - Z[i]) + 
                                        sample(angles, n.null), 
                                      mod = sample(steps, n.null))
    }else{
      Z.new.options <- Z[i] + complex(arg = Phi + sample(angles, n.null), 
                                      mod = sample(steps, n.null))
        
    }
    findZnew <- function(Z.new.options, step.barrier, kappa){
      cross_matrix <- sapply(Z.new.options, function(z)
        apply(step.barrier, 1, function(b) isIntersect(c(Z[i], z), b)))
      # determine which barrier segments were crossed by each step
      which_barrier_crossed <- apply(cross_matrix, 2, function(col) {
        row_idx <- which(col)  # Find row indices where TRUE
        if (length(row_idx) > 0) row_idx[1] else NA  # Return first TRUE row or 0
      })
      # find barrier segment with the most crossings
      barrier_seg_crossed <- as.numeric(names(
        which.max(table(which_barrier_crossed))))
      
      if(length(barrier_seg_crossed)==0){
        #if all 0's (no crossings), assign equal probabilities
        Z.new <-  sample(Z.new.options, 1)
        return(Z.new)
      }else{
        # determine how many null steps cross barrier
        samplestep_cross <- cross_matrix |> 
          apply(2, any)
        # find kappa value for barrier seg crossed the most
        kappa_barrier <- rep(kappa[barrier_seg_crossed],
                             length(samplestep_cross))
        # determine null prob of NOT crossing
        ps0 <- prop.table(table(samplestep_cross))[1]
        # prob of crossing for steps that cross: 1-ps0^kappa
        # prob of crossing for steps that don't cross: 1- (1-ps0^kappa)
        prob_step <- ifelse(samplestep_cross, (1-ps0^kappa_barrier),
                            1-(1-ps0^kappa_barrier))
        # sample for Z.new using probability = prob_step
        Z.new <-  sample(Z.new.options, 1, prob = prob_step)
        return(Z.new)
      }
    }
    
    Z.new <-  findZnew(Z.new.options, step.barrier, kappa)
    
    while(Im(Z.new) < -ymax | Im(Z.new) > ymax | 
          Re(Z.new) > xmax | Re(Z.new) < -xmax){
      
      if(type == "bcrw"){
        Z.new.options <- Z[i] + complex(arg = Arg(target.location - Z[i]) + 
                                          sample(angles, n.null), 
                                        mod = sample(steps, n.null))
      }else{
        Z.new.options <- Z[i] + complex(arg = Phi + sample(angles, n.null), 
                                        mod = sample(steps, n.null))
      }
      Z.new <-  findZnew(Z.new.options, step.barrier, kappa)
    }
    
    # determine if step crosses
    crossed <- apply(step.barrier, 1, 
                               function(b) isIntersect(c(Z[i], 
                                                         Z.new), b))
    if(sum(crossed)>0){
     crossing.index <- c(crossing.index, i)
    }
    
    Z[i+1] <- Z.new
  }
  
  if(plot.track){
    z.cross1 <- Z[crossing.index]
    z.cross2 <- Z[crossing.index+1]
    
    plot(Z0, col="white", xlim = c(-xmax,xmax), ylim =c(-ymax,ymax))
    apply(step.barrier, 1, function(x) lines(x, lwd = 2))
    if("covar" %in% colnames(barrier_df)){
      step1 <- step.barrier[,1]
      step2 <- step.barrier[,2]
      covar <- barrier_df$covar
      cols <- (covar - min(covar))/diff(range(covar))
      
      segments(Re(step1), Im(step1), Re(step2), Im(step2), lwd = 3, 
               col = rgb(cols/2, cols, 0))
    }
    if(type == "bcrw"){
      points(target.location, col = "red", cex = 1, pch = 19)
    }
    lines(Z, asp = 1, type = "o", pch= 19, col = rgb(0,0,0,.2), cex = 0.7)
    segments(Re(z.cross1), Im(z.cross1), Re(z.cross2), Im(z.cross2), lwd = 3, 
             col = "orange")
    legend("bottom", legend=c("steps that cross"), 
           fill=c("orange"), cex = 0.7,inset=c(1,1), xpd=TRUE, horiz = TRUE)
  }
  # format track data to be similar to prepPermeability output
  ## add columns for steplength and time (simulated)
  sim_df <- data.frame(Z = Z) |>
    dplyr::mutate(L = c(NA,Mod(diff(Z))),
                  Time = seq.POSIXt(from = as.POSIXct("2025-02-25 00:00:00"), 
                                    by = "1 hour", length.out = length(Z)))
  
  # convert track to step-wise dataframe
  df.step <- sim_df |> 
    plyr::summarize(Z.start = Z[-length(Z)],
                    Z.end = Z[-1],
                    D.time = round(difftime(Time[-1], 
                                            Time[-length(Time)], 
                                            units = "hours")) |> 
                      as.numeric(),
                    Time = Time[-length(Time)]) |> 
    plyr::mutate(Step = Z.end - Z.start,
                 L = Mod(Step), 
                 Phi = Arg(Step), 
                 Theta = c(NA,diff(Phi))) |>
    plyr::mutate(Step.ID = cumsum(rep(1, length(Step))))
  
  # save movement steps to matrix
  step.move <- cbind(start = df.step$Z.start, end = df.step$Z.end)
  
  # determine if steps are within buffer radius of barrier
  z.barrier <- barrier_df[,"Z1"] # use first loc of each barrier seg
  buffer.m <- max(df.step$L, na.rm = TRUE)*2 # buffer - 2x max step L
  dist_starttobarrier <- getClosest(df.step$Z.start, z.barrier)
  dist_endtobarrier <- getClosest(df.step$Z.end, z.barrier)
  df.step$In.buffer <- dist_starttobarrier < buffer.m & dist_endtobarrier < buffer.m
  
  # new column in step df: distance of first loc in step to barrier
  df.step$Dist_toBarrier <- dist_starttobarrier
  
  # check if any steps in buffer 
  if(nrow(subset(df.step, In.buffer))>0){
    
    step.barrier <- cbind(barrier_df[,1], barrier_df[,2])
    step.move <- cbind(start = df.step$Z.start, end = df.step$Z.end)
    
    # set empty crossing column
    df.step$Crossed <- NA
    
    # loop through steps/rows and determine for each if they cross the barrier
    for(i in c(1:nrow(df.step))){
      my.step <- step.move[i,]
      
      crossing <- any(apply(step.barrier, 1, isIntersect, 
                            z2 = c(my.step[1], my.step[2])))
      if(crossing){
        df.step$Crossed[i] <- TRUE
      }else{
        df.step$Crossed[i] <- FALSE
      }
    }
  }else{
    df.step$Crossed <- NA
  }
  return(df.step)
}


#' @export
reflectStep <- function(z1, z2){
  z.intersect <- findSegmentIntersection(z1,z2)
  z1.intersect <- c(z1[1], z.intersect)
  z1.otherside <- c(z.intersect, z1[2])
  arg.reflect <- 2 * Arg(diff(z2)) - Arg(diff(z1))
  c(z.intersect, 
    z.intersect + complex(arg = arg.reflect, mod = Mod(diff(z1.otherside))))
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

#' @export
getClosest <- function(z1, z2){
  D <- outer(z1, z2, function(z1,z2) Mod(z1-z2))
  apply(D, 1, min)
}

#' @export
findSegmentIntersection <- function(z1, z2){
  theta <- Arg(z1[2]-z1[1])
  z1.t <- (z1-z1[1]) * exp(-1i*theta)
  z2.t <- (z2-z1[1]) * exp(-1i*theta)
  intersect <- (prod(Im(z2.t)) < 0) & 
    (any(Re(z2.t) < Re(z1.t[2])) & any(Re(z2.t) > Re(z1.t[1])))
  crossing.t <-  Re(z2.t[1]) - (diff(Re(z2.t))/diff(Im(z2.t))) * Im(z2.t[1])
  crossing.t * complex(mod = 1, arg = theta)  + z1[1]
}