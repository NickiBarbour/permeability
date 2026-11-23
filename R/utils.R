# Internal helper functions -------------------------------------------------

# distance from each point in z1 to the closest point in z2
getClosest <- function(z1, z2){
  D <- outer(z1, z2, function(z1, z2) Mod(z1 - z2))
  apply(D, 1, min)
}

# does segment z2 intersect segment z1?
isIntersect <- function(z1, z2){
  theta <- Arg(z1[2] - z1[1])
  z1.t <- (z1 - z1[1]) * exp(-1i * theta)
  z2.t <- (z2 - z1[1]) * exp(-1i * theta)
  slope.r <- (Re(z2.t[2]) - Re(z2.t[1])) / (Im(z2.t[2]) - Im(z2.t[1]))
  x.intercept <- Re(z2.t[1]) - slope.r * Im(z2.t[1])
  intersect.x <- x.intercept > 0 & x.intercept < Re(z1.t[2])
  intersect.x & (Im(z2.t[2]) * Im(z2.t[1])) < 0
}

# parallel processing helpers -------------------------------------------------

# stop with a clear message if the (optional) future packages are missing
check_future <- function(){
  if(!requireNamespace("future", quietly = TRUE) ||
     !requireNamespace("future.apply", quietly = TRUE)){
    stop("Packages 'future' and 'future.apply' are required when ",
         "parallel = TRUE. Please install them.", call. = FALSE)
  }
}

# number of workers to use: user-specified, or `frac` of available cores.
# Never more than 2 when running R CMD check (CRAN policy).
set_workers <- function(frac, n.cores = NULL){
  if(is.null(n.cores)){
    n.cores <- max(1L, round(parallel::detectCores() * frac))
  }
  if(identical(Sys.getenv("_R_CHECK_LIMIT_CORES_"), "TRUE")){
    n.cores <- min(n.cores, 2L)
  }
  n.cores
}
