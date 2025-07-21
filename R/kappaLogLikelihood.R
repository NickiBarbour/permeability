#' @export
kappaLogLikelihood_1D <- function(par, data){
  
  crossed <- which(data[,"crossed"] == 1)
  stayed <- which(data[,"crossed"] == 0)
  
  # null probability of STAYING
  P.stayed.null <- data[,"null.stayed"]/(data[,"null.crossed"] + 
                                           data[,"null.stayed"])
  
  # kappa-corrected probability of crossing or staying
  kappa <- exp(par)
  
  P.crossed.kappa <- 1 - P.stayed.null^kappa
  #P.crossed.kappa <- P.crossed.null^kappa
  
  p_cross <- P.crossed.kappa[crossed]
  p_stay <- 1 - P.crossed.kappa[stayed] 
  
  sum(c(log(p_cross), log(p_stay))) # log-likelihood
}


#' @export
kappaLogLikelihood <- function(pars, data, X){
  
  crossed <- which(data[,"crossed"] == 1)
  stayed <- which(data[,"crossed"] == 0)
  
  # null probability of STAYING
  P.stayed.null <- data[,"null.stayed"]/(data[,"null.crossed"] + 
                                           data[,"null.stayed"])
  
  # kappa-corrected probability of crossing or staying
  kappa <- exp(X %*% pars) |> as.vector()
  
  P.crossed.kappa <- 1 - P.stayed.null^kappa
  #P.crossed.kappa <- P.crossed.null^kappa
  
  p_cross <- P.crossed.kappa[crossed]
  p_stay <- 1 - P.crossed.kappa[stayed]      
  
  sum(c(log(p_cross), log(p_stay))) 
}

#' @export
expit <- function(x){
  exp(x)/(1+exp(x))
}
