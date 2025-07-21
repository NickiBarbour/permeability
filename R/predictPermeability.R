#' Prediction method for permeability fit 
#' 
#' @param fit permeability fit object, output of \link{`fitPermeability`}
#' @param new.data an optional data frame if predictions are needed over new 
#' values.
#' @param conf.level confidence level (default 95%) .
#' @param bootstrap whether to get standard errors on estimates using 
#' bootstrapping. Recommended if there are issues getting standard errors with
#' the Hessian with \link{`fitPermeability`}. 
#' @param bootstrap.n number of iterations to run the bootstrapping for, if 
#' selected. Defaults to 1000.
#' @param parallel if bootstrap = TRUE, whether to run in parallel using the
#' future.apply package or not. Recommended (bootstrapping can take 10+ min).
#' @returns data frame with columns `kappa.hat` - point estimate of 
#' permeability, 
#' `ci.low` and `ci.high` - confidence intervals around the kappa estimate
#' @export
predict.permfit <- function(fit, new.data = NULL, se.fit = TRUE, 
                            conf.level = 0.95, bootstrap = FALSE, 
                            bootstrap.n = 1000, parallel = FALSE) {
  formula <- fit$formula
  if(se.fit){
    beta.hat <- suppressWarnings(data.frame(estimate = fit$optim.output$par, 
                           se = sqrt(diag(solve(-fit$optim.output$hessian)))))
    if(sum(is.na(beta.hat$se))> 0 & bootstrap==FALSE){
      warning("Unable to get standard errors on some estimates. You may want to use bootstrap = TRUE.")
    }
  }else{
    beta.hat <- data.frame(estimate = fit$optim.output$par)
  }
  Z <- abs(qnorm((1 - conf.level)/2))
  # check for covariates
  if(length(beta.hat$estimate) > 1){
    allterms <- attributes(terms(formula))$term.labels
    which_spline <- sapply(allterms, function(t) grepl("s(", t, fixed = TRUE))
    
    if(!is.null(new.data)){
      model_data <- new.data
    } else{
      model_data <- fit$data
    }
    if(any(which_spline)){ # check for splined covariate
      spline_text <- names(which_spline)[which_spline]
      spline_formula <- as.formula(call("~", 
                                        parse(text = spline_text)[[1]]))[[2]]
      fit.spline <- smoothCon(eval(spline_formula), 
                              data = model_data)[[1]]$X 
      
      term <- eval(spline_formula)$term
      colnames(fit.spline) <- paste0("s(",term,")", 1:ncol(fit.spline))
      
      if(any(!which_spline)){
        covars <-  names(which_spline)[!which_spline]
        
        fixed_formula <- eval(parse(text = 
                                      paste("formula(~", 
                                            paste(covars, 
                                                  collapse = " + "), ")")))
        fit.fixed <- model.matrix(fixed_formula, data = model_data)
        fit.raw <- cbind(fit.fixed, fit.spline)
        
      }else fit.raw <- cbind(Intercept = rep(1, nrow(fit.spline)), 
                              fit.spline)
    } else fit.raw <- model.matrix(formula, data = model_data)
    
    if (any(is.na(fit.raw))) stop("NA values detected in model matrix. Check covariate data for NA values.")
    
    fit.scaled <- scaleX(fit.raw)
    
    if(se.fit){
      if(bootstrap){
        # function to run bootstrap for the ith iteration
        runBootstrap <- function(i, model_formula, og_data, new_data = NULL){
          if(is.null(new_data)){
            bs_data <- og_data[sample(1:nrow(og_data), replace = TRUE), ]
            myfit <- suppressWarnings(fitPermeability(model_formula, 
                                                      data = bs_data))
            predict_df <- suppressWarnings(predict(myfit, se.fit = FALSE,
                                                   new.data = og_data))
            
            return(predict_df)
          }else{
            bs_data <- og_data[sample(1:nrow(og_data), replace = TRUE), ]
            myfit <- suppressWarnings(fitPermeability(model_formula,
                                                      data = bs_data))
            predict_df <- suppressWarnings(predict(myfit, new.data = new_data, 
                                                   se.fit = FALSE))
            
            return(predict_df)
          }
        }
        iterations <- bootstrap.n
        if(parallel){
          n_cores <- round(parallel::detectCores()*0.80)
          message(paste("Running bootstrap using", iterations, "iterations and",
                        n_cores,"cores (80% of available).", sep=" "))
          
          if(is.null(new.data)){
            plan(multisession, workers = n_cores)
            predict.bs <- future_sapply(
              1:iterations,
              function(i) runBootstrap(i, 
                                       model_formula = formula,
                                       og_data = fit$data), future.seed = TRUE)
            plan(sequential)
          }else{
            plan(multisession, workers = n_cores)
            predict.bs <- future_sapply(
              1:iterations,
              function(i) runBootstrap(i, 
                                       model_formula = formula,
                                       og_data = fit$data,
                                       new_data = new.data), future.seed = TRUE)
            plan(sequential)
          }
        }else{
          warning(paste("Running bootstrap with", iterations,"iterations and without parallel processing. Adjust your expectations accordingly.", sep=" "))
          
          if(is.null(new.data)){
            predict.bs <- sapply(
              1:iterations,
              function(i) runBootstrap(i, 
                                       model_formula = formula,
                                       og_data = fit$data))
          }else{
            predict.bs <- sapply(
              1:iterations,
              function(i) runBootstrap(i, 
                                       model_formula = formula,
                                       og_data = fit$data,
                                       new_data = new.data))
          }
        }
        # get kappa est
        Y.hat <- fit.scaled %*% beta.hat$estimate
        kappa_hat <- exp(Y.hat) 
        predict.bs2 <- do.call("cbind", predict.bs)
        predict.bs2[is.infinite(predict.bs2)] <- NA
        
        return(data.frame(kappa.hat = kappa_hat,
                          ci.low =  
                            apply(predict.bs2, 1, 
                                  function(x) quantile(x,
                                                       probs = (1-conf.level)/2, 
                                                       na.rm = TRUE)),
                          ci.high = 
                            apply(predict.bs2, 1, 
                                  function(x) quantile(x,
                                                       probs = 
                                                         conf.level + 
                                                         ((1-conf.level)/2),
                                                       na.rm = TRUE))))
        
        return(data.frame(kappa.hat = kappa_hat,
                   ci.low =  apply(predict.bs, 1, quantile,
                                (1-conf.level)/2),
                   ci.high = apply(predict.bs, 1, quantile,
                                conf.level + ((1-conf.level)/2))))
      }else{
        if(sum(is.na(beta.hat$se))> 0) warning("Unable to get standard errors on some estimates. Predictions may be unreliable.")
        hessian <- fit$optim.output$hessian
        covariance <- solve(-hessian)
        
        Y.hat <- fit.scaled %*% beta.hat$estimate
        
        kappa.hat <- exp(Y.hat)
        se <- sqrt(rowSums((fit.scaled %*% covariance) * fit.scaled))
        
        ci.low <- exp(Y.hat - Z * se)
        ci.high <- exp(Y.hat + Z * se)
        
        return(data.frame(kappa.hat, ci.low, ci.high))
      }
    }else{
      Y.hat <- fit.scaled %*% beta.hat$estimate
      
      kappa.hat <- exp(Y.hat)
      
      return(data.frame(kappa.hat))
    }
    
  }else{ # null model
    if(se.fit){
      return(data.frame(
        kappa.hat = exp(beta.hat$estimate),
        ci.low = exp(beta.hat$estimate - Z * beta.hat$se),
        ci.high = exp(beta.hat$estimate + Z * beta.hat$se)
      ))
    }else{
      return(data.frame(
        kappa.hat = exp(beta.hat$estimate)
      ))
    }
  }
}



#' @export
scaleX <- function(fit) {
  # Identify the intercept column (usually the first column in 
  ## model.matrix output)
  intercept.col <- grepl("Intercept", colnames(fit))
  
  # Identify columns to scale (all except the intercept)
  cols.to.scale <- setdiff(seq_along(fit[1, ]), intercept.col)
  
  # Compute means and standard deviations
  fit.means <- colMeans(fit[, cols.to.scale, drop = FALSE], na.rm = TRUE)
  fit.sds <- apply(fit[, cols.to.scale, drop = FALSE], 2, sd, na.rm = TRUE)
  
  # Avoid division by zero for constant columns
  fit.sds[fit.sds == 0] <- 1
  
  # Standardize the data
  fit.scaled <- fit
  fit.scaled[, cols.to.scale] <- sweep(fit.scaled[, 
                                                  cols.to.scale, 
                                                  drop = FALSE], 2, 
                                       fit.means, "-")
  fit.scaled[, cols.to.scale] <- sweep(fit.scaled[, 
                                                  cols.to.scale, 
                                                  drop = FALSE], 2, 
                                       fit.sds, "/")
  
  # Store attributes for potential back-transformation
  attr(fit.scaled, "cols.to.scale") <- cols.to.scale
  attr(fit.scaled, "mus") <- fit.means
  attr(fit.scaled, "sds") <- fit.sds
  
  return(fit.scaled)
}

