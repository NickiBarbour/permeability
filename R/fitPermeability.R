#' Estimate permeability 
#' @description Estimates permeability using the crossing table as output from
#' \link{buildCrossingTable}. If no formula is provided, will default to a
#' null model of permeability, returning one value for kappa. Otherwise, if 
#' covariate values are provided in the `formula`, will fit a model that 
#' estimates permeability as a logit-linear function of covariate(s), returning
#' a kappa value for each unique covariate value. Models are fit using
#' maximum likelihood estimation, with the `optim` function.
#' @param formula one-sided formula - uses names of columns in the crossing 
#' table as possible covariates. Defaults to NULL (null model).
#' @param data crossing table object - output of
#'  \code{\link{build_crossingTable}}.
#' @param beta0 optional parameter for initial values for model to be optimized
#' over. Default NULL (see `optim` function for more details).
#' @example examples/examples_fitPermeability.R
#' 

#' @export
fitPermeability <- function(formula = ~1, data, beta0 = NULL){
  # base estimation of kappa
  if(formula == ~1){
    pars.fit <- suppressWarnings(optim(0, kappaLogLikelihood_1D, data = data,
                      control = list(fnscale = -1), hessian = TRUE))
    if(pars.fit$hessian[,1]==0){
      se <- NA
      coef <-  data.frame(estimate = pars.fit$par, 
                          se = se)
      message("Unable to get standard errors on estimate - estimate may be too close to 1.")
    } else {
      se <- suppressWarnings(diag(sqrt(solve(-pars.fit$hessian))))
      coef <-  data.frame(estimate = pars.fit$par, 
                          se = se)
    }
    coef$z_score = coef$estimate/coef$se 
    coef$p_value <- (1 - pnorm(abs(coef$z_score))) * 2
    
    
    logLik <- pars.fit$value
    AIC <- -2*logLik + 2*nrow(coef)
    
    permfit <- list(coef.table = coef, 
                    logLik = logLik, AIC = AIC, 
                    formula = formula, 
                    optim.output = pars.fit,
                    data = data)
    
    class(permfit) <- "permfit"      
    return(permfit)
  } else {
    
    allterms <- attributes(terms(formula))$term.labels
    which_spline <- sapply(allterms, function(t) grepl("s(", t, fixed = TRUE))
    if(sum(which_spline) > 1) 
      warning("permeability does not currently support more than one splined covariate.")
    
    if(any(which_spline)){
      spline_text <- names(which_spline)[which_spline]
      spline_formula <- as.formula(call("~", parse(text = 
                                                     spline_text)[[1]]))[[2]]
      X.spline <- smoothCon(eval(spline_formula), 
                            data = data)[[1]]$X 
      
      term <- eval(spline_formula)$term
      colnames(X.spline) <- paste0("s(",term,")", 1:ncol(X.spline))
      
      if(any(!which_spline)){
        covars <-  names(which_spline)[!which_spline]
        
        fixed_formula <- eval(parse(text = 
                                      paste("formula(~", 
                                            paste(covars, 
                                                  collapse = " + "), ")")))
        X.fixed <- model.matrix(fixed_formula, data = data)
        X.raw <- cbind(X.fixed, X.spline)
        
      } else X.raw <- cbind(Intercept = rep(1, nrow(X.spline)), 
                            X.spline)
    } else X.raw <- model.matrix(formula, data = data)
    
    # check for NA values
    if(any(is.na(X.raw))) stop("NA values detected in model matrix. Check covariate data for NA values.")
    
    # scale design matrix (accounts for large values)
    X.scaled <- scaleX(X.raw)
    beta_init <- X.scaled[1,] * 0
    
    if(!is.null(beta0)) beta_init[1:length(beta0)] <- beta0
    pars.fit <- suppressWarnings(optim(beta_init, kappaLogLikelihood, 
                      data = data, X = X.scaled,
                      control=list(fnscale=-1),
                      hessian = TRUE)) 
    
    se_try <- suppressWarnings(try(diag(sqrt(solve(-pars.fit$hessian))), 
                                   silent = TRUE))
    
    if(!inherits(se_try, "try-error")){
      se <- suppressWarnings(diag(sqrt(solve(-pars.fit$hessian))))
      coef <-  data.frame(estimate = pars.fit$par, se = se)
      if(sum(is.na(coef$se)>0)){
        warning("Unable to get standard errors on some estimates. It's suggested to get use bootstrapping with the predict function if you want standard errors.")
      }
    } else {
      se <- NA
      coef <-  data.frame(estimate = pars.fit$par, 
                          se = se)
      warning("Unable to get standard errors on estimates. It's suggested to get use bootstrapping with the predict function if you want standard errors.")
    }
  }
  
  coef$z_score = coef$estimate/coef$se 
  coef$p_value <- (1 - pnorm(abs(coef$z_score))) * 2

  
  logLik <- pars.fit$value
  AIC <- -2*logLik + 2*nrow(coef)
  
  permfit <- list(coef.table = coef, 
                  logLik = logLik, AIC = AIC, 
                  formula = formula, 
                  optim.output = pars.fit,
                  data = data)
  
  class(permfit) <- "permfit"      
  return(permfit)
}


#' @export
print.permfit <- function(x){
  if(nrow(x$coef.table) == 1){
    kappa.hat <- predict(x)
    cat("\nEstimate of kappa (with 95% confidence intervals):\n")
    print(kappa.hat)
  }else{
    # check for smoother terms
    allterms <- attributes(terms(x$formula))$term.labels
    which_spline <- sapply(allterms, function(t) grepl("s(", t, fixed = TRUE))
    if(any(which_spline)){
      # find linear vs nonlinear terms
      nonspline_text <- names(!which_spline)[!which_spline]
      spline_text <- names(which_spline)[which_spline]
      # check for no linear terms:
      if(length(nonspline_text)==0){
        spline_formula <- as.formula(call("~", 
                                          parse(text = spline_text)[[1]]))[[2]]
        spline_term <- eval(spline_formula)$term
        intercept <- x$coef.table[grepl("Intercept", rownames(x$coef.table)),]
        spline_coef <- x$coef.table[grepl(spline_term, rownames(x$coef.table)),]
        
        # fit null model to compare to
        null_model <- suppressWarnings(fitPermeability(~ 1, 
                                      data = x$data))
        # find deviance for each model
        Dev0 <- -2*null_model$logLik
        Dev1 <- -2*x$logLik
        # find the likelihood ratio test statistic
        LRT <- Dev0 - Dev1
        # find the # of parameters
        k0 <- nrow(null_model$coef.table)
        k1 <- nrow(x$coef.table)
        k <- k1-k0
        
        p_value <- 1-pchisq(LRT, k)
        
        # create dataframe of results
        smoother_results <- data.frame(k = k, LRT = LRT, p_value = p_value)
        rownames(smoother_results) <- paste("s(",spline_term,")", sep = "")
        
        cat("Permeability fit call:", paste(x$formula, collapse = ""), "\n\n")
        cat("Linear coefficients:", "\n\n")
        print(intercept)
        cat("\nApproximate significance of smooth terms:", "\n\n")
        print(smoother_results)
        cat("\nlogLik:", x$logLik, 
            "\n   AIC:", x$AIC, "\n")
      }else{
        nonspline_formula <- as.formula(call("~", 
                                             parse(text = 
                                                     nonspline_text)[[1]]))[[2]]
        spline_formula <- as.formula(call("~", 
                                          parse(text = spline_text)[[1]]))[[2]]
        spline_term <- eval(spline_formula)$term
        intercept <- x$coef.table[grepl("Intercept", rownames(x$coef.table)),]
        linear_coef <- x$coef.table[grepl(nonspline_text, 
                                          rownames(x$coef.table)), ]
        spline_coef <- x$coef.table[grepl(spline_term, rownames(x$coef.table)),]
        
        # fit null/linear model to compare to
        nonspline_formula <- as.formula(paste("~", nonspline_formula))
        null_model <- suppressWarnings(fitPermeability(nonspline_formula, 
                                      data = x$data))
        # find deviance for each model
        Dev0 <- -2*null_model$logLik
        Dev1 <- -2*x$logLik
        # find the likelihood ratio test statistic
        LRT <- Dev0 - Dev1
        # find the # of parameters
        k0 <- nrow(null_model$coef.table)
        k1 <- nrow(x$coef.table)
        k <- k1-k0
        
        p_value <- 1-pchisq(LRT, k)
        
        # create dataframe of results
        smoother_results <- data.frame(k = k, LRT = LRT, p_value = p_value)
        rownames(smoother_results) <- paste("s(",spline_term,")", sep = "")
        
        cat("Permeability fit call:", paste(x$formula, collapse = ""), "\n\n")
        cat("Linear coefficients:", "\n\n")
        print(rbind(intercept, linear_coef))
        cat("\nApproximate significance of smooth terms:", "\n\n")
        print(smoother_results)
        cat("\nlogLik:", x$logLik, 
            "\n   AIC:", x$AIC, "\n")
      }
    }else{
      cat("Permeability fit call:", paste(x$formula, collapse = ""), "\n\n")
      print(x$coef.table)
      cat("\nlogLik:", x$logLik, 
          "\n   AIC:", x$AIC, "\n")
    }
  }
}

#' @export
summary.permfit <- function(x){
  list(coef = x$coef.table, logLik = x$logLik, AIC = x$AIC)
}

#' @export
logLik.permfit <- function(x){
  x$logLik
}


#' @export
AIC.permfit <- function(x){
  x$AIC
}



