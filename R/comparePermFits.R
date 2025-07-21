#' compare Permeability Model Fits
#' @name comparePermFits
#' @description Given a list of models or model formulae, compares model
#' fits using AIC values. Output is ranked by the lowest AIC to highest
#' AIC model.
#'
#' @param formulae list of formulae, with the structure "~ covariate..."
#' @param models optional - a list of fitted models, as from
#'  \link{fitPermeability}
#' @param data crossing table, as output from \link{buildCrossingTable}
#' @return A dataframe with rows for each model, columns for AIC,
#' the delta AIC ("dAIC" or difference in AIC between consecutive models),
#' number of parameters (`k`), log Likelihood value, and 
#' sorted by AIC values.
#' @example examples/examples_comparePermFits.R
#' @export

comparePermFits <- function(formulae, models = NULL, data){

  if(is.null(models)){
    model_fits <- 
      lapply(formulae, 
             function(x) suppressWarnings(fitPermeability(x, data = data)))
  }else{
    model_fits <- models
  }
  
  plyr::ldply(model_fits, 
              function(fit)
                data.frame(k = nrow(fit$coef.table), 
                           logLik = fit$logLik,
                           AIC = fit$AIC) 
              ) |> 
    plyr::mutate(dAIC = AIC - min(AIC)) |> 
    plyr::rename(c(.id = "Model")) |> 
    plyr::arrange(dAIC)

}  


