library(permeability)
data("example_cT")

# NULL MODEL ###########################################
null_model <- fitPermeability(data = example_cT)
print(null_model)

# COVARIATE MODEL ##################################

## fixed effect
traffic_model <- fitPermeability(~ traffic_volume, data = example_cT)

## get predictions with 95% CI
predict_kappa <- predict(traffic_model, se.fit = TRUE)
predict_kappa$traffic_volume <- traffic_model$data$traffic_volume
predict_kappa <- predict_kappa[order(predict_kappa$kappa.hat), ]

plot(predict_kappa$traffic_volume, predict_kappa$kappa.hat, type = "l", 
     xlab = "Daily Traffic Volume", 
     ylab = "Estimated Kappa Values", 
     col ="orange", lwd = 2)

## day of year smoother
example_cT$doy <- as.numeric(format(example_cT$Date, "%j"))
doy_model <- fitPermeability(~ s(doy, k = 6, bs = "cc"), data = example_cT)
doy_model

\donttest{
## bootstrapped confidence intervals (use bootstrap.n = 1000 or more for 
## real analyses; reduced here for speed)
doy_values <- data.frame(doy = seq(min(example_cT$doy), 
                                   max(example_cT$doy), by = 1))
predict_kappa <- predict(doy_model, se.fit = TRUE, bootstrap = TRUE,
                         bootstrap.n = 50, new.data = doy_values)
predict_kappa$doy <- doy_values$doy
}
