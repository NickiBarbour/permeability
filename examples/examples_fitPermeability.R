require(permeability)
data("example_cT")

# NULL MODEL ###########################################
null_model <- fitPermeability(data = example_cT)
print(null_model)

# COVARIATE MODEL ##################################

## fixed effect
traffic_model <- fitPermeability(~ traffic_volume, data = example_cT)

## get predictions with 95% CI using bootstrapping
predict_kappa <- predict(traffic_model, se.fit = TRUE)
predict_kappa$traffic_volume <- traffic_model$data$traffic_volume
predict_kappa <- predict_kappa |>
  dplyr::arrange(kappa.hat)

plot(predict_kappa$traffic_volume, predict_kappa$kappa.hat, type = "l", 
     xlab = "Daily Traffic Volume", 
     ylab = "Estimated Kappa Values", 
     col ="orange", lwd = 2)

## day of year smoother
example_cT$doy <- lubridate::yday(example_cT$Date)
doy_model <- fitPermeability(~ s(doy, k = 6, bs = "cc"), data = example_cT)
doy_model

predict_kappa <- predict(doy_model, se.fit = TRUE, bootstrap = TRUE,
                         new.data = data.frame(doy = seq(min(example_cT$doy),
                                                         max(example_cT$doy),
                                                         by = 1)))
predict_kappa$doy <- seq(min(example_cT$doy), max(example_cT$doy),
                         by = 1)
