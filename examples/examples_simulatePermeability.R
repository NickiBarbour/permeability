require(permeability)

# BASIC EXAMPLE #################

# simulate one track (correlated random walk)
sim_track <- simulatePermeability(xmax = 50, ymax = 50, 
                                  N.steps = 1000, theta.rho = 0.5, 
                                  step.shape = 5, step.scale = 10, 
                                  n.segments = 20,
                                  kappa = 0.7,
                                  plot.track = TRUE) 

# build crossing table
sim_cT <- buildCrossingTable(permdata = sim_track)
# fit null model
null_model <- fitPermeability(data = sim_cT)

null_model

# OTHER EXAMPLES ################

# multiple tracks-----------
library(future.apply)
sim_tracks <- simulatePermeability(xmax = 50, ymax = 50,  
                                   N.steps = 1000, theta.rho = 0.5, 
                                   step.shape = 5, step.scale = 10, 
                                   kappa = 0.2,
                                   n.tracks = 30,
                                   n.segments = 20,
                                   plot.track = FALSE,
                                   parallel = TRUE) 

# build crossing table
sim_cT <- buildCrossingTable(permdata = sim_tracks, 
                             parallel = TRUE)

# estimate kappa (null model)
null_model <- fitPermeability(data = sim_cT)

null_model

# covariate example (simulated covariate data) ---------------
## one track
sim_track <- simulatePermeability(xmax = 50, ymax = 50,  
                                  N.steps = 1000, theta.rho = 0.5, 
                                  step.shape = 5, step.scale = 10, 
                                  n.segments = 20,
                                  beta = c(-2,1),
                                  plot.track = TRUE) 

sim_cT <- buildCrossingTable(permdata = sim_track)

sim_cT <- merge(sim_cT, sim_track$barrier[,c("barrier.id","covar")],
                by = "barrier.id", all.x = TRUE)

covar_model <- fitPermeability(~ covar, data = sim_cT)

covar_model$coef.table

# covariate example (user-provided) -----------------
## multiple tracks
n.segments <- 20
xmax <- ymax <- 50

barrier <- cbind(X = xmax * rnorm(n.segments+1, sd = xmax / 5e3),
                 Y = ymax * seq(-1,1,length = n.segments+1),
                 line_id = 1)

covar_values <- sort(runif(n.segments, 1, 30))

library(future.apply)
sim_tracks <- simulatePermeability(xmax = 50, ymax = 50,  
                                   N.steps = 1000, theta.rho = 0.5, 
                                   step.shape = 5, step.scale = 10, 
                                   barrier = barrier,
                                   covar = covar_values,
                                   beta = c(-2,1),
                                   n.tracks = 30,
                                   n.null = 60,
                                   plot.track = FALSE,
                                   parallel = TRUE) 
## build crossing table
sim_cT <- buildCrossingTable(permdata = sim_tracks,
                             parallel = TRUE)

## add covariate values to crossing table
sim_cT <- merge(sim_cT, sim_tracks$barrier[,c("barrier.id","covar")],
                by = "barrier.id", all.x = TRUE)

## estimate kappa (covar model)
covar_model <- fitPermeability(~ covar, data = sim_cT)

covar_model$coef.table

## predict and plot
predict_kappa <- predict(covar_model, se.fit = TRUE)
predict_kappa$covar <- covar_model$data$covar

predict_kappa <- predict_kappa |>
  dplyr::arrange(kappa.hat)

# true kappa values
predict_kappa <- merge(predict_kappa, 
                       sim_tracks$barrier[,c("kappa","covar")], by = "covar")

plot(predict_kappa$covar, predict_kappa$kappa.hat, type = "l", ylim = c(0,1),
     xlab = "Covariate Values", ylab = "Estimated Kappa Values", col ="orange", lwd = 2)
lines(predict_kappa$covar, predict_kappa$kappa, col = "red", lwd = 2)
lines(predict_kappa$covar, predict_kappa$ci.high, col = "grey")
lines(predict_kappa$covar, predict_kappa$ci.low, col = "grey")

# simulate one track (biased random walk with target location) ----------
target <- data.frame(X = -40, Y = 20) |>
  as.matrix()

sim_track <- simulatePermeability(xmax = 50, ymax = 50,  
                                  N.steps = 500, theta.rho = 0.5, 
                                  step.shape = 5, step.scale = 10, 
                                  n.segments = 10,
                                  kappa = 0.2,
                                  type = "bcrw",
                                  target.location = target,
                                  plot.track = TRUE) 

# user-provided buffer ("grid" example, with multiple lines) -----------
xmax <- 50 
ymax <- 50
n.segments <- 10 

barrier_sep <- (1/2)*xmax # spacing between linear features
y <- ymax-barrier_sep # allow spacing around barrier lines

barrier1 <- cbind(X = (xmax-barrier_sep),
                  Y = y * seq(-1,1,length = n.segments+1),
                  line_id = 1)
barrier2 <- cbind(X = (xmax-barrier_sep*2),
                  Y = y * seq(-1,1,length = n.segments+1),
                  line_id = 2)
barrier3 <- cbind(X = (xmax-barrier_sep*3),
                  Y = y * seq(-1,1,length = n.segments+1),
                  line_id = 3)

grid.barrier <- rbind(barrier1, barrier2, barrier3)

sim_track <- simulatePermeability(xmax = 50, ymax = 50, 
                                  N.steps = 500, theta.rho = 0.5, 
                                  step.shape = 5, step.scale = 10, 
                                  kappa = 0.2,
                                  barrier = grid.barrier,
                                  plot.track = TRUE) 

