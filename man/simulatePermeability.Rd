% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulatePermeability.R
\name{simulatePermeability}
\alias{simulatePermeability}
\title{Simulate track(s) with a barrier}
\usage{
simulatePermeability(
  xmax = 50,
  ymax = 50,
  N.steps = 500,
  theta.rho = 0.5,
  step.shape = 5,
  step.scale = 10,
  barrier = NULL,
  n.segments = NULL,
  n.null = 30,
  kappa = NULL,
  beta = NULL,
  covar = NULL,
  n.tracks = 1,
  type = "crw",
  target.location = NULL,
  plot.track = FALSE,
  parallel = FALSE
)
}
\arguments{
\item{xmax}{maximum x-value of bounding box}

\item{ymax}{maximum y-value of bounding box}

\item{N.steps}{number of steps to simulate (length of simulated track)}

\item{theta.rho}{turning angle concentration parameter, to simulate track
turning angles using a wrapped Cauchy distribution}

\item{step.shape}{shape parameter to simulate track step lengths, using a
Weibull distribution}

\item{step.scale}{scale parameter to simulate track step lengths, using a
Weibull distribution}

\item{barrier}{three column matrix or data frame of barrier vertices
(named \code{X} and \code{Y}) and line ID (named \code{line_id} - should be equal to 1 if
there is only one line segment in the barrier)}

\item{n.segments}{if \code{barrier = NULL}, number of segments for simulated
barriers  - for a barrier with multiple lines, assumes an equal
number of segments per line.}

\item{n.null}{number of null steps to simulate around each actual step
to estimate the null probability of crossing. Default 30.}

\item{kappa}{desired permeability value of the barrier -
values closer to 0 represent a less permeable barrier and values
closer to 1 represent a more permeable barrier.  Should be \code{NULL} if using
a covariate and \code{beta}..}

\item{beta}{a vector of two coefficients (intercept and main effect) for the
impact of the covariate on permeability - default \code{NULL} (no covariate).}

\item{covar}{a numeric vector of covariate values for each barrier segment,
equal in length to \code{n.segments} - default \code{NULL} (if \code{NULL}, function will
simulate covar values on a small-to-big numeric scale)}

\item{n.tracks}{number of tracks to simulate - default 1}

\item{type}{whether to simulate a correlated random walk (\code{crw}) or biased
correlated random walk (\code{bcrw}) - if the latter, will generate a random
target location, unless \code{target.location} is provided. Default type is \code{crw}.}

\item{target.location}{target location (two column matrix of \code{X}, \code{Y}
locations) for biased random walk, if \code{type = "bcrw"}.}

\item{plot.track}{whether to plot resulting track - default \code{FALSE}.}

\item{parallel}{if generating multiple tracks, whether to run in parallel or
not (recommended). Default \code{FALSE}.}
}
\value{
\code{permdata} object, which is a list containing two "stepwise" data
frames object. The \code{track} object has rows corresponding to each observed
step in the data and columns for: \code{ID}, \code{Time}, start/end location of each
step (complex, \code{Z.start}, \code{Z.end}), the movement step (\code{Step}, complex, start
and locations of step), the time step (in hours, \code{D.time}), whether the step
was in the buffer (\code{In.buffer}, TRUE/FALSE), movement metrics
(step length, m - \code{L}, absolute angle, radians - \code{Phi}, relative turning
angle, radians - \code{Theta}), and whether the step crossed the barrier or not
(\code{Crossed} - TRUE/FALSE). The \code{barrier} object has rows
corresponding to each barrier segment and columns for the complex locations
of the start and end of the barrier segments (\code{Z1} and \code{Z2}), the line id
(\code{line_id}), covariate values (\code{covar}), the
id of each barrier segment (\code{barrier.id}), and "true" kappa values
(\code{kappa}).
}
\description{
Simulate track(s), where each track is either a correlated
random walk (default) or biased random walk, with a set permeability with
respect to a linear barrier and the option to have permeability vary
as a function of covariate values associated with barrier segments.
}
\examples{
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

}
