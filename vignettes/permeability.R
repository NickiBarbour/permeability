## ----setup, include = FALSE---------------------------------------------------
require(permeability)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = ">",
  warning = FALSE,
  message = TRUE
)

## ----package, eval = FALSE----------------------------------------------------
# library(devtools)
# install_github("pathprovidedonpublication", vignettes = TRUE)
# library(permeability)

## ----simtrack1, fig.align = "center", fig.width = 6, fig.height=6, eval = FALSE----
# sim_track <- simulatePermeability(xmax = 30, ymax = 30,
#                                   N.steps = 500, theta.rho = 0.5,
#                                   step.shape = 5, step.scale = 10,
#                                   n.segments = 10,
#                                   kappa = 0.2,
#                                   plot.track = TRUE)

## ----loadsimdata, include= FALSE----------------------------------------------
#save(sim_track, file = "./data/sim_example_1.rda")
load("./data/sim_example_1.rda")

## ----simclass-----------------------------------------------------------------
class(sim_track)

## ----simstr-------------------------------------------------------------------
str(sim_track)

## ----simtrackhead-------------------------------------------------------------
head(sim_track$track)

## ----simbarrierhead-----------------------------------------------------------
head(sim_track$barrier)

## ----target-------------------------------------------------------------------
target <- data.frame(X = -40, Y = 20) |>
  as.matrix()

## ----simtrack2, fig.align = "center", fig.width = 6, fig.height=6, eval = FALSE----
# sim_track <- simulatePermeability(xmax = 50, ymax = 50,
#                                   N.steps = 100, theta.rho = 0.5,
#                                   step.shape = 5, step.scale = 10,
#                                   n.segments = 10,
#                                   kappa = 0.2,
#                                   type = "bcrw",
#                                   target.location = target,
#                                   plot.track = TRUE)

## ----gridbarrier--------------------------------------------------------------
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

## ----plotgridbarrier, fig.align = "center", fig.width = 6, fig.height=6, eval = FALSE----
# plot(barrier1, type= "l", ylim = c(-ymax,ymax), xlim = c(-xmax,xmax))
# lines(barrier2)
# lines(barrier3)

## ----headgridbarrier----------------------------------------------------------
grid.barrier <- rbind(barrier1, barrier2, barrier3)

head(grid.barrier)

## ----simtrack3, fig.align = "center", fig.width = 6, fig.height=6, eval = FALSE----
# sim_track <- simulatePermeability(xmax = 50, ymax = 50,
#                               N.steps = 500, theta.rho = 0.5,
#                               step.shape = 5, step.scale = 10,
#                               kappa = 0.2,
#                               barrier = grid.barrier,
#                               plot.track = TRUE)

## ----loadsimdata2, include= FALSE---------------------------------------------
#save(sim_track, file = "./data/sim_example_2.rda")
load("./data/sim_example_2.rda")

## ----simct, message = FALSE, results='hide'-----------------------------------
sim_cT <- buildCrossingTable(permdata = sim_track)

## ----simctstr-----------------------------------------------------------------
str(sim_cT)

## ----simnullmodel, warning = FALSE--------------------------------------------
null_model <- fitPermeability(data = sim_cT)

null_model

## ----summarysimnullmodel------------------------------------------------------
summary(null_model)

## ----predictsimnullmodel------------------------------------------------------
predict(null_model)

## ----simtrack4, eval = FALSE--------------------------------------------------
# library(future.apply)
# 
# sim_tracks <- simulatePermeability(xmax = 50, ymax = 50,
#                               N.steps = 500, theta.rho = 0.5,
#                               step.shape = 5, step.scale = 10,
#                               kappa = 0.2,
#                               barrier = grid.barrier,
#                               n.tracks = 50,
#                               plot.track = FALSE,
#                               parallel = TRUE)

## ----simcovar, fig.align = "center", fig.width = 6, fig.height=6, eval = FALSE----
# sim_track <- simulatePermeability(xmax = 50, ymax = 50,
#                               N.steps = 500, theta.rho = 0.5,
#                               step.shape = 5, step.scale = 10,
#                               n.segments = 15,
#                               beta = c(-2,1),
#                               plot.track = TRUE)

## ----loadsimdata3, include= FALSE---------------------------------------------
#save(sim_track, file = "./data/sim_example_3.rda")
load("./data/sim_example_3.rda")

## ----headsimcovar-------------------------------------------------------------
head(sim_track$barrier[,c("covar","kappa")])

## ----ctsimcovar, message = FALSE, results='hide'------------------------------
sim_cT <- buildCrossingTable(permdata = sim_track)

## ----processctsimcovar--------------------------------------------------------
sim_cT <- merge(sim_cT, sim_track$barrier[,c("barrier.id","covar")],
                by = "barrier.id", all.x = TRUE)

str(sim_cT)

## ----simcovarmodel, warning = FALSE-------------------------------------------
covar_model <- fitPermeability(~ covar, data = sim_cT)

## ----coefsimcovarmodel--------------------------------------------------------
covar_model$coef.table

## ----predictsimcovarmodel2----------------------------------------------------
predict_kappa <- predict(covar_model, se.fit = TRUE)
predict_kappa$covar <- covar_model$data$covar
predict_kappa <- predict_kappa |>
  dplyr::arrange(kappa.hat) # sort by increasing kappa values

head(predict_kappa)

## ----plotsimcovarmodel2, fig.align = "center", fig.width = 6, fig.height=4, eval = FALSE----
# plot(predict_kappa$covar, predict_kappa$kappa.hat, type = "l",
#      xlab = "Covariate Values", ylab = "Estimated Kappa Values", col ="orange", lwd = 2, ylim = c(0, 2))
# lines(predict_kappa$covar, predict_kappa$ci.high, col = "grey")
# lines(predict_kappa$covar, predict_kappa$ci.low, col = "grey")

## ----datasimcovar2, eval = FALSE----------------------------------------------
# n.segments <- 10
# barrier <- cbind(X = xmax * rnorm(n.segments+1, sd = xmax / 5e3),
#                      Y = ymax * seq(-1,1,length = n.segments+1),
#                      line_id = 1)
# covar_values <- sort(runif(n.segments, 1, 5))

## ----simcovar2, eval = FALSE--------------------------------------------------
# sim_tracks <- simulatePermeability(xmax = 50, ymax = 50,
#                               N.steps = 1000, theta.rho = 0.5,
#                               step.shape = 5, step.scale = 10,
#                               barrier = barrier,
#                               covar = covar_values,
#                               beta = c(-2,1),
#                               plot.track = FALSE)

## ----bordata------------------------------------------------------------------
data("boreal_caribou")
data("hwy1_west")
data("traffic_data")

## ----strbordata---------------------------------------------------------------
str(boreal_caribou)

## ----strbarrier---------------------------------------------------------------
str(hwy1_west)

## ----headbarrier--------------------------------------------------------------
hwy1_west_xy <- st_coordinates(hwy1_west)
head(hwy1_west_xy)

## ----checkuniquelines---------------------------------------------------------
unique(hwy1_west_xy[,"L1"])

## ----colnamesbarrier----------------------------------------------------------
colnames(hwy1_west_xy) <- c("X","Y","line_id")

## ----checkxyduplicates--------------------------------------------------------
hwy1_west_xy[566:573,]

## ----plotbarrier, fig.align = "center", fig.width = 7, fig.height=6, eval = FALSE----
# z <- hwy1_west_xy[,1] + 1i*hwy1_west_xy[,2]
# plot(z, type = "l")

## ----timesteps, eval = FALSE--------------------------------------------------
# time_steps <- boreal_caribou |>
#   group_by(ID) |>
#   arrange(Time) |>
#   mutate(dT = c(NA, as.numeric(round(difftime(Time[-1], Time[-length(Time)], units = "hours"))))) |>
#   ungroup() |>
#   data.frame() |>
#   select(dT)

## ----plottimesteps, fig.align = "center", fig.width = 6, fig.height=4, eval = FALSE----
# barplot(table(time_steps))

## ----prepbor1, fig.align = "center", fig.width = 6, fig.height=5, eval = FALSE----
# bor_format <- prepPermeability(move.df = subset(boreal_caribou, ID == "ID_3"),
#                          barrier = hwy1_west_xy,
#                          sample.rate = 8,
#                          plot.track = TRUE)

## ----barriercovariate, eval = FALSE-------------------------------------------
# z <- hwy1_west_xy[,1] + 1i*hwy1_west_xy[,2]
# 
# barrier_covar <- cbind(hwy1_west_xy, z) |>
#   data.frame() |>
#   group_by(line_id) |>
#   reframe(km = Mod(diff(z))) |>
#   select(km) |>
#   mutate(km = cumsum(km)/1000) |>
#   data.frame()

## ----prepborcovar, eval = FALSE-----------------------------------------------
# bor_format <- prepPermeability(move.df = boreal_caribou,
#                          move.covar.names = c("Season"),
#                          barrier.covar.values = barrier_covar,
#                          barrier = hwy1_west,
#                          sample.rate = 8)

## ----saveborformatcovar, include = FALSE--------------------------------------
#save(bor_format, file="./vignettes/data/bor_format.rda")
load(file="./data/bor_format.rda")

## ----classprepbor2------------------------------------------------------------
class(bor_format)

## ----strprepbor2track---------------------------------------------------------
str(bor_format$track)

## ----strprepbor2barrier-------------------------------------------------------
str(bor_format$barrier)

## -----------------------------------------------------------------------------
summary(bor_format)

## ----ctbor, eval = FALSE------------------------------------------------------
# bor_cT <- buildCrossingTable(permdata = bor_format)

## ----savectbor, include = FALSE-----------------------------------------------
#save(bor_cT, file = "./vignettes/data/bor_cT.rda")

load("./data/bor_cT.rda")

## ----parallelpackage, eval = FALSE--------------------------------------------
# require(future.apply)

## ----parallelct, eval = FALSE-------------------------------------------------
# bor_cT <- buildCrossingTable(permdata = bor_format,
#                              parallel = TRUE)

## ----strborct-----------------------------------------------------------------
str(bor_cT)

## ----mergepermdata, eval = FALSE----------------------------------------------
# track_covar <- bor_format$track[,c("ID", "Step.ID","Season")]
# 
# bor_cT <- merge(bor_cT, track_covar, by=c("ID", "Step.ID"),
#                     all.x=TRUE)

## ----mergepermdata2, eval = FALSE---------------------------------------------
# barrier_covar <- bor_format$barrier[,c("barrier.id", "km")]
# 
# bor_cT <- merge(bor_cT, barrier_covar, by=c("barrier.id"),
#                      all.x=TRUE)

## ----annotatetraffic, eval = FALSE--------------------------------------------
# bor_cT$Date <- as.Date(bor_cT$Time)
# 
# bor_cT <- merge(bor_cT, traffic_data, by = "Date", all.x = TRUE) |>
#   unique()

## ----annotatedoy, eval = FALSE------------------------------------------------
# bor_cT$doy <- lubridate::yday(bor_cT$Date)

## ----loadfinalct, include=FALSE-----------------------------------------------
#save(bor_cT, file = "./data/final_borct.rda")
load("./data/final_borct.rda")

## ----finalct------------------------------------------------------------------
str(bor_cT)

## ----plotcrossingsbasic, fig.align = "center", fig.width = 6, fig.height=6, eval = FALSE----
# plotCrossings(permdata = bor_format)

## ----mappingpackages, eval = TRUE---------------------------------------------
library(mapview)

## ----plotcrossings, fig.align = "center", fig.width = 6, fig.height=6, eval = TRUE----
plotCrossings(permdata = bor_format, 
              original.data = boreal_caribou,
              barrier.sf = hwy1_west,
              interactive.map = TRUE)

## ----plotcrossings2, fig.align = "center", fig.width = 8, fig.height=6, eval = FALSE----
# crossing_map <- plotCrossings(permdata = bor_format,
#                               original_data = boreal_caribou,
#                               barrier.sf = hwy1_west,
#                               static.map = TRUE,
#                               add.basemap = TRUE)
# 
# crossing_map

## ----loadmodels, echo = FALSE-------------------------------------------------
#save(null_model, traffic_model, season_model, additive_model, interactive_model, spline_doy_model, spline_km_model, doy_traffic_model, file = "./data/boreal_models.rda")

load("./data/boreal_models.rda")

## ----bornullmodel, eval = FALSE-----------------------------------------------
# null_model <- fitPermeability(data = bor_cT)
# null_model

## ----bornullmodel2, echo = FALSE----------------------------------------------
print(null_model)

## ----borcovarmodel, eval = FALSE----------------------------------------------
# traffic_model <- fitPermeability(~ traffic_volume, data = bor_cT)
# traffic_model

## ----borcovarmodel2, echo = FALSE---------------------------------------------
print(traffic_model)

## ----borcovarmodel3, eval = FALSE---------------------------------------------
# season_model <- fitPermeability(~ Season, data = bor_cT)
# season_model

## ----borcovarmodel4, echo = FALSE---------------------------------------------
print(season_model)

## ----boradditive, eval = FALSE------------------------------------------------
# additive_model <- fitPermeability(~ Season + traffic_volume, data = bor_cT)
# additive_model

## ----boradditive2, echo = FALSE-----------------------------------------------
print(additive_model)

## ----borinteractive, eval = FALSE---------------------------------------------
# interactive_model <- fitPermeability(~ Season*traffic_volume, data = bor_cT)
# interactive_model

## ----borinteractive2, echo = FALSE--------------------------------------------
print(interactive_model)

## ----splinedoy, eval = FALSE--------------------------------------------------
# spline_doy_model <- fitPermeability(~ s(doy, k = 4, bs = "cc"), data = bor_cT)
# spline_doy_model

## ----splinedoy2, echo = FALSE-------------------------------------------------
print(spline_doy_model)

## ----splinekm, eval = FALSE---------------------------------------------------
# spline_km_model <- fitPermeability(~ s(km, k = 6, bs = "cs"), data = bor_cT)
# spline_km_model

## ----splinekm2, echo = FALSE--------------------------------------------------
print(spline_km_model)

## ----borsplinecovar, eval = FALSE---------------------------------------------
# doy_traffic_model <- fitPermeability(~ traffic_volume + s(doy, k = 4, bs = "cc"),
#                                         data = bor_cT)
# doy_traffic_model

## ----borsplinecovar2, echo = FALSE--------------------------------------------
print(doy_traffic_model)

## ----modelcomparison, eval = FALSE--------------------------------------------
# bor_models <- list(Null = null_model,
#                    Traffic = traffic_model,
#                    Season = season_model,
#                    SeasonandTraffic = additive_model,
#                    SeasonxTraffic = interactive_model,
#                    SplineDoy = spline_doy_model,
#                    SplineKm = spline_km_model,
#                    SplineDoyandTraffic = doy_traffic_model)
# compare_models <- comparePermFits(models = bor_models)

## ----modelcomparison2, eval = FALSE-------------------------------------------
# bor_formulas <- list(Null = ~1,
#                      Traffic = ~traffic_volume,
#                      Season = ~Season,
#                      SeasonandTraffic = ~Season + traffic_volume,
#                      SeasonxTraffic = ~Season*traffic_volume,
#                      SplineDoy = ~ s(doy, k = 4, bs = "cc"),
#                      SplineKm = ~ s(km, k = 6, bs = "cs"),
#                      SplineDoyandTraffic = ~ traffic_volume + s(doy, k = 4, bs = "cc"))
# 
# compare_models <- comparePermFits(bor_formulas, data = bor_cT)
# compare_models

## ----loadmodelcomparison, echo = FALSE----------------------------------------
#save(compare_models, file = "./data/compare_models.rda")
load("./data/compare_models.rda")
compare_models

## ----predictborcovarmodel, eval = FALSE---------------------------------------
# predict_kappa <- data.frame(traffic_volume = traffic_model$data$traffic_volume)
# prediction <- predict(traffic_model, se.fit = TRUE)
# predict_kappa <- cbind(predict_kappa, prediction)

## ----predictborcovarmodel2, eval = FALSE--------------------------------------
# predict_kappa <- data.frame(traffic_volume =
#                          seq(0, 100, by = 1))
# prediction <- predict(traffic_model, new.data = predict_kappa, se.fit = TRUE)
# predict_kappa <- cbind(predict_kappa, prediction)

## ----plotborcovarmodel, fig.align = "center", fig.width = 6, fig.height=4, eval = FALSE----
# plot(predict_kappa$traffic_volume, predict_kappa$kappa.hat, type = "l",
#      xlab = "Daily Traffic Volume",
#      ylab = "Estimated Kappa Values",
#      col ="orange", lwd = 2, ylim = c(0, 0.15))
# lines(predict_kappa$traffic_volume, predict_kappa$ci.high, col = "grey")
# lines(predict_kappa$traffic_volume, predict_kappa$ci.low, col = "grey")

## ----bootstrap, eval = FALSE--------------------------------------------------
# prediction <- predict(traffic_model, se.fit = TRUE, bootstrap = TRUE)

## ----bootstrapparallel, eval = FALSE------------------------------------------
# prediction <- predict(traffic_model, se.fit = TRUE, bootstrap = TRUE, parallel = TRUE)

## ----dfpredictions, eval = FALSE----------------------------------------------
# hwy1_west_z <- hwy1_west_xy[,1] + 1i*hwy1_west_xy[,2]
# predictions <- cbind(hwy1_west_xy, hwy1_west_z) |>
#   data.frame() |>
#   group_by(line_id) |>
#   reframe(km = Mod(diff(hwy1_west_z))) |>
#   select(km) |>
#   mutate(km = cumsum(km)/1000) |>
#   mutate(barrier.id = 1:length(km)) |>
#   data.frame() |>
#   subset(barrier.id %in%  # make predictions for barrier segments in crossing table (crossings possible)
#            seq(min(bor_cT$barrier.id), max(bor_cT$barrier.id, by = 1)))

## ----getpredictions, eval = FALSE---------------------------------------------
# predictions$kappa.hat <- predict(spline_km_model, se.fit = FALSE,
#                                  new.data = predictions)[,1]
# str(predictions)

## ----loadkmpredictions, echo = FALSE------------------------------------------
#save(predictions, file = "./data/boreal_km_predictions.rda")
load("./data/boreal_km_predictions.rda")
str(predictions)

## ----mappredictions, eval = FALSE---------------------------------------------
# boreal_caribou_sf <- boreal_caribou |>
#   st_as_sf(coords = c("X","Y"), crs = st_crs(hwy1_west))
# 
# map <- mapPredictions(predictions, barrier.sf = hwy1_west,
#                       add.basemap = TRUE,
#                       tracks.sf = boreal_caribou_sf)
# map

## ----eval = FALSE-------------------------------------------------------------
# citation("permeability")

