#' @keywords internal
"_PACKAGE"

#' @importFrom stats AIC logLik predict as.formula model.matrix na.omit optim
#'   pchisq pnorm qnorm quantile rnorm runif rweibull sd terms
#' @importFrom graphics legend lines points segments title
#' @importFrom grDevices rgb
#' @importFrom sf st_as_sf st_as_sfc st_bbox st_buffer st_cast st_coordinates
#'   st_crop st_crs "st_crs<-" st_linestring st_transform
#' @importFrom dplyr group_by mutate reframe summarize ungroup
#' @importFrom MASS fitdistr
#' @importFrom mgcv s smoothCon
#' @import ggplot2
#' @importFrom ggspatial layer_spatial annotation_scale annotation_north_arrow
#'   annotation_map_tile shadow_spatial
NULL

utils::globalVariables(c(
  "D.time", "In.buffer", "L", "Crossed", "crossed", "null.crossed", "dAIC",
  "Feature.ID", "X", "Y", "Z", "Z1", "X1", "X2", "Y1", "Y2", "line_id", "ID",
  "Time", "row_no", "kappa.hat", "step_to_cross", "Z.end", "Z.start",
  "Step", "Phi"))