library(permeability)
library(ggplot2)
library(ggspatial)

data("example_id")
data("example_permdata")
data("hwy1_west")

# base R version -----------------------------

## simulated data
sim_tracks <- simulatePermeability(xmax = 30, ymax = 30, 
                                   N.steps = 50, theta.rho = 0.5, 
                                   step.shape = 5, step.scale = 10, 
                                   kappa = 0.2,
                                   n.tracks = 10,
                                   n.segments = 20,
                                   plot.track = FALSE) 

plotCrossings(permdata = sim_tracks)

## "real" data
plotCrossings(permdata = example_permdata)

# static map (ggspatial and ggplot packages) --------------------------------

## just ggplot
crossing_map <- plotCrossings(permdata = example_permdata, 
                              original.data = example_id,
                              barrier.sf = hwy1_west,
                              static.map = TRUE)

crossing_map

\donttest{
# interactive map (mapview package)
if(requireNamespace("mapview", quietly = TRUE)){
  crossing_map <- plotCrossings(permdata = example_permdata, 
                                original.data = example_id,
                                barrier.sf = hwy1_west,
                                interactive.map = TRUE)
  
  crossing_map
}

## with basemap (requires internet)
crossing_map <- plotCrossings(permdata = example_permdata, 
                              original.data = example_id,
                              barrier.sf = hwy1_west,
                              static.map = TRUE,
                              add.basemap = TRUE)

crossing_map
}
