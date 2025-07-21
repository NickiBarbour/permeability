require(permeability)
require(ggplot2)
require(ggspatial)
require(mapview)
require(basemaps)

data("example_id")
data("example_cT")
data("hwy1_west")

# fit model
hwy1_w_model <- fitPermeability(~ s(km, k = 4, bs = "cr"),
                                     data = example_cT)

# get predictions with barrier.id as column
hwy1_w_xy <- st_coordinates(hwy1_west)
colnames(hwy1_w_xy) <- c("X","Y","line_id")
hwy1_w_z <- hwy1_w_xy[,1] + 1i*hwy1_w_xy[,2]

predictions <- cbind(hwy1_w_xy, hwy1_w_z) |> 
  data.frame() |>
  group_by(line_id) |>
  reframe(km = Mod(diff(hwy1_w_z))) |>
  select(km) |>
  mutate(km = cumsum(km)/1000) |>
  mutate(barrier.id = 1:length(km)) |>
  data.frame() |>
  subset(barrier.id %in%  # make predictions for barrier segments in cT
           seq(min(example_cT$barrier.id), max(example_cT$barrier.id, by = 1)))

predictions$kappa.hat <- predict(hwy1_w_model, se.fit = FALSE, 
                                 new.data = predictions)[,1]

# static map (no basemap) ----------------------------------------------
map <- mapPredictions(predictions, barrier.sf = hwy1_west)

## add animal tracks
example_sf <- example_id |>
  st_as_sf(coords=c("X","Y"), crs = st_crs(hwy1_west))
map <- mapPredictions(predictions, barrier.sf = hwy1_west,
                      tracks.sf = example_sf)

# static map (with basemap) --------------------------------------------
## default - landscape from ESRI and basemaps package
map <- mapPredictions(predictions, barrier.sf = hwy1_west, 
                      add.basemap = TRUE)
map

## OSM streep map with ggspatial package
map <- mapPredictions(predictions, barrier.sf = hwy1_west, 
                      add.basemap = TRUE, maptype = "osm")
map

# interactive map ------------------------------------------------------
mapPredictions(predictions, barrier.sf = hwy1_west,
               interactive = TRUE)

