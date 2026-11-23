library(sf)
library(dplyr)

# load caribou data
data(boreal_caribou)

# load highway data
data(hwy1_west)

# convert to sf structure
boreal_caribou_sf <- st_as_sf(boreal_caribou, 
                              coords = c("X","Y"), crs = 32611)

# map track highway data
boreal_tracks <- boreal_caribou_sf |> 
  group_by(ID) |>
  summarize(do_union=FALSE) |> 
  st_cast("LINESTRING")

plot(st_geometry(boreal_tracks), col = "grey40")
plot(st_geometry(hwy1_west), col = "red", lwd = 2, add = TRUE)

# interactive map
if(interactive() && requireNamespace("mapview", quietly = TRUE)){
  mapview::mapview(boreal_tracks, zcol = "ID") + 
    mapview::mapview(hwy1_west, color = "red")
}
