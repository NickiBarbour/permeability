library(sf);library(mapview);library(dplyr)

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

mapview(boreal_tracks, zcol = "ID") + 
  mapview(hwy1_west, color = "red")

