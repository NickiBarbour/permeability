#' Map Predictions  
#' 
#' @description Function to map kappa predictions as output from the `predict` 
#' function and a fitted permeability model from `fitPermeability`, to
#' segments of a barrier `sf` object. By default uses the `ggplot2` and 
#' `ggspatial` packages to create a static map (with optional basemap).
#' Users can also create an interactive map, using the `mapview` package.
#' @param predictions dataframe, at minimum, with predicted kappa values and 
#' a `barrier.id` column to identify distinct barrier segments. This column 
#' corresponds to the same `barrier.id` values in the `permdata` `barrier` 
#' object and crossing table object (as output from `buildCrossingTable`). If 
#' there are multiple linear features, these should be identified with an 
#' additional `Feature.ID` column; this column will be used to correspond with
#' a matching column in the provided `barrier.sf` object.
#' @param barrier.sf sf LINESTRING object(s) of the barrier. Should be the same
#' barrier object(s) used as input to the `prepPermeability` function. Can 
#' optionally have a `Feature.ID` column, which should have unique values 
#' that can be used to match with a similar column in the `predictions` object.
#' @param interactive if TRUE, will return an interactive map using the 
#' `mapview` package. Default FALSE.
#' @param tracks.sf optionally, an sf POINT object of animal track locations,
#' with a columns for `ID`. Will be used to map animal tracks behind the
#' barrier and should have a matching CRS to the `barrier.sf` object.
#' @param breaks Default NULL but if specified, should be a vector of values
#' to break up the predicted kappa values in the plot legend. If NULL, will
#' use a log scale to break up data and plot it. Applicable for static maps
#' only.
#' @param add.basemap if TRUE (default FALSE), will add a basemap for a
#' static map, using an ESRI landscape map from the `basemaps` package. 
#' Alternatively you can provide your own basemap (as a raster) with the 
#' `basemap` argument.
#' @param basemap default NULL but you can prove your own basemap for plotting
#' here. Should be a raster with similar spatial extent and CRS as your
#' `barrier.sf` object. Applicable for static maps only. 
#' @param box default NULL but you can provide a different spatial extent for
#' the map using this box parameter, which should be an sf bbox object with a
#' CRS matching the barrier sf object. If NULL, will use the spatial extent of 
#' the provided barrier object to create map. Applicable to static maps only. 
#' @param map.service for static maps with an automatic basemap (provided 
#' through the "basemaps" package), the map service to use - default is "esri"
#' but you can also use other options available through the "basemap_raster"
#' function from basemaps (e.g., "osm").
#' @param map.type for static maps with an automatic basemap (provided 
#' through the "basemaps" package), the map type to use - default is 
#' "world_imagery" from ESRI but you can provide any other options 
#' available through the "get_maptypes" function from the basemaps package 
#' (just provide your map service type).
#' @param map.res for static maps with an automatic basemap (provided 
#' through the "basemaps" package), the resolution of the map (0-1). Default 
#' NULL (will use the basemaps default).
#' @return An interactive or static map of predicted kappa values for barrier 
#' segments.
#' @example examples/examples_mapPredictions.R
#' @export

mapPredictions <- function(predictions, barrier.sf, interactive = FALSE,
                           tracks.sf = NULL, breaks = NULL,
                           add.basemap = FALSE, basemap = NULL,
                           box = NULL, map.service = "esri", 
                           map.type = "world_imagery", map.res = NULL){
  if(! "barrier.id" %in% colnames(predictions)) 
    stop("Predictions dataframe does not have a barrier.id column.")
  # check for multiple linear features
  if("Feature.ID" %in% colnames(barrier.sf)){
    barrier_segment <- plyr::ddply(predictions, "Feature.ID", 
                                   getAnnotatedBarrierSegments,
                                   barrier_sf = barrier.sf) |>
      sf::st_as_sf()
  }else{
    barrier_segment <- getAnnotatedBarrierSegments(predictions,
                                                   barrier_sf = barrier.sf)
  }
  suppressWarnings(
    # map results
    if(interactive){
      barrier_segment$kappa.hat <- signif(barrier_segment$kappa.hat,
                                          digits = 1)
      if(is.null(tracks.sf)){
        map <- mapview::mapview(barrier_segment, 
                                zcol = "kappa.hat", lwd = 3,
                                layer.name = "kappa")
        return(map)
      }else{
        tracks <- suppressWarnings(tracks.sf |> 
          dplyr::group_by(ID) |>
          dplyr::summarize(do_union=FALSE) |>
          sf::st_cast("LINESTRING"))
        map <- mapview::mapview(tracks, color = "grey40", legend = FALSE) +
          mapview::mapview(barrier_segment, 
                           zcol = "kappa.hat", lwd = 3, layer.name = "kappa")
        return(map)
      }
    }else{
      barrier_segment <- barrier_segment |>
        st_transform(4326)
      barrier_segment$kappa.hat <- signif(barrier_segment$kappa.hat,
                                          digits = 1)
      if(is.null(breaks)){
        if(0 %in% unique(barrier_segment$kappa.hat)){
          barrier_segment[which(barrier_segment$kappa.hat==0),]$kappa.hat <- 1e-4
        }
        breaks <- pretty(log10(barrier_segment$kappa.hat))
        breaks_labels <- signif(10^breaks, 1)
        log_kappa <- TRUE
      }else{
        breaks <- breaks
        breaks_labels <- breaks
        log_kappa <- FALSE
      }
      if(is.null(tracks.sf)){
        if(add.basemap){
          if(is.null(basemap)){
            if(is.null(box)){
              box <- barrier.sf |>
                st_bbox()
              box <- barrier.sf |>
                st_buffer(abs(as.numeric(box[1])-as.numeric(box[2]))*.001) |>
                st_transform(4326) |>
                st_bbox() 
            }
            bm <- suppressWarnings(basemaps::basemap_raster(box, 
                                           map_service = map.service, 
                                           map_type = map.type,
                                           map_res = map.res))
            
            box <- box |>
              st_as_sfc() |> 
              st_transform(crs = st_crs(bm)) |>
              st_bbox()
            barrier_segment2 <- barrier_segment |>
              st_transform(st_crs(bm))
            if(log_kappa){
              barrier_segment2$kappa.hat <- log10(barrier_segment2$kappa.hat)
            }
            
            map <- ggplot() +
              layer_spatial(data = bm) +
              guides(fill = "none") +
              geom_sf(data = barrier_segment2, aes(color = kappa.hat), 
                      linewidth = 2) +
              scale_color_viridis_c(option = "B", breaks = breaks,
                                    labels = breaks_labels,
                                    name = expression(hat(kappa))) +
              annotation_scale() +
              annotation_north_arrow(height = unit(0.5,"cm"), 
                                     width = unit(0.5,"cm"), 
                                     pad_y = unit(1,"cm")) +
              xlab("Longitude") + ylab("Latitude")+
              coord_sf(expand=FALSE)+
              theme_classic() + theme(text = element_text(size = 11))
          }else if(!is.null(basemap)){
            bm <- basemap
            if(is.null(box)){
              box <- barrier.sf |>
                st_bbox()
              box <-  suppressWarnings(barrier.sf |>
                st_buffer(abs(as.numeric(box[1])-as.numeric(box[2]))*.001) |>
                st_transform(st_crs(bm)) |>
                st_bbox())
            }else{
              box <- box |>
                st_as_sfc() |> 
                st_transform(crs = st_crs(bm)) |>
                st_bbox()
            }
            barrier_segment2 <- barrier_segment |>
              st_transform(st_crs(bm))
            if(log_kappa){
              barrier_segment2$kappa.hat <- log10(barrier_segment2$kappa.hat)
            }
            
            map <- ggplot() +
              layer_spatial(data = bm) +
              guides(fill = "none") +
              geom_sf(data = barrier_segment2, aes(color = kappa.hat), 
                      linewidth = 2) +
              scale_color_viridis_c(option = "B", breaks = breaks,
                                    labels = breaks_labels,
                                    name = expression(hat(kappa))) +
              annotation_scale() +
              annotation_north_arrow(height = unit(0.5,"cm"), 
                                     width = unit(0.5,"cm"), 
                                     pad_y = unit(1,"cm")) +
              xlab("Longitude") + ylab("Latitude")+
              coord_sf(expand=FALSE)+
              theme_classic() + theme(text = element_text(size = 11))
            
          }
          return(map)
        }else{
          if(is.null(box)){
            box <- barrier.sf |>
              st_transform(4326) |>
              st_bbox()
          }else{
            box <- box |>
              st_as_sfc() |> 
              st_transform(crs = 4326) |>
              st_bbox()
          }
          if(log_kappa){
            barrier_segment$kappa.hat <- log10(barrier_segment$kappa.hat)
          }
          map <- ggplot() +
            geom_sf(data = barrier_segment, aes(color = kappa.hat), 
                    linewidth = 2) +
            scale_color_viridis_c(option = "B", breaks = breaks,
                                  labels = breaks_labels,
                                  name = expression(hat(kappa))) +
            annotation_scale() +
            annotation_north_arrow(height = unit(0.5,"cm"), 
                                   width = unit(0.5,"cm"), 
                                   pad_y = unit(1,"cm")) +
            xlab("Longitude") + ylab("Latitude")+
            shadow_spatial(box) +
            theme_classic() + theme(text = element_text(size = 11))
          return(map)
        }
      }else{
        if(add.basemap){
          if(is.null(basemap)){
            if(is.null(box)){
              box <- barrier.sf |>
                st_bbox()
              box <- barrier.sf |>
                st_buffer(abs(as.numeric(box[1])-as.numeric(box[2]))*.001) |>
                st_transform(4326) |>
                st_bbox()
            }
            bm <- suppressWarnings(basemaps::basemap_raster(box, 
                                                            map_service = 
                                                              map.service, 
                                                            map_type = map.type,
                                                            map_res = map.res))
            box <- box |>
              st_as_sfc() |> 
              st_transform(crs = st_crs(bm)) |>
              st_bbox()
            barrier_segment2 <- barrier_segment |>
              st_transform(st_crs(bm))
            if(log_kappa){
              barrier_segment2$kappa.hat <- log10(barrier_segment2$kappa.hat)
            }
            tracks <- suppressWarnings(tracks.sf |> 
                                         dplyr::group_by(ID) |>
                                         dplyr::summarize(do_union=FALSE) |>
                                         sf::st_cast("LINESTRING") |> 
                                         st_transform(st_crs(bm)) |>
                                         sf::st_crop(box))
            map <- ggplot() +
              layer_spatial(data = bm) +
              guides(fill = "none") +
              geom_sf(data = tracks, color = "grey40", linewidth = 0.7) +
              geom_sf(data = barrier_segment2, aes(color = kappa.hat), 
                      linewidth = 2) +
              scale_color_viridis_c(option = "B", breaks = breaks,
                                    labels = breaks_labels,
                                    name = expression(hat(kappa))) +
              annotation_scale() +
              annotation_north_arrow(height = unit(0.5,"cm"), 
                                     width = unit(0.5,"cm"), 
                                     pad_y = unit(1,"cm")) +
              xlab("Longitude") + ylab("Latitude")+
              coord_sf(expand=FALSE)+
              theme_classic() + theme(text = element_text(size = 11))
          }else if(!is.null(basemap)){
            bm <- basemap
            if(is.null(box)){
              box <- barrier.sf |>
                st_bbox()
              box <- barrier.sf |>
                st_buffer(abs(as.numeric(box[1])-as.numeric(box[2]))*.001) |>
                st_transform(st_crs(bm)) |>
                st_bbox() 
            }
            barrier_segment2 <- barrier_segment |>
              st_transform(st_crs(bm))
            if(log_kappa){
              barrier_segment2$kappa.hat <- log10(barrier_segment2$kappa.hat)
            }
            tracks <- suppressWarnings(tracks.sf |> 
                                         dplyr::group_by(ID) |>
                                         dplyr::summarize(do_union=FALSE) |>
                                         sf::st_cast("LINESTRING") |> 
                                         st_transform(st_crs(bm)) |>
                                         sf::st_crop(box))
            map <- ggplot() +
              layer_spatial(data = bm) +
              guides(fill = "none") +
              geom_sf(data = tracks, color = "grey40", linewidth = 0.7) +
              geom_sf(data = barrier_segment2, aes(color = kappa.hat), 
                      linewidth = 2) +
              scale_color_viridis_c(option = "B", breaks = breaks,
                                    labels = breaks_labels,
                                    name = expression(hat(kappa))) +
              annotation_scale() +
              annotation_north_arrow(height = unit(0.5,"cm"), 
                                     width = unit(0.5,"cm"), 
                                     pad_y = unit(1,"cm")) +
              xlab("Longitude") + ylab("Latitude")+
              coord_sf(expand=FALSE)+
              theme_classic() + theme(text = element_text(size = 11))
          }
          return(map)
        }else{
          if(is.null(box)){
            box <- barrier.sf |>
              st_transform(4326) |>
              st_bbox()
          }else{
            box <- box |>
              st_as_sfc() |> 
              st_transform(crs = 4326) |>
              st_bbox()
          }
          if(log_kappa){
            barrier_segment$kappa.hat <- log10(barrier_segment$kappa.hat)
          }
          tracks <- suppressWarnings(tracks.sf |> 
                                       dplyr::group_by(ID) |>
                                       dplyr::summarize(do_union=FALSE) |>
                                       sf::st_cast("LINESTRING") |> 
                                       st_transform(4326) |>
                                       sf::st_crop(box))
          map <- ggplot() +
            geom_sf(data = tracks, color = "grey40", linewidth = 0.7) +
            geom_sf(data = barrier_segment, aes(color = kappa.hat), 
                    linewidth = 2) +
            scale_color_viridis_c(option = "B", breaks = breaks,
                                  labels = breaks_labels,
                                  name = expression(hat(kappa))) +
            annotation_scale() +
            annotation_north_arrow(height = unit(0.5,"cm"), 
                                   width = unit(0.5,"cm"), 
                                   pad_y = unit(1,"cm")) +
            xlab("Longitude") + ylab("Latitude")+
            shadow_spatial(box) +
            theme_classic() + theme(text = element_text(size = 11))
          return(map)
        }
      }
    }
  )
}

getAnnotatedBarrierSegments <- function(feature_predictions,
                                        barrier_sf){
  if("Feature.ID" %in% colnames(feature_predictions)){
    if(! "Feature.ID" %in% colnames(feature_predictions)) 
      stop("Feature.ID column in predictions dataframe but not barrier sf object.")
    id <- feature_predictions$Feature.ID[1]
    barrier_sf <- subset(barrier_sf, Feature.ID==id)
  }
  barrier_xy <- st_coordinates(barrier_sf)
  colnames(barrier_xy) <- c("X","Y","line_id")
  # for each line_id: convert XY to Z, then convert to stepwise dataframe
  barrier_step_df <- barrier_xy |>
    data.frame() |>
    dplyr::mutate(Z = X + 1i*Y) |>
    dplyr::group_by(line_id) |>
    dplyr::reframe(Z1 = Z[-length(Z)],
                   Z2 = Z[-1]) |>
    dplyr::mutate(barrier.id = 1:length(Z1)) |>
    data.frame()
  # create start pts and end pts
  start_points <- cbind(X1 = Re(barrier_step_df$Z1), 
                        Y1 = Im(barrier_step_df$Z1)) |> data.frame()
  end_points <- cbind(X2 = Re(barrier_step_df$Z2), 
                      Y2 = Im(barrier_step_df$Z2)) |> data.frame()
  
  points <- cbind(start_points, end_points) |>
    dplyr::mutate(line_id = barrier_step_df$line_id,
                  barrier.id = barrier_step_df$barrier.id) |>
    unique()
  # convert points to line
  steps <- points |> 
    dplyr::rowwise() |> 
    dplyr::mutate(
      geometry = list(
        st_linestring(
          matrix(c(X1, Y1, X2, Y2), ncol = 2, byrow = TRUE)))) |>
    ungroup() |> 
    data.frame() |>
    dplyr::select(-c(X1, Y1, X2, Y2)) |>
    sf::st_as_sf()
  st_crs(steps) <- st_crs(barrier_sf)
  barrier_segment <- merge(steps, feature_predictions, by = "barrier.id",
                           all.x = TRUE) |> unique()
  # change any NA's for kappa.hat to 0 
  ## (segments with no actual/potential crossings)
  if(nrow(subset(barrier_segment, is.na(kappa.hat)))>0){
    barrier_segment[which(is.na(
      barrier_segment$kappa.hat)),]$kappa.hat <- 0
  }
  if("Feature.ID" %in% colnames(feature_predictions)){
    barrier_segment$Feature.ID <- id
  }
  return(barrier_segment)
}
