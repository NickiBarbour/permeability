#' Plot Crossings 
#' 
#' @description Function to plot crossings as output from either
#'  \link{prepPermeability} or \link{simulatePermeability}. Defaults to basic 
#'  base R plots but for "real" data (e.g., data not simulated with
#'  `simulatePermeability`) users can load additional packages for interactive 
#'  spatial maps (`mapview` package) or static spatial maps (`ggplot2` package),
#'  optionally with a basemap (`ggspatial` package).
#' @param permdata output from \link{prepPermeability} with `track` and 
#' `barrier` objects, as stepwise dataframes.
#' @param original.data needed if making an interactive or static map with the
#' `mapview` or `ggplot2`/`ggspatial` packages. Should be a dataframe of 
#' original data, containing columns for:
#' `ID`, `Time` (POSIXct format, date and time), and spatial coordinates 
#' (`X`/`Y` - units in meters) - coordinate system should match barrier `sf` 
#' object. Default NULL. 
#' @param barrier.sf needed if making an interactive or static map with the
#' `mapview` or `ggplot2`/`ggspatial` packages. Should be an sf LINESTRING
#' object and coordinate reference system should match `original.data`. Default
#' NULL.
#' @param interactive.map if TRUE, will return interactive map(s) using the 
#' `mapview` package. Default NULL.
#' @param static.map if TRUE, will return static map(s) using the `ggplot2` and
#' optionally the `ggspatial` package (if `add.basemap = TRUE`). Default FALSE.
#' @param individual.plots if TRUE, will return individual plots, for those 
#' individuals that crossed the barrier - default FALSE (plots all tracks).
#' @param add.basemap if TRUE (default FALSE), will add a basemap for a
#' static map, using the `ggspatial` package. If your data has a large 
#' spatial extent, it might take a couple minutes to compile the basemap.
#' @return Plot(s) - if multiple individuals, will be returned as a list of 
#' plots.
#' @example examples/examples_plotCrossings.R
#' @export

plotCrossings <- function(permdata, original.data = NULL, 
                          barrier.sf = NULL, interactive.map = FALSE, 
                          static.map = FALSE,
                          individual.plots = FALSE,
                          add.basemap = FALSE){
  if(sum(permdata$track$Crossed, na.rm=TRUE)==0){
    stop("No crossings in permdata object. Nothing to plot!")
  }
  if(interactive.map == FALSE & static.map == FALSE){
    # track and barrier segments
    ids_in_buff <- unique(subset(permdata$track, In.buffer)$ID)
    track <- subset(permdata$track, ID %in% ids_in_buff)
    track <- track[,1:2] |> as.matrix()
    barrier <- permdata$barrier[,1:2] |> as.matrix()
    
    if(individual.plots == FALSE){
      # extract start/end locs for steps that cross
      z.cross1 <- subset(permdata$track, Crossed==TRUE)$Z.start
      z.cross2 <- subset(permdata$track, Crossed==TRUE)$Z.end
      
      plot(track, col = "darkgrey", xlab = "X", ylab = "Y", cex = 0.5, asp = 1)
      apply(barrier, 1, function(x) lines(x, lwd = 2, col = "black"))
      apply(track, 1, function(x) lines(x, lwd = 1, col = "darkgrey"))
      segments(Re(z.cross1), Im(z.cross1), Re(z.cross2), Im(z.cross2), 
               lwd = 2.5, col = "orange")
      legend("bottom", legend=c("steps that cross"), 
             fill=c("orange"), cex = 0.7,inset=c(1,1), xpd=TRUE, horiz = TRUE)
    }else{
      if(length(unique(subset(permdata$track,Crossed)$ID))<2) stop("You set individual.plots = TRUE, but there is only one individual that crossed. Please keep the default settings (individual.plots = FALSE).")
      # find ids that had tracks that crossed barrier
      ids_to_cross <- unique(subset(permdata$track,Crossed)$ID)
      message(length(ids_to_cross), " tracks crossed the barrier and will be plotted individually.")
      
      # select ids in original data that crossed
      track <- permdata$track
      barrier <- permdata$barrier
      
      track_crossed <- track |>
        subset(ID %in% ids_to_cross)
      
      # split into list of ids
      crossed_list <- split(track_crossed, track_crossed$ID)
      
      lapply(crossed_list,
             plotCrossings_individual,
             barrier = barrier)
      
      
    }
    
  }else{
    # subset only for steps that cross
    track <- subset(permdata$track, Crossed)
    
    track$ID <- as.character(track$ID)
    original.data$ID <- as.character(original.data$ID)
    
    ids_to_cross <- unique(track$ID)
    original.data <- subset(original.data, ID %in% ids_to_cross)
    
    steps_to_cross <- plyr::ddply(track, "ID",
                                  getSteps,
                                  barrier = barrier.sf) |>
      sf::st_as_sf()
    
    if(length(steps_to_cross)==0){
      stop("No steps crossed barrier. Nothing to map!")
    }
    
    if( individual.plots == FALSE){
      
      # subset for tracks in buffer
      ids_in_buff <- unique(subset(permdata$track, In.buffer)$ID)
      
      original.data <- subset(original.data, ID %in% ids_in_buff)
      
      if(interactive.map){
        
          original_tracks <- suppressWarnings(original.data |>
            st_as_sf(coords=c("X","Y"), crs = st_crs(barrier.sf)) |> 
            dplyr::group_by(ID) |>
            summarize(do_union=FALSE) |> 
            st_cast("LINESTRING"))
          
          map <- suppressWarnings(mapview(original_tracks, color = "darkgrey", legend = FALSE) +
            mapview(steps_to_cross, color = "orange") +
            mapview(barrier.sf, color = "black"))
          
          return(map)
        
      }else if(static.map){
        
        box <- suppressWarnings(steps_to_cross |>
          st_transform(4326) |>
          st_bbox())
        
        original_tracks_crop <- suppressWarnings(original.data |>
          group_by(ID) |>
          dplyr::mutate(row_no = length(Time)) |>
          subset(row_no >= 2) |> # drop any tracks w/ < 2 obs
          st_as_sf(coords=c("X","Y"), crs = st_crs(barrier.sf)) |> 
          st_transform(4326) |>
          group_by(ID) |>
          summarize(do_union=FALSE) |> 
          st_cast("LINESTRING") |>
          st_crop(box))
        
        barrier_crop <- suppressWarnings(barrier.sf |>
          st_transform(4326) |>
          st_crop(box))
        
        steps_to_cross2 <- steps_to_cross |>
          st_transform(4326)
        
        if(add.basemap){
          
          map <- suppressMessages(ggplot() +
            annotation_map_tile(type = 'osm', zoomin = -1) +
            geom_sf(data = original_tracks_crop, color = "darkgrey", 
                    linewidth = 0.7) +
            geom_sf(data = barrier_crop, color = "black") +
            geom_sf(data = steps_to_cross2, color = "orange", linewidth = 1) +
            annotation_scale() +
            annotation_north_arrow(height = unit(0.5,"cm"), 
                                   width = unit(0.5,"cm"), 
                                   pad_y = unit(1,"cm")) +
            xlab("Longitude") + ylab("Latitude")+
            shadow_spatial(box) +
            theme_classic() + theme(text = element_text(size = 11)))
          
          return(map)
        }else{
          
          map <- suppressWarnings(ggplot() +
            geom_sf(data = original_tracks_crop, color = "grey", 
                    linewidth= 0.7) +
            geom_sf(data = barrier_crop, color = "black") +
            geom_sf(data = steps_to_cross2, color = "orange", linewidth = 1) +
            theme_classic() +
            #shadow_spatial(box) +
            theme(text = element_text(size = 11)))
          
          return(map)
        }
      }
    }else{
      # subset for tracks in buffer
      ids_in_buff <- unique(subset(permdata$track, In.buffer)$ID)
      
      original.data <- subset(original.data, ID %in% ids_in_buff)
      
      # subset only for steps that cross
      track <- subset(permdata$track, Crossed)
      
      track$ID <- as.character(track$ID)
      original.data$ID <- as.character(original.data$ID)
      
      ids_to_cross <- unique(track$ID)
      og_data <- subset(original.data, ID %in% ids_to_cross)
      
      # get steps that cross
      steps_to_cross <- plyr::ddply(track, "ID",
                                    getSteps,
                                    barrier = barrier.sf) |>
        sf::st_as_sf()
      
      # split into list of ids
      og_data_list <- split(og_data, og_data$ID)
      
      plot_list <- lapply(og_data_list,
                          plotCrossings_individual,
                          barrier = barrier.sf,
                          steps_to_cross = steps_to_cross,
                          interactive.map = interactive.map,
                          static.map = static.map,
                          add.basemap = add.basemap)
      
      return(plot_list)
    }
  }
}

# @rdname plotCrossings
# @export
#' 
plotCrossings_individual <- function(df_ID, barrier,
                                     interactive.map = FALSE, 
                                     static.map = FALSE,
                                     add.basemap = FALSE){
  
  if(interactive.map==FALSE & static.map==FALSE){
    # track and barrier segments
    track <- df_ID[,1:2] |> as.matrix()
    barrier <- barrier[,1:2] |> as.matrix()
    
    z.cross1 <- subset(df_ID, Crossed==TRUE)$Z.start
    z.cross2 <- subset(df_ID, Crossed==TRUE)$Z.end
    
    plot(track, col = "darkgrey", xlab = "X", ylab = "Y", cex = 0.5,
         main = paste("ID", df_ID$ID[1]), asp = 1)
    apply(barrier, 1, function(x) lines(x, lwd = 2, col = "black"))
    apply(track, 1, function(x) lines(x, lwd = 1, col = "darkgrey"))
    segments(Re(z.cross1), Im(z.cross1), Re(z.cross2), Im(z.cross2), lwd = 2.5, 
             col = "orange")
    legend("bottomleft", legend=c("steps that cross"), 
           fill=c("orange"), cex = 0.7,inset=c(0.7,1), xpd=TRUE, horiz = TRUE)
  }else{
    if(interactive.map){
        original_track <- suppressWarnings(df_ID |>
          st_as_sf(coords=c("X","Y"), crs = st_crs(barrier)) |> 
          summarize(do_union=FALSE) |> 
          st_cast("LINESTRING"))
        
        map <- suppressWarnings(mapview(original_track, color = "darkgrey", legend = FALSE) +
          mapview(step_to_cross, color = "orange") +
          mapview(barrier, color = "black"))
        
        return(map)
      
    }else if(static.map){
        box <- suppressWarnings(step_to_cross |>
          st_transform(4326) |>
          st_bbox())
        
        original_track_crop <- suppressWarnings(df_ID |>
          st_as_sf(coords=c("X","Y"), crs = st_crs(barrier)) |> 
          st_transform(4326) |>
          summarize(do_union=FALSE) |> 
          st_cast("LINESTRING") |>
          st_crop(box))
        
        barrier_crop <- suppressWarnings(barrier |>
          st_transform(4326) |>
          st_crop(box))
        
        step_to_cross2 <- step_to_cross |>
          st_transform(4326)
        
        if(add.basemap){
          
          map <- suppressMessages(ggplot() +
            annotation_map_tile(type = 'osm', zoomin = -1) +
            geom_sf(data = original_track_crop, color = "darkgrey", 
                    linewidth = 0.7) +
            geom_sf(data = barrier_crop, color = "black") +
            geom_sf(data = step_to_cross2, color = "orange", linewidth = 1) +
            annotation_scale() +
            annotation_north_arrow(height = unit(0.5,"cm"), 
                                   width = unit(0.5,"cm"), 
                                   pad_y = unit(1,"cm")) +
            xlab("Longitude") + ylab("Latitude")+
            shadow_spatial(box) +
            theme_classic() + theme(text = element_text(size = 11)) +
            ggtitle(paste("ID", unique(df_ID$ID))))
          
          return(map)
        }else{
          
          map <- suppressWarnings(ggplot() +
            geom_sf(data = original_track_crop, color = "grey", linewidth = 0.7) +
            geom_sf(data = barrier_crop, color = "black") +
            geom_sf(data = step_to_cross2, color = "orange", linewidth = 1) +
            theme_classic() +
            theme(text = element_text(size = 11)) +
            #shadow_spatial(box) +
            ggtitle(paste("ID", unique(df_ID$ID))))
          
          return(map)
        }
    }
  }
}



getSteps <- function(track.ID, barrier){
  # function to convert sets of point to line
  points_to_line <- function(row){
    pt1_sf <- row |> 
      dplyr::select(X1, Y1) |>
      sf::st_as_sf(coords = c("X1", "Y1"), crs = st_crs(barrier))
    pt2_sf <- row |> 
      dplyr::select(X2, Y2) |>
      sf::st_as_sf(coords = c("X2", "Y2"), crs = st_crs(barrier))
    
    line <- rbind(pt1_sf, pt2_sf) |>
      sf::st_union() |>
      sf::st_cast("LINESTRING") |>
      data.frame() |>
      sf::st_as_sf()
    
    return(line)
  }
  
  start_points <- cbind(X1 = Re(track.ID[,1]), Y1 = Im(track.ID[,1])) |>
    data.frame()
  end_points <- cbind( X2 = Re(track.ID[,2]), Y2 = Im(track.ID[,2])) |>
    data.frame()
  
  points <- cbind(start_points, end_points) |>
    unique()
  
  step_list <- c()
  for(i in c(1:nrow(points))){
    step_list[[i]] <- points_to_line(points[i,])
  }
  
  steps <- do.call("rbind", step_list)
  
  steps$ID <- track.ID$ID[1]
  
  steps$Step.ID <- track.ID$Step.ID
  
  return(steps)
}
