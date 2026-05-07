
#' Highway 1 West spatial data object
#'
#' @usage
#' data(hwy1_west)
#' @name hwy1_west
#' @format Linear simple feature (LINESTRING)
#' @example examples/datadocumentation_examples_borealcaribou.R
#' @source Extracted and processed from the Road Network File - 2005 - 
#' Northwest Territories 
#' (https://open.canada.ca/data/en/dataset/97adfa76-7eee-477b-b3d3-e430be7bd812)
#' @keywords data
"hwy1_west"

#' Boreal caribou (*Rangifer tarandus*) movement data  
#' 
#' 
#' N= 10 anonymized caribou, 2015-2017, with an 8 hr median sampling rate and 
#' annotated by relevant ecological seasons. 
#' @usage 
#' data(boreal_caribou)
#' @name boreal_caribou
#' @format Data frame with columns for ID, Time (POSIX UTC),
#' X,Y spatial coordinates (ESPG code 32611), and season.
#' @example examples/datadocumentation_examples_borealcaribou.R
#' @source GNWT
#' @keywords data
"boreal_caribou"

#' Traffic data for HWY 1 West
#' 
#' Daily traffic volume data for Highway 1 West, from a GNWT traffic counter.
#'
#' @usage 
#' data(traffic_data)
#' @name traffic_data
#' @format Data frame with columns for Date (POSIX UTC),
#' and traffic volume (number of vehicles),
#' @example examples/datadocumentation_examples_borealcaribou.R
#' @source GNWT
#' @keywords data
"traffic_data"

#' Example caribou individual 
#' 
#' Example caribou individual, for use in running package function examples.
#'
#' @usage 
#' data(example_id)
#' @name example_id
#' @example examples/datadocumentation_examples_borealcaribou.R
#' @keywords data
"example_id"

#' Example permdata 
#' 
#' Example permdata object, as the result of running the examples in the 
#' prepPermeability  function.
#'
#' @usage 
#' data(example_permdata)
#' @name example_permdata
#' @example examples/datadocumentation_examples_borealcaribou.R
#' @keywords data
"example_permdata"

#' Example crossing table
#' 
#' Example crossing table, as the result of running the examples in the 
#' BuildCrossingTable function.
#'
#' @usage 
#' data(example_cT)
#' @name example_cT
#' @example examples/datadocumentation_examples_borealcaribou.R
#' @keywords data
"example_cT"
