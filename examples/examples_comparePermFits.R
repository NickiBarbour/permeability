require(permeability)
require(lubridate)

data("example_cT")

example_cT$doy <- lubridate::yday(example_cT$Date)

model_formulae <- list(null = ~ 1, 
                       Traffic = ~ traffic_volume, 
                       Season = ~ Season, 
                       sDOY = ~ s(doy, k = 5, bs = "cc"),
                       Traffic_Season = ~ traffic_volume + Season)

comparePermFits(model_formulae, data = example_cT)
