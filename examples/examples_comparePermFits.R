library(permeability)

data("example_cT")

# day of year (base R)
example_cT$doy <- as.numeric(format(example_cT$Date, "%j"))

model_formulae <- list(null = ~ 1, 
                       Traffic = ~ traffic_volume, 
                       Season = ~ Season, 
                       sDOY = ~ s(doy, k = 5, bs = "cc"),
                       Traffic_Season = ~ traffic_volume + Season)

comparePermFits(model_formulae, data = example_cT)
