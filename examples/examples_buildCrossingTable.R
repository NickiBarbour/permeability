require(permeability)

data("example_permdata")
data("traffic_data")

# build crossing table
example_cT <- buildCrossingTable(permdata = example_permdata,
                                 verbose = TRUE)

str(example_cT)

## add covariate data from permdata object to crossing table
### covariate specific to movement step: can merge based on step.id
trackcovariate_data <- example_permdata$track[,c("ID", "Step.ID","Season")]
example_cT2 <- merge(example_cT, trackcovariate_data, by=c("ID", "Step.ID"), 
                    all.x=TRUE)

### covariate specific to barrier segment: can merge based on barrier.id
barriercovariate_data <- example_permdata$barrier[,c("barrier.id", "km")]
example_cT3 <- merge(example_cT2, barriercovariate_data, by=c("barrier.id"), 
                     all.x=TRUE)

str(example_cT3)

## covariate representing interaction of barrier and track:
## merge by date
example_cT3$Date <- as.Date(example_cT3$Time)
example_cT4 <- merge(example_cT3, traffic_data, by = "Date", all.x = TRUE) |>
  unique()

str(example_cT4)
