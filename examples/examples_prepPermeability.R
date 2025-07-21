require(permeability)

data("example_id")
data("hwy1_west")

# extract barrier coords (optional)
hwy1_west_xy <- st_coordinates(hwy1_west)
colnames(hwy1_west_xy) <- c("X","Y","line_id")

# no covariates
example_permdata <- prepPermeability(move.df = example_id,
                               barrier = hwy1_west_xy,
                               sample.rate = 8,
                               verbose = TRUE,
                               plot.track = TRUE)

str(example_permdata)

# Covariates example:

# e.g., barrier-specific covariate: km along hwy
z <- hwy1_west_xy[,1] + 1i*hwy1_west_xy[,2]

barrier_covar <- cbind(hwy1_west_xy, z) |> 
  data.frame() |>
  group_by(line_id) |>
  reframe(km = Mod(diff(z))) |>
  select(km) |>
  mutate(km = cumsum(km)/1000) |>
  data.frame()

# and movement track covariate: season
example_permdata <- prepPermeability(move.df = example_id,
                               move.covar.names = c("Season"), 
                               move.covar.values = "start",
                               barrier.covar.values = barrier_covar,
                               barrier = hwy1_west,
                               sample.rate = 8,
                               verbose = TRUE,
                               plot.track = TRUE)

# class permdata, list structure, one for track(s), one for barrier
class(example_permdata)

str(example_permdata)

# print: shows summary of data
print(example_permdata)

