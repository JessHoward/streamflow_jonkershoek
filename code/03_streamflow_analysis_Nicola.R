### Streamflow data analysis

# Load packages
library(tidyverse)
library(rjags)

# Read daily data
ddat <- read_csv("data/data_daily_cleaned.csv")
ddat$rainfall_dayback = ddat$`Rainfall Total`[c(2:length(ddat$`Rainfall Total`),NA)]

# Plot daily streamflow data with horizontal line showing flood threshold
ddat |> ggplot() +
  geom_line(aes(y = `Streamflow Ave`, x = as.Date(Date))) +
  # geom_hline(aes(yintercept = 4.076)) +
  ggtitle("Langrivier daily streamflow") +
  xlab("Date") +
  ylab("Streamflow (Cubic metres per second)")

# We need to withold some data for forecasting and validation
# We will withhold the data from the start of 2024 by making streamflow data NA
cal_ddat <- ddat |> mutate(`Streamflow Ave` = ifelse(Date < "2024-01-01", `Streamflow Ave`, NA)) %>% 
  filter(Date>"2013-03-16")

# #Identify missing values in rainfall
rain <- cal_ddat$rainfall_dayback

# # Track missing rainfall values
is_na_rain <- is.na(rain)
n_missing <- sum(is_na_rain)
missing_idx <- which(is_na_rain)

# ## DAILY data 
# # Format data for model
time <- cal_ddat$Date
y <- cal_ddat$`Streamflow Ave`
z <- cal_ddat$`Streamflow Ave` # for plotting later

#Add a seasonality component
doy <- as.numeric(format(time, "%j")) / 365  # day of year scaled 0–1
season_sin <- sin(2 * pi * doy)
season_cos <- cos(2 * pi * doy)

RandomWalk_rain_missing <- "
model {

  # Observation model
  for(t in 1:n){
    y[t] ~ dnorm(x[t], tau_obs)
  }

  # Process model with rainfall
  for(t in 2:n){
    x[t] ~ dnorm(
  beta_decay * x[t-1] + beta_rain * rain[t] +
  beta_season_sin * season_sin[t] + beta_season_cos * season_cos[t],
  tau_add)
  }

  # Impute missing rain values
  for(i in 1:n_missing){
    rain[missing_idx[i]] ~ dnorm(mu_rain, tau_rain)
  }

  # Priors
  beta_decay ~ dunif(0, 1)
  x[1] ~ dnorm(x_ic, tau_ic)
  tau_obs ~ dgamma(a_obs, r_obs)
  tau_add ~ dgamma(a_add, r_add)
  beta_rain ~ dnorm(0, 0.01)
  beta_season_sin ~ dnorm(0, 0.01)
  beta_season_cos ~ dnorm(0, 0.01)
  
  # Prior for rainfall distribution
  mu_rain ~ dnorm(0, 0.01)
  tau_rain ~ dgamma(1, 1)
}
"


data <- list(
  y = log(y),
  rain = rain,               # vector with NAs
  missing_idx = missing_idx,     # indices to impute
  n_missing = n_missing,         # how many to impute
  n = length(y),
  x_ic = log(0.1),
  tau_ic = 0.1,
  a_obs = 1,
  r_obs = 1,
  a_add = 1,
  r_add = 1, 
  season_sin = season_sin,      
  season_cos = season_cos )

# Run the model
j.model <- jags.model(file = textConnection(RandomWalk_rain_missing),
                      data = data,
                      n.chains = 3)

# First convergence check
jags.out <- coda.samples(model = j.model,
                         variable.names = c("tau_add", "tau_obs", "beta_rain", "beta_decay", "mu_rain", "tau_rain"),
                         n.iter = 2000)
plot(jags.out) # traceplot and density check


# Full posterior sampling
jags.out <- coda.samples(model = j.model,
                         variable.names = c("x", "tau_add", "tau_obs", "beta_rain", "beta_decay", "mu_rain", "tau_rain", "rain"),,
                         n.iter = 5000)


# Remove burn-in
burnin <- 1000
jags.burn <- window(jags.out, start = burnin)

# Convert MCMC output to matrix
out <- as.matrix(jags.out)

# Extract all 'x' columns (state estimates)
x.cols <- grep("^x\\[", colnames(out))

# Compute credible intervals on original (exp) scale
ci <- apply(exp(out[, x.cols]), 2, quantile, c(0.025, 0.5, 0.975))

# Plot setup
plot(time, ci[2,], type = 'n',
     ylim = range(z, na.rm = TRUE),
     ylab = "Streamflow average (m³/s)",
     log = "y",
     xlim = range(time),
     xlab = "Date",
     main = "Forecast with Rainfall Effect")

# Add credible interval ribbon
ecoforecastR::ciEnvelope(time, ci[1,], ci[3,], col = ecoforecastR::col.alpha("lightblue", 0.75))

# Add observed data
# Identify which data was used in the model (not NA)
included <- !is.na(log(y))
heldout <- is.na(log(y))

# Extract 'x' predictions from posterior matrix
x.cols <- grep("^x\\[", colnames(out))
ci <- apply(exp(out[, x.cols]), 2, quantile, c(0.025, 0.5, 0.975))

# Plot setup
plot(time, ci[2,], type = 'n',
     ylim = range(z, ci, na.rm = TRUE),
     ylab = "Streamflow average (m³/s)",
     log = "y",
     xlim = range(time),
     xlab = "Date",
     main = "Forecast with Rainfall + Seasonality")

# Add prediction intervals
ecoforecastR::ciEnvelope(time, ci[1,], ci[3,], col = ecoforecastR::col.alpha("lightblue", 0.75))

# Add observed data
points(time[included], z[included], pch = "+", col = 'black', cex = 0.6)  # model-used
points(time[heldout], z[heldout], pch = 1, col = 'red', cex = 0.8)       # held-out forecast

# Add scaled rainfall for context
par(new = TRUE)
plot(time, rain, type = "l", col = "darkgreen", axes = FALSE, xlab = "", ylab = "")
axis(4)
mtext("Rainfall (mm)", side = 4, line = 2)

