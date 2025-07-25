### Streamflow data analysis

# Load packages
library(tidyverse)
library(rjags)

# Read daily data
ddat <- read_csv("data/data_daily_cleaned.csv")

# Read hourly data
hdat <- read_csv("data/data_hourly_cleaned.csv")

# Plot daily streamflow data with horizontal line showing flood threshold
ddat |> ggplot() +
  geom_line(aes(y = `Streamflow Ave`, x = as.Date(Date))) +
  geom_hline(aes(yintercept = 4.076)) +
  ggtitle("Langrivier daily streamflow") +
  xlab("Date") +
  ylab("Streamflow (Cubic metres per second)")

# Plot hourly streamflow data with horizontal line showing flood threshold
hdat |> ggplot() +
  geom_line(aes(y = `Streamflow`, x = as.Date(Date))) +
  geom_hline(aes(yintercept = 4.076)) +
  ggtitle("Langrivier hourly streamflow") +
  xlab("Date") +
  ylab("Streamflow (Cubic metres per second)")

# We need to withold some data for forecasting and validation
# We will withhold the data from the start of 2024 by making streamflow data NA
cal_ddat <- ddat |> mutate(`Streamflow Ave` = ifelse(Date < "2024-01-01", `Streamflow Ave`, NA))
cal_hdat <- hdat |> mutate(`Streamflow` = ifelse(Date < "2024-01-01", `Streamflow`, NA))

# save the calibration data
# write_csv(cal_ddat, "data/data_daily_calibration.csv")
# write_csv(cal_hdat, "data/data_hourly_calibration.csv")

## We are going to use DAILY data for our forecasts

# We need to withold some data for forecasting and validation
cal_ddat <- cal_ddat %>%
  filter(Date > "2013-03-16") # remove dates before 2013-03-16 where there is no rainfall data

ddat <- ddat %>%
  filter(Date > "2013-03-16") # remove dates before 2013-03-16 where there is no rainfall data - need to do this for plotting

# need to move previous day's rainfall to current day 
cal_ddat <- cal_ddat %>%
  mutate(rainfall_dayback = lag(`Rainfall Total`, 1))

# Format data for model
time <- cal_ddat$Date
y <- cal_ddat$`Streamflow Ave`
z <- ddat$`Streamflow Ave` # for plotting later
y_log <- log(y)
y_log[is.infinite(y_log)] <- NA
#Keep the observed streamflow values before removing the held out portion
y_full <- log(ddat$`Streamflow Ave`)
y_full[is.infinite(y_full)] <- NA


## Start with a null time-series model using JAGS
# Define the model
RandomWalk <- "
model{
  
  #### Data Model
  for(t in 1:n){
    y[t] ~ dnorm(x[t],tau_obs)
  }
  
  #### Process Model (random walk)
  for(t in 2:n){
    x[t] ~ dnorm(x[t-1], tau_add)
  }

  #### Priors
  x[1] ~ dnorm(x_ic, tau_ic)
  tau_obs ~ dgamma(a_obs,r_obs) ## prior on observation error
  tau_add ~ dgamma(a_add,r_add) ## prior on process error
}
"

data <- list(y=y_log,n=length(y),      ## data
             x_ic=log(0.1),tau_ic=0.1,    ## initial condition prior
             a_obs=1,r_obs=1,           ## obs error prior
             a_add=1,r_add=1            ## process error prior
)

nchain = 3
# init <- list()
# for(i in 1:nchain){
#   y.samp = sample(y,length(y),replace=TRUE)
#   init[[i]] <- list(tau_add=1/var(diff(log(y.samp))),  ## initial guess on process precision
#                     tau_obs=5/var(log(y.samp)))        ## initial guess on obs precision

# Run the model
j.model   <- jags.model (file = textConnection(RandomWalk),
                         data = data,
                         # inits = init,
                         n.chains = 3)

# Sample from model without x to check that convergence has happened
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("tau_add","tau_obs"),
                            n.iter = 5000)
# See if convergence has happened
plot(jags.out)

jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add","tau_obs"),
                            n.iter = 10000)

burnin = 1000                                   ## determine convergence
jags.burn <- window(jags.out, start = burnin)  ## remove burn-in

# Plot data and confidence interval
time.rng = c(1,length(time))       ## adjust to zoom in and out
out <- as.matrix(jags.out)         ## convert from coda to matrix  
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(exp(out[,x.cols]),2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

plot(time,ci[2,],type='n',ylim=range(y,na.rm=TRUE),ylab="Streamflow average", log='y', xlim=time[time.rng])
## adjust x-axis label to be monthly if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), format = "%Y-%m")
}
ecoforecastR::ciEnvelope(time,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75)) # add confidence interval 
# add data points
included <- !is.na(y)
heldout <- is.na(y)
# Plot included data points (model saw these)
points(time[included], y[included], pch="+", col='black', cex=0.6)  # filled black dots
# Plot held-out data points (model did NOT see these)
points(time[heldout], z[heldout], pch=1, col='red', cex=0.8)       # open red circles 

## Now let's see if we can add covariates to the model
# We will use the daily rainfall data
rainfall_lag <- cal_ddat$rainfall_dayback 
rainfall_lag[is.na(rainfall_lag)] <- 0 # for now make NA values 0

# Define the model
Rainfall_RandomWalk <- "
model{
  
  #### Data Model
  for(t in 1:n){
    y[t] ~ dnorm(x[t],tau_obs)
  }
  
  #### Process Model (random walk)
  for(t in 2:n){
    x[t] ~ dnorm(x[t-1] + beta * c[t], tau_add)
  }

  #### Priors
  x[1] ~ dnorm(x_ic, tau_ic)
  beta ~ dnorm(0, 0.01) ## prior on the beta for rainfall
  tau_obs ~ dgamma(a_obs,r_obs) ## prior on observation error
  tau_add ~ dgamma(a_add,r_add) ## prior on process error
}
"

data <- list(y=y_log,n=length(y),      ## data
             c = rainfall_lag,  ## rainfall
             x_ic=log(0.1),tau_ic=0.1,    ## initial condition prior
             a_obs=1,r_obs=1,           ## obs error prior
             a_add=1,r_add=1            ## process error prior
)


# Run the model
j.model   <- jags.model (file = textConnection(Rainfall_RandomWalk),
                         data = data,
                         # inits = init,
                         n.chains = 3)

jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add","tau_obs"),
                            n.iter = 10000)

burnin = 1000                                   ## determine convergence
jags.burn <- window(jags.out, start = burnin)  ## remove burn-in

# Plot data and confidence interval
time.rng = c(1,length(time))       ## adjust to zoom in and out
out <- as.matrix(jags.out)         ## convert from coda to matrix  
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(exp(out[,x.cols]),2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

plot(time,ci[2,],type='n',ylim=range(y,na.rm=TRUE),ylab="Streamflow average", log='y', xlim=time[time.rng])
## adjust x-axis label to be monthly if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), format = "%Y-%m")
}
ecoforecastR::ciEnvelope(time,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75)) # add confidence interval 
# add data points
included <- !is.na(y)
heldout <- is.na(y)
# Plot included data points (model saw these)
points(time[included], y[included], pch="+", col='black', cex=0.6)  # filled black dots
# Plot held-out data points (model did NOT see these)
points(time[heldout], z[heldout], pch=1, col='red', cex=0.8)       # open red circles 


# Now let's try model the missing rainfall data and add a seasonality component to the model

#Add a seasonality component
doy <- as.numeric(format(time, "%j")) / 365  # day of year scaled 0–1
season_sin <- sin(2 * pi * doy)
season_cos <- cos(2 * pi * doy)

# #Identify missing values in rainfall
rain <- cal_ddat$rainfall_dayback

# # Track missing rainfall values
is_na_rain <- is.na(rain)
n_missing <- sum(is_na_rain)
missing_idx <- which(is_na_rain)

RandomWalk_rain_decay <- "
model {

  # Observation model
  for(t in 1:n){
    y[t] ~ dnorm(x[t], tau_obs)
  }

  # Process model with autoregressive decay and covariates
  for(t in 2:n){
    mu[t] <- mu0 + beta_decay * (x[t-1] - mu0) + 
             beta_rain * rain[t] + 
             beta_season_sin * season_sin[t] + 
             beta_season_cos * season_cos[t]
    
    x[t] ~ dnorm(mu[t], tau_add)
  }

  # Impute missing rain values
  for(i in 1:n_missing){
    rain[missing_idx[i]] ~ dnorm(mu_rain, tau_rain)
  }

  # Priors
  mu0 ~ dnorm(0, 0.001)                     # Mean log-streamflow level
  x[1] ~ dnorm(mu0, tau_ic)                 # Initial latent state
  
  tau_obs ~ dgamma(a_obs, r_obs)            # Observation error
  tau_add ~ dgamma(a_add, r_add)            # Process error

  beta_decay ~ dunif(0, 1)                  # AR(1) coefficient bounded for stability
  beta_rain ~ dnorm(0, 0.01)
  beta_season_sin ~ dnorm(0, 0.01)
  beta_season_cos ~ dnorm(0, 0.01)

  mu_rain ~ dnorm(0, 0.01)                  # Mean log-rainfall for imputation
  tau_rain ~ dgamma(1, 1)                   # Rainfall imputation variance
}
"

data <- list(
  y = y_log,
  rain = log(rain+1),               # vector with NAs
  missing_idx = missing_idx,     # indices to impute
  n_missing = n_missing,         # how many to impute
  n = length(y),
  #  x_ic = log(0.1),
  tau_ic = 0.1,
  a_obs = 1,
  r_obs = 1,
  a_add = 1,
  r_add = 1, 
  season_sin = season_sin,      
  season_cos = season_cos )

# Run the model
j.model <- jags.model(file = textConnection(RandomWalk_rain_decay),
                      data = data,
                      n.chains = 3)

# First convergence check
jags.out <- coda.samples(model = j.model,
                         variable.names = c("tau_add", "tau_obs", "beta_rain", "beta_decay", "mu_rain", "tau_rain"),
                         n.iter = 1000)
plot(jags.out) # traceplot and density check


# Full posterior sampling
jags.out <- coda.samples(model = j.model,
                         variable.names = c("x", "tau_add", "tau_obs", "beta_rain", "beta_decay", "mu_rain", "tau_rain", "rain", "mu0", "beta_season_sin", "beta_season_cos"),,
                         n.iter = 5000)

# Remove burn-in
burnin <- 1000
jags.burn <- window(jags.out, start = burnin)

# Plot data and confidence interval
time.rng = c(1,length(time))       ## adjust to zoom in and out
out <- as.matrix(jags.out)         ## Convert MCMC output to matrix
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(exp(out[,x.cols]),2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

plot(time,ci[2,],type='n',ylim=range(y,na.rm=TRUE),ylab="Streamflow average (m³/s)", log='y', xlim=time[time.rng], xlab = "Date",
     main = "Forecast with Rainfall + Seasonality + Decay")

## adjust x-axis label to be monthly if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), format = "%Y-%m")
}

ecoforecastR::ciEnvelope(time,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightblue",0.75)) # add confidence interval 

# add line for mean prediction
forecast_period <- time >= as.Date("2024-01-01")
lines(time[forecast_period], ci[2, forecast_period], col = "blue", lwd = 2)

# Plot observed data
included <- !is.na(y)
heldout <- is.na(y) & time >= as.Date("2024-01-01")
# Plot included data points (model saw these)
points(time[included], y[included], pch="+", col='black', cex=0.6)  # filled black dots
# Plot held-out data points (model did NOT see these)
points(time[heldout], z[heldout], pch=1, col='red', cex=0.8)       # open red circles 


# Let's plot just the last year before predicting
forecast_start <- as.Date("2024-01-01") # Define the start of the forecast period
forecast_end <- max(time, na.rm = TRUE)  # Adjust if you want a specific cutoff
plot_start <- forecast_start - 365

# Logical vector to subset the full time range
plot_range <- time >= plot_start & time <= forecast_end

plot(time[plot_range], ci[2, plot_range], type = 'n',
     ylim = range(y, na.rm = TRUE),
     ylab = "Streamflow average (m³/s)",
     log = 'y',
     xlim = c(plot_start, forecast_end),
     xlab = "Date",
     main = "Forecast with Rainfall + Seasonality + Decay")
# Confidence envelope
ecoforecastR::ciEnvelope(time[plot_range], ci[1, plot_range], ci[3, plot_range],
                         col = ecoforecastR::col.alpha("lightblue", 0.75))

# Mean forecast line (for forecast period only)
forecast_period <- time >= forecast_start
lines(time[forecast_period], ci[2, forecast_period], col = "blue", lwd = 2)

# Observed data points
included <- !is.na(y)
points(time[plot_range & included], y[plot_range & included], pch = "+", col = 'black', cex = 0.6)

# Held-out points
heldout <- is.na(y)
points(time[plot_range & heldout], z[plot_range & heldout], pch = 1, col = 'red', cex = 0.8)