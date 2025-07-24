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
    x[t] ~ dnorm(x[t-1],tau_add)
  }

  #### Priors
  x[1] ~ dnorm(x_ic, tau_ic)
  tau_obs ~ dgamma(a_obs,r_obs) ## prior on observation error
  tau_add ~ dgamma(a_add,r_add) ## prior on process error
}
"

## DAILY data 
# Format data for model
time <- cal_ddat$Date
y <- cal_ddat$`Streamflow Ave`
z <- ddat$`Streamflow Ave` # for plotting later

data <- list(y=log(y),n=length(y),      ## data
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
points(time[included], z[included], pch="+", col='black', cex=0.6)  # filled black dots
# Plot held-out data points (model did NOT see these)
points(time[heldout], z[heldout], pch=1, col='red', cex=0.8)       # open red circles 

## Now let's see if we can add covariates to the model
# We will use the daily rainfall data as a covariate but the 
