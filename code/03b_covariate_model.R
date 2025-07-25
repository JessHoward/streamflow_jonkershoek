## Adding rainfall covariate to the model using the ecoforecastR package

# Load packages
library(tidyverse)
library(rjags)


# We will use the daily rainfall data as a covariate but the 
# Define the model
RandomWalk = "
model{
  
  #### Data Model
  for(t in 1:n){
    y[t] ~ dnorm(x[t],tau_obs)
  }
  
  #### Process Model
  for(t in 2:n){
    x[t]~dnorm(x[t-1],tau_add)
  }
  
  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic) ## priors on the initial conditions (since this is unknown)
  tau_obs ~ dgamma(a_obs,r_obs) ## prior on observation error
  tau_add ~ dgamma(a_add,r_add) ## prior on process error
}
"

# Read daily streamflow data
cal_ddat <- read_csv("data/data_daily_calibration.csv")
ddat <- read_csv("data/data_daily_cleaned.csv")
#Choose the covariant as Rainfall total

# Format data for model
rain <- c(cal_ddat$`Rainfall Total`) # today's rainfall 
#rain <- c(NA, cal_ddat$`Rainfall Total`[-length(ddat$`Rainfall Total`)]) # rainfall a day forward 
#rain <- c(cal_ddat$`Rainfall Total`[-1], NA) # rainfall a day back
time <- cal_ddat$Date
y <- cal_ddat$`Streamflow Ave`
z <- ddat$`Streamflow Ave` # for held-out data points (to plot later)
data <- list(y=log(y),n=length(y),      ## data
             x_ic=log(0.1),tau_ic=0.1,    ## initial condition prior
             a_obs=1,r_obs=1,           ## obs error prior
             a_add=1,r_add=1            ## process error prior
)


## fit the model
ef.out <- ecoforecastR::fit_dlm(model=list(obs="y",fixed="~ 1 + X + rain"),data)
names(ef.out)

## parameter diagnostics
params <- window(ef.out$params,start=1000) ## remove burn-in
plot(params)
summary(params)
cor(as.matrix(params))
pairs(as.matrix(params))


##Confidence interval
time.rng = c(1,length(time))       ## adjust to zoom in and out
out <- as.matrix(ef.out$predict)         ## convert from coda to matrix  
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(exp(out),2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

plot(time,ci[2,],type='n',ylim=range(y,na.rm=TRUE),ylab="Streamflow",log='y',xlim=time[time.rng])
## adjust x-axis label to be daily if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='day'), format = "%Y-%m")
}
ecoforecastR::ciEnvelope(time,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
# add data points
included <- !is.na(y)
heldout <- is.na(y)
# Plot included data points (model saw these)
points(time[included], y[included], pch="+", col='black', cex=0.6)  # filled black dots
# Plot held-out data points (model did NOT see these)
points(time[heldout], z[heldout], pch=1, col='red', cex=0.8)       # open red circles 

#JAGS model returned and viewed
strsplit(ef.out$model,"\n",fixed = TRUE)[[1]]

