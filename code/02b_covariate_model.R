## Now let's see if we can add covariates to the model
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

# Load packages
library(tidyverse)

# Read daily streamflow data
ddat <- read_csv("data/data_daily_calibration.csv")

#Choose the covariant as Rainfall total
#Create a one day lag between rain falling and streamflow increase
ddat$rainfall_dayback = ddat$`Rainfall Total`[c(2:length(ddat$`Rainfall Total`),NA)]

#Remove the last row of NA because this would be rain data but no streamflow data
ddat <- ddat[-nrow(ddat), ]


# Format data for model
rain <- ddat$rainfall_dayback
time <- ddat$Date
y <- ddat$`Streamflow Ave`

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
points(time,y,pch="+",cex=0.5)

#JAGS model returned and viewed
strsplit(ef.out$model,"\n",fixed = TRUE)[[1]]
