## Model Assessment
## load libraries
library("plotrix")
library(rpart)
library(randomForest)
library(scoringRules)

## load SSEM output
load("data/SFT_output.RData") # SFT need to update to our outputs
load("data/SFT_output.pf.RData") # SFT need to update to our outputs

## load streamflow data
stream = read.csv("data/data_hourly_2025-07-09.csv",header=TRUE)
#L4[L4==-9999] = NA

## Calculate ensemble means & apply QAQC
qaqc = (stream$qf_NEE_st == 0)
NEE.ens.bar = -apply(NEE.ens[,],1,mean)
NEE.pf.bar  = -apply(NEE.pf[,],1,mean)
E = NEE.ens.bar[qaqc]
P = NEE.pf.bar[qaqc]
O = stream$NEE_st_fMDS[qaqc]

## Model vs obs regressions
NEE.ens.fit = lm(O ~ E)
NEE.pf.fit = lm(O ~ P)

## performance stats
stats = as.data.frame(matrix(NA,4,2))
rownames(stats) <- c("RMSE","Bias","cor","slope")
colnames(stats) <- c("ens","pf")
stats["RMSE",'ens'] = sqrt(mean((E-O)^2))
stats["RMSE",'pf']  = sqrt(mean((P-O)^2))
stats['Bias','ens'] = mean(E-O)
stats['Bias','pf']  = mean(P-O)
stats['cor','ens']  = cor(E,O)
stats['cor','pf']   = cor(P,O)
stats['slope','ens'] = coef(NEE.ens.fit)[2]
stats['slope','pf']  = coef(NEE.pf.fit)[2]
knitr::kable(stats)

## predicted-observed
plot(E,O,pch=".",xlab="ensemble",ylab='observed',main='NEE (umol/m2/sec)')
abline(0,1,col=2,lwd=2)
abline(NEE.ens.fit,col=3,lwd=3,lty=2)
legend("bottomright",legend=c('obs','1:1','reg'),col=1:3,lwd=3)

plot(P,O,pch=".",xlab="particle filter",ylab='observed',main='NEE (umol/m2/sec)')
abline(0,1,col=2,lwd=2)
abline(NEE.pf.fit,col=3,lwd=3,lty=2)
legend("bottomright",legend=c('obs','1:1','reg'),col=1:3,lwd=3)


# Which model performed better? What parameters might need to be fixed, or processeses refined in the model?
#Plot a time-series of observed and modeled **daily mean** streamflow
#Make sure to use the gap-filled streamflow estimates, since data are not missing at random. Hint: you can use `tapply` and the day of year column in the data

## Taylor diagrams: Ensemble
taylor.diagram(ref=O,model=E,normalize=TRUE,ref.sd=TRUE)

## add full ensemble
for(i in 1:ncol(NEE.pf)){
  taylor.diagram(ref=O,model=-NEE.ens[qaqc,i],col=2,pch=".",add=TRUE,normalize=TRUE)
}

# Taylor Diagram: particle filter
taylor.diagram(ref=O,model=P,add=TRUE,normalize=TRUE,col=3)

## add data uncertainty
rlaplace = function(n,mu,b){
  return(mu + ifelse(rbinom(n,1,0.5),1,-1)*rexp(n,b))
}
beta = ifelse(O > 0,0.62+0.63*O,1.42-0.19*O) #Heteroskedasticity, parameters from Richardson et al 2006
for(i in 1:200){
  x = rlaplace(length(O),O,beta)
  taylor.diagram(ref=O,model=x,col=5,add=TRUE,normalize=TRUE)
}
legend("topright",legend=c("ens","PF","obsUncert"),col=2:5,pch=20,cex=0.7)

#What did you learn about model performance from the Taylor diagram? 

#Bayesian p-values
O = stream$NEE_st_fMDS  ## observed
pval.pf = 0
for(i in 1:nrow(NEE.pf)){
  pval.pf[i] = sum(O[i] > -NEE.pf[i,])/ncol(NEE.pf)  ## quantiles of the particle filter
}
plot(pval.pf)   ## quantile 'residuals'
hist(pval.pf,probability=TRUE) ## quantile distribution (should be flat)

crps = scoringRules::crps_sample(y = O, dat = NEE.pf)
time = as.Date(stream$DoY,origin = "2005-01-01")
mean(crps)
plot(time,crps)
hist(crps)

### diurnal cycle of model skill
crps.diurnal = tapply(crps,stream$Hour,mean)
time.of.day = as.numeric(names(crps.diurnal))
plot(time.of.day,crps.diurnal)

## CRPS skill as a function of observed NEE
plot(O,crps)

#In the final section we'll use a few off-the-shelf data mining approaches to look at the model residuals and ask what parts of our input space are associated with the largest model error. Note that we are not limited to just examining the effects of the model inputs, we might also look at other potential drivers that are not included in our model, such as soil moisture, to ask if model error is associated with our failure to include this (or other) drivers. Alternatively, we could have looked at other factors such as the time of day or even other model variables (e.g. is model error higher when LAI is larger or small?)
#Of the many algorithms out there we'll look at two: the Classification and Regression Tree (CART) model and the Random Forest model. For both we'll define our error metric as $(E-O)/beta$, where beta is the parameter equivalent to the variance in Laplace distribution. Specifically, we're using the heteroskedastic observation error to reweight the residuals to account for the fact that large residuals at times of high flux is likely due to high measurement error. Thus the errors can be interpreted as similar to the number of of standard deviations.
#The CART model is a classification algorithm which will build a tree that discretely classifies when the model has high and low error.
#The Random Forest model is more like a response surface. The Random Forest will generate 'partial dependence' plots, which indicate the importance of each factor across its range, as well as an overall estimate of the importance of each factor in the model error. 
#The key thing to remember in all these plots is that we're modelling the RESIDUALS in order to diagnose errors, not modeling the NEE itself.

{r, fig.height=6}
## define error metric and dependent variables
O = stream$NEE_st_fMDS[qaqc]
err = (E-O)/beta
x = cbind(inputs$PAR[qaqc],inputs$temp[qaqc])
colnames(x) = c("PAR","temp")
smp = sample.int(length(err),1000)  ## take a sample of the data since some alg. are slow

### Classification tree
rpb = rpart(err ~ x) ## bias
plot(rpb)
text(rpb)
e2 = err^2
rpe = rpart(e2 ~ x) ## sq error
plot(rpe)
text(rpe)

## Random Forest
rfe = randomForest(x[smp,],abs(err[smp]))
rfe$importance
partialPlot(rfe,x[smp,],"PAR")
partialPlot(rfe,x[smp,],"temp")

#Overall, which driver is most important in explaining model error? What conditions are most associated with model success? With model failure?  Where do these results reinforce conclusions we reached earlier and where do they shine light on new patterns you may have missed earlier?
  
#In this section we look at how well the model performed by assessing the modeled relationships between inputs and outputs and comparing that to the same relationship in the data. The raw relationships are very noisy, as many covariates are changing beyond just the single input variable we are evaluating, so in addition we calculate binned means for both the model and data.
## raw
plot(inputs$temp[qaqc],O,pch=".",ylab="NEE")
points(inputs$temp[qaqc],E,pch=".",col=2)

## binned
nbin = 25
#PAR = inputs$PAR[qaqc]
#x = seq(min(PAR),max(PAR),length=nbin)
Tair = inputs$temp[qaqc]
xd = seq(min(Tair),max(Tair),length=nbin)
xmid = xd[-length(xd)] + diff(xd)
bin = cut(Tair,xd)
Obar = tapply(O,bin,mean,na.rm=TRUE)
Ose  = tapply(O,bin,std.error,na.rm=TRUE)
Ebar = tapply(E,bin,mean,na.rm=TRUE)
Ese  = tapply(E,bin,std.error,na.rm=TRUE)
OCI = -cbind(Obar-1.96*Ose,Obar,Obar+1.96*Ose)
ECI = -cbind(Ebar-1.96*Ese,Ebar,Ebar+1.96*Ese)
rng = range(rbind(OCI,ECI))

col2=ecoforecastR::col.alpha("darkgrey",0.9)
col1=ecoforecastR::col.alpha("lightgrey",0.6)

plot(xmid,Obar,ylim=rng,type='n',xlab="Air Temperature (C)",ylab="NEP (umol/m2/s)",cex.lab=1.3)
ecoforecastR::ciEnvelope(xmid,ECI[,1],ECI[,3],col=col2)
lines(xmid,ECI[,2],col="black",lwd=4)
ecoforecastR::ciEnvelope(xmid,OCI[,1],OCI[,3],col=col1)
lines(xmid,OCI[,2],col="lightgrey",lwd=4)

legend("topleft",legend=c("Model","Data"),lwd=10,col=c(col2,col1),lty=1,cex=1.7)


#Question 8 [A]: Evaluate the model's ability to capture functional responses to both Temperature and PAR.

#Overall 
### other summary figures to go in multi-panel
par(mfrow=c(2,2))

## Time-series visualization, daily means
DoY = floor(stream$DoY-0.02)
uDoY = sort(unique(DoY))
ci.pf  = apply(apply(NEE.pf[,],2,tapply,DoY,mean),1,mean)
NEE = -stream$NEE_st_fMDS
NEEd = tapply(NEE,DoY,mean)
plot(uDoY,ci.pf,xlab="time",ylab="NEE",type='l',ylim=range(c(ci.pf,NEEd)),cex.lab=1.3)
points(uDoY,NEEd,col=2,pch="+")
legend("topright",legend=c("Model","Data"),lty=c(1,NA),pch=c(NA,"+"),col=1:2,cex=1.3)

## predicted vs observed
plot(NEEd,ci.pf,xlab="Model",ylab="Data",cex.lab=1.3)
abline(0,1,lty=2,lwd=4)
abline(lm(ci.pf ~ NEEd),col=2,lwd=3,lty=3)
legend("topleft",legend=c("1:1","Reg"),lty=2:3,lwd=4,col=1:2,cex=1.3)

## Functional response
plot(xmid,Obar,ylim=rng,type='n',xlab="Air Temperature (C)",ylab="NEP (umol/m2/s)",cex.lab=1.3)
ecoforecastR::ciEnvelope(xmid,ECI[,1],ECI[,3],col=col2)
lines(xmid,ECI[,2],col="black",lwd=4)
ecoforecastR::ciEnvelope(xmid,OCI[,1],OCI[,3],col=col1)
lines(xmid,OCI[,2],col="lightgrey",lwd=4)

legend("bottom",legend=c("Model","Data"),lwd=10,col=c(col2,col1),lty=1,cex=1.3)

### Classification tree
par(mar=c(0,0,0,0))
rpe = rpart(e2 ~ PAR+Tair,as.data.frame(x),method="anova") ## sq error
plot(rpe,margin=0.1)
text(rpe,cex=1.5)


# Extra Material

Comparison to flux "climatology"
-------------------------------
  
  In the section below we calculate the long-term average NEE for each 30 min period in the year, excluding the year we modeled (2005) as an alternative model to judge our process model against. We then update our summary statistics and predicted-observed plot

```{r}
## flux "climatology"
fluxfiles = dir("data",pattern="AMF")
fluxfiles = fluxfiles[grep("txt",fluxfiles)]
fluxfiles = fluxfiles[-grep("2005",fluxfiles)]
clim.NEE = clim.doy = NULL
for(f in fluxfiles){
  ff = read.csv(file.path("data",f),header=TRUE,na.strings="-9999")
  ff[ff == -9999] = NA
  clim.NEE = c(clim.NEE,ff$NEE_st_fMDS)
  clim.doy = c(clim.doy,ff$DoY)
}
NEE.clim=tapply(clim.NEE,clim.doy,mean,na.rm=TRUE)[1:length(qaqc)]
C = NEE.clim[qaqc]
NEE.clim.fit = lm(O ~ C)
summary(NEE.clim.fit)
stats["RMSE",3]  = sqrt(mean((C-O)^2))
stats['Bias',3]  = mean(C-O)
stats['cor',3]   = cor(C,O)
stats['slope',3] = coef(NEE.clim.fit)[2]
colnames(stats)[3] <- "clim"
knitr::kable(stats)
plot(C,O,pch=".",xlab="climatology",ylab='observed',main='NEE (umol/m2/sec)')
abline(0,1,col=2,lwd=2)
abline(NEE.clim.fit,col=3,lwd=3,lty=2)
legend("bottomright",legend=c('obs','1:1','reg'),col=1:3,lwd=3)

## example cycle
plot(stream$DoY,-stream$NEE_st_fMDS,xlim=c(200,210),type='l',lwd=2,ylim=c(-10,20),xlab="Day of Year",ylab="NEE")
lines(stream$DoY,-NEE.clim,col=4,lwd=2,lty=2)
legend("topright",legend=c("Obs","clim"),lty=1:2,col=c(1,4),lwd=2)

#How does the process model perform relative to the average flux data? Which statistics showed the largest differences between the model and climatology? 
#Time-scales
#In the next section we look at the average diurnal cycle of the data and models.

## diurnal cycle
NEE.ens.diurnal = tapply(E,stream$Hour[qaqc],mean)
NEE.pf.diurnal  = tapply(P,stream$Hour[qaqc],mean)
NEE.clim.diurnal  = tapply(C,stream$Hour[qaqc],mean)
NEE.obs.diurnal = tapply(O,stream$Hour[qaqc],mean)
ylim=range(c(NEE.ens.diurnal,NEE.pf.diurnal,NEE.obs.diurnal))
tod = sort(unique(stream$Hour))
plot(tod,NEE.ens.diurnal,ylim=ylim,col=2,xlab="Time of Day",ylab='NEE',main="Diurnal Cycle",type='l',lwd=3)
lines(tod,NEE.pf.diurnal,col=3,lwd=3)
lines(tod,NEE.clim.diurnal,col=4,lwd=3)
lines(tod,NEE.obs.diurnal,lwd=3)
legend("bottomright",legend=c("obs","ens","PF","clim"),col=1:4,pch=20,cex=0.75)

