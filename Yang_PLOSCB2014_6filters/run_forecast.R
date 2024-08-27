## Script to run any of the filters for retrospective forecast of influenza
##  More information in 'Comparison of filtering methods for the modeling and retrospective forecasting of influenza epidemics' (PLoS Compute Biol)
##  by Wan Yang, Alicia Karspeck, and Jeffrey Shaman, 2014

## Home directories. CHANGE ACCORDINGLY
dir_home_code="~/code/";
dir_home_data="~/data/";

library("truncnorm"); library("tgp"); # for lhs
library("MASS"); # for multivariate normal distribution
require(plyr); # for function count
source(paste(dir_home_code,"EpiModelsTrackIncidence.R",sep="")); # SIRS model
source(paste(dir_home_code,"Fn_checkxnobounds.R",sep="")); # function to check DA aphysicality
## source of some funtions for initializations
source(paste(dir_home_code,'Fn_initializations.R',sep=""));
## source of different filters
source(paste(dir_home_code,'Fn_PF_rFC.R',sep=""));
source(paste(dir_home_code,'Fn_MIF_rFC.R',sep=""));
source(paste(dir_home_code,'Fn_pMCMC_rFC.R',sep=""));
source(paste(dir_home_code,'Fn_EAKF_rFC.R',sep=""));
source(paste(dir_home_code,'Fn_RHF_rFC.R',sep=""));
source(paste(dir_home_code,'Fn_EnKF_rFC.R',sep=""));

## read list of variables
if (! exists("phi")){
  phi=read.csv(paste(dir_home_data,'phi.csv',sep=""),header=F);
  params=phi[,1:4];
  susceps=infects=NULL; 
  for (i in 5:31){
    susceps=append(susceps,phi[,i])
    infects=append(infects,phi[,27+i])
  }
}
# global variables
N=1e5; # population
dt=1; # time step for the SIRS integration
tmstep=7; # wkly data
wk_start=40; # start from Week 40 of the year
## parameter boundaries:
D_low=1.5; L_low=1*365; Rmx_low=1.3; Rmn_low=0.8;
D_up=7; L_up=3650; Rmx_up=4; Rmn_up=1.2;
theta_low=c(L_low, D_low, Rmx_low, Rmn_low);
theta_up=c(L_up, D_up, Rmx_up, Rmn_up);
param.bound=cbind(theta_low,theta_up);
## parameters for the filters
disrete=T; # run the SIRS model discretely
metricsonly=F; # save all outputs
lambda=1.03; # inflation factor for the ensemble filters
num_ens=10000;  # if to run the particle filters, used 10000 particles, 300 for the ensemble filters

## read in data
# ah data for New York City 
ah=read.csv(paste(dir_home_data,'nyc_qclim79_02.csv',sep=""),header=F,sep=",")
AH=matrix(c(ah[,1],ah[,1]),365*2,1); # ah for NYC
## Read ILI+ data for New York City from files
iliiso=read.csv(paste(dir_home_data,"iliiso_nyc_2003wk40_2012wk47.csv",sep=""),header=F)
iliiso=as.matrix(iliiso,dim(iliiso)[1],dim(iliiso)[2]);
gamma=2.5; # use the same scaling for all cities & all seasons.

# try any epideimc season from 2003-04 ('03-04') to 2011-12 ('11-12'), exclude the pandemic
season="03-04"; # try the 2003-04 season
# get the observation data for that season
tmp=Fn_dates(season);
weeks=tmp$weeks; 
start_date=tmp$start_date; end_date=tmp$end_date;
obs_i=iliiso[weeks,1]*gamma 
# variance of ILI+ data
tmp=rep(0,length(obs_i))
for (i in 4:length(obs_i)){
  tmp[i]=mean(obs_i[(i-3):(i-1)]);
}
obs_vars=(1e4+(tmp^2)/50);

## get the first and last date of the simuliaton
clim_start=as.numeric(start_date-as.Date(paste("20",
           substr(season,gregexpr("-",season)[[1]][1]-2,gregexpr("-",season)[[1]][1]-1),"-01-01",sep="")))+1-6; 
# number of days in the year at the beginning of the week
clim_end=as.numeric(end_date-as.Date(paste("20",
            substr(season,gregexpr("-",season)[[1]][1]-2,gregexpr("-",season)[[1]][1]-1),"-01-01",sep="")))+1;

tm.ini=clim_start-1; # the end of the former week
tm.range=clim_start:clim_end;

# names of all functions
fn.names=c("PF_rFC","MIF_rFC","pMCMC_rFC","EnKF_rFC","EAKF_rFC","RHF_rFC")
## 
fn.id=1; # assign function ID, if fn.id=1, run the PF, ect..
fn=get(fn.names[fn.id]);
ntrn=15; # number of observation for the training period, could be 1 to 51 weeks
res=fn(num_ens, tmstep, param.bound, obs_i=obs_i,ntrn, 
       obs_vars,tm.ini, tm.range)

## plot results
out=res$train[,2:9]; out[,'newI']=out[,'newI'];
fcast=res$fcast[,2:5]; fcast[,'newI']=fcast[,'newI'];
par(mfrow=c(3,2), mar = c(4, 4, 1, 0), oma = c(0, 1, 4, 1),cex=.8, mgp=c(1.8,.5,0))
plot(out[,1],out[,2],xlim=c(tm.ini,clim_end),ylim=c(min(out[,2],fcast[,2])*.9,max(out[,2],fcast[,2])*1.1),
     ylab="# S per 100,000 population",xlab="Time",type="l")
lines(c(tail(out[,1],1),fcast[,1]),c(tail(out[,2],1),fcast[,2]),lty=5);
plot(out[,1],out[,3],xlim=c(tm.ini,clim_end),
     ylim=c(min(out[,3],fcast[,3])*.9,max(out[,3],fcast[,3])*1.1),
     ylab="# I per 100,000 population",xlab="Time",type="l")
lines(c(tail(out[,1],1),fcast[,1]),c(tail(out[,3],1),fcast[,3]),lty=5);
plot(out[,1],out[,8],xlim=c(tm.ini,clim_end),
     ylim=c(0,max(out[,8],obs_i,fcast[,4])*1.1),
     ylab="# New cases per 100,000 patient visits",xlab="Time",type="l")
points(c(out[,1],fcast[,1]),obs_i,pch="x",cex=.9);
lines(c(tail(out[,1],1),fcast[,1]),c(tail(out[,8],1),fcast[,4]),lty=5);
legend("topright",c("modeled","observed"),lty=c(1,NA),pch=c(NA,"x"),cex=.6,bty='n')
plot(out[,1],out[,4],xlim=c(tm.ini,clim_end),ylab="L",xlab="Time",type="l")
lines(fcast[,1],rep(tail(out[,4],1),length(fcast[,1])),lty=5);
plot(out[,1],out[,5],xlim=c(tm.ini,clim_end),ylab="D",xlab="Time",type="l")
lines(fcast[,1],rep(tail(out[,5],1),length(fcast[,1])),lty=5);
plot(out[,1],out[,6],xlim=c(tm.ini,clim_end),ylab="R",ylim=c(0,max(out[,6])+.5),
     xlab="Time",type="l",col="blue")
lines(fcast[,1],rep(tail(out[,6],1),length(fcast[,1])),lty=5,col="blue");
lines(out[,1],out[,7],col="green");
lines(fcast[,1],rep(tail(out[,7],1),length(fcast[,1])),lty=5,col="green");
legend("topright",c("R0max","R0min"),lty=c(1,1),col=c("blue","green"),cex=.6,bty='n')
mtext(bquote(.(fn.names[fn.id])*","~.(num_ens)~"particles; ILI+:"~.(season)~"season;"~
               "Forecast from Wk"~.(40+ntrn-1)),side = 3, outer = TRUE, cex = 0.8, line = 0, col = "black")
