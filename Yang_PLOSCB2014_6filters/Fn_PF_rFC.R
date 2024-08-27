##  R Function to paricile filter, with resampling and regularization, and then forecast
##  Resampling is initiated when Neff<N/2.
##  Based on Arulampalam MS, Maskell S, Gordon N, Clapp T (2002) IEEE Trans Signal Process 50: 174-188.
##           and Shaman J, Karspeck A (2012) Proc Natl Acad Sci USA 109: 20425-20430.
##  use Eqn 63 in Arulampalam et al. 2002, where the importance density is simply taken to be the prior.
##  More information in 'Comparison of filtering methods for the modeling and retrospective forecasting of influenza epidemics' (PLoS Compute Biol)
##  by Wan Yang, Alicia Karspeck, and Jeffrey Shaman, 2014

PF_rFC<-function(num_ens, tmstep, param.bound, obs_i=obs_i, ntrn=1,
                          obs_vars,tm.ini=273, tm.range=273:500,regul=TRUE){
  ## PF to run retrospecitve forecast of influenza
  ## Inputs: 
  ##        num_ens: number of ensemble members
  ##        tmstep: time step (e.g., 7 days if weekly observations are used)
  ##        param.bound: bounds for the parameters 
  ##          (a 2-column matrix prescribed the lower (col#1) and upper (col#2 of the bounds))
  ##        obs_i: observations (the time series);
  ##        ntrn: number of observations for training the filter
  ##        obs_vars: variances of the obseravations
  ##        tm.ini: the initial time (the number of day in the year)
  ##        tm.range: the sequence of day from the first to the last day of simulation 
  ##         from the current year to the next year of a season (used to match the time of climatological humidity)
  ## Outputs:
  ##        train: modeled time series of the variables/parameters (posteriors) for the training period
  ##        trainsd: standard deviation of the variables/parameters during the training
  ##        fcast: forecasted time series of the variables (i.e., S, I, and newI)
  ##        fcsd: standard deviation of the variables of the forecast
  ##        metrics: root mean squared error (rms), correlation with the observations, ect.
  ##        Note: if only wish to save the metrics, set metricsonly=TRUE
  
  library(tgp);
  source(paste(dir_home_code,"EpiModelsTrackIncidence.R",sep="")); # SIR model
  source(paste(dir_home_code,"Fn_checkxnobounds.R",sep="")); # function to check DA aphysicality
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
  
  num_times=length(obs_i);  
  num_var=2+4+1; # include S, I, L, D, R0max, R0min, newI
  nsn=length(obs_i);
  nfc=nsn-ntrn; # number of weeks for forecasting;
  fcast=array(0,c(3,num_ens,nfc)); # fcast: S, I, newI;
  
  So=matrix(0,6,num_ens);  
  
  x=array(0,c(num_var,num_ens,ntrn+1));
  # y=matrix(0,num_ens,ntrn+1); # new I
  # Y=rep(0,ntrn);  # mean new I
  
  theta_low=param.bound[,1];
  theta_up=param.bound[,2];
  
  paramsEns=lhs(num_ens,param.bound)
  
  # Initialization:
  
  rnd=cbind(ceiling(27*1e4*runif(num_ens)), ceiling(27*1e4*runif(num_ens)));
  So[1,]=susceps[rnd[,1]]/5;
  So[2,]=infects[rnd[,2]]/5;
  So[3:6,]=t(paramsEns);
  
  # integrate forward 1 step parallelly
  beta.range=tm.range[1]:(tail(tm.range,1)+2*tmstep); 
  # make sure to include enought beta, b/c the prior integrate 1 time step more than needed
  AHpt=AH[beta.range,]; # only use AH during the time of interest, to minimize space needed
  ones1=matrix(1,length(AHpt),1); # matrix of 1 to make the dimension match
  b=log(So[5,]-So[6,]); b=ones1%*%b; # expand b to a matrix with each col for each particle
  a=-180;
  ones2=matrix(1,1,num_ens); AHpt=as.matrix(AHpt,length(AHpt),1)
  AHpt=AHpt%*%ones2
  BT1=exp(a*AHpt+b)+ones1%*%So[6,];
  beta=BT1/(ones1%*%So[4,]);
  tcurrent=tm.ini;
  Sr_tmp=propagateParSIR(tcurrent+dt,tcurrent+tmstep,dt,S0=So[1,], I0=So[2,], N, D=So[4,], L=So[3,], beta, realdata=T)
  x[1,,1]=tail(Sr_tmp$S,1) 
  x[2,,1]=tail(Sr_tmp$I,1)
  x[3:6,,1]=So[3:6,];
  x[num_var,,1]=tail(Sr_tmp$newI,1);
  # y[,1]=tail(Sr_tmp$newI,1);
  
  ### This is for the regularization
  #  Getting the regularized component (random draw off a Epanechnikov kernel) is a bit complex.  
  # If n is the dimension of the state vector (could include the parameters as these are also adjusted), 
  # then the Epanechnikov kernel is K = (n+2)/(2C) * (1 - x^2) if |x|<1
  # and 0 otherwise.  C is the volume of a unit hypersphere in dimensions n.  
  # So for n=2 it is simply pi (i.e. pi*r^2, where r=1).  For n=6, it is pi^3/6.
  # here, for n=8 (even number), Vn(R)=Cn*R^n; Cn=pi^(n/2)/(n/2)!=pi^4/4!=pi^4/24.
  
  #  The optimal bandwidth for the kernel draw is
  #  h_opt=A/(N^(1/n+4)), where N is the number of particles
  #  and A = [(8/(C^(n+4)))*(2*sqrt(pi)^n)]^(1/n+4)
  # if use SEIR model:
  # C=pi^4/24;
  # A=((8/C^(8+4))*(2*sqrt(pi)^8))^(1/(8+4));
  # hopt=2*A/num_ens^(1/(8+4));
  # SIR model, 6 variables
  ## same famular as in Arulampalam et al. 2002
  # C=pi^3/6;
  # A=((8/C^(6+4))*(2*sqrt(pi)^6))^(1/(6+4)); # typo?
  # hopt=2*A/num_ens^(1/(6+4));  # typo?
  nn=num_var; # number of parameters
  C=pi^(nn/2)/gamma(nn/2+1);
  A=((8/C*(nn+4))*(2*sqrt(pi))^nn)^(-1/(nn+4)); # # change the exponent '(1/(6+4))' to (-1/(6+4))
  # hopt=A*num_ens^(-1/(6+4));  # checked right, from Density Estimation by Silverman
  hopt=2*A*num_ens^(-1/(nn+4));  # increase the band width
  xK=seq(-1,1,by=.001)
  dK=((nn+2)/2/C)*(1-xK^2); # density of the Epanechnickov Kernel
  dK=dK/sum(dK); # normolized
  # construct the cumulative density function
  cK=dK;
  for (i in 2:length(dK)){
    cK[i]=cK[i-1]+dK[i];
  }
  
  # Particle Filtering
  SD1=apply(x[,,1],1,sd);
  # print(SD1);
  wts0=rep(1,num_ens)/num_ens;
  wts=nwts=matrix(0,num_ens,ntrn);
  
  for (tt in 1:ntrn){
    # wts[,tt]=dpois(round(x[2,,tt],0),as.numeric(obs_i[tmstep*tt+1])); # use a Poisson dist. all wts[,tt-1]=1/num_ens, b/c of resampling
    # Poison dist. does not work!
    if (tt==1){
      wts[,tt]=wts0*dnorm(x=x[num_var,,tt],mean=obs_i[tt],sd=sqrt(obs_vars[tt]));
    } else{
      # ccumulative weights
      wts[,tt]=wts[,tt-1]*dnorm(x=x[num_var,,tt],mean=obs_i[tt],sd=sqrt(obs_vars[tt]));
    }
    
    # Normalizing the weights here
    # sometimes obs_i=0, if so, density=0 for all particles, can nornalize (divide by 0)
    if (sum(wts[,tt])==0 | any(is.na(wts[,tt]))){
      nwts[,tt]=1/num_ens; # assign equal weights if get no information from the likelihood
    } else {
      nwts[,tt]=wts[,tt]/sum(wts[,tt]);  # Normalizing here
    }
    # Y[tt]=sum(y[,tt]*nwts[,tt]); # weighted mean newI
    # Resampling with regulation
    # resample only when Neff<num_ens/2;
    neff=1/(sum(nwts[,tt]^2));
    if(neff<num_ens/2){
      # use R's sample function
      currx=x[,,tt];   #  Getting current state and parameters
      ind_new=sample(x=1:num_ens, size=num_ens, replace=T, prob=nwts[,tt])
      nwts[,tt]=1/num_ens; # reset the weights to equal after resampling.
      wts[,tt]=1/num_ens;
      currx=currx[,ind_new];
      if(regul==TRUE){
        ## add regularization noise
        SD=apply(currx,1,sd);
        if (length(unique(currx[2,]))<max(.01*num_ens,20) 
            #  |obs_i[tmstep*tt+1]>max(currx[2,]) 
            # |obs_i[tmstep*tt+1]<min(currx[2,])  # don't use, keep it simple!
        ){ 
          # | abs(mean(currx[2,])-obs_i[tmstep*tt+1])>obs_i[tmstep*tt+1]*.2
          SD=apply(x[,,max(tt-2,1)],1,sd);  # go back 2 steps before it degenerates.
          SD[2]=max(SD[2],obs_i[tt]/5); # in case ~ the peak
          # abs(median(currx[2,])-obs_i[tmstep*tt+1])
          # scale SD1 in case it happened in the end.
          # print(c(tt,length(unique(currx[2,]))));
        }
        # print(c(tt,SD));
        
        nze=matrix(0,num_var,num_ens)
        for (i in 1:num_var){
          nze[i,]=approx(cK,xK,runif(num_ens))$y;  # Regularize noise ?
        }
        
        # currx=currx[,ind_new]+hopt*SD*nze; # standard regularization
        # adjust the band width according to obs_i such at hopt.adj->1 for small obs_i, 
        # and hopt.adj >1 for big obs_i
        # hopt.adj=
        currx=currx+hopt*diag(SD,length(SD))%*%nze; # small jitter
      }
      
      ## check DA aphyiscality
      currx=Fn_checkDA(currx,bound.low=c(rep(0,2),theta_low,0),  # note: newI included
                       bound.up=c(rep(N,2),theta_up,N));
      
      # x[1:2,,tt]=round(currx[1:2,],0);
      x[,,tt]=currx;  
    }
    
    # integrate forward 1 step parallelly
    b=log(x[5,,tt]-x[6,,tt]); b=ones1%*%b; # expand b to a matrix with each col for each particle
    a=-180;
    BT1=exp(a*AHpt+b)+ones1%*%x[6,,tt];
    beta=BT1/(ones1%*%x[4,,tt]); 
    tcurrent = tm.ini+tmstep*tt;
    Sr_tmp=propagateParSIR(tcurrent+dt,tcurrent+tmstep,dt,x[1,,tt],x[2,,tt],N,D=x[4,,tt],L=x[3,,tt],beta,realdata=T)
    x[1,,tt+1]=tail(Sr_tmp$S,1);
    x[2,,tt+1]=tail(Sr_tmp$I,1);
    x[3:6,,tt+1]=x[3:6,,tt];
    x[num_var,,tt+1]=tail(Sr_tmp$newI,1);
    # y[,tt+1]=tail(Sr_tmp$newI,1);
  } # end for-loop
  #  Now getting the mean of the distribution at each time step
  xmean=xsd=matrix(0,num_var,ntrn)
  for (k in 1:num_var){
    for (t in 1:ntrn){
      xmean[k,t]=sum(x[k,,t]*nwts[,t]);  # after resampling, weights are equal for all particles
      xsd[k,t]=sqrt(sum(nwts[,t]*(x[k,,t]-xmean[k,t])^2));  # standard deviation
    }
  }
  #### Forecast
  ## resample according to the final nwts, so don't need to adjust for nwts later
  currx=x[,,ntrn];   #  Getting current state and parameters
  # use R's sample function
  # resample using the normalized weight
  # ind_new=sample(x=1:num_ens, size=num_ens, replace=T, prob=nwts[,ntrn]);
  # currx=currx[,ind_new];
  b=log(currx[5,]-currx[6,]); b=ones1%*%b; # expand b to a matrix with each col for each particle
  a=-180;
  BT1=exp(a*AHpt+b)+ones1%*%currx[6,];
  beta=BT1/(ones1%*%currx[4,]); 
  tcurrent = tm.ini+tmstep*ntrn;
  Sr_tmp=propagateParSIR(tcurrent+dt,tcurrent+tmstep*nfc,dt,currx[1,],currx[2,],N,
                         D=currx[4,],L=currx[3,],beta,realdata=T)
  
  fcast[1,,]=t(Sr_tmp$S[tmstep*(1:nfc)+1,]); # num_ens, time 
  fcast[2,,]=t(Sr_tmp$I[tmstep*(1:nfc)+1,]);
  fcast[3,,]=t(Sr_tmp$newI[tmstep*(1:nfc)+1,]-Sr_tmp$newI[tmstep*(0:(nfc-1))+1,]);
  
  fcast_mean=fcast_sd=matrix(0,3,nfc)
  for (tt in 1:nfc){
    # fcast_mean[,tt]=apply(fcast[,,tt],1,mean, na.rm=T);
    # fcast_sd[,tt]=apply(fcast[,,tt],1,sd, na.rm=T);
    for(ii in 1:3){
      fcast_mean[ii,tt]=sum(fcast[ii,,tt]*nwts[,ntrn]);
      fcast_sd[ii,tt]=sqrt(sum(nwts[,ntrn]*(fcast[ii,,tt]-fcast_mean[ii,tt])^2));
    }
  }
  # Getting the mean predictions
  #  First normalize the weights
  # nwts=wts/(ones(num_ens,1)*sum(wts));
  
  ## output data
  pkwks=rep(0,num_ens);
  obs_pkwk=which.max(obs_i);
  for (i in 1:num_ens){
    pkwks[i]=which.max(c(x[num_var,i,1:ntrn],fcast[3,i,1:nfc]));
  }
  pkwk_mode=MODE(pkwks)[1];
  pkwk_mode_perc=MODE(pkwks)[2]/num_ens;
  leadpkwk_mode=pkwk_mode-ntrn;
  ## check the mean as well
  Y=xmean[num_var,]; # newI
  y.all=c(Y,fcast_mean[3,]);  # newI
  pkwk_mean=which.max(y.all);
  delta_pkwk_mean=pkwk_mean-obs_pkwk;
  leadpkwk_mean=pkwk_mean-ntrn;
  ## weighted var
  if(pkwk_mean<=ntrn){ # peak before the forecast
    peak_nwts=nwts[,pkwk_mean];
  } else {
    peak_nwts=nwts[,ntrn]; # b/c resmp before fcast
  }
  pkwk.ens.mean=sum(pkwks*peak_nwts);
  pkwk_var=sum(peak_nwts*(pkwks-pkwk.ens.mean)^2) # note the weighted variance is biased
  # Note: the mean and nwts at the time of pkwk_mean according to the mean trajectory are not matched
  # the mean of ens of pkwks is sum(peak_nwts*pkwks)
  # pkwk_var=var(pkwks-obs_pkwk); # Wrong if forecast is made after the predicted peak!!
  
  ## output data
  ## output the prediction from the last iteration
  tstep=seq(tm.ini+tmstep,num_times*tmstep+tm.ini,by=tmstep); # time
  fc_start=ntrn+wk_start-1;
  
  out1=cbind(rep(fc_start,ntrn),tstep[1:ntrn],t(xmean)); # NOTE: Y here only include the training
  
  # metrics for comparison
  corr=cor(y.all,head(obs_i,nsn)); rms=sqrt(mean((y.all-head(obs_i,nsn))^2));
  delta_sum_newI=sum(y.all)-sum(head(obs_i,nsn));
  delta_pkwk=pkwk_mode-obs_pkwk;
  
  out2=cbind(rep(fc_start,ntrn),tstep[1:ntrn],t(xsd));
  out3=cbind(rep(fc_start,nfc),tstep[(ntrn+1):nsn],t(fcast_mean));
  out4=cbind(rep(fc_start,nfc),tstep[(ntrn+1):nsn],t(fcast_sd));
  out5=t(c(fc_start,rms,corr,delta_sum_newI,delta_pkwk,leadpkwk_mode,
           pkwk_var,pkwk_mode_perc,pkwk_mean,leadpkwk_mean,delta_pkwk_mean));
  
  colnames(out1)=colnames(out2)=c("fc_start","time","S","I","L","D","R0max","R0min","newI");
  # colnames(out2)=c("fc_start","time","S","I","L","D","R0max","R0min");
  colnames(out3)=colnames(out4)=c("fc_start","time","S","I","newI");
  colnames(out5)=c("fc_start","rms","corr","delta_sum_newI","delta_pkwk","leadpkwk_mode",
                   "pkwk_var","pkwk_mode_perc","mn_pkwk","mn_leadpk","delta_mn_pkwk");
  
  if (metricsonly==F){  
    out=list(train=out1,trainsd=out2,fcast=out3,fcsd=out4,metrics=out5);
  } else {
    out=list(metrics=out5);
  }
}
