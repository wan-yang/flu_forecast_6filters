## Function to run Maximum Likelihood Extimation with Particle Filtering (MIF), Then forecast
## Algorithm based on Ionides et al. (2006) PNAS 103: 18438-18443.
## Procedure 1. (MIF), Supplementary S3 MIF via sequential Monte Carlo
## 1. Select starting values theta[1], a discount factor 0<alpha<1, (alpha=0.95 in the paper),
##     an initial variance multiplier c_sq (=20 in the paper), and the number of iterations N.
## 2. For n in 1,..., N
##   (i) Set sig[n]=alpha^(n-1), for t=1,...,T, evaluate theta[t,n]=theta_t(theta[n],sig[n]) 
##       and V[t,n]=V_t(theta[n],sig[n]);
##       use particle filtering to estimate theta[t,n] and V[t,n]
##   (ii) Set theta[n+1]=theta[n]+V[1,n]*sum(t=1:T, V[t,n]^(-1)* (theta[t,n]-theta[t-1,n])),
##       where theta[0,n]=theta[n]
## 3. Take theta[N+1] to be a maximum likelihood estimate of the parameter theta for the fixed parameter model.
##     E[theta_t|theta_t-1]=theta_t-1, Var(theta_t|theta_t-1)=sig^2*SIG (eqn.1)
##     E[theta_0]=theta, Var(theta_0)=sig^2*c_sq*SIG
## SIG: typically a diagonal matrix giving the respective scales of each component of theta
## theta_t=theta_t(theta,sig)=E[theta_t|y1:t]  (eqn.2)
## V_t=V_t(theta,sig)=Var(theta_t|y1:t-1)
##  More information in 'Comparison of filtering methods for the modeling and retrospective forecasting of influenza epidemics' (PLoS Compute Biol)
##  by Wan Yang, Alicia Karspeck, and Jeffrey Shaman, 2014

MIF_rFC<-function(num_ens, tmstep, param.bound, obs_i=obs_i, ntrn=1,
                   obs_vars,tm.ini=273, tm.range=273:500,
                   IN=30,regul=T){
  ## MIF to run retrospecitve forecast of influenza
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
  ##        IN: number of iterations
  ## Outputs:
  ##        train: modeled time series of the variables/parameters (posteriors) for the training period
  ##        trainsd: standard deviation of the variables/parameters during the training
  ##        fcast: forecasted time series of the variables (i.e., S, I, and newI)
  ##        fcsd: standard deviation of the variables of the forecast
  ##        metrics: root mean squared error (rms), correlation with the observations, ect.
  ##        trainparams: parameter estimates in the iteration
  ##        Note: if only wish to save the metrics, set metricsonly=TRUE
  
  library("truncnorm"); library("tgp"); # for lhs
  library("MASS"); # for multivariate normal distribution
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
  nsn=length(obs_i);
  num_var=2+4+1; # include S, I, L, D, R0max, R0min, newI
  nfc=nsn-ntrn; # number of weeks for forecasting;
  fcast=array(0,c(3,num_ens,nfc)); # fcast: S, I, newI;
  
  theta_low=param.bound[,1];
  theta_up=param.bound[,2];
  
  # 1. set initial values: use x for variable matrix (including state variables and parameters)
  theta=matrix(0,4,IN+1); # matrix storing parameters at each iteration
  xmean=matrix(0,num_var,ntrn); # matirx storing state variables and parameters at each time step
  sig=rep(0,IN);
  alp=0.95; # alpha
  c_sq=20;
  
  SIG=(theta_up-theta_low)^2/4/num_times; # variance of theta
  # Initial theta:
  rnd=ceiling(1e4*runif(1))
  
  theta[,1]=c(params[rnd,1]*365,params[rnd,2]*365,params[rnd,3],params[rnd,4])
  # log likelihood
  l=rep(0,IN)
  
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
  ## same famular as in Arulampalam et al. 2002
  nn=2; # number of parameters
  C=pi;
  # A=((8/C^(6+4))*(2*sqrt(pi)^6))^(1/(6+4)); # typo?
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
  flag=NULL; # record err
  # run N iteration
  for (n in 1:IN){
    # print(c("interation #:",n));
    sig[n]=alp^(n-1);
    # initialize particles
    So=matrix(0,6,num_ens);
    x=array(0,c(num_var,num_ens,ntrn+1)); # include newI
    # y=matrix(0,num_ens,ntrn+1);  # newI
    # Y=rep(0,ntrn); # weighted mean newI
    # Like=rep(0,num_times); # likelihood
    rnd=cbind(ceiling(27*1e4*runif(num_ens)), ceiling(27*1e4*runif(num_ens)));
    So[1,]=susceps[rnd[,1]]/5;  # 100,000 population
    So[2,]=infects[rnd[,2]]/5;
    
    # generate new particles 
    # generate new particles using multivariate normal distribution
    if (n==1){
      So[3:6,]=t(mvrnorm(num_ens,mu=theta[,n],Sigma=diag(sig[n]^2*c_sq*SIG,4))) 
    } else{
      So[3:6,]=t(mvrnorm(num_ens,mu=theta[,n],Sigma=diag(sig[n]^2*SIG,4)))
    }
    # correct lower/upper bounds of the proposals
    ug=min(So[3:6,])
    if (ug<=0){
      for (ii in 3:6){
        So[ii,]=pmax(So[ii,],theta_low[ii-2])
      }
    }
    So[5,]=pmax(So[5,],So[6,]+0.01); # make sure Rmx > Rmn
    
    # integrate forward 1 step parallelly
    ## Calculate the reproductive number at time t BT1 and the transmission rate 
    ## according to Shaman J, Karspeck A (2012) Proc Natl Acad Sci USA 109: 20425-20430.
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
    x[2,,1]=tail(Sr_tmp$I,1);
    x[3:6,,1]=So[3:6,];
    x[num_var,,1]=tail(Sr_tmp$newI,1);
    # y[,1]=tail(Sr_tmp$newI,1);
    
    wts0=rep(1,num_ens)/num_ens;
    wts=nwts=matrix(0,num_ens,num_times);
    # resmp=matrix(0,num_ens,num_times);
    # tm.resmp=NULL;
    cumlike=NULL;
    # SD1=apply(x[,,1],1,sd); # set SD fixed, all the same as the initial
    SD1=apply(x[1:2,,1],1,sd); # add noise to S and I only
    # print(c(1,SD1));
    # Particle Filtering
    for (tt in 1:ntrn){
      # wts[,tt]=dpois(round(x[2,,tt],0),obs_i[tmstep*tt+1]); # use a Poisson dist. all wts[,tt-1]=1/num_ens, b/c of resampling
      # Poison dist. does not work!
      if (tt==1){
        wts[,tt]=wts0*dnorm(x=x[num_var,,tt],mean=obs_i[tt],sd=sqrt(obs_vars[tt]));
      } else{
        # ccumulative weights
        wts[,tt]=wts[,tt-1]*dnorm(x=x[num_var,,tt],mean=obs_i[tt],sd=sqrt(obs_vars[tt]));
      }
      
      # Like[tt]=mean(wts[,tt]); # cumulative likelihood
      # Normalizing the weights here
      # sometimes obs_i=0, if so, density=0 for all particles, can nornalize (divide by 0)
      if (sum(wts[,tt])==0 | any(is.na(wts[,tt]))){
        nwts[,tt]=1/num_ens; # assign equal weights if get no information from the likelihood
        cumlike=append(cumlike,mean(wts[,tt-1]));
        # record err: iteration, tm step
        flag=rbind(flag,c("degeneracy",n,tt));
        # once degenerate, no point to resample again (only 1 particle exists)
      } else {
        nwts[,tt]=wts[,tt]/sum(wts[,tt]);  # Normalizing here
      }
      # Y[tt]=sum(y[,tt]*nwts[,tt]); # weighted mean newI
      
      # Resampling with or without regulation 
      neff=1/(sum(nwts[,tt]^2));
      #  Resample if Neff is too small
      if (neff<num_ens/2){
        # print(c("resmp at time step:",tt));
        currx=x[,,tt];   #  Getting current state and parameters
        # use R's sample function
        # resample using the normalized weight
        ind_new=sample(x=1:num_ens, size=num_ens, replace=T, prob=nwts[,tt]);
        # record cummulative likelihood upto resampling
        # should record before resetting the weights!
        cumlike=append(cumlike,mean(wts[,tt]));
        nwts[,tt]=1/num_ens; # reset the weights to equal after resampling.
        wts[,tt]=1/num_ens;
        # tm.resmp=append(tm.resmp,tt); # record time of resampling
        currx=currx[,ind_new];
        if (regul==TRUE){
          ## add regularization noise
          # SD=apply(currx,1,sd);
          SD=apply(currx[1:2,],1,sd);
          SD=pmax(SD,c(mean(x[1,,tt])/100,obs_i[tt]/20)); # in case ~ the peak
          ## check how many distinct particles exist, if too few, probably degenerate, widden SD.
          ## not to check SD, b/c SD could still be huge 
          ## when there are only a few very different particles (e.g., bimodal) 
          if (length(unique(currx[2,]))<max(.01*num_ens,20)){
            # SD=apply(x[,,max(tt-2,1)],1,sd);  # go back 2 steps before it degenerates.
            # SD[2]=max(SD[2],obs_i[tmstep*tt+1]/5); # in case ~ the peak
            SD=apply(x[1:2,,max(tt-2,1)],1,sd);  # go back 2 steps before it degenerates.
            SD=pmax(SD,c(mean(x[1,,tt]),obs_i[tt])/10); # in case ~ the peak
            # print(c(tt,length(unique(currx[2,]))));
          }
          # print(c(tt,SD));
          # draw epsilon ~ K from the Epanechnickov Kernel
          # nze=matrix(0,6,num_ens)
          nze=matrix(0,2,num_ens)
          for (i in 1:2){
            nze[i,]=approx(cK,xK,runif(num_ens))$y;  # Regularize noise ?
          }
          
          # currx=currx[,ind_new]+hopt*SD*nze; # standard regularization
          currx[1:2,]=currx[1:2,]+hopt*diag(SD,length(SD))%*%nze; # small jitter
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
    } # end for-loop (all time steps)
    
    # get theta[t,n] and V[t,n], that is the mean and variance at each time step
    xmean=xsd=matrix(0,num_var,ntrn)
    for (k in 1:num_var){
      for (t in 1:ntrn){
        xmean[k,t]=sum(x[k,,t]*nwts[,t]);  # after resampling, weights are equal for all particles
        xsd[k,t]=sqrt(sum(nwts[,t]*(x[k,,t]-xmean[k,t])^2));  # standard deviation
      }
    }
    # estimate theta[n+1]
    # theta[n+1]=theta[n]+V[1,n]*sum(t=1:T, V[t,n]^(-1)* (theta[t,n]-theta[t-1,n])),
    
    # update theta using Eqn.17 in the Supplementary Material
    # should use this for the first 5 iterations and then switch to the eq. in the main text
    # but more stable if use Eqn.17 for the whole process 
    # Note: this does not guarantee maximum likelihood.
    # theta[n+1]=1/T*sum(theta[t,n])
    theta[,n+1]=rowMeans(xmean[3:6,]);
    # l[n]=sum(log(Like)); # total likelihood (marginal likelihood, i.e., p(y1:T))
    if (is.null(cumlike)){
      l[n]=log(mean(wts[,tt])); # no resample done
    } else{
      l[n]=sum(log(cumlike),log(mean(wts[,tt])))
    }
  } # end iteration
  
  #### Forecast
  ## resample according to the final nwts, so don't need to adjust for nwts later
  currx=x[,,ntrn];   #  Getting current state and parameters
  # use R's sample function
  # resample using the normalized weight
  ind_new=sample(x=1:num_ens, size=num_ens, replace=T, prob=nwts[,ntrn]);
  currx=currx[,ind_new];
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
    fcast_mean[,tt]=apply(fcast[,,tt],1,mean, na.rm=T);
    fcast_sd[,tt]=apply(fcast[,,tt],1,sd, na.rm=T);
  }
  
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
    peak_nwts=rep(1/num_ens,num_ens); # b/c resmp before fcast
  }
  pkwk.ens.mean=sum(pkwks*peak_nwts);
  leadpkwk_var=sum(peak_nwts*(pkwks-pkwk.ens.mean)^2) # note the weighted variance is biased
  # Note: the mean and nwts at the time of pkwk_mean according to the mean trajectory are not matched
  # the mean of ens of pkwks is sum(peak_nwts*pkwks)
  # leadpkwk_var=var(pkwks-obs_pkwk); # Wrong if forecast is made after the predicted peak!!
  
  ## output data
  ## output the prediction from the last iteration
  tstep=seq(tm.ini+tmstep,num_times*tmstep+tm.ini,by=tmstep); # time
  fc_start=ntrn+wk_start-1;
  
  out1=cbind(rep(fc_start,ntrn),tstep[1:ntrn],t(xmean)); # NOTE: Y here only include the training
  # metrics for comparison
  # Y=c(Y,fcast_mean[3,]);  # newI
  corr=cor(y.all,head(obs_i,nsn)); rms=sqrt(mean((y.all-head(obs_i,nsn))^2));
  delta_sum_newI=sum(y.all)-sum(head(obs_i,nsn));
  delta_pkwk=pkwk_mode-obs_pkwk;
  
  out2=cbind(rep(fc_start,ntrn),tstep[1:ntrn],t(xsd));
  out3=cbind(rep(fc_start,nfc),tstep[(ntrn+1):nsn],t(fcast_mean));
  out4=cbind(rep(fc_start,nfc),tstep[(ntrn+1):nsn],t(fcast_sd));
  out5=t(c(fc_start,rms,corr,delta_sum_newI,delta_pkwk,leadpkwk_mode,
           leadpkwk_var,pkwk_mode_perc,pkwk_mean,leadpkwk_mean,delta_pkwk_mean));
  out6=cbind(rep(fc_start,IN+1),0:IN,t(theta),c(NA,l));
  
  colnames(out1)=c("fc_start","time","S","I","L","D","R0max","R0min","newI");
  colnames(out2)=c("fc_start","time","S","I","L","D","R0max","R0min","newI");
  colnames(out3)=colnames(out4)=c("fc_start","time","S","I","newI");
  colnames(out5)=c("fc_start","rms","corr","delta_sum_newI","delta_pkwk","leadpkwk_mode",
                   "leakpkwk_var","pkwk_mode_perc","mn_pkwk","mn_leadpk","delta_mn_pkwk");
  colnames(out6)=c("fc_start","iteration","L","D","R0max","R0min","log_like");
  
  if (metricsonly==F){  
    out=list(train=out1,trainsd=out2,fcast=out3,fcsd=out4,metrics=out5,trainparams=out6);
  } else {
    out=list(metrics=out5);
  }
}
