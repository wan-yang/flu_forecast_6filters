## Function to run particle MCMC, then forecast
## Based on Andrieu C, Doucet A, Holenstein R (2010) J. R. Statist. Soc. B 72: 269-342.
##          and Rasmussen et al. 2011 PLoS Comput Biol 7: e1002136.
##          and Shaman J, Karspeck A (2012) Proc Natl Acad Sci USA 109: 20425-20430.
##  More information in 'Comparison of filtering methods for the modeling and retrospective forecasting of influenza epidemics' (PLoS Compute Biol)
##  by Wan Yang, Alicia Karspeck, and Jeffrey Shaman, 2014

pMCMC_rFC<-function(num_ens, tmstep, param.bound, obs_i=obs_i,ntrn=1,
                      obs_vars=obs_vars,tm.ini=273, tm.range=273:500,
                      M=80,regul=TRUE){
  ## pMCMC to run retrospecitve forecast of influenza
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
  ##        M: number of sets of parameter to test in the first round (coarse resolution)
  ## Outputs:
  ##        train: modeled time series of the variables/parameters (posteriors) for the training period
  ##        trainsd: standard deviation of the variables/parameters during the training
  ##        fcast: forecasted time series of the variables (i.e., S, I, and newI)
  ##        fcsd: standard deviation of the variables of the forecast
  ##        metrics: root mean squared error (rms), correlation with the observations, ect.
  ##        trainparams: parameter estimates in each MCMC step
  ##        Note: if only wish to save the metrics, set metricsonly=TRUE
  # library to use multi-norm sampling and latin hypercube sampling
  library(MASS); library(sampling); library(tgp);
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
  
  num_times=length(obs_i);  num_var=3; # S, I, newI, for the PF
  nsn=length(obs_i);
  nfc=nsn-ntrn; # number of weeks for forecasting;
  
  theta_low=param.bound[,1];
  theta_up=param.bound[,2];
  D_mean=3; L_mean=4*365; Rmn_mean=0.8; Rmx_mean=3;
  theta_mean=c(L_mean, D_mean, Rmx_mean, Rmn_mean );
  tuning=1/10
  L_len=tuning*(theta_up[1]-theta_low[1]); D_len=tuning*(theta_up[2]-theta_low[2]);   
  Rmx_len=tuning*(theta_up[3]-theta_low[3]); 
  Rmn_len=tuning*(theta_up[4]-theta_low[4]);  
  # theta_len=c(L_len, D_len, Rmx_len, Rmn_len);
  ## AH for the whole function, not restrict to PF_regul
  beta.range=tm.range[1]:(tail(tm.range,1)+tmstep);
  AHpt=AH[beta.range,]; # only use AH during the time of interest, to minimize space needed
  AHpt=as.matrix(AHpt,length(AHpt),1)
  AHpt=AHpt%*%matrix(1,1,num_ens);
  ## Algorithm 2: the particle filter component of particle MCMC
  ## In the following algorithm, theta is the vector containing all model parameters,
  ## X1:T are the latent variables and Z1:T are the observed data. t = 1 to T are the observation times. 
  ## j= 1 to N are the particle indices. For example, xt^j  represents the state of particle j at time t.  
  ## The notation At^j  is used to track the ancestry of particles backward in time, 
  ## such that At-1^j  represents the parent index of particle j at time t.
  ## run PF with regularization
  # function to run particle filter
  PF_regul<-function(theta,num_ens){
    ### This is for the regularization
    #  Getting the regularized component (random draw off a Epanechnikov kernel) is a bit complex.  
    # If n is the dimension of the state vector (could include the parameters as these are also adjusted), 
    # then the Epanechnikov kernel is K = (n+2)/(2C) * (1 - x^2) if |x|<1
    # and 0 otherwise.  C is the volume of a unit hypersphere in dimensions n.  
    # So for n=2 it is simply pi (i.e. pi*r^2, where r=1).  For n=6, it is pi^3/6.
    # here, for n=8 (even number), Vn(R)=Cn*R^n; Cn=pi^(n/2)/(n/2)!=pi^4/4!=pi^4/24.
    # for n=7 (odd number), Vn=16/105*pi^3
    
    #  The optimal bandwidth for the kernel draw is
    #  h_opt=A/(N^(1/n+4)), where N is the number of particles
    #  and A = [(8/(C^(n+4)))*(2*sqrt(pi)^n)]^(1/n+4)
    # if use SEIR model:
    # C=pi^4/24;
    # A=((8/C^(8+4))*(2*sqrt(pi)^8))^(1/(8+4));
    # hopt=2*A/num_ens^(1/(8+4));
    # SIR model, 6 variables
    ## same famular as in Arulampalam et al. 2002
    # C=pi;
    # C=pi^3/6;
    # C=16/105*pi^3; # for n=7
    # A=((8/C^(6+4))*(2*sqrt(pi)^6))^(1/(6+4)); # typo?
    # hopt=2*A/num_ens^(1/(6+4));  # typo?
    
    ### UPDATED HERE
    ### NOTE: FIX THE NUMBER OF VARIABLES, WHICH COULD CHANGE DUE TO MODEL FORMULATION
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
    
    # step 1: initialize particle filter at time t=1
    # (a) set x1^j to initial values for all particles.
    # rnd=cbind(ceiling(1e4*runif(num_ens)), ceiling(27*1e4*runif(num_ens)), ceiling(27*1e4*runif(num_ens)));
    # rnd=cbind(ceiling(27*1e4*runif(num_ens)), ceiling(27*1e4*runif(num_ens)));
    # num_times=floor((tail(tm.range,1)-tm.range[1])/tmstep)
    So=matrix(0,2,num_ens); # initials
    
    ### UPDATED HERE
    ### NOTE: INCLUDE THE OBSERVED VARIABLE (EARLY STORED SEPEARTELY AS LITTLE Y)
    ### IN THE COMBINED LITTLE X ARRAY, SO IT WILL BE UPDATED SIMULTANEOUSLY AFTER RESAMPLING
    x=array(0,c(num_var,num_ens,ntrn+1));
    # y=matrix(0,num_ens,ntrn); # newI
    # Y=rep(0,ntrn); # weighted mean newI
    w=W=matrix(0,num_ens,ntrn) # w: unnormalized weights; W: normalized weights
    cumlike=NULL;  # cummulative likelihood
    # p_hat=rep(0,num_times+1)   # estimated marginal likelihood
    flag=NULL; # record err
    A=matrix(0,num_ens,ntrn)  # index
    rnd=cbind(ceiling(27*1e4*runif(num_ens)), ceiling(27*1e4*runif(num_ens)));
    So[1,]=susceps[rnd[,1]]/5;  # per 100,000 population
    So[2,]=infects[rnd[,2]]/5; 
    L=theta[1]; D=theta[2]; # alpha=1/theta[5];
    R0max=theta[3]; R0min=theta[4];
    # beta.range=tm.range[1]:(tail(tm.range,1)+tmstep);
    # make sure to include enought beta, b/c the prior integrate 1 time step more than needed
    # AHpt=AH[beta.range,]; # only use AH during the time of interest, to minimize space needed
    ones1=matrix(1,dim(AHpt)[1],num_ens); # matrix of 1 to make the dimension match
    b=log(R0max-R0min);  # expand b to a matrix with each col for each particle
    a=-180;
    # ones2=matrix(1,1,num_ens); 
    # AHpt=as.matrix(AHpt,length(AHpt),1)
    # AHpt=AHpt%*%ones2
    BT1=exp(a*AHpt+b*ones1)+R0min*ones1;
    beta=BT1/D; # for all particles
    tcurrent=tm.ini;
    Sr_tmp=propagateParSIR(tcurrent+dt,tcurrent+tmstep,dt,S0=So[1,], I0=So[2,], N, D, L, beta, realdata=T)
    x[1,,1]=tail(Sr_tmp$S,1) 
    x[2,,1]=tail(Sr_tmp$I,1);
    
    ### UPDATED HERE
    ### NOTE: INCLUDE THE OBSERVED VARIABLE (EARLY STORED SEPEARTELY AS LITTLE Y)
    ### IN THE COMBINED LITTLE X ARRAY, SO IT WILL BE UPDATED SIMULTANEOUSLY AFTER RESAMPLING
    x[num_var,,1]=tail(Sr_tmp$newI,1);  # at the end of 1st week
    # y[,1]=tail(Sr_tmp$newI,1);  # at the end of 1st week
    
    # x[3,,1]=infects[rnd[,3]]; # new I;
    # x[1,,1]=rnorm(num_ens,S0,S0/50);
    # x[1,x[1,,1]<0,1]=mean(x[1,,1])
    # x[2,,1]=rpois(num_ens,10);
    # x[3,,1]=rpois(num_ens,1); 
    SD1=apply(x[,,1],1,sd); # set SD fixed, all the same as the initial
    # print(SD1);
    # (b) assign unormalized particle weights: w1^j=p(z1|x1^j,theta)
    # w[,1]=dnorm(y[,1],obs_i[1],sd=sqrt(obs_vars[1]))  # this is actually p(xt|zt), but normal dist. is symetric, so p(zt|xt)=p(xt|zt)
    
    ### UPDATED HERE
    ### NOTE: INCLUDE THE OBSERVED VARIABLE (EARLY STORED SEPEARTELY AS LITTLE Y)
    ### IN THE COMBINED LITTLE X ARRAY, SO IT WILL BE UPDATED SIMULTANEOUSLY AFTER RESAMPLING
    ### x[num_var,,1] IS NOW THE OBSERVED VARIABLE
    w[,1]=dnorm(x[num_var,,1],obs_i[1],sd=sqrt(obs_vars[1]))
    # (c) normalize the particle weights
    if (sum(w[,1])==0){
      W[,1]=1/num_ens; # assign equal weights if get no information from the likelihood
    } else {
      W[,1]=w[,1]/sum(w[,1]);  # Normalizing here
    }
    # Y[1]=sum(y[,1]*W[,1]);
    # p_hat[1]=mean(w[,1])
    neff=1/(sum(W[,1]^2));
    # A[,1]=1:num_ens
    ## Step 2: Run particle filter
    for (t in 2:(ntrn)){
      if (neff<num_ens/2){
        #(a) resample particles by sampling parent particle indices k according to their weights, such that
        # p(At-1^j=k)=Wt-1^k
        # A[,t]=sample(x=A[,t-1],size=num_ens, replace=T, prob=W[,t-1]); # every round, indices are reset to 1:num_ens
        A[,t-1]=sample(x=1:num_ens,size=num_ens, replace=T, prob=W[,t-1]);
        # (b) set xt-1^j=xt-1^(At-1^j) for all particles, with regularization
        currx=x[,A[,t-1],t-1]; # resampled parents
        # Y[t-1]=sum(y[A[,t-1],t-1]*W[A[,t-1],t-1]); # recalculate newI
        # NOTE: DON'T DO IT!
        # record cummulative likelihood upto resampling
        # should record before resetting the weights!
        cumlike=append(cumlike,mean(w[,t-1]));
        W[,t-1]=1/num_ens; # reset the weights to equal after resampling.
        w[,t-1]=1/num_ens;
        if (regul==TRUE){
          ## add regularization noise
          SD=apply(currx,1,sd);
          SD=pmax(SD,c(mean(currx[1,])/100,mean(currx[2,])/20,mean(currx[3,])/20))
          ## if allow # suscept to change freely, 
          ## the filter is prone to settle at high # suscept combined with low R0max
          if (length(unique(currx[2,]))<max(.01*num_ens,20)){
            SD=apply(x[,,max(t-2,1)],1,sd);  # go back 2 steps before it degenerates.
            SD[2]=max(SD[2],mean(currx[2,])/5); # in case ~ the peak
          }
          # print(c(t, SD));
          # draw epsilon ~ K from the Epanechnickov Kernel
          nze=matrix(0,num_var,num_ens); # x only contains 2 state parameters, i.e., S and I
          for (i in 1:num_var){
            nze[i,]=approx(cK,xK,runif(num_ens))$y;  # Regularize noise ?
          }
          currx=currx+hopt*diag(SD,length(SD))%*%nze; # small jitter
          ## check DA aphyiscality
          currx=Fn_checkDA(currx,bound.low=rep(0,3),  # note: newI included
                           bound.up=rep(N,3));
          
          # x[,,t-1]=round(currx[1:2,],0);  # Auxiliary PF
        } # end of regularization
        x[,,t-1]=currx;  # new parents
      } # end of resmp;
      
      # (c) propagate particles by simulating from the process model: 
      # xt^j ~ p(xt^j|xt-1^j,theta) to the next obs. time t
      tcurrent = tm.ini+tmstep*(t-1)  # go back one time step 
      Sr_tmp=propagateParSIR(tcurrent+dt,tcurrent+tmstep,dt,x[1,,t-1],x[2,,t-1],N,D,L,beta,realdata=T)
      # (d) set x1:t^j =(x1:t-1^j, xt^j) for all particles
      x[1,,t]=tail(Sr_tmp$S,1);
      x[2,,t]=tail(Sr_tmp$I,1);
      ### UPDATED HERE
      ### NOTE: INCLUDE THE OBSERVED VARIABLE (EARLY STORED SEPEARTELY AS LITTLE Y)
      ### IN THE COMBINED LITTLE X ARRAY, SO IT WILL BE UPDATED SIMULTANEOUSLY AFTER RESAMPLING
      ### x[num_var=3,,t] IS NOW THE OBSERVED VARIABLE
      x[3,,t]=tail(Sr_tmp$newI,1);
      # y[,t]=tail(Sr_tmp$newI,1);
      
      # (e) assign incremental, unnormalized particle weights
      # cummulative weights
      # w[,t]=w[,t-1]*dnorm(y[,t],obs_i[t],sd=sqrt(obs_vars[t]))
      
      ### UPDATED HERE
      ### NOTE: INCLUDE THE OBSERVED VARIABLE (EARLY STORED SEPEARTELY AS LITTLE Y)
      ### IN THE COMBINED LITTLE X ARRAY, SO IT WILL BE UPDATED SIMULTANEOUSLY AFTER RESAMPLING
      ### x[num_var=3,,t] IS NOW THE OBSERVED VARIABLE
      w[,t]=w[,t-1]*dnorm(x[3,,t],obs_i[t],sd=sqrt(obs_vars[t]))
      # (f) normalize the particle weights
      if (sum(w[,t])==0| any(is.na(w[,t]))){
        cumlike=append(cumlike,mean(w[,t])); # if degenerate, affect likelihood.
        W[,t]=1/num_ens; # assign equal weights if get no information from the likelihood
        w[,t]=1/num_ens;
        # record err: iteration, tm step
        flag=rbind(flag,c("degeneracy",t));
        # print(paste("Step",t," degenerate!"));
      } else {
        W[,t]=w[,t]/sum(w[,t]);  # Normalizing here
      }
      # Y[t]=sum(y[,t]*W[,t]); # weighted newI
      neff=1/(sum(W[,t]^2));
    } # end of times
    
    # resample and adjust the parents the last time
    # A[,ntrn]=sample(x=1:num_ens,size=num_ens, replace=T, prob=W[,ntrn]);
    # x[,,ntrn]=x[,A[,ntrn],ntrn]; # resampled parents
    # W[,ntrn]=1/num_ens; # reset normalized weights after resampling;
    # w[,ntrn]=w[A[,ntrn],ntrn]; # match the un-normalized weights
    
    ## Step 3: Estimate marginal likelihood: p_hat(z1:T|theta)=prod_1:T(p_hat(zt|z1:t-1,theta))
    ## where p_hat(zt|z1:t-1,theta)=1/N*sum(wt^j)
    if (is.null(cumlike)){
      p_hat_log=log(mean(w[,t])); # no resampling is done
    } else{
      p_hat_log=sum(log(cumlike),log(mean(w[,t])))
    }
    # p_hat_log=sum(log(p_hat))
    
    ## Step 4: sample x1:T* from p_hat(x1:T|theta, z1:T) by 
    ## tracing the lineage of one particle trajectory backwards through time
    
    ## it is probably not a good idea to sample at the end 
    ## and then use the samples and their ancestors to represent the trajectory
    ## becuase information at the end of the epi is not very accurate
    ## try just use the mean at each time step:
    x_star=x_sd=matrix(0,num_var,ntrn)
    for (k in 1:num_var){
      for (t in 1:ntrn){
        x_star[k,t]=sum(x[k,,t]*W[,t]);
        # x_sd[k,t]=sd(x[k,,t]*W[,t]);  # need to resample every time step to get the weights equal for all particles!!
        x_sd[k,t]=sqrt(sum(W[,t]*(x[k,,t]-x_star[k,t])^2)); # biased weighted sd
      }
    }
    
    ## UPDATED HERE
    ### NOTE: PASS THE WEITHS AS WELL: nwts
    ## also pass the particles (S&I) from the final training step: xtrn
    list(x_star=x_star,xtrn=x[,,ntrn],y=x[num_var,,1:ntrn],nwts=W,p_hat_log=p_hat_log,x_sd=x_sd,err=flag)  
  } # end of PF
  
  ## Algorithm 1: the MCMC component of particle MCMC
  ## In the following algorithm, theta is the vector containing all model parameters,   
  ## x1:T are the latent variables and z1:T are the observed data. m indexes the MCMC iterations from 1 to M.
  # step 1: initialize MCMC
  # (a) set m=0
  # (b) set theta(0) arbitrorily
  # (c) rum particle filter (Alg. 2) to sample x1:T(0) from p_hat(x1:T|theta(0),z1:T) 
  # and obtain the marginal likelihood estimate p_hat(z1:T|theta(0)
  
  MM=round(M*1.25,0);  # MM: M plus addtional 'zoom in' steps
  
  theta=matrix(0,4,MM)
  X=array(0,c(3,ntrn,MM)) # x1:T
  y=matrix(0,num_ens,ntrn);
  nwts=matrix(0,num_ens,ntrn); # normalized weights
  xtrn=matrix(0,2,num_ens); # particles (S&I) from the final training step
  xsd=matrix(0,2,ntrn); # standard deviation, only save the best run
  p_hat=rep(0,MM)  # p_hat(z1:T|theta)
  paramsMC=lhs(M,cbind(theta_low,theta_up))
  theta[,1]=paramsMC[1,]; 
  
  pf=PF_regul(theta[,1],num_ens)
  X[,,1]=pf$x_star; y=pf$y; nwts=pf$nwts;
  p_hat[1]=pf$p_hat_log
  xtrn=pf$xtrn;
  xsd=pf$x_sd; 
  flag=pf$err;  # only need the err record for the best run
  # Step 2: run MCMC
  err_cnt=0;  # count errors
  acc_cum=0   # index to track acceptance rate
  ind_acc=1   # records for unique parameter sets and X's
  m=1
  ## use 'repeat' to do for-loop, in order to control loop counter once error occurs
  repeat {
    m=m+1;
    # print(c("m",m));
    # (a) sample theta_star from a proposal density q(theta_star|theta(m-1))
    # use a ramdom walk for the proposal density, such that q(theta_star|theta(m-1))=q(theta(m-1)|theta_star)
    # theta_star=theta[,m-1]
    # propose new theta by a multipnormial dist.
    # theta_star=theta.prop.multinorm(theta[,m-1],theta_low,theta_up,Sig)
    # propose new theta from the lhs
    theta_star=paramsMC[m,]
    # theta_star=paramsMC[rndMC[m],]; # if more than M combinations
    # propose new theta by independent uniform dist.
    # for (i in 1:5){
    #  theta_star[i]=theta.prop.unif(theta[i,m-1],theta_low[i],theta_up[i],theta_len[i])
    # }
    
    # (b) run Alg. 2 to sample x1:T* from p_hat(x1:T|theta*,z1:T) and 
    # obtain the marginal likelihood estimate p_hat(z1:T|theta*)
    pf=PF_regul(theta_star,num_ens)
    p_hat_star=pf$p_hat_log
    # test if Particle Filter runs properly; if not, try again
    if (is.null(p_hat_star)){
      m=m-1;   
      err_cnt=err_cnt+1;
      next;
    }
    # prior dist: assume the 4 prames (theta vector) is independent normal
    p_old=dnorm(theta[1,m-1],L_mean,L_len)*dnorm(theta[2,m-1],D_mean,D_len)*dnorm(theta[3,m-1],Rmx_mean,Rmx_len)*
      dnorm(theta[4,m-1],Rmn_mean,Rmn_len); #*dnorm(theta[5,m-1],alpha_inv_mean,alpha_inv_len)
    p_star=dnorm(theta_star[1],L_mean,L_len)*dnorm(theta_star[2],D_mean,D_len)*dnorm(theta_star[3],Rmx_mean,Rmx_len)*
      dnorm(theta_star[4],Rmn_mean,Rmn_len); #*dnorm(theta_star[5],alpha_inv_mean,alpha_inv_len)
    # (c) with probability 
    # min(p_hat(z1:T|theta*)*p(theta*)*q(theta(m-1)|theta*)/[p_hat(z1:T|theta(m-1))*p(theta(m-1))*q(theta*|theta(m-1))],1)
    # set theta(m)=theta*, x1:T(m)=x1:T* and p_hat(z1:T|theta(m))=p_hat(z1:T|theta*);
    # else set theta(m)=theta(m-1), x1:T(m)=x1:T(m-1) and p_hat(z1:T|theta(m))=p_hat(z1:T|theta(m-1));
    
    p_accept=min(p_hat_star+log(p_star)-p_hat[m-1]-log(p_old),0)
    u=log(runif(1))
    # if the both the old and new p_hat are -Inf, then p_accept=NaN, accept new parameters
    if (is.na(p_accept)) p_accept=0;
    if (u<p_accept){
      # accept proposal
      acc_cum=acc_cum+1;
      theta[,m]=theta_star;
      
      ## UPDATED HERE
      ### NOTE: PASS THE WEITHS AS WELL: nwts
      X[,,m]=pf$x_star; y=pf$y; nwts=pf$nwts;
      p_hat[m]=pf$p_hat_log;
      xtrn=pf$xtrn;
      xsd=pf$x_sd; # changed only new proposal is accepted.
      flag=pf$err;  # only need the err record for the best run
    } else {
      # keep old parameters
      theta[,m]=theta[,m-1];
      X[,,m]=X[,,m-1];
      p_hat[m]=p_hat[m-1];
    }
    #test
    # print(c(m, theta_star, acc_cum));
    if (m==M) break;
  } # end repeat
  ## zoom in to finer parameter space, centered at theta[,M]
  # paramsMCzi=lhs(MM-M,cbind(theta[,M]-10/M*(theta_up-theta_low),theta[,M]+10/M*(theta_up-theta_low)))
  # use prior to inform the width of parameter space [?? is this kind of cheating?? maybe not, if we have the info, why not us it]
  theta_zi_low=theta[,M]-1/2*abs(theta_mean-theta[,M])
  theta_zi_low=apply(cbind(theta_zi_low,theta_low),1,max)
  theta_zi_up=theta[,M]+1/2*abs(theta_mean-theta[,M])
  theta_zi_up=apply(cbind(theta_zi_up,theta_up),1,min)
  paramsMCzi=lhs(MM-M,cbind(theta_zi_low,theta_zi_up))
  
  # MCMC Step 2
  ## use 'repeat' to do for-loop, in order to control loop counter once error occurs
  repeat {
    m=m+1;
    # print(c("m",m));
    # (a) sample theta_star from a proposal density q(theta_star|theta(m-1))
    # use a ramdom walk for the proposal density, such that q(theta_star|theta(m-1))=q(theta(m-1)|theta_star)
    # theta_star=theta[,m-1]
    # propose new theta by a multipnormial dist.
    # theta_star=theta.prop.multinorm(theta[,m-1],theta_low,theta_up,Sig)
    # propose new theta from the lhs
    theta_star=paramsMCzi[m-M,]; # if more than M combinations: theta_star=paramsMC[rndMC[m],]
    # propose new theta by independent uniform dist.
    # for (i in 1:5){
    #  theta_star[i]=theta.prop.unif(theta[i,m-1],theta_low[i],theta_up[i],theta_len[i])
    # }
    
    # (b) run Alg. 2 to sample x1:T* from p_hat(x1:T|theta*,z1:T) and 
    # obtain the marginal likelihood estimate p_hat(z1:T|theta*)
    pf=PF_regul(theta_star,num_ens)
    p_hat_star=pf$p_hat_log
    # test if Particle Filter runs properly; if not, try again
    if (is.null(p_hat_star)){
      m=m-1;   
      err_cnt=err_cnt+1;
      next;
    }
    # prior dist: assume the 5 prames (theta vector) is independent normal
    p_old=dnorm(theta[1,m-1],L_mean,L_len)*dnorm(theta[2,m-1],D_mean,D_len)*dnorm(theta[3,m-1],Rmx_mean,Rmx_len)*
      dnorm(theta[4,m-1],Rmn_mean,Rmn_len); #*dnorm(theta[5,m-1],alpha_inv_mean,alpha_inv_len)
    p_star=dnorm(theta_star[1],L_mean,L_len)*dnorm(theta_star[2],D_mean,D_len)*dnorm(theta_star[3],Rmx_mean,Rmx_len)*
      dnorm(theta_star[4],Rmn_mean,Rmn_len); #*dnorm(theta_star[5],alpha_inv_mean,alpha_inv_len)
    # (c) with probability 
    # min(p_hat(z1:T|theta*)*p(theta*)*q(theta(m-1)|theta*)/[p_hat(z1:T|theta(m-1))*p(theta(m-1))*q(theta*|theta(m-1))],1)
    # set theta(m)=theta*, x1:T(m)=x1:T* and p_hat(z1:T|theta(m))=p_hat(z1:T|theta*);
    # else set theta(m)=theta(m-1), x1:T(m)=x1:T(m-1) and p_hat(z1:T|theta(m))=p_hat(z1:T|theta(m-1));
    
    p_accept=min(p_hat_star+log(p_star)-p_hat[m-1]-log(p_old),0)
    u=log(runif(1))
    # if the both the old and new p_hat are -Inf, then p_accept=NaN, accept new parameters
    if (is.na(p_accept)) p_accept=0;
    if (u<p_accept){
      # accept proposal
      acc_cum=acc_cum+1;
      theta[,m]=theta_star;
      
      ### UPDATED HERE
      ### NOTE: PASS THE WEITHS AS WELL: nwts
      X[,,m]=pf$x_star; y=pf$y; nwts=pf$nwts;
      p_hat[m]=pf$p_hat_log;
      xtrn=pf$xtrn;
      xsd=pf$x_sd; # changed only new proposal is accepted.
      flag=pf$err;  # only need the err record for the best run
    } else {
      # keep old parameters
      theta[,m]=theta[,m-1];
      X[,,m]=X[,,m-1];
      p_hat[m]=p_hat[m-1];
    }
    
    if (m==MM) break;
  }  
  acc_rate=acc_cum/MM
  
  #### Forecast
  ## resample according to the final nwts, so don't need to adjust for nwts later
  currx=xtrn;   #  Getting current state and parameters
  L=theta[1,MM]; D=theta[2,MM]; # alpha=1/theta[5];
  R0max=theta[3,MM]; R0min=theta[4,MM];
  
  ## Calculate the reproductive number at time t BT1 and the transmission rate 
  ## according to Shaman J, Karspeck A (2012) Proc Natl Acad Sci USA 109: 20425-20430.
  beta.range=tm.range[1]:(tail(tm.range,1)+2*tmstep);
  # make sure to include enought beta, b/c the prior integrate 1 time step more than needed
  # AHpt=AH[beta.range,]; # only use AH during the time of interest, to minimize space needed
  ones1=matrix(1,dim(AHpt)[1],num_ens); # matrix of 1 to make the dimension match
  b=log(R0max-R0min);  # expand b to a matrix with each col for each particle
  a=-180;
  # ones2=matrix(1,1,num_ens); 
  # AHpt=as.matrix(AHpt,length(AHpt),1)
  # AHpt=AHpt%*%ones2
  BT1=exp(a*AHpt+b*ones1)+R0min*ones1;
  beta=BT1/D; # for all particles 
  tcurrent = tm.ini+tmstep*ntrn;
  Sr_tmp=propagateParSIR(tcurrent+dt,tcurrent+tmstep*nfc,dt,currx[1,],currx[2,],N,
                         D=theta[2,MM],L=theta[1,MM],beta,realdata=T)
  
  fcast=array(0,c(3,num_ens,nfc)); # fcast: S, I, newI; 
  fcast[1,,]=t(Sr_tmp$S[tmstep*(1:nfc)+1,]); # num_ens, time 
  fcast[2,,]=t(Sr_tmp$I[tmstep*(1:nfc)+1,]);
  fcast[3,,]=t(Sr_tmp$newI[tmstep*(1:nfc)+1,]-Sr_tmp$newI[tmstep*(0:(nfc-1))+1,]);
  
  fcast_mean=fcast_sd=matrix(0,num_var,nfc)
  for (tt in 1:nfc){
    # fcast_mean[,tt]=apply(fcast[,,tt],1,mean, na.rm=T);
    # fcast_sd[,tt]=apply(fcast[,,tt],1,sd, na.rm=T);
    for(ii in 1:num_var){
      fcast_mean[ii,tt]=sum(fcast[ii,,tt]*nwts[,ntrn]);
      fcast_sd[ii,tt]=sqrt(sum(nwts[,ntrn]*(fcast[ii,,tt]-fcast_mean[ii,tt])^2));
    }
  }
  
  pkwks=rep(0,num_ens);
  obs_pkwk=which.max(obs_i);
  for (i in 1:num_ens){
    pkwks[i]=which.max(c(y[i,1:ntrn],fcast[3,i,1:nfc]));
  }
  pkwk_mode=MODE(pkwks)[1];
  pkwk_mode_perc=MODE(pkwks)[2]/num_ens;
  leadpkwk_mode=pkwk_mode-ntrn;
  ## check the mean as well
  y.all=c(X[3,,MM],fcast_mean[3,]);  # newI 
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
  leadpkwk_var=sum(peak_nwts*(pkwks-pkwk.ens.mean)^2) # note the weighted variance is biased
  # Note: the mean and nwts at the time of pkwk_mean according to the mean trajectory are not matched
  # the mean of ens of pkwks is sum(peak_nwts*pkwks)
  # leadpkwk_var=var(pkwks-obs_pkwk); # Wrong if forecast is made after the predicted peak!!
  
  ones=matrix(1,1,ntrn);
  xsd=rbind(xsd[1:2,],matrix(0,4,1)%*%ones,xsd[num_var,]);
  
  # metrics for comparison
  # Y=c(X[3,,MM],fcast_mean[3,]);  # newI 
  corr=cor(y.all,head(obs_i,nsn)); rms=sqrt(mean((y.all-head(obs_i,nsn))^2));
  delta_sum_newI=sum(y.all)-sum(head(obs_i,nsn));
  delta_pkwk=pkwk_mode-obs_pkwk;
  
  ## output data
  ## output the prediction from the last iteration
  tstep=seq(tm.ini+tmstep,num_times*tmstep+tm.ini,by=tmstep); # time
  fc_start=ntrn+wk_start-1;
  
  out1=cbind(rep(fc_start,ntrn),tstep[1:ntrn],t(X[1:2,,MM]),t(theta[,MM]%*%ones),X[3,,MM]);
  out2=cbind(rep(fc_start,ntrn),tstep[1:ntrn],t(xsd));
  out3=cbind(rep(fc_start,nfc),tstep[(ntrn+1):nsn],t(fcast_mean));
  out4=cbind(rep(fc_start,nfc),tstep[(ntrn+1):nsn],t(fcast_sd));
  out5=t(c(fc_start,rms,corr,delta_sum_newI,delta_pkwk,leadpkwk_mode,leadpkwk_var,pkwk_mode_perc,
           pkwk_mean,leadpkwk_mean,delta_pkwk_mean));
  out6=cbind(rep(fc_start,MM),1:MM,t(theta),p_hat);
  
  colnames(out1)=colnames(out2)=c("fc_start","time","S","I","L","D","R0max","R0min","newI");
  # colnames(out2)=c("fc_start","time","S","I","L","D","R0max","R0min");
  colnames(out3)=colnames(out4)=c("fc_start","time","S","I","newI");
  colnames(out5)=c("fc_start","rms","corr","delta_sum_newI","delta_pkwk",
                   "leadpkwk_mode","leakpkwk_var",'pkwk_mode_perc',"mn_pkwk","mn_leadpk","delta_mn_pkwk");
  colnames(out6)=c("fc_start","iteration","L","D","R0max","R0min","log_like");
  
  if (metricsonly==F){  
    out=list(train=out1,trainsd=out2,fcast=out3,fcsd=out4,metrics=out5,trainparams=out6);
  } else {
    out=list(metrics=out5);
  }
}


