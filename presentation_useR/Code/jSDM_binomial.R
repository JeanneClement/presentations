## ==============================================================================
## author          :Ghislain Vieilledent, Jeanne Clément
## email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

##################################################################################
## betadens #####################################################################
betadens <- function (beta_k, pos_beta, beta_run, mubeta, Vbeta, Y, T, X) {
  # Indicating the rank of the parameter of interest
  k = pos_beta;
  NOBS = nrow(X);
  NP = ncol(X);
  # logLikelihood
  logL = 0.0;
  for (n in 1:NOBS) {
    # theta 
    Xpart_theta = 0.0;
    for ( p in 1:NP) {
      if ( p != k ) {
        Xpart_theta <- Xpart_theta + X[n,p] * beta_run[p];
      }
    }
    Xpart_theta <- Xpart_theta + X[n,k] * beta_k;
    theta <- inv.logit(Xpart_theta);
    # log Likelihood 
    logL <- logL + dbinom(Y[n], T[n], theta, 1);
  }
  # logPosterior = logL + logPrior
  logP <- logL + dnorm(beta_k, mubeta[k], sqrt(Vbeta[k]), 1);
  
  return (logP)
}

###########################################################################################################################
#################################### Gibbs sampler #######################################################################
jSDM_binomial <- function(ngibbs, nthin, nburn,
                          Y, T, X,
                          beta_start, mubeta, Vbeta,
                          seed, ropt, verbose) {
  
  #########################################################
  ##### Defining and initializing objects ################
  
  ### Initialize random number generator 
  set.seed(seed)
  
  ### Redefining constants 
  NGIBBS = ngibbs;
  NTHIN = nthin;
  NBURN = nburn;
  NSAMP = (NGIBBS-NBURN)/NTHIN;
  NOBS = nrow(X);
  NP = ncol(X);
  
  #########################################################
  #### Declaring new objects to store results #############
  
  ### Parameters
  beta <- matrix(0,NSAMP,NP);
  beta_run <- t(beta_start);
  ### Latent variable
  theta_run <- rep(0,NOBS);
  theta_latent <- rep(0,NOBS);
  Deviance <- rep(0,NSAMP);
  
  ############################################################
  ## Proposal variance and acceptance for adaptive sampling ##
  
  ### beta
  sigmap_beta <- rep(1,NP);
  nA_beta <- rep(0,NP);
  Ar_beta <- rep(0,NP); ## Acceptance rate
  
  ### Message 
  print("Running the Gibbs sampler. It may be long, please keep cool :)", quote =F);
  
  ##############################################################
  ###### Gibbs sampler #########################################
  
  for ( g in 1:NGIBBS) {
    
    ####################################################
    ####### beta
    
    for ( p in  1:NP ) {
      pos_beta <- p; 
      x_now <- beta_run[p];
      x_prop = x_now + rnorm(1,0,sigmap_beta[p]);
      p_now = betadens(x_now, pos_beta, beta_run, mubeta,Vbeta, Y, T, X);
      p_prop = betadens(x_prop, pos_beta, beta_run, mubeta, Vbeta, Y, T, X);
      ratio = exp(p_prop - p_now); # ratio
      z = runif(1);
      ## Actualization
      if ( z < ratio ) {
        beta_run[p] <- x_prop;
        nA_beta[p] <- nA_beta[p]+1;
      }
    }
    
    
    #########################################################################
    ######### Deviance ######################################################
    
    # logLikelihood
    logL = 0.0;
    for ( n in 1:NOBS) {
      # theta
      Xpart_theta <- 0.0;
      for ( p in 1:NP) {
        Xpart_theta <- Xpart_theta + X[n, p] * beta_run[p];
      }
      theta_run[n] <- inv.logit(Xpart_theta);
      # log Likelihood 
      logL <- logL + dbinom(Y[n], T[n], theta_run[n], 1);
    }
    
    # Deviance
    Deviance_run <- -2 * logL;
    
    
    #########################################################################
    ######### Output ######################################################
    if ( ((g+1) > NBURN) & (((g+1) %% NTHIN) == 0) ) {
      isamp = ((g+1)-NBURN) / NTHIN;
      beta[isamp-1,] = beta_run;
      Deviance[isamp-1] = Deviance_run;
      # We compute the mean of NSAMP values
      for ( n in 1:NOBS) {
        theta_latent[n] = theta_latent[n] + theta_run[n] / NSAMP;
      }
    }
    
    
    #########################################################################
    ## Adaptive sampling (on the burnin period) ############################
    ROPT = ropt;
    DIV = 0;
    if ( NGIBBS >= 1000 ) {
      DIV=100;
    }
    else {
      DIV = NGIBBS / 10;
    }
    ## During the burnin period 
    if ( ((g+1) %% DIV == 0) & ((g+1) <= NBURN) ) {
      for (p in 1:NP) {
        Ar_beta[p] = nA_beta[p] / DIV;
        if ( Ar_beta[p] >= ROPT ){
          sigmap_beta[p] = sigmap_beta[p] * (2-(1-Ar_beta[p]) / (1-ROPT)); 
        } 
        else {
          sigmap_beta[p] = sigmap_beta[p]  / (2-Ar_beta[p]  / ROPT);}
        # We reinitialize the number of acceptance to zero
        nA_beta[p] = 0.0; 
      }
    }
    ## After the burnin period
    if ( ((g+1) %% DIV == 0) & ((g+1) > NBURN) ) {
      for (p in 1:NP) {
        Ar_beta[p] = nA_beta[p] / DIV;
        # We reinitialize the number of acceptance to zero
        nA_beta[p] = 0.0;
      }
    }
    
    
    #########################################################################
    ## Progress bar ########################################################
    Perc = 100 * (g+1) / (NGIBBS);
    if ( ((g+1) %% (NGIBBS/100) == 0) & (verbose == 1)) {
      cat("*", sep ="");
      if( (g+1) %% (NGIBBS/10) == 0 ) {
        # Mean acceptance rate
        mAr_beta=0; 
        for ( p in 1:NP) {
          mAr_beta = mAr_beta + Ar_beta[p] / NP;
        }
        print(paste(Perc, "%", "mean accept. rates = beta:", mAr_beta), quote =F);
      }
    }
    
    
  } # Gibbs sampler
  
  # Return results as a Rcpp::List
  results = list(beta = beta, Deviance = Deviance,
                 theta_latent = theta_latent);
  
  return (results) 
} # end jSDM_binomial function

library(coda)

inv.logit <- function(x, min=0, max=1) {
  p <- exp(x)/(1+exp(x))
  p <- ifelse( is.na(p) & !is.na(x), 1, p ) # fix problems with +Inf
  p * (max-min) + min
}

nsite <- 200
seed <- 1234
set.seed(seed)
visits<- rpois(nsite,3)
visits[visits==0] <- 1

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- cbind(rep(1,nsite),x1,x2)
beta.target <- c(-1,1,-1)
logit.theta <- X %*% beta.target
theta <- inv.logit(logit.theta)
Y <- rbinom(nsite,visits,theta)

# Data-sets
data.obs <- data.frame(Y,visits,x1,x2)

# Iterations
nsamp <- 1000
nburn <- 1000
nthin <- 1
ngibbs <- nsamp+nburn

mf.suit <- model.frame(formula=~x1+x2, data=data.obs)
X <- model.matrix(attr(mf.suit,"terms"), data=mf.suit)

# Call to C++ function
mod <- jSDM_binomial(
  ngibbs=ngibbs, nthin=nthin, nburn=nburn,
  Y=data.obs$Y,
  T=data.obs$visits,
  X=X,
  beta_start=rep(0, 3),
  mubeta=rep(0, 3), Vbeta=rep(1.0E6, 3),
  seed=1234, ropt=0.44, verbose=1)

# Parameter estimates
MCMC <- mcmc(mod$beta, start=nburn+1, end=ngibbs, thin=nthin)
summary(MCMC)
mean(mod$Deviance)
plot(MCMC)

# GLM resolution to compare
mod.glm <- glm(cbind(Y,visits-Y)~x1+x2,family="binomial",data=data.obs)
summary(mod.glm)

# Predictions
plot(theta,mod$theta_latent, main="theta",
     xlab="obs", ylab="fitted")
abline(a=0,b=1,col="red")
