## remove previous list and set working directory
rm(list=ls())

## libraries used
library(statmod)
library(MASS)
library(coda)
library(GeneralizedHyperbolic)
library(rmutil)
library(MCMCpack)
library(mclust)
library(class)
library(teigen)
library(MixGHD)
library(truncnorm)

######################################
###                            #######
### functions to be used later #######
###                            #######
######################################
### Likelihood function
### likelihood calculation function
lik <- function(mu,beta,gamma,delta,Delt,alpha,sq_phi,lam,y){
  first <- delta*(1/(2^((p-1)/2)))*exp(gamma*delta+(y-mu)%*%(beta))
  sec <- (alpha/(pi*delta*sq_phi))^(p+1)/2
  third <- besselK(alpha*delta*sq_phi,lam)
  product <- first * sec * third
  if(third == 0){
    product <- 0
  }
  return(product)
}

### sample group labels
sample_z <- function(pi_g,mu,beta,gamma,delta,Delt,G){
  z <- c()
  alpha_g <- c()
  phi_y_g <- matrix(ncol=G, nrow=n)
  forz <- matrix(ncol=G, nrow=n)
  for (grp in 1:G) {
    alpha_g[grp] = sqrt(gamma[[grp]]^2+t(beta[[grp]])%*%Delt[[grp]]%*%beta[[grp]])
    lambda = -(p+1)/2
    for (i in 1:n) {
      phi_y_g[i,grp]<-1+(1/delta[[grp]]^2)*t(Y[i,]-mu[[grp]])%*%solve(Delt[[grp]])%*%(Y[i,]-mu[[grp]])
      likelihood<-lik(mu=mu[[grp]],beta=beta[[grp]],gamma=gamma[[grp]],delta=delta[[grp]],Delt=Delt[[grp]], 
                      alpha=alpha_g[grp],sq_phi=sqrt(phi_y_g[i,grp]),lam=lambda,y=Y[i,])
      forz[i,grp] <- pi_g[grp]*likelihood
      # print(i)
    }
  }
  if(G == 1){
    z <- as.matrix(rep(1,n))
    # forz[which(forz == 0),] <- 1e-300
  }else{
    forclass<-forz/rowSums(forz)
    if(any(rowSums(forz) == 0)){
      # print(forclass)
      forclass[which(rowSums(forz)==0),] <- rep(0,G)
      forclass[which(rowSums(forz)==0),1] <- 1
    }
    # put constrains when G = 1
    # print(forclass)##sample probability for class labels
    for(i in 1:n){
      z[i]<-sample(x=1:G,size=1,prob=forclass[i,])
    }
  }
  return(list(z,forz))
}

### update all parameters
estm <- function(delta,gamma,mu,beta,Delt,Y,z_g){
  #### Default values given for hyperparameters
  nu_0 = p + 2
  L_0 = diag(p)
  a_0_old = 1
  a_1_old = 0
  a_2_old = 0
  a_3_old = 1
  a_4_old = 1 
  
  alpha = sqrt(gamma^2+t(beta)%*%Delt%*%beta)
  phi_y <- rep(0, n)
  for (i in 1:n) {
    phi_y[i]<-1+(1/delta^2)*t(Y[i,]-mu)%*%solve(Delt)%*%(Y[i,]-mu)
  }
  ###########################################################################
  #### update the latent variable U and U^(-1)
  #### using the GIG distribution sampling trick
  #### U follows a GIG(lambda,chi,psi) distribution has a density function of 
  #### u^(lambda-1)*exp{-(1/2)*[chi/u+psi*u]}
  lambda = -(p + 1)/2
  U_gig_chi = delta^2*phi_y # GIG chi parameter
  U_gig_psi = alpha^2 # GIG psi parameter
  #### check MCMC convergence
  U_check <- c()
  U_inv <- rgamma(n,1,1) # start from a gamma distribution
  maxit <- 200
  i <- 1
  check <- 0
  while(check < 1){
    if(i < (maxit+1)){
      S <- rgamma(n,shape = -lambda,rate = U_gig_chi/2)
      V <- rgamma(n,shape = -lambda,rate = U_gig_psi/2)
      U_inv <- S+1/(V+1/U_inv)
      U <- 1/U_inv
    }
    U_check[i] <- mean(U)
    if(i > maxit){
      conv<-heidel.diag(U_check)
      if (conv[,1]==1){
        if (conv[,4]==1) check<-1 else check<-0} else check<-0
        maxit <- maxit+50
    }
    i <- i+1
    if(maxit > 5000){check <- 1}
  }
  
  ###############################
  #### update the hyperparameters
  a_0 = a_0_old + sum(z_g)
  a_1 = a_1_old + colSums(z_g * Y)
  a_2 = a_2_old + t(z_g * Y) %*% U_inv
  a_3 = a_3_old + sum(z_g * U)/2
  a_4 = a_4_old + sum(z_g * U_inv)/2
  
  #### update delta,gamma
  delta <- sqrt(rgamma(1,shape=(a_0/2+1),rate=(a_4-(a_0^2/(4*a_3)))))
  
  require(truncnorm)
  mu_ga <- (a_0* delta)/(2*a_3)
  sigma_ga <- sqrt(1/(2*a_3))
  # aa = - mu_ga/sigma_ga
  gamma <-rtruncnorm(1, a=0, b=Inf, mean = mu_ga, sd = sigma_ga)
  
  # gamma <- sqrt(rgamma(1,shape = (a_0/2),rate = (-a_0 + a_3 + a_4)))
  # delta <- gamma
  
  # #### update beta
  # mean_beta <- (solve(Delt)/(2*a_3))%*%(a_1-a_0*mu)
  # sigm_beta <- solve(Delt)/(2*a_3)
  # beta <- mvrnorm(1,mean_beta,sigm_beta)
  # # print(beta)
  # 
  # #### update mu
  # mean_mu <- (Delt/(2*a_4))%*%(solve(Delt)%*%a_2 - a_0*beta)
  # sigm_mu <- Delt/(2*a_4)
  # mu <- mvrnorm(1,mean_mu,sigm_mu)
  # # print(mu)
  
  #### update mubeta jointly
  rho <- 4*a_3*a_4-a_0*a_0
  for_mean <- c((2*a_2*a_3-a_0*a_1)/rho, solve(Delt)%*%(2*a_1*a_4-a_0*a_2)/rho)
  sig_m <- 2*a_3*Delt/rho
  sig_b <- 2*a_4*solve(Delt)/rho
  sig_mb <- -a_0*diag(1,p)/rho
  for_sigma <- rbind(cbind(sig_m,sig_mb),cbind(sig_mb,sig_b))
  sample_mubeta <- mvrnorm(1,for_mean,for_sigma)
  mu <- sample_mubeta[1:p]
  beta <- sample_mubeta[(p+1):(2*p)]
  
  #### update Delt 
  S_0 <- matrix(rep(0,p^2),nrow=p,ncol=p)
  S <- matrix(nrow =p,ncol=p)
  for (i in 1:n) {
    S = (Y[i,]-mu)%*%t(Y[i,]-mu)*z_g[i]*U_inv[i]
    S_0 = S_0 + S
  }
  lambda_gig <- -(nu_0+sum(z_g))/2
  A = (L_0+S_0)/2
  B = a_3
  z <- beta
  gig_chi <-  2*B
  gig_psi <- 2*t(z)%*%A%*%z
  param_gig = c(gig_chi,gig_psi,lambda_gig-(1-p)/2)
  ### include a try(silent = TRUE)
  X_gig <- GeneralizedHyperbolic::rgig(1, param = param_gig)
  Y_gig <- rWishart(1,-2*lambda_gig,solve(2*A))
  k <- solve(z %*% as.matrix(X_gig) %*% t(z) + Y_gig[, , 1])
  Delt <- k/det(k)^(1/p)
  # Delt <- k
  # print(Delt)
  
  return(list(delta,gamma,mu,beta,Delt))
}

### update all parameters and class labels
update_all_parameters <- function(parameters,z_clus,pi_g,G){
  #### order everything from the previous step based on the first variate of mu
  mu.temp <- parameters[[3]]
  # print(mu.temp)
  orders <- rank(-sapply(mu.temp,'[[',1))
  pi_g.temp <- pi_g
  z_clus.temp <- z_clus
  # True_Para.temp <- True_Para
  
  delta_g <- list()
  gamma_g <- list()
  mu_g <- list()
  beta_g <- list()
  Delt_g <- list()
  p_ig<-numeric(G)
  z_clus <- matrix(nrow = n, ncol = G)
  # True_Para <- list()
  
  for(grp in 1:G){
    order_ind <- which(orders == grp)
    pi_g[grp] <- pi_g.temp[order_ind]
    z_clus[,grp] <- z_clus.temp[,order_ind]
    delta_g[[grp]] <- parameters[[1]][[order_ind]]
    gamma_g[[grp]] <- parameters[[2]][[order_ind]]
    mu_g[[grp]] <- parameters[[3]][[order_ind]]
    beta_g[[grp]] <- parameters[[4]][[order_ind]]
    Delt_g[[grp]] <- parameters[[5]][[order_ind]]
    # True_Para[[grp]] <- True_Para.temp[[order_ind]]
  }
  
  ##### update parameters (including hyperparameters and U,U_inv)
  for (grp in 1:G) {
    z <- z_clus[, grp]
    delta <- delta_g[[grp]]
    gamma <- gamma_g[[grp]]
    mu <- mu_g[[grp]]
    beta <- beta_g[[grp]]
    Delt <- Delt_g[[grp]]
    
    est <- estm(delta=delta,gamma=gamma,mu=mu,beta=beta,Delt=Delt,Y=Y,z_g=z)
    
    delta_g[[grp]] <- est[[1]]
    gamma_g[[grp]] <- est[[2]]
    mu_g[[grp]] <- est[[3]]
    beta_g[[grp]] <- est[[4]]
    Delt_g[[grp]] <- est[[5]]
  }
  
  parameters <- list()
  parameters[[1]] <- delta_g
  parameters[[2]] <- gamma_g
  parameters[[3]] <- mu_g
  parameters[[4]] <- beta_g
  parameters[[5]] <- Delt_g
  
  para_output <- parameters
  
  #### update class labels and loglikelihood, and pi_g
  samplez<-sample_z(pi_g=pi_g,mu=mu_g,beta=beta_g,gamma=gamma_g,delta=delta_g,Delt=Delt_g,G=G)
  class<-samplez[[1]]
  forz <- samplez[[2]]
  
  z_clus<-matrix(nrow = n,ncol = G)
  for(grp in 1:G){
    z_clus[,grp] <- as.numeric(class == grp)
  }
  loglik<-sum(log(rowSums(forz)))
  
  forpi <- colSums(z_clus)
  # pi_g <- numeric(G)
  pi_g <- rdirichlet(1, alpha=(forpi+1))
  
  updatedallresult <- list(para_output,pi_g,class,z_clus,loglik)
  return(updatedallresult)
}

### remove small grps 
remove_small_group <- function(parameters,rem,G){
  delta_g <- parameters[[1]]
  gamma_g <- parameters[[2]]
  mu_g <- parameters[[3]]
  beta_g <- parameters[[4]]
  Delt_g <- parameters[[5]]
  
  G<-G-length(rem)
  
  delta_g <- delta_g[-rem]
  gamma_g <- gamma_g[-rem]
  mu_g <- mu_g[-rem]
  beta_g <- beta_g[-rem]
  Delt_g <- Delt_g[-rem]
  
  parameters[[1]] <- delta_g
  parameters[[2]] <- gamma_g
  parameters[[3]] <- mu_g
  parameters[[4]] <- beta_g
  parameters[[5]] <- Delt_g
  
  return(list(parameters,G))
}

### return max BIC index
arrayIndex <- function(i, dim) { 
  ndim <- length(dim);      # number of dimension 
  v <- cumprod(c(1,dim));  # base 
  
  # Allocate return matrix 
  j <- matrix(NA, nrow=length(i), ncol=ndim); 
  
  i <- (i-1);     # one-based indices 
  for (kk in 1:ndim) 
    j[,kk] <- (i %% v[kk+1])/v[kk]; 
  1 + floor(j);  # one-based indices 
} 

###################################################################################
##########################
####                  ####
#### Simulation Study ####
####                  ####
##########################
seed.no=2019
set.seed(seed.no)
# seed.no=1995
# for(seed.no in 1242:1251)
# set.seed(seed.no)

###################################################################
## generate two groups of observtions from the MNIG distribution ##
###################################################################
p = 2
### grp 1
n = 500
t_delta <- 0.9 # true delta grp 1
t_gamma <- 0.9 # true gamma grp 1
t_mu <- c(2, 5) # true mu
t_beta <- c(0.5, 0.5) #true beta
# t_beta <- c(0,0)
# t_Delta <- matrix(c(1, 0, 0, 1), p, p)
t_Delta <- matrix(c(2, -1, -1, 1), p, p)
# latent variable
u <- statmod::rinvgauss(n, mean = t_delta/t_gamma, shape = t_delta^2)
Y_1 <- matrix(nrow = n, ncol = p)
for (i in 1:n) {
  Y_1[i, ] = mvrnorm(1, (t_mu + u[i] * t_Delta %*% t_beta), (u[i] * t_Delta))
}
### grp2
n_2 = 500
t_delta_2 <- 1.2 #true delta
t_gamma_2 <- 1.2 #true gamma
t_mu_2 <- c(-4, 3) #true mu
t_beta_2 <- c(-0.5, -0.1) #true beta
# t_beta_2 <- c(0,0)
# t_Delta_2 <- matrix(c(2, 1, 1, 1), p, p)
t_Delta_2 <- matrix(c(1, 0, 0, 1), p, p) #true Delt determinant = 1
u_2 <- statmod::rinvgauss(n_2, mean = t_delta_2/t_gamma_2, shape = t_delta_2^2)
Y_2 <- matrix(nrow = n_2, ncol = p)
for (i in 1:n_2) {
  Y_2[i, ] = mvrnorm(1, (t_mu_2 + u_2[i] * t_delta_2 %*% t_beta_2), (u_2[i] * t_Delta_2))
}
### true group labels of the observations
true_lab <- c(rep(1, n), rep(2, n_2))
### combine two grps and plot the data
Y <- rbind(Y_1, Y_2)
n <- nrow(Y)
p <- ncol(Y)
plot(Y,col=true_lab)


########################
#### Gibbs Sampling ####
########################
Gibbs_output <- list()#store the parameter updated for each component and each iteration
allG_out <- list() #store things need for more iterations at the end of first 1000 iteration
allG_class <- list() #store class labels from each iteration
store_all <- list() #store everything
collectG <- list()

for(G_act in 1:5){# run through 4 possible number of components
  it_max<-1000
  # it_max <- 500
  all_loglik <- list()
  Gibbs_output[[G_act]] <- list()
  allG_out[[G_act]] <- list()
  allG_class[[G_act]] <- list()
  store_all[[G_act]] <- list()
  collectG[[G_act]] <- numeric(3)
  
  for(chain.no in 1:3){
    all_loglik[[chain.no]] <- numeric(it_max)
    Gibbs_output[[G_act]][[chain.no]] <- list()
    allG_out[[G_act]][[chain.no]] <- list()
    allG_class[[G_act]][[chain.no]] <- matrix(nrow=n,ncol=it_max)
    store_all[[G_act]][[chain.no]] <- list()
  }
  
  for(chain.no in 1:3){
    # set.seed(chain.no*seed.no)
    G <- G_act
    k_means <- kmeans(Y,centers=G)
    while(any(k_means$size == 1)){
      k_means <- kmeans(Y,centers=G)
    }
    ### Crude Estimate (i.e Starting values)
    delta_g <- list()
    gamma_g <- list()
    mu_g <- list()
    beta_g <- list()
    Delt_g <- list()
    pi_g<-numeric(G)
    for (grp in 1:G){
      obs<-which(k_means$cluster==grp)
      delta_g[[grp]] <- 1
      gamma_g[[grp]] <- 1
      mu_g[[grp]] <- colMeans(Y[obs,])
      beta_g[[grp]] <- rep(0.05,p)
      Delt_g[[grp]] <- cov(Y[obs,])/(det(cov(Y[obs,]))^(1/p))
      if(rcond(Delt_g[[grp]]) <= 1e-12){Delt_g[[grp]] <- diag(1,p)}
      pi_g[grp]<-length(obs)/n
    }
    
    ### get starting values for class label
    samplez<-sample_z(pi_g=pi_g,mu=mu_g,beta=beta_g,gamma=gamma_g,delta=delta_g,Delt=Delt_g,G=G)
    class<-samplez[[1]]
    while(!all(table(class)>1)|length(table(class))!=G){
      class<-sample_z(pi_g=pi_g,mu=mu_g,beta=beta_g,gamma=gamma_g,delta=delta_g,Delt=Delt_g,G=G)[[1]]
    }
    z_clus<-matrix(nrow = n, ncol = G)
    for(grp in 1:G){
      z_clus[,grp] <- as.numeric(class == grp)
    }
    z_clus.init <- z_clus
    pi_g.init <- pi_g
    
    ### iterations
    it <- 1
    while(it<(it_max+1)){
      if(it != 1){
        para_g <- Gibbs_output[[G_act]][[chain.no]][[it-1]]
        z_clus <- store_all[[G_act]][[chain.no]][[it-1]][[4]]
        pi_g <- store_all[[G_act]][[chain.no]][[it-1]][[2]]
      }else{
        para_g <- list(delta_g,gamma_g,mu_g,beta_g,Delt_g)
        z_clus <- z_clus.init
        pi_g <- pi_g.init
      }
      parameter_updated<-update_all_parameters(parameters=para_g,z_clus=z_clus,pi_g=pi_g,G=G)
      pi_g <- parameter_updated[[2]]
      # print(pi_g)
      z_clus <- parameter_updated[[4]]
      all_loglik[[chain.no]][it] <- parameter_updated[[5]]
      
      if(any(colSums(z_clus)<2)){
        rem <- which(colSums(z_clus)<2)
        z_clus <- as.matrix(z_clus[,-rem])
        pi_g <- pi_g[-rem]
        
        afterRem <- remove_small_group(parameters=parameter_updated[[1]],
                                       rem=rem,G=G)
        parameter_updated[[1]] <- afterRem[[1]]
        parameter_updated[[2]] <- pi_g
        parameter_updated[[3]] <- mclust::map(z_clus)
        parameter_updated[[4]] <- z_clus
        G <- afterRem[[2]]
      }
      
      class <- parameter_updated[[3]]
      allG_class[[G_act]][[chain.no]][,it] <- class
      Gibbs_output[[G_act]][[chain.no]][[it]] <- parameter_updated[[1]]
      store_all[[G_act]][[chain.no]][[it]] <- parameter_updated
      
      cat("Seed",seed.no,"\n")
      cat("Group",G,"\n")
      cat("Chain",chain.no,"\n")
      cat("Iteration", it, "\n")
      
      if(it == it_max){
        allG_out[[G_act]][[chain.no]][[1]] <- Gibbs_output[[G_act]][[chain.no]][[it]]
        allG_out[[G_act]][[chain.no]][[2]] <- all_loglik[[chain.no]]
        allG_out[[G_act]][[chain.no]][[3]] <- z_clus
      }
      it <- it+1
    }
    collectG[[G_act]][chain.no] <- G
  }
  
  # plot likelihood for all chains
  plot(all_loglik[[1]],type = 'l')
  lines(all_loglik[[2]],col='red')
  lines(all_loglik[[3]],col='blue')

  ### check convergence
  # burnin <- 1:(it_max-500)
  burnin <- 1:(0.5*it_max)
  object <- list()
  negInfLik <- 0
  for(chain.no in 1:3){
    if(any(allG_out[[G_act]][[chain.no]][[2]] == -Inf)){negInfLik <- 1}
    object[[chain.no]]=mcmc(allG_out[[G_act]][[chain.no]][[2]][-burnin])
  }
  obj <- mcmc.list(object)
  conv <- gelman.diag(obj)
  check.thresh <- 1.1
  if(negInfLik == 1){
    check.conv <- 1
  }
  else{
    check.conv <- conv[[1]][1]
  }
  check.thresh <- 1.1 
  
  ### add 100 more iteration if gelman.rubin factor larger than 1.1
  while(check.conv > check.thresh){
    it_old <- it
    it_max <- it_max+100
    
    for(chain.no in 1:3){
      # set.seed(chain.no*seed.no)
      G<-collectG[[G_act]][chain.no]
      z_clus <- allG_out[[G_act]][[chain.no]][[3]]
      all_loglik[[chain.no]] <- allG_out[[G_act]][[chain.no]][[2]]
      
      it <- it_old
      while(it<(it_max+1)){
        para_g <- Gibbs_output[[G_act]][[chain.no]][[it-1]]
        z_clus <- store_all[[G_act]][[chain.no]][[it-1]][[4]]
        pi_g <- store_all[[G_act]][[chain.no]][[it-1]][[2]]
        parameter_updated<-update_all_parameters(parameters=para_g,z_clus=z_clus,pi_g=pi_g,G=G)
        pi_g <- parameter_updated[[2]]
        # print(pi_g)
        z_clus <- parameter_updated[[4]]
        all_loglik[[chain.no]]<-append(all_loglik[[chain.no]],parameter_updated[[5]])
        
        if(any(colSums(z_clus)<2)){
          rem <- which(colSums(z_clus)<2)
          z_clus <- as.matrix(z_clus[,-rem])
          pi_g <- pi_g[-rem]
          
          afterRem <- remove_small_group(parameters=parameter_updated[[1]],
                                         rem=rem,G=G)
          parameter_updated[[1]] <- afterRem[[1]]
          parameter_updated[[2]] <- pi_g
          parameter_updated[[3]] <- mclust::map(z_clus)
          parameter_updated[[4]] <- z_clus
          G <- afterRem[[2]]
        }
        
        class <- parameter_updated[[3]]
        allG_class[[G_act]][[chain.no]]<-cbind(allG_class[[G_act]][[chain.no]],class)
        Gibbs_output[[G_act]][[chain.no]][[it]] <- parameter_updated[[1]]
        store_all[[G_act]][[chain.no]][[it]] <- parameter_updated
        
        # cat("Seed",seed.no,"\n")
        cat("Group",G,"\n")
        cat("Chain",chain.no,"\n")
        cat("Iteration", it, "\n")
        
        if(it == it_max){
          allG_out[[G_act]][[chain.no]][[1]] <- Gibbs_output[[G_act]][[chain.no]][[it]]
          allG_out[[G_act]][[chain.no]][[2]] <- all_loglik[[chain.no]]
          allG_out[[G_act]][[chain.no]][[3]] <- z_clus
        }
        it <- it+1
      }
      collectG[[G_act]][chain.no] <- G
    }
    
    # plot likelihood for all chains
    plot(all_loglik[[1]],type = 'l')
    lines(all_loglik[[2]],col='red')
    lines(all_loglik[[3]],col='blue')
    
    ### check convergence
    # burnin <- 1:(it_max-500)
    burnin <- 1:(0.5*it_max)
    object <- list()
    negInfLik <- 0
    for(chain.no in 1:3){
      if(any(allG_out[[G_act]][[chain.no]][[2]] == -Inf)){negInfLik <- 1}
      object[[chain.no]]=mcmc(allG_out[[G_act]][[chain.no]][[2]][-burnin])
    }
    obj <- mcmc.list(object)
    conv <- gelman.diag(obj)
    if(it_max == 10000|negInfLik == 1){check.conv <- 1}
    else{
      check.conv <- conv[[1]][1]
    }
    check.thresh <- 1.1
  }
  
  it_old <- it
  it_max <- it_max+500
  
  for(chain.no in 1:3){
    # set.seed(chain.no*seed.no)
    G<-collectG[[G_act]][chain.no]
    z_clus <- allG_out[[G_act]][[chain.no]][[3]]
    all_loglik[[chain.no]] <- allG_out[[G_act]][[chain.no]][[2]]
    
    it <- it_old
    while(it<(it_max+1)){
      para_g <- Gibbs_output[[G_act]][[chain.no]][[it-1]]
      z_clus <- store_all[[G_act]][[chain.no]][[it-1]][[4]]
      pi_g <- store_all[[G_act]][[chain.no]][[it-1]][[2]]
      parameter_updated<-update_all_parameters(parameters=para_g,z_clus=z_clus,pi_g=pi_g,G=G)
      pi_g <- parameter_updated[[2]]
      # print(pi_g)
      z_clus <- parameter_updated[[4]]
      all_loglik[[chain.no]]<-append(all_loglik[[chain.no]],parameter_updated[[5]])
      
      if(any(colSums(z_clus)<2)){
        rem <- which(colSums(z_clus)<2)
        z_clus <- as.matrix(z_clus[,-rem])
        pi_g <- pi_g[-rem]
        
        afterRem <- remove_small_group(parameters=parameter_updated[[1]],
                                       rem=rem,G=G)
        parameter_updated[[1]] <- afterRem[[1]]
        parameter_updated[[2]] <- pi_g
        parameter_updated[[3]] <- mclust::map(z_clus)
        parameter_updated[[4]] <- z_clus
        G <- afterRem[[2]]
      }
      
      class <- parameter_updated[[3]]
      allG_class[[G_act]][[chain.no]]<-cbind(allG_class[[G_act]][[chain.no]],class)
      Gibbs_output[[G_act]][[chain.no]][[it]] <- parameter_updated[[1]]
      store_all[[G_act]][[chain.no]][[it]] <- parameter_updated
      
      # cat("Seed",seed.no,"\n")
      cat("Group",G,"\n")
      cat("Chain",chain.no,"\n")
      cat("Iteration", it, "\n")
      
      if(it == it_max){
        allG_out[[G_act]][[chain.no]][[1]] <- Gibbs_output[[G_act]][[chain.no]][[it]]
        allG_out[[G_act]][[chain.no]][[2]] <- all_loglik[[chain.no]]
        allG_out[[G_act]][[chain.no]][[3]] <- z_clus
      }
      it <- it+1
    }
    collectG[[G_act]][chain.no] <- G
  }
  
  ### export traceplot for the whole chain for each G
  pdf(file=paste0("./Traceplot_",seed.no,"_G=",G_act,".pdf",sep=""))
  plot(all_loglik[[1]],type = 'l')
  lines(all_loglik[[2]],col='red')
  lines(all_loglik[[3]],col='blue')
  dev.off()
  
  rm(all_loglik)
  rm(class)
  rm(z_clus)
}

#####################
### calculate BIC ###
#####################
BIC <- matrix(nrow = 3,ncol = 5)

for(G in 1:5){
  delta_exp <- list()
  gamma_exp <-list()
  mu_exp <- list()
  beta_exp <- list()
  Delt_exp <- list()
  class_labs <- list()
  
  ### get parameters and class labels from last 400 interations of all three chains
  for(chain.no in 1:3){
    zz <- collectG[[G]][chain.no]
    len <- length(allG_out[[G]][[chain.no]][[2]])
    allParameters <- Gibbs_output[[G]][[chain.no]][(len-399):len]
    class_labs[[chain.no]] <- allG_class[[G]][[chain.no]][,c((len-399):len)]
    
    class.temp <- numeric(n)
    for(i in 1:n){
      class.temp[i] <- as.numeric(names(which.max(table(class_labs[[chain.no]][i,]))))
    }
    z_clus.temp <- mclust::unmap(class.temp)
    forpi <- colSums(z_clus.temp)
    pi_g <- rdirichlet(1, alpha=forpi+1)
    
    delta_exp[[chain.no]]<-matrix(nrow=400,ncol=zz)
    gamma_exp[[chain.no]]<-matrix(nrow=400,ncol=zz)
    mu_exp[[chain.no]]<-matrix(nrow=400,ncol=p*zz)
    beta_exp[[chain.no]]<-matrix(nrow=400,ncol=p*zz)
    Delt_exp[[chain.no]]<-matrix(nrow=400,ncol=p*p*zz)
    
    for(i in 1:400){
      for(grp in 1:zz){
        delta_exp[[chain.no]][i,grp]<-allParameters[[i]][[1]][[grp]]
        gamma_exp[[chain.no]][i,grp]<-allParameters[[i]][[2]][[grp]]
        mu_exp[[chain.no]][i,(1+p*(grp-1)):(p*grp)]<-allParameters[[i]][[3]][[grp]]
        beta_exp[[chain.no]][i,(1+p*(grp-1)):(p*grp)]<-allParameters[[i]][[4]][[grp]]
        Delt_exp[[chain.no]][i,(1+p*p*(grp-1)):(p*p*grp)]<-as.vector(allParameters[[i]][[5]][[grp]])
      }
    }
    delta <- colMeans(delta_exp[[chain.no]])
    gamma <- colMeans(gamma_exp[[chain.no]])
    mu <- colMeans(mu_exp[[chain.no]])
    beta <- colMeans(beta_exp[[chain.no]])
    Delta <- colMeans(Delt_exp[[chain.no]])
    
    delta_g <- list()
    gamma_g <- list()
    mu_g <- list()
    beta_g <- list()
    Delt_g <- list()
    for (grp in 1:zz){
      delta_g[[grp]] <- delta[grp]
      gamma_g[[grp]] <- gamma[grp]
      mu_g[[grp]] <- mu[(1+p*(grp-1)):(p*grp)]
      beta_g[[grp]] <- beta[(1+p*(grp-1)):(p*grp)]
      k <- matrix(Delta[(1+p*p*(grp-1)):(p*p*grp)],nrow=p,ncol=p)
      Delt_g[[grp]] <- k/det(k)^(1/p)
    }
    
    samplez<-sample_z(pi_g=pi_g,mu=mu_g,beta=beta_g,gamma=gamma_g,delta=delta_g,Delt=Delt_g,G=zz)
    forz <- samplez[[2]]
    loglik<-sum(log(rowSums(forz)))
    
    par = 2*p*zz+p*p*zz+zz+(zz-1)
    BIC[chain.no,G] <- 2*loglik-par*log(n)
  }
}

############################
### parameter estimation ###
############################
maxBICind = arrayIndex(which.max(BIC), dim=dim(BIC)) 
zzB = maxBICind[2]
ccB = maxBICind[1]
zz=zzB
cc=ccB

GlobalMaxBIC <- collectG[[zz]][cc]

delta_exp<-matrix(nrow=400,ncol=collectG[[zz]][cc])
gamma_exp<-matrix(nrow=400,ncol=collectG[[zz]][cc])
mu_exp<-matrix(nrow=400,ncol=p*collectG[[zz]][cc])
beta_exp<-matrix(nrow=400,ncol=p*collectG[[zz]][cc])
Delt_exp<-matrix(nrow=400,ncol=p*p*collectG[[zz]][cc])

len <- length(allG_out[[zz]][[1]][[2]])

### get parameters and class labels from last 400 interations of all three chains
for(chain.no in cc:cc){
  allParameters <- Gibbs_output[[zz]][[chain.no]][(len-399):len]
  class_labs[[chain.no]] <- allG_class[[zz]][[chain.no]][,c((len-399):len)]

  for(i in 1:400){
    for(grp in 1:collectG[[zz]][cc]){
      delta_exp[i,grp]<-allParameters[[i]][[1]][[grp]]
      gamma_exp[i,grp]<-allParameters[[i]][[2]][[grp]]
      mu_exp[i,(1+p*(grp-1)):(p*grp)]<-allParameters[[i]][[3]][[grp]]
      beta_exp[i,(1+p*(grp-1)):(p*grp)]<-allParameters[[i]][[4]][[grp]]
      Delt_exp[i,(1+p*p*(grp-1)):(p*p*grp)]<-as.vector(allParameters[[i]][[5]][[grp]])
    }
  }
}

delta <- colMeans(delta_exp)
gamma <- colMeans(gamma_exp)
mu <- colMeans(mu_exp)
beta <- colMeans(beta_exp)
Delta <- colMeans(Delt_exp)

### get the empirical 0.025-qurtile and .975-quartile
del_l_var<-numeric(collectG[[zz]][cc])
del_u_var<-numeric(collectG[[zz]][cc])
gam_l_var<-numeric(collectG[[zz]][cc])
gam_u_var<-numeric(collectG[[zz]][cc])
mu_l_var<-numeric(p*collectG[[zz]][cc])
mu_u_var<-numeric(p*collectG[[zz]][cc])
beta_l_var<-numeric(p*collectG[[zz]][cc])
beta_u_var<-numeric(p*collectG[[zz]][cc])
Delt_l_var<-numeric(p*p*collectG[[zz]][cc])
Delt_u_var<-numeric(p*p*collectG[[zz]][cc])

for(j in 1:collectG[[zz]][cc]){
  del_l_var[j] <- delta_exp[which(rank(delta_exp[,j]) == 10),j]
  gam_l_var[j] <- gamma_exp[which(rank(gamma_exp[,j]) == 10),j]
  del_u_var[j] <- delta_exp[which(rank(delta_exp[,j]) == 390),j]
  gam_u_var[j] <- gamma_exp[which(rank(gamma_exp[,j]) == 390),j]
}

for(j in 1:(p*collectG[[zz]][cc])){
  mu_l_var[j] <- mu_exp[which(rank(mu_exp[,j]) == 10),j]
  mu_u_var[j] <- mu_exp[which(rank(mu_exp[,j]) == 390),j]
  beta_l_var[j] <- beta_exp[which(rank(beta_exp[,j]) == 10),j]
  beta_u_var[j] <- beta_exp[which(rank(beta_exp[,j]) == 390),j]
}

for(j in 1:(p*p*collectG[[zz]][cc])){
  Delt_l_var[j] <- Delt_exp[which(rank(Delt_exp[,j]) == 10),j]
  Delt_u_var[j] <- Delt_exp[which(rank(Delt_exp[,j]) == 390),j]
}

delta_g <- list()
gamma_g <- list()
mu_g <- list()
beta_g <- list()
Delt_g <- list()
for (grp in 1:collectG[[zz]][cc]){
  delta_g[[grp]] <- delta[grp]
  gamma_g[[grp]] <- gamma[grp]
  mu_g[[grp]] <- mu[(1+p*(grp-1)):(p*grp)]
  beta_g[[grp]] <- beta[(1+p*(grp-1)):(p*grp)]
  k <- matrix(Delta[(1+p*p*(grp-1)):(p*p*grp)],nrow=p,ncol=p)
  Delt_g[[grp]] <- k/det(k)^(1/p)
}
Delta_fin <- unlist(Delt_g)


classAllChains <- Reduce(cbind,class_labs)
class <- numeric(n)
for(i in 1:n){
  class[i] <- as.numeric(names(which.max(table(classAllChains[i,]))))
}

ARI <- mclust::adjustedRandIndex(class,true_lab)
BIC.choose <- as.integer(GlobalMaxBIC)
# BIC.choose <- as.integer(max.avg.BIC)

if(BIC.choose == 2){
  #############################
  ### generate contour plot ###
  #############################
  den <- function(mu, beta, gamma, delta, Delt, lam, y) {
    alpha <- sqrt(gamma^2 + t(beta) %*% Delt %*% beta)
    sq_phi <- sqrt(1 + (1/delta^2) * t(y - mu) %*% solve(Delt) %*% (y - mu))
    
    first <- delta*(1/(2^((p-1)/2)))*exp(gamma*delta+(y-mu)%*%(beta))
    sec <- (alpha/(pi*delta*sq_phi))^(p+1)/2
    third <- besselK(alpha*delta*sq_phi,lam)
    product <- first * sec * third
    return(product)
    # sq_phi <- sqrt(1 + (1/gamma^2) * t(y - mu) %*% solve(Delt) %*% (y - mu))
    # first <- gamma^(1/2) * (1/(2^((p - 1)/2))) * exp(gamma^2 + (y - mu) %*% (beta))
    # sec <- (alpha/(pi * gamma * sq_phi))^(p + 1)/2
    # third <- besselK(alpha * gamma * sq_phi, lam)
    # return(first * sec * third)
  }
  
  z_fun <- function(a, b, mu, beta, gamma, delta, Delt, lam = lambda,pi_g) {
    z1 <- matrix(0, nrow = length(a), ncol = length(b))
    colnames(z1) <- b
    rownames(z1) <- a
    for (i in 1:length(a)) {
      for (j in 1:length(b)) {
        y <- c(a[i], b[j]) #current coordinates
        tmp <- den(mu, beta, gamma, delta, Delt, lam = lambda, y)
        z1[i, j] <- pi_g*tmp #Put density here...
      }
    }
    return(list(z1 = z1))
  }
  
  lambda = -(p+1)/2
  dat <- Y
  a = seq(min(Y[, 1])-1.5*sd(Y[,1]) , max(Y[, 1]) + 1.5*sd(Y[,1]), length.out = 100)
  b = seq(min(Y[, 2])-1.5 *sd(Y[,2]), max(Y[, 2]) + 1.5*sd(Y[,2]), length.out = 100)
  
  zvals1 <- z_fun(a, b, mu = mu_g[[1]], beta = beta_g[[1]], gamma = gamma_g[[1]], delta = delta_g[[1]], Delt = Delt_g[[1]], lam = lambda,pi_g[1])
  zvals2 <- z_fun(a, b, mu = mu_g[[2]], beta = beta_g[[2]], gamma = gamma_g[[2]], delta = delta_g[[2]], Delt = Delt_g[[2]], lam = lambda,pi_g[2])
  
  forclist <- list()
  forclist[[1]] <- zvals1$z1
  forclist[[2]] <- zvals2$z1
  forcall <- Reduce("+", forclist)
  
  summary(c(zvals1$z1))
  summary(c(zvals2$z1))
  
  Predict <- as.factor(class)
  # pch_v <- PredictLabel
  # col_d<- PredictLabel
  # col_d[col_d==1]<-5
  
  ## creat eps plot
  # setEPS()
  # postscript(file=paste("./run/ContorPlot_",seed.no,".eps",sep=""))
  # contour(a, b,forcall,xlab="Y[,1]",ylab="Y[,2]",xlim=c(-10,12),ylim=c(-1,10),
  #         levels = c(1e-7,5e-7,1e-6,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,0.001,0.005,0.01,0.1,0.5,1))
  # points(dat, col = col_d, pch = pch_v)
  # dev.off()
  library(reshape2)
  library(ggplot2)
  mtrx.melt <- melt(forcall)
  
  cp <- ggplot()+
    stat_contour(data = mtrx.melt,aes(x=Var1, y=Var2, z=value), colour = "black",
      breaks=c(1e-7,5e-7,1e-6,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,0.001,0.005,0.01,0.1,0.5,1))+
    xlab(expression(Y[1])) +  ylab(expression(Y[1]))+
    geom_point(data=data.frame(dat),size=3,aes(x=X1,y=X2,color=Predict,shape=Predict))
  cp
  ggsave(file=paste("./ContorPlot_2var_",seed.no,".eps",sep=""),cp,device=cairo_ps,
         fallback_resolution = 600)
}

# fit GMM
require(mclust)
fit_gaussian <- Mclust(Y,modelNames = "VVV")
GMMchoose <- as.integer(fit_gaussian$G)
GMM_ARI <- adjustedRandIndex(fit_gaussian$classification,true_lab)

# fit tMM
require(teigen)
fit_tMix <- teigen(Y, models="UUUU", Gs=1:5, init="hard")
tMixchoose <- as.integer(fit_tMix$G)
tMix_ARI <- adjustedRandIndex(fit_tMix$classification,true_lab)
# pdf(file=paste("./run/tMixcontour_",seed.no,".pdf",sep=""))
# plot(fit_tMix,what='contour')
# pdf(file=paste("./run/tMixPairs_",seed.no,".pdf",sep=""))
# pairs(Y,col=fit_tMix$classification)
# dev.off()

# fit MixGHD
require(MixGHD)
fit_MGHD <- MGHD(data=Y,G=1:5,gpar0=NULL,method="kmeans",scale=FALSE,modelSel="BIC")
GMGHD <- summary(fit_MGHD)$Cluster
MGHDchoose <- as.integer(length(GMGHD))
MixGHD_ARI <- adjustedRandIndex(fit_MGHD@map,true_lab)
# pdf(file=paste("./run/MGHDcontour_",seed.no,".pdf",sep=""))
# contourpl(fit_MGHD)
# pdf(file=paste("./run/MGHDPairs_",seed.no,".pdf",sep=""))
# pairs(Y,col=fit_MGHD@map)
# dev.off()

all_results <- c(gamma,delta,mu,beta,Delta_fin,BIC.choose,ARI,GMMchoose,GMM_ARI,tMixchoose,
                 tMix_ARI,MGHDchoose,MixGHD_ARI)
BIC_results <- c(as.vector(BIC),collectG[[zzB]][[ccB]])
gam_results <- c(gam_l_var,gam_u_var,gamma)
del_results <- c(del_l_var,del_u_var,delta)
mu_results <- c(mu_l_var,mu_u_var,mu)
beta_results <- c(beta_l_var,beta_u_var,beta)
Delta_results <- c(Delt_l_var,Delt_u_var,Delta)


##############################
### save plots and results ###
##############################
# outputbic.filename = paste("./BIC",".csv",sep='')
# write.table(t(BIC_results),file=outputbic.filename,sep=",",
            # row.names=paste("simulation",seed.no),col.names=F,append = TRUE)
# output.filename = paste("./simulationOutput",".csv",sep='')
# write.table(t(all_results),file=output.filename,sep=",",
            # row.names=paste("simulation",seed.no),col.names=F,append = TRUE)
# if(BIC.choose == 2){
#   # get ready for credible intervals
#   outputbic.filename = paste("./GammaCredibleInterval",".csv",sep='')
#   write.table(t(gam_results),file=outputbic.filename,sep=",",
#               row.names=paste("simulation",seed.no),col.names=F,append = TRUE)
#   outputbic.filename = paste("./DeltaCredibleInterval",".csv",sep='')
#   write.table(t(del_results),file=outputbic.filename,sep=",",
#               row.names=paste("simulation",seed.no),col.names=F,append = TRUE)
#   outputbic.filename = paste("./MuCredibleInterval",".csv",sep='')
#   write.table(t(mu_results),file=outputbic.filename,sep=",",
#               row.names=paste("simulation",seed.no),col.names=F,append = TRUE)
#   outputbic.filename = paste("./BetaCredibleInterval",".csv",sep='')
#   write.table(t(beta_results),file=outputbic.filename,sep=",",
#               row.names=paste("simulation",seed.no),col.names=F,append = TRUE)
#   outputbic.filename = paste("./BigDeltaCredibleInterval",".csv",sep='')
#   write.table(t(Delta_results),file=outputbic.filename,sep=",",
#               row.names=paste("simulation",seed.no),col.names=F,append = TRUE)
# }


# pdf(file=paste("./run/Yplot_",seed.no,".pdf",sep=""))
# pairs(Y, col = true_lab)
# dev.off()

########################
### pairwise plots #####
########################
library(GGally)
Predict <- as.factor(class)
True <- as.factor(true_lab)
dataY <- data.frame(Y)
colnames(dataY)<-c("Y1","Y2")

## Y plot with true label
p_t <- ggplot()+
  geom_point(data=dataY,size=3,aes(x=Y1,y=Y2,color=True))+
  xlab(expression(Y[1]))+ylab(expression(Y[2]))
p_t
ggsave(file=paste("./YPlot_2var_",seed.no,".eps",sep=""),p_t,device=cairo_ps,
       fallback_resolution = 600)
## Pairwise plot with true label
# pp_t<-ggpairs(dataY,upper=list(continuous="points",combo="facethist",
#                                discrete="facetbar",na="na"),aes(colour=True,alpha=0.4),
#               legend=grab_legend(p_t))
# pp_t
# ggsave(file=paste("./YPairwise_",seed.no,".eps",sep=""),pp_t,device=cairo_ps,
#        fallback_resolution = 600)

## Pairwise plot wt predicted label
# pp<-ggpairs(dataY,upper=list(continuous="points",combo="facethist",
#             discrete="facetbar",na="na"),
#             aes(colour=Predict,shape=Predict,alpha=0.4),
#             legend=grab_legend(cp))
# pp
# ggsave(file=paste("./PairwisePlot_",seed.no,".eps",sep=""),pp,device=cairo_ps,
#        fallback_resolution = 600)

image.filename = paste("./Simulation_2var_seed",seed.no,".RData",sep='')
save.image(file = image.filename, version = NULL, safe = TRUE)