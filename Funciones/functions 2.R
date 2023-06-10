#' Forward Filtering Backward Sampling algorithm
#' simulates the states parameters in Dynamic linear model
#' 
#'   yt = FF mu_t + et, et ~ N(0,W)
#'   mu_t = G mu_{t-1} + vt, vt ~ N(0,V)
#'   y0 ~ N(m0,C0)
#' The function computes the distribution of the states parameters
#' p(mu_t | y1,y2,...,yt)
#' 
#' @param m0 is the mean of the initial value y0
#' @param C0 is the covariance matrix of the initial value y0
#' @param FF is the link matrix between states (mu_t) and observations (yt)
#' @param G the transition matrix in the states equation.
#' @param n is the length of the data set n = length(y_t)
#' @param m is the dimension of the data yt in R^m
#' @param p is the dimension of the state parameters mu_t in R^p.
#' @param data a matrix of dimensions m rows and n cols containing the data set.
#' @param samples an integer with the amount of samples for approximating the
#' latent distribution of mu_t | y_t.
#' 
#' @author Elvis Arrazola y Cristian Cruz
#' 
#' @return an array mu of dimension n X samples x p
#' 
FFBS <- function(m0 = 0, C0 = 0.6, FF = 1, G = 1, V = 1, 
                 W = 1, n, m = 1, p = 1, data, samples){
  
  at <- mt <- mts <- Mt <-  theta_t <- matrix(ncol = n,nrow = p)
  Rt <- Ct <- Cts <- Bts <- array(dim = c(n,p,p))
  ft <- et <-h_t <- H_t <- matrix(ncol = n,nrow = m)
  Qt <- matrix(ncol = n,nrow = m)
  At <- array(dim = c(n,p,n))
  Theta_t_mayus <- theta_t <- matrix(0,ncol = n,nrow = p)
  
  ##	ECUACIONES ACTUALIZACION PARA t = 1, ... ,length(Go).
  at[,1] = G%*%m0
  Rt[1,,] = (G%*%C0%*%t(G))+W
  ft[,1] = FF%*%at[,1]
  Qt[,1] = FF%*%Rt[1,,]%*%t(FF) + V
  At[1,,] = Rt[1,,]%*%t(FF)%*%solve(Qt[,1])
  et[,1] = data[,1] - ft[,1]  
  mt[,1] = at[,1] + At[1,,]%*%t(et[,1])
  Ct[1,,] = Rt[1,,] - At[1,,]%*%t(Qt[1,1])%*%t(At[1,,])
  
  
  #forward Filter
  for(t in 2:n){
    at[,t] = G%*%mt[,t-1]    
    Rt[t,,] = (G%*%Ct[t-1,,]%*%t(G)) + W
    ft[,t] = FF%*%at[,t]
    Qt[,t] = FF%*%Rt[t,,]%*%t(FF) + V
    At[t,,] = Rt[t,,]%*%t(FF)%*%solve(Qt[,t])
    et[,t] = data[,t] - ft[,t]  
    mt[,t] = at[,t] + At[t,,]%*%t(et[,t])
    Ct[t,,] = Rt[t,,] - At[t,,]%*%t(Qt[,t])%*%t(At[t,,])
    
    
    theta_t[,t] = rnorm(1 , mean =at[,t], sd = sqrt(Rt[t,,]))
    M_t[,t] = rnorm(1, mean = ft[,t] , sd = sqrt(Qt[,t]))
    Theta_t_mayus[,t] = rnorm(1 , mean = M_t[,t] , sd = sqrt(Ct[t,,]))
    
    H_t[,t] = solve(solve(Ct[t,,])+t(G)%*%solve(W)%*%G)
    h_t[,t] = H_t[,t]%*%(solve(Ct[t,,])%*%mt[,t] + t(G)%*%solve(W)%*%M_t[,t])
  }
  # BACKWARD SAMPLING
  if(samples == 1){
    mu = rep(0,n)
    mu[n] = rnorm(1,at[,n],sqrt(Ct[n,,]))
    for (t in (n-1):1) {
      mu[t] = rnorm(1,h_t[,(t+1)],sqrt(H_t[,(t+1)]))
    }
  }else{
    #mus = matrix(0,samples,T)
    mu <- array(dim = c(n,samples,p))
    mu[n,,] = rnorm(samples,at[,n],sqrt(Ct[n,,]))
    
    for (t in (n-1):1) {
      mu[t,,] = rnorm(samples,h_t[,(t+1)],sqrt(H_t[,(t+1)]))
    }
  }

  return(mu)
} 

#' Metropolis Sampler algorithm
#' 
#' Generates different Markov Chains with the simulated marginal posteriors
#' 
#' @param y matrix with the used data for estimating the posterior
#' @param chains integer with the total amount of chains 
#' @param iter integer with the total iterations per chain after warm up
#' @param scale scale matrix for the proposed Jump distribution
#' @param warm_up integer for the first iterations to be burned in warm-up
#' iterations.
#' @param thin integer with the amount of lags to be burned for avoid stucked
#' jumps.
#' @param Hastings performs a Metropolis-Hastings step with an asymmetric Jump
#' function.
#' @param diff_adjust Boolean that indicates a differential adjustment For a
#' Jump candidate in the Metropolis step.
#' @param mala Boolean that indicates a MALA step for the Jump candidate in
#' the Metropolis step.
#' 
#' @author Asael Alonzo Matamoros
#' 
#' @return an data frame with the simulated posteriors, a vector `.chain`
#' containing the chain which the sample belongs, and a vector `reject`
#' which stats if the current jump was rejected.
#' 
#' 
sampling <- function(y,chains = 4, iter = 5000, scale = 0.5, thin = 5,
                     warm_up = round(iter/2,0), Hastings = TRUE,
                     diff_adjust = FALSE) {
  
  init = inits()
  d = length(init)
  mu = rep(0,d)
  
  if(length(diag(scale)) != d)
    stop("The dimension of the scale matrix and parameters don't match")
  
  post = NULL
  
  for(k in 1:chains){
    results = NULL
    current_state = inits() 
    # burn-in
    for(i in 1:warm_up) {
      out = step(y,current_state+mu,scale,Hastings)
      current_state = out$value
    }
    # sample
    for(i in 1:iter) {
      for(j in 1:thin) {
        out = step(y,current_state+mu,scale,Hastings)
        current_state = out$value
      }
      results = rbind(results, c(out$value,out$accepted))
    }
   post = rbind(post,results) 
  }
  row.names(post) = NULL
  post = as.data.frame(post)
  post$.chain = sort(rep(1:chains,iter))
  
  return(post)
}

#' MCMC Steps
#' 
#' Generates a Jump for a MCMC algorithm 
#' 
#' @param y matrix with the used data for estimating the posterior.
#' @param prop vector with the previous step.
#' @param scale scale matrix for the proposed Jump distribution.
#' @param Hastings Boolean which performs a Metropolis-Hastings step with an
#' asymmetric Jump function.
#' @param diff_adjust Boolean that indicates a differential adjustment For a
#' Jump candidate in the Metropolis step.
#' @param mala Boolean that indicates a MALA step for the Jump candidate in
#' the Metropolis step.
#' 
#' This function generates a  MCMC step using one of the three algorithms:
#'  - Metropolis Hastings
#'  - Differential Adjusted Metropolis
#'  - Metropolis-adjusted Langevin Algorithm
#' 
#' @author Asael Alonzo Matamoros
#' 
#' @return a vector with the proposed jump
#' 
step <- function(y, prev, scale, Hastings = TRUE) {
  
  u  =  runif(1)
  lstep = metropolis(y, prev, scale, Hastings)
  
  if(u <= lstep$ap){
    value = lstep$proposed
    accepted = TRUE
  }else{
    value = prev
    accepted = FALSE
  }
  
  out = list(value = value, accepted = accepted)
  return(out)
}
#' Metropolis-Hastings Steps
#' 
#' Generates the acceptance probability for a Metropolis step 
#' 
#' @param y matrix with the used data for estimating the posterior.
#' @param prev vector with the previous step.
#' @param scale scale matrix for the proposed Jump distribution.
#' @param Hastings Boolean which performs a Metropolis-Hastings step with an
#' 
#' This function use a symmetric multivariate normal distribution
#' to generates the Jumps, in case of a predefined asymmetric jump, use the
#' the Metropolis-Hastings jump.
#' 
#' @author Asael Alonzo Matamoros
#' 
#' @return a list with a real value with the acceptance probability,
#' and a vector with the proposed value
#' 
metropolis <- function(y, prev, scale, Hastings = TRUE){
  
  proposed = as.numeric(mvtnorm::rmvnorm(1,mean = prev,sigma = scale))
  
  t1 = target(y, proposed) + mvtnorm::dmvnorm(prev, proposed, scale) 
  t2 = target(y, prev) + mvtnorm::dmvnorm(proposed, prev, scale) 
  et = exp(t1 - t2)
  et = ifelse(is.na(et),0,et)
  
  lstep = list(ap = min(1,et),proposed = proposed)
  
  return(lstep)
}

#' Differential adjusted Metropolis algorithm
#' 
#' Generates the acceptance probability for a differential adjusted
#' Metropolis step. 
#' 
#' @param y matrix with the used data for estimating the posterior.
#' @param prev vector with the previous step.
#' @param scale scale matrix for the proposed Jump distribution.
#' 
#' This function use a symmetric multivariate normal distribution
#' to generates the Jumps, where the mean and covariance are:
#' 
#'   - opt = arg_max(target)
#'   - loc  = opt + Sigma* GRAD(target(opt))
#'   - scale^{-1} = -Hessian(target(opt))
#'   - target(y,par) = -loglik(y | par) - log_prior(par)
#' 
#' where GRAD is the gradient operator, pt is the local optimum 
#' of the target distribution at the current step.
#' 
#' @author Asael Alonzo Matamoros
#' 
#' @return a real value with the acceptance probability
#' 
adjusted_metropolis <- function(y, prev, mu, scale){
  
  proposed = as.numeric(mvtnorm::rmvnorm(1,mean = prev+mu,sigma = scale))
  
  t1 = target(y, proposed) + mvtnorm::dmvnorm(prev + mu, proposed, scale) 
  t2 = target(y, prev) + mvtnorm::dmvnorm(proposed, prev + mu, scale) 
  et = exp(t1 - t2)
  (et = ifelse(is.na(et) || is.infinite(et),0,et))
  
  lstep = list(ap = min(1,et),proposed = proposed)
  
  return(lstep)
}

#' Generate a proposed value from a Jump distribution
#' 
#' @param loc a vector with the proposed parameters from previous
#' iteration.
#' @param scale a covariance matrix used in the Gaussian Jump distribution
#' 
#' @author Asael Alonzo Matamoros
#' 
#' @return a vector following a Gaussian distribution of dimension d,
#' with mean = prop, and covariance matrix C = scale.
#' 
rjump <- function(loc,scale){
  d = length(loc)
  x = as.numeric(scale%*%rnorm(n = d))
  return(x + loc)
}

#' Evaluate the density of a Gaussian distribution used as Jump density
#' 
#' @param proposed a vector to evaluate the multivariate density.
#' @param loc a vector with the proposed parameters from previous
#' iteration.
#' @param scale a covariance matrix used in the Gaussian Jump distribution
#' 
#' @author Asael Alonzo Matamoros
#' 
#' @return a vector following a Gaussian distribution of dimension d,
#' with mean = prop, and covariance matrix C = scale.
#' 
djump <- function(proposed, loc, scale){

  z = solve(scale)%*%(proposed - loc)
  d2 = dnorm(z,log = TRUE)
  
  return(sum(d2))
}

#' Evaluates the target function in a metropolis step
#' 
#'  target(y,theta) = log_lik(y,theta) + log_prior(theta)
#'  
#' 
#' @param y a matrx with the used data for estimating the posterior.
#' @param theta a vector of the unknown parameters to be evaluated in 
#' the prior density.
#' 
#' @author Asael Alonzo Matamoros
#' 
#' @return a real value with the image of the target distribution
#' 
target <- function(y,theta){
  
  t = loglik(y,theta) + log_prior(theta)
  return(t)
}

#' Computes the Log_likelihood matrix 
#' 
#' @param y a matrx with the used data for estimating the posterior.
#' @param post a matrix with the posterior parameters, with ncols = d,
#' and nrows = iter.
#' 
#' @author Asael Alonzo Matamoros
#' 
#' @return a matrix of dimension ncols n = length(y) and nrows = iter,
#' with the evaluated log_likelihood.
#' 
log_lik <- function(y,post){
  
  LL = apply(post, 1, function(z){
    z = as.numeric(z)
    gev_lpdf(y,mu = z[1],sigma = z[2],k = z[3])
  })
  t(LL)
}
#' Defaults log_prior density
#' 
#' @param theta a matrix with the posterior parameters, with ncols = d,
#' and nrows = iter.
#' 
#' @author Asael Alonzo Matamoros
#' 
#' @return a real value with the image of the log_prior density
#' 
log_prior <- function(theta){
  return(0)
}
