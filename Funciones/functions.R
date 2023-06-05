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
#' Generates different MARKOV Chains with the simulated marginal posteriors
#' 
#' @param y a matrx with the used data for estimating the posterior
#' @param chains integer with the total amount of chains 
#' @param iter an integer with the total iterations per chain after warm up
#' @param scale scale matrix for the proposed Jump distribution
#' @param warm_up an integer for the first iterations to be burned in warm-up
#' iterations.
#' @param thin the amount of lags to be burned for avoid stuck jumps
#' @param Hastings performs a Metropolis-Hastings step with an asymmetric Jump
#' function.
#' @param dif_adjust performs a differential adjustment to perform a Jump
#' candidate in the Metropolis step.
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
                               dif_adjust = FALSE,mala = FALSE) {
  
  if(dif_adjust) {
    mala = FALSE
    Hastings = TRUE
  }
  
  if(mala) {
    dif_adjust = FALSE
    Hastings = TRUE
    ds = diag(scale)
    tau = mean(ds)
    scale1 = abs(tau)*diag(length(ds))
  }
  
  d1 = length(inits())
  if(length(diag(scale)) != d1)
    stop("The dimension of the scale matrix and parameters don't match")
  
  
  post = NULL
  for(k in 1:chains){
    results = NULL
    current_state = inits()
    # burn-in
    for(i in 1:warm_up) {
      out = step(y,current_state,scale,Hastings,dif_adjust,mala)
      current_state = out$value
    }
    # sample
    for(i in 1:iter) {
      for(j in 1:thin) {
        out = step(y,current_state,scale,Hastings,dif_adjust,mala)
        current_state = out$value
      }
      results = rbind(results,c(out$value, out$reject))
    }
   post = rbind(post,results) 
  }
  row.names(post) = NULL
  post = as.data.frame(post)
  post$.chain = sort(rep(1:chains,iter))
  
  return(post)
}

#' Metropolis-Hastings Steps
#' 
#' Generates a Jump for the Metropolis Algorithm 
#' 
#' @param y a matrx with the used data for estimating the posterior.
#' @param prop a vector with the previous step
#' @param scale scale matrix for the proposed Jump distribution
#' @param Hastings performs a Metropolis-Hastings step with an asymmetric Jump
#' function.
#' @param dif_adjust performs a differential adjustment to perform a Jump
#' candidate in the Metropolis step.
#' 
#' This function use a symmetric multivariate normal distribution
#' to generates the Jumps, in case of an asymmetric jump, use the
#' the Metropolis-Hastings jump.
#' 
#' @author Asael Alonzo Matamoros
#' 
#' @return a vector with the proposed jump
#' 
step <- function(y,prop,scale,Hastings = TRUE,
                 dif_adjust = FALSE,mala = FALSE) {
  reject = 1
  ls = diff_adjustment(prop, scale, dif_adjust,mala)
  proposed = rjump(prop = ls$prop, scale = ls$scale)

  u  =  runif(1)
  t1 = target(y, proposed) + Hastings*djump(prop, ls$prop, ls$scale) 
  t2 = target(y, prop) + Hastings*djump(proposed, ls$prop, ls$scale) 
  et = exp(t1 - t2)
  et = ifelse(is.na(et),0,et)
  accept_prob = min(1,et)
  
  if(u <= accept_prob){
    value = proposed
    reject = 0
  }else{
    value = prop
  }
  
  out = list(value = value,reject = reject)
  return(out)
}

#' Generate a proposed value from a Jump distribution
#' 
#' @param prop a vector with the proposed parameters from previous
#' iteration.
#' @param scale a covariance matrix used in the Gaussian Jump distribution
#' 
#' @author Asael Alonzo Matamoros
#' 
#' @return a vector following a Gaussian distribution of dimension d,
#' with mean = prop, and covariance matrix C = scale.
#' 
rjump <- function(prop,scale){
  d = length(prop)
  x = as.numeric(scale%*%rnorm(n = d))
  return(x + prop)
}

#' Evaluate the density of a Gaussian distribution used as Jump density
#' 
#' @param proposed a vector to evaluate the multivariate density.
#' @param prop a vector with the proposed parameters from previous
#' iteration.
#' @param scale a covariance matrix used in the Gaussian Jump distribution
#' 
#' @author Asael Alonzo Matamoros
#' 
#' @return a vector following a Gaussian distribution of dimension d,
#' with mean = prop, and covariance matrix C = scale.
#' 
djump <- function(proposed,prop,scale){
  d = length(proposed)
  
  z = solve(scale)%*%(proposed - prop)
  d2 = dnorm(z,log = TRUE)
  
  return(sum(d2))
}
#' Estimate the mean and covariance using a differential adjustment
#' 
#' The values are obtained by applying:
#'   
#'   opt = arg_max(target)
#'   mu  = opt + Sigma* GRAD(target(opt))
#'   Sigma^{-1} = -Hessian(target(opt))
#' 
#' @param prop a vector with the proposed parameters from previous
#' iteration.
#' @param scale a covariance matrix used in the Gaussian Jump distribution
#' @param dif_adjust performs a differential adjustment to perform a Jump
#' candidate in the Metropolis step.
#' 
#' @return if `dif_adjut` is `TRUE` returns a list with the adjusted mean 
#' and covariance, if `FALSE` returns a list with the original values.
#' 
diff_adjustment <- function(prop,scale,dif_adjust = FALSE, mala = FALSE){
  t_prop = function(par) target(y,par)
  
  if(dif_adjust){
    opt = optim(prop,fn =  function(par) -target(y,par))$par
    
    scale = chol2inv(
      -numDeriv::hessian(t_prop,x = opt)
    )
    
    prop = scale %*% numDeriv::grad(t_prop,x = opt)
    prop = opt + as.numeric(prop)  
  }
  if(mala){
    prop = prop + as.numeric(scale %*% numDeriv::grad(t_prop,x = prop))
    scale = 2*scale
  }
  return(list(prop = prop,scale = chol(scale) ))
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
