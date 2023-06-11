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
sampling <- function(y,chains = 4, iter = 5000, scale = 0.5, thin = 0,
                     warm_up = round(iter/2,0), Hastings = TRUE,
                     MALA = FALSE, h = 1,seed = NULL) {
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  post = NULL
  init = inits()
  d = length(init)
  
  if(length(diag(scale)) != d)
    stop("The dimension of the scale matrix and parameters don't match")
  
  scale = sqrt(h)*chol(scale)
  
  for(k in 1:chains){
    results = NULL
    current_state = inits()
    # burn-in
    for(i in 1:warm_up) {
      out = step(y, current_state, h, scale, Hastings)
      current_state = out$value
    }
    current_state = current_state
    # sample
    for(i in 1:iter) {
      for(j in 1:thin) {
        out = step(y, current_state, h, scale, Hastings, MALA)
        current_state = out$value
      }
      results = rbind(results, c(out$value, out$accepted))
    }
   post = rbind(post,results) 
  }
  row.names(post) = NULL
  post = as.data.frame(post)
  post$.chain = sort(rep(1:chains, iter))
  
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
step <- function(y, prev, h = 1, scale, Hastings = TRUE, MALA = FALSE) {
  
  u  =  runif(1)
  
  ml1 = mala_step(x = prev, h = h, scale = scale, MALA = MALA)
  
  prop = rjump(prev + ml1, scale)
  
  ml2 = mala_step(x = prop, h = h, scale = scale, MALA = MALA)
  
  t1 = target(y, prop) + Hastings*djump(prev, prop + ml1, scale)
  t2 = target(y, prev) + Hastings*djump(prop, prev + ml2, scale)
  accept_prob = exp(t1 - t2)
  accept_prob = ifelse(is.na(accept_prob),0,accept_prob)
  
  if(u <= accept_prob){
    value = prop
    accepted = TRUE
  }else{
    value = prev
    accepted = FALSE
  }
  
  out = list(value = value, accepted = accepted)
  return(out)
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
rjump <- function(loc, scale){
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
#' Evaluates the target function in a metropolis step
#' for a constant data set y
#' 
#'  target(y,theta) = log_lik(y,theta) + log_prior(theta)
#'  
#' @param x a a vector of the unknown parameters to be evaluated in 
#' the prior density.
#' 
#' @author Asael Alonzo Matamoros
#' 
#' @return a real value with the image of the target distribution
#' 
t_prop = function(x) {
  target(y,x)
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

#' MALA step 
#' 
#' Performs a MAetropolis Adjusted Lavenging algorithm steo
#' 
#' @param x vector with the current metropolis step.
#' @param h scale correction for a scale adjusted Metropolis.
#' @param scale scale matrix for the proposed Jump distribution
#' @param mala Boolean that indicates a MALA step for the Jump candidate in
#' the Metropolis step.
#' 
#' This function generates a  MCMC step using one of the three algorithms:
#'  - if MALA == TRUE, returns thr h GRAD (log f(x))/2 MALA factor
#'  - else returns 0
#' 
#' @author Asael Alonzo Matamoros
#' 
#' @return a vector with the proposed MALA step
#' 
mala_step <-function(x, h = 1, scale, MALA = FALSE){
  d = length(x)
  r = rep(0,d)
  if(MALA){
    r = try(h*scale%*%numDeriv::grad(func = t_prop,x)/2,silent = TRUE)
    if(inherits(r, "try-error")){
      r = rep(0,d)
    }else{
      r = as.numeric(r) 
    }
  }
  return(r)
}

#' Kalman Filtering update equations
#' 
#' Obtain the mean and covariance parameters of the posterior of mu
#' when applying a Vanilla Kalman Filter for a constan DLM model
#' 
#'   yt = FF mu_t + et, et ~ N(0,W)
#'   mu_t = G mu_{t-1} + vt, vt ~ N(0,V)
#'   y0 ~ N(m0,C0)
#'   
#' The function computes the poserior of the states parameters
#' p(mu_t | y1,y2,...,yt)
#' 
#' @param y a matrix of dimensions of k columns and n rows containing 
#' a sample of  `n` elements `y` of dimension `m`.
#' @param m0 is the mean vector [k] with mean for the prior of the y0.
#' @param C0 is the covariance matrix [k,k] with the the initial value y0
#' @param FF is the link matrix between states (mu_t) and observations (yt)
#' @param G the transition matrix of dimension `kxk` for the states equation.
#' @param samples an integer with the amount of samples for approximating the
#' latent distribution of mu_t | y_t.
#' 
#' @author Asael Alonzo Matamoros
#' 
#' @return Returns a list with parameters of the posterior of mu_t at every step
#' 
#' the list lt = list(mt, Ct), where: 
#'  
#'   - mt is a data frame of n row and k columns.
#'   - Ct is an array of matrix with dimensions c(n,k,k)
#' 
#' Such that for every time-step the posterior of mu_t is:
#'    
#'    mu_t | y1,y2,...,yt ~ N(mt, C_t)
#' 
forward_filter<-function(y, G = 1, FF = 1, V = 1, W = 1, m0 = 0, C0 = 1){
  
  # filter dimensions
  if(is.matrix(y)){
    n = nrow(y)
    m = ncol(y)
    
    if(any(dim(V) != m))
      stop("Formato invalido para W, debe ser una matriz de m x m")
    
  }else{
    n = length(y)
    m = 1
    y = matrix(y,ncol = 1)
    
    if(!is.numeric(V) || length(V) != 1)
      stop("Formato invalido para V, debe ser un numero real positivo")
  }
  
  if(is.matrix(G)){
    k = ncol(G)
    
    if(any(length(m0) != k))
      stop("Formato invalido para m0, debe ser un vector de tamanio k")
    
    if(any(dim(C0) != k))
      stop("Formato invalido para C0, debe ser una matriz de k x k")
    
    if(any(dim(FF) != c(k,m) ))
      stop("Formato invalido para FF, debe ser una matriz de k x m")
    
    if(any(dim(W) != k))
      stop("Formato invalido para W, debe ser una matriz de k x k")
    
  }else{
    if(is.numeric(G) && length(G) == 1){
      k = 1
      
      if(!is.numeric(m0) || length(m0) != 1)
        stop("Formato invalido para m0, debe ser un numero real")
      
      if(!is.numeric(C0) || length(C0) != 1)
        stop("Formato invalido para C0, debe ser un numero real positivo")
      
      if(!is.numeric(FF) || length(FF) != 1)
        stop("Formato invalido para FF, debe ser un numero real")
      
      if(!is.numeric(W) || length(W) != 1)
        stop("Formato invalido para W, debe ser un numero real positivo")
      
    }else{
      stop("Formato invalido para G, debe ser una matriz o un numero real") 
    }
  }
  
  # initial values
  mt = m_temp = m0
  C_temp = C0
  Ct = list(C0)
  
  for(t in 1:n){
    # prior at t
    a = as.numeric(G%*%m_temp)
    R = (G%*%C_temp %*%t(G)) + W
    
    # Marginal likelihood: 
    f = as.numeric( t(FF) %*% a)
    Q = (t(FF) %*%R%*%FF) + V
    
    # additional computations at t
    e = as.numeric(y[t,] - f)
    A = R %*% FF %*% solve(Q)
    
    #  Posterior at t
    m_temp = a + as.numeric(A%*%e)
    C_temp = R - (A %*%Q%*%t(A))
    
    # Append othe historical data
    mt = rbind(mt,m_temp)
    row.names(mt) = NULL
    
    Ct[[t+1]] = C_temp
  }
  pp = list(mt = mt, Ct = Ct)
  return(pp)
}


