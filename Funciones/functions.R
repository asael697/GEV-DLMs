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
#' @param p is the dimension of the state parameters mu_t in R^p
#' @param data a matrix of dimensions m rows and n cols containing the data set.
#' @param samples an integer with the amount of samples for approximating the
#' latent distribution of mu_t | y_t.
#' 
#' @return an array mu of dimension n X samples x p

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

