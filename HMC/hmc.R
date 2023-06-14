library(MASS)
library(lubridate)
library(readr)
library(dplyr)
library(readxl)
library(tibble)
library(WriteXLS)
library(tidyverse)
library(datos)
library(profvis)
library(timeDate)
library(timeSeries)
library(ggfortify)
library(dlm)
library(coda)
library(LaplacesDemon)
library(fExtremes)
library(extRemes)
library(car)
library(carData)
library(evdbayes)
library(coda)
library(compiler)
library(evir)
library(ggplot2)
library(openxlsx)


DLMG=function(m0,C0,FF,G,V,W,T,n,p,Xt,muestra){
  #m0=0;C0=0.6;FF=1;G=1;V=1;W=1;T=length(Go);n=1;p=1;Xt=t(Go);
  
  at=matrix(ncol=T,nrow=p)
  Rt=array(dim=c(T,p,p))
  ft=matrix(ncol=T,nrow=n)
  Qt=matrix(ncol=T,nrow=n)
  Ct=array(dim=c(T,p,p))
  mt=matrix(ncol=T,nrow=p)
  et=matrix(ncol=T,nrow=n)
  At=array(dim=c(T,p,n))
  
  mts=matrix(ncol=T,nrow=p)
  Cts=array(dim=c(T,p,p))
  Bts=array(dim=c(T,p,p))
  
  #  theta=matrix(ncol=T,nrow=p)
  ##	ECUACIONES ACTUALIZACION PARA t=1, ... ,length(Go).
  ##t=1
  
  M_t=matrix(0,ncol=T,nrow=p)
  theta_t=matrix(0,ncol=T,nrow=p)
  Theta_t_mayus=matrix(0,ncol=T,nrow=p)
  
  h_t=matrix(0,ncol=T,nrow=n)
  H_t=matrix(0,ncol=T,nrow=n)
  #Theta_t_mayus=matrix(0,ncol=T,nrow=p)
  
  
  at[,1] = G%*%m0
  Rt[1,,] = (G%*%C0%*%t(G))+W
  ft[,1] = FF%*%at[,1]
  Qt[,1] = FF%*%Rt[1,,]%*%t(FF) + V
  At[1,,] = Rt[1,,]%*%t(FF)%*%solve(Qt[,1])
  et[,1] = Xt[,1] - ft[,1]  
  mt[,1] = at[,1] + At[1,,]%*%t(et[,1])
  Ct[1,,] = Rt[1,,] - At[1,,]%*%t(Qt[1,1])%*%t(At[1,,])
  
  
  #forward Filter
  for(t in 2:T){
    at[,t] = G%*%mt[,t-1]    
    Rt[t,,] = (G%*%Ct[t-1,,]%*%t(G)) + W
    ft[,t] = FF%*%at[,t]
    Qt[,t] = FF%*%Rt[t,,]%*%t(FF) + V
    At[t,,] = Rt[t,,]%*%t(FF)%*%solve(Qt[,t])
    et[,t] = Xt[,t] - ft[,t]  
    mt[,t] = at[,t] + At[t,,]%*%t(et[,t])
    Ct[t,,] = Rt[t,,] - At[t,,]%*%t(Qt[,t])%*%t(At[t,,])
    
    
    theta_t[,t] = rnorm(1 , mean =at[,t], sd = sqrt(Rt[t,,]))
    M_t[,t] = rnorm(1, mean = ft[,t] , sd = sqrt(Qt[,t]))
    Theta_t_mayus[,t] = rnorm(1 , mean = M_t[,t] , sd = sqrt(Ct[t,,]))
    
    H_t[,t] = solve(solve(Ct[t,,])+t(G)%*%solve(W)%*%G)
    h_t[,t] = H_t[,t]%*%(solve(Ct[t,,])%*%mt[,t] + t(G)%*%solve(W)%*%M_t[,t])
  }
  #BACKWAR SAMPLIN
  if(muestra==1){
    mus = rep(0,T)
    mus[T] = rnorm(1,at[,T],sqrt(Ct[T,,]))
    for (t in (T-1):1) {
      mus[t]=rnorm(1,h_t[,(t+1)],sqrt(H_t[,(t+1)]))
    }
  }else{
    #mus = matrix(0,muestra,T)
    mus=array(dim=c(T,muestra,p))
    mus[T,,]=rnorm(muestra,at[,T],sqrt(Ct[T,,]))
    for (t in (T-1):1) {
      mus[t,,]=rnorm(muestra,h_t[,(t+1)],sqrt(H_t[,(t+1)]))
    }
  }
  #list(mus=mus,M_t=M_t)  
  return(mus)
  ##list(theta_t = theta_t , M_t = M_t , Theta_t_mayus = Theta_t_mayus)
  ##list(ft=ft,Qt=Qt,mt=mt,Ct=Ct,At=At,et=et,at=at,Rt=Rt)
} 



require(compiler) # speed up R
enableJIT(3)


##################
J <- -digamma(1) 

# Densidad GEV
########################
f <- function(x, mu, sig, xi) {
  z = 1 + xi*(x-mu)/sig
  ifelse(z > 0, z^(-1/xi-1)/sig * exp(-z^(-1/xi)), 0)
}

#Distribucion GEV
#############################
# F <- function(x, mu, sig, xi) {
#      z = 1 + xi*(x-mu)/sig
#      ifelse(z > 0, exp(-z^(-1/xi)), 0)
# }

# Valores aleatorios GEV
#######################################################################
rgev <- function(n, mu, sig, xi) {
  u <- runif(n)
  x <- ( mu - (sig/xi) * (1-(-log(u))^(-xi)) ); x
}

# y <- rgev(25, 2, .5, -.5); hist(y, prob = T)
# n <- length(y)

#######################################################################
### log-posterior para GEV
#######################################################################

ll <- function(theta, x) {
  mu    = theta[1]
  delta = theta[2]
  xi    = theta[3]
  
  ll <- sum(log(f(x, mu, exp(delta), xi))) + sum(dnorm(theta, mean = 0, sd = 5, log = T))
  ifelse(ll != "-Inf", ll, -100000)
}


### gradiente log-posterior para GEV
######################################################################

nll <- function(theta, x) {
  D <- rep(0, 3)
  mu    = theta[1]
  delta = theta[2]
  xi    = theta[3]
  n     = length(x)
  
  sig = exp(delta)
  d   = x - mu
  z   = 1 + xi*d/sig
  w   = d/sig^2
  q   = d/sig
  
  D[1] = 1/sig*sum(z^(-1)*((1+xi) - z^(-1/xi))) - mu/25
  D[2] = (sum((1+xi)*w*z^(-1) - w*z^(-1-1/xi))*sig - n) - delta/25
  D[3] = suppressWarnings(1/xi^(2)*sum(log(z)) - (1/xi+1)*sum(q*z^(-1)) + 
                            1/xi*sum(q*z^(-1-1/xi)) - 1/xi^(2)*sum(log(z)*z^(-1/xi)) - xi/25)
  
  if (suppressWarnings(any(z < 0 | D == "NaN" | D == "Inf" | D == "-Inf"))) rep(-100000, 3)
  else D
}




# Matriz de información esperada de Fisher y sus derivadas
################################################################################

G <- function(theta, x) {
  G     = matrix(, 3, 3)
  mu    = theta[1]
  delta = theta[2]
  sig   = exp(theta[2])
  xi    = theta[3]
  
  n = length(x)
  
  pg <- gamma(2 + xi)
  p  <- (1 + xi)^(2) * gamma(1 + 2*xi)
  q  <- pg * (psigamma(1 + xi) + 1 + 1/xi)
  
  G[1, 1] =  p/sig^2 + 1/25
  G[1, 2] = -1/(sig*xi)*(p - pg)
  G[1, 3] = -1/(sig*xi)*(q - p/xi)
  G[2, 1] = G[1, 2]
  G[2, 2] =  1/(xi^2) * (1 - 2*pg + p) + 1/25
  G[2, 3] = -1/(xi^2) * (1 - J + (1 - pg)/xi - q + p/xi)
  G[3, 1] = G[1, 3]
  G[3, 2] = G[2, 3]
  G[3, 3] = 1/(xi^2) * ((pi^2)/6 + (1 - J + 1/xi)^2 - 2*q/xi + p/(xi^2)) + 1/25
  
  return(n*G)
} 




leapfrog_step = function(theta, r, eps, M_diag , x){
  r_tilde <- r + 0.5 * eps * nll(theta , x)
  theta_tilde <- theta + eps * r_tilde / M_diag
  r_tilde <- r_tilde + 0.5 * eps * nll(theta , x)
  list(theta = theta_tilde, r = r_tilde)
}




joint_log_density = function(theta, r, M_diag , x){
  ll(theta , x) - 0.5*sum(r**2 / M_diag)
}



find_reasonable_epsilon = function(theta, M_diag, eps = 1, verbose = TRUE, x){
  r <- rnorm(length(theta), 0, sqrt(M_diag))
  
  proposed <- leapfrog_step(theta, r, eps, M_diag , x)
  
  log_ratio <- joint_log_density(proposed$theta, proposed$r, M_diag , x) - 
    joint_log_density(theta, r, M_diag , x)
  
  alpha <- ifelse(exp(log_ratio) > 0.5, 1, -1)
  if(is.nan(alpha)) alpha <- -1
  count <- 1
  while(is.nan(log_ratio) || alpha * log_ratio > (-alpha)*log(2)){
    eps <- 2**alpha * eps
    proposed <- leapfrog_step(theta, r, eps, M_diag,x)
    
    log_ratio <- joint_log_density(proposed$theta, proposed$r, M_diag , x) -
      joint_log_density(theta, r, M_diag , x)
    
    count <- count + 1
    if(count > 100) {
      stop("no se encontro epsilon en 100 iteraciones")
    }
  }
  if(verbose)
    #message("epsilon = ", eps, " se encontro despues de ", count, " pasos")
    eps
  #return(eps)
}





# Hamiltonian Monte Carlo

hmc <- function(theta.current, SS, burn, lag, epsilon, LF, M, y)  {
  
  ratio = 0
  
  SS. = (SS + burn) * lag
  idx = seq(burn * lag + 1, SS., by = lag)
  
  ## Matrices para las cadenas
  D.         <- length(theta.current)
  mD.        <- rep(0, D.)
  G.         <- M
  invG       <- solve(G.)
  theta      <- matrix( , SS., D.)
  theta[1, ] <- theta.current  
  
  start.time <- Sys.time()
  
  ## generando las cadenas
  for(i in 2:SS.) {
    
    pn     <- mvrnorm(1, mD., G.) ## G.
    thetan <- theta[i-1, ]
    
    ## hamilton corriente
    H.current <- - ll(thetan, y) + t(pn) %*% invG %*% pn/2 ## invG
    
    ## saltos de rana
    pn <- pn + epsilon/2 * nll(thetan, y)
    for (l in 1:LF) {
      thetan <- thetan + epsilon * invG %*% pn  ## invG
      if (l != LF)  pn <- pn + epsilon * nll(thetan, y)
    }
    pn <- pn + epsilon/2 * nll(thetan, y)
    pn <- -pn
    
    ## hamilton propuetso y hamilton logaritmico
    H.prop <- - ll(thetan, y) + t(pn) %*% invG %*% pn/2 ## invG
    logratio <- - log(H.prop) + log(H.current)
    
    ## applying metropolis rule
    if ( (logratio > 0 | logratio > log(runif(1))) & logratio != "NaN" ) {
      theta[i, ] <- thetan
      ratio  <- ratio + 1
    }
    else theta[i, ] <- theta[i-1, ]
    
  }
  print(Sys.time()-start.time)
  theta = theta[idx, ] 

  
  return(list(theta = theta, r = ratio/SS))
}






fgfg <- read_xlsx("C:/Users/Elvis/Documents/SU_100/TODAS_1000.xlsx")
f1 <- nrow(fgfg) 
c1 <- ncol(fgfg)



for (i in 1:1) {
  temp_1 = fgfg[1,i]
  j = 2
  contador = 2
  #print(i)
  while( (j != 0)  && (contador <= f1) ){
    temp_1 = rbind(temp_1, fgfg[contador,i])   
    contador = contador + 1
    j=fgfg[contador,i]
    
  }
  
  
  Go<-temp_1
  tam<-nrow(Go)
  
  landa<-0.6
  
  mu_tep<-mean(DLMG(m0=0,C0=0.6,FF=1,G=1,V=1,W=1,T=tam,n=1,p=1,Xt=t(Go),muestra=1000))
  
  q <- c(mu_tep , 1.3 , 0.4)
  
  epsilons<-find_reasonable_epsilon( q , 1 , eps=1,verbose = TRUE, as.numeric(unlist(Go)))
  
  
  LF <- max(1, round((0.9+runif(1)/5)*landa/epsilons))
  
  
  n=10000  
  
  BB = hmc(q, n, 10, 1, epsilons, LF , diag(3) , as.numeric(unlist(Go)))
  
  
}


