library(MASS)
library(mcmc)
library(MCMCpack)
library(extRemes)
library(MCMCglmm)
library(nsRFA)




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




enableJIT(3)


##################
### constante euler
##################
J <- -digamma(1) 

########################
### GEV densidad
########################
f <- function(x, mu, sig, xi) {
  z = 1 + xi*(x-mu)/sig
  ifelse(z > 0, z^(-1/xi-1)/sig * exp(-z^(-1/xi)), 0)
}

#############################
### GEV Distribucio
#############################
F <- function(x, mu, sig, xi) {
  z = 1 + xi*(x-mu)/sig
  ifelse(z > 0, exp(-z^(-1/xi)), 0)
}

#######################################################################
### GEV numeros aleatorios
#######################################################################
rgev <- function(n, mu, sig, xi) {
  u <- runif(n)
  x <- ( mu - (sig/xi) * (1-(-log(u))^(-xi)) );
  x
}

# y <- rgev(25, 2, .5, -.5); hist(y, prob = T)
# n <- length(y)


ll2 <- function(theta) {
  mu    = theta[1]
  delta = theta[2]
  xi    = theta[3]
  
  ll <- sum(log(f(xx, mu, exp(delta), xi))) + sum(dnorm(theta, mean = 0, sd = 2, log = T))
  ifelse(ll != "-Inf", ll, -100000)
  
  
  #return(exp(ll))  
}

MCMC01 <- function (x, N=1000, theta0=c(1,0,-.5), pseudo_var=c(1,1,1), burnin=100) {
  #N = tamaño final de la muestra (es decir, excluyendo la duración del rodaje)
  #theta0 = punto de partida de su cadena Metropolis que contiene (mu0, log(sigma0), xi0)
  #pseudo_var = varianza de la normal que se usa como distribución propuesta para el paseo aleatorio
  #Metrópolis (muestreo independiente)
  #burnin = el número especificado será el número de muestras iniciales arrojadas
  
  thetas <- theta0
  
  for (i in 2:(burnin+N)) {
    #loglikelihood0 <- sum(log(dGEV(x, mu=theta0[1],
    #                               sigma=exp(theta0[2]), xi=theta0[3])))
    loglikelihood0<- ll2(theta0)
    #loglikelihood0<-sum(log(f(x,mu=theta0[1],sig=exp(theta0[2]), xi=theta0[3])))
    logprior0 <- 0 # because 1/sigma corresponds to uniform distr of the log(sigma)
    logtarget0 <- loglikelihood0 + logprior0
    
    if(is.nan(logtarget0)) logtarget0 <- -10000000
    prop <- mvrnorm(n=1, mu=theta0, Sigma=diag(pseudo_var))
    #loglikelihood1 <- sum(log(dGEV(x, mu=prop[1],
    #sigma=exp(prop[2]), xi=prop[3])))
    loglikelihood1 <- ll2(prop)
    #loglikelihood1 <- sum(log(f(x,mu=prop[1],sig=exp(prop[2]), xi=prop[3])))
    logprior1 <- 0
    logtarget1 <- loglikelihood1 + logprior1
    if(is.nan(logtarget1)) logtarget1 <- -10000000
    if (runif(1) < min(1, exp(logtarget1 - logtarget0))) {
      theta0 <- prop
    }
    thetas <- rbind(thetas, theta0)
  }
  thetas[(burnin+1):(N+burnin),]
}

xx<-rnorm(928,mean = 3,sd=0.2)

mDLM<-DLMG(m0=0,C0=0.6,FF=1,G=1,V=1,W=1,T=length(xx),n=1,p=1,Xt=t(xx),muestra=1000)


yy<-MCMC01(xx,N=10000,theta0 = c(mean(mDLM),0.3,-0.1),pseudo_var = c(1,1,1),burnin = 500)


