#library(MASS)
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
library(MARSS)
library(MASS)


datosds<-read_xlsx("C:/Users/Elvis/Documents/estaciones_ubicacion.xlsx")

data<-select(datosds , 2:3)


#aqui leemos las longitudes y latitudes de los datos

x<-c(data$latitud)
y<-c(data$longitud)


pas<-cbind(x,y)

# definimos las matrices para alamcenar el nucelo
Kij<-matrix(,nrow = 18,ncol=18)
v<-matrix(,nrow = 18,ncol=18)
w<-matrix(,nrow = 18,ncol=18)
h<-matrix(,nrow = 18,ncol=18)


#----calculamos el nucleo comparando cada estacion 
#--- una a una mediante formual
#---- exp([-d||S_i - S_j||^2]/2 )
#---- tomamos el valor de d=1
#---- S_i es la longitud y latitud de las estacion i
#---- S_j es la longitud y latitud de las estacion j
#----- usamos norma euclidea para verificar la distancia entre estaciones 


for (k in 1:18) {
  for (i in 1:18) {
    j=1
    while(j<2){
      v[k,i]<-(pas[k,j]-pas[i,j])*(pas[k,j]-pas[i,j])
      w[k,i]<-(pas[k,(j+1)]-pas[i,(j+1)])*(pas[k,(j+1)]-pas[i,(j+1)])
      h[k,i]<-(v[k,i]+w[k,i])*(v[k,i]+w[k,i])
      Kij[k,i] = exp(-0.5*(h[k,i]))
      j <- j+1
    }
    
  }
  
}



#----paso dlm con FFBS 
#---- Recordemos que aqui tomamos luego la media del vector obtenido
#---- esta es la parte que se debe mejorar para que genera un vector por cada iteracion


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
  
  
  H_t=array(dim=c(T,p,p))
  
  h_t=array(dim=c(T,p,n))
  
  
  M_t= array(dim=c(T,p,n))
  theta_t= array(dim=c(T,p,p))
  Theta_t_mayus= array(dim=c(T,p,p))
  
  at[,1] = G%*%m0
  Rt[1,,] = (G%*%C0%*%t(G))+W
  ft[,1] = FF%*%at[,1]
  Qt[,1] = FF%*%Rt[1,,]%*%t(FF) + V
  At[1,,] = Rt[1,,]%*%t(FF)%*%solve(Qt[,1])
  et[,1] = Xt[,1] - ft[,1]  
  mt[,1] = at[,1] + At[1,,]%*%t(et[,1])
  Ct[1,,] = Rt[1,,] - At[1,,]%*%t(Qt[1,1])%*%t(At[1,,])
  
  
  
  M_t[1,,] = mvrnorm(1,ft[,1],Qt[,1])
  
  
  H_t[1,,] = solve( (solve(Ct[1,,]))+t(G)%*%solve(W)%*%G )
  
  
  h_t[1,,] = H_t[1,,]%*%(solve(Ct[1,,])%*%mt[,1] + t(G)%*%solve(W)%*%M_t[1,,])
  
  #forward Filter
  for(t in 2:T){
    #t=10
    at[,t] = G%*%mt[,t-1]    
    Rt[t,,] = (G%*%Ct[t-1,,]%*%t(G)) + W
    #Rt
    ft[,t] = FF%*%at[,t]
    Qt[,t] = FF%*%Rt[t,,]%*%t(FF) + V
    At[t,,] = Rt[t,,]%*%t(FF)%*%solve(Qt[,t])
    et[,t] = Xt[,t] - ft[,t]  
    mt[,t] = at[,t] + At[t,,]%*%t(et[,t])
    Ct[t,,] = Rt[t,,] - At[t,,]%*%t(Qt[,t])%*%t(At[t,,])
    
    
    #H_t[,t] = solve( (solve(Ct[t,,]))+t(G)%*%solve(W)%*%G )
    
    #theta_t[,t] = rnorm(1 , mean =at[,t], sd = sqrt(Rt[t,,]))
    #M_t[,t] = rnorm(1, mean = ft[,t] , sd = sqrt(Qt[,t]))
    #Theta_t_mayus[,t] = rnorm(1 , mean = M_t[,t] , sd = sqrt(Ct[t,,]))
    
    #print(t)
    
    theta_t[t,,] = mvrnorm(1,at[,t],Rt[t,,])
    M_t[t,,] = mvrnorm(1,ft[,t],Qt[,t])
    Theta_t_mayus[t,,] = mvrnorm(1,M_t[t,,],Ct[t,,])
    
    
    H_t[t,,] = solve( (solve(Ct[t,,]))+t(G)%*%solve(W)%*%G )
    
    #h_t[,t] = H_t[t,,]%*%(solve(Ct[t,,])%*%mt[,t] + t(G)%*%solve(W)%*%M_t[t,,])
    
    h_t[t,,] = H_t[t,,]%*%(solve(Ct[t,,])%*%mt[,t] + t(G)%*%solve(W)%*%M_t[t,,])
    #print(t)
  }
  
  #BACKWAR SAMPLIN
  

  
  if(muestra==1){
    
    #muestra=1
    mus=array(dim=c(T,1,p))
    #mus = rep(0,T)
    #mus[T] = rnorm(1,at[,T],sqrt(Ct[T,,]))
    mus[T,,] = mvrnorm(1,at[,T],Ct[T,,])
    for (t in (T-1):1) {
      #mus[t]=rnorm(1,h_t[,(t+1)],sqrt(H_t[(t+1),,]))
      #t= T-1 
      mus[t,,] = mvrnorm(1,h_t[(t+1),,],H_t[(t+1),,])
    }
  }else{
    #mus = matrix(0,muestra,T)
    mus=array(dim=c(T,muestra,p))
    #mus[T,,]=rnorm(muestra,at[,T],sqrt(Ct[T,,]))
    mus[T,,] = mvrnorm(muestra,at[,T],Ct[T,,])
    for (t in (T-1):1) {
      #mus[t,,]=rnorm(muestra,h_t[,(t+1)],sqrt(H_t[(t+1),,]))
      
      mus[t,,]=mvrnorm(muestra,h_t[(t+1),,],H_t[(t+1),,])
    }
  }
return(mus)
}  






#----- lectura de los datos que superaron el umbral por estacio

fgfg <- read_xlsx("C:/Users/Elvis/Documents/SU_100/TODAS_1000.xlsx")
f1 <- nrow(fgfg) 
c1 <- ncol(fgfg)


#medias_2 = matrix(0 , nrow=18 , ncol = 1)
I<- diag(18)
#0.6*I



#----- matrices de almacenamiento de sampler von metropolisy paquete evdbayes

mu_col <- matrix(0,10001,152)
sigma_col <- matrix(0,10001,152)
xi_col <- matrix(0,10001,152)

acum<-0


m0<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
FFF<-c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)


#--- codigo de metropolis espacial
#-----  F toma el valor de la casilla de la matriz i , j
#------ recorddemos que es estcion i con estacion j (no tendria sentido hacer j , i porque seria lo mismo por la norma)
#------ aqui lo hacemos para cada estacion
#------ F seria un vector donde la primer cordenada sera el elemnto de la columna Kij
#------ El primer ciclo for es para ubcar los datos de la estacion 1 en un solo vector y asi 
#------ solo tener los datos correspondientes a la estacion
#------ llamanos a la funcion dlm para obtener la media del parametro \mu
#------ condicionamos a que cada la estacion 1 sea comparada en espacio mediante Kij
#------ con el resto de las estaciones.
#------ almacenamos 152 muestreos por cada parametro ya que ese es el resultado de evaluar estacion 
#------ por estacion


for (i in 1:18) {
  temp_1 = na.exclude(fgfg[,i])
  
  mat <- diag(c(1000, 50, 100))  #matriz diagonal
  
  pn <- prior.norm(mean = c(0,0,0), cov = mat) 
  
  Go<-temp_1
  tam<-nrow(Go)
  
  
  for (hh in (1:18)) {
    
  if(i<hh){
  
    mu_tep<-mean(DLMG(m0=m0,C0=0.6*I,FF=t(Kij[i,hh]*FFF),G=I,V=1,W=I,T=tam,n=1,p=18,Xt=t(Go),muestra=1000))
  
    
    
    n <- 10000 ; t0 <- c(mu_tep , 1.3 , 0.4) ; #s <- c(.02,.01,.1)
    s<- c(.1,.05,.04)
  
   ptpmc <- posterior(n, t0, prior = pn, lh = "gev", data = temp_1, psd = s)
  
   acum<-acum + 1
  
   mu_col[,acum] <- ptpmc[,1]
   sigma_col[,acum] <- ptpmc[,2]
   xi_col[,acum] <- ptpmc[,3]
   
  }
 }
}


