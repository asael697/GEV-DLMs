

library(cmdstanr)
library(cowplot)
library(posterior)
library(bayesplot)
library(ggplot2)
library(dlm)


load("Datos/K_matrix.Rdata")
Estaciones <- readxl::read_excel("Datos/Estaciones.xlsx")
Estaciones <- as.matrix(Estaciones)
mod <- cmdstan_model("Espacial/gev.stan")
mod1<- cmdstan_model("Espacial/gev1.stan")

loglik <- function(y,theta){
  d = SpatialExtremes::dgev(y, loc = 0, scale = theta[1],
                            shape = theta[2], log = TRUE)
  return(sum(d))
}

# Calculo de la log prior logP(theta)
log_prior <- function(theta){
  d2 = dt(theta[1],df = 3,log = TRUE) + log(2)
  d3 = dnorm(theta[2],sd = 1,log = TRUE)
  
  return(d2+d3)
}

# Generacion de los valores iniciales
inits <- function(){
  c(rlnorm(1),rnorm(1))
}

I <- diag(18)