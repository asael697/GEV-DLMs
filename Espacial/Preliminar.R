
source("funciones/functions.R")

library(ggplot2)
library(bayesplot)
library(posterior)

load("Datos/K_matrix.Rdata")
Estaciones <- readxl::read_excel("Datos/Estaciones.xlsx")

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