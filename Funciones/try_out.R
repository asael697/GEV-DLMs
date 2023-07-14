
library(cmdstanr)
library(cowplot)
library(posterior)
library(bayesplot)
library(ggplot2)

setwd("/Users/asael_am/Documents/GEV-DLMs/Funciones")
source("functions.R")

mod <- cmdstan_model("gev.stan")
mod1<- cmdstan_model("gev1.stan")

n = 559

preal <- c(0,rlnorm(1,meanlog = 0.1),rnorm(1))

# datos
y = SpatialExtremes::rgev(n,preal[1],scale = preal[2],shape = preal[3])

ggplot(data.frame(y),aes(x = y))+geom_density(fill = "darkred")+
  labs(x = "y",y = "conteo",title = "Histograma de los datos simulados y")


data_list1 = list(n = n, y = y)

fit <- mod1$sample(data = data_list1, chains = 4, parallel_chains = 4, max_treedepth = 12, 
                   refresh = 500, iter_sampling = 2000, iter_warmup = 2000,adapt_delta = 0.9)

fit$summary(variables = c("scale",'k'))
preal

mcmc_combo(
   fit$draws(variables = c("scale",'k'))
)

ppc_ribbon(y = as.numeric(y),yrep = mu_post)
ppc_ribbon(y = eta,yrep = eta_post)

post = fit$draws(variables = c("k","scale","mu"),format = "matrix")

y_rep = apply(post,1,FUN = function(x){
  x = as.numeric(x)
  SpatialExtremes::rgev(n,loc = 0,scale = x[2],shape = x[1]) + x[4:262]
})

ppc_ribbon(y = y, yrep = t(y_rep))
ppc_dens_overlay(y = y, yrep = t(y_rep))
