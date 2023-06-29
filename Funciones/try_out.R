
library(cmdstanr)
library(cowplot)
library(posterior)
library(bayesplot)
library(ggplot2)

setwd("/Users/asael_am/Documents/GEV-DLMs/Funciones")
source("functions.R")

mod <- cmdstan_model("gev_dlm.stan")
mod1 <- cmdstan_model("GEV-DLM.stan")

n = 259

preal <- c(rlnorm(1),rnorm(1),rlnorm(1,meanlog = -2))
y0 <- rdlm(n, G = 1, FF = 1, V = preal[3], W = 0.6*preal[3])

plot_ts(y0$mt[-1])
preal

eta = SpatialExtremes::rgev(n)
y = preal[1]*(exp(preal[2]*eta) - 1)/preal[2] + rnorm(n = n,mean = y0$mt[-1],sd = preal[3])
plot_ts(y, name = "y",color = "darkblue")

#data_list = list(n = n, y = y, d = 1, G = matrix(1,nrow = 1,ncol = 1), 
#                 m0 = 0, C0 = matrix(100,nrow = 1,ncol = 1))

data_list = list(n = n, y = y, d = 1, phi = 1, m0 = 0, C0 = 100)

fit <- mod$sample(data = data_list, chains = 4, parallel_chains = 4, max_treedepth = 11, 
                  refresh = 500, iter_sampling = 2000, iter_warmup = 4000,adapt_delta = 0.9)

fit$summary(variables = c('k',"scale","tau","sigma"))
preal

mcmc_combo(
   fit$draws(variables = c('k',"scale","tau","sigma"))
)

mu_post = fit$draws(variables = "mu",format = "matrix")
eta_post = fit$draws(variables = "eta",format = "matrix")

ppc_ribbon(y = as.numeric(y0$mt),yrep = mu_post)
ppc_ribbon(y = eta,yrep = eta_post)

post = fit$draws(variables = c("k","scale","mu"),format = "matrix")

y_rep = apply(post,1,FUN = function(x){
  x = as.numeric(x)
  SpatialExtremes::rgev(n,loc = 0,scale = x[2],shape = x[1]) + x[4:262]
})

ppc_ribbon(y = y, yrep = t(y_rep))
ppc_dens_overlay(y = y, yrep = t(y_rep))
