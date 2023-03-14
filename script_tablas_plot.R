

library(posterior)
library(bayesplot)
library(xtable)

mu = read.csv("~/Downloads/datos dos metodos/datos hamilton/mu_hamilto.csv")
sigma = read.csv("~/Downloads/datos dos metodos/datos hamilton/sigma_hamilto.csv")
x = read.csv("~/Downloads/datos dos metodos/datos hamilton/xi_hamilto.csv")

dat = draws_df(mu = mu[,2],sigma = sigma[,2],x = x[,2])
##Burn-in
dat1 = subset_draws(dat, iteration =1000:18000)

# Tabla resumen
xtable(summarise_draws(dat1))
# Grafico
mcmc_combo(dat1)

#grafico ACF

library(forecast)
ggAcf(dat1[,1])
