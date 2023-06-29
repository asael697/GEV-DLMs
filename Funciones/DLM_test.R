
library(dlm)
library(bayesplot)
library(cowplot)
library(cmdstanr)
source("functions.R")

mod1 <- cmdstan_model("DLMs.stan")

n = 250
d = 1
k = 2

rm1 <- rdlm(n = n, G = diag(c(-0.7,0.9)), FF = c(1,1), m0 = c(0,0), 
            C0 = 100*diag(2), V = 0.6, W = 3*diag(2))

ts.plot(rm1$y)

data_list = list(n = n, d = d, k = k, y = rm1$y,
                 FF = matrix(1,ncol = d,nrow = k),
                 G = diag(c(-0.7,0.9)), 
                 m0 = c(0,0), C0 = 100*diag(2))

fit <- mod1$sample(data = data_list, chains = 4, parallel_chains = 4)
fit$summary(variables = c("tau","sigma","W","V"))

mu =  fit$draws(variables = "mu",format = "matrix")
y_pred = fit$draws(variables = "y_pred",format = "matrix")

g1 = ppc_ribbon(rm1$mt[,1],yrep = mu[,1:251])
g2 = ppc_ribbon(rm1$mt[,2],yrep = mu[,252:502])
plot_grid(g1,g2,ncol = 1)

ppc_ribbon(as.numeric(rm1$y),yrep = y_pred)

