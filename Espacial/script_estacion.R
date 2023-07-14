source("Espacial/Preliminar.R")

est = 6
S_1 =  na.exclude(Estaciones[,est])
V = c(0.1,0.2,0.1,0.1,0.1,0.1)

plot.ts(S_1)
hist(S_1,breaks = 50)

n = length(S_1)
iter = 3000

## Construccion del DLM 
buildFun <- function(x){ 
  dlm::dlm(GG = 0.5*I, FF = t(Kij[est,]), V = V[est], W = exp(x[1])*I, m0 = rep(0,18), 
           C0 = 100*I)
 }

dlm_MLE <- dlmMLE(S_1, parm = c(0,0), build = buildFun)
dlm_MLE$par

dlm_temp <- dlm::dlmFilter(y = S_1,mod =  buildFun(dlm_MLE$par))

## Aplicacion de FFBS
ft = matrix(0,nrow = 500, ncol = n+1)

for(i in 1:500){
  Mu = dlm::dlmBSample(dlm_temp)
  ft[i, ] = Mu %*% Kij[est, ]
}

# Media a posterior del DLM 
ft_mean = apply(ft, 2, mean)
plot.ts(ft_mean)

ft = ft[,1:n +1]
colnames(ft) = paste0("mu",1:n)

# Estimacion de los parametros del GEV
z = as.numeric(S_1) - ft_mean[-1]
plot.ts(z)
hist(z,breaks = 50)

data_list1 = list(n = n, y = S_1)

fit <- mod1$sample(data = data_list1, chains = 4, parallel_chains = 4, 
                     max_treedepth = 12, refresh = 500, iter_sampling = iter,
                     iter_warmup = 2000,adapt_delta = 0.9)

# Revisamos las posterior
fit$summary()
mcmc_combo(fit$draws(variables = c("scale",'k')))

#Agregamos las cadenas del FFBS al Metropolis
post2 = fit$draws(variables = c("scale",'k'),format = "matrix")
post2 = cbind(post2[sample(size = 500,x = 1:(4*iter)),],ft)

# Posterior predictive checks
y_rep = apply(post2,1,FUN = function(x){
  x = as.numeric(x)
  SpatialExtremes::rgev(n,loc = 0,scale = x[1],shape = x[2]) + x[1:n + 2]
})

g1 = ppc_dens_overlay(y = as.numeric(S_1),yrep = t(y_rep),trim = FALSE)
g2 = ppc_ribbon(y = as.numeric(S_1),yrep = t(y_rep)) + 
  labs(title = "Revisión de la predictiva",subtitle = paste("Estación",est),
       x = "tiempo",y = "Precipitaciones")

cowplot::plot_grid(g1,g2,ncol = 1)
