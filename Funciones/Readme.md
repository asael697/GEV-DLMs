    source("functions.R")

    library(posterior)
    library(bayesplot)
    library(ggplot2)
    library(loo)

## Pruebas para validar los MCMC en distribuciones univariadas

Este cuaderno hace un prueba de la funcionalidad del algortimo de
Metropolis implementado en el archivo `functions.R`, en este experimento
se simularan datos de una normal univariada mediante las siguientes
ecuaciones

$$
y = \mu + \sigma \varepsilon, \quad \varepsilon \sim N(0,1),\\
\mu \sim N(0,1), \quad \sigma \sim logN(1,1).
$$

Los datos se simulan mediante el siguiente codigo:

    #priors
    mu = rnorm(1)
    sigma = rlnorm(1,meanlog = 1)
    # datos
    y = rnorm(n = 100,mean = mu,sd = sigma)

    ggplot(data.frame(y),aes(x = y))+geom_histogram(bins = 20)+
      labs(x = "y",y = "conteo",title = "Histograma de los datos simulados y")

![](Figures/unnamed-chunk-2-1.png)

Ahora bien, previo a la estimacion de los parametros, actualizamos las
funciones `loglik()`, `log_prior()` y `inits()` para que calculen la
verosimilitud, el logaritmo de la prior y valores iniciales, de forma
correcta. Estas funciones son necesarias para el salto de metropolis
(`metropolis_step`) y para la iniciacion correcta del algoritmo
(`metropolis_sampler()`).

    # Estimacion de la log-verosimilitud
    loglik <- function(y,theta){
        d = sum(dnorm(x = y,mean = theta[1],sd = theta[2],log = TRUE))
        d = ifelse(is.na(d),-10^(64),d)
       return(sum(d))
    }

    # Calculo de la log prior logP(theta)
    log_prior <- function(theta){
      d1 = dnorm(theta[1],log = TRUE)
      d2 = dlnorm(theta[2],meanlog = 1,log = TRUE)
      return(d1 + d2)
    }

    # Generacion de los valores iniciales
    inits <- function(y,d = 2){
      c(rnorm(1),rlnorm(1,meanlog = 1))
    }

Ahora bien, procedemos a generar nuestras cadenas de Markov usando el
algoritmo de Metropolis. En este caso la funcion de salto es una normal
multivariada simetrica $proposed \sim N_2(prop, dia(1,1))$,
con matriz de covarianza siendo la matriz identidad. Se generaron 4
cadenas de Markov de un total de 10,000 iteraciones por cadena
(`iter = 10,000`), donde se eliminaron las primeras 2,000 iteraciones
(`burn-in = 2500`), y se acortaron las cadenas con un `thinning` de cada
5 saltos (`thin = 5`) para evitar que las cadenas se quedaran
estancadas.

    scale_mat  = diag(c(1,1))

    start = Sys.time()
    post1 = sampling(y,iter = 5000, scale = scale_mat, thin = 3)
    print( Sys.time() - start)
    ## Time difference of 6.675022 secs

Los resultados obtenidos muestran convergencia en las cadenas $\hat R approx 1$,
pero los tamaños de muestra efectivos son lo muy bajos, como para
aceptar las simulaciones obtenidas. Notemos que las posteriors si
recuperaron los valores reales de $\mu$ (-1.2886447) y $\sigma$ (1.0481244),
pero dada su poca convergencia no quantifican bien la incertidumbre de
los parametros.

    colnames(post1) = c("mu","sigma","accpetence_rate",".chain")
    post_df = as_draws_df(post1)

    summarise_draws(post_df)
    ## # A tibble: 3 × 10
    ##   variable         mean median     sd    mad     q5    q95  rhat ess_b…¹ ess_t…²
    ##   <chr>           <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <dbl>   <dbl>   <dbl>
    ## 1 mu            -1.14    -1.14 0.107  0.0995 -1.33  -0.969  1.02    543.    654.
    ## 2 sigma          1.05     1.05 0.0777 0.0740  0.924  1.17   1.01    637.    694.
    ## 3 accpetence_r…  0.0144   0    0.119  0       0      0      1.00  19329.  19329.
    ## # … with abbreviated variable names ¹​ess_bulk, ²​ess_tail

La siguiente visualización muestra las posterios y cadenas para ambos
parametros, los cadenas quedan en puntos de acumulacion debido al bajo
ratio de rechazo y las posteriors son deformes, indicando que es
necesario realizar mas iteraciones.

    color_scheme_set("mix-gray-brightblue")
    mcmc_combo(post_df,pars = c("mu","sigma"))

![](Figures/unnamed-chunk-6-1.png)

## Adjusted scale metropolis

Una alternativa a aumentar el numero de iteraciones es controlar el
factor de escala de la cadena, para eso haremos un algoritmo de
Metropolis-Hastings con factor de escala ajustado, en este caso la
distribution de salto es:

$$x\_k = x\_{k-1} + \sqrt h \Sigma^{1./2} \varepsilon, \quad \varepsilon \sim N\_d(0,I).$$
Donde $H$ es el factor de escala que corrige la matriz de covarianza
$\Sigma$, y d es la dimension de los parametros a obtener ($x$). Un mejor
control de $h$ aumentaria la tasa de aceptacion del MCMC. para este caso
elegimos $h = 0.01$ y reducimos el acortamiento a 1

    scale_mat  = diag(c(1,1))

    start = Sys.time()
    post1 = sampling(y,iter = 5000, scale = scale_mat, thin = 1, h = 0.1)
    print( Sys.time() - start)
    ## Time difference of 2.879897 secs

El tiempo de ejecucion se ve reducido, ademas que la convergencia y
tamaño de las cadenas es mejor con respecto al tiempo anterior, y la
tasa de acceptacion incremento.

    colnames(post1) = c("mu","sigma","accpetence_rate",".chain")
    post_df = as_draws_df(post1)

    summarise_draws(post_df)
    ## # A tibble: 3 × 10
    ##   variable         mean median     sd    mad     q5    q95  rhat ess_b…¹ ess_t…²
    ##   <chr>           <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <dbl>   <dbl>   <dbl>
    ## 1 mu             -1.15   -1.15 0.105  0.110  -1.32  -0.971  1.00   1397.   1471.
    ## 2 sigma           1.04    1.04 0.0735 0.0743  0.927  1.16   1.00   1472.   1756.
    ## 3 accpetence_ra…  0.124   0    0.329  0       0      1      1.00  14966.     NA 
    ## # … with abbreviated variable names ¹​ess_bulk, ²​ess_tail

La siguiente visualización muestra las posterios y cadenas para ambos
parametros, las cadenas se visualizan estacionaria y las posteriors
menos deformes, indicando convergencia.

    color_scheme_set("viridisC")
    mcmc_combo(post_df,pars = c("mu","sigma"))

![](Figures/unnamed-chunk-9-1.png)

Para una mayor tasa de rechazo se puede incrementar el tamanio de
muesta, mejorar el factor de escala, o intentar un MCMC adaptativo.

## Metropolis adjusted Lavengin algorithm

El siguiente algoritmo, hace el salto de la cadena mediante un descenso
en gradiente, de tal forma que el la funcion de salto sea normal
simetrica, el salto se describe mediante la siguiente ecuación

$$x\_k = x\_{k-1} + h \Sigma \nabla \log f(x\_{k-1})/2 +  \sqrt{h} \Sigma^{1/2}\varepsilon, \quad \varepsilon \sim N\_d(0,I)$$

En este caso $\Sigma$ es la matriz de covarianza de la funcion de salto, h es
el factor de escala, $\log f$ es el logaritmo de la funcion objetivo tal
que $f(x) \propto likelihood(x)\pi(x)$, y ∇ es el
gradiente de una funcion. Se sabe que cuando estos algoritmos usan un
factor de escala optimo, el salto es ergodico y la tasa de aceptación es
de alrededor de $r \approx 0.56$. Empiricamente, se sabe que el factor de
escala optimo esta en la region $0 < h << 1$.

El siguiente codigo genera 2000 iteraciones de un MALA con un warm-up de
1000 iteraciones, sin acortamiento, y una matriz de covarianza con
diagonal 1, y 0.1 de covarianza. El factor de escala en este caso es
bastante cercano a 0, $h = 0.01$.

    scale_mat = matrix(c(1,0.1,0.1,1),nrow = 2)

    start = Sys.time()
    post2 = sampling(y,iter = 2000, scale = scale_mat, thin = 1, h = 0.01,MALA = TRUE)
    print( Sys.time() - start)
    ## Time difference of 9.080563 secs

Aunque los tiempos de ejecución incrementaron debido al computo de los
gradientes, el metodo solo duplico el tiempo de ejecución alcanzando la
convergencia con menos iteraciones que los primeros dos. Donde los
factores de convergencia y tamños de muestra efectivos se indican
convergencia y la tasa de aceptación es mucho mayor que la de los
metodos anteriores y cercana a la optima.

    colnames(post2) = c("mu","sigma","accpetence_rate",".chain")
    post_df = as_draws_df(post2)

    summarise_draws(post_df)
    ## # A tibble: 3 × 10
    ##   variable         mean median     sd    mad     q5    q95  rhat ess_b…¹ ess_t…²
    ##   <chr>           <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <dbl>   <dbl>   <dbl>
    ## 1 mu             -1.15   -1.15 0.100  0.101  -1.31  -0.978  1.00    898.   1325.
    ## 2 sigma           1.04    1.04 0.0750 0.0736  0.925  1.17   1.00   1275.   1722.
    ## 3 accpetence_ra…  0.528   1    0.499  0       0      1      1.00   7447.     NA 
    ## # … with abbreviated variable names ¹​ess_bulk, ²​ess_tail

Los graficos indican que las cadenas son ergodicas y las posteriors son
unimodales y simetricas indicando convergencia y una buena
cuantificacion de la incertidumbre de los parametros.

    color_scheme_set("mix-gray-green")
    mcmc_combo(post_df,pars = c("mu","sigma"))

![](Figures/unnamed-chunk-12-1.png)

## Posterior predictive checks

Revisar el ajuste de los datos es una buena praxis para validar los
supuestos del modelo, para eso es necesario realizar posterior
predictive checks. En este caso, extraemos 500 simulaciones aleatorias,
y las comparamos con los datos reales. El grafico siguiente muestra que
el modelo provee un buen ajuste de los datos.

    color_scheme_set("viridisC")

    y_rep = apply(post1,1,FUN = function(x){
      x = as.numeric(x)
      rnorm(n = 100,mean = x[1],sd = x[2])
    })

    s = sample(1:2000,size = 500)
    ppc_dens_overlay(y = y,yrep = t(y_rep)[s,])

![](Figures/unnamed-chunk-13-1.png)

Finalmente, realizamos validacion cruzada para determinar el ajuste del
modelo, para eso es necesario computar la matriz de de
log-verosimilitud, al aproximar la elpd (Expected log-predictive
density) notamos que no hay advertencias de un mal ajusto de los valores
de pareto, por lo tanto aceptamos el modelo propuesto.

    log_lik <- function(y,theta) {
      t(apply(theta,1,function(theta)
       dnorm(x = y,mean = theta[1],sd = theta[2],log = TRUE)
      ))
    }
    LL = log_lik(y,post1)

    loo1 = loo(LL, relative_eff(exp(LL)))
    loo1
    ## 
    ## Computed from 20000 by 100 log-likelihood matrix
    ## 
    ##          Estimate   SE
    ## elpd_loo   -146.1  5.9
    ## p_loo         1.7  0.2
    ## looic       292.2 11.8
    ## ------
    ## Monte Carlo SE of elpd_loo is 0.0.
    ## 
    ## All Pareto k estimates are good (k < 0.5).
    ## See help('pareto-k-diagnostic') for details.

## Pruebas para validar el filtro de Kalman hacia adelante

Para validar el algoritmo se simula una serie de tiempo multivariada de
dimension $m = 2$ y tamaño $n = 250$ observaciones. Los parámetros
latentes $\mu \in \mathbb R^3$, Usando un DLM constante. En
este caso, las matrices $G, FF, V$ y $W$ se simulan
aleatoriamente, y los parametros latentes iniciales son:

$$m_0 = [0, 0, 0] \quad y \quead C_0 = 100*diag([1, 1, 1]).$$

Dicha serie de tiempo se simula mediante el siguiente codigo:

    library(dlm)
    ## 
    ## Attaching package: 'dlm'
    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     %+%

    rm = dlmRandom(m = 2,p = 3,nobs = 250)

Los parametros simulados del DLM constante son

    list(FF = rm$mod$FF, G = rm$mod$GG, V = rm$mod$V, W = rm$mod$W)
    ## $FF
    ##             [,1]      [,2]       [,3]
    ## [1,] -0.33063053 1.1097491 -0.1600778
    ## [2,] -0.09691218 0.9632287 -0.2066727
    ## 
    ## $G
    ##              [,1]        [,2]         [,3]
    ## [1,]  0.024024948  0.04277632 -0.009977307
    ## [2,]  0.106508870  0.04660221 -0.001959571
    ## [3,] -0.001793041 -0.11177610 -0.045605856
    ## 
    ## $V
    ##           [,1]      [,2]
    ## [1,]  7.356003 -4.008054
    ## [2,] -4.008054  2.381135
    ## 
    ## $W
    ##           [,1]     [,2]      [,3]
    ## [1,] 2.6824080 1.123869 0.1234766
    ## [2,] 1.1238690 2.430207 1.9838683
    ## [3,] 0.1234766 1.983868 2.8324780

Ajustamos el filtro de Kalman hacia delante con nuestra implementación y
la realizada por el paquete `DLM`. Comparamos los resultados de
$m_t$ y $C_t$ obtenidos, con la norma 2 para el
vector de medias y la norma 1 para las matrices de covarianza.

    # Nuestra implementacion 
    kf1 = forward_filter(y = rm$y, G = rm$mod$GG, FF = t(rm$mod$FF),
                   V = rm$mod$V, W = rm$mod$W, m0 = rm$mod$m0,
                   C0 = rm$mod$C0)

    # Implementacion paquete DLM
    filt1 = dlmFilter(y = rm$y,mod = rm$mod)

$$e_{mt} = || \text{OUR} - \text{dlm}||_2.$$

    e_mt = kf1$mt - filt1$m
    apply(e_mt,2,function(x) sqrt(sum(x^2)) )
    ## [1] 1.554387e-13 4.337018e-14 3.150671e-14

y

$$e_{Ct} = || \text{OUR} - \text{dlm}||_2.$$


    e_Ct = rep(0,250)

    for (i in 1:250) {
      C_temp = filt1$U.C[[i]] %*% diag(filt1$D.C[i,]^2) %*% t(filt1$U.C[[i]])
      e_Ct[i] = norm(kf1$Ct[[i]] - C_temp,type = "1")
    }

    summary(e_Ct)
    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## 0.000e+00 2.498e-15 2.498e-15 2.545e-15 2.498e-15 7.161e-15

Observamos que todos los errores son aproximadamente cero y por lo tanto
validamos nuestra implementacion del Filtro de Kalman hacia adelante.
