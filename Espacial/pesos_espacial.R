library(posterior)
library(bayesplot)
library(xtable)
library(MASS)
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

#cargamos las cordenadas de las estaciones meteorologicas

cordenadas_reales = read.xlsx("C:/Users/elvis/Downloads/descargas todo/datos parte espacial/estaciones_ubicacion.xlsx")


#------en este caso hicimos una aproximacion mediante una cuadricula en el mapa
#------ del estado de whashinton para obterner las cordenadas de los pesos
#------ como suuiere el articulo de huerta, podemos usar paquetes para poder obtener 
#------ los pesos tambien.
#---- peso-long es la longitud
#---- peso_lat es la latitud, recodemos que una ubicacion corresnde de longitud y latitud
#-----

peso_lon<-c(-123.633351, -121.733351, -119.733351, -118.033351, 
            -123.633351, -121.733351, -119.733351, -118.033351,
            -123.633351, -121.733351, -119.733351, -118.033351,
            -123.633351, -121.733351, -119.733351, -118.033351)

peso_lat<-c( 47.936488 ,  47.936488 ,  47.936488 ,  47.936488 ,
             47.436488 ,  47.436488 ,  47.436488 ,  47.436488 ,
             46.936488 ,  46.936488 ,  46.936488 ,  46.936488 ,
             46.436488 ,  46.436488 ,  46.436488 ,  46.436488)


#--- aqui solo calculamos la primer parte de la norma euclidea que es el producto
#--- de los ((S_i(lat) - peso_lat)^2 + S_i(long) - peso_long)^2)
#--- ojo S_i(lat) es la latitud de la estacion i, recordemos que este esta almacenado
#--- en la variable cordenadas reales
#--- aqui se entiende que es la estacion con el peso mas cercano (ver mapa del estado)
#---

e1<- (((cordenadas_reales[1,2] - peso_lat[1])*(cordenadas_reales[1,2] - peso_lat[1]))+
        ((cordenadas_reales[1,3] - peso_lon[1])*(cordenadas_reales[1,3] - peso_lon[1])))

e2<- (((cordenadas_reales[2,2] - peso_lat[7])*(cordenadas_reales[2,2] - peso_lat[7]))+
        ((cordenadas_reales[2,3] - peso_lon[7])*(cordenadas_reales[2,3] - peso_lon[7])))

e3<- (((cordenadas_reales[3,2] - peso_lat[5])*(cordenadas_reales[3,2] - peso_lat[5]))+
        ((cordenadas_reales[3,3] - peso_lon[5])*(cordenadas_reales[3,3] - peso_lon[5])))

e4<- (((cordenadas_reales[4,2] - peso_lat[9])*(cordenadas_reales[4,2] - peso_lat[9]))+
        ((cordenadas_reales[4,3] - peso_lon[9])*(cordenadas_reales[4,3] - peso_lon[9])))

e5<- (((cordenadas_reales[5,2] - peso_lat[2])*(cordenadas_reales[5,2] - peso_lat[2]))+
        ((cordenadas_reales[5,3] - peso_lon[2])*(cordenadas_reales[5,3] - peso_lon[2])))

e6<- (((cordenadas_reales[6,2] - peso_lat[15])*(cordenadas_reales[6,2] - peso_lat[15]))+
        ((cordenadas_reales[6,3] - peso_lon[15])*(cordenadas_reales[6,3] - peso_lon[15])))

e7<- (((cordenadas_reales[5,2] - peso_lat[1])*(cordenadas_reales[5,2] - peso_lat[1]))+
        ((cordenadas_reales[5,3] - peso_lon[1])*(cordenadas_reales[5,3] - peso_lon[1])))

e8<- (((cordenadas_reales[8,2] - peso_lat[6])*(cordenadas_reales[8,2] - peso_lat[6]))+
        ((cordenadas_reales[8,3] - peso_lon[6])*(cordenadas_reales[8,3] - peso_lon[6])))

e9<- (((cordenadas_reales[9,2] - peso_lat[6])*(cordenadas_reales[9,2] - peso_lat[16]))+
        ((cordenadas_reales[9,3] - peso_lon[6])*(cordenadas_reales[9,3] - peso_lon[6])))

e10<- (((cordenadas_reales[10,2] - peso_lat[9])*(cordenadas_reales[10,2] - peso_lat[9]))+
         ((cordenadas_reales[10,3] - peso_lon[9])*(cordenadas_reales[10,3] - peso_lon[9])))

e11<- (((cordenadas_reales[11,2] - peso_lat[11])*(cordenadas_reales[11,2] - peso_lat[11]))+
         ((cordenadas_reales[11,3] - peso_lon[11])*(cordenadas_reales[11,3] - peso_lon[11])))

e12<- (((cordenadas_reales[12,2] - peso_lat[1])*(cordenadas_reales[12,2] - peso_lat[1]))+
         ((cordenadas_reales[12,3] - peso_lon[1])*(cordenadas_reales[12,3] - peso_lon[1])))

e13<- (((cordenadas_reales[13,2] - peso_lat[14])*(cordenadas_reales[13,2] - peso_lat[14]))+
         ((cordenadas_reales[13,3] - peso_lon[14])*(cordenadas_reales[13,3] - peso_lon[14])))


e14<- (((cordenadas_reales[14,2] - peso_lat[15])*(cordenadas_reales[14,2] - peso_lat[15]))+
         ((cordenadas_reales[14,3] - peso_lon[15])*(cordenadas_reales[14,3] - peso_lon[15])))

e15<- (((cordenadas_reales[15,2] - peso_lat[16])*(cordenadas_reales[15,2] - peso_lat[16]))+
         ((cordenadas_reales[15,3] - peso_lon[16])*(cordenadas_reales[15,3] - peso_lon[16])))

e16<- (((cordenadas_reales[16,2] - peso_lat[8])*(cordenadas_reales[16,2] - peso_lat[8]))+
         ((cordenadas_reales[16,3] - peso_lon[8])*(cordenadas_reales[16,3] - peso_lon[8])))

e17<- (((cordenadas_reales[17,2] - peso_lat[7])*(cordenadas_reales[17,2] - peso_lat[7]))+
         ((cordenadas_reales[17,3] - peso_lon[7])*(cordenadas_reales[17,3] - peso_lon[7])))

e18<- (((cordenadas_reales[18,2] - peso_lat[7])*(cordenadas_reales[18,2] - peso_lat[7]))+
         ((cordenadas_reales[18,3] - peso_lon[7])*(cordenadas_reales[18,3] - peso_lon[7])))


#--- almacenamos los las diferencias en un vector que llamamos nucleso
#---
nucleos<-c(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18)



para_espa<-matrix(NA,nrow = 1,ncol = 18)

#----- almacenamiento de los K^t en un vector, para tomarlos 1 a 1.
for (i in 1:18) {
  para_espa[,i]<-exp(-1*((dis/2)*sqrt(nucleos[i])))
}

