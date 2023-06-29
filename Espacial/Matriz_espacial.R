
#####################################################################
### Script para Calcular la matriz de distancias entre estaciones
#####################################################################

datosds = read_xlsx("estaciones_ubicacion.xlsx")

pas = cbind(datosds$latitud,datosds$longitud)

Kij<-matrix(,nrow = 18,ncol=18)
v<-matrix(0,nrow = 18,ncol=18)
w<-matrix(0,nrow = 18,ncol=18)
h<-matrix(0,nrow = 18,ncol=18)

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

rm(i,j,k,h,pas,v,w,datosds)

save.image("K_matrix.Rdata")
