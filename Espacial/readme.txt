1. El archivo pesos_espaciales.r contiene el calculo de los nucleos para la parte espacial
   aqui usamos los pesos tomando una regilla de 16 pesos (cada peso tiene su longitud y latitud 
   correspondiente) esto no contiene el desarrollo en metropolis-hasting.

2. El archivo metropolis_espacial contiene lo que es los nucleos calculados tomado ubicacion por ubicacion
   es decir la estacion 1 tiene longitud y latitud al igual que las demas estaciones, entonces se 
   genero una combinacion de 18x18 nucleos ya que son 18 estaciones.
   
   Aqui si hacemos el calculo de la estimacion de los parametros lo cual genera una estimacion de 
   152 muestreos por parametros es decir 3X152 ya que son 3 parametros.

   obs: recordemos que no se toma la estacion i con la misma estacion i puesto no tendria centido.

3. el archivo estaciones_ubicaciones contiene la longitud y latitud mas la altura de cada estacion 
   solo se tomo la longitud y la latitud. 