###################################################
###################INTERPOLACIÓN###################
##################POR GABRIEL ROA##################
###################################################
###################CI: 25919459####################
#########Asignatura: Programación Numérica#########
############Profesor: Abelardo Monsalve############
####################UCLA - DCYT####################
###################################################

InterpolacionCanonica<-function(X,fx,a){ ### Interpolación de forma canónica.
	if (length(X)!=length(fx)){ ### Como ingresa dos vectores distintos para x y fx, verifica que estos tengan igual número de elementos.
		return("ERROR: Debe ingresar un número equivalente de puntos (x,fx)")
	}

	k = length(X) 

	A = matrix(1,nrow=k,ncol=k) ### Crea la matriz de Vandermonde vacía.

	print("Para la interpolación canónica, haremos uso de la llamada Matriz de Vandermonde.")
	print("Al resolver el sistema de ecuaciones asociado a esta matriz, obtendremos los coeficientes del polinomio interpolador.")
	print("La matriz de Vandermonde en este caso de interpolación es la siguiente:")

	for (j in 1 : k){
		for (i in 2 : k){
			A[j,i]=X[j]^(i-1) ### Llena los coeficientes de la matriz de Vandermonde.
		}
	}

	print(A)

	Ab=cbind(A,fx) ### Crea la matriz a la cual se aplicará eliminación gaussiana.

	print("Mientras que el sistema asociado que será resuelto mediante eliminación gaussiana es el siguiente:")
	print(Ab)

	for (i in 1:(k-1))
	{ ### Eliminación Gaussiana.
		for (j in (i+1):k){
		factor=Ab[j,i]/Ab[i,i]
			if (factor!=0){
				Ab[j,]=Ab[j,]-factor*Ab[i,]
			}
		}
	}

	x=numeric(k)
	x[k]=Ab[k,(k+1)]/Ab[k,k]

	for (i in (k-1) : 1){ ### Sustitución. Tras esto, en x se almacenan los coeficientes del polinomio.
		suma = sum((Ab[i,1:k])*x[1:k])
		x[i]=(Ab[i,(k+1)]-suma)/Ab[i,i]
	}

	print("Resolviendo el sistema, obtenemos los siguientes coeficientes del polinomio interpolador:")
	print(x)

	interp = x[1] ### Inicializa con el primer coeficiente, puesto que este no se multiplica ni se modifica.

	print("Finalmente, hallemos la aproximación del valor dado:")

	for (i in k : 2){
		interp = interp + x[i] * a ^ (i-1) ### Calcula el valor de la aproximación.
	}

	print("Con x =")
	print(a)
	print("El resultado es: ")
	print(interp)
}

InterpolacionNewton<-function(X,fx,a){ ### Interpolación de forma de Newton.
	if (length(X)!=length(fx)){ ### Misma validación anterior.
		return("ERROR: Debe ingresar un número equivalente de puntos (x,fx)")
	}

	k = length(X)

	### Para programar esto, consideré más sencillo separar el cálculo de las bases del cálculo de los coeficientes.
	### Para el cálculo de los coeficientes, también consideré más fácil hacer uso de una matriz-cuadro similar al dado en clase.
	### Esta matriz-cuadro va almacenando todos los valores de las diferencias divididas a medida que las calcula.
	### Se puede prescindir de él, sin embargo, considero más sencilla esa forma de realizar el cálculo.

	base = numeric(k) ### Inicia calculando las bases.
	base[1]=1 

	print("Para resolver la interpolación mediante la forma de Newton, hallemos las bases y sus diferencias divididas.")
	print("Primero, hallando las bases tomando el valor de x a aproximar, tenemos que B será igual a:")

	for (i in 2 : k){
		base[i]=(a-X[i-1])*base[i-1] ### {1,(x-x0),(x-x0)(x-x1),(x-x0)(x-x1)(x-x2),...}
	}

	print(base)
	n = 1 ### Indica la "distancia" entre los índices de Xi que se están manipulando.
		  ### También puede considerarse como el número de la iteración actual. Al alcanzar la iteración k-1, se detiene.
		  ### Se utiliza para restar los coeficientes adecuados a la hora de calcular el cociente.

	print("Ahora, procedamos a realizar el proceso para hallar las distintas diferencias divididas.")
	print("Estas diferencias divididas nos servirán como coeficientes del polinomio interpolador.")

	difdiv = matrix(nrow=k,ncol=k) ### Crea la matriz.
	difdiv[,1]=fx ### La matriz-cuadro guarda los valores de fx en la primera columna

	for (j in 2 : k){ ### A partir de la segunda columna...
		for (i in 1 : (k-j+1)){ ### Y a partir de la primera fila, hasta la fila que debe realizar el cálculo
			difdiv[i,j]=(difdiv[i,(j-1)]-difdiv[(i+1),(j-1)])/(X[i]-X[i+n])
			### En este espacio de la matriz, se almacena la diferencia dividida hallada nueva.
			### A destacar la presencia del n en el cociente (X[i]-X[i+n]).
		}
		n=n+1
	}

	print("La matriz estudiada para las diferencias divididas es la siguiente:")
	print(difdiv)

	print("De la matriz anterior, extraemos los coeficientes del polinomio interpolador.")
	print("Estos coeficientes son los siguientes:")
	print(difdiv[1,]) ### Realmente es lo que se deseaba hallar.

	print("Finalmente, hallemos el valor de la aproximación deseada, dado por la sumatoria del producto de las bases y sus coeficientes.")

	interp = 0
	for (i in 1 : k){ ### Realiza la aproximación. Sumatoria de la base por el coeficiente.
		interp = interp + base[i] * difdiv[1,i]
	}

	print("Así, la aproximación de x = ")
	print(a)
	print("Es la siguiente:")
	print(interp)
}

InterpolacionLagrange<-function(X,fx,a){ ### Interpolación en forma de Lagrange.
	if (length(X)!=length(fx)){ ### Validación.
		return("ERROR: Debe ingresar un número equivalente de puntos (x,fx)")
	}

	k = length(X)
	L = numeric(k) ### Vector para almacenar los valores de Li.

	for (i in 1 : k){
		L[i]=1 ### Como se va a utilizar una multiplicación recursiva, se inicializan en 1.
	}

	print("Para resolver la interpolación mediante la forma de Lagrange, utilizaremos los coeficientes L[i]")
	print("Hallaremos estos coeficientes mediante la operación del productorio, que se ejecutará a continuación.")
	print("Esta operación toma en consideración el valor de x a aproximar.")

	for (i in 1 : k){ ### Recorre los Li
		for (j in 1 : k){ ### Recorre los elementos realizando las multiplicaciones convenientes.
			if (i!=j){
				L[i]=L[i]*(a-X[j])/(X[i]-X[j]) ### Toma el valor de x aproximado de inmediato.
			}
		}
	}

	print("Los distintos valores de L[i] hallados son los siguientes:")
	print(L)

	interp = 0

	for (i in 1 : k){ ### Calcula la aproximación como la sumatoria de los Li por sus f(x) correspondientes.
		interp = interp + L[i] * fx[i]
	}

	print("Ahora, hallemos la aproximación de x = ")
	print(a)
	print("La cual será igual a:")
	return(interp)
}
