densite  = function(x) {
  return(dbeta(x, 1.7, 5.1))
}

# Les deux fonctions suivantes sont écrites différements dans le projet, mais voici ce qu'elles donnaient au TP
dividifTP=function(x,y){
  ##  Newton's Divided differences
  n = length(x) - 1 ## degree of the Lagrange polynomial
  d  = y 
  for (j in 2:(n+1) ) {
    d[j : (n+1) ] = (d[j : (n+1)] - d[(j-1):n])/(x[j:(n+1)] - x[1:(n+2-j)])
  }
  return(d)    
}


hornerNewtonTP = function(a,x,z){
  ## Horner's method: Evaluates  a polynom P at points z, given
  n = length(x);
  f  = a[n] * (z-x[n-1]) + a[n-1]
  for( i in 2:(n-1)){
    f = a[n-i] + f  * (z-x[n-i])
  }
  return(f)
}



interpolDividif=function(x,y,z){
  ## Efficient Lagrange interpolation using Horner's method with  
  ## Newton basis for evaluation
  coeff = dividif(x,y)
  return (hornerNewton(coeff, x, z))
}

NTcheby = function(n, a, b){
  T = seq(0, n-1);
  for (k in (1:n)){
    T[k]= ((a+b)/2) + ((a-b)/2)*cos(((T[k]+(1/2))*pi)/n)}
  return (T)
}

# Je l'ai modifiée par rapport au tp afin de pouvoir superposer des courbes dans le projet
interpolLagrange =function(n, a, b, neval, nodes = 'equi', FUN){ 
  ## Generic Lagrange interpolation, with equidistant or Chebyshev nodes. 
  ## @param n : the degree of the interpolating polynomial on each
  ## subinterval
  ## @param a : left end-point of the interval
  ## @param b : right end-point of the interval
  ## @param neval :number of evaluation points (a regular grid will be
  ## used on [a,b]
  ## @param nodes :string, either "equi" (default) for equidistant
  ## Lagrange interpolation (on each subinterval) or "cheby" for
  ## using Chebyshev nodes.
  ## @param FUN: the function to be interpolated 
  ## @return : vector of size 'neval': the values of the Lagrange
  ## polynomial on an equi-distant grid.
  
  if (nodes == 'equi'){
    x =  seq(a,b,length.out=n)  
  }
  else if (nodes == 'cheby'){
    x = NTcheby(n,a,b)
  }
  else{stop("the nodes must be either 'equi' or 'cheby'") }
  
  z = seq(a,b,length.out=neval)
  y = FUN(x)
  f = interpolDividif(x, y, z)
  
  return(f)              
}


dividif=function(x,y){ #
  ##  Newton's Divided differences
  n = length(x) - 1 ## degree of the Lagrange polynomial
  d  = y 
  if(n==0){
    d = d/x
  } else {
    for (j in 2:(n+1) ) {d[j : (n+1) ] = (d[j : (n+1)] - d[(j-1):n])/(x[j:(n+1)] - x[1:(n+2-j)])}
  }
  return(d)    
}

hornerNewton = function(a,x,z){
  ## Horner's method: Evaluates a polynom P at points z, given ## nodes x and the coefficients a of P in Newton's basis
  ##
  ## @param a : vector: the coefficients of the polynomial in ## Newton's basis
  ## @param x : the interpolation nodes.
  ## @param z : vector of points where the polynom needs to be ## evaluated.
  ## @return : a vector of same size as z: the value of the
  ## polynomial at points z.
  ##
  n = length(x)
  f = a[n]*rep(1,length(z))
  if(n==0){stop('at least one interpolating point is needed')} 
  if(n==1){
    f = f*(z-x[1]) 
  }
  if(n >= 2){
    for( i in 1:(n-1)){
      f = f * (z-x[n-i]) + a[n-i]
    } }
  return(f) }


piecewiseInterpol=function(n,nInt,a,b,neval, nodes = 'equi', FUN, Plot) {
  ## @param n : the degree of the interpolating polynomial on each
  ## subinterval
  ## @param nInt :  the number of sub-intervals
  ## @param a, b : endpoints of the interval
  ## @param neval : the number of points on the interpolating grid (on
  ## each subinterval)
  ## @param nodes : string, either "equi" (default) for equidistant
  ## Lagrange interpolation (on each subinterval) or "cheby" for
  ## chebyshev nodes.
  ## @param FUN the function to be interpolated
  ## @param Plot : logical. Should the result be plotted ?
  ## @return : a matrix with 2 rows and neval * nInt -neval + 1:
  ## values of the interpolated funtion on a regular grid (first row)
  ## and the corresponding abscissas (second row).
  
  intEndPoints = seq(a,b,length.out = nInt+1)
  f = c()
  z = c()
  for (m in 1:nInt){
    A = intEndPoints[m]; B = intEndPoints[m+1] 
    
    fm = interpolLagrange(n, A, B, neval, nodes, FUN)
    zm = seq(A,B,length.out=length(fm)) 
    
    if( m >= 2){
      zm = zm[2:length(fm)];
      fm = fm[2:length(fm)];
      z = c(z,zm);
      f = c(f,fm);
    }
  } 
  if (Plot == 0) {
    return(f)
  } else {
    if (nodes == "equi") {methodName = " equidistant "}
    else  {methodName = " Chebyshev "}
    
    
    plot(z, sapply(z,FUN),type="l");
    title(main = paste("Piecewise Lagrange interpolation with",
                       toString(n+1), methodName, "nodes on",
                       toString(nInt), "Intervals", sep=" "));
    lines(z,f, col='red', lwd=2);
    legend('topright', legend = c('function','interpolation'),
           lwd=c(1,2), col=c('black','red'));
  }}


erreur = function(f,FUN){
  z = seq(0,1,length.out = 1000)
  y = FUN(z)
  erreur = 0
  for (k in(length(z))){
    erreur = c(erreur, abs(y[k]-f[k]))
  }
  return (max(erreur))
}


funHisto = function(x){ #On transforme hist_value en une fonction de x afin de pouvoir réutiliser le code précédent.
  return(hist_value(histog,x))
}


trapezeInt =function(FUN,a,b,M){
  ##' TRAPEZOIDAL INTEGRATION RULE (COMPOSITE)
  ##' @param FUN : the function to be integrated
  ##' @param a, b : interval end points 
  ##' @param M : number of intervals (each of size (b-a)/M)
  ##' @return: the value of the composite trapezoidal quadrature. 
  x = seq(a,b, length.out= M+1)
  y = sapply(x, FUN)
  h = (b-a)/M 
  q = h*sum(y)-(h/2)*(y[1]+y[M+1]) 
  return(q)
}


refineTrapeze=function(FUN,a,b,M,q){
  ##' refinement of the subdivision step: incremental method
  ##' @param FUN : the function to be integrated
  ##' @param a, b : interval end points 
  ##' @param M : initial number of intervals (each of size (b-a)/M)
  ##'  having been used to compute q
  ##' @param  q : the value of the trapezoidal  quadrature method
  ##'  of stepsize (b-a)/M
  ##' @return : the value of the quadrature for a stepsize h' = h/2
  h = (b-a)/M
  x = seq(a+h/2,b-h/2, length.out = M)
  ##  x : a vector of size M :
  ##     the additional abscissas where 'fun' must be evaluated.
  y = sapply(x, FUN)
  Q = q/2 + (h/2)*sum(y) 
  return(Q)
}

simpsonInt = function(FUN,a,b,M){
  ##' Simpson integration via trapeze rule
  ##' uses the fact that 
  ##' simpson(h) = 4/3(trapeze(h/2) - 1/4 trapeze(h))
  h = (b-a)/M;
  qtrapeze = trapezeInt(FUN,a,b,M) 
  qrefined = refineTrapeze(FUN,a,b,M,qtrapeze)
  q =  (4/3)*(qrefined) - ((1/3)*qtrapeze)
  return(q)
}

# version projet 
evalErrSimpson=function(FUN,a,b,M){
  ## Computes an approximation E of the error 
  ## for the composite Simpson rule of step h=(b-a)/(2M). 
  ##This requires computing I_M and I_{2M}. 
  ##The value  q = I_{2M} is also returned. 
  qth = trapezeInt(FUN,a,b,M)   ## M +1 evaluations
  qth2 = refineTrapeze ( FUN,a,b,M,qth )  ## M evaluations
  qth4 = refineTrapeze ( FUN,a,b,2*M,qth2 )   ## 2M evaluations
  simps_h =   4/3*(qth2 - 1/4* qth ) 
  simps_h2 =  4/3*(qth4 - 1/4* qth2 )
  q = simps_h2  
  E = (simps_h - q)/15
  return(E)
}

#version TP
evalErrSimpsonTP=function(FUN,a,b,M){
  ## Computes an approximation E of the error 
  ## for the composite Simpson rule of step h=(b-a)/(2M). 
  ##This requires computing I_M and I_{2M}. 
  ##The value  q = I_{2M} is also returned. 
  qth = trapezeInt(FUN,a,b,M)   ## M +1 evaluations
  qth2 = refineTrapeze ( FUN,a,b,M,qth )  ## M evaluations
  qth4 = refineTrapeze ( FUN,a,b,2*M,qth2 )   ## 2M evaluations
  simps_h =   4/3*(qth2 - 1/4* qth ) 
  simps_h2 =  4/3*(qth4 - 1/4* qth2 )
  q = simps_h2  
  E = (simps_h - q)/15
  return(q, E)
}


quadrature = function(FUN,a,b,M){
  qth = trapezeInt(FUN,a,b,M)   ## M +1 evaluations
  qth2 = refineTrapeze ( FUN,a,b,M,qth )  ## M evaluations
  qth4 = refineTrapeze ( FUN,a,b,2*M,qth2 )   ## 2M evaluations
  simps_h =   4/3*(qth2 - 1/4* qth ) 
  simps_h2 =  4/3*(qth4 - 1/4* qth2 )
  return (simps_h2)
}


choixDuM = function(FUN, a, b, tolerance){
  M = 10 ; erreur = 1 ; Mvecteur = 0 ; Qvecteur = 0 ; Evecteur = 0 ;
  while ( 2*erreur > tolerance) {
    a = evalErrSimpson(FUN,a,b,2*M)
    quadrature = quadrature(FUN, a, b, 2*M)
    Mvecteur = c(Mvecteur, M) ; Qvecteur = c(Qvecteur, quadrature) ;   Evecteur = c(Evecteur, erreur)
    erreur = a
    M = 2 * M
  }
  return (M)
}


densite_bruitee = function(x) {
  return (funHisto(x))
}


romberg =function(FUN,n,a,b,M){## methode de Romberg avec n etapes
  ## appliquee sur la fonction FUN sur l'intervalle (a,b), avec un
  ## pas initial h = (b-a)/M
  h= (b-a)/M 
  A = rep(0, n+1)
  A[1] = trapezeInt(FUN,a,b,M);
  Mc = M
  ## initialisation des differences divisees
  for( i in 2:(n+1)){
    A[i] = refineTrapeze( FUN,a,b, Mc, q= A[i-1])
    Mc = 2*Mc 
  }
  delta = 1/4;
  for (j in 2:(n+1)){
    A[j : (n+1) ] = (A[j:(n+1)] - (delta^(j-1))*A[(j-1):n]) / (1 - delta^(j-1))
  }
  return(A)
}

# richardson (seule fonction non utiliée dans le projet)
richardson = function(FUN,n,t,delta){
  ## Calcule le tableau des differences  divisees en 0 du 
  ## polynome d'interpolation en t,delta t, ... delta^n t
  ## renvoie un vecteur de taille n+1:
  ## le vecteur des A_{k,k}, k= 0 .. n 
  ## (pas la matrice).   
  ## La meilleure approximation est le dernier element A[n+1].
  ##
  lx = log(t)  +  log(delta) *(0:n)
  x = exp(lx) 
  A = sapply(x,FUN) 
  for( j in 2:(n+1)){
    A[j : (n+1) ] =  (A[j:(n+1)] - delta^(j)*A[(j-1):(n)])/(1-delta^(j))## Completer le code 
  }
  return(A)
}
