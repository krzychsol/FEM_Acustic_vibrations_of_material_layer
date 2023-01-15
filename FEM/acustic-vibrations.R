### METODA ELEMENTÓW SKOŃCZONYCH ###
### RÓWNANIE WIBRACJI FAL AKUSTYCZNYCH WARSTY MATERIAŁU ###

start_time <- Sys.time()

# liczba elementów
n <- 50

# funkcje bazowe
e <- function(x, i){
  h <- 2/n
  
  if (x >= h*(i-1) && x <= h*i){
    return ((x-h*(i-1))/h)
  }
  else if (x >= h*i && x <= h*(i+1)){
    return ((h*(i+1) - x)/h)
  }
  else {
    return (0)
  }
}

# pochodna funkcji bazowej
e_prim <- function(x, i){
  h <- 2/n
  
  if (x >= h*(i-1) && x < h*i){
    return (1/h)
  }
  else if (x >= h*i && x < h*(i+1)){
    return (-1/h)
  }
  else {
    return (0)
  }
}

# funkcja rysująca funkcje bazowe
plot_basis <- function(){
  plot(seq(0, 2, 1/(100*n)), mapply(e, seq(0, 2, 1/(100*n)), 1), 
       main = 'Funkcje bazowe',
       xlab='',
       ylab='',
       type='l')
  for (i in 1:n-1){
    lines(seq(0, 2, 1/(10*n)), mapply(e, seq(0, 2, 1/(10*n)), i))
  }
}

#plot_basis()

# całkowanie kwadraturą Gaussa-Legendre’a 
# 2 punkty kwadratury
integrate_B <- function(u,v){
  h <- 2/n
  #u <- u-1
  #v <- v-1
  b = -e(2,u)*e(2,v)
  
  #nad przekatna
  if (u == v-1){
    factor_1 = h/2 * (1/sqrt(3))+(2*u*h+h)/2
    factor_2 = h/2 * (-1/sqrt(3))+(2*u*h+h)/2
    b = b + h/2 * (e_prim(factor_1,u)*e_prim(factor_1,v)
                    +e_prim(factor_2,u)*e_prim(factor_2,v))
    b = b - h/2 * (e(factor_1,u)*e(factor_1,v)
                    +e(factor_2,u)*e(factor_2,v))
  }
  
  #pod przekatna
  else if(u == v+1){
    factor_1 = h/2 * (1/sqrt(3)) + (2*u*h-h)/2;
    factor_2 = h/2 * (-1/sqrt(3)) + (2*u*h-h)/2;
    b = b + h/2 * (e_prim(factor_1,u)*e_prim(factor_1,v)
                    +e_prim(factor_2,u)*e_prim(factor_2,v))
    b = b - h/2 * (e(factor_1,u)*e(factor_1,v)
                    +e(factor_2,u)*e(factor_2,v))
  }
  
  #glowna przekatna
  else if(u == n){
    factor_1 = h/2 * (1/sqrt(3)) + (2*u*h-h)/2;
    factor_2 = h/2 * (-1/sqrt(3)) + (2*u*h-h)/2;
    b = b + h/2 * (e_prim(factor_1,u)^2+e_prim(factor_2,u)^2)
    b = b - h/2 * (e(factor_1,u)^2+e(factor_2,u)^2)
  }
  
  else{
    factor_1 = h * (1/sqrt(3)) + u*h;
    factor_2 = h * (-1/sqrt(3)) + u*h;
    b = b + h * (e_prim(factor_1,u)^2+e_prim(factor_2,u)^2)
    b = b - h * (e(factor_1,u)^2+e(factor_2,u)^2)
  }
  
  return (b)
}

integrate_L <- function(v){
  #v <- v-1
  h <- 2/n
  factor_1 = h*(1/sqrt(3))+v*h
  factor_2 = h*(-1/sqrt(3))+v*h
  l = h*(sin(factor_1)*e(factor_1,v)+sin(factor_2)*e(factor_2,v))
  
  return (l)
}

# glowna funkcja 
solution <- function(){
  b_matrix <- matrix(0,nrow=n,ncol=n)
  
  #z warunku Dirichleta
  b_matrix[1,1] <- 1
  
  l_matrix <- matrix(0,nrow=n,ncol=1)
  
  for (i in 2:n){
    b_matrix[i,i] <- integrate_B(i,i)
    if (i < n){
      b_matrix[i,i+1] <- integrate_B(i,i+1)
    }
    if (i > 2){
      b_matrix[i,i-1] <- integrate_B(i,i-1)
    }
    l_matrix[i,1] <- integrate_L(i)
  }
  
  #print(b_matrix)
  #print(l_matrix)
  
  #rozwiazanie ukladu
  vs = solve(b_matrix,l_matrix)
  
  #kombinacja liniowa funkcji bazowych
  result_function <- function(x,v=vs){
    result = 0
    for (i in 1:n){
      result = result + v[i]*e(x,i)
    }
    return(result)
  }
  
  return (result_function)
}

#funcja rysujaca rozwiazanie 
plot_result <- function(){
  u = solution()
  plot(seq(0,2,1/(100*n)),
       mapply(u,seq(0,2,1/(100*n))),
       main='Rozwiązanie równania wibracji fal akustycznych warsty materiału',
       xlab=' ',
       ylab=' ',
       type='l')
}

plot_result()
end_time <- Sys.time()
print(paste("Czas wykonania: ",round(end_time-start_time,digits=2)," s"))