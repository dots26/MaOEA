#' The WFG1 test function.
#' @param nObj The number of objective
#' @param individual The individual to be evaluated
#' @param k Number of distance related parameters. The reference suggests a positive integer multiplied by (nObj-1). Default to nObj-1
#' @return A matrix of size nObjective, containing the objective values.
#' @examples
#' individual <- runif(14)
#' nObj <- 4
#' WFG1(individual,nObj)
#' @references Huband, S., Hingston, P., Barone, L., While, L.: A review of multiobjective test problems and a scalable test problem toolkit. Trans. Evol. Comp 10 (5), 477–506 (2006)
#' @export

WFG1 <- function(individual, nObj,k = nObj-1){
  M <- nObj
  n <- length(individual) # number of variables
  l <- n-k

  individual1 <- individual
  individual2 <- individual
  individual3 <- individual
  x <- rep(0,M)
  h <- x
  #individual <- individual/seq(2,2*n,2)
  # first transformation
  for(i in 1:k){
    individual1[i] <- individual[i]
  }
  for(i in (k+1):n){
    individual1[i] <- s_linear(individual[i],0.35)
  }

  # second transform
  for(i in 1:k){
    individual2[i] <- individual1[i]
  }
  for(i in (k+1):n){
    individual2[i] <- b_flat(individual1[i],0.8,0.75,0.85)
  }

  # third transform
  for(i in 1:n){
    individual3[i] <- b_poly(individual2[i],0.02)
  }

  # fourth transform
  for(i in 1:(M-1)){
    rsumMinIndex <- (i-1)*k/(M-1)+1
    rsumMaxIndex <- (i*k)/(M-1)
    weightMin <- 2*rsumMinIndex
    weightMax <- 2*rsumMaxIndex

    weightVector <- seq(weightMin,weightMax,2)

    x[i] <- r_sum(individual3[rsumMinIndex:rsumMaxIndex],weightVector)
  }
  rsumMinIndex <- k+1
  rsumMaxIndex <- n
  weightMin <- 2*rsumMinIndex
  weightMax <- 2*rsumMaxIndex

  weightVector <- seq(weightMin,weightMax,2)

  x[M] <- r_sum(individual3[rsumMinIndex:rsumMaxIndex],weightVector)


  # shape function
  for(i in 1:(M-1)){
    h[i] <- shape_convex(M,i,x)
  }
  h[M] <- shape_mixed(M,x,1,5)

  S <- seq(2,2*M,2)

  obj_val <- x + h*S
  return(obj_val)
}

#' The WFG2 test function.
#' @param nObj The number of objective
#' @param individual The individual to be evaluated
#' @param k Number of distance related parameters. The reference suggests a positive integer multiplied by (nObj-1). Default to nObj-1
#' @return A matrix of size nObjective, containing the objective values.
#' @examples
#' individual <- runif(14)
#' nObj <- 4
#' WFG2(individual,nObj)
#'
#' @references Huband, S., Hingston, P., Barone, L., While, L.: A review of multiobjective test problems and a scalable test problem toolkit. Trans. Evol. Comp 10 (5), 477–506 (2006)
#' @export
WFG2 <- function(individual, nObj,k = nObj-1){
  M <- nObj
  n <- length(individual) # number of variables
  l <- n-k


  individual1 <- individual
  individual2 <- individual
  individual3 <- individual
  x <- rep(0,M)
  h <- x

  #individual <- individual/seq(2,2*n,2)
  # first transformation
  for(i in 1:k){
    individual1[i] <- individual[i]
  }
  for(i in (k+1):n){
    individual1[i] <- s_linear(individual[i],0.35)
  }

  # second transform
  for(i in 1:k){
    individual2[i] <- individual1[i]
  }
  for(i in (k+1):(k+ (l/2) ) ){
    individual2[i] <- r_nonsep(c(individual1[k+2*(i-k)-1],individual1[k+2*(i-k)]),2)
  }

  # third transform
  for(i in 1:(M-1)){
    rsumMinIndex <- (i-1)*k/(M-1)+1
    rsumMaxIndex <- (i*k)/(M-1)
    weightVector <- rep(1,rsumMaxIndex-rsumMinIndex+1)

    x[i] <- r_sum(individual2[rsumMinIndex:rsumMaxIndex],weightVector)
  }
  rsumMinIndex <- k+1
  rsumMaxIndex <- k+l/2

  weightVector <- rep(1,rsumMaxIndex-rsumMinIndex+1)
  x[M] <- r_sum(individual2[rsumMinIndex:rsumMaxIndex],weightVector)


  # shape function
  for(i in 1:(M-1)){
    h[i] <- shape_convex(M,i,x)
  }
  h[M] <- shape_disconnected(M,x,1,1,5)

  S <- seq(2,2*M,2)

  obj_val <- x + h*S
  return(obj_val)
}

#' The WFG4 test function.
#' @param nObj The number of objective
#' @param individual The individual to be evaluated
#' @param k Number of distance related parameters. The reference suggests a positive integer multiplied by (nObj-1). Default to nObj-1
#' @return A matrix of size nObjective, containing the objective values.
#' @examples
#' individual <- runif(14)
#' nObj <- 4
#' WFG4(individual,nObj)
#' @references Huband, S., Hingston, P., Barone, L., While, L.: A review of multiobjective test problems and a scalable test problem toolkit. Trans. Evol. Comp 10 (5), 477–506 (2006)
#' @export
WFG4 <- function(individual, nObj, k = nObj-1){
  M <- nObj
  n <- length(individual) # number of variables
  l <- n-k

  individual1 <- individual
  x <- rep(0,M)
  h <- x
  # first transformation
  for(i in 1:n){
    individual1[i] <- s_multi(individual[i],30,10,0.35)
  }

  # second transform
  for ( i in 1:(M-1)){
    rsumMinIndex <- (i-1)*k/(M-1)+1
    rsumMaxIndex <- (i*k)/(M-1)

    weightVector <- rep(1,rsumMaxIndex-rsumMinIndex+1)

    x[i] <- r_sum(individual1[rsumMinIndex:rsumMaxIndex],weightVector)
  }
  rsumMinIndex <- k+1
  rsumMaxIndex <- n

  weightVector <- rep(1,rsumMaxIndex-rsumMinIndex+1)

  x[M] <- r_sum(individual1[rsumMinIndex:rsumMaxIndex],weightVector)
  # shape function
  for(i in 1:M){
    h[i] <- shape_concave(M,i,x)
  }
  S <- seq(2,2*M,2)

  obj_val <- x[M] + h*S
  return(obj_val)
}

#' The WFG5 test function.
#' @param nObj The number of objective
#' @param individual The individual to be evaluated
#' @param k Number of distance related parameters. The reference suggests a positive integer multiplied by (nObj-1). Default to nObj-1
#' @return A matrix of size nObjective, containing the objective values.
#' @examples
#' individual <- runif(14)
#' nObj <- 4
#' WFG5(individual,nObj)
#' @references Huband, S., Hingston, P., Barone, L., While, L.: A review of multiobjective test problems and a scalable test problem toolkit. Trans. Evol. Comp 10 (5), 477–506 (2006)
#' @export
WFG5 <- function(individual, nObj,k = nObj-1){
  M <- nObj
  n <- length(individual) # number of variables
  l <- n-k
 # individual <- individual/seq(2,2*n,2)
  individual1 <- individual
  x <- rep(0,M)
  h <- x
  # first transformation
  for(i in 1:n){
    individual1[i] <- s_deceptive(individual[i],0.35,0.001,0.05)
  }

  # second transform
  for ( i in 1:(M-1)){
    rsumMinIndex <- (i-1)*k/(M-1)+1
    rsumMaxIndex <- (i*k)/(M-1)

    weightVector <- rep(1,rsumMaxIndex-rsumMinIndex+1)

    x[i] <- r_sum(individual1[rsumMinIndex:rsumMaxIndex],weightVector)
  }

  rsumMinIndex <- k+1
  rsumMaxIndex <- n

  weightVector <- rep(1,rsumMaxIndex-rsumMinIndex+1)

  x[M] <- r_sum(individual1[rsumMinIndex:rsumMaxIndex],weightVector)
  # shape function
  for(i in 1:M){
    h[i] <- shape_concave(M,i,x)
  }
  S <- seq(2,2*M,2)

  obj_val <- x[M] + h*S
  return(obj_val)
}

#' The WFG6 test function.
#' @param nObj The number of objective
#' @param individual The individual to be evaluated
#' @param k Number of distance related parameters. The reference suggests a positive integer multiplied by (nObj-1). Default to nObj-1
#' @return A matrix of size nObjective, containing the objective values.
#' @examples
#' individual <- runif(14)
#' nObj <- 4
#' WFG6(individual,nObj)
#' @references Huband, S., Hingston, P., Barone, L., While, L.: A review of multiobjective test problems and a scalable test problem toolkit. Trans. Evol. Comp 10 (5), 477–506 (2006)
#' @export
WFG6 <- function(individual, nObj,k = nObj-1){
  M <- nObj
  n <- length(individual) # number of variables
  l <- n-k

 # individual <- individual/seq(2,2*n,2)
  individual1 <- individual
  x <- rep(0,M)
  h <- x
  # first transformation
  for(i in 1:k){
    individual1[i] <- individual[i]
  }
  for(i in (k+1):n){
    individual1[i] <- s_linear(individual[i],0.35)
  }

  # second transform
  for(i in 1:(M-1)){
    rsumMinIndex <- (i-1)*k/(M-1)+1
    rsumMaxIndex <- (i*k)/(M-1)

    x[i] <- r_nonsep(individual1[rsumMinIndex:rsumMaxIndex],k/(M-1))
  }
  rsumMinIndex <- k+1
  rsumMaxIndex <- n

  x[M] <- r_nonsep(individual1[rsumMinIndex:rsumMaxIndex],l)
  # shape function
  for(i in 1:M){
    h[i] <- shape_concave(M,i,x)
  }
  S <- seq(2,2*M,2)

  obj_val <- x[M] + h*S
  return(obj_val)
}

#' The WFG7 test function.
#' @param nObj The number of objective
#' @param individual The individual to be evaluated
#' @param k Number of distance related parameters. The reference suggests a positive integer multiplied by (nObj-1). Default to nObj-1
#' @return A matrix of size nObjective, containing the objective values.
#' @examples
#' individual <- runif(14)
#' nObj <- 4
#' WFG7(individual,nObj)
#' @references Huband, S., Hingston, P., Barone, L., While, L.: A review of multiobjective test problems and a scalable test problem toolkit. Trans. Evol. Comp 10 (5), 477–506 (2006)
#' @export
WFG7 <- function(individual, nObj,k = nObj-1){
  M <- nObj
  n <- length(individual) # number of variables
  l <- n-k
  #individual <- individual/seq(2,2*n,2)
  individual1 <- individual
  individual2 <- individual

  x <- rep(0,M)
  h <- x
  # first transformation
  for(i in 1:k){
    individual1[i] <- b_param(individual[i],r_sum(individual[(i+1):n],rep(1,(n-i))),0.98/49.98,0.02,50)
  }
  for(i in (k+1):n){
    individual1[i] <- (individual[i])
  }

  # second transform
  for(i in 1:k){
    individual2[i] <- individual1[i]
  }
  for(i in (k+1):n){
    individual2[i] <- s_linear(individual1[i],0.35)
  }


  for ( i in 1:(M-1)){
    rsumMinIndex <- (i-1)*k/(M-1)+1
    rsumMaxIndex <- (i*k)/(M-1)

    weightVector <- rep(1,rsumMaxIndex-rsumMinIndex+1)

    x[i] <- r_sum(individual2[rsumMinIndex:rsumMaxIndex],weightVector)
  }

  rsumMinIndex <- k+1
  rsumMaxIndex <- n

  weightVector <- rep(1,rsumMaxIndex-rsumMinIndex+1)

  x[M] <- r_sum(individual2[rsumMinIndex:rsumMaxIndex],weightVector)
  # shape function
  for(i in 1:M){
    h[i] <- shape_concave(M,i,x)
 }
  S <- seq(2,2*M,2)


  obj_val <- x[M] + h*S
  return(obj_val)
}

#' The WFG8 test function.
#' @param nObj The number of objective
#' @param individual The individual to be evaluated
#' @param k Number of distance related parameters. The reference suggests a positive integer multiplied by (nObj-1). Default to nObj-1
#' @return A matrix of size nObjective, containing the objective values.
#' @examples
#' individual <- runif(14)
#' nObj <- 4
#' WFG8(individual,nObj)
#' @references Huband, S., Hingston, P., Barone, L., While, L.: A review of multiobjective test problems and a scalable test problem toolkit. Trans. Evol. Comp 10 (5), 477–506 (2006)
#' @export
WFG8 <- function(individual, nObj,k = nObj-1){
  M <- nObj
  n <- length(individual) # number of variables
  l <- n-k
#  individual <- individual/seq(2,2*n,2)
  individual1 <- individual
  x <- rep(0,M)
  h <- x
  # first transformation
  for(i in 1:k){
    individual1[i] <- individual[i]
  }
  for(i in (k+1):n){
    individual1[i] <- s_linear(individual[i],0.35)
  }


  # second transform
  for ( i in 1:(M-1)){
    rsumMinIndex <- (i-1)*k/(M-1)+1
    rsumMaxIndex <- (i*k)/(M-1)

    weightVector <- rep(1,rsumMaxIndex-rsumMinIndex+1)

    x[i] <- r_sum(individual1[rsumMinIndex:rsumMaxIndex],weightVector)
  }

  rsumMinIndex <- k+1
  rsumMaxIndex <- n

  weightVector <- rep(1,rsumMaxIndex-rsumMinIndex+1)

  x[M] <- r_sum(individual1[rsumMinIndex:rsumMaxIndex],weightVector)
  # shape function
  for(i in 1:M){
    h[i] <- shape_concave(M,i,x)
  }
  S <- seq(2,2*M,2)

  obj_val <- x[M] + h*S
  return(obj_val)
}

#' The WFG9 test function.
#' @param nObj The number of objective
#' @param individual The individual to be evaluated
#' @param k Number of distance related parameters. The reference suggests a positive integer multiplied by (nObj-1). Default to nObj-1
#' @return A matrix of size nObjective, containing the objective values.
#' @examples
#' individual <- runif(14)
#' nObj <- 4
#' WFG9(individual,nObj)
#' @references Huband, S., Hingston, P., Barone, L., While, L.: A review of multiobjective test problems and a scalable test problem toolkit. Trans. Evol. Comp 10 (5), 477–506 (2006)
#' @export
WFG9 <- function(individual, nObj,k = nObj-1){
  M <- nObj
  n <- length(individual) # number of variables
  l <- n-k
#  individual <- individual/seq(2,2*n,2)
  individual1 <- individual
  x <- rep(0,M)
  h <- x
  # first transformation
  for(i in 1:(n-1)){
    individual1[i] <- b_param(individual[i],r_sum(individual[(i+1):n],rep(1,(n-i))),0.98/49.98,0.02,50)
  }
  individual1[n] <- (individual[n])


  # second transform
  for ( i in 1:(k)){
    s_deceptive(individual1[i],0.35,0.001,0.05)
  }
  for ( i in (k+1):(n)){
    s_multi(individual1[i],30,95,0.35)
  }

  # third transform
  for(i in 1:(M-1)){
    rsumMinIndex <- (i-1)*k/(M-1)+1
    rsumMaxIndex <- (i*k)/(M-1)

    x[i] <- r_nonsep(individual1[rsumMinIndex:rsumMaxIndex],k/(M-1))
  }
  rsumMinIndex <- k+1
  rsumMaxIndex <- n

  x[M] <- r_nonsep(individual1[rsumMinIndex:rsumMaxIndex],l)
  # shape function
  for(i in 1:M){
    h[i] <- shape_concave(M,i,x)
  }
  S <- seq(2,2*M,2)

  obj_val <- x[M] + h*S
  return(obj_val)
}
