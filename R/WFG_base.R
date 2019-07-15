shape_linear <- function(nObj, objectiveIndex, x){
  multiplication <- 1

  for(varIndex in 1:(nObj-1)){
    if(objectiveIndex == 1){
      multiplication <- multiplication * x[varIndex]
    }

    if((objectiveIndex < nObj) && (objectiveIndex > 1)){
      if(varIndex <= (nObj-objectiveIndex))
        multiplication <- multiplication * x[varIndex]
      if(varIndex == (nObj-1))
        multiplication <- multiplication * (1-x[nObj-objectiveIndex+1])
    }

    if(objectiveIndex == nObj){
      multiplication <- multiplication * (1-x[1])
    }
  }

  return(multiplication)
}

shape_convex <- function(nObj, objectiveIndex, x){
  multiplication <- 1

  for(varIndex in 1:(nObj-1)){
    if(objectiveIndex == 1){
      multiplication <- multiplication * (1-cos(x[varIndex]*pi/2))
    }

    if((objectiveIndex < nObj) && (objectiveIndex > 1)){
      if(varIndex <= (nObj-objectiveIndex))
        multiplication <- multiplication * (1-cos(x[varIndex]*pi/2))
      if(varIndex == (nObj-1))
        multiplication <- multiplication * (1-sin(x[nObj-objectiveIndex+1]*pi/2))
    }

    if(objectiveIndex == nObj){
      multiplication <- multiplication * (1-sin(x[1]*pi/2))
    }
  }

  return(multiplication)
}

shape_concave <- function(nObj, objectiveIndex, x){
  multiplication <- 1

  for(varIndex in 1:(nObj-1)){
    if(objectiveIndex == 1){
      multiplication <- multiplication * (sin(x[varIndex]*pi/2))
    }

    if((objectiveIndex < nObj) && (objectiveIndex >1)){
      if(varIndex <= (nObj-objectiveIndex))
        multiplication <- multiplication * (sin(x[varIndex]*pi/2))
      if(varIndex == (nObj-1))
        multiplication <- multiplication * (cos(x[nObj-objectiveIndex+1]*pi/2))
    }

    if(objectiveIndex == nObj){
      multiplication <- multiplication * (cos(x[1]*pi/2))
      break
    }
  }

  return(multiplication)
}

shape_mixed<- function(nObj, x, alpha, A){
  obj <- (1 - x[1] - ((cos((2*A*pi*x[1])+(pi/2)))/(2*pi*A)))^alpha
}

shape_disconnected <- function(nObj, x, alpha, beta, A){
  obj <- 1 - (x[1]^alpha) * cos(A*(x[1]^beta)*pi)*cos(A*(x[1]^beta)*pi)
}

b_poly <- function(y,alpha){
  x <- y^alpha
  return(x)
}

b_flat <- function(y,A,B,C){
  x <- A + min(c(0,floor(y-B))*A*(B-y)/B) + min(c(0,floor(C-y))*(1-A)*(y-C)/(1-C))
  return(x)
}

b_param <- function(y,secondary_y,A,B,C){
  u <- secondary_y
  v <- A - (1 - 2*u) * abs(floor(0.5-u) + A)
  x <- y^(B+(C-B)*v)
  return(x)
}

s_linear <- function(y,A){
  x <- abs(y-A)/abs(floor(A-y)+A)
  return(x)
}

s_deceptive <- function(y,A,B,C){
  bracket_term1 <- floor(y-A+B)*(1-C+((A-B)/B))/(A-B)
  bracket_term2 <- floor(A+B-y)*(1-C+((1-A-B)/B))/(1-A-B)
  brackets <- bracket_term1 + bracket_term2 + (1/B)

  x <- 1 + ((abs(y-A)-B) * brackets)
}

s_multi <- function(y,A,B,C){
  first_term <- cos((4*A)+2*pi*(0.5-(abs(y-C)/(2*(floor(C-y)+C)))))
  second_term <- 4*B * (abs(y-C)/(2*(floor(C-y)+C))) * (abs(y-C)/(2*(floor(C-y)+C)))
  nominator <- 1 + first_term + second_term
  denominator <- B + 2
  x <- nominator/denominator
}

r_sum <- function(y, weight){
  x <- sum(y*weight)/sum (weight)
  return(x)
}

r_nonsep <- function(y, A){
  nVar <- length(y)
  nominator <- 0
  for(varIndex in 1:nVar){
    subsum <- 0
    for(k in 0:(A-2)){
      if((A-2) > 0){
        subsum <- subsum + abs(y[varIndex]-y[1+((varIndex+k)%%nVar)])
      }
    }
    nominator <- nominator + y[varIndex] + subsum
  }
  denominator <- nVar/A * ceiling(A/2) * (1+ 2*A - 2*ceiling(A/2))
  x <- nominator/denominator
  return(x)
}
