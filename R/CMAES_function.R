UpdateStepSize <- function(a , successProbability,cp,d,ps_target,minimumStep =0.001){ 
  # a is a 5-tuple, in this case a list with 5 members
  # members of a:
  # a[1] <- a design point x, dim = nVar, column vector
  # a[2] <- average success rate, dim = 1
  # a[3] <- step size, dim = 1, column vector
  # a[4] <- evolution path, dim = nVar, column vector
  # a[5] <- covariance matrix, dim = nVar x nVar
  
  # success probability is a scalar
  newAverageSuccessRate <- (1-cp)*a[[2]] + cp*successProbability
  newStepSize <- a[[3]]* exp( (1/d) * ( a[[2]] - ( ps_target / (1-ps_target) ) * (1-a[[2]]) ) )
  
  if(newStepSize < minimumStep){
    newStepSize <- minimumStep
  }
  a_new <- list(x = a[[1]], 
                averageSuccessRate = newAverageSuccessRate, 
                stepSize = newStepSize, 
                evoPath = a[[4]],
                covarianceMatrix = a[[5]])
  
  return(a_new)
} 

UpdateCovarianceMatrix<- function(a, x_step, ps_threshold,cc,ccov){
  p_success <- a[[2]]
  pc <- a[[4]]
  Cov <- a [[5]]
  if(p_success < ps_threshold){
    newEvoPath <- (1-cc) * pc + ( ( cc * (2-cc) ) ^ 0.5 ) * x_step
    newCov <- (1-ccov) * Cov + ccov * (pc %*% t(pc))
  }else{
    newEvoPath <- (1-cc) * pc 
    newCov <- (1-ccov) * Cov + ccov * ( (pc %*% t(pc)) + (cc * (2-cc) * Cov) )
  }
  
  a_new <- list(x = a[[1]], 
                averageSuccessRate = a[[2]], 
                stepSize = a[[3]], 
                evoPath = newEvoPath,
                covarianceMatrix = newCov)
  
  return(a_new)
}


