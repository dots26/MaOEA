#' Do an iteration of population based Multi-Objective Covariance Matrix Adaptation Evolution Strategy (MO-CMA-ES). The variation is using simulated binary crossover (SBX) and mutation following the CMA. The original MO-CMA-ES does not use crossover, to do this simply set crossoverProbability to zero.
#' @title Multi-Objective CMA-ES
#' @param parent The parent generation, an object of class cmaes_gen. The MO-CMA-ES parent is a 5 tuple: x (the design point, length = number of variable),averageSuccessRate (scalar),stepSize (scalar), evoPath (evolution path, vector, length = number of variable ),covarianceMatrix (square matrix with ncol = nrow = number of variable). The parent is then should be a vector of lists (see example).
#' @param nObjective The number of objective functions. A scalar value.
#' @param fun Objective function being solved.
#' @param control List of parameters for CMA-ES. Available control are as follows:
#' \code{successProbTarget} Target success probability
#' \code{successProbThreshold} The threshold for success probability. If the average success probability is higher than this value, the success rate growth is slowed.
#' \code{crossoverProbability} The probability of doing crossover. Should be between 0-1. Negative value will behave like a zero, and values larger than 1 will behave like 1. Default to 1.
#' \code{crossoverDistribution} The distribution index for SBX. Larger index makes the distribution sharper around each parent.
#' @param ... Further arguments to be passed to \code{fun}
#' @return Returns a list for the next generation. It contains list$new_generation (class: cmaes_gen), list$population (basically a copy of list$new_generation[[]]$x), and list$populationObjective
#' @references Voß, T., Hansen, N., Igel, C.: Improved step size adaptation for the MO-CMA-ES. In: Genetic and Evolutionary Computation (GECCO). pp. 487–494. ACM, New York, NY (2010)
#' @examples
#' \donttest{
#' nVar <- 14
#' nObjective <- 5
#' nIndividual <- 100
#' crossoverProbability <- 1
#' ps_target <- 1 / (5 + ( 1 / 2  )^0.5 )
#' pop <- matrix(stats::runif(nIndividual*nVar), nrow = nVar) # create the population
#' a_list <- cmaes_gen(pop)
#' control <- list(successProbTarget=ps_target,crossoverProbability=crossoverProbability)
#'
#' # run a generation of MO-CMA-ES with standard WFG8 test function.
#' numpyready <- reticulate::py_module_available('numpy')
#' pygmoready <- reticulate::py_module_available('pygmo')
#' py_module_ready <- numpyready && pygmoready
#' if(py_module_ready) # prevent error on testing the example
#' newGeneration <- MOCMAES(a_list,nObjective,WFG8,control,nObjective)
#' }
#' @export
MOCMAES <- function( parent, nObjective, fun, control=list(), ...){
  if (!pkg.globals$have_pygmo)
    pkg.globals$have_pygmo <- reticulate::py_module_available("pygmo")
  if (!pkg.globals$have_pygmo)
    stop("MOCMAES requires PyGMO to compute hypervolume")

  con <- list(successProbTarget = 1 / ( 5 + ( 1 / 2  ) ^ 0.5 ),
              successProbThreshold = 0.44,
              crossoverProbability=1,
              crossoverDistribution=30,
              hypervolumeMethod="exact",
              hypervolumeMethodParam=list(),
              referencePoint = NULL)
  con[names(control)] <- control
  control <- con
  if (identical(fun, DTLZ4))
    message("DTLZ4 may suffer from floating-point inaccuracy.")

  ps_threshold <- control$successProbThreshold
  ps_target <- control$successProbTarget

  chromosomeLength <- length(parent[[1]]$x)
  populationSize <- length(parent)

  d <- 1 + chromosomeLength / 2
  cc <- ( 2 / (chromosomeLength + 2))
  ccov <-  2 / ( (chromosomeLength ^ 2) + 6)
  cp <- ps_target * 1 / (2 + ps_target * 1)

  # object of class cmaes_gen
  a_list <- parent

  population<-matrix(, nrow=chromosomeLength, ncol=0);
  for (populationIndex in 1:length(parent)){
    population <- cbind(population, parent[[populationIndex]]$x)
  }

  populationObjective <- matrix(, nrow=nObjective,ncol=0);
  #evaluation of parents
  for (parentIndex in 1:populationSize){
    individualObjectiveValue <- EvaluateIndividual(population[, parentIndex], fun, ...)
    populationObjective <- cbind(populationObjective, individualObjectiveValue)
  }

  new_a_list <- a_list
  new_population <- population
  newObjectiveValue <- populationObjective

  # not used in the reference
  if (stats::runif(1) < control$crossoverProbability){
    parentCount <- 0
    parent_x <- population
    parent_stepSizeAverage <- matrix(, nrow=1)

    while (parentCount<ncol(population)){
      parent_index <- sample.int(ncol(population), 2, replace = FALSE)
      parent_x[,(parentCount+1)] <- a_list[[parent_index[1]]]$x
      parent_x[,(parentCount+2)] <- a_list[[parent_index[2]]]$x

      # step size
      parent_stepSize1 <- a_list[[parent_index[1]]]$stepSize
      parent_stepSize2 <- a_list[[parent_index[2]]]$stepSize

      new_a_list[[parent_index[1]]]$stepSize <- (parent_stepSize1+parent_stepSize2)/2
      new_a_list[[parent_index[2]]]$stepSize <- (parent_stepSize1+parent_stepSize2)/2

      # covariance matrix
      parent_covmatrix1 <- a_list[[parent_index[1]]]$covarianceMatrix
      parent_covmatrix2 <- a_list[[parent_index[2]]]$covarianceMatrix

      new_a_list[[parent_index[1]]]$covarianceMatrix <- (parent_covmatrix1+parent_covmatrix2)/2
      new_a_list[[parent_index[2]]]$covarianceMatrix <- (parent_covmatrix1+parent_covmatrix2)/2

      # evoPath
      parent_evoPath1 <- a_list[[parent_index[1]]]$evoPath
      parent_evoPath2 <- a_list[[parent_index[2]]]$evoPath

      new_a_list[[parent_index[1]]]$evoPath<- (parent_evoPath1+parent_evoPath2)/2
      new_a_list[[parent_index[2]]]$evoPath<- (parent_evoPath1+parent_evoPath2)/2

      # successRate
      parent_sR1 <- a_list[[parent_index[1]]]$averageSuccessRate
      parent_sR2 <- a_list[[parent_index[2]]]$averageSuccessRate

      new_a_list[[parent_index[1]]]$averageSuccessRate <- (parent_sR1 +parent_sR2 )/2
      new_a_list[[parent_index[2]]]$averageSuccessRate <- (parent_sR1 +parent_sR2 )/2

      parentCount<- parentCount +2
    }

    #recombination
    offspring <- nsga2R::boundedSBXover(t(parent_x),
                                        rep(0,chromosomeLength),
                                        rep(1,chromosomeLength),
                                        control$crossoverProbability,
                                        control$crossoverDistribution)
    offspring <- t(offspring)
    # update the population"s x
    for (k in 1:populationSize){
      new_a_list[[k]]$x <- offspring[,k]
    }
  }

  # variation based on multivariate normal distribution (CMA mutation)
  for (k in 1:populationSize){
    individual <- new_a_list[[k]]$x
    covariance <- (new_a_list[[k]]$stepSize^2) *  new_a_list[[k]]$covarianceMatrix
    new_a_list[[k]]$x <- MASS::mvrnorm(1, mu = individual, Sigma = covariance)

    new_population[,k] <- new_a_list[[k]]$x
    newObjectiveValue[,k] <- EvaluateIndividual(new_a_list[[k]]$x,fun,...)
  }

  # non-dominated sorting
  combinedPopulation <- cbind(population,new_population)
  combinedObjectiveValue <- cbind(populationObjective,newObjectiveValue)
  rnk <- nsga2R::fastNonDominatedSorting(t(combinedObjectiveValue))
  rnkIndex <- integer(ncol(combinedPopulation))
  sortingIndex <- 1;
  while (sortingIndex <= length(rnk)) {
    rnkIndex[rnk[[sortingIndex]]] <- sortingIndex;
    sortingIndex <- sortingIndex + 1;
  }

  # selection
  newPopulation <- matrix(, nrow = chromosomeLength,ncol=0)
  newPopulationObjective <- matrix(, nrow = nObjective,ncol=0)
  newPopulationIndex <- matrix(, nrow = 1,ncol=0)
  newPopulationSize <- integer(1)

  frontIndex <- 1
  while (newPopulationSize< populationSize) #new population is filled by every next front
  {
    newPopulationIndex <- append(newPopulationIndex,rnk[[frontIndex]])
    newPopulation <- cbind( newPopulation,combinedPopulation[ , rnk[[frontIndex]] ] )
    newPopulationObjective <- cbind(newPopulationObjective,
                                    combinedObjectiveValue[ , rnk[[frontIndex]] ])
    newPopulationSize <- newPopulationSize + length( rnk[[frontIndex]] )

    frontIndex <- frontIndex+1
  }

  while (newPopulationSize > populationSize) #new population is overcapacity remove least contributor
  {
    leastContrib <- GetLeastContributor(newPopulationObjective,
                                        control$referencePoint,
                                        control$hypervolumeMethod,
                                        control$hypervolumeMethodParam)
    newPopulationIndex <- newPopulationIndex[-leastContrib]
    newPopulation <- newPopulation[,-leastContrib]
    newPopulationObjective <- newPopulationObjective[,-leastContrib]
    newPopulationSize <- newPopulationSize -1
  }


  #update step size and cov matrix
  for (k in 1:populationSize){
    successfullUpdate <- max( as.integer ( newPopulationIndex==(k+populationSize) ) )

    # update step size 1
    a_list[[k]] <- UpdateStepSize(a_list[[k]],successfullUpdate,cp,d,ps_target)
    # update step size 2
    new_a_list[[k]] <- UpdateStepSize(new_a_list[[k]],successfullUpdate,cp,d,ps_target)

    # update cov matrix
    x_step <- (new_a_list[[k]]$x-a_list[[k]]$x) / a_list[[k]]$stepSize
    new_a_list[[k]] <- UpdateCovarianceMatrix(new_a_list[[k]],x_step,ps_threshold,cc,ccov)
  }

  # set new parent
  combined_a_list <- append(a_list,new_a_list)
  a_list <- combined_a_list[newPopulationIndex]

  population <- newPopulation
  populationObjective <- newPopulationObjective

  return(list(new_generation=a_list,
              population=population,
              populationObjective=populationObjective))
}
