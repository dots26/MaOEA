#' Do an iteration of population bases steady state MO-CMA-ES. The variation is using SBX and CMA mutation. The original MO-CMA-ES does not use crossover, to do this simply set crossoverProbability to zero.
#' @title Steady-state Multi-Objective CMA-ES
#' @param parent The parent generation, an object of class cmaes_gen. The MO-CMA-ES parent is a 5 tuple: x (the design point, length = number of variable),averageSuccessRate (scalar),stepSize (scalar), evoPath (evolution path, vector, length = number of variable ),covarianceMatrix (square matrix with ncol = nrow = number of variable). The parent is then should be a vector of lists (see example).
#' @param nObjective The number of objective functions. A scalar value.
#' @param problemType A string containing which problem are being solved. Currently available DTLZ1-DTLZ4, WFG4-WFG9. Default to the "easy" DTLZ2.
#' @param successProbTarget Target success probability
#' @param successProbThreshold The threshold for success probability. If the average success probability is higher than this value, the success rate growth is slowed.
#' @param crossoverProbability The probability of doing crossover. Should be between 0-1. Negative value will behave like a zero, and values larger than 1 will behave like 1. Default to 1.
#' @param WFGScaling The use of scaling factor in WFG. Will be ignored in DTLZ problems. Without the scaling, the Pareto front would be on the all-positive portion of hypersphere with radius 1.
#' @param crossoverDistribution The distribution index for SBX. Larger index makes the distribution sharper around each parent.
#' @return Returns a list for the next generation. It contains It contains list$new_generation (class: cmaes_gen), list$population (basically a copy of list$new_generation[[]]$x), and list$populationObjective
#' @examples
#' nVar <- 14
#' nObjective <- 5
#' nIndividual <- 100
#' crossoverProbability <- 1
#' ps_target <- 1 / (5 + ( 1 / 2  ) )
#' pop <- matrix(runif(nIndividual*nVar), nrow = nVar) # create the population
#' a_list <- cmaes_gen(pop)
#'
#' newGeneration <- SCMAES(a_list,nObjective ,"WFG8",crossoverProbability,TRUE) # will run a generation of steady state MOCMAES with standard WFG8 test function.
#' @export
SCMAES <- function(parent,...){
  print(class(parent))
 # UseMethod(SCMAES,parent, nObjective,problem="DTLZ2",successProbTarget,successProbThreshold = 0.44,crossoverProbability=1,WFGScaling=TRUE,crossoverDistribution=30)
  UseMethod('SCMAES')
}

#' @export
SCMAES.cmaes_gen <- function(parent,nObjective ,problem="DTLZ2",successProbTarget,successProbThreshold = 0.44,crossoverProbability=1,WFGScaling=TRUE,crossoverDistribution=30,...){
  ps_threshold <- successProbThreshold
  ps_target <- successProbTarget
  if(!is.na(stringr::str_match(problem,"DTLZ"))){
    if(!is.na(stringr::str_match(problem,"DTLZ1")))
      problemType <- 1
    else
      problemType <- 2
  }
  if(!is.na(stringr::str_match(problem,"WFG"))){
    problemType <- 3
  }
  chromosomeLength <- length(parent[[1]]$x)
  populationSize <- length(parent)

  d <- 1 + chromosomeLength / 2
  cc <- ( 2 / (chromosomeLength + 2))
  ccov <-  2 / ((chromosomeLength^2) + 6)
  cp <- ps_target * 1 / (2 + ps_target*1)

  a_list <- parent

  population<-matrix(,nrow=chromosomeLength,ncol=0);
  for(populationIndex in 1:length(parent)){
    population <- cbind(population,parent[[populationIndex]]$x)
  }

  populationObjective<-matrix(,nrow=nObjective,ncol=0);
  #evaluation of parents
  for(parentIndex in 1:populationSize){
    #individualObjectiveValue<-EvaluateIndividual(population[,parentIndex])
    individualObjectiveValue<-EvaluateIndividual(population[,parentIndex],nObjective)
    populationObjective<-cbind(populationObjective,individualObjectiveValue)
  }

  k <- sample(1:populationSize,1)

  new_a_list <- a_list[[k]]
  new_population <- population[,k]
  new_objectiveValue <- populationObjective[,k]

  if(runif(1) < crossoverProbability)
  {
    parentCount <- 0
    parent_x <- population[,1:2]
    parent_stepSizeAverage <- matrix(,nrow=1)

    while(parentCount<2){
      parent_index <- sample.int(ncol(population),2,replace = FALSE)
      parent_x[,(parentCount+1)] <- a_list[[parent_index[1]]]$x
      parent_x[,(parentCount+2)] <- a_list[[parent_index[2]]]$x

      # step size
      parent_stepSize1 <- a_list[[parent_index[1]]]$stepSize
      parent_stepSize2 <- a_list[[parent_index[2]]]$stepSize

      new_a_list$stepSize <- (parent_stepSize1+parent_stepSize2)/2

      # covariance matrix
      parent_covmatrix1 <- a_list[[parent_index[1]]]$covarianceMatrix
      parent_covmatrix2 <- a_list[[parent_index[2]]]$covarianceMatrix
      new_a_list$covarianceMatrix <- (parent_covmatrix1+parent_covmatrix2)/2

      # evoPath
      parent_evoPath1 <- a_list[[parent_index[1]]]$evoPath
      parent_evoPath2 <- a_list[[parent_index[2]]]$evoPath

      new_a_list$evoPath<- (parent_evoPath1+parent_evoPath2)/2

      # successRate
      new_a_list$averageSuccessRate <- a_list[[parent_index[1]]]$averageSuccessRate

      parentCount<- parentCount +2
    }

    #recombination
    if((problemType==3) && (WFGScaling==TRUE)){
      offspring <- nsga2R::boundedSBXover(t(parent_x),rep(0,chromosomeLength),seq(2,2*chromosomeLength,2),crossoverProbability,30)
    }else{
      offspring <- nsga2R::boundedSBXover(t(parent_x),rep(0,chromosomeLength),rep(1,chromosomeLength),crossoverProbability,30)
    }
    offspring <- t(offspring)

    # update the population's x
    new_a_list$x <- offspring[,1]
  }


  # variation based on multivariate normal distribution
  individual <- new_a_list$x
  covariance <- (new_a_list$stepSize^2) *  new_a_list$covarianceMatrix
  new_a_list$x <- MASS::mvrnorm(1, mu = individual, Sigma = covariance)
  #check feasibility

  newObjectiveValue <- EvaluateIndividual(new_a_list$x,nObjective)
  new_population <- new_a_list$x
  notFeasible <- FALSE

  # non-dominated sorting
  combinedPopulation <- cbind(population,new_population)
  combinedObjectiveValue <- cbind(populationObjective,newObjectiveValue)
  rnk <- nsga2R::fastNonDominatedSorting(t(combinedObjectiveValue))
  rnkIndex <- integer(ncol(combinedObjectiveValue))
  sortingIndex <- 1;
  while (sortingIndex <= length(rnk)) {
    rnkIndex[rnk[[sortingIndex]]] <- sortingIndex;
    sortingIndex <- sortingIndex + 1;
  }

  newPopulation <- matrix(, nrow = chromosomeLength,ncol=0)
  newPopulationIndex <- matrix(, nrow = 1,ncol=0)
  newPopulationObjective <- matrix(, nrow = nObjective,ncol=0)
  newPopulationSize <- integer(1)

  frontIndex <- 1
  while (newPopulationSize< populationSize) #new population is filled by every next front
  {
    newPopulationIndex <- append(newPopulationIndex,rnk[[frontIndex]])
    newPopulation <- cbind( newPopulation,combinedPopulation[ , rnk[[frontIndex]] ] )
    newPopulationObjective <- cbind(newPopulationObjective,combinedObjectiveValue[ , rnk[[frontIndex]] ])
    newPopulationSize <- newPopulationSize + length( rnk[[frontIndex]] )

    frontIndex <- frontIndex+1
  }
  while (newPopulationSize > populationSize) #new population is overcapacity, remove least contributor
  {
    leastContrib <- GetLeastContributor(newPopulationObjective,,getOption("hypervolumeMethod"))
    newPopulationIndex <- newPopulationIndex[-leastContrib]
    newPopulation <- newPopulation[,-leastContrib]
    newPopulationObjective <- newPopulationObjective[,-leastContrib]
    newPopulationSize <- newPopulationSize -1
  }



  # check if the offspring is a success
  successfulUpdate <- max( as.integer ( newPopulationIndex==(1+populationSize) ) )

  # update step size 1
  a_list[[k]] <- UpdateStepSize(a_list[[k]],successfulUpdate,cp,d,ps_target)
  # update step size 2
  new_a_list <- UpdateStepSize(new_a_list,successfulUpdate,cp,d,ps_target)

  # update cov matrix
  x_step <- (new_a_list$x-a_list[[k]]$x) / a_list[[k]]$stepSize
  new_a_list <- UpdateCovarianceMatrix(new_a_list,x_step,ps_threshold,cc,ccov)

  # set new parent
  combined_a_list <- vector('list',populationSize+1)
  for(i in 1:populationSize)
    combined_a_list[i] <- a_list[i]

  combined_a_list[[populationSize+1]] <- new_a_list
  a_list <- combined_a_list[newPopulationIndex]

  population <- newPopulation
  populationObjective <- newPopulationObjective

  for(i in 1:populationSize){
    if(length(a_list[[i]]) != 5) print(c('error',evalIndex,length(a_list[[i]]),i,leastContrib))
  }

  return(list(new_generation=a_list,population=population,populationObjective=populationObjective))
}
