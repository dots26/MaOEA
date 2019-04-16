#' Do an iteration of SMS-EMOA on DTLZ and WFG test functions. THe variation used is SBX and polynomial mutation.
#' @title S-Metric Selection EMOA on DTLZ and WFG test functions
#'
#' @param nObjective The number of objective functions
#' @param population The parent generation. One individual per column.
#' @param problemType A string containing which problem are being solved. Currently available DTLZ1-DTLZ4, WFG4-WFG9. Default to the "easy" DTLZ2.
#' @param crossoverProbability The probability of doing crossover. Should be between 0-1. Negative value will behave like a zero, and values larger than 1 will behave like 1. Default to 1.
#' @param mutationProbability The probability of doing mutation. Should be between 0-1. Negative value will behave like a zero, and values larger than 1 will behave like 1. Default to 1
#' @param WFGScaling The use of scaling factor in WFG. Will be ignored in DTLZ problems. Without the scaling, the Pareto front would be on the all-positive portion of hypersphere with radius 1.
#' @param mutationDistribution The distribution index for polynomial mutation. Larger index makes the distribution sharper around the parent.
#' @param crossoverDistribution The distribution index for SBX. Larger index makes the distribution sharper around each parent.
#' @return Returns a list for the next generation. It contains list$population (the new generation), list$successfulOffspring (normally used for adative schemes), and list$populationObjective
#' @examples
#' nVar <- 14
#' nObjective <- 5
#' nIndividual <- 100
#' crossoverProbability <- 1
#' mutationProbability <- 1/nVar
#' population <- matrix(runif(nIndividual*nVar), nrow = nVar)
#' SMSEMOA(population,nObjective,"WFG8",crossoverProbability,mutationProbability,TRUE) # will run a generation of SMS-EMOA with standard WFG8 test function.
#' @export

SMSEMOA.test <- function( population,nObjective,problem="DTLZ2",crossoverProbability=1,mutationProbability=1,WFGScaling=TRUE,mutationDistribution=20,crossoverDistribution=30,hypervolumeMethod = 'exact'){
  if(!is.na(stringr::str_match(problem,"DTLZ"))){
    if(!is.na(stringr::str_match(problem,"DTLZ1")))
      problemType <- 1
    else
      problemType <- 2
  }
  if(!is.na(stringr::str_match(problem,"WFG"))){
    problemType <- 3
  }
  chromosomeLength <- nrow(population)
  populationSize <- ncol(population)

  newPointSurvives <- TRUE
  #evaluation of parents
  populationObjective<-matrix(,nrow=nObjective,ncol=0);
  for(parentIndex in 1:populationSize){
    #individualObjectiveValue<-EvaluateIndividual(population[,parentIndex])
    individualObjectiveValue<-EvaluateIndividual(population[,parentIndex],nObjective,problem)
    populationObjective<-cbind(populationObjective,individualObjectiveValue)
  }

  # create offspring
  offspring<-matrix(,nrow = chromosomeLength,ncol = 0)
  offspringObjectiveValue<-matrix(,nrow = nObjective,ncol=0)

  parentIndex <- sample(1:populationSize,2,replace = FALSE)

  #Crossover
  if((problemType==3) && (WFGScaling==TRUE)){
    offspring <- nsga2R::boundedSBXover(t(population[,parentIndex]),rep(0,chromosomeLength),seq(2,2*chromosomeLength,2),crossoverProbability,crossoverDistribution)
  }else{
    offspring <- nsga2R::boundedSBXover(t(population[,parentIndex]),rep(0,chromosomeLength),rep(1,chromosomeLength),crossoverProbability,crossoverDistribution)
  }
  offspring <- matrix(offspring[sample(1:2,1),],nrow=1, ncol=chromosomeLength)

  #Mutation
  if((problemType==3) && (WFGScaling==TRUE)){
    offspring <- myBoundedPolyMutation(offspring,rep(0,chromosomeLength),seq(2,2*chromosomeLength,2),mutationProbability,mutationDistribution)
  }else{
    offspring <- myBoundedPolyMutation(offspring,rep(0,chromosomeLength),rep(1,chromosomeLength),mutationProbability,mutationDistribution)
  }
  offspring <- t(offspring)

  # evaluate objective
  offspringObjectiveValue <- EvaluateIndividual(offspring[,1],nObjective,problem)

  combinedPopulation <- cbind(population,offspring[,1])
  combinedObjectiveValue <- cbind(populationObjective,offspringObjectiveValue)

  # nondominated sorting
  rnk <- nsga2R::fastNonDominatedSorting(t(combinedObjectiveValue))

  rnkIndex <- integer(ncol(combinedPopulation))
  sortingIndex <- 1;
  while (sortingIndex <= length(rnk)) {
    rnkIndex[rnk[[sortingIndex]]] <- sortingIndex;
    sortingIndex <- sortingIndex + 1;
  }

  # fill new population
  newPopulation <- matrix(, nrow = chromosomeLength,ncol=0)
  newPopulationObjective <- matrix(, nrow = nObjective,ncol=0)
  newPopulationSize <- integer(1)

  frontIndex <- 1


  while (newPopulationSize < populationSize) #new population is filled by every next front
  {
    newPopulation <- cbind(newPopulation,combinedPopulation[,rnk[[frontIndex]]])
    newPopulationObjective <- cbind(newPopulationObjective,combinedObjectiveValue[,rnk[[frontIndex]]])
    newPopulationSize <- newPopulationSize + length(rnk[[frontIndex]])

    frontIndex <- frontIndex+1
  }

  worstFrontIndex <- rnk[[frontIndex-1]]
  individualIndexTobeChecked <- NULL

  # if there are too many individual in new population, do secondary selection (to remove 1 individual), otherwise return.
  # to do this, we check the last front
  # remove one according to the secondary selection method.
  if(newPopulationSize>populationSize){
    for(individualIndex in worstFrontIndex){
      individualIndexTobeChecked <- c(individualIndexTobeChecked,individualIndex)
    }
    #print("Getting least contributor")
    smallestContributor <- GetLeastContributor(combinedObjectiveValue[,worstFrontIndex], ,hypervolumeMethod)

    #  hypervolumeContribution <- emoa::hypervolume_contribution(combinedObjectiveValue[,worstFrontIndex], rep(20,nObjective))
    #  smallestContributor <- which.is.max(-hypervolumeContribution)

    removedPoint <- individualIndexTobeChecked[smallestContributor]

    newPopulation <- matrix(,nrow = chromosomeLength,ncol = 0)
    newPopulationObjective <- matrix(,nrow = nObjective,ncol = 0)
    for(populationIndex in 1:ncol(combinedPopulation)){
      if(populationIndex != removedPoint){
        newPopulation <- cbind(newPopulation,combinedPopulation[,populationIndex])
        newPopulationObjective <-cbind(newPopulationObjective,combinedObjectiveValue[,populationIndex])
      }
    }
    if(removedPoint == ncol(combinedPopulation))
      newPointSurvives <- FALSE
  }else{
    if(rnk[[length(rnk)]][1] == ncol(combinedPopulation))
      newPointSurvives <- FALSE
  }

  population <- newPopulation
  populationObjective <- newPopulationObjective

  generation <- list(population=population,
                     populationObjective=populationObjective,
                     successfulOffspring = newPointSurvives)
  return(generation)
  #print(c('Rindex',Rindex))
}
