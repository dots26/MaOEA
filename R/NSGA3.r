#' Do an iteration of Elitist Non-dominated Sorting Genetic Algorithm version III (NSGA-III). THe variation is using SBX and polynomial mutation.
#' @title Elitist Non-dominated Sorting Genetic Algorithm version III
#' @param population The parent generation. One individual per column. nrow = number of variable, ncol = number of individuals in the population.
#' @param fun Objective function being solved. Currently available in the package DTLZ1-DTLZ4, WFG4-WFG9.
#' @param nObjective The number of objective functions. A scalar value.
#' @param control A list, containing the following:
#' \code{weightVector} NSGA-III require a set of reference points defined a priori. The reference can be any point. If not supplied, 5*nObjective points are generated from a sobol sequence. Column major: nrow = nObjective, ncol = number of reference points
#' \code{crossoverProbability} The probability of doing crossover. Should be between 0-1. Negative value will behave like a zero, and values larger than 1 will behave like 1. Default to 1.
#' \code{mutationProbability} The probability of doing mutation. Should be between 0-1. Negative value will behave like a zero, and values larger than 1 will behave like 1. Default to 1
#' \code{mutationDistribution} The distribution index for polynomial mutation. Larger index makes the distribution sharper around the parent.
#' \code{crossoverDistribution} The distribution index for SBX. Larger index makes the distribution sharper around each parent.
#' @param ... Further arguments to be passed to \code{fun}
#' @return #' @return Returns a list for the next generation
#' \code{population} The new generation design points. Column major.
#' \code{populationObjective} The new generation's objective values. Column major.
#' @references Deb, K., Jain, H.: An evolutionary many-objective optimization algorithm using
#' reference-point-based  nondominated  sorting  approach,  part  I:  Solving  problems with box constraints.
#' Trans. Evol. Comp. 18 (4), 577–601 (2014)
#'
#' @examples
#' nVar <- 14
#' nObjective <- 5
#' nIndividual <- 100
#' #control for NSGA3
#' ctrl <- list(crossoverProbability = 1,
#'              mutationProbability = 1/nVar)
#' #Initial population
#' population <- matrix(runif(nIndividual*nVar), nrow = nVar)
#'
#' # run a generation of NSGA-III with standard WFG8 test function.
#' NSGA3(population, WFG8,nObjective,ctrl,nObjective)
#' @export
NSGA3 <- function(population,fun,nObjective,control=list(),...){
  con <- list(crossoverProbability=1,
              mutationProbability=1,
              mutationDistribution=20,
              crossoverDistribution=30,
              weightVector=NULL)
  con[names(control)] <- control

  if(identical(fun,DTLZ4))
    message('DTLZ4 may suffer from floating-point inaccuracy.')

  if(is.null(con$weightVector)){
    con$weightVector <- createWeightsSobol(nDim=nObjective,nWeights = 5*nObjective)
  }
  nRefPoint <- ncol(con$weightVector)

  if(nObjective!= nrow(con$weightVector)){
    stop('Number of rows in reference set must be the same as number of objective!' )
  }
  chromosomeLength <- nrow(population)
  populationSize <- ncol(population)

  keepClass <- class(population)

  populationObjective<-matrix(,nrow=nObjective,ncol=0);
  for(parentIndex in 1:populationSize){
    #individualObjectiveValue<-EvaluateIndividual(population[,parentIndex])
    ind <- population[,parentIndex,drop=F]
    class(ind) <- class(population)
    individualObjectiveValue<-EvaluateIndividual(ind,fun,...)
    populationObjective<-cbind(populationObjective,individualObjectiveValue)
  }

  offspring<-matrix(,nrow = chromosomeLength,ncol = 0)
  offspringObjectiveValue<-matrix(,nrow = nObjective,ncol=0);
  parentCount <- 0
  parent <- population

  while(parentCount<ncol(population)){
    parent[,(parentCount+1):(parentCount+2)] <- population[,sample.int(ncol(population),2,replace = FALSE)]
    parentCount<- parentCount +2
  }

  #recombination
  offspring <- nsga2R::boundedSBXover(t(population),rep(0,chromosomeLength),rep(1,chromosomeLength),con$crossoverProbability,con$crossoverDistribution)

  #Mutation
  offspring <- nsga2R::boundedPolyMutation(offspring,rep(0,chromosomeLength),rep(1,chromosomeLength),con$mutationProbability,con$mutationDistribution)
  offspring <- t(offspring)


  #evaluation of offsprings
  for(offspringIndex in 1:populationSize){
    #  offspringObjectiveValue <- cbind(offspringObjectiveValue,EvaluateIndividual(offspring[,offspringIndex]))
    offspr <- offspring[,offspringIndex,drop=F]
    class(offspr) <- keepClass
    offspringObjectiveValue <- cbind(offspringObjectiveValue,EvaluateIndividual(offspr,fun,...))
  }

  # start selection procedure of NSGA-III
  combinedPopulation <- cbind(population,offspring)
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

  # if there are too many individual in new population, do secondary selection, otherwise return.
  # to do this, the new population is reverted by one front.
  # Then, one by one the members of the last front are added according to the secondary selection method.
  if(newPopulationSize>populationSize){
    # NSGA-III
    lastAccomodatedFront <- frontIndex - 2
    if(lastAccomodatedFront>0){
      lastAccomodatedSize <- newPopulationSize - length(rnk[[frontIndex-1]])

      newPopulationWithoutLastFront <- newPopulation[,1:lastAccomodatedSize]
      newPopulationObjectiveWithoutLastFront <- newPopulationObjective[,1:lastAccomodatedSize]
    }
    else{
      newPopulationWithoutLastFront <- newPopulation[0]
      newPopulationObjectiveWithoutLastFront <- newPopulationObjective[0]
      lastAccomodatedSize <- 0
    }
    populationSizeToBeFilled <- populationSize - lastAccomodatedSize

    # normalize objectives to 0-1, minimization problem, so ideal point will be at (0,0)
    normalizationList <- AdaptiveNormalization(newPopulationObjective)
    normalizedObjective <- normalizationList$normalizedObjective
    idealPoint <- normalizationList$idealPoint

    # calculate distances from reference lines
    distFromLine <- CalcDistRefLineToSolution(normalizedObjective,con$weightVector)

    # associate solution with refline
    associatedLine <- AssociateLine(distFromLine)

    solutionsAddedFromWorstFront <- matrix(,nrow = chromosomeLength,ncol = 0)
    objectiveAddedFromWorstFront<- matrix(,nrow = nObjective,ncol = 0)

    lineNiching <- matrix(rep(0,nRefPoint),nrow = nRefPoint,ncol=1)

    # count associated solution for each refline (not including the worst front)
    if(lastAccomodatedSize>0){
      for (reflineIndex in 1:nRefPoint) {
        lineNiching[reflineIndex] <- sum(as.integer(associatedLine[1:lastAccomodatedSize] == reflineIndex)) # do not include the worst front
      }
    }
    addedSolutionCount <- 0

    while (addedSolutionCount < populationSizeToBeFilled) {
      # get which refline has the smallest number of associated solution
      sadRefline <- nnet::which.is.max(-lineNiching) #sad because it has the least member

      # get which solution can be associated with the sadRefline
      probableChoice <- as.integer(associatedLine[(lastAccomodatedSize+1):newPopulationSize]==sadRefline)
      probableChoiceIndex <- matrix(,nrow = 1,ncol = 0)

      for (index in (1:length(worstFrontIndex))) {
        if(probableChoice[index]==1){
          probableChoiceIndex  <- cbind(probableChoiceIndex,index)
        }
      }

      # if sadRefline has no chance of getting any more associated solution, do not consider it again (make it infinite)
      if(sum(probableChoice)==0){
        lineNiching[sadRefline] <- Inf
      }
      # if sadRefline gets any associated solution, add to the niche and increment addedSolutionCount
      else{
        # get the minimum distance to the sadRefline
        addedSolution <- nnet::which.is.max (-distFromLine[sadRefline,probableChoiceIndex])
        #add solution with smallest distance to the sadRefline
        solutionsAddedFromWorstFront <- cbind(solutionsAddedFromWorstFront,newPopulation[,probableChoiceIndex[addedSolution]])
        objectiveAddedFromWorstFront <- cbind(objectiveAddedFromWorstFront,newPopulationObjective[,probableChoiceIndex[addedSolution]])

        distFromLine[,probableChoiceIndex[addedSolution]] <- Inf

        addedSolutionCount <- addedSolutionCount + 1
        lineNiching[sadRefline] <- lineNiching[sadRefline]+1
      }
    } # endwhile: adding enough solution too fill the new population
    newPopulation <- cbind(newPopulationWithoutLastFront, solutionsAddedFromWorstFront)
    newPopulationObjective <- cbind(newPopulationObjectiveWithoutLastFront, objectiveAddedFromWorstFront)
  } #endif newpopulation too large
  population <- newPopulation
  populationObjective <- newPopulationObjective

  generation <- list(population=(population), objective=(populationObjective))
  return(generation)
}
