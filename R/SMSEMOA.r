#' Do an iteration of  S-Metric Selection (SMS)-EMOA. The variation used is simulated binary crossover (SBX) and polynomial mutation.
#' @title S-Metric Selection EMOA
#'
#' @param population The parent generation. One individual per column.
#' @param nObjective Number of objective. Ignored as of version 0.6.1; number of row from fun is used instead.
#' @param fun Objective function being solved. Currently available in the package DTLZ1-DTLZ4, WFG4-WFG9.
#' @param control (list) Options to control the SMS-EMOA:
#' \code{mutationProbability} The probability of doing mutation. Should be between 0-1. Negative value will behave like a zero, and values larger than 1 will behave like 1. Default to 1
#' \code{mutationDistribution} The distribution index for polynomial mutation. Larger index makes the distribution sharper around the parent.
#' \code{crossoverDistribution} The distribution index for SBX. Larger index makes the distribution sharper around each parent.
#' \code{referencePoint} The reference point for HV computation on normalized objective space, i.e. (1,...,1) is the nadir point. If not supplied, the ref_multiplier is used instead.
#' \code{ref_multiplier} In case that a reference point is not supplied, the reference is set as a multiply of the current nadir. Default to 1.1.
#' \code{lbound} A vector containing the lower bound for each gene
#' \code{ubound} A vector containing the upper bound for each gene
#' \code{scaleinput} Whether the input should be scaled to 0-1.
#' @return Returns a list for the next generation
#' \code{population} The new generation. Column major, each row contain 1 set of objectives.
#' \code{successfulOffspring} Binary, 1 if the offspring is kept in the new generation. Used in some adaptive schemes.
#' \code{populationObjective} The new generation's objective values.
#' @param ... Further arguments to be passed to \code{fun}
#' @references Beume,  N.,  Naujoks,  B.,  Emmerich,  M.:  SMS-EMOA:  Multiobjective  selection
#' based on dominated hypervolume. Eur. J. Oper. Res. 181 (3), 1653 – 1669 (2007)
#'
#' @examples
#'  \donttest{
#' nVar <- 14
#' nObjective <- 5
#' nIndividual <- 100
#' crossoverProbability <- 1
#' mutationProbability <- 1/nVar
#' population <- matrix(runif(nIndividual*nVar), nrow = nVar)
#'
#' # run a generation of SMS-EMOA with standard WFG6 test function.
#' numpyready <- reticulate::py_module_available('numpy')
#' pygmoready <- reticulate::py_module_available('pygmo')
#' py_module_ready <- numpyready && pygmoready
#' if(py_module_ready) # prevent error on testing the example
#' SMSEMOA(population,WFG6,nObjective,list(crossoverProbability = crossoverProbability,
#'                                           mutationProbability = mutationProbability),nObjective)
#' }
#' @export
SMSEMOA <- function(population,fun,nObjective,control=list(),...){
  if (!pkg.globals$have_pygmo)
    pkg.globals$have_pygmo <- reticulate::py_module_available("pygmo")
  if (!pkg.globals$have_pygmo)
    stop("SMS-EMOA requires PyGMO to compute hypervolume")

  chromosomeLength <- nrow(population)
  populationSize <- ncol(population)

  if(identical(fun,DTLZ4))
    message('DTLZ4 may suffer from floating-point inaccuracy.')

  con <- list(crossoverProbability=1,
              mutationProbability=1,
              mutationDistribution=20,
              crossoverDistribution=30,
              hypervolumeMethod='exact',
              hypervolumeMethodParam=list(),
              referencePoint = NULL,
              scaleinput=T,
              lower=rep(0,chromosomeLength),
              upper=rep(1,chromosomeLength),
              ref_multiplier=1.1)

  con[names(control)] <- control
  control <- con

  ubound <- control$upper
  lbound <- control$lower

  scale_multip <- 1
  scale_shift <- 0
  if(control$scaleinput){
    scale_shift <- -lbound
    scale_multip <- (ubound-lbound)

    population <- ((population) - scale_shift) / scale_multip # scale available population

    ubound <- rep(1,chromosomeLength)
    lbound <- rep(0,chromosomeLength)

    # if(!is.null(population))
      # population <- t((t(population) - scale_shift) / scale_multip) # scale available population
  }

  newPointSurvives <- TRUE
  #evaluation of parents
  populationObjective<-NULL
  for(parentIndex in 1:populationSize){
    ind <- population[,parentIndex,drop=F]*scale_multip + scale_shift
    class(ind) <- class(population)
    individualObjectiveValue<-fun(ind,...)
    populationObjective<-cbind(populationObjective,individualObjectiveValue)
  }

  # create offspring
  offspring<-matrix(,nrow = chromosomeLength,ncol = 0)
  offspringObjectiveValue<-NULL

  parentIndex <- sample(1:populationSize,2,replace = FALSE)
  #Crossover
  offspring <- nsga2R::boundedSBXover(parent_chromosome = t(population[,parentIndex]),
                                      lowerBounds = lbound,
                                      upperBounds = ubound,
                                      cprob = con$crossoverProbability,
                                      mu = con$crossoverDistribution)
  offspring <- matrix(offspring[sample(1:2,1),],nrow=1, ncol=chromosomeLength)
  #Mutation
  offspring <- nsga2R::boundedPolyMutation(parent_chromosome = offspring,
                                           lowerBounds = lbound,
                                           upperBounds = ubound,
                                           mprob = con$mutationProbability,
                                           mum = con$mutationDistribution)
  offspring <- t(offspring)

   # evaluate objective
  this_offspring <- offspring[,1,drop=FALSE]


  class(this_offspring) <- class(population)
  offspringObjectiveValue <- fun(this_offspring,...)

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
  newPopulationObjective <- NULL
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
    # normalize objectives to 0-1, minimization problem, ideal point is at (0,0)
    normalizationList <- Normalize(newPopulationObjective)
    normalizedObjective <- normalizationList$normalizedObjective
    idealPoint <- normalizationList$idealPoint

    for(individualIndex in worstFrontIndex){
      individualIndexTobeChecked <- c(individualIndexTobeChecked,individualIndex)
    }
    smallestContributor <- GetLeastContributor(combinedObjectiveValue[,worstFrontIndex],
                                               reference=control$referencePoint,
                                               method=control$hypervolumeMethod,
                                               hypervolumeMethodParam=control$hypervolumeMethodParam,
                                               ref_multiplier=control$ref_multiplier)
    removedPoint <- individualIndexTobeChecked[smallestContributor]

    newPopulation <- newPopulation[,-removedPoint]
    newPopulationObjective <- newPopulationObjective[,-removedPoint]
    if(removedPoint == ncol(combinedPopulation))
      newPointSurvives <- FALSE
  }else{
    if(rnk[[length(rnk)]][1] == ncol(combinedPopulation))
      newPointSurvives <- FALSE
  }

  keep_class <- class(population)
  population <- newPopulation*scale_multip + scale_shift
  class(population) <- keep_class
  populationObjective <- newPopulationObjective

  generation <- list(population=(population),
                     objective=(populationObjective),
                     successfulOffspring = newPointSurvives)

  return(generation)
}

