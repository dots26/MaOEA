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
#' \code{crossoverMethod} Either "sbx" (simulated binary) or "uniform"
#' \code{mutationMethod} Either "polymut" (polynomial mutation), "truncgauss" (truncated gaussian), "ortho" (orthogonal sampling)
#' \code{searchDir} The search direction used for orthogonal sampling. Will be orthonormalized using the Gram-Schmidt process.
#' \code{nDirection} Used in "ortho". Number of orthogonal search direction used.
#' \code{checkMirror} If TRUE, and "ortho" is used, then the mirrorred mutation is checked for one additional evaluation for each offspring.
#' \code{scaleinput} Whether the input should be scaled to 0-1.
#' \code{stepsize} Stepsize for gaussian mutation. Default to 0.3.
#' @return Returns a list for the next generation
#' \code{population} The new generation. Column major, each row contain 1 set of objectives.
#' \code{successfulOffspring} Binary, 1 if the offspring is kept in the new generation. Used in some adaptive schemes.
#' \code{populationObjective} The new generation's objective values.
#' \code{searchDir} The search direction used. Not NULL if orthogonal sampling is used. Can be used to do weighted optimization.
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
#' SMSEMOA(population = population,
#'         fun = WFG6,
#'         nObjective = nObjective,
#'         populationObjective = NULL,
#'         control = list(crossoverProbability = crossoverProbability,
#'                        mutationProbability = mutationProbability,
#'                        mutationMethod = "poly"),
#'         nObj=nObjective)
#' }
#' @export
SMSEMOA <- function(population,fun,
                    nObjective=nrow(populationObjective),
                    populationObjective=NULL,
                    control=list(),...){
  nEval<-0
  searchDir <- NULL
  if (!pkg.globals$have_pygmo)
    pkg.globals$have_pygmo <- reticulate::py_module_available("pygmo")
  if (!pkg.globals$have_pygmo)
    stop("SMS-EMOA requires PyGMO to compute hypervolume")


  chromosomeLength <- nrow(population)
  populationSize <- ncol(population)

  if(identical(fun,DTLZ4))
    message('DTLZ4 may suffer from floating-point inaccuracy.')

  con <- list(crossoverProbability=0.5,
              mutationProbability=1/chromosomeLength,
              mutationDistribution=20,
              crossoverDistribution=30,
              hypervolumeMethod='exact',
              hypervolumeMethodParam=list(),
              referencePoint = NULL,
              stepsize=0.3,
              scaleinput=T,
              lower=rep(0,chromosomeLength),
              upper=rep(1,chromosomeLength),
              ref_multiplier=1.1,
              crossoverMethod =c("uniform","sbx"),
              mutationMethod = c("truncgauss","poly","ortho"),
              nDirection = 1,
              checkMirror = T)
  con[names(control)] <- control
  control <- con
  stepsize <- con$stepsize
  checkMirror <- con$checkMirror
  ubound <- control$upper
  lbound <- control$lower
  crossovermethod <-control$crossoverMethod[1]
  mutationmethod <- control$mutationMethod[1]

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

  if(is.null(populationObjective))
  {
    for(parentIndex in 1:populationSize){
      ind <- population[,parentIndex,drop=F]*scale_multip + scale_shift
      class(ind) <- class(population)
      individualObjectiveValue<-fun(ind,...)
      populationObjective<-cbind(populationObjective,individualObjectiveValue)
      nEval <- nEval +1
    }
  }
  # print((populationObjective))
  # create offspring
  offspring<-matrix(,nrow = chromosomeLength,ncol = 0)
  offspringObjectiveValue<-NULL

  parentIndex <- sample(1:populationSize,2,replace = FALSE)

  #Crossover

  if(crossovermethod=="sbx"){
    offspring <- nsga2R::boundedSBXover(parent_chromosome = t(population[,parentIndex]),
                                        lowerBounds = lbound,
                                        upperBounds = ubound,
                                        cprob = con$crossoverProbability,
                                        mu = con$crossoverDistribution)
  }else{
    offspring <- uniformXover(parent_chromosome = t(population[,parentIndex]),
                              lowerBounds = lbound,
                              upperBounds =  ubound,
                              cprob = con$crossoverProbability)
  }
  offspring <- matrix(offspring[sample(1:2,1),],nrow=1, ncol=chromosomeLength)
  #Mutation
  eff_stepsize <- stepsize#[parentIndex[mainParent]]
  eff_pm <- con$mutationProbability
  if(mutationmethod=="poly"){
    offspring <- nsga2R::boundedPolyMutation(parent_chromosome = offspring,
                                             lowerBounds = lbound,
                                             upperBounds = ubound,
                                             mprob = con$mutationProbability,
                                             mum = con$mutationDistribution)
  }else if(mutationmethod=="ortho"){
    ortho <- orthogonal_sampling_mutation(parent_chromosome = offspring,
                                          lowerBounds = lbound,
                                          upperBounds =  ubound,
                                          mprob = eff_pm,
                                          sigma = eff_stepsize,
                                          nOffspring = 1,
                                          nDirection = nDirection,
                                          searchDir = searchDir)
    searchDir <- ortho$searchDir
    offspring <- ortho$offspring

    if(checkMirror)
      offspring <- cbind(offspring,ortho$mirror)
  }else{
    offspring <- truncnormMutation(parent_chromosome = offspring,
                                   lowerBounds = lbound,
                                   upperBounds =  ubound,
                                   mprob = eff_pm,
                                   sigma = eff_stepsize)
  }
  offspring <- t(offspring)

  # evaluate objective
  offspringHV <- rep(0,ncol(offspring))
  offspringObjectiveValueCol <- NULL
  # print(paste("noffspring", ncol(offspring)))
  for(offspringIndex in 1:ncol(offspring)){
    this_offspring <- offspring[,offspringIndex,drop=FALSE]
    class(this_offspring) <- class(population)
    offspringObjectiveValue <- fun(this_offspring*scale_multip + scale_shift,...)
    offspringObjectiveValueCol <- cbind(offspringObjectiveValueCol,offspringObjectiveValue)
    nEval <- nEval + 1
    combinedPopulation <- cbind(population,offspring[,offspringIndex])
    combinedObjectiveValue <- cbind(populationObjective,offspringObjectiveValue)

    # nondominated sorting
    rnk <- nsga2R::fastNonDominatedSorting(t(combinedObjectiveValue))

    rnkIndex <- integer(ncol(combinedPopulation))
    sortingIndex <- 1;
    while (sortingIndex <= length(rnk)) {
      rnkIndex[rnk[[sortingIndex]]] <- sortingIndex;
      sortingIndex <- sortingIndex + 1;
    }

    # get offspring HV contrib, one at a time
    worstFront <- length(rnk)
    worstFrontIndex <- rnk[[worstFront]]
    if(length(worstFrontIndex)>1){
      smallestContributor <- GetLeastContributor(combinedObjectiveValue[,worstFrontIndex],
                                                 reference=control$referencePoint,
                                                 method=control$hypervolumeMethod,
                                                 hypervolumeMethodParam=control$hypervolumeMethodParam,
                                                 ref_multiplier=control$ref_multiplier)
      removedPoint <- worstFrontIndex[smallestContributor]
      if(removedPoint == ncol(combinedPopulation)){
        newPointSurvives <- FALSE
        offspringHV[offspringIndex] <- 0
      }else{
        for(frontId in 1:length(rnk)){
          if(any(rnk[[frontId]]==ncol(combinedPopulation))){
            offspringFrontIndex <- frontId
            break
          }
        }
        offspringFrontLength <- length(offspringFrontIndex)
        offspringHV[offspringIndex] <- GetHVContribution(combinedObjectiveValue[,offspringFrontIndex],
                                                         reference=control$referencePoint,
                                                         method=control$hypervolumeMethod,
                                                         ref_multiplier=control$ref_multiplier)[offspringFrontLength]
      }
    }else{
      removedPoint <- worstFrontIndex
      if(worstFrontIndex==ncol(combinedPopulation)){
        newPointSurvives <- FALSE
        offspringHV[offspringIndex] <- 0
      }else{
        for(frontId in 1:length(rnk)){
          if(any(rnk[[frontId]]==ncol(combinedPopulation))){
            offspringFrontIndex <- frontId
            break
          }
        }
        offspringFrontLength <- length(offspringFrontIndex)
        offspringHV[offspringIndex] <- GetHVContribution(combinedObjectiveValue[,offspringFrontIndex],
                                                         reference=control$referencePoint,
                                                         method=control$hypervolumeMethod,
                                                         ref_multiplier=control$ref_multiplier)[offspringFrontLength]
      }
    }
  }
  # pick one offspring
  # #TODO get RECORDD rmovd point
  chosenOffspring <- which.max(offspringHV)
  newPopulation <- cbind(population,offspring[,chosenOffspring])[,-removedPoint]
  newPopulationObjective <- cbind(populationObjective,
                                  offspringObjectiveValueCol[,chosenOffspring])[,-removedPoint]
  keep_class <- class(population)
  population <- newPopulation*scale_multip + scale_shift
  class(population) <- keep_class
  populationObjective <- newPopulationObjective

  generation <- list(population=(population),
                     objective=(populationObjective),
                     successfulOffspring = newPointSurvives,
                     searchDirection = searchDir)

  return(generation)
}


