#' Create initial sample using Latin Hypercube Sampling (LHS) method. The variables will be ranged between 0-1
#' @title Initialize population with Latin Hypercube Sampling
#' @param numberOfIndividuals The number of individual in the population. Integer > 0.
#' @param chromosomeLength The number of variables per individual
#' @param minVal Minimum value of the resulting sample
#' @param maxVal Maximum value of the resulting sample
#' @param samplingMethod Not used
#' @return A matrix of size chromosomeLength x nIndividual.
#' @examples
#' nVar <- 14
#' nIndividual <- 100
#' InitializePopulationLHS(nIndividual,nVar,FALSE)
#' @export
#'
InitializePopulationLHS <- function(numberOfIndividuals,chromosomeLength, minVal=0,maxVal=1,samplingMethod=0) {
  #population<-optimumLHS(n=numberOfIndividuals,k=chromosomeLength,maxSweeps=10,eps=.1,verbose=FALSE)
  population<-lhs::randomLHS(n=numberOfIndividuals,k=chromosomeLength)

  population<-t(population)

#  if(binaryEncoding==TRUE){
#    population<-round(population)
#  }else{
    population<-population * (maxVal-minVal) + minVal
#  }

  return(population)
}

#' Evaluate individual with the specified test function. Non-feasible solution are given Inf as objective values.
#' @title Evaluate objective values of a single individual
#' @param individual The individual to be evaluated
#' @param fun A string containing which problem is being solved. Currently available DTLZ1-DTLZ4, WFG4-WFG9.
#' @param ... Further parameters used by \code{fun}
#' @return A matrix of size nObjective, containing the objective values.
#' @examples
#' individual <- stats::runif(8)
#' EvaluateIndividual(individual,WFG4,3) # the 3 is passed to WFG4 nObj
#' @export
EvaluateIndividual <- function(individual,fun,...){
  nVar <- length(individual)

  objective <- fun(individual,...)
  return(objective)
}

#' Evaluate a population with the specified test function. Non-feasible solution are given Inf as objective values.
#' @title Evaluate objective value of a set of individuals
#' @param pop The population to be evaluated
#' @param fun A string containing which problem is being solved. Currently available in the package: DTLZ1-DTLZ4, WFG4-WFG9.
#' @param ... Further parameters used by \code{fun}
#' @return A matrix of size nObjective, containing the objective values.
#' @examples
#' pop <- matrix(runif(8*50),nrow=8) # 8 variables, 50 individuals
#' EvaluatePopulation(pop,WFG4,3) # the 3 is passed to WFG4 nObj
#' @export
EvaluatePopulation <- function(pop,fun,...){
  popSize <- ncol(pop)
  popObjective <- NULL
  for(individualIndex in 1:popSize){
    indObjective <- EvaluateIndividual(pop[,individualIndex],fun,...)
    popObjective <- cbind(popObjective,indObjective)
  }

  return(popObjective)
}

#' The DTLZ1 test function.
#' @param individual The individual to be evaluated
#' @param nObj The number of objective
#' @return A matrix of size nObjective, containing the objective values.
#' @references Deb,  K.,  Thiele,  L.,  Laumanns,  M.,  Zitzler,  E.:  Scalable  Multi-Objective  Optimization Test Problems. In: Congress on Evolutionary Computation (CEC). pp. 825–830. IEEE Press, Piscataway, NJ (2002)
#' @examples
#' individual <- stats::runif(14)
#' nObj <- 4
#' DTLZ1(individual,nObj)
#' @export
DTLZ1 <- function(individual,nObj){
  nVar <- length(individual)
  obj <- matrix(rep(1,nObj),nrow = nObj,ncol = 1)
  gSigma <- 0
  for (subIndex in nObj:nVar) {
    gSigma <- gSigma + ((individual[subIndex] - 0.5)^2 - cos(20*pi*((individual[subIndex] - 0.5))))
  }
  g <- 100*(nVar- nObj + 1 + gSigma)

  for(objectiveIndex in 1:nObj){
    obj[objectiveIndex] <- 0.5* (1+g)

    if( (nObj-objectiveIndex) > 0){
      for(cosIndex in 1:(nObj-objectiveIndex)){
        obj[objectiveIndex] <- obj[objectiveIndex] * individual[cosIndex]
      }
    }
    if(objectiveIndex > 1)
      obj[objectiveIndex] <- obj[objectiveIndex] *  (1 - individual[nObj - objectiveIndex + 1])
  }

  return(obj)
}

#' The DTLZ2 test function.
#' @param individual The individual to be evaluated
#' @param nObj The number of objective
#' @return A matrix of size nObjective, containing the objective values.
#' @references Deb,  K.,  Thiele,  L.,  Laumanns,  M.,  Zitzler,  E.:  Scalable  Multi-Objective  Optimization Test Problems. In: Congress on Evolutionary Computation (CEC). pp. 825–830. IEEE Press, Piscataway, NJ (2002)
#' @examples
#' individual <- stats::runif(14)
#' nObj <- 4
#' DTLZ2(individual,nObj)
#' @export
DTLZ2 <- function(individual,nObj){
  nVar <- length(individual)
  obj <- matrix(rep(1,nObj),nrow = nObj,ncol = 1)
  g <- 0
  for (subIndex in nObj:nVar) {
    g <- g + (individual[subIndex] - 0.5)^2
  }
  for(objectiveIndex in 1:nObj){
    obj[objectiveIndex] <- (1+g)

    if( (nObj-objectiveIndex) > 0){
      for(cosIndex in 1:(nObj-objectiveIndex)){
        obj[objectiveIndex] <- obj[objectiveIndex] * cos(individual[cosIndex] * pi / 2)
      }
    }
    if(objectiveIndex > 1)
      obj[objectiveIndex] <- obj[objectiveIndex] *  sin(individual[nObj - objectiveIndex + 1] * pi / 2)
  }
  return(obj)
}

#' The DTLZ3 test function.
#' @param individual The individual to be evaluated
#' @param nObj The number of objective
#' @return A matrix of size nObjective, containing the objective values.
#'
#' @references Deb,  K.,  Thiele,  L.,  Laumanns,  M.,  Zitzler,  E.:  Scalable  Multi-Objective  Optimization Test Problems. In: Congress on Evolutionary Computation (CEC). pp. 825–830. IEEE Press, Piscataway, NJ (2002)
#' @examples
#' individual <- stats::runif(14)
#' nObj <- 4
#' DTLZ3(individual,nObj)
#' @export
DTLZ3 <- function(individual,nObj){
  nVar <- length(individual)
  obj <- matrix(rep(1,nObj),nrow = nObj,ncol = 1)
  gSigma <- 0
  for (subIndex in nObj:nVar) {
    gSigma <- gSigma + ((individual[subIndex] - 0.5)^2 - cos(20*pi*((individual[subIndex] - 0.5))))
  }
  g <- 100*(nVar- nObj + 1 + gSigma)
  for(objectiveIndex in 1:nObj){
    obj[objectiveIndex] <- (1+g)

    if( (nObj-objectiveIndex) > 0){
      for(cosIndex in 1:(nObj-objectiveIndex)){
        obj[objectiveIndex] <- obj[objectiveIndex] * cos(individual[cosIndex] * pi / 2)
      }
    }
    if(objectiveIndex > 1)
      obj[objectiveIndex] <- obj[objectiveIndex] *  sin(individual[nObj - objectiveIndex + 1] * pi / 2)
  }
  return(obj)
}

#' The DTLZ4 test function.
#' @param individual The individual to be evaluated
#' @param nObj The number of objective
#' @param alpha Alpha value of DTLZ4 function.
#' @return A matrix of size nObjective, containing the objective values.
#' @references Deb,  K.,  Thiele,  L.,  Laumanns,  M.,  Zitzler,  E.:  Scalable  Multi-Objective  Optimization Test Problems. In: Congress on Evolutionary Computation (CEC). pp. 825–830. IEEE Press, Piscataway, NJ (2002)
#' @examples
#' individual <- stats::runif(14)
#' nObj <- 4
#' DTLZ4(individual,nObj)
#' @export
DTLZ4 <- function(individual,nObj,alpha=100){
  nVar <- length(individual)
  obj <- matrix(rep(1,nObj),nrow = nObj,ncol = 1)
  g <- 0
  for (subIndex in nObj:nVar) {
    g <- g + (individual[subIndex] - 0.5)^2
  }
  for(objectiveIndex in 1:nObj){
    obj[objectiveIndex] <- (1+g)

    if( (nObj-objectiveIndex) > 0){
      for(cosIndex in 1:(nObj-objectiveIndex)){
        obj[objectiveIndex] <- (obj[objectiveIndex]) * cos((individual[cosIndex]^alpha) * pi / 2)
      }
    }
    if(objectiveIndex > 1)
      obj[objectiveIndex] <- obj[objectiveIndex] *  sin(individual[nObj - objectiveIndex + 1] * pi / 2)
  }
  return(obj)
}



#' Normalize the objectives to 0-1. The origin is the ideal point. (1,...,1) is not the nadir point.
#' The normalization is done by using adaptive normalization used in NSGA-III.
#'
#' @title Objective space normalization.
#' @param objectiveValue Set of objective vectors to normalize
#'
#' @return A list containing the following:
#' \code{normalizedObjective} The normalized values
#' \code{idealPoint} The ideal point corresponding to the origin
#' \code{nadirPoint} The location of nadir point in the normalized Space
#' @examples
#' nObj <- 5
#' nIndividual <- 100
#' nVar <- 10
#' population <- InitializePopulationLHS(nIndividual,nVar,FALSE)
#' objective <- matrix(,nrow=nObj,ncol=nIndividual)
#' for(individual in 1:nIndividual){
#'    objective[,individual] <- WFG4(population[,individual],nObj)
#' }
#' AdaptiveNormalization(objective)
#' @export
#'
AdaptiveNormalization <- function(objectiveValue){
  minObjVal <- matrix(,nrow = nrow(objectiveValue),ncol = 1)
  nObjective <- nrow(objectiveValue)
  idealPoint <- rep(Inf,nObjective)
  nadirPoint <- rep(-Inf,nObjective)

  indexOfMax <- matrix(,nrow = nObjective,ncol = 1)
  for (i in 1:nObjective) {
    indexOfMax[i] <- nnet::which.is.max(objectiveValue[i,]) # get the individual with max objective for objective i
    minObjVal[i] <- min(objectiveValue[i,])
    if(minObjVal[i]>idealPoint[i])
      minObjVal[i] <- idealPoint[i]
    else
      idealPoint[i] <- minObjVal[i]
  }

  # get the normalization minimum value to zero
  normalizedObjective <- objectiveValue - rep(minObjVal,ncol(objectiveValue))
  #get the intercept, set it as maximum value
  for (i in 1:nObjective) {
    summedMaxRow <- sum(normalizedObjective[,indexOfMax[i]])
    normalizedObjective[i,] <- normalizedObjective[i,] / summedMaxRow
    nadirPoint[i] <- max(normalizedObjective[i,])
    #   normalizedRefPoint[i,] <- normalizedRefPoint[i,] / summedMaxRow
  }


  newList <- list(normalizedObjective = normalizedObjective,idealPoint = idealPoint, nadirPoint = nadirPoint)# "numeric" = normalizedRefPoint)

  return(newList)
}

#' Normalize the objectives AND reference (combined) to 0-1. The origin is the ideal point. (1,...,1) is the nadir.
#'
#' @title Objective space normalization.
#' @param objectiveValue Set of objective vectors to normalize
#' @param referencePoints Set of reference points to transform following the objective vector normalization
#'
#' @return A list containing the following:
#' \code{normalizedObjective} The normalized values
#' \code{idealPoint} The ideal point corresponding to the origin
#' \code{transformedReference} The location of reference points in the normalized Space
#' @examples
#' nObj <- 5
#' nVar <- 10
#' nIndividual <- 100
#' population <- InitializePopulationLHS(nIndividual,nVar,FALSE)
#' objective <- matrix(,nrow=nObj,ncol=nIndividual)
#' for(individual in 1:nIndividual){
#'    objective[,individual] <- WFG4(population[,individual],nObj)
#' }
#' Normalize(objective)
#' @export
#'
Normalize <- function(objectiveValue,referencePoints=NULL){
  minObjVal <- matrix(,nrow = nrow(objectiveValue),ncol = 1)
  nObjective <- nrow(objectiveValue)
  idealPoint <- rep(Inf,nObjective)
  popSize <- ncol(objectiveValue)
  if(!is.null(referencePoints)) {
    referencePoints <- matrix(referencePoints,nrow = nObjective)
    objectiveValue <- cbind(objectiveValue,referencePoints)
  }
  indexOfMax <- matrix(,nrow = nObjective,ncol = 1)

  for (i in 1:nObjective) {
    indexOfMax[i] <- nnet::which.is.max(objectiveValue[i,]) # get the individual with max objective for objective i
    minObjVal[i] <- min(objectiveValue[i,])
    if(minObjVal[i]>idealPoint[i])
      minObjVal[i] <- idealPoint[i]
    else
      idealPoint[i] <- minObjVal[i]
  }
  # get the normalization minimum value to zero
  normalizedObjective <- objectiveValue - rep(minObjVal,ncol(objectiveValue))

  #get the intercept, set it as maximum value
  for (i in 1:nObjective) {
    maxObj <- (normalizedObjective[i,indexOfMax[i]])
    normalizedObjective[i,] <- normalizedObjective[i,] / maxObj
  }

  newList <- list(normalizedObjective = normalizedObjective[,1:popSize], transformedReference =normalizedObjective[,-(1:popSize)],idealPoint = idealPoint)# "numeric" = normalizedRefPoint)

  return(newList)
}

# Used in NSGA-III line-association
CalcDistanceFromLine <- function(pointAInLine, pointBInLine, pointC){
  pa <- pointC - pointAInLine;
  ba <- pointBInLine - pointAInLine;
  t  <- (pa %*% ba)/(ba %*% ba);

  distance = norm(matrix(pa - (ba*c(t)),nrow=length(pointAInLine),ncol=1),type="f")
  return(distance)
}

# Used in NSGA-III line-association
CalcDistRefLineToSolution <- function(normalizedObjective,referencePoint){
  nRefPoint <- ncol(referencePoint)
  populationSize <- ncol(normalizedObjective)
  origin <- rep(0,nrow(normalizedObjective))
  distFromLine <- matrix(,nrow = nRefPoint,ncol = 0)

  for(objectiveIndex in 1:(populationSize)){
    currentPointDistanceFromLine <- matrix(,nrow = nRefPoint,ncol = 1)
    for (refPointIndex in 1:nRefPoint) {
      distanceToLine <- CalcDistanceFromLine(origin, referencePoint[,refPointIndex],normalizedObjective[,objectiveIndex])
      currentPointDistanceFromLine[refPointIndex] <- distanceToLine
    }
    distFromLine <- cbind(distFromLine,currentPointDistanceFromLine)

  } # end for calculate distances of all P+Q individuals from reference lines, we obtain distFromLine

  return(distFromLine) #each column represent different solution, each row different refline
}

CalcIGDDistanceToObjective <- function(objective,referencePoint){
  nRefPoint <- ncol(referencePoint)
  populationSize <- ncol(objective)
  nObjective <- nrow(objective)
  distances <- matrix(,nrow = nRefPoint,ncol = populationSize)

  intersectionPoint <- referencePoint
  # get the intersection of the refline and sigma(x^2)=1
  # x^2  + y^2  + z^2  = 1
  # the refline is the line from (0,0,0) to some {(x,y,z) | x+y+z=1} (the ref point)
  # which means if for example the ref point is at (x*,y*,z*)
  # the line has length = a = root(x*^2+y*^2+z*^2)
  # to make the line to have length = 1, we simply multiply the initial ref point by 1/a
  #

  for (refPointIndex in 1:nRefPoint){
    a <- (sum(referencePoint[,refPointIndex]*referencePoint[,refPointIndex]))^0.5
    intersectionPoint[,refPointIndex] <- intersectionPoint[,refPointIndex]/a
    for(objectiveIndex in 1:populationSize){
      distances[refPointIndex,objectiveIndex] <- norm(matrix(intersectionPoint[,refPointIndex]-objective[,objectiveIndex],nrow=nObjective,ncol = 1),type = "F")
    }
  }
  return(distances) #each column represent different solution, each row different refline
}

# Used in NSGA-III
AssociateLine <- function(distanceFromLines){
  populationSize <- ncol(distanceFromLines)
  associatedLine <- matrix(,nrow = 1,ncol = populationSize)
  for(populationIndex in 1:populationSize){
    associatedLine[populationIndex] <- nnet::which.is.max(-distanceFromLines[,populationIndex])
  }
  return(associatedLine)
}

#' Get Inverted Generational Distance (IGD) value of the population objective w.r.t. a matrix of reference set (each row contain 1 point).
#' @title Get IGD value
#' @param populationObjective The objective value of the corresponding individual
#' @param referenceSet The reference points for computing IGD
#' @return The IGD metric. A Scalar value.
#' @export
GetIGD <- function(populationObjective, referenceSet){
  nRef <- ncol(referenceSet)

  IGDval <- 0
  distances <- CalcIGDDistanceToObjective(populationObjective,referenceSet)
  for(refPointIndex in 1: nRef){
    IGDval <- IGDval + min(distances[refPointIndex,])
  }
  IGDval <- IGDval/nRef
  return(IGDval)
}


ApproximateHypervolumeContribution <- function(populationObjective,referencePoint,numberOfSamples){
  nObjective <- nrow(populationObjective)
  nRef <- length(referencePoint)
  nPop <- ncol(populationObjective)

  nSuccess <- 0
  nUnique <- integer(nPop)

  objectiveAxisMinimum <- matrix(,nrow= nObjective,ncol = 1)
  samplePoint <- rep(0,nObjective)

  for(objectiveIndex in 1:nObjective){
    objectiveAxisMinimum[objectiveIndex] <- min(populationObjective[objectiveIndex,])
  }

  for(sampleIndex in 1:numberOfSamples){
    for (objectiveIndex in 1:nObjective) {
      samplePoint[objectiveIndex] <- stats::runif(1)*(referencePoint[objectiveIndex]-objectiveAxisMinimum[objectiveIndex]) + objectiveAxisMinimum[objectiveIndex]
    }
    successCount <- 0
    for(popIndex in 1:nPop){
      if(min(as.integer(samplePoint >= populationObjective[,popIndex])) == 1 ){
        successCount <- successCount + 1
        successfulPop <- popIndex
        if(successCount==1)
          nSuccess <- nSuccess+1
      }

      if(popIndex==nPop){
        if(successCount == 1){
          nUnique[successfulPop] <- nUnique[successfulPop]+1
        }
      }
    }
  }

  boundingBoxLengths <- referencePoint-objectiveAxisMinimum
  boundingBoxVolume <- 1

  for(objectiveIndex in 1:nObjective){
    boundingBoxVolume <- boundingBoxVolume * boundingBoxLengths[objectiveIndex]
  }
  contribution <- nUnique / numberOfSamples * boundingBoxVolume
  totalHypervolume <- nSuccess / numberOfSamples * boundingBoxVolume

  return(list(contribution,totalHypervolume))
}

#' Get index of the individual with least hypervolume (HV) contribution. For the contribution itself, use GetLeastContribution()
#' @title Get least HV contributor
#' @param populationObjective The objective value of the corresponding individual
#' @param reference The reference point for computing HV
#' @param method the HV computation method
#' @param hypervolumeMethodParam A list of parameters to be passed to the hypervolumeMethod
#' @return The index of the least contributor, an integer.
#' @examples
#' \donttest{
#' nObjective <- 5 # the number of objectives
#' nPoint <- 10 # the number of points that will form the hypervolume
#' objective <- matrix(stats::runif(nObjective*nPoint), nrow = nObjective, ncol = nPoint)
#' GetHypervolume(objective,,"exact") # no reference supplied
#'
#' reference <- rep(2,nObjective) # create a reference point at (2,2,2,2,2)
#' GetLeastContributor(objective,reference,"exact")
#' }
#' @export
GetLeastContributor<- function(populationObjective,reference=NULL,method="exact",hypervolumeMethodParam=list()){
  if(is.null(reference)){
    for(objectiveIndex in 1:nrow(populationObjective))
      reference <- append(reference,max(populationObjective[objectiveIndex,])*1.1)
  }

  if(method=="exact"){
    #hypervolumeContribution <- HVContrib_WFG(populationObjective,reference)
    smallestContributor <- LeastContributorExact(populationObjective,reference)#nnet::which.is.max(-hypervolumeContribution)
  }else{
    string <- paste('LeastContributor',method,'(populationObjective,reference,hypervolumeMethodParam)',sep='')
    smallestContributor <- eval(parse(text=string))
  }

  return(smallestContributor)
}

#' Get the hypervolume (HV) contribution of the individual with least HV contribution.
#' @title Get least HV contribution
#' @param populationObjective The objective value of the corresponding individual
#' @param reference The reference point for computing HV
#' @param method the HV computation method

#' @return The HV contribution value of the least contributor.
#' @examples
#'  \donttest{
#' nObjective <- 5 # the number of objectives
#' nPoint <- 10 # the number of points that will form the hypervolume
#' objective <- matrix(stats::runif(nObjective*nPoint), nrow = nObjective, ncol = nPoint)
#' GetHypervolume(objective,,"exact") # no reference supplied
#'
#' reference <- rep(2,nObjective) # create a reference point at (2,2,2,2,2)
#' GetLeastContribution(objective,reference,"exact")
#' }
#' @export
GetLeastContribution<- function(populationObjective,reference=NULL,method="exact"){
  if(method=="exact"){
    if(is.null(reference))
      hypervolumeContribution <- HVContrib_WFG(populationObjective)
    else
      hypervolumeContribution <- HVContrib_WFG(populationObjective, reference)

    smallestContribution <- min(hypervolumeContribution)
  }
  if(method=="approx"){
    smallestContribution <- LeastContributionApprox(populationObjective,reference)
  }

  return(smallestContribution)
}

#' Get the hypervolume (HV) contribution of the population. Dominated front will give 0 contribution.
#' @title Get HV contribution of all points.
#' @param populationObjective The objective value of the corresponding individual
#' @param reference The reference point for computing HV
#' @param method the HV computation method. Currently ignored and uses the WFG exact method.
#' @return A vector of length ncol(populationObjective)
#' @examples
#' \donttest{
#' nObjective <- 5 # the number of objectives
#' nPoint <- 10 # the number of points that will form the hypervolume
#' objective <- matrix(stats::runif(nObjective*nPoint), nrow = nObjective, ncol = nPoint)
#' GetHypervolume(objective,,"exact") # no reference supplied
#'
#' reference <- rep(2,nObjective) # create a reference point at (2,2,2,2,2)
#' GetHVContribution(objective,reference)
#' }
#' @export
GetHVContribution<- function(populationObjective,reference=NULL,method="exact"){
  #  if(method=="exact"){
  if (!pkg.globals$have_pygmo)
    stop("HV computation requires PyGMO")

  hypervolumeContribution <- HVContrib_WFG(populationObjective, reference)

  return(hypervolumeContribution)
}

#' Compute the hypervolume formed by the points w.r.t. a reference point. If no reference is supplied, use the nadir point*(1.1,...,1.1).
#' @title Compute hypervolume
#' @param objective The set of points in the objective space (The objective values). A single column should contain one point, so the size would be numberOfObjective x nPoint, e.g. in 5 objective problem, it is 5 x n.
#' @param reference The reference points. Each column represent one point. Size: numberOfObjective x nPoint, e.g. in 5 objective problem, it is 5 x n.
#' @param method Exact using WFG method or approximate HV using the method by Bringmann and Friedrich. Default to "exact".
#'
#' @return Hypervolume size, a scalar value.
#' @examples
#' \donttest{
#' nObjective <- 5 # the number of objectives
#' nPoint <- 10 # the number of points that will form the hypervolume
#' objective <- matrix(stats::runif(nObjective*nPoint), nrow = nObjective, ncol = nPoint)
#' GetHypervolume(objective,,"exact") # no reference supplied
#'
#' reference <- rep(2,nObjective) # create a reference point at (2,2,2,2,2)
#' GetHypervolume(objective,reference,"exact") # using reference point
#' }
#' @export
#'
GetHypervolume <- function(objective,reference=NULL,method="exact"){
  if (!pkg.globals$have_pygmo)
    stop("HV computation requires PyGMO")

  if(method=="exact"){
    if(is.null(reference))
      hypervolume <- HypervolumeExact(objective)
    else{
      hypervolume <- HypervolumeExact(objective, reference)
    }
  }
  if(method=="approx"){
    hypervolume <- HypervolumeBFapprox(objective,reference)
  }
  return(hypervolume)
}

