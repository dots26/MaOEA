#' Compute the R2-mtch indicator from Shang et al.
#' @title Modified tchebyscheff R2-indicator
#' @param dataPoints The Points coordinate. Each column contains a single point (column major).
#' @param reference The reference point for computing R2-mtch (similar as reference for HV)
#' @param weights The weights/direction to be used to compute the achievement scalarization. Each column contains a single weight vector. If no weight is supplied, weights are generated using Sobol sequences.
#' @param nWeight Used only when no weights are supplied. An input for the sobol weight generation. This defines how many points to be generated.
#' @references Ke Shang, Hisao Ishibuchi, Min-Ling Zhang, and Yiping Liu. 2018. A new R2 indicator for better hypervolume approximation. In Proceedings of the Genetic and Evolutionary Computation Conference (GECCO '18), Hernan Aguirre (Ed.). ACM, New York, NY, USA, 745-752. DOI: https://doi.org/10.1145/3205455.3205543
#' @return The function return the R2-indicator of the set.
#' @examples
#' nPointToSample <- 100
#' nObjective <- 3
#' points <- matrix(runif(nPointToSample*nObjective), nrow = nObjective) # sample the points
#' ranks <- nsga2R::fastNonDominatedSorting(t(points)) # non-dominated sorting
#' points <- points[,ranks[[1]],drop=FALSE] # take only the non-dominated front
#' nPoints <- ncol(points) # check how many points are on the non-dominated front
#' reference <- rep(2,nObjective)
#'
#' compute_R2mtch(points,reference)
#' @export
compute_R2mtch <- function(dataPoints,reference,weights=NULL,nWeight = 100){
  nObj <- nrow(dataPoints)
  if(is.null(weights)) {
    weights <- createWeightsSobol(nDim=nObj,nWeights = nWeight)
  }
  #  if(is.null(weights)) {
  #    weights <- createWeights(nObj,axisDivision,noZero = TRUE)
  #  }
  sumR2 <- 0

  nWeight <- ncol(weights)
  nPoints <- ncol(dataPoints)
  for(weightIndex in 1:nWeight){
    pointAchievementVector <- NULL
    for(pointIndex in 1:nPoints){
      pointAchievementVector <- append(pointAchievementVector,
                                       gmtch(dataPoints[,pointIndex],reference,weights[,weightIndex])
      )
    }
    sumR2 <- sumR2+(max(pointAchievementVector))
  }
  return(sumR2/nWeight)
}

#' Compute the R2-HV from Shang et al.
#' @title Modified powered tchebyscheff R2-indicator designed to approximate HV
#' @param dataPoints The Points coordinate. Each column contains a single point (column major).
#' @param reference The reference point for computing R2-mtch (similar as reference for HV)
#' @param weights The weights/direction to be used to compute the achievement scalarization. Each column contains a single weight vector. If no weight is supplied, weights are generated using Sobol sequences.
#' @param nPoints Used only when no weights are supplied. An input for the weight generator (sobol sequences). This defines how many points are created.
#' @references Ke Shang, Hisao Ishibuchi, Min-Ling Zhang, and Yiping Liu. 2018. A new R2 indicator for better hypervolume approximation. In Proceedings of the Genetic and Evolutionary Computation Conference (GECCO '18), Hernan Aguirre (Ed.). ACM, New York, NY, USA, 745-752. DOI: https://doi.org/10.1145/3205455.3205543
#' @return The function return the powered R2-indicator of the set.
#' @examples
#' nPointToSample <- 100
#' nObjective <- 3
#' points <- matrix(runif(nPointToSample*nObjective), nrow = nObjective) # sample the points
#' ranks <- nsga2R::fastNonDominatedSorting(t(points)) # non-dominated sorting
#' points <- points[,ranks[[1]],drop=FALSE] # take only the non-dominated front
#' nPoints <- ncol(points) # check how many points are on the non-dominated front
#' reference <- rep(2,nObjective)
#'
#' compute_R2HV(points,reference)
#' @export
compute_R2HV <- function(dataPoints,reference,weights=NULL,nPoints = 100){
  nObj <- nrow(dataPoints)
  if(is.null(weights)) {
    weights <- createWeightsSobol(nDim=nObj,nWeights = nPoints)
  }
  #if(is.null(weights)) {
  #  weights <- createWeights(nObj,axisDivision,noZero = TRUE)
  #}
  sumR2 <- 0

  nWeight <- ncol(weights)
  nPoints <- ncol(dataPoints)
  for(weightIndex in 1:nWeight){
    pointAchievementBest <- -Inf
    for(pointIndex in 1:nPoints){
      this.gmtch <- gmtch(dataPoints[,pointIndex],reference,weights[,weightIndex])^nObj
      if(pointAchievementBest< this.gmtch)
        pointAchievementBest <- this.gmtch
    }
    sumR2 <- sumR2+pointAchievementBest
  }
  return(sumR2/nWeight)
}

#' Compute the R2-HVC from Shang et al.
#' @title Modified tchebyscheff R2-indicator contribution designed to approximate HV
#' @param dataPoints The Points coordinate. Each column contains a single point (column major).
#' @param reference The reference point for computing R2-mtch (similar as reference for HV)
#' @param weights The weights/direction to be used to compute the achievement scalarization. Each column contains a single weight vector. If no weight is supplied, weights are generated using Sobol sequences
#' @param nWeight Used only when no weights are supplied. The number of weights generated by sobol sequence.
#' @param alpha Power factor on the gmtch and g*2tch utility functions.
#' @param indexOfInterest individuals to be evaluated. The R2 values will only be reported/returned for these individuals.
#' @references K. Shang, H. Ishibuchi and X. Ni, "R2-based Hypervolume Contribution Approximation," in IEEE Transactions on Evolutionary Computation. doi: 10.1109/TEVC.2019.2909271
#' @return The function return R2-indicator contribution of each point.
#' @examples
#' nPointToSample <- 100
#' nObjective <- 3
#' points <- matrix(runif(nPointToSample*nObjective), nrow = nObjective) # sample the points
#' ranks <- nsga2R::fastNonDominatedSorting(t(points)) # non-dominated sorting
#' points <- points[,ranks[[1]],drop=FALSE] # take only the non-dominated front
#' nPoints <- ncol(points) # check how many points are on the non-dominated front
#' reference <- rep(2,nObjective)
#'
#' compute_R2HVC(points,reference)
#' @export
compute_R2HVC <- function(dataPoints,reference,weights=NULL,alpha=1,nWeight = 300,indexOfInterest = 1:ncol(dataPoints)){
  nObj <- nrow(dataPoints)
  if(is.null(weights)) {
    weights <- createWeightsSobol(nDim=nObj,nWeights = nWeight)
  }
  sumR2 <- 0
  R2contrib <- NULL
  log_R2 <- NULL
  logsd_R2 <- NULL
  skew_R2 <- NULL
  nWeight <- ncol(weights)
  nPoints <- ncol(dataPoints)
  for(sIndex in indexOfInterest){
    sumR2 <- 0
    minRset <- NULL
    for(weightIndex in 1:nWeight){
      minimumStar <- Inf
      pointAchievementToBoundary <- gmtch(reference,dataPoints[,sIndex],weights[,weightIndex])

      for(secondaryPointIndex in 1:nPoints){
        if( secondaryPointIndex != sIndex){
          new_g2tch <- g2tch_star(dataPoints[,secondaryPointIndex],dataPoints[,sIndex],weights[,weightIndex])

          if(minimumStar > new_g2tch ){
            minimumStar <- new_g2tch
          }
        }
      }
      minR <- min(c(minimumStar,pointAchievementToBoundary))
      minR <- minR^alpha
      minRset <- append(minRset,minR)
    }
    skew_R2 <- append(skew_R2,skewness(log(minRset)))
    R2contrib <- append(R2contrib,mean(minRset))
    log_R2 <- append(log_R2,mean(log(minRset)))
    logsd_R2 <- append(logsd_R2,stats::sd(log(minRset)))
  }
  return(list(R2=R2contrib,log_R2=log_R2,logsd_R2=logsd_R2,skew_R2=skew_R2))
}

gmtch <- function(point, reference, weight){
  nObj <- nrow(point)
  achievement_vector <- abs(point-reference)/weight

  return(min(achievement_vector))
}


g2tch <- function(point, ref, weight){
  nObj <- nrow(point)
  achievement_vector <- abs(point-ref)/weight
  return(max(achievement_vector))
}

g2tch_star <- function(a, s, weight){
  nObj <- nrow(a)
  achievement_vector <- (a-s)/weight


  return(max(achievement_vector))
}

#' Generate a set of weights following Das and Dennis's method. Each column returned is a weight vector.
#' @title Das and Dennis's structured weight generation, normal boundary intersection (NBI).
#' @param nDim The dimensionality of the problem. In EA, usually this is used in the objective space, hence nDim = nObjective
#' @param axisDivision Used only when no weights are supplied. An input for the structured weight distribution. This defines how many division are created in each axis.
#' @param noZero Default to false. If set to TRUE, reference vector containing zero, e.g. (1,0,0) will be removed. Used to generate weight in modified tch method.
#' @references Indraneel Das and J. E. Dennis. 1998. Normal-Boundary Intersection: A New Method for Generating the Pareto Surface in Nonlinear Multicriteria Optimization Problems. SIAM Journal on Optimization 1998 8:3, 631-657.
#' @return The function return a set of weight vectors.
#' @examples
#' nObjective <- 3
#' axisDiv <- 6
#'
#' createWeights(nObjective,axisDiv)
#' @export
createWeights <- function(nDim,axisDivision = nDim+2,noZero=FALSE){
  #REFERENCE POINT DEFINITION
  nRefPoint <- choose(nDim+axisDivision-1,axisDivision)

  referencePoint <- matrix(,nrow = nDim,ncol = nRefPoint)
  referenceCombination <- gtools::combinations(nDim,axisDivision,1:nDim,repeats.allowed = TRUE)
  for(objectiveIndex in 1:nDim){
    this.objective <- as.integer(referenceCombination==objectiveIndex)
    this.objective <- matrix(this.objective,nrow = nRefPoint,ncol = axisDivision)
    for(referencePointIndex in 1:nRefPoint)
      referencePoint[objectiveIndex,referencePointIndex] <- sum(this.objective[referencePointIndex,])
  }

  if(noZero){
    if(axisDivision <= nDim){
      stop('Number of axis division must be larger than problem dimension, otherwise no internal point will be created.')
    }
    removeCol <- NULL
    for(refIndex in 1:nRefPoint){
      if(min(referencePoint[,refIndex])==0){
        removeCol <- append(removeCol,refIndex)
      }
    }
    referencePoint <- referencePoint[,-removeCol]
  }

  nRefPoint <- ncol(referencePoint)
  for(refIndex in 1:nRefPoint){
    referencePoint[,refIndex] <- referencePoint[,refIndex]/norm(referencePoint[,refIndex,drop=FALSE],'F')
  }

  return(referencePoint)
}

#' Generate a set of weights following Sobol sequence generator
#' @title Sobol sequence weights
#' @param nWeights Number of weights to generate.
#' @param nDim The dimensionality of the problem. In EA, usually this is used in the objective space, hence nDim = nObjective
#' @param seed Seed for scrambling
#' @return The function return a set of weight vectors.
#' @examples
#' nObjective <- 3
#' nPoint <- 1000
#'
#' createWeightsSobol(nPoint,nObjective)
#' @export
createWeightsSobol <- function(nWeights, nDim,seed=4177){
  weights <- randtoolbox::sobol(nWeights, dim = nDim,scrambling=3,seed = seed)
  weights <- t(weights)
  for(pointIndex in 1:nWeights){
    weights[,pointIndex] <- weights[,pointIndex]/norm(weights[,pointIndex,drop=FALSE],'F')
  }

  return(weights)
}
