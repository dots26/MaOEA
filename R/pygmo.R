LeastContributorExact <- function(populationObjective,reference=NULL,ref_multiplier=1.1){
  if(is.vector(populationObjective))
    populationObjective <- matrix(populationObjective)
  if(is.null(reference)){ # use dynamic reference
    hv <- pkg.globals$pygmo$hypervolume(t(populationObjective))
    hv_cont_alg <- pkg.globals$pygmo$hvwfg()
    hv_forall <- hv$compute(reference,hv_cont_alg)

    popSize <- ncol(populationObjective)
    contrib <- NULL
    for(i in 1:popSize){
      tempPop <- populationObjective[,-i]
      for(objectiveIndex in 1:nrow(populationObjective))
        reference<-append(reference,max(tempPop[objectiveIndex,])*ref_multiplier)
      subhv_object <- pkg.globals$pygmo$hypervolume(t(tempPop))
      hv_subhv <- subhv_object$compute(reference,hv_cont_alg)

      contrib <- append(contrib,hv_forall-hv_subhv)
    }
    leastContributor <- which.min(contrib)
  }else{
    # check reference is dominated by all points
    rmIndex <- NULL
    for(pointIndex in 1:ncol(populationObjective)){
      if(any(populationObjective[,pointIndex]>reference)){
        rmIndex <- append(rmIndex,pointIndex)
        populationObjective[,pointIndex] <- reference - .Machine$double.eps*reference*10*sign(reference)
      }
    }

    if(!is.null(rmIndex)){
      warning("Some points are dominated by the reference and ignored")
    }
    hv <- pkg.globals$pygmo$hypervolume(t(populationObjective))
    hv_cont_alg <- pkg.globals$pygmo$hvwfg()
    if(is.matrix(reference)){
      reference <- reference[,]
    }
    leastContributor <- hv$least_contributor(reference,hv_cont_alg)
  }

  return(leastContributor+1)
}

LeastContributorBFapprox <- function(populationObjective,reference=NULL,ref_multiplier=1.1){
  if(is.vector(populationObjective))
    populationObjective <- matrix(populationObjective)

  if(is.null(reference)){
    hv <- pkg.globals$pygmo$hypervolume(t(populationObjective))
    hv_cont_alg <- pkg.globals$pygmo$bf_fpras()
    hv_forall <- hv$compute(reference,hv_cont_alg)

    popSize <- ncol(populationObjective)
    contrib <- NULL
    for(i in 1:popSize){
      tempPop <- populationObjective[,-i]
      for(objectiveIndex in 1:nrow(populationObjective))
        reference<-append(reference,max(tempPop[objectiveIndex,])*ref_multiplier)
      subhv_object <- pkg.globals$pygmo$hypervolume(t(tempPop))
      hv_subhv <- subhv_object$compute(reference,hv_cont_alg)

      contrib <- append(contrib,hv_forall-hv_subhv)
    }
    leastContributor <- which.min(contrib)
  }else{
    # check reference is dominated by all points
    rmIndex <- NULL
    for(pointIndex in 1:ncol(populationObjective)){
      if(any(populationObjective[,pointIndex]>reference)){
        rmIndex <- append(rmIndex,pointIndex)
        populationObjective[,pointIndex] <- reference - .Machine$double.eps*reference*10*sign(reference)
      }
    }
    if(!is.null(rmIndex)){
      warning("Some points are dominated by the reference and ignored")
    }
    hv <- pkg.globals$pygmo$hypervolume(t(populationObjective))
    hv_cont_alg <- pkg.globals$pygmo$bf_approx(TRUE,1L,0.01,0.05)

    leastContributor <- hv$least_contributor(reference,hv_cont_alg)
  }

  return(leastContributor+1)
}

LeastContributionApprox <- function(populationObjective,reference=NULL,ref_multiplier=1.1){
  if(is.vector(populationObjective))
    populationObjective <- matrix(populationObjective)
  if(is.null(reference)){
    for(objectiveIndex in 1:nrow(populationObjective))
      reference<-append(reference,max(populationObjective[objectiveIndex,])*ref_multiplier)
  }else{
    # check reference is dominated by all points
    rmIndex <- NULL
    for(pointIndex in 1:ncol(populationObjective)){
      if(any(populationObjective[,pointIndex]>reference)){
        rmIndex <- append(rmIndex,pointIndex)
        populationObjective[,pointIndex] <- reference - .Machine$double.eps*reference*10*sign(reference)
      }
    }
    if(!is.null(rmIndex)){
      warning("Some points are dominated by the reference and ignored")
    }
    # populationObjective <- populationObjective[,-rmIndex]
  }

  hv <- pkg.globals$pygmo$hypervolume(t(populationObjective))
  contribution <- hv$contributions(reference)

  return(min(contribution))
}

HypervolumeBFapprox <- function(populationObjective,reference=NULL,ref_multiplier=1.1){
  if(is.vector(populationObjective))
    populationObjective <- matrix(populationObjective)

  if(ncol(populationObjective)==0)
    return(0)

  if(is.null(reference)){
    for(objectiveIndex in 1:nrow(populationObjective))
      append(reference,max(populationObjective[objectiveIndex,])*ref_multiplier)
  }else{
    # check reference is dominated by all points
    rmIndex <- NULL
    for(pointIndex in 1:ncol(populationObjective)){
      if(any(populationObjective[,pointIndex]>reference)){
        rmIndex <- append(rmIndex,pointIndex)
        populationObjective[,pointIndex] <- reference - .Machine$double.eps*reference*10*sign(reference)
      }
    }
    if(!is.null(rmIndex)){
      warning("Some points are dominated by the reference and ignored")
    }
    # populationObjective <- populationObjective[,-rmIndex]
  }

  hv <- pkg.globals$pygmo$hypervolume(t(populationObjective))
  hv_alg <- pkg.globals$pygmo$bf_fpras(0.05,0.05,as.integer(1000))

  approxHv <- hv$compute(reference,hv_alg)
}

HypervolumeExact <- function(populationObjective,reference=NULL,ref_multiplier=1.1){
  if(is.vector(populationObjective))
    populationObjective <- matrix(populationObjective)

  if(ncol(populationObjective)==0)
    return(0)

  if(is.null(reference)){
    for(objectiveIndex in 1:nrow(populationObjective)){
      reference <- append(reference,max(populationObjective[objectiveIndex,])*ref_multiplier)
    }
  }else{
    # check reference is dominated by all points
    rmIndex <- NULL
    for(pointIndex in 1:ncol(populationObjective)){
      if(any(populationObjective[,pointIndex]>reference)){
        rmIndex <- append(rmIndex,pointIndex)
        populationObjective[,pointIndex] <- reference - .Machine$double.eps*reference*10*sign(reference)
      }
    }
    if(!is.null(rmIndex)){
      warning("Some points are dominated by the reference and ignored")
    }
    # populationObjective <- populationObjective[,-rmIndex]
  }


  hv <- pkg.globals$pygmo$hypervolume(t(populationObjective))
  algo <- pkg.globals$pygmo$hvwfg()

  hypervolume <- hv$compute(reference,algo)
}


HVContrib_WFG <- function(populationObjective,reference=NULL,ref_multiplier=1.1){
  if(is.vector(populationObjective))
    populationObjective <- matrix(populationObjective)

  if(is.null(reference)){
    for(objectiveIndex in 1:nrow(populationObjective))
      reference <- append(reference,max(populationObjective[objectiveIndex,])*ref_multiplier)
  }else{
    # check reference is dominated by all points
    rmIndex <- NULL
    for(pointIndex in 1:ncol(populationObjective)){
      if(any(populationObjective[,pointIndex]>reference)){
        rmIndex <- append(rmIndex,pointIndex)
        populationObjective[,pointIndex] <- reference - .Machine$double.eps*reference*10*sign(reference)
      }
    }
    if(!is.null(rmIndex)){
      warning("Some points are dominated by the reference and ignored")
    }
  }

  hv <- pkg.globals$pygmo$hypervolume(t(populationObjective))
  algo <- pkg.globals$pygmo$hvwfg()
  if(is.matrix(reference)){
    reference <- reference[,]
  }
  hvContrib <- hv$contributions(reference,algo)
  return(hvContrib)
}
