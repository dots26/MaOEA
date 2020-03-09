LeastContributorExact <- function(populationObjective,reference=NULL){
  if(is.vector(populationObjective))
    populationObjective <- matrix(populationObjective)
  if(is.null(reference)){
    for(objectiveIndex in 1:nrow(populationObjective))
      reference<-append(reference,max(populationObjective[objectiveIndex,])*1.1)
  }else{
    # check reference is dominated by all points
    rmIndex <- NULL
    for(pointIndex in 1:ncol(populationObjective)){
      if(any(populationObjective[,pointIndex]>reference)){
        rmIndex <- append(rmIndex,pointIndex)
        populationObjective[,pointIndex] <- reference - .Machine$double.eps*reference
      }
    }
    if(!is.null(rmIndex)){
      warning("Some points are dominated by the reference and ignored")
    }
    # populationObjective <- populationObjective[,-rmIndex]
  }

  hv <- pkg.globals$pygmo$hypervolume(t(populationObjective))
  hv_cont_alg <- pkg.globals$pygmo$bf_approx(TRUE,1L,0.01,0.05)


  leastContributor <- hv$least_contributor(reference,hv_cont_alg)
  return(leastContributor+1)
}

LeastContributorBFapprox <- function(populationObjective,reference=NULL){
  if(is.vector(populationObjective))
    populationObjective <- matrix(populationObjective)
  if(is.null(reference)){
    for(objectiveIndex in 1:nrow(populationObjective))
      reference<-append(reference,max(populationObjective[objectiveIndex,])*1.1)
  }else{
    # check reference is dominated by all points
    rmIndex <- NULL
    for(pointIndex in 1:ncol(populationObjective)){
      if(any(populationObjective[,pointIndex]>reference)){
        rmIndex <- append(rmIndex,pointIndex)
        populationObjective[,pointIndex] <- reference - .Machine$double.eps*reference
      }
    }
    if(!is.null(rmIndex)){
      warning("Some points are dominated by the reference and ignored")
    }
    # populationObjective <- populationObjective[,-rmIndex]
  }

  hv <- pkg.globals$pygmo$hypervolume(t(populationObjective))
  hv_cont_alg <- pkg.globals$pygmo$bf_approx(TRUE,1L,0.01,0.05)


  leastContributor <- hv$least_contributor(reference,hv_cont_alg)
  return(leastContributor+1)
}

LeastContributionApprox <- function(populationObjective,reference=NULL){
  if(is.vector(populationObjective))
    populationObjective <- matrix(populationObjective)
  if(is.null(reference)){
    for(objectiveIndex in 1:nrow(populationObjective))
      reference<-append(reference,max(populationObjective[objectiveIndex,])*1.1)
  }else{
    # check reference is dominated by all points
    rmIndex <- NULL
    for(pointIndex in 1:ncol(populationObjective)){
      if(any(populationObjective[,pointIndex]>reference)){
        rmIndex <- append(rmIndex,pointIndex)
        populationObjective[,pointIndex] <- reference - .Machine$double.eps*reference
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

HypervolumeBFapprox <- function(populationObjective,reference=NULL){
  if(is.vector(populationObjective))
    populationObjective <- matrix(populationObjective)
  if(is.null(reference)){
    for(objectiveIndex in 1:nrow(populationObjective))
      append(reference,max(populationObjective[objectiveIndex,])*1.1)
  }else{
    # check reference is dominated by all points
    rmIndex <- NULL
    for(pointIndex in 1:ncol(populationObjective)){
      if(any(populationObjective[,pointIndex]>reference)){
        rmIndex <- append(rmIndex,pointIndex)
        populationObjective[,pointIndex] <- reference - .Machine$double.eps*reference
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

HypervolumeExact <- function(populationObjective,reference=NULL){
  if(is.vector(populationObjective))
    populationObjective <- matrix(populationObjective)
  if(is.null(reference)){
    for(objectiveIndex in 1:nrow(populationObjective)){
      reference <- append(reference,max(populationObjective[objectiveIndex,])*1.1)
    }
  }else{
    # check reference is dominated by all points
    rmIndex <- NULL
    for(pointIndex in 1:ncol(populationObjective)){
      if(any(populationObjective[,pointIndex]>reference)){
        rmIndex <- append(rmIndex,pointIndex)
        populationObjective[,pointIndex] <- reference - .Machine$double.eps*reference
      }
    }
    if(!is.null(rmIndex)){
      warning("Some points are dominated by the reference and ignored")
    }
  }


  hv <- pkg.globals$pygmo$hypervolume(t(populationObjective))
  algo <- pkg.globals$pygmo$hvwfg()

  hypervolume <- hv$compute(reference,algo)
}


HVContrib_WFG <- function(populationObjective,reference=NULL,outsideReferenceHandling=c('zero','penalty')){
  if(is.vector(populationObjective))
    populationObjective <- matrix(populationObjective)
  if(is.null(reference)){
    for(objectiveIndex in 1:nrow(populationObjective))
      reference <- append(reference,max(populationObjective[objectiveIndex,])*1.1)
  }else{
    # check reference is dominated by all points
    rmIndex <- NULL
    for(pointIndex in 1:ncol(populationObjective)){
      if(any(populationObjective[,pointIndex]>reference)){
        rmIndex <- append(rmIndex,pointIndex)
        if(outsideReferenceHandling=='zero')
          populationObjective[,pointIndex] <- reference - .Machine$double.eps*reference
      }
    }
    if(!is.null(rmIndex)){
      warning("Some points are dominated by the reference and ignored")
    }
  }
  if(outsideReferenceHandling=='zero'){
    hv <- pkg.globals$pygmo$hypervolume(t(populationObjective))
    algo <- pkg.globals$pygmo$hvwfg()
    if(is.matrix(reference)){
      reference <- reference[,]
    }
    hvContrib <- hv$contributions(reference,algo)
  }else if(outsideReferenceHandling='penalty'){
    hv <- pkg.globals$pygmo$hypervolume(t(populationObjective[,-rmIndex]))
    algo <- pkg.globals$pygmo$hvwfg()
    if(is.matrix(reference)){
      reference <- reference[,]
    }
    # original HV when all objs are used
    hv_pre <- hv$compute(reference,algo)
    hvContrib <- vector()

    for(i in 1:ncol(populationObjective)){
      if(i==any(rmIndex)){
        new_ref <- vector()
        for(objIndex in 1:nrow(populationObjective)){
          new_ref[objIndex] <- max(populationObjective[objIndex,i],reference[objIndex])
        }
        hv_new <- hv$compute(new_ref,algo)
        hvContrib[i] <- hv_pre - hv_new # penalized contribution
      }else{
        hv_obj_new <- pkg.globals$pygmo$hypervolume(t(populationObjective[,c(-rmIndex,-i)]))
        hv_new <- hv$compute(reference,algo)
        hvContrib[i] <- hv_pre - hv_new  # normal contribution
      }
    }
  }
  return(hvContrib)
}
