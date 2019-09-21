LeastContributorExact <- function(populationObjective,reference=NULL){
  if(is.null(reference)){
    for(objectiveIndex in 1:nrow(populationObjective))
      reference<-append(reference,max(populationObjective[objectiveIndex,])*1.1)
  }

  hv <- pkg.globals$pygmo$hypervolume(t(populationObjective))
  hv_cont_alg <- pkg.globals$pygmo$bf_approx(TRUE,1L,0.01,0.05)


  leastContributor <- hv$least_contributor(reference,hv_cont_alg)
  return(leastContributor+1)
}

LeastContributorBFapprox <- function(populationObjective,reference=NULL){
  if(is.null(reference)){
    for(objectiveIndex in 1:nrow(populationObjective))
      reference<-append(reference,max(populationObjective[objectiveIndex,])*1.1)
  }

  hv <- pkg.globals$pygmo$hypervolume(t(populationObjective))
  hv_cont_alg <- pkg.globals$pygmo$bf_approx(TRUE,1L,0.01,0.05)


  leastContributor <- hv$least_contributor(reference,hv_cont_alg)
  return(leastContributor+1)
}

LeastContributionApprox <- function(populationObjective,reference=NULL){
  if(is.null(reference)){
    for(objectiveIndex in 1:nrow(populationObjective))
      reference<-append(reference,max(populationObjective[objectiveIndex,])*1.1)
  }

  hv <- pkg.globals$pygmo$hypervolume(t(populationObjective))
  contribution <- hv$contributions(reference)

  return(min(contribution))
}

HypervolumeBFapprox <- function(populationObjective,reference=NULL){
  if(is.null(reference)){
    for(objectiveIndex in 1:nrow(populationObjective))
      append(reference,max(populationObjective[objectiveIndex,])*1.1)
  }

  hv <- pkg.globals$pygmo$hypervolume(t(populationObjective))
  hv_alg <- pkg.globals$pygmo$bf_fpras(0.05,0.05,as.integer(1000))

  approxHv <- hv$compute(reference,hv_alg)
}

HypervolumeExact <- function(populationObjective,reference=NULL){
  if(is.null(reference)){
    for(objectiveIndex in 1:nrow(populationObjective)){
      reference <- append(reference,max(populationObjective[objectiveIndex,])*1.1)
    }
  }
  hv <- pkg.globals$pygmo$hypervolume(t(populationObjective))
  algo <- pkg.globals$pygmo$hvwfg()

  hypervolume <- hv$compute(reference,algo)
}


HVContrib_WFG <- function(populationObjective,reference=NULL){
  if(is.null(reference)){
    for(objectiveIndex in 1:nrow(populationObjective))
      reference <- append(reference,max(populationObjective[objectiveIndex,])*1.1)
  }

  hv <- pkg.globals$pygmo$hypervolume(t(populationObjective))
  algo <- pkg.globals$pygmo$hvwfg()
  if(is.matrix(reference)){
    reference <- reference[,]
  }
  hvContrib <- hv$contributions(reference,algo)
  return(hvContrib)
}
