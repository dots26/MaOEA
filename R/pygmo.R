#pygmo <- import("pygmo")
#rndGen <- import("numpy")
#garbage <- import("gc")
#p<-1
#rndGen$random$seed(as.integer(p*1000))


LeastContributorBFapprox <- function(populationObjective,reference=NULL){
  if(is.null(reference)){
    for(objectiveIndex in 1:nrow(populationObjective))
      reference<-append(reference,max(populationObjective[objectiveIndex,])+1)
  }

  hv <- pygmo$hypervolume(t(populationObjective))
  hv_cont_alg <- pygmo$bf_approx(TRUE,1L,0.01,0.05)


  leastContributor <- hv$least_contributor(reference,hv_cont_alg)
  return(leastContributor+1)
}

LeastContributionApprox <- function(populationObjective,reference=NULL){
  if(is.null(reference)){
    for(objectiveIndex in 1:nrow(populationObjective))
      reference<-append(reference,max(populationObjective[objectiveIndex,])+1)
  }

  hv <- pygmo$hypervolume(t(populationObjective))
  contribution <- hv$contributions(reference)

  return(min(contribution))
}

HypervolumeBFapprox <- function(populationObjective,reference=NULL){
  if(is.null(reference)){
    for(objectiveIndex in 1:nrow(populationObjective))
      append(reference,max(populationObjective[objectiveIndex,])+1)
  }

  hv <- pygmo$hypervolume(t(populationObjective))
  hv_alg <- pygmo$bf_fpras(0.05,0.05,as.integer(p*1000))

  approxHv <- hv$compute(reference,hv_alg)
}

HypervolumeExact <- function(populationObjective,reference=NULL){
  if(is.null(reference)){
    for(objectiveIndex in 1:nrow(populationObjective)){
      reference <- append(reference,max(populationObjective[objectiveIndex,])+1)
    }
  }
  hv <- pygmo$hypervolume(t(populationObjective))
  hypervolume <- hv$compute(reference)
}


HVContrib_WFG <- function(populationObjective,reference=NULL){
  if(is.null(reference)){
    for(objectiveIndex in 1:nrow(populationObjective))
      reference <- append(reference,max(populationObjective[objectiveIndex,])+1)
  }
  #totalHV <- HypervolumeExact(populationObjective,reference[,,drop=TRUE])
  #hvContrib <- matrix(,ncol=0,nrow = 1)
  #for(i in 1:ncol(populationObjective)){
  #  remHV <- HypervolumeExact(populationObjective[,-i],reference[,,drop=TRUE])
  #  hvContrib <- cbind(hvContrib,totalHV-remHV)
  #}
  hv <- pygmo$hypervolume(t(populationObjective))
  algo <- pygmo$hvwfg()
  if(is.matrix(reference)){
    reference <- reference[,]
  }
  hvContrib <- hv$contributions(reference,algo)
  #print(hvContrib)
  return(hvContrib)
}
