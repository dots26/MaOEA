optimMaOEA <- function(problem,
                       optimizer,
                       nObjective,
                       nVar,
                       populationSize,
                       referencePoint,
                       hypervolumeTarget=0.99,
                       crossoverProbability=1,
                       crossoverDistribution=30,
                       mutationProbability=1,
                       mutationDistribution=20,
                       successProbThreshold = 0.44,
                       successProbTarget = 0.5){
  # check which problem are being solved
  if(!is.na(stringr::str_match(problem,"DTLZ"))){
    if(!is.na(stringr::str_match(problem,"DTLZ1")))
      problemType <- 1
    else
      problemType <- 2
  }
  if(!is.na(stringr::str_match(problem,"WFG"))){
    problemType <- 3
  }
  #initialize population
  if((problemType==3) && (WFGScaling==TRUE)){
    population<-InitializePopulationWFG(populationSize,nVar,FALSE,0,1)
  }else{
    population<-InitializePopulationLHS(populationSize,nVar,FALSE,0,1)
  }

  intersectionPoint <- referencePoint
  # get the intersection of the refline and sigma(x^2)=1 (for the concave problems, slightly different for DTLZ1)
  # x^2  + y^2  + z^2  = 1
  # the refline is the line from (0,0,0) to some {(x,y,z) | x+y+z=1} (the ref point)
  # which means if for example the ref point is at (x*,y*,z*)
  # the line has length = a = root(x*^2+y*^2+z*^2)
  # to make the line to have length = 1, we simply multiply the initial ref point by 1/a
  #
  intersectionPoint <- referencePoint
  nRefPoint <- ncol(referencePoint)
  for (refPointIndex in 1:nRefPoint){
    if(problemType == 1) #DTLZ1
      a <- (sum(referencePoint[,refPointIndex]))*2 # for DTLZ1
    else
      a <- (sum(referencePoint[,refPointIndex]*referencePoint[,refPointIndex]))^0.5 # for other DTLZ and WFG
    intersectionPoint[,refPointIndex] <- intersectionPoint[,refPointIndex]/a
  }
  origin <- integer(nObjective)

  targetSmetric <- hypervolumeTarget * GetHypervolume(intersectionPoint,rep(2,nObjective),"exact")
  Smetric <- 0 # s-metric is the hypervolume
  parent_list <- cmaes_gen(population)

  while(Smetric < (targetSmetric)){ # use this if there is a Hypervolume target
    populationObjective<-matrix(,nrow=nObjective,ncol=0)
    if((optimizer=="NSGA"))
      NSGA3(nObjective, population,referencePoint,problem,crossoverProbability,mutationProbability,TRUE,mutationDistribution,crossoverDistribution)
    else if((optimizer=="SCMAES") ||(optimizer=="CMAES")) {

      if(optimizer=="SCMAES")
        newGeneration <- SCMAES(nObjective, parent_list,problem,successProbTarget,crossoverProbability,TRUE)
      else
        newGeneration <- CMAES(nObjective, parent_list,problem,successProbTarget,crossoverProbability,TRUE)


      parent_list <- newGeneration[[1]]
      population <- newGeneration[[2]]
    }else{
      newGeneration <- SMSEMOA(nObjective, population,problem,crossoverProbability,mutationProbability,TRUE,mutationDistribution,crossoverDistribution)
      population <- newGeneration$population
    }


    #evaluate the new generation
    for(popIndex in 1:populationSize){
      individualObjectiveValue<-EvaluateIndividual(population[,popIndex],nObjective)
      populationObjective<-cbind(populationObjective,individualObjectiveValue)
    }
    # scale the objective (in case of WFG problem)
    if((problemType == 3) && (WFGScaling == TRUE))
      scaledObj <- populationObjective/(seq(2,2*nObjective,2))
    else
      scaledObj <- populationObjective

    #compute the hypervolume
    if(max(scaledObj)<2) # if the farthest point is smaller than 2, no shift needed
      Smetric <- GetHypervolume(scaledObj,reference = rep(2,nObjective),method=getOption('hypervolumeMethod'))
    else{ # if any objective exceed 2, shift it to a number slightly less than 2 so it will give nearly zero hypervolume
      scaledObjMaxed <- (scaledObj>2)*1.999999+(scaledObj<2)*scaledObj
      Smetric <- GetHypervolume(scaledObjMaxed,reference = rep(2,nObjective),method=getOption('hypervolumeMethod'))
    }
  }
  return(list(population,Smetric))
}
