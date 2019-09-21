#' Main interface for the many-objective optimization evolutionary algorithm (MaOEA) package.
#'
#' @title Elitist Non-dominated Sorting Genetic Algorithm version III
#' @param x The initial population. If not supplied, will be generated using LHS. Column major, each column contain one entry.
#' @param fun Objective function being solved.
#' @param nObjective The number of objective functions. A scalar value.
#' @param solver Function name of the solver. Currently available: SMSEMOA, MOCMAES, SMOCMAES, and NSGA3.
#' @param nGeneration Optional, the number of generation the solver should run.
#' @param nVar Number of variables, will be used if \code{x} is not given.
#' @param populationSize Number of individuals in the population, will be used if \code{x} is not given.
#' @param seed random number seed for reproduction of code
#' @param control A list, containing the following:
#' \code{weightVectorSet} A set of weight vector for the optimizer. The weight vector can be any point in the objective space. If not supplied, 5*nObjective points are generated from a sobol sequence. Size: nrow = nObjective,ncol = number of weight vectors
#' \code{crossoverProbability} The probability of doing crossover. Should be between 0-1. Negative value will behave like a zero, and values larger than 1 will behave like 1. Default to 1.
#' \code{mutationProbability} The probability of doing mutation. Should be between 0-1. Negative value will behave like a zero, and values larger than 1 will behave like 1. Default to 1
#' \code{WFGScaling} The use of scaling factor in WFG. Will be ignored in DTLZ problems. Without the scaling, the Pareto front would be on the all-positive portion of hypersphere with radius 1.
#' \code{mutationDistribution} The distribution index for polynomial mutation. Larger index makes the distribution sharper around the parent.
#' \code{crossoverDistribution} The distribution index for SBX. Larger index makes the distribution sharper around each parent.
#' @param ... Further arguments to be passed to \code{fun}
#' @return Returns a list for the next generation
#' \code{population} The new generation design points.
#' \code{populationObjective} The new generation's objective values.
#' @examples
#' \donttest{
#' nVar <- 14
#' nObjective <- 5
#' nIndividual <- 100
#' #control for NSGA3
#' ctrl <- list(crossoverProbability = 1,
#'              mutationProbability = 1/nVar)
#' #Initial population can be supplied, like below but for this example, we skip it
#' #population <- matrix(runif(nIndividual*nVar), nrow = nVar)
#'
#' # Hybrid NSGA-III and SMSEMOA example
#' # 2 calls for nObjective. 1 for optimMaOEA, 1 for WFG8
#' # generate initial population and run 10 gen. NSGA-III with standard WFG8 test function.
#' newPop <- optimMaOEA( , WFG8,NSGA3,nObjective,10,nVar,nIndividual,,ctrl,nObjective)$x
#'
#' # run 5 generations of SMSEMOA with standard WFG8 test function starting with newPop.
#' result <- optimMaOEA( newPop, WFG8,SMSEMOA,nObjective,5,,,1000,ctrl,nObjective)
#' finalPop <- result$x
#' finalObjective <- result$y
#'
#' }
#' @export
optimMaOEA <- function(x=NULL,
                       fun,
                       solver=NSGA3,
                       nObjective,
                       nGeneration = 1,
                       nVar=nrow(x),
                       populationSize=ncol(x),
                       seed =2000,
                       control=list(),...){
  set.seed(seed)
  #initialize population if not specified
  if(is.null(x)){
    if(is.null(nVar)) stop('Initial population and number of variable is not specified. Need either one to run.')
    x<-InitializePopulationLHS(populationSize,nVar,FALSE,0,1)
  }
  if(identical(fun,DTLZ4))
    message('DTLZ4 may suffer from floating-point inaccuracy.')
  con <- list( hypervolumeReferencePoint=NULL,
               hypervolumeTarget=0.99,
               weightVectorSet=NULL,
               crossoverProbability=1,
               crossoverDistribution=30,
               mutationProbability=1,
               mutationDistribution=20,
               successProbThreshold = 0.44,
               successProbTarget = 0.5)

  con[names(control)] <- control

  for(i in 1:nGeneration){
    newGeneration <- solver(x,fun,nObjective,con,...)
    x <- newGeneration$population
    y <- newGeneration$objective
  }

  return(list(x=x,y=y))
}

pkg.globals <- new.env()

pkg.globals$pygmo <- NULL
pkg.globals$rndGen <- NULL
pkg.globals$have_numpy <- F
pkg.globals$have_pygmo <- F

.onLoad <- function(libname, pkgname){
  pkg.globals$have_numpy <- reticulate::py_module_available("numpy")
  pkg.globals$have_pygmo <- reticulate::py_module_available("pygmo")

  if(pkg.globals$have_pygmo && pkg.globals$have_numpy){
    pkg.globals$pygmo <- reticulate::import("pygmo", delay_load = TRUE)
    pkg.globals$rndGen <- reticulate::import("numpy", delay_load = TRUE)
  }
}

setLoadAction(function(ns){
  #have_numpy <<- reticulate::py_module_available("numpy")
  if (!pkg.globals$have_numpy)
    packageStartupMessage("Numpy not available")

  #have_pygmo <<- reticulate::py_module_available("pygmo")
  if (!pkg.globals$have_pygmo)
    packageStartupMessage("PyGMO not available")

  if(!pkg.globals$have_numpy || !pkg.globals$have_pygmo)
    packageStartupMessage("Missing required python modules.
Try using MaOEA::install_python_dependencies()
or follow the instructions in https://esa.github.io/pagmo2/install.html
and call MaOEA::load_python_dependencies().")
})

#' Install the required python package via conda.
#' @title Install python modules required by MaOEA: numpy and PyGMO
#' @param conda Default: auto
#' @param envname Python virtual environment where the modules will be installed, default to 'r-reticulate'
#' @param ... Further argument to pass to reticulate::py_install
#' @return 0 if dependencies installed and loaded successfully, 1 if fails.
#' @export
install_python_dependencies <- function(conda = "auto",envname=NULL,...) {
  if(is.null(envname)){
    envname <- 'r-reticulate'
  }

  will_install <- utils::askYesNo(paste0('This will install numpy and pygmo and all their dependencies to the \'', envname,'\' environment.
           Do you wish to continue?'),default=F,prompts = gettext(c("Y","N","Cancel")))

  if(will_install){
    if (!pkg.globals$have_numpy)
      reticulate::py_install("numpy", method = 'conda', conda = conda,envname = envname)

    if (!pkg.globals$have_pygmo)
      reticulate::py_install("pygmo", method = 'conda', conda = conda,envname = envname)
  }
  success <- load_python_dependencies()
  return(success)
}

#' Import the required python package if it fails onLoad.
#' @title Install python modules required by MaOEA: numpy and PyGMO
#' @return 0 if dependencies loaded successfully, 1 if fails.
#' @export
load_python_dependencies <- function(){
  pkg.globals$have_numpy <- reticulate::py_module_available("numpy")
  pkg.globals$have_pygmo <- reticulate::py_module_available("pygmo")

  if(pkg.globals$have_pygmo && pkg.globals$have_numpy){
    pkg.globals$pygmo <- reticulate::import("pygmo", delay_load = TRUE)
    pkg.globals$rndGen <- reticulate::import("numpy", delay_load = TRUE)

    return(0)
  }else{
    return(1)
  }
}
