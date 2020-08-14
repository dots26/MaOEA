#' Create a list with cmaes_gen class. Basically, the function transform the population into a class that is accepted by the MOCMAES and SMOCMAES function.
#' @title Generator for cmaes_gen class.
#' @param population The number of objective functions. A scalar value.
#' @param ps_target The target success rate. Used to initialize cmaes_gen$averageSuccessRate.
#' @param stepSize The initial step size.
#' @param evoPath A vector of numbers indicating evolution path of each variable.
#' @param covarianceMatrix Covariance matrix of the variables.
#' @return An object of cmaes_gen class. It can be used as MO-CMA-ES parent. It is a 5 tuple: x (the design point, length = number of variable),averageSuccessRate (scalar),stepSize (scalar), evoPath (evolution path, vector, length = number of variable ),covarianceMatrix (square matrix with ncol = nrow = number of variable).
#' @examples
#' nVar <- 14
#' nObjective <- 5
#' nIndividual <- 100
#' crossoverProbability <- 1
#' ps_target <- 1 / (5 + ( 1 / 2  )^0.5 )
#' pop <- matrix(stats::runif(nIndividual*nVar), nrow = nVar) # create the population
#' a_list <- cmaes_gen(pop)
#' control <- list(successProbTarget=ps_target,crossoverProbability=crossoverProbability)
#'  \donttest{
#' # run a generation of MO-CMA-ES with standard WFG8 test function.
#' numpyready <- reticulate::py_module_available('numpy')
#' pygmoready <- reticulate::py_module_available('pygmo')
#' py_module_ready <- numpyready && pygmoready
#' if(py_module_ready) # prevent error on testing the example
#' newGeneration <- MOCMAES(a_list,nObjective,WFG8,control,nObjective)
#' }
#' @export
cmaes_gen <- function(population,ps_target= (1 / (5 + ( 1 / 2  )^0.5)),stepSize=0.5,evoPath = rep(0,nrow(population)),
                      covarianceMatrix = diag(nrow(population)) ){
  a_list <- vector('list') # container for the parent list
  populationSize <- ncol(population)
  for (populationIndex in 1:populationSize) {
    a <- list(x = population[,populationIndex],
              averageSuccessRate = ps_target,
              stepSize = stepSize,
              evoPath = evoPath,
              covarianceMatrix = covarianceMatrix)
    a_list[[populationIndex]] <- a
  }
  class(a_list) <- 'cmaes_gen'
  return(a_list)
}
