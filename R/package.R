# Package Description for Roxygen:
#' MaOEA contains several algorithms for solving many-objective optimization problems.
#' The algorithms are provided as a sequence of operators used in a single iteration.
#' For example, the SMSEMOA function calls the recombination (SBX) and mutation operator (polynomial mutation) to produce 1 offspring, and perform the S-metric selection.
#' The function then returns a list containing the population and population objective after the procedure is conducted once.
#' The purpose of only doing a single iteration is to support users if they wish to formulate hybrid algorithms.
#'
#' Alternatively, users can use the optimMaOEA function to solve an optimization problem with their chosen algorithm.
#' This function is a simple wrapper to call the algorithms listed above for several iterations.
#' Using this function, users can simply supply the initial population, objective function, the chosen algorithm, and the number of iterations. If number of iteration is not supplied, then only a single iteration is conducted.
#'
#' Note: This package uses column-major ordering, i.e. an individual should be contained in a single column, each row represents different variable.
#' All optimization variable should be scaled to 0-1.
#' \tabular{ll}{
#' Package: \tab MaOEA\cr
#' Type: \tab Package\cr
#' Version: \tab 0.4.1\cr
#' Date: \tab 2019-07-12\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' @name MaOEA-package
#' @aliases MaOEA
#' @docType package
#' @title Many-Objective Evolutionary Algorithm
#' @author Dani Irawan \email{irawan_dani@@yahoo.com}
#' @keywords package
#' @seealso Main interface function is \code{\link{optimMaOEA}}.
#' @import reticulate
#' @import nsga2R
#' @import lhs
#' @import nnet
#' @import stringr
#' @import randtoolbox
#' @import e1071
#'
#' @section Acknowledgments:
#' This work is funded by the European Commission's H2020 programme through
#' the UTOPIAE  Marie  Curie  Innovative  Training Network, H2020-MSCA-ITN-2016,
#' under Grant Agreement No. 722734 as well as  through the Twinning project SYNERGY
#' under Grant Agreement No. 692286.
#'
#' @section Maintainer:
#' Dani Irawan \email{irawan_dani@@yahoo.com}
#End of Package Description
NA
