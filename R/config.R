#file:config.R
options(hypervolumeMethod = "approx") # no caps: approx or exact.
#' Set hypervolume method
#' @param method The number of individual in the population. Integer > 0.
#'
#' @return Nothing
#' @examples
#' setHypervolumeMethod("approx") # available: exact or approx
#' @export
setHypervolumeMethod<- function(method){
  options(hypervolumeMethod= method)
}

#' Set hypervolume method
#' @param method The number of individual in the population. Integer > 0.
#'
#' @return The hypervolume method set by either setHypervolumeMethod() or options(hypervolumeMethod = ...)
#' @examples
#' setHypervolumeMethod("approx") # this is also the default value
#' getHypervolumeMethod() # returns NSGA
#' @export
getHypervolumeMethod<- function(){
  return( getOptions(hypervolumeMethod))
}
