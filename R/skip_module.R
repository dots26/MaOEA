skip_if_missing_module <- function() {
  pkg.globals$have_numpy <<- reticulate::py_module_available("numpy")
  pkg.globals$have_pygmo <<- reticulate::py_module_available("pygmo")
  if (!pkg.globals$have_numpy)
    testthat::skip("numpy not available for testing")

  if (!pkg.globals$have_pygmo)
    testthat::skip("pygmo not available for testing")

  if(!pkg.globals$have_numpy || !pkg.globals$have_pygmo)
    testthat::skip("Missing required python modules, skipping test.")
}
