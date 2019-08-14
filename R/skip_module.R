skip_if_missing_module <- function() {
  have_numpy <- reticulate::py_module_available("numpy")
  if (!have_numpy)
    testthat::skip("numpy not available for testing")

  have_pygmo <- reticulate::py_module_available("pygmo")
  if (!have_pygmo)
    testthat::skip("pygmo not available for testing")

  if(!have_numpy || !have_pygmo)
    testthat::skip("Missing required python modules, skipping test.")
}
