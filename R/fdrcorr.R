#' fdrcorr
## this function is used to calculate for FDR based on the list of P val.
#' @export
fdrcorr <- function(test_list) {
  # Use the p.adjust function from the stats package in R to adjust the p-values
  # in the test_list using the false discovery rate method.
  p.adjust(test_list, method = 'fdr')
}
#' err_handle
## this function will be call with the libary
#' @export
err_handle<-function(x){ tryCatch(x, error=function(e){NA}) }
