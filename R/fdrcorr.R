## this function is used to calculate for FDR based on the list of P val.
fdrcorr <- function(test_list) {
  # Use the p.adjust function from the stats package in R to adjust the p-values
  # in the test_list using the false discovery rate method.
  p.adjust(test_list, method = 'fdr')
}

## this function will be call with the libary
err_handle<-function(x){ tryCatch(x, error=function(e){NA}) }
