#Custom Error handling: The err_handle() function is a simple function that tries to evaluate an expression x, and if an error occurs it returns NA instead. This is known as error handling, and it allows the code to continue running even if an error is encountered.
err_handle<-function(x){ tryCatch(x, error=function(e){NA}) }
