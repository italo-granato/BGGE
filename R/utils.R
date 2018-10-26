#' Print BGGE information object
#'
#' @param x BGGE object
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Displays the most relevant model fit information.
#' @export
#'
print.BGGE <- function(x, ...){
  cat('Model Fitted with: \n',
      x$ite, ' Iterations, burning the first ', x$burn, ' and thining every ', x$thin, '\n\n',
      'Some predicted Values: \n')
  
  print.default(format(head(x$yHat, 10), digits = 3), print.gap = 2L, quote = FALSE)
  
  cat('\n Use str() function to found more datailed information.\n')
  invisible(x)
}

#' @title Comparative plot
#'
#' @description Simple plot of the predicted values versus observed values
#' @usage plot(BGGE_Object, ...)
#' 
#' @param x \code{BGGE object}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @importFrom graphics plot abline
#' @export
plot.BGGE <- function(x, ...){
  ### Check that object is compatible
  if (!inherits(x, "BGGE")) stop("This function only works for objects of class 'BGGE'")
  response <- x$y
  predictions <- x$yHat
  limits <- range(c(response, predictions), na.rm = TRUE)
  plot(response, predictions, xlim = limits, ylim = limits, xlab = 'Observed values', ylab = 'Predicted values', ...);
  abline(a = 0, b = 1, lty = 3)
}
