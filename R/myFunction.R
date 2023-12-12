#' myFunction
#'
#' @description
#' adds two numbers
#'
#'
#' @param x a numeric vector
#' @param y a numeric vector
#'
#' @return a numeric value
#'
#' @export
#'
#' @examples
#' myFunction(x = 1, y = 2)
myFunction <- function(x, z) {
  sumResult <- x + z

  myHiddenFunction(sumResult)

  return(sumResult)
}

myHiddenFunction <- function(sumResult) {
  message("The sum is: ", sumResult)
}
