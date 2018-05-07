#' An example graphNEL
#'
#' Toy example of graphNEL object for pathway the survival test
#'
#' @format A graphNEL object of a pathway with 308 nodes.
#'  
"graph"

#' Expression dataset
#'
#' Toy example of an expression dataset for survival test
#'
#' @format A matrix with 246 genes (rows) measured across 73 patients (columns)
#'  
"exp"

#' Anotations
#'
#' Survival annotation for toy example dataset
#'
#' #' @format A data frame with 73 observations on the following 2 variables.
#' \describe{
#'   \item{\code{days}}{a character vector, with the days to Death/last follow up after 1st surgical event.}
#'   \item{\code{status}}{dead = 1, alive = 0}
#' }
#'
#'  
"survAnnot"
