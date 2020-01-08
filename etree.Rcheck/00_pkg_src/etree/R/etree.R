#' etree: Energy Trees for Mixed Type Data
#'
#' The etree package allows to build regression/classification conditional inference trees.
#' The splitting procedure is carried out by means of energy tests, thus allowing for covariates of mixed type.
#' Current implementation covers the case of functions, graphs and persistence diagrams,
#' as well as standard categorical and numerical values, as covariates.
#'
#'
#' @docType package
#' @name etree
#' @importFrom methods as new
#' @importFrom stats delete.response fitted model.matrix na.omit napredict pnorm predict terms
NULL
