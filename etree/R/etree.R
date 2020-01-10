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
#' @importFrom stats model.frame terms approxfun density ecdf getCall model.response model.weights predict quantile weighted.mean
#'             dist var
#' @importFrom fda.usc optim.basis kmeans.fd metric.lp
#' @importFrom utils capture.output combn
#' @importFrom grid depth
#' @importFrom cluster pam
NULL
