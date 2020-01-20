#' etree: Energy Trees for Mixed Type Data
#'
#' The etree package allows to build regression/classification conditional
#' inference trees. The splitting procedure is carried out by means of energy
#' tests, thus allowing for covariates of mixed type. Current implementation
#' covers the case of functions, graphs and persistence diagrams, as well as
#' standard categorical and numerical values, as covariates.
#'
#'
#' @import grid
#' @import partykit
#' @import graphics
#' @importFrom stats model.frame terms approxfun density ecdf getCall dist var
#' model.response model.weights predict quantile weighted.mean knots na.exclude
#' @importFrom fda.usc optim.basis kmeans.fd metric.lp
#' @importFrom utils capture.output combn
#' @importFrom igraph sample_gnp ecount as_adjacency_matrix coreness
#' @importFrom NetworkDistance nd.csd
#' @importFrom grDevices gray.colors
#'
#' @docType package
#' @name etree

NULL
