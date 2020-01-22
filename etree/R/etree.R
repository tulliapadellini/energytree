#' etree: Energy Trees for Mixed Type Data
#'
#' The etree package allows to build regression/classification conditional
#' inference trees. The splitting procedure is carried out by means of energy
#' tests, thus allowing for covariates of mixed type. Current implementation
#' covers the case of functions, graphs and persistence diagrams, as well as
#' standard categorical and numerical values, as covariates.
#'
#'
#' @rawNamespace import(partykit, except = c(partynode, as.partynode.list,
#'   kidids_node, fitted_node, party, predict.party, predict_party.default,
#'   predict_party.constparty, data_party.default, edge_simple, partysplit,
#'   kidids_split))
#' @import graphics
#' @import grid
#' @importFrom stats model.frame terms approxfun density ecdf getCall dist var
#' model.response model.weights predict quantile weighted.mean knots na.exclude
#' @importFrom utils capture.output
#' @importFrom grDevices gray.colors
#' @importFrom survival survfit
#'
#' @docType package
#' @name etree

NULL
