#' Energy Tree
#'
#' Fits an energy tree for classification/regression using mixed type data.
#'
#' @param response response variable (either numeric or factor).
#' @param covariates covariates. Must be provided as a list, where each element of the list is a different variable.
#' @param case.weights an optional numeric vector of weights to be used in the fitting process.
#' @param minbucket minimum number of observations that each terminal node must contain. Default is 1.
#' @param alpha significance level for the global test of association and, if \code{split.type = "coeff"} and \code{coef.split.type = "test"}, for the test used in each split. Default is 0.05.
#' @param R number of replicates for the global test of association and, if \code{split.type = "coeff"} and \code{coef.split.type = "test"}, for the test used in each split. Deafult is 1000.
#' @param split.type type of the split when covariates are "complex" (i.e. they are not numeric or factor). This variable can be se to either "coeff" (only available when there is a coeffient representation) or "cluster"
#' @param coef.split.type either "variance" or "test"
#' @param nb number of basis to use for fdata covariates if \code{split.type = "coeff"}.
#'
#' @export
#'
#' @examples
#'  ## returns 3
#'

etree <- function(response, covariates, case.weights = NULL, minbucket = 1, alpha = 0.05, R = 1000, split.type = 'coeff', coef.split.type = 'test', nb = 5) {

  # Check whether covariates is a list
  if(!is.list(covariates)) stop("Argument 'covariates' must be provided as a list")

  # Number of covariates
  n.var = length(covariates)

  # If the case weights are not provided, they are all initialized as 1
  if(is.null(case.weights))
    case.weights <- rep(1L, as.numeric(length(response)))

  # New list of covariates (needed here to build the df used by party)
  newcovariates = lapply(covariates, function(j){
    if(class(j) == 'fdata'){

      if(split.type == "coeff"){

        foo <- fda.usc::optim.basis(j, numbasis = nb)
        fd3 <- fda.usc::fdata2fd(foo$fdata.est,
                                 type.basis = "bspline",
                                 nbasis = foo$numbasis.opt)
        foo <- t(fd3$coefs)

      } else if(split.type == "cluster"){

        foo <- as.factor(1:length(response))

      }

      return(foo)

    } else if(class(j) == 'list' &
              all(sapply(j, class) == 'igraph')){

      if(split.type == "coeff"){
        foo <- graph.shell(j)
      } else if(split.type == "cluster"){
        foo <- as.factor(1:length(response))
      }

      return(foo)

    } else {

      return(j)

    }
  }
  )

  names(newcovariates) <- 1:length(newcovariates)

  # Distances
  cov.distance <- lapply(covariates, compute.dissimilarity)

  # Large list with covariates, newcovariates and distances
  covariates.large = list('cov' = covariates, 'newcov' = newcovariates, 'dist' = cov.distance)

  # Growing the tree (finds the split rules)
  nodes <- growtree(id = 1L,
                    response = response,
                    covariates = covariates.large,
                    case.weights = case.weights,
                    minbucket = minbucket,
                    alpha = alpha,
                    R = R,
                    n.var = n.var,
                    split.type = split.type,
                    coef.split.type = coef.split.type,
                    nb = nb)
  print(c('NODES', nodes))

  # Actually performing the splits
  fitted.obs <- fitted_node(nodes, data = newcovariates)

  # Returning a rich constparty object
  ret <- party(nodes,
               data = newcovariates,
               fitted = data.frame("(fitted)" = fitted.obs,
                                   "(response)" = response,
                                   check.names = FALSE),
               terms = terms(response ~ ., data = newcovariates))

  return(etree = as.constparty(ret))

}



#' Energy Tree Predictions
#'
#' Compute predictions based on an Energy Tree Fit.
#'
#' @param object object of class party.
#' @param newdata an optional list of variables used to make predictions. If omitted, the fitted values are used.
#' @param perm an optional character vector of variable names. Splits of nodes with a primary split in any of these variables will be permuted (after dealing with surrogates). Note that surrogate split in the \code{perm} variables will no be permuted.
#' @param ... additional arguments.
#'
#' @export
#'
#' @examples
#'  ## returns 3
#'


predict.party <- function(object, newdata = NULL, split.type, nb, perm = NULL, ...)
{

  if(!is.null(newdata)){

    newdata = lapply(newdata, function(j){
      if(class(j) == 'fdata' && split.type == "coeff"){

        foo <- fda.usc::optim.basis(j, numbasis = nb)
        fd3 <- fda.usc::fdata2fd(foo$fdata.est,
                                 type.basis = "bspline",
                                 nbasis = foo$numbasis.opt)
        foo <- t(fd3$coefs)
        return(foo)

      } else if(class(j) == 'list' &
                all(sapply(j, class) == 'igraph') & split.type == "coeff"){
        foo <- graph.shell(j)
        return(foo)
      } else {

        return(j)

      }
    }
    )
  }

  ### compute fitted node ids first
  fitted <- if(is.null(newdata) && is.null(perm)) {
    object$fitted[["(fitted)"]]
  } else {
    if (is.null(newdata)) newdata <- model.frame(object)
    ### make sure all the elements in newdata have the same number of rows
    stopifnot(length(unique(sapply(newdata, NROW))) == 1L)

    terminal <- nodeids(object, terminal = TRUE)

    if(max(terminal) == 1L) {
      rep.int(1L, unique(sapply(newdata, NROW)))
    } else {

      inner <- 1L:max(terminal)
      inner <- inner[-terminal]

      primary_vars <- nodeapply(object, ids = inner, by_node = TRUE, FUN = function(node) {
        varid_split(split_node(node))
      })
      surrogate_vars <- nodeapply(object, ids = inner, by_node = TRUE, FUN = function(node) {
        surr <- surrogates_node(node)
        if(is.null(surr)) return(NULL) else return(sapply(surr, varid_split))
      })
      vnames <- names(object$data)

      ### the splits of nodes with a primary split in perm
      ### will be permuted
      if (!is.null(perm)) {
        if (is.character(perm)) {
          stopifnot(all(perm %in% vnames))
          perm <- match(perm, vnames)
        } else {
          ### perm is a named list of factors coding strata
          ### (for varimp(..., conditional = TRUE)
          stopifnot(all(names(perm) %in% vnames))
          stopifnot(all(sapply(perm, is.factor)))
          tmp <- vector(mode = "list", length = length(vnames))
          tmp[match(names(perm), vnames)] <- perm
          perm <- tmp
        }
      }

      ## ## FIXME: the is.na() call takes loooong on large data sets
      ## unames <- if(any(sapply(newdata, is.na)))
      ##     vnames[unique(unlist(c(primary_vars, surrogate_vars)))]
      ## else
      ##     vnames[unique(unlist(primary_vars))]
      unames <- vnames[unique(unlist(c(primary_vars, surrogate_vars)))]

      vclass <- structure(lapply(object$data, class), .Names = vnames)
      ndnames <- names(newdata)
      ndclass <- structure(lapply(newdata, class), .Names = ndnames)
      checkclass <- all(sapply(unames, function(x)
        isTRUE(all.equal(vclass[[x]], ndclass[[x]]))))
      factors <- sapply(unames, function(x) inherits(object$data[[x]], "factor"))
      checkfactors <- all(sapply(unames[factors], function(x)
        isTRUE(all.equal(levels(object$data[[x]]), levels(newdata[[x]])))))
      ## FIXME: inform about wrong classes / factor levels?
      if(all(unames %in% ndnames) && checkclass && checkfactors) {
        vmatch <- match(vnames, ndnames)
        fitted_node_predict(node_party(object), data = newdata,
                            vmatch = vmatch, perm = perm)
      } else {
        if (!is.null(object$terms)) {
          ### <FIXME> this won't work for multivariate responses
          ### </FIXME>
          xlev <- lapply(unames[factors],
                         function(x) levels(object$data[[x]]))
          names(xlev) <- unames[factors]
          #         mf <- model.frame(delete.response(object$terms), newdata,
          #                          xlev = xlev)
          # fitted_node_predict(node_party(object), data = newdata,
          #             vmatch = match(vnames, names(mf)), perm = perm)
          fitted_node_predict(node_party(object), data = newdata,
                              perm = perm)
        } else
          stop("") ## FIXME: write error message
      }
    }
  }
  ### compute predictions
  predict_party(object, fitted, newdata, ...)
}

## Versione Riccardo
# predict.party <- function(object, newdata = NULL, perm = NULL, ...)
# {
#
#   ### compute fitted node ids first
#   fitted <- if(is.null(newdata) && is.null(perm)) {
#     object$fitted[["(fitted)"]]
#   } else {
#     if (is.null(newdata)) newdata <- model.frame(object)
#     ### make sure all the elements in newdata have the same number of rows
#     stopifnot(length(unique(sapply(newdata, NROW))) == 1L)
#
#     terminal <- nodeids(object, terminal = TRUE)
#
#     if(max(terminal) == 1L) {
#       rep.int(1L, unique(sapply(newdata, NROW)))
#     } else {
#
#       inner <- 1L:max(terminal)
#       inner <- inner[-terminal]
#
#       primary_vars <- nodeapply(object, ids = inner, by_node = TRUE, FUN = function(node) {
#         varid_split(split_node(node))
#       })
#       surrogate_vars <- nodeapply(object, ids = inner, by_node = TRUE, FUN = function(node) {
#         surr <- surrogates_node(node)
#         if(is.null(surr)) return(NULL) else return(sapply(surr, varid_split))
#       })
#       vnames <- names(object$data)
#
#       ### the splits of nodes with a primary split in perm
#       ### will be permuted
#       if (!is.null(perm)) {
#         if (is.character(perm)) {
#           stopifnot(all(perm %in% vnames))
#           perm <- match(perm, vnames)
#         } else {
#           ### perm is a named list of factors coding strata
#           ### (for varimp(..., conditional = TRUE)
#           stopifnot(all(names(perm) %in% vnames))
#           stopifnot(all(sapply(perm, is.factor)))
#           tmp <- vector(mode = "list", length = length(vnames))
#           tmp[match(names(perm), vnames)] <- perm
#           perm <- tmp
#         }
#       }
#
#       ## ## FIXME: the is.na() call takes loooong on large data sets
#       ## unames <- if(any(sapply(newdata, is.na)))
#       ##     vnames[unique(unlist(c(primary_vars, surrogate_vars)))]
#       ## else
#       ##     vnames[unique(unlist(primary_vars))]
#       unames <- vnames[unique(unlist(c(primary_vars, surrogate_vars)))]
#
#       vclass <- structure(lapply(object$data, class), .Names = vnames)
#       ndnames <- names(newdata)
#       ndclass <- structure(lapply(newdata, class), .Names = ndnames)
#       checkclass <- all(sapply(unames, function(x)
#         isTRUE(all.equal(vclass[[x]], ndclass[[x]]))))
#       factors <- sapply(unames, function(x) inherits(object$data[[x]], "factor"))
#       checkfactors <- all(sapply(unames[factors], function(x)
#         isTRUE(all.equal(levels(object$data[[x]]), levels(newdata[[x]])))))
#       ## FIXME: inform about wrong classes / factor levels?
#       if(all(unames %in% ndnames) && checkclass && checkfactors) {
#         vmatch <- match(vnames, ndnames)
#         fitted_node_predict(node_party(object), data = newdata,
#                             vmatch = vmatch, perm = perm)
#       } else {
#         if (!is.null(object$terms)) {
#           ### <FIXME> this won't work for multivariate responses
#           ### </FIXME>
#           xlev <- lapply(unames[factors],
#                          function(x) levels(object$data[[x]]))
#           names(xlev) <- unames[factors]
#           #         mf <- model.frame(delete.response(object$terms), newdata,
#           #                          xlev = xlev)
#           # fitted_node_predict(node_party(object), data = newdata,
#           #             vmatch = match(vnames, names(mf)), perm = perm)
#           fitted_node_predict(node_party(object), data = newdata,
#                               perm = perm)
#         } else
#           stop("") ## FIXME: write error message
#       }
#     }
#   }
#   ### compute predictions
#   predict_party(object, fitted, newdata, ...)
# }



#' Visualization of Energy Trees
#'
#' \code{plot} method for \code{party} objects with extended facilities for plugging in panel functions.
#'
#' @param x	an object of class party or constparty.
#' @param main an optional title for the plot.
#' @param type a character specifying the complexity of the plot: extended tries to visualize the distribution of the response variable in each terminal node whereas simple only gives some summary information.
#' @param terminal_panel an optional panel function of the form function(node) plotting the terminal nodes. Alternatively, a panel generating function of class "grapcon_generator" that is called with arguments x and tp_args to set up a panel function. By default, an appropriate panel function is chosen depending on the scale of the dependent variable.
#' @param tp_args	a list of arguments passed to terminal_panel if this is a "grapcon_generator" object.
#' @param inner_panel	an optional panel function of the form function(node) plotting the inner nodes. Alternatively, a panel generating function of class "grapcon_generator" that is called with arguments x and ip_args to set up a panel function.
#' @param ip_args	a list of arguments passed to inner_panel if this is a "grapcon_generator" object.
#' @param edge_panel an optional panel function of the form function(split, ordered = FALSE, left = TRUE) plotting the edges. Alternatively, a panel generating function of class "grapcon_generator" that is called with arguments x and ip_args to set up a panel function.
#' @param ep_args	a list of arguments passed to edge_panel if this is a "grapcon_generator" object.
#' @param drop_terminal	a logical indicating whether all terminal nodes should be plotted at the bottom.
#' @param tnex a numeric value giving the terminal node extension in relation to the inner nodes.
#' @param newpage	a logical indicating whether grid.newpage() should be called.
#' @param pop	a logical whether the viewport tree should be popped before return.
#' @param gp graphical parameters.
#' @param margins	numeric vector of margin sizes.
#' @param digits number of digits to be printed.
#' @param ... additional arguments.
#'
#' @export
#'
#' @examples
#' ## returns 3
#'

plot.party <- function(x, main = NULL,
                       terminal_panel = node_terminal, tp_args = list(),
                       inner_panel = node_inner, ip_args = list(),
                       edge_panel = edge_simple, ep_args = list(),
                       drop_terminal = FALSE, tnex = 1,
                       newpage = TRUE, pop = TRUE, gp = gpar(),
                       margins = NULL, ...)
{

  ### extract tree
  node <- node_party(x)
  ### total number of terminal nodes
  nx <- width(node)
  ### maximal depth of the tree
  ny <- depth(node, root = TRUE)

  ## setup newpage
  if (newpage) grid.newpage()

  ## setup root viewport
  margins <- if(is.null(margins)) {
    c(1, 1, if(is.null(main)) 0 else 3, 1)
  } else {
    rep_len(margins, 4L)
  }
  root_vp <- viewport(layout = grid.layout(3, 3,
                                           heights = unit(c(margins[3L], 1, margins[1L]),
                                                          c("lines", "null", "lines")),
                                           widths = unit(c(margins[2L], 1, margins[4L]),
                                                         c("lines", "null", "lines"))),
                      name = "root",
                      gp = gp)
  pushViewport(root_vp)

  ## viewport for main title (if any)
  if (!is.null(main)) {
    main_vp <- viewport(layout.pos.col = 2, layout.pos.row = 1,
                        name = "main")
    pushViewport(main_vp)
    grid.text(y = unit(1, "lines"), main, just = "center")
    upViewport()
  }

  ## setup viewport for tree
  tree_vp <- viewport(layout.pos.col = 2, layout.pos.row = 2,
                      xscale = c(0, nx), yscale = c(0, ny + (tnex - 1)),
                      name = "tree")
  pushViewport(tree_vp)

  ### setup panel functions (if necessary)
  if(inherits(terminal_panel, "grapcon_generator"))
    terminal_panel <- do.call("terminal_panel", c(list(x), as.list(tp_args)))
  if(inherits(inner_panel, "grapcon_generator"))
    inner_panel <- do.call("inner_panel", c(list(x), as.list(ip_args)))
  if(inherits(edge_panel, "grapcon_generator"))
    edge_panel <- do.call("edge_panel", c(list(x), as.list(ep_args)))


  if((nx <= 1 & ny <= 1)) {
    if(is.null(margins)) margins <- rep.int(1.5, 4)
    pushViewport(plotViewport(margins = margins, name = paste("Node", id_node(node), sep = "")))
    terminal_panel(node)
  } else {
    ## call the workhorse
    .plot_node(node,
               xlim = c(0, nx), ylim = c(0, ny - 0.5 + (tnex - 1)),
               nx = nx, ny = ny,
               terminal_panel = terminal_panel,
               inner_panel = inner_panel,
               edge_panel = edge_panel,
               tnex = tnex,
               drop_terminal = drop_terminal,
               debug = FALSE)
  }
  upViewport()
  if (pop) popViewport() else upViewport()
}



growtree <- function(id = 1L,
                     response,
                     covariates,
                     case.weights,
                     minbucket,
                     alpha,
                     R,
                     n.var,
                     split.type = 'coeff',
                     coef.split.type = 'test',
                     nb) {

  # For less than <minbucket> observations, stop here
  if (sum(case.weights) < minbucket)
    return(partynode(id = id))

  # Finding the best split (variable selection & split point search)
  split <- findsplit(response = response,
                     covariates = covariates,
                     alpha = alpha,
                     R = R,
                     lp = rep(2, 2),
                     split.type = split.type,
                     coef.split.type = coef.split.type,
                     nb = nb)

  # If no split is found, stop here
  if (is.null(split))
    return(partynode(id = id))

  # Selected variable index and possibly selected basis index
  varid <- split$varid
  if(!is.null(split$basid)){
    basid <- split$basid
  }

  breaks <- split$breaks
  index <- split$index

  # Assigning the ids to the observations
  kidids <- c()
  switch(class(covariates$cov[[varid]]),

         fdata = {

           if(split.type == 'coeff'){

             # observations before the split point are assigned to node 1
             kidids[which(covariates$newcov[[varid]][, basid] <= breaks)] <- 1
             #  observations before the split point are assigned to node 2
             kidids[which(covariates$newcov[[varid]][, basid] > breaks)] <- 2

           } else if (split.type == 'cluster') {

             kidids <- na.exclude(index)

           }
         },

         numeric = {

           kidids[(which(covariates$cov[[varid]] <= breaks))] <- 1
           kidids[(which(covariates$cov[[varid]] > breaks))] <- 2

         },

         integer = {

           kidids[(which(covariates$newcov[[varid]] <= breaks))] <- 1
           kidids[(which(covariates$newcov[[varid]] > breaks))] <- 2

         },

         factor = {

           kidids <- na.exclude(index)

         },

         list = if(FALSE){
           #attributes(x[[1]])$names == 'diagram'
         } else if(all(sapply(covariates$cov[[varid]], class) == 'igraph')){

           if(split.type == 'coeff'){

             kidids[which(covariates$newcov[[varid]][, basid] <= breaks)] <- 1
             kidids[which(covariates$newcov[[varid]][, basid] > breaks)] <- 2

           } else if(split.type == 'cluster') {

             kidids <- na.exclude(index)

           }
         }
  )

  # If all the observations belong to the same node, no split is done
  if (all(kidids == 1) | all(kidids == 2))
    return(partynode(id = id))

  # Initialization of the kid nodes
  kids <- vector(mode = "list", length = max(kidids, na.rm = TRUE))

  # Giving birth to the kid nodes
  for (kidid in 1:length(kids)) {
    # selecting observations for the current node
    w <- case.weights
    w[kidids != kidid] <- 0

    # getting next node id
    if (kidid > 1) {
      myid <- max(nodeids(kids[[kidid - 1]]))
    } else{
      myid <- id
    }

    # starting recursion on this kid node
    covariates.updated <- list()
    covariates.updated$cov <- lapply(covariates$cov, function(cov) subset(cov, as.logical(w)))
    covariates.updated$newcov <- lapply(covariates$newcov, function(cov) subset(cov, as.logical(w)))
    covariates.updated$dist <- lapply(covariates$dist, function(cov) subset(cov, subset = as.logical(w), select = which(w == 1)))

    kids[[kidid]] <-
      growtree(
        id = as.integer(myid + 1),
        response = subset(response, as.logical(w)),
        covariates = covariates.updated,
        case.weights = rep(1L, sum(w, na.rm = TRUE)),
        minbucket,
        alpha,
        R,
        n.var = n.var,
        split.type = split.type,
        coef.split.type = coef.split.type,
        nb = nb)
  }

  # Return the nodes (i.e. the split rules)
  return(partynode(id = as.integer(id),
                   split = split,
                   kids = kids,
                   info = list(p.value = min(info_split(split)$p.value, na.rm = TRUE))
  ))
}



# Find split --------------------------------------------------------------

findsplit <- function(response,
                      covariates,
                      alpha,
                      R,
                      lp = rep(2,2),
                      split.type = 'coeff',
                      coef.split.type = 'test',
                      nb) {

  # Number of original covariates
  n.cov = length(covariates$cov)

  print('one round again')
  # Performing an independence test between the response and each covariate
  p = lapply(covariates$dist,
             function(cov.dist) {
               #set.seed(12345)
               ct <- energy::dcor.test(cov.dist, compute.dissimilarity(response), R = R)
               if (!is.na(ct$statistic)) {
                 return(c(ct$statistic, ct$p.value))
               } else{
                 c(NA, NA)
               }
             }
  )

  p = t(matrix(unlist(p), ncol = 2, byrow = T))
  rownames(p) <- c("statistic", "p-value")
  if (all(is.na(p[2,]))) return(NULL)

  # Bonferroni correction
  minp <- min(p[2,], na.rm = TRUE)
  minp <- 1 - (1 - minp) ^ sum(!is.na(p[2,]))
  if (minp > alpha) return(NULL)

  # Variable selection
  if (length(which(p[2,] == min(p[2,], na.rm = T))) > 1) {
    xselect <- which.max(p[1,])    # in case of multiple minima, take that with the highest test statistic
  } else{
    xselect <- which.min(p[2,])
  }

  # Selected covariates
  x <- covariates$cov[[xselect]]
  newx <- covariates$newcov[[xselect]]
  if(split.type == 'cluster'){
    xdist <- covariates$dist[[xselect]]
  }

  # Split point search
  split.objs = split.opt(y = response,
                         x = x,
                         newx = newx,
                         xdist = xdist,
                         split.type = split.type,
                         coef.split.type = coef.split.type,
                         nb = nb)

  # Separately saving split.objs outputs
  splitindex <- split.objs$splitindex
  bselect <- split.objs$bselect
  centroids <- split.objs$centroids

  # Returning the split point
  switch(class(x),

         numeric = {

           return(sp = partysplit(varid = as.integer(xselect),
                                  breaks = splitindex,
                                  info = list(p.value = 1-(1-p)^sum(!is.na(p)))))

         },

         integer = {

           return(sp = partysplit(varid = as.integer(xselect),
                                  breaks = splitindex,
                                  info = list(p.value = 1-(1-p)^sum(!is.na(p)))))

         },

         factor = {

           return(sp = partysplit(varid = as.integer(xselect),
                                  index = splitindex,
                                  info = list(p.value = 1-(1-p)^sum(!is.na(p)))))

         },

         fdata = {

           if(split.type == 'coeff'){
             return(sp = partysplit(varid = as.integer(xselect),
                                    basid = as.integer(bselect),
                                    breaks = splitindex,
                                    info = list(p.value = 1-(1-p[2,])^sum(!is.na(p[2,])))))
           } else if(split.type == 'cluster'){
             return(sp = partysplit(varid = as.integer(xselect),
                                    centroids = centroids,
                                    index = as.integer(splitindex),
                                    info = list(p.value = 1-(1-p[2,])^sum(!is.na(p[2,])))))
           }

         },

         list = if(FALSE){
           #attributes(v[[1]])$names == 'diagram'
           return(sp = partysplit(varid = as.integer(xselect),
                                  index = as.integer(splitindex),
                                  info = list(p.value = 1-(1-p[2,])^sum(!is.na(p[2,])))))

         } else if(all(sapply(x, class) == 'igraph')){

           if(split.type == 'coeff'){

             return(sp = partysplit(varid = as.integer(xselect),
                                    basid = as.integer(bselect),
                                    breaks = splitindex,
                                    info = list(p.value = 1-(1-p[2,])^sum(!is.na(p[2,])))))

           } else if(split.type == 'cluster') {

             return(sp = partysplit(varid = as.integer(xselect),
                                    centroids = centroids,
                                    index = as.integer(splitindex),
                                    info = list(p.value = 1-(1-p[2,])^sum(!is.na(p[2,])))))

           }
         }

  )
}



# Split point search ------------------------------------------------------

#' Find Split Value
#'
#' Computes optimal split value
#'
#' @param y response variable
#' @param x selected covariate
#'
#' @export
#'
#' @examples
#' add_numbers(1, 2) ## returns 3
#'

split.opt <- function(y,
                      x,
                      newx,
                      xdist,
                      split.type = 'coeff',
                      coef.split.type = 'test',
                      nb,
                      R=1000,
                      wass.dist = NULL){

  switch(class(x),

         factor     = {

           lev <- levels(x[drop = TRUE])
           if (length(lev) == 2) {
             splitpoint <- lev[1]
           } else{
             comb <- do.call("c", lapply(1:(length(lev) - 1),
                                         ### TBC: isn't this just floor(length(lev)/2) ??
                                         function(x)combn(lev,
                                                          x,
                                                          simplify = FALSE)))
             xlogp <- sapply(comb, function(q) mychisqtest(x %in% q, y))
             splitpoint <- comb[[which.min(xlogp)]]
           }

           # split into two groups (setting groups that do not occur to NA)
           splitindex <- !(levels(x) %in% splitpoint)
           splitindex[!(levels(x) %in% lev)] <- NA_integer_
           splitindex <- splitindex - min(splitindex, na.rm = TRUE) + 1L

         },

         numeric    = {

           s  <- sort(x)
           comb = sapply(s[2:(length(s)-1)], function(j) x<j)
           #first and last one are excluded (trivial partitions)
           xp.value <- apply(comb, 2, function(q) independence.test(x = q, y = y))
           if (length(which(xp.value[2,] == min(xp.value[2,], na.rm = T))) > 1) {
             splitindex <- s[which.max(xp.value[1,])]
           } else {
             splitindex <- s[which.min(xp.value[2,])]
           }

         },

         integer    = {

           s  <- sort(x)
           comb = sapply(s[2:(length(s)-1)], function(j) x<j)
           xp.value <- apply(comb, 2, function(q) independence.test(x = q, y = y))
           if (length(which(xp.value[2,] == min(xp.value[2,], na.rm = T))) > 1) {
             splitindex <- s[which.max(xp.value[1,])]
           } else {
             splitindex <- s[which.min(xp.value[2,])]
           }

         },

         fdata      = {

           if(split.type == 'coeff'){
             x1 = newx
             bselect <- 1:dim(x1)[2]
             p1 <- c()
             p1 <- sapply(bselect, function(i) independence.test(x1[, i], y, R = R))
             colnames(p1) <- colnames(x1)
             if (length(which(p1[2,] == min(p1[2,], na.rm = T))) > 1) {
               bselect <- as.integer(which.max(p1[1,]))
             } else{
               bselect <- as.integer(which.min(p1[2,]))
             }
             sel.coeff = x1[,bselect]
             s  <- sort(sel.coeff)
             comb = sapply(s[1:(length(s)-1)], function(j) sel.coeff<=j)

             if(coef.split.type == 'variance'){

               obj <- apply(comb, 2, function(c){
                 data1 <- y[c]
                 data2 <- y[!c]
                 v1 <- var(data1)
                 v2 <- var(data2)
                 n1 <- length(data1)
                 n2 <- length(data2)
                 n <- n1+n2
                 obj_c <- (n1*v1+n2*v2)/n
                 return(obj_c)})
               splitindex <- s[which.min(obj)]

             } else if (coef.split.type == 'test'){

               xp.value <- apply(comb, 2, function(q) independence.test(x = q, y = y))
               if (length(which(xp.value[2,] == min(xp.value[2,], na.rm = T))) > 1) {
                 splitindex <- s[which.max(xp.value[1,])]
               } else {
                 splitindex <- s[which.min(xp.value[2,])]
               }

             }

           } else if(split.type == 'cluster') {
             cl.fdata = kmeans.fd(x, ncl=2, draw = FALSE, par.ini=list(method="exact"), cluster.size = 1)
             clindex <- cl.fdata$cluster
             lev = levels(newx)
             splitindex = rep(NA, length(lev))
             splitindex[lev %in% newx[clindex==1]]<- 1
             splitindex[lev %in% newx[clindex==2]]<- 2

             c1 <- cl.fdata$centers[1]
             c2 <- cl.fdata$centers[2]
             centroids <- list(c1 = c1, c2 = c2)
           }

         },

         list = if(FALSE){
           #attributes(x[[1]])$names == "diagram"
           cl.diagrams = cluster::pam(wass.dist, k = 2, diss = TRUE)
           splitindex <- cl.diagrams$clustering

         } else if(all(sapply(x, class) == 'igraph')){

           if(split.type == 'coeff'){
             x1 = newx
             bselect <- 1:dim(x1)[2]
             p1 <- c()
             p1 <- sapply(bselect, function(i) independence.test(x1[, i], y, R = R))
             colnames(p1) <- colnames(x1)
             if (length(which(p1[2,] == min(p1[2,], na.rm = T))) > 1) {
               bselect <- as.integer(which.max(p1[1,]))
             } else{
               bselect <- as.integer(which.min(p1[2,]))
             }
             sel.coeff = x1[,bselect]
             s  <- sort(sel.coeff)
             comb = sapply(s[1:(length(s)-1)], function(j) sel.coeff<=j)

             if(coef.split.type == 'variance'){

               obj <- apply(comb, 2, function(c){
                 data1 <- y[c]
                 data2 <- y[!c]
                 v1 <- var(data1)
                 v2 <- var(data2)
                 n1 <- length(data1)
                 n2 <- length(data2)
                 n <- n1+n2
                 obj_c <- (n1*v1+n2*v2)/n
                 return(obj_c)})
               splitindex <- s[which.min(obj)]

             } else if (coef.split.type == 'test'){

               xp.value <- apply(comb, 2, function(q) independence.test(x = q, y = y))
               if (length(which(xp.value[2,] == min(xp.value[2,], na.rm = T))) > 1) {
                 splitindex <- s[which.max(xp.value[1,])]
               } else {
                 splitindex <- s[which.min(xp.value[2,])]
               }

             }

           } else if(split.type == 'cluster') {
             cl.graph <- cluster::pam(xdist, k = 2, diss = TRUE)
             clindex <- cl.graph$clustering
             lev = levels(newx)
             splitindex = rep(NA, length(lev))
             splitindex[lev %in% newx[clindex==1]]<- 1
             splitindex[lev %in% newx[clindex==2]]<- 2

             ceindex1 <- as.integer(cl.graph$medoids[1])
             c1 <- x[[which(newx == ceindex1)]]
             ceindex2 <- as.integer(cl.graph$medoids[2])
             c2 <- x[[which(newx == ceindex2)]]
             centroids <- list(c1 = c1, c2 = c2)
             #the which part is necessary since ceindex (pam medoids indices) go from 1 to the TOTAL number of observations
           }


         }
  )

  out <- list('splitindex' = splitindex)
  if(exists('bselect')) out$bselect <- bselect
  if(exists('centroids')) out$centroids <- centroids
  return(out)

}



# Independence (dcor) test ------------------------------------------------

independence.test <- function(x,
                              y,
                              R = 1000,
                              lp = c(2,2)) {

  # Computing the dissimilarities within x and y
  d1 = compute.dissimilarity(x, lp = lp[1])
  d2 = compute.dissimilarity(y, lp = lp[2])

  # Distance correlation test
  #set.seed(12345)
  ct <- energy::dcor.test(d1, d2, R = R)
  if (!is.na(ct$statistic)) {
    return(c(ct$statistic, ct$p.value))
  } else{
    c(NA, NA)
  }
}


# Distances ---------------------------------------------------------------

compute.dissimilarity <- function(x,
                                  lp = 2){

  # Computing the dissimilarities
  switch(class(x),
         logical    = as.matrix(dist(x)),
         factor     = as.matrix(cluster::daisy(as.data.frame(x))),
         numeric    = as.matrix(dist(x)),
         integer    = as.matrix(dist(x)),
         matrix     = as.matrix(dist(x)),
         data.frame = as.matrix(dist(x)),
         fdata      = metric.lp(x, lp=lp),
         list       = {
           if(all(sapply(x, class) == 'igraph')){
             adj_matrices <- lapply(x, as_adjacency_matrix)
             d <- nd.csd(adj_matrices) #continuous spectral density for the moment
             return(as.matrix(d$D))
           } else if(all(sapply(x, function(x) attributes(x)$names) == 'diagram')){
             k.fun = function(i,j) TDA::wasserstein(x[[i]], x[[j]])
             k.fun = Vectorize(k.fun)
             d.idx = seq_along(x)
             outer(d.idx,d.idx, k.fun)
           }
         })

}




compute.dissimilarity.cl <- function(centroid, x,
                                     lp = 2){

  switch(class(x),
         fdata      = metric.lp(fdata1 = x, fdata2 = centroid, lp=lp),
         list       = {
           if(all(sapply(x, class) == 'igraph')){
             adj_matrices <- lapply(x, as_adjacency_matrix)
             d <- nd.csd(adj_matrices) #continuous spectral density for the moment
             return(as.matrix(d$D))
           } else if (all(sapply(x, function(x) attributes(x)$names) == 'diagram')){
             k.fun = function(i, centroid) TDA::wasserstein(x[[i]], centroid)
             k.fun = Vectorize(k.fun)
             d.idx = seq_along(x)
             outer(d.idx, centroid, k.fun)
           }
         })
  # list       = {
  #   if(!is.null(attributes(x[[1]]))){
  #   if(attributes(x[[1]])$names == "diagram"){
  #     d1 = x[case.weights]
  #     k.fun = function(i, j) TDA::wasserstein(d1[[i]], d1[[j]])
  #     k.fun = Vectorize(k.fun)
  #     d.idx = seq_along(d1)
  #     outer(d.idx,d.idx, k.fun)
  #   }}
  #})

}



# Graphs ------------------------------------------------------------------

graph.shell <- function(graph.list, shell.limit = NULL){

  # Number of observations (graphs)
  n.graphs <- length(graph.list)

  # Shell distribution for each graph
  table.shell <- lapply(graph.list, function(g){table(coreness(g))})

  # Maximum shell index
  max.shell <- do.call(max, lapply(table.shell,
                                   function(s){
                                     as.integer(names(s))
                                   }))

  # Column names for the shell df
  col.names = as.character(seq(1, max.shell, 1))
  #starting from 1 since we presumably only deal with connected graphs

  # Shell df inizialization
  all.shell.df = data.frame(matrix(
    data = 0L,
    nrow = n.graphs,
    ncol = length(col.names)))
  colnames(all.shell.df) <- col.names

  # Fill in with the actual shell distibutions
  invisible(sapply(1:n.graphs, function(i){
    cols <- names(table.shell[[i]])
    all.shell.df[i, cols] <<- table.shell[[i]][cols] # <<- for global environment assignment
  }))
  # better a for cycle?
  # for(i in 1:n.graphs){
  #   cols <- names(table.shell[[i]])
  #   all.shell.df[i, cols] = table.shell[[i]][cols]
  # }

  # No more than 'shell.limit' indices for each graph
  if(!is.null(shell.limit) && max.shell > shell.limit){
    all.shell.df <- all.shell.df[,as.character(seq(1, shell.limit, 1))]
  }

  # Return the final shell df
  return(all.shell.df)

}
