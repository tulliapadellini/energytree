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
#' add_numbers(1, 2) ## returns 3
#'

predict.party <- function(object, newdata = NULL, perm = NULL, ...)
{

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
#' add_numbers(1, 2) ## returns 3
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

