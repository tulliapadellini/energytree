
# node --------------------------------------------------------------------

partynode <- function(id, split = NULL, kids = NULL, surrogates = NULL, info = NULL, centroids = NULL) {

  if (!is.integer(id) || length(id) != 1) {
    id <- as.integer(id0 <- id)
    if (any(is.na(id)) || !isTRUE(all.equal(id0, id)) || length(id) != 1)
      stop(sQuote("id"), " ", "must be a single integer")
  }

  if (is.null(split) != is.null(kids)) {
    stop(sQuote("split"), " ", "and", " ", sQuote("kids"), " ",
         "must either both be specified or unspecified")
  }

  if (!is.null(kids)) {
    if (!(is.list(kids) && all(sapply(kids, inherits, "partynode")))
        || length(kids) < 2)
      stop(sQuote("kids"), " ", "must be an integer vector or a list of",
           " ", sQuote("partynode"), " ", "objects")
  }

  if (!is.null(surrogates)) {
    if (!is.list(surrogates) || any(!sapply(surrogates, inherits, "partysplit")))
      stop(sQuote("split"), " ", "is not a list of", " ", sQuote("partysplit"),
           " ", "objects")
  }

  node <- list(id = id, split = split, kids = kids, surrogates = surrogates, info = info, centroids = centroids)
  class(node) <- "partynode"
  return(node)
}

as.partynode.list <- function(x, ...) {

  if (!all(sapply(x, inherits, what = "list")))
    stop("'x' has to be a list of lists")

  if (!all(sapply(x, function(x) "id" %in% names(x))))
    stop("each list in 'x' has to define a node 'id'")

  ok <- sapply(x, function(x)
    all(names(x) %in% c("id", "split", "kids", "surrogates", "info", "centroids")))
  if (any(!ok))
    sapply(which(!ok), function(i)
      warning(paste("list element", i, "defines additional elements:",
                    paste(names(x[[i]])[!(names(x[[i]]) %in%
                                            c("id", "split", "kids", "surrogates", "info", "centroids"))],
                          collapse = ", "))))

  ids <- as.integer(sapply(x, function(node) node$id))
  if(any(duplicated(ids))) stop("nodeids must be unique integers")
  x <- x[order(ids)]
  ids <- ids[order(ids)]

  new_recnode <- function(i) {
    x_i <- x[[which(ids == i)]]
    if (is.null(x_i$kids))
      partynode(id = x_i$id, info = x_i$info)
    else
      partynode(id = x_i$id, split = x_i$split,
                kids = lapply(x_i$kids, new_recnode),
                surrogates = x_i$surrogates,
                info = x_i$info)
  }

  ret <- new_recnode(ids[1L])
  ### <FIXME> duplicates recursion but makes sure
  ###    that the ids are in pre-order notation with
  ###    from defined in as.partynode.partynode
  ### </FIXME>
  as.partynode(ret, ...)
}

kidids_node <- function(node, data, vmatch = 1:length(data), obs = NULL,
                        perm = NULL) {

  primary <- split_node(node)
  surrogates <- surrogates_node(node)

  ### perform primary split
  x <- kidids_split(primary, data, vmatch, obs)

  ### surrogate / random splits if needed
  if (any(is.na(x))) {
    ### surrogate splits
    if (length(surrogates) >= 1) {
      for (surr in surrogates) {
        nax <- is.na(x)
        if (!any(nax)) break;
        x[nax] <- kidids_split(surr, data, vmatch, obs = obs)[nax]
      }
    }
    nax <- is.na(x)
    ### random splits
    if (any(nax)) {
      prob <- prob_split(primary)
      x[nax] <- sample(1:length(prob), sum(nax), prob = prob,
                       replace = TRUE)
    }
  }

  ### permute variable `perm' _after_ dealing with surrogates etc.
  if (!is.null(perm)) {
    if (is.integer(perm)) {
      if (varid_split(primary) %in% perm)
        x <- .resample(x)
    } else {
      if (is.null(obs)) obs <- 1:length(data)
      strata <- perm[[varid_split(primary)]]
      if (!is.null(strata)) {
        strata <- strata[obs, drop = TRUE]
        for (s in levels(strata))
          x[strata == s] <- .resample(x[strata == s])
      }
    }
  }
  return(x)
}

kidids_node_predict <- function(node, data, vmatch = 1:length(data), obs = NULL,
                                perm = NULL) {

  primary <- split_node(node)
  surrogates <- surrogates_node(node)

  ### perform primary split
  x <- kidids_split_predict(primary, data, vmatch, obs)

  ### surrogate / random splits if needed
  if (any(is.na(x))) {
    ### surrogate splits
    if (length(surrogates) >= 1) {
      for (surr in surrogates) {
        nax <- is.na(x)
        if (!any(nax)) break;
        x[nax] <- kidids_split_predict(surr, data, vmatch, obs = obs)[nax]
      }
    }
    nax <- is.na(x)
    ### random splits
    if (any(nax)) {
      prob <- prob_split(primary)
      x[nax] <- sample(1:length(prob), sum(nax), prob = prob,
                       replace = TRUE)
    }
  }

  ### permute variable `perm' _after_ dealing with surrogates etc.
  if (!is.null(perm)) {
    if (is.integer(perm)) {
      if (varid_split(primary) %in% perm)
        x <- .resample(x)
    } else {
      if (is.null(obs)) obs <- 1:length(data)
      strata <- perm[[varid_split(primary)]]
      if (!is.null(strata)) {
        strata <- strata[obs, drop = TRUE]
        for (s in levels(strata))
          x[strata == s] <- .resample(x[strata == s])
      }
    }
  }
  return(x)
}

fitted_node <- function(node, data, vmatch = 1:length(data),
                        obs = 1:unique(sapply(data, NROW)), perm = NULL) {

  if (is.logical(obs)) obs <- which(obs)
  if (is.terminal(node))
    return(rep(id_node(node), length(obs)))
  retid <- nextid <- kidids_node(node, data, vmatch, obs, perm)

  for (i in unique(nextid)) {
    indx <- nextid == i
    retid[indx] <- fitted_node(kids_node(node)[[i]], data,
                               vmatch, obs[indx], perm)
  }
  return(retid)
}

fitted_node_predict <- function(node, data, vmatch = 1:length(data),
                                obs = 1:unique(sapply(data, NROW)), perm = NULL) {

  if (is.logical(obs)) obs <- which(obs)
  if (is.terminal(node))
    return(rep(id_node(node), length(obs)))
  retid <- nextid <- kidids_node_predict(node, data, vmatch, obs, perm)

  for (i in unique(nextid)) {
    indx <- nextid == i
    retid[indx] <- fitted_node_predict(kids_node(node)[[i]], data,
                                       vmatch, obs[indx], perm)
  }
  return(retid)
}


# split -------------------------------------------------------------------



partysplit <- function(varid, breaks = NULL, index = NULL, right = TRUE,
                       prob = NULL, info = NULL, centroids = NULL, basid = NULL) {

  ### informal class for splits
  split <- vector(mode = "list", length = 8)
  names(split) <- c("varid", "breaks", "index", "right", "prob", "info", "centroids", "basid")

  ### split is an id referring to a variable
  if (!is.integer(varid))
    stop(sQuote("varid"), " ", "is not integer")
  split$varid <- varid

  if (is.null(breaks) && is.null(index))
    stop("either", " ", sQuote("breaks"), " ", "or", " ",
         sQuote("index"), " ", "must be given")

  ### vec
  if (!is.null(breaks)) {
    if (is.numeric(breaks) && (length(breaks) >= 1)) {
      ### FIXME: I think we need to make sure breaks are double in C
      split$breaks <- as.double(breaks)
    } else {
      stop(sQuote("break"), " ",
           "should be a numeric vector containing at least one element")
    }
  }

  if (!is.null(index)) {
    if (is.integer(index)) {
      if (!(length(index) >= 2))
        stop(sQuote("index"), " ", "has less than two elements")
      if (!(min(index, na.rm = TRUE) == 1))
        stop("minimum of", " ", sQuote("index"), " ", "is not equal to 1")
      if (!all.equal(diff(sort(unique(index))), rep(1, max(index, na.rm = TRUE) - 1)))
        stop(sQuote("index"), " ", "is not a contiguous sequence")
      split$index <- index
    } else {
      stop(sQuote("index"), " ", "is not a class", " ", sQuote("integer"))
    }
    if (!is.null(breaks)) {
      if (length(breaks) != (length(index) - 1))
        stop("length of", " ", sQuote("breaks"), " ",
             "does not match length of", " ", sQuote("index"))
    }
  }

  if (is.logical(right) & !is.na(right))
    split$right <- right
  else
    stop(sQuote("right"), " ", "is not a logical")

  if (!is.null(prob)) {
    if (!is.double(prob) ||
        (any(prob < 0) | any(prob > 1) | !isTRUE(all.equal(sum(prob), 1))))
      stop(sQuote("prob"), " ", "is not a vector of probabilities")
    if (!is.null(index)) {
      if (!(max(index, na.rm = TRUE) == length(prob)))
        stop("incorrect", " ", sQuote("index"))
    }
    if (!is.null(breaks) && is.null(index)) {
      if (!(length(breaks) == (length(prob) - 1)))
        stop("incorrect", " ", sQuote("breaks"))
    }
    split$prob <- prob
  }

  if (!is.null(info))
    split$info <- info

  if (!is.null(centroids))
    split$centroids <- centroids

  if (!is.null(basid))
    if (!is.integer(varid)){
      stop(sQuote("varid"), " ", "is not integer")
    } else {
      split$basid <- basid
    }

  class(split) <- "partysplit"

  return(split)
}

centroids_split <- function(split) {
  if (!(inherits(split, "partysplit")))
    stop(sQuote("split"), " ", "is not an object of class",
         " ", sQuote("partysplit"))
  split$centroids
}

basid_split <- function(split) {
  if (!(inherits(split, "partysplit")))
    stop(sQuote("split"), " ", "is not an object of class",
         " ", sQuote("partysplit"))
  split$basid
}

kidids_split <- function(split, data, vmatch = 1:length(data), obs = NULL) {

  varid <- varid_split(split)
  basid <- basid_split(split)
  class(data) <- "list" ### speed up
  if(!is.null(basid)){ #means we are in the coeff case
    x <- data[[vmatch[varid]]][,basid]
  } else { #means we are in the cluster case
    x <- data[[vmatch[varid]]]
  }
  if (!is.null(obs)) x <- x[obs]

  if (is.null(breaks_split(split))) {
    if (storage.mode(x) != "integer")
      stop("variable", " ", vmatch[varid], " ", "is not integer")
  } else {
    ### labels = FALSE returns integers and is faster
    ### <FIXME> use findInterval instead of cut?
    #        x <- cut.default(as.numeric(x), labels = FALSE,
    #                 breaks = unique(c(-Inf, breaks_split(split), Inf)),  ### breaks_split(split) = Inf possible (MIA)
    #                 right = right_split(split))
    x <- .bincode(as.numeric(x), # labels = FALSE,
                  breaks = unique(c(-Inf, breaks_split(split), Inf)),  ### breaks_split(split) = Inf possible (MIA)
                  right = right_split(split))
    ### </FIXME>
  }
  index <- index_split(split)
  ### empty factor levels correspond to NA and return NA here
  ### and thus the corresponding observations will be treated
  ### as missing values (surrogate or random splits):
  if (!is.null(index))
    x <- index[x]
  return(x)
}

kidids_split_predict <- function(split, data, vmatch = 1:length(data), obs = NULL) {

  varid <- varid_split(split)
  basid <- basid_split(split)
  class(data) <- "list" ### speed up
  if(!is.null(basid)){ #means we are in the coeff case
    x <- data[[vmatch[varid]]][,basid]

  } else { #means we are in the cluster case
    x <- data[[vmatch[varid]]]
  }
  if (!is.null(obs)) x <- x[obs]

  if (is.null(breaks_split(split))) {

    if(!is.null(centroids_split(split))){

      cl.idx = lapply(centroids_split(split), compute.dissimilarity.cl, x = x, lp = 2)
      x <- apply(matrix(unlist(cl.idx),ncol = 2),1, which.min)
    }

    if (storage.mode(x) != "integer")
      stop("variable", " ", vmatch[varid], " ", "is not integer")
  } else {
    ### labels = FALSE returns integers and is faster
    ### <FIXME> use findInterval instead of cut?
    #        x <- cut.default(as.numeric(x), labels = FALSE,
    #                 breaks = unique(c(-Inf, breaks_split(split), Inf)),  ### breaks_split(split) = Inf possible (MIA)
    #                 right = right_split(split))
    x <- .bincode(as.numeric(x), # labels = FALSE,
                  breaks = unique(c(-Inf, breaks_split(split), Inf)),  ### breaks_split(split) = Inf possible (MIA)
                  right = right_split(split))
    ### </FIXME>
  }
  index <- index_split(split)
  ### empty factor levels correspond to NA and return NA here
  ### and thus the corresponding observations will be treated
  ### as missing values (surrogate or random splits):
  if (!is.null(index) & is.null(centroids_split(split)))
    x <- index[x]
  return(x)
}


# party -------------------------------------------------------------------

## FIXME: data in party
##   - currently assumed to be a data.frame
##   - potentially empty
##   - the following are all assumed to work:
##     dim(data), names(data)
##     sapply(data, class), lapply(data, levels)
##   - potentially these need to be modified if data/terms
##     should be able to deal with data bases

party <- function(node, data, fitted = NULL, terms = NULL, names = NULL, info = NULL) {

  stopifnot(inherits(node, "partynode"))
  #stopifnot(inherits(data, "list"))   #give rise to problems for classif plots
  ### make sure all split variables are there
  ids <- nodeids(node)[!nodeids(node) %in% nodeids(node, terminal = TRUE)]
  varids <- unique(unlist(nodeapply(node, ids = ids, FUN = function(x)
    varid_split(split_node(x)))))
  #stopifnot(varids %in% 1:length(data))

  if(!is.null(fitted)) {
    stopifnot(inherits(fitted, "data.frame"))
    stopifnot(all(sapply(data, NROW) == 0L) | all(sapply(data, NROW) == NROW(fitted)))

    # try to provide default variable "(fitted)"
    if(all(sapply(data, NROW) > 0L)) {
      if(!("(fitted)" %in% names(fitted)))
        fitted[["(fitted)"]] <- fitted_node(node, data = data)
    } else {
      stopifnot("(fitted)" == names(fitted)[1L])
    }

    nt <- nodeids(node, terminal = TRUE)
    stopifnot(all(fitted[["(fitted)"]] %in% nt))

    node <- as.partynode(node, from = 1L)
    nt2 <- nodeids(node, terminal = TRUE)
    fitted[["(fitted)"]] <- nt2[match(fitted[["(fitted)"]], nt)]
  } else {
    node <- as.partynode(node, from = 1L)
    # default "(fitted)"
    if(all(sapply(data, NROW) > 0L) & missing(fitted))
      fitted <- data.frame("(fitted)" = fitted_node(node,
                                                    data = data), check.names = FALSE)
  }

  party <- list(node = node, data = data, fitted = fitted,
                terms = NULL, names = NULL, info = info)
  class(party) <- "party"

  if(!is.null(terms)) {
    stopifnot(inherits(terms, "terms"))
    party$terms <- terms
  }

  if (!is.null(names)) {
    n <- length(nodeids(party, terminal = FALSE))
    if (length(names) != n)
      stop("invalid", " ", sQuote("names"), " ", "argument")
    party$names <- names
  }

  party
}

predict.party <- function(object, newdata = NULL, nb = 10, perm = NULL, ...)
{

  split.type <- det_split.type(object)

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

### do nothing expect returning the fitted ids
predict_party.default <- function(party, id, newdata = NULL, FUN = NULL, ...) {

  if (length(list(...)) > 1)
    warning("argument(s)", " ", sQuote(names(list(...))), " ", "have been ignored")

  ## get observation names: either node names or
  ## observation names from newdata
  nam <- if(is.null(newdata)) {
    if(is.null(rnam <- rownames(data_party(party)))) names(party)[id] else rnam
  } else {
    rownames(newdata[[1]])
  }
  if(length(nam) != length(id)) nam <- NULL

  if (!is.null(FUN))
    return(.simplify_pred(nodeapply(party,
                                    nodeids(party, terminal = TRUE), FUN, by_node = TRUE), id, nam))

  ## special case: fitted ids
  return(structure(id, .Names = nam))
}

predict_party.constparty <- function(party, id, newdata = NULL,
                                     type = c("response", "prob", "quantile", "density", "node"),
                                     at = if (type == "quantile") c(0.1, 0.5, 0.9),
                                     FUN = NULL, simplify = TRUE, ...)
{
  ## extract fitted information
  response <- party$fitted[["(response)"]]
  weights <- party$fitted[["(weights)"]]
  fitted <- party$fitted[["(fitted)"]]
  if (is.null(weights)) weights <- rep(1, NROW(response))

  ## get observation names: either node names or
  ## observation names from newdata
  nam <- if(is.null(newdata)) names(party)[id] else rownames(newdata[[1]])
  if(length(nam) != length(id)) nam <- NULL

  ## match type
  type <- match.arg(type)

  ## special case: fitted ids
  if(type == "node")
    return(structure(id, .Names = nam))

  ### multivariate response
  if (is.data.frame(response)) {
    ret <- lapply(response, function(r) {
      ret <- .predict_party_constparty(node_party(party), fitted = fitted,
                                       response = r, weights, id = id, type = type, at = at, FUN = FUN, ...)
      if (simplify) .simplify_pred(ret, id, nam) else ret
    })
    if (all(sapply(ret, is.atomic)))
      ret <- as.data.frame(ret)
    names(ret) <- colnames(response)
    return(ret)
  }

  ### univariate response
  ret <- .predict_party_constparty(node_party(party), fitted = fitted, response = response,
                                   weights = weights, id = id, type = type, at = at, FUN = FUN, ...)
  if (simplify) .simplify_pred(ret, id, nam) else ret[as.character(id)]
}

data_party.default <- function(party, id = 1L) {

  extract <- function(id) {
    if(is.null(party$fitted))
      if(length(party$data) == 0) return(NULL)
    else
      stop("cannot subset data without fitted ids")

    ### which terminal nodes follow node number id?
    nt <- nodeids(party, id, terminal = TRUE)
    wi <- party$fitted[["(fitted)"]] %in% nt

    ret <- if (length(party$data) == 0)
      subset(party$fitted, wi)
    else
      subset(cbind(party$data, party$fitted), wi)
    ret
  }
  if (length(id) > 1)
    return(lapply(id, extract))
  else
    return(extract(id))
}


# plot --------------------------------------------------------------------

edge_simple <- function(obj, digits = 3, abbreviate = FALSE,
                        justmin = Inf, just = c("alternate", "increasing", "decreasing", "equal"),
                        fill = "white")
{
  meta <- obj$data

  split.type <- det_split.type(obj)

  justfun <- function(i, split) {
    myjust <- if(mean(nchar(split)) > justmin) {
      match.arg(just, c("alternate", "increasing", "decreasing", "equal"))
    } else {
      "equal"
    }
    k <- length(split)
    rval <- switch(myjust,
                   "equal" = rep.int(0, k),
                   "alternate" = rep(c(0.5, -0.5), length.out = k),
                   "increasing" = seq(from = -k/2, to =  k/2, by = 1),
                   "decreasing" = seq(from =  k/2, to = -k/2, by = -1)
    )
    unit(0.5, "npc") + unit(rval[i], "lines")
  }

  ### panel function for simple edge labelling
  function(node, i) {
    split <- character_split(split_node(node), meta, digits = digits)$levels
    y <- justfun(i, split)
    split <- split[i]
    # try() because the following won't work for split = "< 10 Euro", for example.
    if(any(grep(">", split) > 0) | any(grep("<", split) > 0)) {
      tr <- suppressWarnings(try(parse(text = paste("phantom(0)", split)), silent = TRUE))
      if(!inherits(tr, "try-error")) split <- tr
    }
    if (split.type == 'coeff'){
      grid.rect(y = y, gp = gpar(fill = fill, col = 0), width = unit(1, "strwidth", split))
      grid.text(split, y = y, just = "center")
    } else {
      # the number of obs in each kid node is calculated as the number of commas
      # appearing in split (which is a string where the levels are separated by
      # commas), plus one
      n_kid <- as.character(lengths(regmatches(split, gregexpr(",", split))) + 1)
      n_kid <- paste('n =', n_kid)
      grid.rect(y = y, gp = gpar(fill = fill, col = 0), width = unit(1, "strwidth", n_kid))
      grid.text(n_kid, y = y, just = "center")
    }
  }
}
class(edge_simple) <- "grapcon_generator"
