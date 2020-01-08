
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

is.partynode <- function(x) {
  if (!inherits(x, "partynode")) return(FALSE)
  rval <- diff(nodeids(x, terminal = FALSE))
  isTRUE(all.equal(unique(rval), 1))
}

as.partynode <- function(x, ...)
  UseMethod("as.partynode")

as.partynode.partynode <- function(x, from = NULL, recursive = TRUE, ...) {
  if(is.null(from)) from <- id_node(x)
  from <- as.integer(from)
  if (!recursive) {
    if(is.partynode(x) &
       id_node(x) == from) return(x)
  }
  id <- from - 1L

  new_node <- function(x) {
    id <<- id + 1L
    if(is.terminal(x)) return(partynode(id, info = info_node(x)))
    partynode(id,
              split = split_node(x),
              kids = lapply(kids_node(x), new_node),
              surrogates = surrogates_node(x),
              info = info_node(x))
  }

  return(new_node(x))
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

as.list.partynode <- function(x, ...)
{
  ids <- nodeids(x)
  obj <- vector(mode = "list", length = length(ids))
  thisnode <- NULL

  nodelist <- function(node) {
    if (is.terminal(node))
      obj[[which(ids == id_node(node))]] <<- list(id = id_node(node), info = info_node(node))
    else {
      thisnode <<- list(id = id_node(node), split = split_node(node),
                        kids = sapply(kids_node(node), function(k) id_node(k)))
      if (!is.null(surrogates_node(node)))
        thisnode$surrogates <- surrogates_node(node)
      if (!is.null(info_node(node)))
        thisnode$info <- info_node(node)
      obj[[which(ids == id_node(node))]] <<- thisnode
      lapply(kids_node(node), nodelist)
    }
  }
  nodelist(x)
  return(obj)
}


id_node <- function(node) {
  if (!(inherits(node, "partynode")))
    stop(sQuote("node"), " ", "is not an object of class", " ", sQuote("node"))
  node$id
}

kids_node <- function(node) {
  if (!(inherits(node, "partynode")))
    stop(sQuote("node"), " ", "is not an object of class", " ", sQuote("node"))
  node$kids
}

info_node <- function(node) {
  if (!(inherits(node, "partynode")))
    stop(sQuote("node"), " ", "is not an object of class", " ", sQuote("node"))
  node$info
}

formatinfo_node <- function(node, FUN = NULL, default = "", prefix = NULL, ...) {
  info <- info_node(node)

  ## FIXME: better dispatch to workhorse FUN probably needed in the future, e.g.:
  ## (1) formatinfo() generic with formatinfo.default() as below,
  ## (2) supply default FUN from party$info$formatinfo() or similar.
  if(is.null(FUN)) FUN <- function(x, ...) {
    if(is.null(x)) x <- ""
    if(!is.object(x) & is.atomic(x)) x <- as.character(x)
    if(!is.character(x)) x <- capture.output(print(x), ...)
    x
  }

  info <- if(is.null(info)) default else FUN(info, ...)
  if(!is.null(prefix)) {
    info <- if(length(info) > 1L) c(prefix, info) else paste(prefix, info, sep = "")
  }
  info
}

### FIXME: permutation and surrogate splits: is only the primary
### variable permuted?
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


length.partynode <- function(x)
  length(kids_node(x))

"[.partynode" <- "[[.partynode" <- function(x, i, ...) {
  if (!(length(i) == 1 && is.numeric(i)))
    stop(sQuote("x"), " ", "is incorrect node")
  kids_node(x)[[i]]
}

split_node <- function(node) {
  if (!(inherits(node, "partynode")))
    stop(sQuote("node"), " ", "is not an object of class", " ", sQuote("node"))
  node$split
}

surrogates_node <- function(node) {
  if (!(inherits(node, "partynode")))
    stop(sQuote("node"), " ", "is not an object of class", " ", sQuote("node"))
  node$surrogates
}

is.terminal <- function(x, ...)
  UseMethod("is.terminal")

is.terminal.partynode <- function(x, ...) {
  kids <- is.null(kids_node(x))
  split <- is.null(split_node(x))
  if (kids != split)
    stop("x", " ", "is incorrect node")
  kids
}

## ## depth generic now taken from package 'grid'
## depth <- function(x, ...)
##     UseMethod("depth")

depth.partynode <- function(x, root = FALSE, ...) {
  if (is.terminal(x)) return(as.integer(root))
  max(sapply(kids_node(x), depth, root = root)) + 1L
}

width <- function(x, ...)
  UseMethod("width")

width.partynode <- function(x, ...) {
  if (is.terminal(x)) return(1)
  sum(sapply(kids_node(x), width.partynode))
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

varid_split <- function(split) {
  if (!(inherits(split, "partysplit")))
    stop(sQuote("split"), " ", "is not an object of class",
         " ", sQuote("partysplit"))
  split$varid
}

breaks_split <- function(split) {
  if (!(inherits(split, "partysplit")))
    stop(sQuote("split"), " ", "is not an object of class",
         " ", sQuote("partysplit"))
  split$breaks
}

index_split <- function(split) {
  if (!(inherits(split, "partysplit")))
    stop(sQuote("split"), " ", "is not an object of class",
         " ", sQuote("partysplit"))
  split$index
}

right_split <- function(split) {
  if (!(inherits(split, "partysplit")))
    stop(sQuote("split"), " ", "is not an object of class",
         " ", sQuote("partysplit"))
  split$right
}

prob_split <- function(split) {
  if (!(inherits(split, "partysplit")))
    stop(sQuote("split"), " ", "is not an object of class",
         " ", sQuote("partysplit"))
  prob <- split$prob
  if (!is.null(prob)) return(prob)

  ### either breaks or index must be there
  if (is.null(index <- index_split(split))) {
    if (is.null(breaks <- breaks_split(split)))
      stop("neither", " ", sQuote("prob"), " ", "nor", " ",
           sQuote("index"), " ", "or", sQuote("breaks"), " ",
           "given for", " ", sQuote("split"))
    nkids <- length(breaks) + 1
  } else {
    nkids <- max(index, na.rm = TRUE)
  }
  prob <- rep(1, nkids) / nkids
  return(prob)
}

info_split <- function(split) {
  if (!(inherits(split, "partysplit")))
    stop(sQuote("split"), " ", "is not an object of class",
         " ", sQuote("partysplit"))
  split$info
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
  if (!is.null(index))
    x <- index[x]
  return(x)
}

character_split <- function(split, data = NULL, digits = getOption("digits") - 2) {

  varid <- varid_split(split)

  if (!is.null(data)) {
    ## names and labels
    lev <- lapply(data, levels)[[varid]]
    mlab <- names(data)[varid]

    ## determine split type
    type <- sapply(data, function(x) class(x)[1])[varid_split(split)]
    type[!(type %in% c("factor", "ordered"))] <- "numeric"
  } else {
    ## (bad) default names and labels
    lev <- NULL
    mlab <- paste("V", varid, sep = "")
    type <- "numeric"
  }

  ## process defaults for breaks and index
  breaks <- breaks_split(split)
  index <- index_split(split)
  right <- right_split(split)

  if (is.null(breaks)) breaks <- 1:(length(index) - 1)
  if (is.null(index)) index <- 1:(length(breaks) + 1)

  ## check whether ordered are really ordered
  if (type == "ordered") {
    if (length(breaks) > 1)
      type <- "factor"
  }
  ### <FIXME> format ordered multiway splits? </FIXME>

  switch(type,
         "factor" = {
           nindex <- index[cut(seq_along(lev), c(-Inf, breaks, Inf), right = right)]
           dlab <- as.vector(tapply(lev, nindex, paste, collapse = ", "))
         },
         "ordered" = {
           if (length(breaks) == 1) {
             if (right)
               dlab <- paste(c("<=", ">"), lev[breaks], sep = " ")
             else
               dlab <- paste(c("<", ">="), lev[breaks], sep = " ")
           } else {
             stop("") ### see above
           }
           dlab <- dlab[index]
         },
         "numeric" = {
           breaks <- round(breaks, digits)
           if (length(breaks) == 1) {
             if (right)
               dlab <- paste(c("<=", ">"), breaks, sep = " ")
             else
               dlab <- paste(c("<", ">="), breaks, sep = " ")
           } else {
             dlab <- levels(cut(0, breaks = c(-Inf, breaks, Inf),
                                right = right))
           }
           dlab <- as.vector(tapply(dlab, index, paste, collapse = " | "))
         }
  )

  return(list(name = mlab, levels = dlab))
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
  stopifnot(inherits(data, "list"))
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

length.party <- function(x)
  length(nodeids(x))

names.party <- function(x)
  .names_party(x)

"names<-.party" <- function(x, value) {
  n <- length(nodeids(x, terminal = FALSE))
  if (!is.null(value) && length(value) != n)
    stop("invalid", " ", sQuote("names"), " ", "argument")
  x$names <- value
  x
}

.names_party <- function(party) {
  names <- party$names
  if (is.null(names))
    names <- as.character(nodeids(party, terminal = FALSE))
  names
}

node_party <- function(party) {
  stopifnot(inherits(party, "party"))
  party$node
}

is.constparty <- function(party) {
  stopifnot(inherits(party, "party"))
  if (!is.null(party$fitted))
    return(all(c("(fitted)", "(response)") %in% names(party$fitted)))
  return(FALSE)
}

as.constparty <- function(obj, ...) {
  if(!inherits(obj, "party")) obj <- as.party(obj)
  if (!is.constparty(obj)) {
    if(is.null(obj$fitted))
      obj$fitted <- data.frame("(fitted)" = predict(obj, type = "node"), check.names = FALSE)
    if(!("(fitted)" %in% names(obj$fitted)))
      obj$fitted["(fitted)"] <- predict(obj, type = "node")
    if(!("(response)" %in% names(obj$fitted)))
      obj$fitted["(response)"] <- model.response(model.frame(obj))
    if(!("(weights)" %in% names(obj$fitted))) {
      w <- model.weights(model.frame(obj))
      if(is.null(w) && any(w != 1L)) obj$fitted["(weights)"] <- w
    }
  }
  if (is.constparty(obj)) {
    ret <- obj
    class(ret) <- c("constparty", class(obj))
    return(ret)
  }
  stop("cannot coerce object of class", " ", sQuote(class(obj)),
       " ", "to", " ", sQuote("constparty"))
}

"[.party" <- "[[.party" <- function(x, i, ...) {
  if (is.character(i) && !is.null(names(x)))
    i <- which(names(x) %in% i)
  stopifnot(length(i) == 1 & is.numeric(i))
  stopifnot(i <= length(x) & i >= 1)
  i <- as.integer(i)
  dat <- data_party(x, i)
  if (!is.null(x$fitted)) {
    findx <- which("(fitted)" == names(dat))[1]
    fit <- dat[,findx:ncol(dat), drop = FALSE]
    dat <- dat[,-(findx:ncol(dat)), drop = FALSE]
    if (ncol(dat) == 0)
      dat <- x$data
  } else {
    fit <- NULL
    dat <- x$data
  }
  nam <- names(x)[nodeids(x, from = i, terminal = FALSE)]

  recFun <- function(node) {
    if (id_node(node) == i) return(node)
    kid <- sapply(kids_node(node), id_node)
    return(recFun(node[[max(which(kid <= i))]]))
  }
  node <- recFun(node_party(x))

  ret <- party(node = node, data = dat, fitted = fit,
               terms = x$terms, names = nam, info = x$info)
  class(ret) <- class(x)
  ret
}

nodeids <- function(obj, ...)
  UseMethod("nodeids")

nodeids.partynode <- function(obj, from = NULL, terminal = FALSE, ...) {

  if(is.null(from)) from <- id_node(obj)

  id <- function(node, record = TRUE, terminal = FALSE) {
    if(!record) return(NULL)
    if(!terminal)
      return(id_node(node))
    else
      if(is.terminal(node)) return(id_node(node)) else return(NULL)
  }

  rid <- function(node, record = TRUE, terminal = FALSE) {
    myid <- id(node, record = record, terminal = terminal)
    if(is.terminal(node)) return(myid)
    kids <- kids_node(node)
    kids_record <- if(record)
      rep(TRUE, length(kids))
    else
      sapply(kids, id_node) == from
    return(c(myid,
             unlist(lapply(1:length(kids), function(i)
               rid(kids[[i]], record = kids_record[i], terminal = terminal)))
    ))
  }

  return(rid(obj, from == id_node(obj), terminal))
}

nodeids.party <- function(obj, from = NULL, terminal = FALSE, ...)
  nodeids(node_party(obj), from = from, terminal = terminal, ...)

nodeapply <- function(obj, ids = 1, FUN = NULL, ...)
  UseMethod("nodeapply")

nodeapply.party <- function(obj, ids = 1, FUN = NULL, by_node = TRUE, ...) {

  stopifnot(isTRUE(all.equal(ids, round(ids))))
  ids <- as.integer(ids)

  if(is.null(FUN)) FUN <- function(x, ...) x

  if (length(ids) == 0)
    return(NULL)

  if (by_node) {
    rval <- nodeapply(node_party(obj), ids = ids, FUN = FUN, ...)
  } else {
    rval <- lapply(ids, function(i) FUN(obj[[i]], ...))
  }

  names(rval) <- names(obj)[ids]
  return(rval)
}

nodeapply.partynode <- function(obj, ids = 1, FUN = NULL, ...) {

  stopifnot(isTRUE(all.equal(ids, round(ids))))
  ids <- as.integer(ids)

  if(is.null(FUN)) FUN <- function(x, ...) x

  if (length(ids) == 0)
    return(NULL)

  rval <- vector(mode = "list", length = length(ids))
  rval_id <- rep(0, length(ids))
  i <- 1

  recFUN <- function(node, ...) {
    if(id_node(node) %in% ids) {
      rval_id[i] <<- id_node(node)
      rval[[i]] <<- FUN(node, ...)
      i <<- i + 1
    }
    kids <- kids_node(node)
    if(length(kids) > 0) {
      for(j in 1:length(kids)) recFUN(kids[[j]])
    }
    invisible(TRUE)
  }
  foo <- recFUN(obj)
  rval <- rval[match(ids, rval_id)]
  return(rval)
}

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

predict_party <- function(party, id, newdata = NULL, ...)
  UseMethod("predict_party")

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

### functions for node prediction based on fitted / response
.pred_Surv <- function(y, w) {
  if (length(y) == 0) return(NA)
  survfit(y ~ 1, weights = w, subset = w > 0)
}

.pred_Surv_response <- function(y, w) {
  if (length(y) == 0) return(NA)
  .median_survival_time(.pred_Surv(y, w))
}

.pred_factor <- function(y, w) {
  lev <- levels(y)
  sumw <- tapply(w, y, sum)
  sumw[is.na(sumw)] <- 0
  prob <- sumw / sum(w)
  names(prob) <- lev
  return(prob)
}

.pred_factor_response <- function(y, w) {
  prob <- .pred_factor(y, w)
  return(factor(which.max(prob), levels = 1:nlevels(y),
                labels = levels(y),
                ordered = is.ordered(y)))
  return(prob)
}

.pred_numeric_response <- function(y, w)
  weighted.mean(y, w, na.rm = TRUE)

.pred_ecdf <- function(y, w) {
  if (length(y) == 0) return(NA)
  iw <- as.integer(round(w))
  if (max(abs(w - iw)) < sqrt(.Machine$double.eps)) {
    y <- rep(y, w)
    return(ecdf(y))
  } else {
    stop("cannot compute empirical distribution function with non-integer weights")
  }
}

.pred_quantile <- function(y, w) {
  y <- rep(y, w)
  function(p, ...) quantile(y, probs = p, ...)
}

.pred_density <- function(y, w) {
  d <- density(y, weights = w / sum(w))
  approxfun(d[1:2], rule = 2)
}

### workhorse: compute predictions based on fitted / response data
.predict_party_constparty <- function(node, fitted, response, weights,
                                      id = id, type = c("response", "prob", "quantile", "density"),
                                      at = if (type == "quantile") c(0.1, 0.5, 0.9), FUN = NULL, ...) {

  type <- match.arg(type)
  if (is.null(FUN)) {

    rtype <- class(response)[1]
    if (rtype == "ordered") rtype <- "factor"
    if (rtype == "integer") rtype <- "numeric"
    if (rtype == "AsIs") rtype <- "numeric"

    if (type %in% c("quantile", "density") && rtype != "numeric")
      stop("quantile and density estimation currently only implemented for numeric responses")

    FUN <- switch(rtype,
                  "Surv" = if (type == "response") .pred_Surv_response else .pred_Surv,
                  "factor" = if (type == "response") .pred_factor_response else .pred_factor,
                  "numeric" = switch(type,
                                     "response" = .pred_numeric_response,
                                     "prob" = .pred_ecdf,
                                     "quantile" = .pred_quantile,
                                     "density" = .pred_density)
    )
  }

  ## empirical distribution in each leaf
  if (all(id %in% fitted)) {
    tab <- tapply(1:NROW(response), fitted,
                  function(i) FUN(response[i], weights[i]), simplify = FALSE)
  } else {
    ### id may also refer to inner nodes
    tab <- as.array(lapply(sort(unique(id)), function(i) {
      index <- fitted %in% nodeids(node, i, terminal = TRUE)
      ret <- FUN(response[index], weights[index])
      ### no information about i in fitted
      if (all(!index)) ret[1] <- NA
      return(ret)
    }))
    names(tab) <- as.character(sort(unique(id)))
  }
  if (inherits(tab[[1]], "function") && !is.null(at))
    tab <- lapply(tab, function(f) f(at))
  tn <- names(tab)
  dim(tab) <- NULL
  names(tab) <- tn

  tab
}


### simplify structure of predictions
.simplify_pred <- function(tab, id, nam) {

  if (all(sapply(tab, length) == 1) & all(sapply(tab, is.atomic))) {
    ret <- do.call("c", tab)
    names(ret) <- names(tab)
    ret <- if (is.factor(tab[[1]]))
      factor(ret[as.character(id)], levels = 1:length(levels(tab[[1]])),
             labels = levels(tab[[1]]), ordered = is.ordered(tab[[1]]))
    else
      ret[as.character(id)]
    names(ret) <- nam
  } else if (length(unique(sapply(tab, length))) == 1 &
             all(sapply(tab, is.numeric))) {
    ret <- matrix(unlist(tab), nrow = length(tab), byrow = TRUE)
    colnames(ret) <- names(tab[[1]])
    rownames(ret) <- names(tab)
    ret <- ret[as.character(id),, drop = FALSE]
    rownames(ret) <- nam
  } else {
    ret <- tab[as.character(id)]
    names(ret) <- nam
  }
  ret
}

data_party <- function(party, id = 1L)
  UseMethod("data_party")

data_party.default <- function(party, id = 1L) {

  extract <- function(id) {
    if(is.null(party$fitted))
      if(nrow(party$data) == 0) return(NULL)
    else
      stop("cannot subset data without fitted ids")

    ### which terminal nodes follow node number id?
    nt <- nodeids(party, id, terminal = TRUE)
    wi <- party$fitted[["(fitted)"]] %in% nt

    ret <- if (nrow(party$data) == 0)
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

width.party <- function(x, ...) {
  width(node_party(x), ...)
}

depth.party <- function(x, root = FALSE, ...) {
  depth(node_party(x), root = root, ...)
}

getCall.party <- function(x, ...) {
  x$info$call
}

getCall.constparties <- function(x, ...) {
  x$info$call
}

formula.party <- function(x, ...) {
  x <- terms(x)
  NextMethod()
}

model.frame.party <- function(formula, ...)
{
  mf <- formula$data
  if(nrow(mf) > 0L) return(mf)

  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0L)]
  mf <- getCall(formula)
  mf <- mf[c(1L, match(c("formula", "data", "subset", "na.action"), names(mf), 0L))]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf[names(nargs)] <- nargs
  if(is.null(env <- environment(terms(formula)))) env <- parent.frame()
  eval(mf, env)
}


nodeprune <- function(x, ids, ...)
  UseMethod("nodeprune")



nodeprune.partynode <- function(x, ids, ...) {

  stopifnot(ids %in% nodeids(x))

  ### compute indices path to each node
  ### to be pruned off
  idxs <- lapply(ids, .get_path, obj = x)

  ### [[.partynode is NOT [[.list
  cls <- class(x)
  x <- unclass(x)

  for (i in 1:length(idxs)) {
    ## path to be pruned
    idx <- idxs[[i]]
    if(!is.null(idx)) {
      ### check if we already pruned-off this node
      tmp <- try(x[[idx]], silent = TRUE)
      if(inherits(tmp, "try-error")) next()
      ### prune node by introducing a "new" terminal node
      x[[idx]] <- partynode(id = id_node(tmp), info = info_node(tmp))
    } else {
      ## if idx path is NULL prune everything
      x[2L:4L] <- NULL
    }
  }

  class(x) <- cls
  return(as.partynode(x, from = 1L))
}

nodeprune.default <- function(x, ids, ...)
  stop("No", sQuote("nodeprune"), "method for class", class(x), "implemented")

.list.rules.party <- function(x, i = NULL, ...) {
  if (is.null(i)) i <- nodeids(x, terminal = TRUE)
  if (length(i) > 1) {
    ret <- sapply(i, .list.rules.party, x = x)
    names(ret) <- if (is.character(i)) i else names(x)[i]
    return(ret)
  }
  if (is.character(i) && !is.null(names(x)))
    i <- which(names(x) %in% i)
  stopifnot(length(i) == 1 & is.numeric(i))
  stopifnot(i <= length(x) & i >= 1)
  i <- as.integer(i)
  dat <- data_party(x, i)
  if (!is.null(x$fitted)) {
    findx <- which("(fitted)" == names(dat))[1]
    fit <- dat[,findx:ncol(dat), drop = FALSE]
    dat <- dat[,-(findx:ncol(dat)), drop = FALSE]
    if (ncol(dat) == 0)
      dat <- x$data
  } else {
    fit <- NULL
    dat <- x$data
  }

  rule <- c()

  recFun <- function(node) {
    if (id_node(node) == i) return(NULL)
    kid <- sapply(kids_node(node), id_node)
    whichkid <- max(which(kid <= i))
    split <- split_node(node)
    ivar <- varid_split(split)
    svar <- names(dat)[ivar]
    index <- index_split(split)
    if (is.factor(dat[, svar])) {
      if (is.null(index))
        index <- ((1:nlevels(dat[, svar])) > breaks_split(split)) + 1
      slevels <- levels(dat[, svar])[index == whichkid]
      srule <- paste(svar, " %in% c(\"",
                     paste(slevels, collapse = "\", \"", sep = ""), "\")",
                     sep = "")
    } else {
      if (is.null(index)) index <- 1:length(kid)
      breaks <- cbind(c(-Inf, breaks_split(split)),
                      c(breaks_split(split), Inf))
      sbreak <- breaks[index == whichkid,]
      right <- right_split(split)
      srule <- c()
      if (is.finite(sbreak[1]))
        srule <- c(srule,
                   paste(svar, ifelse(right, ">", ">="), sbreak[1]))
      if (is.finite(sbreak[2]))
        srule <- c(srule,
                   paste(svar, ifelse(right, "<=", "<"), sbreak[2]))
      srule <- paste(srule, collapse = " & ")
    }
    rule <<- c(rule, srule)
    return(recFun(node[[whichkid]]))
  }
  node <- recFun(node_party(x))
  paste(rule, collapse = " & ")
}

