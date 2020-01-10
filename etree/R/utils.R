
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
  if (!is.null(index) & is.null(centroids_split(split)))
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


# plot --------------------------------------------------------------------

.nobs_party <- function(party, id = 1L) {
  dat <- data_party(party, id = id)
  if("(weights)" %in% names(dat)) sum(dat[["(weights)"]]) else NROW(dat)
}

node_inner <- function(obj, id = TRUE, pval = TRUE, abbreviate = FALSE, fill = "white", gp = gpar())
{
  meta <- obj$data
  nam <- names(obj)

  extract_label <- function(node) {
    if(is.terminal(node)) return(rep.int("", 2L))

    varlab <- character_split(split_node(node), meta)$name
    if(abbreviate > 0L) varlab <- abbreviate(varlab, as.integer(abbreviate))

    ## FIXME: make more flexible rather than special-casing p-value
    if(pval) {
      nullna <- function(x) is.null(x) || is.na(x)
      pval <- suppressWarnings(try(!nullna(info_node(node)$p.value), silent = TRUE))
      pval <- if(inherits(pval, "try-error")) FALSE else pval
    }
    if(pval) {
      pvalue <- node$info$p.value
      plab <- ifelse(pvalue < 10^(-3L),
                     paste("p <", 10^(-3L)),
                     paste("p =", round(pvalue, digits = 3L)))
    } else {
      plab <- ""
    }
    return(c(varlab, plab))
  }

  maxstr <- function(node) {
    lab <- extract_label(node)
    klab <- if(is.terminal(node)) "" else unlist(lapply(kids_node(node), maxstr))
    lab <- c(lab, klab)
    lab <- unlist(lapply(lab, function(x) strsplit(x, "\n")))
    lab <- lab[which.max(nchar(lab))]
    if(length(lab) < 1L) lab <- ""
    return(lab)
  }

  nstr <- maxstr(node_party(obj))
  if(nchar(nstr) < 6) nstr <- "aAAAAa"

  ### panel function for the inner nodes
  rval <- function(node) {
    pushViewport(viewport(gp = gp, name = paste("node_inner", id_node(node), "_gpar", sep = "")))
    node_vp <- viewport(
      x = unit(0.5, "npc"),
      y = unit(0.5, "npc"),
      width = unit(1, "strwidth", nstr) * 1.3,
      height = unit(3, "lines"),
      name = paste("node_inner", id_node(node), sep = ""),
      gp = gp
    )
    pushViewport(node_vp)

    xell <- c(seq(0, 0.2, by = 0.01),
              seq(0.2, 0.8, by = 0.05),
              seq(0.8, 1, by = 0.01))
    yell <- sqrt(xell * (1-xell))

    lab <- extract_label(node)
    fill <- rep(fill, length.out = 2L)

    grid.polygon(x = unit(c(xell, rev(xell)), "npc"),
                 y = unit(c(yell, -yell)+0.5, "npc"),
                 gp = gpar(fill = fill[1]))

    ## FIXME: something more general instead of pval ?
    grid.text(lab[1L], y = unit(1.5 + 0.5 * (lab[2L] != ""), "lines"))
    if(lab[2L] != "") grid.text(lab[2L], y = unit(1, "lines"))

    if(id) {
      nodeIDvp <- viewport(x = unit(0.5, "npc"), y = unit(1, "npc"),
                           width = max(unit(1, "lines"), unit(1.3, "strwidth", nam[id_node(node)])),
                           height = max(unit(1, "lines"), unit(1.3, "strheight", nam[id_node(node)])))
      pushViewport(nodeIDvp)
      grid.rect(gp = gpar(fill = fill[2]))
      grid.text(nam[id_node(node)])
      popViewport()
    }
    upViewport(2)
  }

  return(rval)
}
class(node_inner) <- "grapcon_generator"

node_terminal <- function(obj,
                          digits = 3,
                          abbreviate = FALSE,
                          fill = c("lightgray", "white"),
                          id = TRUE,
                          just = c("center", "top"),
                          top = 0.85,
                          align = c("center", "left", "right"),
                          gp = NULL,
                          FUN = NULL,
                          height = NULL,
                          width = NULL)
{
  nam <- names(obj)

  extract_label <- function(node) formatinfo_node(node, FUN = FUN, default = c("terminal", "node"))

  maxstr <- function(node) {
    lab <- extract_label(node)
    klab <- if(is.terminal(node)) "" else unlist(lapply(kids_node(node), maxstr))
    lab <- c(lab, klab)
    lab <- try(unlist(lapply(lab, function(x) strsplit(x, "\n"))), silent = TRUE)
    if(inherits(lab, "try-error")) {
      paste(rep("a", 9L), collapse = "") ## FIXME: completely ad-hoc: possibly throw warning?
    } else {
      return(lab[which.max(nchar(lab))])
    }
  }

  nstr <- if(is.null(width)) maxstr(node_party(obj)) else paste(rep("a", width), collapse = "")

  just <- match.arg(just[1L], c("center", "centre", "top"))
  if(just == "centre") just <- "center"
  align <- match.arg(align[1L], c("center", "centre", "left", "right"))
  if(align == "centre") align <- "center"

  ### panel function for simple n, Y terminal node labeling
  rval <- function(node) {
    fill <- rep(fill, length.out = 2)

    lab <- extract_label(node)

    ## if gp is set, then an additional viewport may be
    ## required to appropriately evaluate strwidth unit
    if(!is.null(gp)) {
      outer_vp <- viewport(gp = gp)
      pushViewport(outer_vp)
    }

    if(is.null(height)) height <- length(lab) + 1L

    node_vp <- viewport(x = unit(0.5, "npc"),
                        y = unit(if(just == "top") top else 0.5, "npc"),
                        just = c("center", just),
                        width = unit(1, "strwidth", nstr) * 1.1,
                        height = unit(height, "lines"),
                        name = paste("node_terminal", id_node(node), sep = ""),
                        gp = if(is.null(gp)) gpar() else gp
    )
    pushViewport(node_vp)

    grid.rect(gp = gpar(fill = fill[1]))

    for(i in seq_along(lab)) grid.text(
      x = switch(align,
                 "center" = unit(0.5, "npc"),
                 "left"   = unit(1, "strwidth", "a"),
                 "right"  = unit(1, "npc") - unit(1, "strwidth", "a")),
      y = unit(length(lab) - i + 1, "lines"), lab[i], just = align)

    if(id) {
      nodeIDvp <- viewport(x = unit(0.5, "npc"), y = unit(1, "npc"),
                           width = max(unit(1, "lines"), unit(1.3, "strwidth", nam[id_node(node)])),
                           height = max(unit(1, "lines"), unit(1.3, "strheight", nam[id_node(node)])))
      pushViewport(nodeIDvp)
      grid.rect(gp = gpar(fill = fill[2], lty = "solid"))
      grid.text(nam[id_node(node)])
      popViewport()
    }

    if(is.null(gp)) upViewport() else upViewport(2)
  }
  return(rval)
}
class(node_terminal) <- "grapcon_generator"

edge_simple <- function(obj, digits = 3, abbreviate = FALSE,
                        justmin = Inf, just = c("alternate", "increasing", "decreasing", "equal"),
                        fill = "white")
{
  meta <- obj$data

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
    grid.rect(y = y, gp = gpar(fill = fill, col = 0), width = unit(1, "strwidth", split))
    grid.text(split, y = y, just = "center")
  }
}
class(edge_simple) <- "grapcon_generator"

.plot_node <- function(node, xlim, ylim, nx, ny,
                       terminal_panel, inner_panel, edge_panel,
                       tnex = 2, drop_terminal = TRUE, debug = FALSE) {

  ### the workhorse for plotting trees

  ### set up viewport for terminal node
  if (is.terminal(node)) {
    x <- xlim[1] + diff(xlim)/2
    y <- ylim[1] + 0.5

    tn_vp <- viewport(x = unit(x, "native"),
                      y = unit(y, "native") - unit(0.5, "lines"),
                      width = unit(1, "native"),
                      height = unit(tnex, "native") - unit(1, "lines"),
                      just = c("center", "top"),
                      name = paste("Node", id_node(node), sep = ""))
    pushViewport(tn_vp)
    if (debug)
      grid.rect(gp = gpar(lty = "dotted", col = 4))
    terminal_panel(node)
    upViewport()
    return(NULL)
  }

  ## convenience function for computing relative position of splitting node
  pos_frac <- function(node) {
    if(is.terminal(node)) 0.5 else {
      width_kids <- sapply(kids_node(node), width)
      nk <- length(width_kids)
      rval <- if(nk %% 2 == 0) sum(width_kids[1:(nk/2)]) else
        mean(cumsum(width_kids)[nk/2 + c(-0.5, 0.5)])
      rval/sum(width_kids)
    }
  }

  ## extract information
  split <- split_node(node)
  kids <- kids_node(node)
  width_kids <- sapply(kids, width)
  nk <- length(width_kids)

  ### position of inner node
  x0 <- xlim[1] + pos_frac(node) * diff(xlim)
  y0 <- max(ylim)

  ### relative positions of kids
  xfrac <- sapply(kids, pos_frac)
  x1lim <- xlim[1] + cumsum(c(0, width_kids))/sum(width_kids) * diff(xlim)
  x1 <- x1lim[1:nk] + xfrac * diff(x1lim)
  if (!drop_terminal) {
    y1 <- rep(y0 - 1, nk)
  } else {
    y1 <- ifelse(sapply(kids, is.terminal), tnex - 0.5, y0 - 1)
  }

  ### draw edges
  for(i in 1:nk) grid.lines(x = unit(c(x0, x1[i]), "native"), y = unit(c(y0, y1[i]), "native"))

  ### create viewport for inner node
  in_vp <- viewport(x = unit(x0, "native"),
                    y = unit(y0, "native"),
                    width = unit(1, "native"),
                    height = unit(1, "native") - unit(1, "lines"),
                    name = paste("Node", id_node(node), sep = ""))
  pushViewport(in_vp)
  if(debug) grid.rect(gp = gpar(lty = "dotted"))
  inner_panel(node)
  upViewport()

  ### position of labels
  y1max <- max(y1)
  ypos <- y0 - (y0 - y1max) * 0.5
  xpos <- x0 - (x0 - x1) * 0.5 * (y0 - y1max)/(y0 - y1)

  ### setup labels
  for(i in 1:nk) {
    sp_vp <- viewport(x = unit(xpos[i], "native"),
                      y = unit(ypos, "native"),
                      width = unit(diff(x1lim)[i], "native"),
                      height = unit(1, "lines"),
                      name =  paste("edge", id_node(node), "-", i, sep = ""))
    pushViewport(sp_vp)
    if(debug) grid.rect(gp = gpar(lty = "dotted", col = 2))
    edge_simple(node, i)
    upViewport()
  }

  ## call workhorse for kids
  for(i in 1:nk) .plot_node(kids[[i]],
                            c(x1lim[i], x1lim[i+1]), c(y1[i], 1), nx, ny,
                            terminal_panel, inner_panel, edge_panel,
                            tnex = tnex, drop_terminal = drop_terminal, debug = debug)
}


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

node_barplot <- function(obj,
                         col = "black",
                         fill = NULL,
                         bg = "white",
                         beside = NULL,
                         ymax = NULL,
                         ylines = NULL,
                         widths = 1,
                         gap = NULL,
                         reverse = NULL,
                         rot = 0,
                         just = c("center", "top"),
                         id = TRUE,
                         mainlab = NULL,
                         text = c("none", "horizontal", "vertical"),
                         gp = gpar())
{
  ## extract response
  y <- obj$fitted[["(response)"]]
  stopifnot(is.factor(y) || isTRUE(all.equal(round(y), y)) || is.data.frame(y))

  ## FIXME: This could be avoided by
  ##   predict_party(obj, nodeids(obj, terminal = TRUE), type = "prob")
  ## but only for terminal nodes                  ^^^^
  probs_and_n <- function(x) {
    y1 <- x$fitted[["(response)"]]
    if(!is.factor(y1)) {
      if(is.data.frame(y1)) {
        y1 <- t(as.matrix(y1))
      } else {
        y1 <- factor(y1, levels = min(y):max(y))
      }
    }
    w <- x$fitted[["(weights)"]]
    if(is.null(w)) w <- rep.int(1L, length(y1))
    sumw <- if(is.factor(y1)) tapply(w, y1, sum) else drop(y1 %*% w)
    sumw[is.na(sumw)] <- 0
    prob <- c(sumw/sum(w), sum(w))
    names(prob) <- c(if(is.factor(y1)) levels(y1) else rownames(y1), "nobs")
    prob
  }
  probs <- do.call("rbind", nodeapply(obj, nodeids(obj), probs_and_n, by_node = FALSE))
  nobs <- probs[, "nobs"]
  probs <- probs[, -ncol(probs), drop = FALSE]

  if(is.factor(y)) {
    ylevels <- levels(y)
    if(is.null(beside)) beside <- if(length(ylevels) < 3L) FALSE else TRUE
    if(is.null(ymax)) ymax <- if(beside) 1.1 else 1
    if(is.null(gap)) gap <- if(beside) 0.1 else 0
  } else {
    if(is.null(beside)) beside <- TRUE
    if(is.null(ymax)) ymax <- if(beside) max(probs) * 1.1 else max(probs)
    ylevels <- colnames(probs)
    if(length(ylevels) < 2) ylevels <- ""
    if(is.null(gap)) gap <- if(beside) 0.1 else 0
  }
  if(is.null(reverse)) reverse <- !beside
  if(is.null(fill)) fill <- gray.colors(length(ylevels))
  if(is.null(ylines)) ylines <- if(beside) c(3, 2) else c(1.5, 2.5)

  ## text labels?
  if(isTRUE(text)) text <- "horizontal"
  if(!is.character(text)) text <- "none"
  text <- match.arg(text, c("none", "horizontal", "vertical"))

  ### panel function for barplots in nodes
  rval <- function(node) {

    ## id
    nid <- id_node(node)

    ## parameter setup
    pred <- probs[nid,]
    if(reverse) {
      pred <- rev(pred)
      ylevels <- rev(ylevels)
    }
    np <- length(pred)
    nc <- if(beside) np else 1

    fill <- rep(fill, length.out = np)
    widths <- rep(widths, length.out = nc)
    col <- rep(col, length.out = nc)
    ylines <- rep(ylines, length.out = 2)

    gap <- gap * sum(widths)
    yscale <- c(0, ymax)
    xscale <- c(0, sum(widths) + (nc+1)*gap)

    top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                                            widths = unit(c(ylines[1], 1, ylines[2]), c("lines", "null", "lines")),
                                            heights = unit(c(1, 1), c("lines", "null"))),
                       width = unit(1, "npc"),
                       height = unit(1, "npc") - unit(2, "lines"),
                       name = paste0("node_barplot", nid),
                       gp = gp)

    pushViewport(top_vp)
    grid.rect(gp = gpar(fill = bg, col = 0))

    ## main title
    top <- viewport(layout.pos.col=2, layout.pos.row=1)
    pushViewport(top)
    if (is.null(mainlab)) {
      mainlab <- if(id) {
        function(id, nobs) sprintf("Node %s (n = %s)", id, nobs)
      } else {
        function(id, nobs) sprintf("n = %s", nobs)
      }
    }
    if (is.function(mainlab)) {
      mainlab <- mainlab(names(obj)[nid], nobs[nid])
    }
    grid.text(mainlab)
    popViewport()

    plot <- viewport(layout.pos.col=2, layout.pos.row=2,
                     xscale=xscale, yscale=yscale,
                     name = paste0("node_barplot", node$nodeID, "plot"),
                     clip = FALSE)

    pushViewport(plot)

    if(beside) {
      xcenter <- cumsum(widths+gap) - widths/2
      if(length(xcenter) > 1) grid.xaxis(at = xcenter, label = FALSE)
      grid.text(ylevels, x = xcenter, y = unit(-1, "lines"),
                just = just, rot = rot,
                default.units = "native", check.overlap = TRUE)
      grid.yaxis()
      grid.rect(gp = gpar(fill = "transparent"))
      grid.clip()
      for (i in 1:np) {
        grid.rect(x = xcenter[i], y = 0, height = pred[i],
                  width = widths[i],
                  just = c("center", "bottom"), default.units = "native",
                  gp = gpar(col = col[i], fill = fill[i]))
        if(text != "none") {
          grid.text(x = xcenter[i], y = pred[i] + 0.025,
                    label = paste(format(round(100 * pred[i], 1), nsmall = 1), "%", sep = ""),
                    just = if(text == "horizontal") c("center", "bottom") else c("left", "center"),
                    rot = if(text == "horizontal") 0 else 90,
                    default.units = "native")
        }
      }
    } else {
      ycenter <- cumsum(pred) - pred

      if(np > 1) {
        grid.text(ylevels[1], x = unit(-1, "lines"), y = 0,
                  just = c("left", "center"), rot = 90,
                  default.units = "native", check.overlap = TRUE)
        grid.text(ylevels[np], x = unit(-1, "lines"), y = ymax,
                  just = c("right", "center"), rot = 90,
                  default.units = "native", check.overlap = TRUE)
      }
      if(np > 2) {
        grid.text(ylevels[-c(1,np)], x = unit(-1, "lines"), y = ycenter[-c(1,np)],
                  just = "center", rot = 90,
                  default.units = "native", check.overlap = TRUE)
      }
      grid.yaxis(main = FALSE)

      grid.clip()
      grid.rect(gp = gpar(fill = "transparent"))
      for (i in 1:np) {
        grid.rect(x = xscale[2]/2, y = ycenter[i], height = min(pred[i], ymax - ycenter[i]),
                  width = widths[1],
                  just = c("center", "bottom"), default.units = "native",
                  gp = gpar(col = col[i], fill = fill[i]))
      }
    }
    grid.rect(gp = gpar(fill = "transparent"))


    upViewport(2)
  }

  return(rval)
}
class(node_barplot) <- "grapcon_generator"

node_boxplot <- function(obj,
                         col = "black",
                         fill = "lightgray",
                         bg = "white",
                         width = 0.5,
                         yscale = NULL,
                         ylines = 3,
                         cex = 0.5,
                         id = TRUE,
                         mainlab = NULL,
                         gp = gpar())
{
  y <- obj$fitted[["(response)"]]
  stopifnot(is.numeric(y))

  if (is.null(yscale))
    yscale <- range(y) + c(-0.1, 0.1) * diff(range(y))

  ### panel function for boxplots in nodes
  rval <- function(node) {

    ## extract data
    nid <- id_node(node)
    dat <- data_party(obj, nid)
    yn <- dat[["(response)"]]
    wn <- dat[["(weights)"]]
    if(is.null(wn)) wn <- rep(1, length(yn))

    ## parameter setup
    x <- boxplot(rep.int(yn, wn), plot = FALSE)

    top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                                            widths = unit(c(ylines, 1, 1),
                                                          c("lines", "null", "lines")),
                                            heights = unit(c(1, 1), c("lines", "null"))),
                       width = unit(1, "npc"),
                       height = unit(1, "npc") - unit(2, "lines"),
                       name = paste("node_boxplot", nid, sep = ""),
                       gp = gp)

    pushViewport(top_vp)
    grid.rect(gp = gpar(fill = bg, col = 0))

    ## main title
    top <- viewport(layout.pos.col=2, layout.pos.row=1)
    pushViewport(top)
    if (is.null(mainlab)) {
      mainlab <- if(id) {
        function(id, nobs) sprintf("Node %s (n = %s)", id, nobs)
      } else {
        function(id, nobs) sprintf("n = %s", nobs)
      }
    }
    if (is.function(mainlab)) {
      mainlab <- mainlab(names(obj)[nid], sum(wn))
    }
    grid.text(mainlab)
    popViewport()

    plot <- viewport(layout.pos.col = 2, layout.pos.row = 2,
                     xscale = c(0, 1), yscale = yscale,
                     name = paste0("node_boxplot", nid, "plot"),
                     clip = FALSE)

    pushViewport(plot)

    grid.yaxis()
    grid.rect(gp = gpar(fill = "transparent"))
    grid.clip()

    xl <- 0.5 - width/4
    xr <- 0.5 + width/4

    ## box & whiskers
    grid.lines(unit(c(xl, xr), "npc"),
               unit(x$stats[1], "native"), gp = gpar(col = col))
    grid.lines(unit(0.5, "npc"),
               unit(x$stats[1:2], "native"), gp = gpar(col = col, lty = 2))
    grid.rect(unit(0.5, "npc"), unit(x$stats[2], "native"),
              width = unit(width, "npc"), height = unit(diff(x$stats[c(2, 4)]), "native"),
              just = c("center", "bottom"),
              gp = gpar(col = col, fill = fill))
    grid.lines(unit(c(0.5 - width/2, 0.5+width/2), "npc"),
               unit(x$stats[3], "native"), gp = gpar(col = col, lwd = 2))
    grid.lines(unit(0.5, "npc"), unit(x$stats[4:5], "native"),
               gp = gpar(col = col, lty = 2))
    grid.lines(unit(c(xl, xr), "npc"), unit(x$stats[5], "native"),
               gp = gpar(col = col))

    ## outlier
    n <- length(x$out)
    if (n > 0) {
      index <- 1:n ## which(x$out > yscale[1] & x$out < yscale[2])
      if (length(index) > 0)
        grid.points(unit(rep.int(0.5, length(index)), "npc"),
                    unit(x$out[index], "native"),
                    size = unit(cex, "char"), gp = gpar(col = col))
    }

    upViewport(2)
  }

  return(rval)
}
class(node_boxplot) <- "grapcon_generator"

node_surv <- function(obj, col = "black", bg = "white", yscale = c(0, 1), ylines = 2,
                      id = TRUE, mainlab = NULL, gp = gpar(), ...)
{
  ## extract response
  y <- obj$fitted[["(response)"]]
  stopifnot(inherits(y, "Surv"))

  ## helper functions
  mysurvfit <- function(y, weights, ...)
    survfit(y ~ 1, weights = weights)
  ### structure(
  ###   survival:::survfitKM(x = gl(1, NROW(y)), y = y, casewt = weights, ...),
  ### class = "survfit")

  dostep <- function(x, y) {
    ### create a step function based on x, y coordinates
    ### modified from `survival:print.survfit'
    if (is.na(x[1] + y[1])) {
      x <- x[-1]
      y <- y[-1]
    }
    n <- length(x)
    if (n > 2) {
      # replace verbose horizonal sequences like
      # (1, .2), (1.4, .2), (1.8, .2), (2.3, .2), (2.9, .2), (3, .1)
      # with (1, .2), (3, .1).  They are slow, and can smear the looks
      # of the line type.
      dupy <- c(TRUE, diff(y[-n]) !=0, TRUE)
      n2 <- sum(dupy)

      #create a step function
      xrep <- rep(x[dupy], c(1, rep(2, n2-1)))
      yrep <- rep(y[dupy], c(rep(2, n2-1), 1))
      RET <- list(x = xrep, y = yrep)
    } else {
      if (n == 1) {
        RET <- list(x = x, y = y)
      } else {
        RET <- list(x = x[c(1,2,2)], y = y[c(1,1,2)])
      }
    }
    return(RET)
  }

  ### panel function for Kaplan-Meier curves in nodes
  rval <- function(node) {

    ## extract data
    nid <- id_node(node)
    dat <- data_party(obj, nid)
    yn <- dat[["(response)"]]
    wn <- dat[["(weights)"]]
    if(is.null(wn)) wn <- rep(1, NROW(yn))

    ## get Kaplan-Meier curver in node
    km <- mysurvfit(yn, weights = wn, ...)
    a <- dostep(km$time, km$surv)

    ## set up plot
    yscale <- yscale
    xscale <- c(0, max(y[,1]))

    top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                                            widths = unit(c(ylines, 1, 1),
                                                          c("lines", "null", "lines")),
                                            heights = unit(c(1, 1), c("lines", "null"))),
                       width = unit(1, "npc"),
                       height = unit(1, "npc") - unit(2, "lines"),
                       name = paste("node_surv", nid, sep = ""), gp = gp)

    pushViewport(top_vp)
    grid.rect(gp = gpar(fill = bg, col = 0))

    ## main title
    top <- viewport(layout.pos.col=2, layout.pos.row=1)
    pushViewport(top)
    if (is.null(mainlab)) {
      mainlab <- if(id) {
        function(id, nobs) sprintf("Node %s (n = %s)", id, nobs)
      } else {
        function(id, nobs) sprintf("n = %s", nobs)
      }
    }
    if (is.function(mainlab)) {
      mainlab <- mainlab(nid, sum(wn))
    }
    grid.text(mainlab)
    popViewport()

    plot <- viewport(layout.pos.col=2, layout.pos.row=2,
                     xscale=xscale, yscale = yscale,
                     name = paste0("node_surv", nid, "plot"),
                     clip = FALSE)

    pushViewport(plot)
    grid.xaxis()
    grid.yaxis()
    grid.rect(gp = gpar(fill = "transparent"))
    grid.clip()
    grid.lines(unit(a$x, "native"), unit(a$y, "native"), gp = gpar(col = col))
    upViewport(2)
  }

  return(rval)
}
class(node_surv) <- "grapcon_generator"

node_ecdf <- function(obj, col = "black", bg = "white", ylines = 2,
                      id = TRUE, mainlab = NULL, gp = gpar(), ...)
{
  ## extract response
  y <- obj$fitted[["(response)"]]
  stopifnot(inherits(y, "numeric") || inherits(y, "integer"))

  dostep <- function(f) {
    x <- knots(f)
    y <- f(x)
    ### create a step function based on x, y coordinates
    ### modified from `survival:print.survfit'
    if (is.na(x[1] + y[1])) {
      x <- x[-1]
      y <- y[-1]
    }
    n <- length(x)
    if (n > 2) {
      # replace verbose horizonal sequences like
      # (1, .2), (1.4, .2), (1.8, .2), (2.3, .2), (2.9, .2), (3, .1)
      # with (1, .2), (3, .1).  They are slow, and can smear the looks
      # of the line type.
      dupy <- c(TRUE, diff(y[-n]) !=0, TRUE)
      n2 <- sum(dupy)

      #create a step function
      xrep <- rep(x[dupy], c(1, rep(2, n2-1)))
      yrep <- rep(y[dupy], c(rep(2, n2-1), 1))
      RET <- list(x = xrep, y = yrep)
    } else {
      if (n == 1) {
        RET <- list(x = x, y = y)
      } else {
        RET <- list(x = x[c(1,2,2)], y = y[c(1,1,2)])
      }
    }
    return(RET)
  }

  ### panel function for ecdf in nodes
  rval <- function(node) {

    ## extract data
    nid <- id_node(node)
    dat <- data_party(obj, nid)
    yn <- dat[["(response)"]]
    wn <- dat[["(weights)"]]
    if(is.null(wn)) wn <- rep(1, NROW(yn))

    ## get ecdf in node
    f <- .pred_ecdf(yn, wn)
    a <- dostep(f)

    ## set up plot
    yscale <- c(0, 1)
    xscale <- range(y)
    a$x <- c(xscale[1], a$x[1], a$x, xscale[2])
    a$x <- a$x - min(a$x)
    a$x <- a$x / max(a$x)
    a$y <- c(0, 0, a$y, 1)

    top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                                            widths = unit(c(ylines, 1, 1),
                                                          c("lines", "null", "lines")),
                                            heights = unit(c(1, 1), c("lines", "null"))),
                       width = unit(1, "npc"),
                       height = unit(1, "npc") - unit(2, "lines"),
                       name = paste("node_ecdf", nid, sep = ""), gp = gp)

    pushViewport(top_vp)
    grid.rect(gp = gpar(fill = bg, col = 0))

    ## main title
    top <- viewport(layout.pos.col=2, layout.pos.row=1)
    pushViewport(top)
    if (is.null(mainlab)) {
      mainlab <- if(id) {
        function(id, nobs) sprintf("Node %s (n = %s)", id, nobs)
      } else {
        function(id, nobs) sprintf("n = %s", nobs)
      }
    }
    if (is.function(mainlab)) {
      mainlab <- mainlab(nid, sum(wn))
    }
    grid.text(mainlab)
    popViewport()

    plot <- viewport(layout.pos.col=2, layout.pos.row=2,
                     xscale=xscale, yscale=yscale,
                     name = paste0("node_surv", nid, "plot"),
                     clip = FALSE)

    pushViewport(plot)
    grid.xaxis()
    grid.yaxis()
    grid.rect(gp = gpar(fill = "transparent"))
    grid.clip()
    grid.lines(a$x, a$y, gp = gpar(col = col))
    upViewport(2)
  }

  return(rval)
}
class(node_ecdf) <- "grapcon_generator"



node_mvar <- function(obj, which = NULL, id = TRUE, pop = TRUE, ylines = NULL, mainlab = NULL, varlab = TRUE, bg = "white", ...)
{
  ## obtain dependent variables
  y <- obj$fitted[["(response)"]]

  ## fitted node ids
  fitted <- obj$fitted[["(fitted)"]]

  ## number of panels needed
  if(is.null(which)) which <- 1L:NCOL(y)
  k <- length(which)

  rval <- function(node) {

    tid <- id_node(node)
    nobs <- .nobs_party(obj, id = tid)

    ## set up top viewport
    top_vp <- viewport(layout = grid.layout(nrow = k, ncol = 2,
                                            widths = unit(c(ylines, 1), c("lines", "null")), heights = unit(k, "null")),
                       width = unit(1, "npc"), height = unit(1, "npc") - unit(2, "lines"),
                       name = paste("node_mvar", tid, sep = ""))
    pushViewport(top_vp)
    grid.rect(gp = gpar(fill = bg, col = 0))

    ## main title
    if (is.null(mainlab)) {
      mainlab <- if(id) {
        function(id, nobs) sprintf("Node %s (n = %s)", id, nobs)
      } else {
        function(id, nobs) sprintf("n = %s", nobs)
      }
    }
    if (is.function(mainlab)) {
      mainlab <- mainlab(tid, nobs)
    }

    for(i in 1L:k) {
      tmp <- obj
      tmp$fitted[["(response)"]] <- y[,which[i]]
      if(varlab) {
        nm <- names(y)[which[i]]
        if(i == 1L) nm <- paste(mainlab, nm, sep = ": ")
      } else {
        nm <- if(i == 1L) mainlab else ""
      }
      pfun <- switch(sapply(y, class)[which[i]],
                     "Surv" = node_surv(tmp, id = id, mainlab = nm, ...),
                     "factor" = node_barplot(tmp, id = id, mainlab = nm,  ...),
                     "ordered" = node_barplot(tmp, id = id, mainlab = nm, ...),
                     node_boxplot(tmp, id = id, mainlab = nm, ...))
      ## select panel
      plot_vpi <- viewport(layout.pos.col = 2L, layout.pos.row = i)
      pushViewport(plot_vpi)

      ## call panel function
      pfun(node)

      if(pop) popViewport() else upViewport()
    }
    if(pop) popViewport() else upViewport()
  }

  return(rval)
}
class(node_mvar) <- "grapcon_generator"

