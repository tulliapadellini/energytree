

# Main function -----------------------------------------------------------

etree <- function(response,
                  covariates,
                  weights = NULL,
                  minbucket = 5,
                  alpha = 0.05,
                  R = 1000,
                  split_type = 'coeff',
                  coeff_split_type = 'test',
                  p_adjust_method = 'fdr',
                  supervised = TRUE,
                  random_covs = NULL) {

  # Check whether covariates is a list
  if (!is.list(covariates)) {
    stop("Argument 'covariates' must be provided as a list")
  }

  # If the case weights are not provided, they are all initialized as 1
  if (is.null(weights)) {
    if (isTRUE(identical(names(response), c('response', 'response_dist')))) {
      weights <- rep(1L, as.numeric(length(response$response)))
    } else {
      weights <- rep(1L, as.numeric(length(response)))
    }
  }

  # Large list both for covariates and for response
  #(get them, if already computed in eforest; otherwise, compute them)
  if (isTRUE(identical(names(covariates), c('cov', 'newcov', 'dist')))) {

    # Large list with covariates, newcovariates and distances
    covariates_large <- covariates

    # Large list with response and the corresponding distances
    response_large <- response

    # New covariates and response (needed divided to build the df used by party)
    newcovariates <- covariates_large$newcov
    response <- response_large$response

  } else {

    # New covariates
    newcovariates <- .create_newcov(covariates = covariates,
                                    response = response,
                                    split_type = split_type)

    # Distances
    cov_distance <- lapply(covariates, dist_comp)

    # Large list with covariates, newcovariates and distances
    covariates_large <- list('cov' = covariates,
                             'newcov' = newcovariates,
                             'dist' = cov_distance)

    # Large list with response and the corresponding distances
    response_large <- list('response' = response,
                           'response_dist' = dist_comp(response))

  }

  # Check that response is factor or numeric
  if (!is.factor(response_large$response) &&
      !is.numeric(response_large$response)) {
    stop("Argument 'response' must be provided either as a factor or as an
         object of mode 'numeric'")
  }


  # Grow the tree (find the split rules)
  nodes <- growtree(id = 1L,
                    response = response_large,
                    covariates = covariates_large,
                    weights = weights,
                    minbucket = minbucket,
                    alpha = alpha,
                    R = R,
                    split_type = split_type,
                    coeff_split_type = coeff_split_type,
                    p_adjust_method = p_adjust_method,
                    random_covs = random_covs)

  # Actually perform the splits
  fitted_obs <- fitted_node(nodes, data = newcovariates)

  # Return a rich constparty object
  obj <- party(nodes,
               data = newcovariates,
               fitted = data.frame("(fitted)" = fitted_obs,
                                   "(response)" = if (supervised) response
                                   else fitted_obs,
                                   check.names = FALSE),
               terms = terms(response ~ ., data = newcovariates))
  etree_obj <- as.constparty(obj)
  attr(etree_obj, 'split_type') <- split_type     #used in predict.party

  return(etree_obj)

}


.create_newcov <- function(covariates, response, split_type) {

  newcovariates <- lapply(covariates, function(j) {

    switch(class(j),

           fdata = {

             if (split_type == "coeff") {

               fdata_est <- fda.usc::optim.basis(j, numbasis =
                                                   floor(seq(4,
                                                             ncol(j) / 2,
                                                             len = 10)))
               #seq starts from 4 as it is the smallest ok number (norder is 4,
               #and nbasis has to be >= norder -- cf. fda::create.bspline.basis)
               coefs <- fda.usc::fdata2fd(fdata_est$fdata.est,
                                          type.basis = "bspline",
                                          nbasis = fdata_est$numbasis.opt)$coefs
               newcov <- data.frame(t(coefs))
               names(newcov) <- seq_along(names(newcov))

             } else if (split_type == "cluster") {

               newcov <- as.factor(seq_along(response))

             }

             attr(newcov, 'cov_type') <- 'fdata'
             return(newcov)

           },

           list = if (all(sapply(j, function(x) attributes(x)$names)
                          == 'diagram')) {

             newcov <- as.factor(seq_along(response))
             attr(newcov, 'cov_type') <- 'diagram'
             return(newcov)

           } else if (all(sapply(j, class) == 'igraph')) {

             if (split_type == "coeff") {

               newcov <- graph_shell(j)

             } else if (split_type == "cluster") {

               newcov <- as.factor(seq_along(response))

             }

             attr(newcov, 'cov_type') <- 'graph'
             return(newcov)

           },

           # In all other cases:
           return(j)

    )
  })

  # Covariates name
  if (!is.null(names(covariates))) {
    names(newcovariates) <- names(covariates)
    #control if any name is void, i.e. if it is ''
    no_name <- which(sapply(names(covariates), function(n) n == '',
                            USE.NAMES = FALSE))
    names(newcovariates) <- replace(names(newcovariates),
                                    no_name,
                                    as.factor(seq_along(no_name)))
  } else {
    warning('No names available for covariates. Numbers are used instead.')
    names(newcovariates) <- seq_along(newcovariates)
  }

  # Return newcovariates
  return(newcovariates)

}



# Grow the tree ---------------------------------------------------------------

growtree <- function(id = 1L,
                     response,
                     covariates,
                     weights,
                     minbucket,
                     alpha,
                     R,
                     split_type,
                     coeff_split_type,
                     p_adjust_method,
                     random_covs) {

  # If the observations are less than <2*minbucket>, stop here (otherwise, there
  #would not be enough obs to have at least <minbucket> obs in each kid node)
  if (sum(weights) < 2 * minbucket) {
    return(partynode(id = id))
  }

  # Find the best split (variable selection & split point search)
  split <- findsplit(response = response,
                     covariates = covariates,
                     alpha = alpha,
                     R = R,
                     split_type = split_type,
                     coeff_split_type = coeff_split_type,
                     p_adjust_method = p_adjust_method,
                     random_covs = random_covs)

  # If no split is found, stop here
  if (is.null(split)) {
    return(partynode(id = id))
  }

  # Selected variable index and possibly selected basis index
  varid <- split$varid
  if (!is.null(split$basid)) {
    basid <- split$basid
  }

  breaks <- split$breaks
  index <- split$index

  # Assign the ids to the observations
  kidids <- c()
  switch(class(covariates$cov[[varid]]),

         integer = {

           kidids[(which(covariates$cov[[varid]] <= breaks))] <- 1
           kidids[(which(covariates$cov[[varid]] > breaks))] <- 2

         },

         numeric = {

           kidids[(which(covariates$cov[[varid]] <= breaks))] <- 1
           kidids[(which(covariates$cov[[varid]] > breaks))] <- 2

         },

         factor = {

           kidids <- index[covariates$cov[[varid]]]
           #replicate 1 in index for each level in splitpoint; 2 otherwise
           #no need for na.exclude()

         },

         fdata = {

           if (split_type == 'coeff') {

             kidids[which(covariates$newcov[[varid]][, basid] <= breaks)] <- 1
             kidids[which(covariates$newcov[[varid]][, basid] > breaks)] <- 2

           } else if (split_type == 'cluster') {

             kidids <- na.exclude(index)

           }
         },

         list = if (all(sapply(covariates$cov[[varid]],
                               function(x) attributes(x)$names) == 'diagram')) {

           kidids <- na.exclude(index)

         } else if (all(sapply(covariates$cov[[varid]], class) == 'igraph')) {

           if (split_type == 'coeff') {

             kidids[which(covariates$newcov[[varid]][as.character(basid)]
                          <= breaks)] <- 1
             kidids[which(covariates$newcov[[varid]][as.character(basid)]
                          > breaks)] <- 2
             #recall that basid is the column's name, so [as.character(basid)]

           } else if (split_type == 'cluster') {

             kidids <- na.exclude(index)

           }
         }
  )

  # Initialization of the kid nodes
  kids <- vector(mode = "list", length = max(kidids, na.rm = TRUE))

  # Give birth to the kid nodes
  for (kidid in seq_along(kids)) {

    # Select observations for the current kid node
    w <- weights
    w[kidids != kidid] <- 0

    # For less than <minbucket> observations, stop here
    if (sum(w) < minbucket) {
      return(partynode(id = id))
    }

    # Previous maximum id (to later set the id for the current kid node)
    if (kidid > 1) {
      prev_id <- max(nodeids(kids[[kidid - 1]]))
    } else {
      prev_id <- id
    }

    # Recursion on this kid node: update covariates and response
    covariates_updated <- list()
    covariates_updated$cov <- lapply(covariates$cov,
                                     function(cov) subset(cov, as.logical(w)))
    covariates_updated$newcov <- lapply(covariates$newcov, function(newcov)
      subset(newcov, as.logical(w)))
    covariates_updated$dist <- lapply(covariates$dist, function(cov_dist)
      usedist::dist_subset(cov_dist, which(w == 1)))
    response_updated <- list()
    response_updated$response <- subset(response$response, as.logical(w))
    response_updated$response_dist <-
      usedist::dist_subset(response$response_dist, which(w == 1))

    # Actual recursion
    kids[[kidid]] <- growtree(id = as.integer(prev_id + 1),
                              response = response_updated,
                              covariates = covariates_updated,
                              weights = rep(1L, sum(w, na.rm = TRUE)),
                              minbucket = minbucket,
                              alpha = alpha,
                              R = R,
                              split_type = split_type,
                              coeff_split_type = coeff_split_type,
                              p_adjust_method = p_adjust_method,
                              random_covs = random_covs)
  }

  # Return the nodes (i.e. the split rules)
  return(partynode(id = as.integer(id),
                   split = split,
                   kids = kids,
                   info = list(pvalue = min(info_split(split)$pvalue,
                                            na.rm = TRUE))
  ))
}



# Find the split --------------------------------------------------------------

findsplit <- function(response,
                      covariates,
                      alpha,
                      R,
                      split_type,
                      coeff_split_type,
                      p_adjust_method,
                      random_covs) {

  # Subset of covariates to consider for splitting (possibly all of them)
  cov_subset <- as.integer(seq_along(covariates$cov))
  if (!is.null(random_covs)) {
    cov_subset <- sample(cov_subset, random_covs)
  }

  # Independence test between the response and each covariate
  resp_dist <- response$response_dist
  stat_pval <- sapply(covariates$dist[cov_subset],
                      function(cov_dist) {
                        indep_test(x_dist = cov_dist, y_dist = resp_dist, R = R)
                      }
  )
  if (all(is.na(stat_pval['Pvalue', ]))) {
    return(NULL)
  }

  # Multiple testing correction
  adj_p <- p.adjust(stat_pval['Pvalue', ], method = p_adjust_method)

  # Stop criterion
  if (min(adj_p, na.rm = TRUE) > alpha) {
    return(NULL)
  }

  # Variable selection (based on original pvalues)
  xselect <- .select_variable(statistic_pvalue = stat_pval)
  #.select_variable is just a copy of .select_component

  # Selected covariates
  xselect <- cov_subset[xselect] #useful only if !is.null(random_covs)
  x <- covariates$cov[[xselect]]
  newx <- covariates$newcov[[xselect]]
  if (split_type == 'cluster') {
    xdist <- covariates$dist[[xselect]]
  }

  # Control for newx not being void
  if ((split_type == 'coeff') &&
      (class(x) == 'list') && all(sapply(x, class) == 'igraph') &&
      (dim(newx)[2] == 0)) {

    # Throw warning: if x is selected (through distance) and newx is void because
    # it had all-equal columns (after coeff exp), it means that either distance
    # or coeff exp are not appropriate, or at least that they are uncompatible
    warning('The selected covariate has non-informative coefficient expansion.
    Please reconsider the choice for at least one between distance and
    coefficient expansion. The selected covariate is ignored for the time
            being.')

    # Ignore the selected covariate and re-run findsplit
    # in order to keep the original indices (i.e. not to create confusion), the
    # selected covariate is not dropped, but it is instead replaced with a
    # vector of 0s, so that it is selected no more
    covariates$dist[[xselect]][] <- 0L
    split <- findsplit(response = response,
                       covariates = covariates,
                       alpha = alpha,
                       R = R,
                       split_type = split_type,
                       coeff_split_type = coeff_split_type,
                       p_adjust_method = p_adjust_method,
                       random_covs = random_covs)
    return(split)
  }

  # Split point search
  split_objs <- split_opt(resp = response,
                          cov = x,
                          new_cov = newx,
                          cov_dist = xdist,
                          R = R,
                          split_type = split_type,
                          coeff_split_type = coeff_split_type)

  # If split_objs is VOID, ignore the selected covariate and re-run findsplit
  if (isTRUE(split_objs$void)) {
    covariates$dist[[xselect]][] <- 0L
    split <- findsplit(response = response,
                       covariates = covariates,
                       alpha = alpha,
                       R = R,
                       split_type = split_type,
                       coeff_split_type = coeff_split_type,
                       p_adjust_method = p_adjust_method,
                       random_covs = random_covs)
    return(split)
  }

  # Separately save split_objs outputs
  splitpoint <- split_objs$splitpoint
  splitindex <- split_objs$splitindex
  bselect <- split_objs$bselect
  centroids <- split_objs$centroids

  # Return the split point
  switch(class(x),

         integer = {

           return(sp = partysplit(varid = as.integer(xselect),
                                  breaks = splitpoint,
                                  info = list(pvalue = stat_pval),
                                  right = TRUE))

         },

         numeric = {

           return(sp = partysplit(varid = as.integer(xselect),
                                  breaks = splitpoint,
                                  info = list(pvalue = stat_pval),
                                  right = TRUE))

         },

         factor = {

           return(sp = partysplit(varid = as.integer(xselect),
                                  index = splitindex,
                                  info = list(pvalue = stat_pval)))

         },

         fdata = {

           if (split_type == 'coeff') {
             return(sp = partysplit(varid = as.integer(xselect),
                                    basid = as.integer(bselect),
                                    breaks = splitpoint,
                                    right = TRUE,
                                    info = list(pvalue = stat_pval)))

           } else if (split_type == 'cluster') {

             sp <- partysplit(varid = as.integer(xselect),
                              centroids = centroids,
                              index = as.integer(splitindex),
                              info = list(pvalue = stat_pval))
             attr(sp, 'curr_split_type') <- 'cluster'   #used in edge.simple
             return(sp)

           }
         },

         list = if (all(sapply(x, function(x) attributes(x)$names)
                        == 'diagram')) {

           #only cluster
           sp <- partysplit(varid = as.integer(xselect),
                            centroids = centroids,
                            index = as.integer(splitindex),
                            info = list(pvalue = stat_pval))
           attr(sp, 'curr_split_type') <- 'cluster'   #used in edge.simple
           return(sp)

         } else if (all(sapply(x, class) == 'igraph')) {

           if (split_type == 'coeff') {

             return(sp = partysplit(varid = as.integer(xselect),
                                    basid = as.integer(bselect),
                                    breaks = splitpoint,
                                    right = TRUE,
                                    info = list(pvalue = stat_pval)))

           } else if (split_type == 'cluster') {

             sp <- partysplit(varid = as.integer(xselect),
                              centroids = centroids,
                              index = as.integer(splitindex),
                              info = list(pvalue = stat_pval))
             attr(sp, 'curr_split_type') <- 'cluster'   #used in edge.simple
             return(sp)

           }
         })
}



# Split point search ------------------------------------------------------


split_opt <- function(resp,
                      cov,
                      new_cov,
                      cov_dist,
                      R = 1000,
                      split_type = 'coeff',
                      coeff_split_type = 'test') {

  # Retrieve response_dist and response
  resp_dist <- resp$response_dist
  resp <- resp$response

  switch(class(cov),

         integer    = {

           s  <- sort(cov)
           comb <- sapply(s[-length(s)], function(j) cov <= j)

           if (coeff_split_type == 'traditional') {

             splitpoint <- .select_traditional(values = s,
                                               comb_logical = comb,
                                               resp = resp)

           } else if (coeff_split_type == 'test') {

             stat_pval <- apply(comb, 2,
                                function(q) indep_test(x = q,
                                                       y_dist = resp_dist,
                                                       R = R))
             splitpoint <- .select_splitpoint(values = s,
                                              statistic_pvalue = stat_pval)

           }
         },

         numeric    = {

           s  <- sort(cov)
           comb <- sapply(s[-length(s)], function(j) cov <= j)

           if (coeff_split_type == 'traditional') {

             splitpoint <- .select_traditional(values = s,
                                               comb_logical = comb,
                                               resp = resp)

           } else if (coeff_split_type == 'test') {

             stat_pval <- apply(comb, 2,
                                function(q) indep_test(x = q,
                                                       y_dist = resp_dist,
                                                       R = R))
             splitpoint <- .select_splitpoint(values = s,
                                              statistic_pvalue = stat_pval)

           }
         },

         factor     = {

           # Drop unused levels
           lev <- levels(cov[drop = TRUE])
           n_lev <- length(lev)

           if (n_lev == 2) {

             splitpoint <- lev[1]  #split point simply given by the first level

           } else {

             # Combination of all the levels
             lev_cmb <- do.call("c",
                                lapply(seq_len(n_lev / 2), function(ntaken) {
                                  #1 for n_lev=2,3; 1,2 for n_lev=4,5; 1,2,3 for
                                  #n_lev=6,7;... since [(1:x) := (1:floor(x))]
                                  combs_ntaken <- combn(x = lev,
                                                        m = ntaken,
                                                        simplify = FALSE)
                                  #fine, when n_lev is odd; for even n_lev, skip
                                  #half the combs for the last ntaken
                                  #e.g. n_lev=6 => last_ntaken=3, but '1,5,6' is
                                  #equivalent to '2,3,4'
                                  if (ntaken == (n_lev / 2)) {
                                    combs_ntaken <-
                                      combs_ntaken[seq_len(
                                        length(combs_ntaken) / 2)]
                                  }
                                  return(combs_ntaken)

                                }))
             comb <- sapply(lev_cmb, function(lc) cov %in% lc)

             if (coeff_split_type == 'traditional') {

               splitpoint <- .select_traditional(values = lev_cmb,
                                                 comb_logical = comb,
                                                 resp = resp)

             } else if (coeff_split_type == 'test') {

               stat_pval <- apply(comb, 2,
                                  function(q) indep_test(x = q,
                                                         y_dist = resp_dist,
                                                         R = R))
               splitpoint <- .select_splitpoint(values = lev_cmb,
                                                statistic_pvalue = stat_pval)

             }

           }

           # Label levels with 1 if they are in splitpoint, 2 otherwise
           # (and with NA if they do not occur)
           #needed in growtree to split observations using their level
           splitindex <- !(levels(cov) %in% splitpoint)
           splitindex[!(levels(cov) %in% lev)] <- NA_integer_
           splitindex <- splitindex - min(splitindex, na.rm = TRUE) + 1L

         },

         fdata      = {

           if (split_type == 'coeff') {

             comp_idx <- seq_len(dim(new_cov)[2])
             stat_pval <- sapply(comp_idx,
                                 function(i) indep_test(x = new_cov[, i],
                                                        y_dist = resp_dist,
                                                        R = R))
             bselect <- .select_component(statistic_pvalue = stat_pval)

             sel_comp <- new_cov[, bselect]
             s  <- sort(sel_comp)
             comb <- sapply(s[-length(s)], function(j) sel_comp <= j)

             if (coeff_split_type == 'traditional') {

               splitpoint <- .select_traditional(values = s,
                                                 comb_logical = comb,
                                                 resp = resp)

             } else if (coeff_split_type == 'test') {

               stat_pval <- apply(comb, 2,
                                  function(q) indep_test(x = q,
                                                         y_dist = resp_dist,
                                                         R = R))
               splitpoint <- .select_splitpoint(values = s,
                                                statistic_pvalue = stat_pval)

             }

           } else if (split_type == 'cluster') {

             clustering_objs <- .select_clustering(cov = cov,
                                                   new_cov = new_cov,
                                                   cov_dist = cov_dist)
             splitindex <- clustering_objs$splitindex
             centroids <- clustering_objs$centroids

           }

         },

         list = if (all(sapply(cov, function(obj) attributes(obj)$names)
                        == 'diagram')) {

           clustering_objs <- .select_clustering(cov = cov,
                                                 new_cov = new_cov,
                                                 cov_dist = cov_dist)
           splitindex <- clustering_objs$splitindex
           centroids <- clustering_objs$centroids

         } else if (all(sapply(cov, class) == 'igraph')) {

           if (split_type == 'coeff') {

             # Drop non-informative (i.e. all-equal) columns
             new_cov <- new_cov[, !as.logical(apply(new_cov, 2, .zero_range))]

             # Control if the df is now void; if so, return 'void'
             if (dim(new_cov)[2] == 0) {
               return(list('void' = TRUE))
             }

             comp_idx <- seq_len(dim(new_cov)[2])
             stat_pval <- sapply(comp_idx,
                                 function(i) indep_test(x = new_cov[, i],
                                                        y_dist = resp_dist,
                                                        R = R))
             bselect <- .select_component(statistic_pvalue = stat_pval)

             # graph_shell may drop columns, so switch from index to name
             bselect <- as.integer(names(new_cov)[bselect])

             sel_comp <- new_cov[[as.character(bselect)]]
             s  <- sort(sel_comp)
             comb <- sapply(s[-length(s)], function(j) sel_comp <= j)

             # Check if all columns of comb are equal
             if (all(apply(comb, 2, identical, comb[, 1]))) {
               # if TRUE, the omitted column [position length(s)] is different;
               # when this is the case, set last included column as splitpoint,
               # as it means that breaks is set before last column
               splitpoint <- s[(length(s) - 1)]
             } else if (coeff_split_type == 'traditional') {

               splitpoint <- .select_traditional(values = s,
                                                 comb_logical = comb,
                                                 resp = resp)

             } else if (coeff_split_type == 'test') {

               stat_pval <- apply(comb, 2,
                                  function(q) indep_test(x = q,
                                                         y_dist = resp_dist,
                                                         R = R))
               splitpoint <- .select_splitpoint(values = s,
                                                statistic_pvalue = stat_pval)

             }

           } else if (split_type == 'cluster') {

             clustering_objs <- .select_clustering(cov = cov,
                                                   new_cov = new_cov,
                                                   cov_dist = cov_dist)
             splitindex <- clustering_objs$splitindex
             centroids <- clustering_objs$centroids

           }
         }
  )

  split_out <- list()
  if (exists('splitpoint')) split_out$splitpoint <- splitpoint
  if (exists('splitindex')) split_out$splitindex <- splitindex
  if (exists('bselect')) split_out$bselect <- bselect
  if (exists('centroids')) split_out$centroids <- centroids
  return(split_out)

}


.select_splitpoint <- function(values, statistic_pvalue) {

  stopifnot(identical(rownames(statistic_pvalue), c('Statistic', 'Pvalue')))

  if (all(is.na(statistic_pvalue['Pvalue', ])) ||
      length(which(statistic_pvalue['Pvalue', ] ==
                   min(statistic_pvalue['Pvalue', ], na.rm = TRUE))) > 1) {

    splitpoint <- values[which.max(statistic_pvalue['Statistic', ])]

  } else {

    splitpoint <- values[which.min(statistic_pvalue['Pvalue', ])]

  }

  return(splitpoint)

}


.select_component <- function(statistic_pvalue) {

  stopifnot(identical(rownames(statistic_pvalue), c('Statistic', 'Pvalue')))

  if (all(is.na(statistic_pvalue['Pvalue', ])) ||
      length(which(statistic_pvalue['Pvalue', ] ==
                   min(statistic_pvalue['Pvalue', ], na.rm = TRUE))) > 1) {

    bselect <- as.integer(which.max(statistic_pvalue['Statistic', ]))

  } else {

    bselect <- as.integer(which.min(statistic_pvalue['Pvalue', ]))

  }

  return(bselect)

}
.select_variable <- .select_component


.select_clustering <- function(cov, new_cov, cov_dist) {

  stopifnot(identical(typeof(cov), 'list'))
  stopifnot(inherits(cov_dist, 'dist'))

  # Initialize splitindex
  N_obs <- length(levels(new_cov))
  splitindex <- rep(NA, N_obs)

  # Only two observations: split is trivial
  if (length(cov) == 2) {

    splitindex[new_cov] <- c(1, 2)

    centroids <- if (class(cov) == 'fdata') {
      list(c1 = cov[1, ], c2 = cov[2, ])
    } else {
      list(c1 = cov[[1]], c2 = cov[[2]])
    }

  } else {  # More than two observations: split via PAM

    pam_obj <- cluster::pam(cov_dist, k = 2, diss = TRUE)
    cl_index <- pam_obj$clustering
    splitindex[as.integer(names(cl_index))] <- cl_index

    medindex1 <- pam_obj$id.med[1]
    c1 <- if (class(cov) == 'fdata') cov[medindex1, ] else cov[[medindex1]]
    medindex2 <- pam_obj$id.med[2]
    c2 <- if (class(cov) == 'fdata') cov[medindex2, ] else cov[[medindex2]]
    centroids <- list(c1 = c1, c2 = c2)

  }

  return(list(splitindex = splitindex, centroids = centroids))

}


.select_traditional <- function(values, comb_logical, resp) {

  if (class(resp) == 'factor') {

    total_gini <- apply(comb_logical, 2,
                        function(c) {

                          y1 <- resp[c]
                          y2 <- resp[!c]
                          # Gini (alt. def.): G(B_r) = 1 - \sum_k p_{kr}^2
                          g1 <- 1 - sum(prop.table(table(y1)) ^ 2)
                          g2 <- 1 - sum(prop.table(table(y2)) ^ 2)

                          total_gini_c <- g1 + g2
                          return(total_gini_c)

                        })

    splitpoint <- values[which.min(total_gini)]

  } else {

    total_rss <- apply(comb_logical, 2,
                       function(c) {

                         y1 <- resp[c]
                         y2 <- resp[!c]
                         v1 <- if (length(y1) == 1) 0 else var(y1)
                         v2 <- if (length(y2) == 1) 0 else var(y2)
                         n1 <- length(y1)
                         n2 <- length(y2)

                         total_rss_c <- (n1 - 1) * v1 + (n2 - 1) * v2
                         return(total_rss_c)

                       })

    splitpoint <- values[which.min(total_rss)]

  }

  return(splitpoint)

}



# Independence test -----------------------------------------------------------

indep_test <- function(x,
                       y,
                       x_dist = NULL,
                       y_dist = NULL,
                       R = 1000) {

  # If distances are not provided, compute them
  if (is.null(x_dist)) x_dist <- dist_comp(x)
  if (is.null(y_dist)) y_dist <- dist_comp(y)

  # Distance correlation test
  dct <- energy::dcor.test(x_dist, y_dist, R = R)

  # Retrieve and return test statistic and p-value
  stat_pval <- if (!is.na(dct$statistic)) {
    c(dct$statistic, dct$p.value)
  } else{
    c(NA, NA)
  }
  names(stat_pval) <- c('Statistic', 'Pvalue')
  return(stat_pval)

}


# Distances ---------------------------------------------------------------

dist_comp <- function(x,
                      lp = 2) {

  # Compute the distance/dissimilarity objects
  dist_obj <- switch(class(x),

                logical    = dist(x),
                #needed for split point search when split_type = 'coeff'

                integer    = dist(x),
                #objects of class integer are not of class numeric

                numeric    = dist(x),

                factor     = cluster::daisy(as.data.frame(x)),

                fdata      = as.dist(fda.usc::metric.lp(x, lp = lp)),

                list       = {

                  # Graphs
                  if (all(sapply(x, class) == 'igraph')) {

                    # Unweighted
                    if (all(sapply(x, function(i)
                      is.null(igraph::edge.attributes(i)$weight)))) {
                      adj_data <- lapply(x, igraph::as_adjacency_matrix)
                      } else {  # Weighted
                      adj_data <- lapply(x, function(i) {
                        igraph::as_adjacency_matrix(i, attr = 'weight')
                      })
                    }

                    #d is obtained in the same way in the two cases:
                    d <- NetworkDistance::nd.edd(adj_data)
                    return(d$D)

                    # Persistence diagrams
                  } else if (all(sapply(x, function(x) attributes(x)$names)
                                 == 'diagram')) {
                    wass_fun <- function(i, j) TDA::wasserstein(x[[i]]$diagram,
                                                                x[[j]]$diagram)
                    vec_wass_fun <- Vectorize(wass_fun)
                    d_idx <- seq_along(x)
                    return(outer(d_idx, d_idx, vec_wass_fun))
                  }

                }
                )

  # Set observations' names (needed for eforest())
  dist_obj <- usedist::dist_setNames(dist_obj, seq_along(x))

  # Return
  return(dist_obj)

}



dist_comp_cl <- function(centroid,
                         x,
                         lp = 2) {

  switch(class(x),

         fdata  = fda.usc::metric.lp(fdata1 = x, fdata2 = centroid, lp = lp),

         list   = {

           # Graphs
           if (all(sapply(x, class) == 'igraph')) {

             # Unweighted
             if (all(sapply(x, function(i)
               is.null(igraph::edge.attributes(i)$weight)))) {
               adj_data <- lapply(x, igraph::as_adjacency_matrix)
               adj_centroid <- igraph::as_adjacency_matrix(centroid)
             } else {  # Weighted
               adj_data <- lapply(x, function(i) {
                 igraph::as_adjacency_matrix(i, attr = 'weight')
               })
               adj_centroid <- igraph::as_adjacency_matrix(centroid,
                                                           attr = 'weight')
             }

             #dist_centroid is obtained in the same way in the two cases:
             dist_centroid <- sapply(adj_data, function(i) {
               d <- NetworkDistance::nd.edd(list(i, adj_centroid))
               return(d$D)
               #since the distances are computed pairwise, d$D is numeric (i.e.,
               #not a 'dist' object anymore)
             })
             return(dist_centroid)

             # Persistence diagrams
           } else if (all(sapply(x, function(x) attributes(x)$names)
                          == 'diagram')) {
             wass_fun <- function(x, centroid)
               TDA::wasserstein(x$diagram, centroid$diagram)
             vec_wass_fun <- Vectorize(wass_fun, vectorize.args = 'x')
             return(vec_wass_fun(x, centroid))
           }

         }
  )
}




# Graphs ------------------------------------------------------------------

graph_shell <- function(graph_list,
                        max_shell = NULL,
                        predicting = FALSE) {

  # Number of observations (graphs)
  n_graphs <- length(graph_list)

  # Shell distribution for each graph
  table_shell <- lapply(graph_list, function(g) table(brainGraph::s_core(g)))

  # Observed maximum shell index
  obs_max_shell <- do.call(max, lapply(table_shell,
                                       function(s) as.integer(names(s))))

  # If max_shell is given, go for it; otherwise, set obs_max_shell as max_shell
  if (is.null(max_shell)) max_shell <- obs_max_shell

  # Column names for the shell df
  col_names <- as.character(seq(0, max_shell, 1))

  # Shell df inizialization
  all_shell_df <- data.frame(matrix(data = 0L,
                                    nrow = n_graphs,
                                    ncol = length(col_names)))
  colnames(all_shell_df) <- col_names

  # Fill in with the actual shell distibutions
  invisible(sapply(seq_len(n_graphs), function(i) {
    shells <- names(table_shell[[i]])
    cols <- intersect(col_names, shells)
    all_shell_df[i, cols] <<- table_shell[[i]][cols]
  }))
  # better a for cycle?
  # for(i in 1:n.graphs){
  #   cols <- names(table.shell[[i]])
  #   all.shell.df[i, cols] = table.shell[[i]][cols]
  # }

  # Ignore non-informative columns only when not in predict
  if (isFALSE(predicting)) {
    # Update all_shell_df ignoring non-informative columns
    all_shell_df <- all_shell_df[, !as.logical(apply(all_shell_df,
                                                     2,
                                                     .zero_range))]
  }

  # Return the final shell df
  return(all_shell_df)

}


# Function to test if all elements of a given vector are equal for tol provided
.zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  # if only one obs, equality cannot be tested -> return FALSE
  if (length(x) == 1) {
    return(FALSE)
  }
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}



# Selected features analysis --------------------------------------------------

# Determine the variables selected for splitting
sel.features_det <- function(object) {

  # Check that object has class party
  stopifnot(inherits(object, 'party'))

  # Find *internal* nodes' ids
  ids <- nodeids(object)[!nodeids(object) %in% nodeids(object, terminal = TRUE)]

  # Retrieve variables' id for each split
  varids <- unique(unlist(nodeapply(object, ids = ids,
                                    FUN = function(x) {
                                      varid_split(split_node(x))
                                    })))

  # Return variables' ids
  return(varids)
}


# Analysis
sel.features_analyze <- function(object_list) {

  # Determine selected features for each object
  lapply(object_list,
         FUN = function(object) {
           sel.features_det(object)
         })

  # Plot

  # Table?

}
