

# Main function -----------------------------------------------------------

etree <- function(response,
                  covariates,
                  case.weights = NULL,
                  minbucket = 5,
                  alpha = 0.05,
                  R = 1000,
                  split.type = 'coeff',
                  coef.split.type = 'test',
                  p.adjust.method = 'fdr',
                  nb = 5) {

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
        foo <- data.frame(t(fd3$coefs))
        names(foo) <- 1:length(names(foo))

      } else if(split.type == "cluster"){

        foo <- as.factor(1:length(response))

      }

      attr(foo, 'cov.type') <- 'fdata'
      return(foo)

    } else if(class(j) == 'list' & all(sapply(j, function(x) attributes(x)$names) == 'diagram')){

      foo <- as.factor(1:length(response))
      attr(foo, 'cov.type') <- 'diagram'
      return(foo)

    } else if(class(j) == 'list' &
              all(sapply(j, class) == 'igraph')){

      if(split.type == "coeff"){

        foo <- graph.shell(j)

      } else if(split.type == "cluster"){

        foo <- as.factor(1:length(response))

      }

      attr(foo, 'cov.type') <- 'graph'
      return(foo)

    } else {

      return(j)

    }
  })

  # Covariates name
  if(!is.null(names(covariates))){
    names(newcovariates) <- names(covariates)
    #control if any name is void, i.e. if it is ''
    no_name <- which(sapply(names(covariates), function(n) n == '', USE.NAMES = FALSE))
    names(newcovariates) <- replace(names(newcovariates),
                                    no_name,
                                    as.factor(1:length(no_name)))
  } else {
    warning('No names available for covariates. Numbers are used instead.')
    names(newcovariates) <- 1:length(newcovariates)
  }

  # Distances
  cov.distance <- lapply(covariates, compute.dissimilarity)

  # Large list with covariates, newcovariates and distances
  covariates.large = list('cov' = covariates, 'newcov' = newcovariates, 'dist' = cov.distance)

  # Grow the tree (finds the split rules)
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
                    p.adjust.method = p.adjust.method,
                    nb = nb)
  print(c('NODES', nodes))

  # Actually perform the splits
  fitted.obs <- fitted_node(nodes, data = newcovariates)

  # Return a rich constparty object
  obj <- party(nodes,
               data = newcovariates,
               fitted = data.frame("(fitted)" = fitted.obs,
                                   "(response)" = response,
                                   check.names = FALSE),
               terms = terms(response ~ ., data = newcovariates))
  etree_obj <- as.constparty(obj)
  attr(etree_obj, 'split.type') <- split.type     #used in predict.party

  return(etree_obj)

}



# Grow the tree ---------------------------------------------------------------

growtree <- function(id = 1L,
                     response,
                     covariates,
                     case.weights,
                     minbucket,
                     alpha,
                     R,
                     n.var,
                     split.type,
                     coef.split.type,
                     p.adjust.method,
                     nb) {

  # Find the best split (variable selection & split point search)
  split <- findsplit(response = response,
                     covariates = covariates,
                     alpha = alpha,
                     R = R,
                     lp = rep(2, 2),
                     split.type = split.type,
                     coef.split.type = coef.split.type,
                     p.adjust.method = p.adjust.method,
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

  # Assign the ids to the observations
  kidids <- c()
  switch(class(covariates$cov[[varid]]),

         integer = {

           kidids[(which(covariates$newcov[[varid]] <= breaks))] <- 1
           kidids[(which(covariates$newcov[[varid]] > breaks))] <- 2

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

           if(split.type == 'coeff'){

             # observations before the split point are assigned to node 1
             kidids[which(covariates$newcov[[varid]][, basid] <= breaks)] <- 1
             #  observations before the split point are assigned to node 2
             kidids[which(covariates$newcov[[varid]][, basid] > breaks)] <- 2

           } else if (split.type == 'cluster') {

             kidids <- na.exclude(index)

           }
         },

         list = if(all(sapply(covariates$cov[[varid]], function(x) attributes(x)$names) == 'diagram')){

           kidids <- na.exclude(index)

         } else if(all(sapply(covariates$cov[[varid]], class) == 'igraph')){

           if(split.type == 'coeff'){

             kidids[which(covariates$newcov[[varid]][as.character(basid)]
                          <= breaks)] <- 1
             kidids[which(covariates$newcov[[varid]][as.character(basid)]
                          > breaks)] <- 2
             #recall that basid is the column's name, so [as.character(basid)]

           } else if(split.type == 'cluster') {

             kidids <- na.exclude(index)

           }
         }
  )

  # Initialization of the kid nodes
  kids <- vector(mode = "list", length = max(kidids, na.rm = TRUE))

  # Give birth to the kid nodes
  for (kidid in 1:length(kids)) {

    # Select observations for the current kid node
    w <- case.weights
    w[kidids != kidid] <- 0

    # For less than <minbucket> observations, stop here
    if (sum(w) < minbucket)
      return(partynode(id = id))

    # Previous maximum id (to later set the id for the current kid node)
    if (kidid > 1) {
      prev_id <- max(nodeids(kids[[kidid - 1]]))
    } else{
      prev_id <- id
    }

    # Start recursion on this kid node: initialization
    covariates.updated <- list()
    covariates.updated$cov <- lapply(covariates$cov, function(cov) subset(cov, as.logical(w)))
    covariates.updated$newcov <- lapply(covariates$newcov, function(cov) subset(cov, as.logical(w)))
    covariates.updated$dist <- lapply(covariates$dist, function(cov) subset(cov, subset = as.logical(w), select = which(w == 1)))

    # Recursion
    kids[[kidid]] <-
      growtree(
        id = as.integer(prev_id + 1),
        response = subset(response, as.logical(w)),
        covariates = covariates.updated,
        case.weights = rep(1L, sum(w, na.rm = TRUE)),
        minbucket = minbucket,
        alpha = alpha,
        R = R,
        n.var = n.var,
        split.type = split.type,
        coef.split.type = coef.split.type,
        p.adjust.method = p.adjust.method,
        nb = nb)
  }

  # Return the nodes (i.e. the split rules)
  return(partynode(id = as.integer(id),
                   split = split,
                   kids = kids,
                   info = list(p.value = min(info_split(split)$p.value, na.rm = TRUE))
  ))
}



# Find the split --------------------------------------------------------------

findsplit <- function(response,
                      covariates,
                      alpha,
                      R,
                      lp,
                      split.type,
                      coef.split.type,
                      p.adjust.method,
                      nb) {

  # Number of original covariates
  n.cov = length(covariates$cov)

  print('one round again')

  # Independence test between the response and each covariate
  p = lapply(covariates$dist,
             function(cov.dist) {
               ct <- energy::dcor.test(cov.dist,
                                       compute.dissimilarity(response),
                                       R = R)
               if (!is.na(ct$statistic)) {
                 return(c(ct$statistic, ct$p.value))
               } else{
                 c(NA, NA)
               }
             }
  )

  # Create matrix with test stats and p-values; check if p-values are all NULL
  p = t(matrix(unlist(p), ncol = 2, byrow = T))
  rownames(p) <- c("statistic", "p-value")
  if (all(is.na(p[2,]))) return(NULL)

  # Multiple testing correction
  adj_p <- p.adjust(p[2,], method = p.adjust.method)

  # Stop criterion
  if (min(adj_p, na.rm = TRUE) > alpha) return(NULL)

  # Variable selection (based on original p-values)
  if (length(which(p[2,] == min(p[2,], na.rm = TRUE))) > 1) {
    xselect <- which.max(p[1,])
    #in case of multiple minima, take that with the highest test statistic
  } else{
    xselect <- which.min(p[2,])
  }

  # Selected covariates
  x <- covariates$cov[[xselect]]
  newx <- covariates$newcov[[xselect]]
  if(split.type == 'cluster'){
    xdist <- covariates$dist[[xselect]]
  }

  # Control for newx not being void
  if((split.type == 'coeff') &&
     (class(x) == 'list') && all(sapply(x, class) == 'igraph') &&
     (dim(newx)[2] == 0)){

    # Throw warning: if x is selected (through distance) and newx is void because
    # it had all-equal columns (after coeff exp), it means that either distance
    # or coeff exp are not appropriate, or at least that they are uncompatible
    warning('The selected covariate has non-informative coefficient expansion.
            Please reconsider the choice for at least one between distance and
            coefficient expansion. The selected covariate is ignored for the
            time being.')

    # Ignore the selected covariate and re-run findsplit
    # in order to keep the original indices (i.e. not to create confusion), the
    # selected covariate is not dropped, but it is instead replaced with a
    # vector of 0s, so that it is selected no more
    covariates$dist[[xselect]][] <- 0L
    split <- findsplit(response = response,
                       covariates = covariates,
                       alpha = alpha,
                       R = R,
                       lp = rep(2, 2),
                       split.type = split.type,
                       coef.split.type = coef.split.type,
                       p.adjust.method = p.adjust.method,
                       nb = nb)
    return(split)
  }

  # Split point search
  split.objs = split.opt(y = response,
                         x = x,
                         newx = newx,
                         xdist = xdist,
                         split.type = split.type,
                         coef.split.type = coef.split.type,
                         nb = nb)

  # If split.objs is VOID, ignore the selected covariate and re-run findsplit
  if(isTRUE(split.objs$void)){
    covariates$dist[[xselect]][] <- 0L
    split <- findsplit(response = response,
                       covariates = covariates,
                       alpha = alpha,
                       R = R,
                       lp = rep(2, 2),
                       split.type = split.type,
                       coef.split.type = coef.split.type,
                       p.adjust.method = p.adjust.method,
                       nb = nb)
    return(split)
  }

  # Separately save split.objs outputs
  splitindex <- split.objs$splitindex
  bselect <- split.objs$bselect
  centroids <- split.objs$centroids

  # Return the split point
  switch(class(x),

         integer = {

           return(sp = partysplit(varid = as.integer(xselect),
                                  breaks = splitindex,
                                  info = list(p.value = p),
                                  right = TRUE))

         },

         numeric = {

           return(sp = partysplit(varid = as.integer(xselect),
                                  breaks = splitindex,
                                  info = list(p.value = p),
                                  right = TRUE))

         },

         factor = {

           return(sp = partysplit(varid = as.integer(xselect),
                                  index = splitindex,
                                  info = list(p.value = p)))

         },

         fdata = {

           if(split.type == 'coeff'){
             return(sp = partysplit(varid = as.integer(xselect),
                                    basid = as.integer(bselect),
                                    breaks = splitindex,
                                    right = TRUE,
                                    info = list(p.value = p)))

           } else if(split.type == 'cluster'){

             sp = partysplit(varid = as.integer(xselect),
                             centroids = centroids,
                             index = as.integer(splitindex),
                             info = list(p.value = p))
             attr(sp, 'this.split.type') <- 'cluster'   #used in edge.simple
             return(sp)

           }
         },

         list = if(all(sapply(x, function(x) attributes(x)$names) == 'diagram')){

           #only cluster
           sp = partysplit(varid = as.integer(xselect),
                           centroids = centroids,
                           index = as.integer(splitindex),
                           info = list(p.value = p))
           attr(sp, 'this.split.type') <- 'cluster'   #used in edge.simple
           return(sp)

         } else if(all(sapply(x, class) == 'igraph')){

           if(split.type == 'coeff'){

             return(sp = partysplit(varid = as.integer(xselect),
                                    basid = as.integer(bselect),
                                    breaks = splitindex,
                                    right = TRUE,
                                    info = list(p.value = p)))

           } else if(split.type == 'cluster') {

             sp = partysplit(varid = as.integer(xselect),
                             centroids = centroids,
                             index = as.integer(splitindex),
                             info = list(p.value = p))
             attr(sp, 'this.split.type') <- 'cluster'   #used in edge.simple
             return(sp)

           }
         })
}



# Split point search ------------------------------------------------------


split.opt <- function(y,
                      x,
                      newx,
                      xdist,
                      split.type = 'coeff',
                      coef.split.type = 'test',
                      nb,
                      R = 500,
                      wass.dist = NULL){

  switch(class(x),

         integer    = {

           s  <- sort(x)
           comb = sapply(s[-length(s)], function(j) x <= j)
           xp.value <- apply(comb, 2, function(q) independence.test(x = q, y = y))
           if (length(which(xp.value[2,] ==
                            min(xp.value[2,], na.rm = T))) > 1 ||
               all(is.na(xp.value[2,]))) {
             splitindex <- s[which.max(xp.value[1,])]
           } else {
             splitindex <- s[which.min(xp.value[2,])]
           }

         },

         numeric    = {

           s  <- sort(x)
           comb = sapply(s[-length(s)], function(j) x <= j)
           #first one is excluded since it only return FALSEs
           xp.value <- apply(comb, 2, function(q) independence.test(x = q, y = y))
           if (length(which(xp.value[2,] ==
                            min(xp.value[2,], na.rm = T))) > 1 ||
               all(is.na(xp.value[2,]))) {
             splitindex <- s[which.max(xp.value[1,])]
           } else {
             splitindex <- s[which.min(xp.value[2,])]
           }

         },

         factor     = {

           # Drop unused levels
           lev <- levels(x[drop = TRUE])

           if (length(lev) == 2) {
             splitpoint <- lev[1]
             #the split point is simply given by the first level
           } else{
             # Combination of all the levels
             comb <- do.call("c",
                             lapply(1:(length(lev) - 2),
                                    function(ntaken) combn(x = lev,
                                                           m = ntaken,
                                                           simplify = FALSE)))
             #todo: take only first length(lev)/2 - 1, the other are complements!
             xp.value <- sapply(comb,
                                function(q) independence.test(x %in% q, y))

             if (length(which(xp.value[2,] ==
                              min(xp.value[2,], na.rm = T))) > 1 ||
                 all(is.na(xp.value[2,]))) {
               splitpoint <- comb[[which.max(xp.value[1,])]]
             } else {
               splitpoint <- comb[[which.min(xp.value[2,])]]
             }
           }

           # Label levels with 1 if they are in splitpoint, 2 otherwise
           # (and with NA if they do not occur)
           #needed in growtree to split observations using their level
           splitindex <- !(levels(x) %in% splitpoint)
           splitindex[!(levels(x) %in% lev)] <- NA_integer_
           splitindex <- splitindex - min(splitindex, na.rm = TRUE) + 1L

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
             comb = sapply(s[-length(s)], function(j) sel.coeff <= j)

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
               if (length(which(xp.value[2,] ==
                                min(xp.value[2,], na.rm = T))) > 1 ||
                   all(is.na(xp.value[2,]))) {
                 splitindex <- s[which.max(xp.value[1,])]
               } else {
                 splitindex <- s[which.min(xp.value[2,])]
               }

             }

           } else if(split.type == 'cluster') {

             cl.fdata <- cluster::pam(xdist, k = 2, diss = TRUE)
             clindex <- cl.fdata$clustering
             lev = levels(newx)
             splitindex = rep(NA, length(lev))
             splitindex[lev %in% newx[clindex==1]]<- 1
             splitindex[lev %in% newx[clindex==2]]<- 2

             medindex1 <- cl.fdata$id.med[1]
             c1 <- x[medindex1,]
             medindex2 <- cl.fdata$id.med[2]
             c2 <- x[medindex2,]
             centroids <- list(c1 = c1, c2 = c2)

           }

         },

         list = if(all(sapply(x, function(x) attributes(x)$names) == 'diagram')){

           cl.diag <- cluster::pam(xdist, k = 2, diss = TRUE)
           clindex <- cl.diag$clustering
           lev = levels(newx)
           splitindex = rep(NA, length(lev))
           splitindex[lev %in% newx[clindex==1]]<- 1
           splitindex[lev %in% newx[clindex==2]]<- 2

           medindex1 <- cl.diag$id.med[1]
           c1 <- x[[medindex1]]
           medindex2 <- cl.diag$id.med[2]
           c2 <- x[[medindex2]]
           centroids <- list(c1 = c1, c2 = c2)


         } else if(all(sapply(x, class) == 'igraph')){

           if(split.type == 'coeff'){
             x1 = newx
             # Drop non-informative (i.e. all-equal) columns
             x1 <- x1[, !as.logical(apply(x1, 2, zero_range))]
             # Control if the df is now void; if so, return 'void'
             if(dim(x1)[2] == 0) return(list('void' = TRUE))
             bselect <- 1:dim(x1)[2]
             p1 <- c()
             p1 <- sapply(bselect, function(i) independence.test(x1[, i], y, R = R))
             colnames(p1) <- colnames(x1)
             if (length(which(p1[2,] ==
                              min(p1[2,], na.rm = T))) > 1 ||
                 all(is.na(p1[2,]))) {
               bselect <- as.integer(which.max(p1[1,]))
             } else{
               bselect <- as.integer(which.min(p1[2,]))
             }
             # graph.shell may drop columns, so switch from index to name
             bselect <- as.integer(names(x1)[bselect])
             sel.coeff = x1[[as.character(bselect)]]
             s  <- sort(sel.coeff)
             comb = sapply(s[-length(s)], function(j) sel.coeff <= j)

             # Check if all columns of comb are equal
             if(all(apply(comb, 2, identical, comb[,1]))){
               # if TRUE, the omitted column [position length(s)] is different;
               # when this is the case, set last included column as splitindex,
               # as it means that breaks is set before last column
               splitindex <- s[(length(s) - 1)]
             } else if(coef.split.type == 'variance'){

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
               if (length(which(xp.value[2,] ==
                                min(xp.value[2,], na.rm = T))) > 1 ||
                   all(is.na(xp.value[2,]))) {
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

             medindex1 <- cl.graph$id.med[1]
             c1 <- x[[medindex1]]
             medindex2 <- cl.graph$id.med[2]
             c2 <- x[[medindex2]]
             centroids <- list(c1 = c1, c2 = c2)
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

  # Compute the dissimilarities within x and y
  d1 = compute.dissimilarity(x, lp = lp[1])
  d2 = compute.dissimilarity(y, lp = lp[2])

  # Distance correlation test
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
         integer    = as.matrix(dist(x)),
         numeric    = as.matrix(dist(x)),
         factor     = as.matrix(cluster::daisy(as.data.frame(x))),
         fdata      = metric.lp(x, lp=lp),
         list       = {
           if(all(sapply(x, class) == 'igraph')){
             if(all(sapply(x,
                           function(i) {
                             is.null(edge.attributes(i)$weight)
                             #if attribute weight is null for all the graphs,
                             #the graph covariate is not weighted
                             })
                    )){
               adj_data <- lapply(x, igraph::as_adjacency_matrix)
             } else { #otherwise, it is weighted
               adj_data <- lapply(x, function(i) {
                 igraph::as_adjacency_matrix(i, attr = 'weight')
               })
             }
             #d is obtained in the same way in the two cases:
             d <- NetworkDistance::nd.edd(adj_data)
             return(as.matrix(d$D))
           } else if(all(sapply(x, function(x) attributes(x)$names) == 'diagram')){
             k.fun = function(i,j) TDA::wasserstein(x[[i]]$diagram, x[[j]]$diagram)
             k.fun = Vectorize(k.fun)
             d.idx = seq_along(x)
             return(outer(d.idx,d.idx, k.fun))
           }
         })

}




compute.dissimilarity.cl <- function(centroid, x,
                                     lp = 2){

  switch(class(x),
         fdata      = metric.lp(fdata1 = x, fdata2 = centroid, lp=lp),
         list       = {
           if(all(sapply(x, class) == 'igraph')){
             if(all(sapply(x, function(i) {
               is.null(edge.attributes(i)$weight)
               #if attribute weight is null for all the graphs, the graph
               #covariate is not weighted
             }))){
               adj_data <- lapply(x, igraph::as_adjacency_matrix)
               adj_centroid <- igraph::as_adjacency_matrix(centroid)
             } else { #otherwise, it is weighted
               adj_data <- lapply(x, function(i) {
                 igraph::as_adjacency_matrix(i, attr = 'weight')
               })
               adj_centroid <- igraph::as_adjacency_matrix(centroid, attr = 'weight')
             }
             #dist_centroid is obtained in the same way in the two cases:
             dist_centroid <- sapply(adj_data, function(i){
               d <- NetworkDistance::nd.centrality(list(i, adj_centroid), mode = 'Degree', directed = TRUE)
               d$D
             })
             return(dist_centroid)
           } else if (all(sapply(x, function(x) attributes(x)$names) == 'diagram')){
             k.fun = function(x, centroid) TDA::wasserstein(x$diagram, centroid$diagram)
             k.fun = Vectorize(k.fun, vectorize.args = 'x')
             return(k.fun(x, centroid))
           }
         })


}




# Graphs ------------------------------------------------------------------

graph.shell <- function(graph.list, max.shell = NULL, predicting = FALSE){

  # Number of observations (graphs)
  n.graphs <- length(graph.list)

  # Shell distribution for each graph
  table.shell <- lapply(graph.list, function(g){table(brainGraph::s_core(g))})

  # Observed maximum shell index
  obs.max.shell <- do.call(max, lapply(table.shell,
                                       function(s){
                                         as.integer(names(s))
                                       }))

  # If max.shell is provided, go for it; otherwise, set obs.max.shell as max.shell
  if(is.null(max.shell)) max.shell <- obs.max.shell

  # Column names for the shell df
  col.names = as.character(seq(0, max.shell, 1))

  # Shell df inizialization
  all.shell.df = data.frame(matrix(
    data = 0L,
    nrow = n.graphs,
    ncol = length(col.names)))
  colnames(all.shell.df) <- col.names

  # Fill in with the actual shell distibutions
  invisible(sapply(1:n.graphs, function(i){
    shells <- names(table.shell[[i]])
    cols <- intersect(col.names, shells)
    all.shell.df[i, cols] <<- table.shell[[i]][cols]
    # <<- for global environment assignment
  }))
  # better a for cycle?
  # for(i in 1:n.graphs){
  #   cols <- names(table.shell[[i]])
  #   all.shell.df[i, cols] = table.shell[[i]][cols]
  # }

  # Ignore non-informative columns only when not in predict
  if(isFALSE(predicting)){
    # Update all.shell.df ignoring non-informative columns
    all.shell.df <- all.shell.df[, !as.logical(apply(all.shell.df,
                                                     2,
                                                     zero_range))]
  }

  # Return the final shell df
  return(all.shell.df)

}


# Function to test if all elements of a given vector are equal for tol provided
zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  # if only one obs, equality cannot be tested -> return FALSE
  if(length(x) == 1) return(FALSE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}



# Selected features analysis --------------------------------------------------

# Determine the variables selected for splitting
sel.features_det <- function(object){

  # Check that object has class party
  stopifnot(inherits(object, 'party'))

  # Find *internal* nodes' ids
  ids <- nodeids(object)[!nodeids(object) %in% nodeids(object, terminal = TRUE)]

  # Retrieve variables' id for each split
  varids <- unique(unlist(nodeapply(object, ids = ids,
                                    FUN = function(x){
                                      varid_split(split_node(x))
                                    })))

  # Return variables' ids
  return(varids)
}


# Analysis
sel.features_analyze <- function(object_list){

  # Determine selected features for each object
  lapply(object_list,
         FUN = function(object){
           sel.features_det(object)
         })

  # Plot

  # Table?

}



