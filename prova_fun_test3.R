
# Loading the libraries
library(cluster)
library(fda.usc)



# Main function -----------------------------------------------------------

etree <- function(response,
                  covariates,
                  case.weights = NULL,
                  minbucket = 1,
                  alpha = 0.05,
                  R = 1000,
                  split.type = 'coeff',
                  coef.split.type = 'test',
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

  # Building a df with all the new 'variables'
  newcovariates.df <- as.data.frame(newcovariates)
  names(newcovariates.df) <- 1:ncol(newcovariates.df)

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
  fitted.obs <- fitted_node(nodes, data = newcovariates.df)

  # Returning a rich constparty object
  ret <- party(nodes,
               data = newcovariates.df,
               fitted = data.frame("(fitted)" = fitted.obs,
                                   "(response)" = response,
                                   check.names = FALSE),
               terms = terms(response ~ ., data = data.frame(response = response, newcovariates.df)))

  return(etree = as.constparty(ret))

}



# growtree ----------------------------------------------------------------

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
  res_splt <- findsplit(response = response,
                        covariates = covariates,
                        alpha = alpha,
                        R = R,
                        lp = rep(2, 2),
                        split.type = split.type,
                        coef.split.type = coef.split.type,
                        nb = nb)

  # Separately saving res_splt outputs
  sp <- res_splt$sp
  varselect <- res_splt$varselect

  # If no split is found, stop here
  if (is.null(sp))
    return(partynode(id = id))

  # Assigning the ids to the observations
  kidids <- c()
  switch(class(covariates$cov[[varselect]]),

         fdata = {

           if(split.type == 'coeff'){

             # observations before the split point are assigned to node 1
             kidids[which(covariates$newcov[[varselect]][, sp$varid] <= sp$breaks)] <- 1
             #  observations before the split point are assigned to node 2
             kidids[which(covariates$newcov[[varselect]][, sp$varid] > sp$breaks)] <- 2

           } else if (split.type == 'cluster') {

             kidids <- na.exclude(sp$index)

           }
         },

         numeric = {

           kidids[(which(covariates$cov[[varselect]] <= sp$breaks))] <- 1
           kidids[(which(covariates$cov[[varselect]] > sp$breaks))] <- 2

         },

         integer = {

           kidids[(which(covariates$newcov[[varselect]] <= sp$breaks))] <- 1
           kidids[(which(covariates$newcov[[varselect]] > sp$breaks))] <- 2

         },

         factor = {

           kidids <- na.exclude(sp$index)

         },

         list = if(FALSE){
           #attributes(x[[1]])$names == 'diagram'
         } else if(all(sapply(covariates$cov[[varselect]], class) == 'igraph')){

           if(split.type == 'coeff'){

             kidids[which(covariates$newcov[[varselect]][, sp$varid] <= sp$breaks)] <- 1
             kidids[which(covariates$newcov[[varselect]][, sp$varid] > sp$breaks)] <- 2

           } else if(split.type == 'cluster') {

             kidids <- na.exclude(sp$index)

           }
         }
  )

  # Total number of features for each covariate
  total_features <- lapply(covariates$cov,
                           function(v) {
                             switch(
                               class(v),
                               logical    = 1,
                               factor     = 1,
                               numeric    = 1,
                               integer    = 1,
                               matrix     = ncol(v),
                               fdata      = {
                                 if(split.type == "coeff"){
                                   foo <- fda.usc::optim.basis(v, numbasis = nb)
                                   foo$numbasis.opt
                                 } else if(split.type == "cluster"){
                                   1
                                 }
                               },
                               list       = if(FALSE){
                                 #attributes(v[[1]])$names == 'diagram'
                               } else if(all(sapply(v, class) == 'igraph')){
                                 if(split.type == 'coeff'){
                                   idx <- which(sapply(covariates$cov, function(c){identical(c,v)}))
                                   ncol(covariates$newcov[[idx]])
                                 } else if(split.type == 'cluster') {
                                   1
                                 }
                               }
                             )
                           })

  # Total number of features before the selected variable
  step <- if(length(total_features) > 1) {
    sum_feat <- do.call(sum, total_features[which(1:n.var < varselect)])
    as.integer(sum_feat)
  } else {
    sum_feat <- 0L
  }

  # Shifting the varid by the number of the previous features
  if(class(covariates$cov[[varselect]]) == 'fdata' &&
     split.type == "coeff"){
    sp$varid = step + sp$varid #since here sp$varid is bselect
  } else if(class(covariates$cov[[varselect]]) == 'list' &&
            all(sapply(covariates$cov[[varselect]], class) == 'igraph') &&
            split.type == "coeff") {
    sp$varid = step + sp$varid #since here sp$varid is bselect
  } else {
    sp$varid = step + 1 #since sp$varid is xselect!
  }

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
                   split = sp,
                   kids = kids,
                   info = list(p.value = min(info_split(sp)$p.value, na.rm = TRUE))
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
  centers <- split.objs$centers

  # Returning the split point
  switch(class(x),

         numeric = {

           return(list(sp = partysplit(varid = as.integer(xselect),
                                       breaks = splitindex,
                                       info = list(p.value = 1-(1-p)^sum(!is.na(p)))),
                       varselect = xselect))

         },

         integer = {

           return(list(sp = partysplit(varid = as.integer(xselect),
                                       breaks = splitindex,
                                       info = list(p.value = 1-(1-p)^sum(!is.na(p)))),
                       varselect = xselect))

         },

         factor = {

           return(list(sp = partysplit(varid = as.integer(xselect),
                                       index = splitindex,
                                       info = list(p.value = 1-(1-p)^sum(!is.na(p)))),
                       varselect = xselect))

         },

         fdata = {

           if(split.type == 'coeff'){
             return(list(sp = partysplit(varid = as.integer(bselect),
                                         breaks = splitindex,
                                         info = list(p.value = 1-(1-p[2,])^sum(!is.na(p[2,])))),
                         varselect = xselect))
           } else if(split.type == 'cluster'){
             return(list(sp = partysplit(varid = as.integer(xselect),
                                         index = as.integer(splitindex),
                                         info = list(p.value = 1-(1-p[2,])^sum(!is.na(p[2,])),
                                                     centers = centers)),
                         varselect = xselect))
           }

         },

         list = if(FALSE){
           #attributes(v[[1]])$names == 'diagram'
           return(list(sp = partysplit(varid = as.integer(xselect),
                                       index = as.integer(splitindex),
                                       info = list(p.value = 1-(1-p[2,])^sum(!is.na(p[2,])))),
                       varselect = xselect))

         } else if(all(sapply(x, class) == 'igraph')){

           if(split.type == 'coeff'){

             return(list(sp = partysplit(varid = as.integer(bselect),
                                         breaks = splitindex,
                                         info = list(p.value = 1-(1-p[2,])^sum(!is.na(p[2,])))),
                         varselect = xselect))

           } else if(split.type == 'cluster') {

             return(list(sp = partysplit(varid = as.integer(xselect),
                                         index = as.integer(splitindex),
                                         info = list(p.value = 1-(1-p[2,])^sum(!is.na(p[2,])),
                                                     centers = centers)),
                         varselect = xselect))

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
             cl.fdata = kmeans.fd(x, ncl=2, draw = FALSE, par.ini=list(method="exact"))
             clindex <- cl.fdata$cluster
             lev = levels(newx)
             splitindex = rep(NA, length(lev))
             splitindex[lev %in% newx[clindex==1]]<- 1
             splitindex[lev %in% newx[clindex==2]]<- 2

             c1 <- cl.fdata$centers[1]
             c2 <- cl.fdata$centers[2]
             centers <- list(c1 = c1, c2 = c2)
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
             c1 <- x[[ceindex1]]
             ceindex2 <- as.integer(cl.graph$medoids[2])
             c2 <- x[[ceindex2]]
             centers <- list(c1 = c1, c2 = c2)
           }


         }
  )

  out <- list('splitindex' = splitindex)
  if(exists('bselect')) out$bselect <- bselect
  if(exists('centers')) out$centers <- centers
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
