

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

        foo <- as.factor(paste0('cluster',1:length(response)))

      }

      return(foo)

    } else if(class(j) == 'list' &
              all(sapply(j, class) == 'igraph')){

      if(split.type == "coeff"){
        foo <- graph.shell(j)
      } else if(split.type == "cluster"){
        foo <- as.factor(paste0('cluster',1:length(response)))
      }

      return(foo)

    } else if(class(j) == 'list' & all(sapply(j, function(x) attributes(x)$names) == 'diagram')){

      foo <- as.factor(paste0('cluster',1:length(response)))

      return(foo)

    }
    else {

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

         list = if(all(sapply(covariates$cov[[varid]], function(x) attributes(x)$names) == 'diagram')){

           kidids <- na.exclude(index)

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

         list = if(all(sapply(x, function(x) attributes(x)$names) == 'diagram')){

           return(sp = partysplit(varid = as.integer(xselect),
                                  centroids = centroids,
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
           comb = sapply(s[2:length(s)], function(j) x<j)
           #first one is excluded since it only return FALSEs
           xp.value <- apply(comb, 2, function(q) independence.test(x = q, y = y))
           if (length(which(xp.value[2,] == min(xp.value[2,], na.rm = T))) > 1) {
             splitindex <- s[which.max(xp.value[1,])]
           } else {
             splitindex <- s[which.min(xp.value[2,])]
           }

         },

         integer    = {

           s  <- sort(x)
           comb = sapply(s[2:length(s)], function(j) x<j)
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

             cl.fdata <- cluster::pam(xdist, k = 2, diss = TRUE)
             clindex <- cl.fdata$clustering
             lev = levels(newx)
             splitindex = rep(NA, length(lev))
             splitindex[lev %in% newx[clindex==1]]<- 1
             splitindex[lev %in% newx[clindex==2]]<- 2

             ceindex1 <- as.integer(cl.fdata$medoids[1])
             c1 <- x[ceindex1,]
             ceindex2 <- as.integer(cl.fdata$medoids[2])
             c2 <- x[ceindex2,]
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

           ceindex1 <- as.integer(cl.diag$medoids[1])
           c1 <- x[[ceindex1]]
           ceindex2 <- as.integer(cl.diag$medoids[2])
           c2 <- x[[ceindex2]]
           centroids <- list(c1 = c1, c2 = c2)


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
             c1 <- x[[which(as.integer(sub('cluster','',newx)) == ceindex1)]]
             ceindex2 <- as.integer(cl.graph$medoids[2])
             c2 <- x[[which(as.integer(sub('cluster','',newx)) == ceindex2)]]
             centroids <- list(c1 = c1, c2 = c2)
             #the which part is necessary since ceindex (pam medoids indices) go
             #from 1 to the TOTAL number of observations, while newx contains
             #only the observations that belong to this node
             #(sub is necessary to retrieve the index 'N' from 'clusterN')
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
             if(all(sapply(x, function(i) {
               is.null(edge.attributes(i)$weight)
               #if attribute weight is null for all the graphs, the graph
               #covariate is not weighted
             }))){
               adj_data <- lapply(x, igraph::as_adjacency_matrix)
             } else { #otherwise, it is weighted
               adj_data <- lapply(x, function(i) {
                 igraph::as_adjacency_matrix(i, attr = 'weight')
               })
             }
             #d is obtained in the same way in the two cases:
             d <- NetworkDistance::nd.centrality(adj_data, mode = 'Degree', directed = TRUE)
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


# Detect split.type -------------------------------------------------------

split.type_det <- function(object){

  # Check that object has class party
  stopifnot(inherits(object, 'party'))

  # Extract basid from the first node (which is necessarily present)
  basid_list <- nodeapply(object, by_node = TRUE, ids = 1,
                          FUN = function(node) basid_split(split_node(node)))

  # If basid is not null, it means we are in the coeff case
  if (!is.null(unlist(basid_list))){
    return(split.type = 'coeff')
  } else {
    return(split.type = 'other') #other means cluster or 'traditional' covariate
  }

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



