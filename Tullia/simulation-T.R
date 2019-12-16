
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

      foo <- fda.usc::min.basis(j, numbasis = nb)
      fd3 <- fda.usc::fdata2fd(foo$fdata.est,
                               type.basis = "bspline",
                               nbasis = foo$numbasis.opt)
      foo$coef <- t(fd3$coefs)
      return(foo)

    } else if(class(j) == 'list' &
              all(sapply(j, class) == 'igraph')){

      shell <- graph.to.shellness.distr.df(j)
      return(shell)

    } else {

      return(j)

    }
  }
  )

  # Building a df with all the new 'variables'
  newcovariates.onlybasis <- newcovariates
  for(i in 1:n.var){
    if(class(covariates[[i]]) == 'fdata') {
      newcovariates.onlybasis[[i]] <- newcovariates.onlybasis[[i]]$coef
    }
  }
  newcovariates.df <- as.data.frame(do.call(cbind, newcovariates.onlybasis))
  names(newcovariates.df) <- 1:ncol(newcovariates.df)

  # Growing the tree (finds the split rules)
  nodes <- growtree(id = 1L,
                    response = response,
                    covariates = covariates,
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
  data1 = cbind('response' = as.data.frame(response), newcovariates.df)
  names(data1) <- c('response', 1:(ncol(data1)-1))
  ret <- party(nodes,
               data = newcovariates.df,
               fitted = data.frame("(fitted)" = fitted.obs,
                                   "(response)" = data1$response,
                                   check.names = FALSE),
               terms = terms(response ~ ., data = data1))

  return(as.constparty(ret))

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

  # New list of covariates (here again, since it must be done at each split)
  newcovariates = lapply(covariates, function(j){
    if(class(j) == 'fdata'){

      foo <- fda.usc::min.basis(j, numbasis = nb)
      fd3 <- fda.usc::fdata2fd(foo$fdata.est,
                               type.basis = "bspline",
                               nbasis = foo$numbasis.opt)
      foo$coef <- t(fd3$coefs)
      return(foo)

    } else if(class(j) == 'list' &
              all(sapply(j, class) == 'igraph')){

      shell <- graph.to.shellness.distr.df(j)
      return(shell)

    } else {

      return(j)

    }
  }
  )

  # Finding the best split (variable selection & split point search)
  res_splt <- findsplit(response = response,
                        covariates = covariates,
                        newcovariates = newcovariates,
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
  switch(class(covariates[[varselect]]),

         fdata = {

           if(split.type == 'coeff'){

             # observations before the split point are assigned to node 1
             kidids[which(newcovariates[[varselect]]$coef[, sp$varid] <= sp$breaks)] <- 1
             #  observations before the split point are assigned to node 2
             kidids[which(newcovariates[[varselect]]$coef[, sp$varid] > sp$breaks)] <- 2

           } else if (split.type == 'cluster') {

             kidids[sp$index == 1] <- 1
             kidids[sp$index == 2] <- 2

           }
         },

         numeric = {

           kidids[(which(covariates[[varselect]] <= sp$breaks))] <- 1
           kidids[(which(covariates[[varselect]] > sp$breaks))] <- 2

         },

         integer = {

           kidids[(which(covariates[[varselect]] <= sp$breaks))] <- 1
           kidids[(which(covariates[[varselect]] > sp$breaks))] <- 2

         },

         factor = {

           kidids[sp$index == 1] <- 1
           kidids[sp$index == 2] <- 2

         }
  )

  # Total number of features for each covariate
  total_features <- lapply(covariates,
                           function(v) {
                             switch(
                               class(v),
                               logical    = 1,
                               factor     = 1,
                               numeric    = 1,
                               integer    = 1,
                               matrix     = ncol(v),
                               fdata      = {
                                 foo <- fda.usc::min.basis(v, numbasis = nb)
                                 foo$numbasis.opt                               }
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
  if(class(covariates[[varselect]]) == 'fdata'){
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
    kids[[kidid]] <-
      growtree(
        id = as.integer(myid + 1),
        response = subset(response, as.logical(w)),
        covariates = lapply(covariates, function(cov) subset(cov, as.logical(w))),
        case.weights = rep(1L, sum(w)),
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
                      newcovariates,
                      alpha,
                      R,
                      lp = rep(2,2),
                      split.type = 'coeff',
                      coef.split.type = 'test',
                      nb) {

  # Performing an independence test between the response and each covariate
  p = lapply(covariates, function(sel.cov) mytestREG(x = sel.cov,
                                                     y = response,
                                                     R = R,
                                                     lp = lp))
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

  # Selected covariate
  x <-  covariates[[xselect]]
  newx <- newcovariates[[xselect]]

  # Split point search
  split.objs = split.opt(y = response,
                         x = x,
                         newx = newx,
                         split.type = split.type,
                         coef.split.type = coef.split.type,
                         nb = nb)

  # Separately saving split.objs outputs
  splitindex <- split.objs$splitindex
  bselect <- split.objs$bselect

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
                                         index = splitindex,
                                         info = list(p.value = 1-(1-p[2,])^sum(!is.na(p[2,])))),
                         varselect = xselect))
           }

         },

         list = if(attributes(x[[1]])$names == 'diagram'){
           return(list(sp = partysplit(varid = as.integer(xselect),
                                       index = splitindex,
                                       info = list(p.value = 1-(1-p[2,])^sum(!is.na(p[2,])))),
                       varselect = xselect))

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
           xp.value <- apply(comb, 2, function(q) mytestREG(x = q, y = y))
           if (length(which(xp.value[2,] == min(xp.value[2,], na.rm = T))) > 1) {
             splitindex <- s[which.max(xp.value[1,])]
           } else {
             splitindex <- s[which.min(xp.value[2,])]
           }

         },

         integer    = {

           s  <- sort(x)
           comb = sapply(s[2:(length(s)-1)], function(j) x<j)
           xp.value <- apply(comb, 2, function(q) mytestREG(x = q, y = y))
           if (length(which(xp.value[2,] == min(xp.value[2,], na.rm = T))) > 1) {
             splitindex <- s[which.max(xp.value[1,])]
           } else {
             splitindex <- s[which.min(xp.value[2,])]
           }

         },

         fdata      = {

           if(split.type == 'coeff'){
             x1 = newx$coef
             bselect <- 1:dim(x1)[2]
             p1 <- c()
             p1 <- sapply(bselect, function(i) mytestREG(x1[, i], y, R = R))
             colnames(p1) <- colnames(x1)
             if (length(which(p1[2,] == min(p1[2,], na.rm = T))) > 1) {
               bselect <- as.integer(which.max(p1[1,]))
             } else{
               bselect <- as.integer(which.min(p1[2,]))
             }
             sel.coeff = x1[,bselect]
             s  <- sort(sel.coeff)
             comb = sapply(s[2:(length(s)-1)], function(j) sel.coeff<j)

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

               xp.value <- apply(comb, 2, function(q) mytestREG(x = q, y = y))
               if (length(which(xp.value[2,] == min(xp.value[2,], na.rm = T))) > 1) {
                 splitindex <- s[which.max(xp.value[1,])]
               } else {
                 splitindex <- s[which.min(xp.value[2,])]
               }

             }

           } else if(split.type == 'cluster') {
             cl.fdata = kmeans.fd(x, ncl=2, draw = FALSE, par.ini=list(method="exact"))
             splitindex <- cl.fdata$cluster
           }

         },

         list = if(attributes(x[[1]])$names == "diagram"){

           cl.diagrams = cluster::pam(wass.dist, k = 2, diss = TRUE)
           splitindex <- cl.diagrams$clustering

         }
  )

  out <- list('splitindex' = splitindex)
  if(class(x) == 'fdata') out$bselect <- bselect
  return(out)

}



# Independence (dcor) test ------------------------------------------------

mytestREG <- function(x,
                      y,
                      R = 1000,
                      lp = c(2,2)) {

  # Computing the dissimilarities within x and y
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
         logical    = dist(x),
         factor     = as.matrix(cluster::daisy(as.data.frame(x))),
         numeric    = dist(x),
         integer    = dist(x),
         matrix     = dist(x),
         fdata      = metric.lp(x, lp=lp))
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

graph.to.shellness.distr.df <- function(data, shell.limit = NULL) {
  tot.graphs = length(data)

  list.df <- list()
  max.shellness = 0

  for (i in 1:tot.graphs) {
    g = data[[i]]
    coreness.distr = count(coreness(g)) # aggr. by count
    rownames(coreness.distr) <-
      coreness.distr$x # re-index the df by the shellness number

    # keep just the frequency column
    coreness.distr = coreness.distr[c('freq')]

    # transpose the df. Convert the column-df into row-df.
    #This will ease the join with df.shellness.distr
    coreness.distr = t(coreness.distr)
    list.df[[i]] <- coreness.distr

    this.max.shellness = colnames(coreness.distr)[
      ncol(coreness.distr)]

    # update the maximum shellness found so far (used to build
    #the df of shellness distr)
    if (this.max.shellness > max.shellness) {
      max.shellness = this.max.shellness
    }
  }

  if (!is.null(shell.limit)) {
    # calculates the max shellness between the number of
    # predictors used in train set and the one calculated in
    # test set
    max.shellness <-
      if (as.numeric(max.shellness) < shell.limit - 1)
        shell.limit - 1
    else
      max.shellness
  }

  col.names = seq(0, max.shellness, 1)
  col.names = lapply(col.names, function(x)
    as.character(x))  # convert to char


# Inizialization ----------------------------------------------------------

  # Loading the libraries
  library(fda.usc)
  library(energy)
  library(entropy)
  library(partykit)
  library(mlr)
  library(fastDummies)
  source("functions.R")

  # Loading the dataset
  data('iris')

  # Errors
  ACC_etree <- c()
  ACC_mlc <- c()


  # Response and covariates lists construction ------------------------------

  # Response
  resp <- iris$Species

  # Covariates
  cov.list <- list('Sepal.Length' = iris$Sepal.Length, 'Sepal.Width' = iris$Sepal.Width, 'Petal.Length' = iris$Petal.Length, 'Petal.Width' = iris$Petal.Width)


  # Model fitting -----------------------------------------------------------

  ### CLASSIFICATION ENERGY TREE ###
  etree_fit <- etree(response = resp,
                     covariates = cov.list,
                     case.weights = NULL,
                     minbucket = 1,
                     alpha = 0.05,
                     R = 1000)
  plot(etree_fit)

  ### MULTILABEL CLASSIFICATION VIA CLASSIFIER CHAINS (MLR PACKAGE) ###
  # Transforming response into logical dummy variables (required by makeMultilabelTask)
  resp_dummy <- dummy_cols(resp)[,-1]
  resp_dummy <- sapply(resp_dummy, as.logical)
  colnames(resp_dummy) <- c('setosa', 'versicolor', 'virginica')
  # All data (coariates and dummies for the response) together
  all_data <- cbind(as.data.frame(do.call(cbind, cov.list)), resp_dummy)
  # Multilabel infrastucture (i.e. setting data and target)
  iris_task <- makeMultilabelTask(data = all_data, target = colnames(resp_dummy))
  # Base Learner
  binary.learner = makeLearner("classif.rpart")
  # Multilabel learner (wrapper of repetitions of the base one)
  mcc = makeMultilabelClassifierChainsWrapper(binary.learner)
  # Train and test sets
  n = getTaskSize(iris_task)
  train_set = sample(1:n, round(n*0.8))
  test_set = (1:n)[!(1:n) %in% train_set]
  # Train the multilabel learner
  iris_train = train(mcc, iris_task, subset = train_set)
  # Prediction
  iris_pred = predict(iris_train, task = iris_task, subset = test_set)
  # Accuracy
  performance(iris_pred, measures = list(multilabel.acc))

  # Prediction --------------------------------------------------------------

  ### CLASSIFICATION ENERGY TREE PREDICTION ###

  # New covariates

  new.cov.list = lapply(cov.list, function(j){
    if(class(j) == 'fdata'){

      foo <- fda.usc::min.basis(j, numbasis = n.bas)
      fd3 <- fda.usc::fdata2fd(foo$fdata.est,
                               type.basis = "bspline",
                               nbasis = foo$numbasis.opt)
      foo$coef <- t(fd3$coefs)
      return(foo$coef)

    } else if(class(j) == 'list' &
              all(sapply(j, class) == 'igraph')){

      shell <- graph.to.shellness.distr.df(j)
      return(shell)

    } else {

      return(j)

    }
  }
  )

  # New covariates dataframe
  new.cov.df <- as.data.frame(do.call(cbind, new.cov.list))
  names(new.cov.df) <- 1:ncol(new.cov.df)

  # Prediction
  y_pred <- predict(etree_fit, newdata = new.cov.df)

  # Error
  y <- resp
  t <- table(y_pred, y)
  ACC_etree <- sum(diag(t))/(length(y))


  ### MULTILABEL CLASSIFICATION VIA CLASSIFIER CHAINS PREDICTION ###
  ACC_mlc <- performance(iris_pred, measures = list(multilabel.acc))


  # Storing the results
  save(ACC_etree, ACC_mlc, file = "results.RData")




  df.shellness.distr = data.frame(matrix(
    data = NA_integer_,
    nrow = tot.graphs,
    ncol = length(col.names)
  )) #df with all graphs shellness distribution
  colnames(df.shellness.distr) <- col.names

  # fill in the df with the shellness distribution of each graph
  for (i in 1:tot.graphs) {
    updated.cols = colnames(list.df[[i]])

    for (x in updated.cols) {
      df.shellness.distr[i, x] = list.df[[i]][, x]
    }
  }

  df.shellness.distr[is.na(df.shellness.distr)] <-
    0 # replace NA by 0

  # converted the df columns to integer
  df.shellness.distr[, seq(1, ncol(df.shellness.distr))] <-
    sapply(df.shellness.distr[, seq(1, ncol(df.shellness.distr))],
           as.integer)

  return(df.shellness.distr[,1])
}
