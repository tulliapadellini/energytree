library(cluster)
library(fda.usc)




# split value -------------------------------------------------------------

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
split.opt <- function(y, x, split.type = "coeff", R=1000, wass.dist = NULL){

      switch(class(x),
           factor     = { #com'era prima ma overall na mezza merda

             lev <- levels(x[drop = TRUE])
             if (length(lev) == 2) {
               splitpoint <- lev[1]
             } else{
               comb <- do.call("c", lapply(1:(length(lev) - 1), ### TBC: isn't this just floor(length(lev)/2) ??
                                           function(x)
                                             combn(lev,
                                                   x,
                                                   simplify = FALSE)))
               xlogp <- sapply(comb, function(q)
                 mychisqtest(x %in% q, y))
               splitpoint <- comb[[which.min(xlogp)]]
             }

             # split into two groups (setting groups that do not occur
             # to NA)
             splitindex <- !(levels(x) %in% splitpoint)
             splitindex[!(levels(x) %in% lev)] <- NA_integer_
             splitindex <- splitindex - min(splitindex, na.rm = TRUE) + 1L
             },

           numeric    = {
             s  <- sort(x)
             comb = sapply(s, function(j) x<j)
             xp.value <- apply(comb, 2, function(q) mytestREG(x = q, y = y)[2])
             splitindex <- s[which.min((xp.value))]
             },

           integer    = {
             s  <- sort(x)
             comb = sapply(s, function(j) x<j)
             xp.value <- apply(comb, 2, function(q) mytestREG(x = q, y = y)[2])
             splitindex <- s[which.min((xp.value))]
           },
           fdata      = {
             if(split.type == "coeff"){
                 foo <- fda.usc::min.basis(x, numbasis = nb)
                 fd3 <- fda.usc::fdata2fd(foo$fdata.est,
                                        type.basis = "bspline",
                                        nbasis = foo$numbasis.opt)
                 foo$coef <- t(fd3$coefs)
                 x1 = foo$coef
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
                 comb = sapply(s, function(j) sel.coeff<j)
                 xp.value <- apply(comb, 2, function(q) mytestREG(x = q, y = y)[2])
                 splitindex <- s[which.min((xp.value))]
             }
             else if(split.type == "cluster"){
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
  if(class(x) == 'fdata') {out$bselect <- bselect}
  return(out)
}


# distances ---------------------------------------------------------------
# CHECK
compute.dissimilarity <- function(x, lp = 2, case.weights, ...){



    switch(class(x),
           logical    = dist((x[case.weights])),
           factor     = daisy(as.data.frame(x[case.weights,])),
           numeric    = dist(x),  # TBC: controlla se si possono accorpare condizioni sullo switch
           integer    = dist(x[case.weights]),
           data.frame = dist(x[case.weights]),
           matrix     = dist(x[case.weights]),
           fdata      = metric.lp(x, lp=lp))
           # list       = {
           #   if(!is.null(attributes(x[[1]]))){
           #   if(attributes(x[[1]])$names == "diagram"){
           #     d1 = x[case.weights]
           #     k.fun = function(i, j) TDA::wasserstein(d1[[i]], d1[[j]], ...)
           #     k.fun = Vectorize(k.fun)
           #     d.idx = seq_along(d1)
           #     outer(d.idx,d.idx, k.fun)
           #   }}
           #})

}




# Test --------------------------------------------------------------------
# CHECK
mytestREG <- function(x, y, R = 1000, lp = c(2,2), case.weights, ...) {

  d1 = compute.dissimilarity(x, lp = lp[1], case.weights = case.weights)
  if(is.null(d1)) print("OPSIE")
  d2 = compute.dissimilarity(y, lp = lp[2], case.weights = case.weights)
  if(is.null(d2)) print("OPSIE")


  ct <- energy::dcor.test(d1, d2, R = R)
  if (!is.na(ct$statistic)) {
    return(c(ct$statistic, ct$p.value))
  } else{
    c(NA, NA)
  }
}


# Find split --------------------------------------------------------------
#CHECK
findsplit <- function(response,
                      covariates,
                      case.weights,
                      alpha,
                      R,
                      rnd.sel,
                      rnd.splt,
                      lp = rep(2,2),
                      split.type = "coeff", ...) {

  if(!is.list(covariates)) stop("Argument 'covariates' must be provided as a list")

  n.var = length(covariates) #TBC: switch for type of covariates?

  # selects one covariate and performs an independece test
  p = lapply(covariates, function(sel.cov) mytestREG(x = sel.cov,
                                                     y = response,
                                                     R = R,
                                                     dist.types = dist.types,
                                                     lp = lp,
                                                     case.weights = case.weights))

   p = t(matrix(unlist(p), ncol = 2, byrow = T))

   rownames(p) <- c("statistic", "p-value")


   # Bonferroni correction
   if (all(is.na(p[2,]))) return(NULL)

   minp <- min(p[2,], na.rm = TRUE)
   minp <- 1 - (1 - minp) ^ sum(!is.na(p[2,]))
   if (minp > alpha) return(NULL)

   if (length(which(p[2,] == min(p[2,], na.rm = T))) > 1) {
     xselect <- which.max(p[1,])    # in case of multiple minima, take that with the highest test statistic
   } else{
     xselect <- which.min(p[2,])
   }

   x <-  covariates[[xselect]]

   split.objs = split.opt(y = response, x = x, split.type = "coeff")
   splitindex <- split.objs$splitindex
   bselect <- split.objs$bselect

   switch(class(x),
          numeric = {
            return(list(
              sp = partysplit(
              varid = as.integer(xselect),
              breaks = splitindex,
              info = list(p.value = 1 - (1 - p) ^ sum(!is.na(p)))
            ),
              varselect = xselect))
          },
          integer = {
            return(list(
              sp = partysplit(
                varid = as.integer(xselect),
                breaks = splitindex,
                info = list(p.value = 1 - (1 - p) ^ sum(!is.na(p)))),
              varselect = xselect))
          },
          factor = {
            return(list(
              sp = partysplit(
               varid = as.integer(xselect),
               index = splitindex,
               info = list(p.value = 1 - (1 - p) ^ sum(!is.na(p)))),
              varselect = xselect))
          },
          fdata = {
            if(split.type == "coeff"){
            return(list(
              sp = partysplit(
                varid = as.integer(bselect),
                breaks = splitindex,
                info = list(p.value = 1 - (1 - p[2,]) ^
                              sum(!is.na(p[2,])))),
                varselect = xselect
            ))} else if(split.type == "cluster"){
              return(list(
                sp = partysplit(
                  varid = as.integer(xselect),
                  index = splitindex,
                  info = list(p.value = 1 - (1 - p[2,]) ^
                                sum(!is.na(p[2,])))),
                varselect = xselect
              ))
            }
            },
          list = if(attributes(x[[1]])$names == 'diagram'){
            return(list(
              sp = partysplit(
                varid = as.integer(xselect),
                index = splitindex,
                info = list(p.value = 1 - (1 - p[2,]) ^
                              sum(!is.na(p[2,])))),
              varselect = xselect
              ))

   }

   )

}


# growtree ----------------------------------------------------------------

growtree <- function(id = 1L,
                     response,
                     covariates,
                     case.weights,
                     minbucket,
                     alpha,
                     R,
                     rnd.sel,
                     rnd.splt,
                     n.var,
                     split.type = 'coeff') {
  # for less than <minbucket> observations stop here
  if (sum(case.weights) < minbucket)
    return(partynode(id = id))

  # find best split
  res_splt <- findsplit(
    response,
    covariates,
    case.weights,
    alpha,
    R,
    rnd.sel,
    rnd.splt,
    dist.types = rep("default", 2),
    lp = rep(2, 2)
  )

  # separately saving res_splt outputs
  sp <- res_splt$sp
  varselect <- res_splt$varselect
  print(c('sp is', sp))
  print(c('varselect is', varselect))

  # no split found, stop here
  if (is.null(sp))
    return(partynode(id = id))

  #kidids_split suppressed since data may be arbirarily complex
  ### isn't it possible to use kidids_split anyway?
  kidids <- c()

  switch(class(covariates[[varselect]]),
         fdata = {
           if(split.type == "coeff"){

             basis_covariates = lapply(covariates, function(j){
               foo <- fda.usc::min.basis(j, numbasis = nb)
               fd3 <- fda.usc::fdata2fd(foo$fdata.est,
                                        type.basis = "bspline",
                                        nbasis = foo$numbasis.opt)
               foo$coef <- t(fd3$coefs)
               return(foo)
             }
             )

             # observations before the split point are assigned to node 1
             kidids[which(basis_covariates[[varselect]]$coef[, sp$varid] <= sp$breaks)] <-
               1
             #  observations before the split point are assigned to node 2
             kidids[which(basis_covariates[[varselect]]$coef[, sp$varid] > sp$breaks)] <-
               2

             # vector containing the optimal number of basis for each covariate
             nb = sapply(1:n.var, function(i) basis_covariates[[i]]$numbasis.opt)
             ###partiva da 0! dove viene utilizzato?

             # shift the varid of the tree based on the quantity of the
             # previous features/basis
             # Ex: if variable 3 is selected for splitting (variable 1
             # is the response, it's ignored), then shift varid by the
             # number of basis of variable 2 (if it's functional) or the
             # maximum k_core found in the graphs (if it's a graph)
             total_features <- lapply(basis_covariates,
                                      function(v) {
                                        switch(
                                          class(v),
                                          logical    = 1,
                                          factor     = 1,
                                          numeric    = 1,
                                          integer    = 1,
                                          data.frame = ncol(v),
                                          matrix     = ncol(v),
                                          fdata      = v$numbasis.opt
                                        )
                                      })

             step <- if(length(total_features) > 1) {
               sum_feat <- sum(total_features[[which(1:n.var <= varselect)]], na.rm = T)
               as.integer(sum_feat)
             } else {
               sum_feat <- 0L
             }

             sp$varid = sp$varid + as.integer(step)

           } else if (split.type == "cluster") {

             kidids[sp$index == 1] <- 1
             kidids[sp$index == 2] <- 2

         }
           },

         numeric = {

           kidids[(which(covariates[[varselect]][, sp$varid] <= sp$breaks))] <- 1
           kidids[(which(covariates[[varselect]][, sp$varid] > sp$breaks))] <- 2

         },

         integer = {

           kidids[(which(covariates[[varselect]][, sp$varid] <= sp$breaks))] <- 1
           kidids[(which(covariates[[varselect]][, sp$varid] > sp$breaks))] <- 2

         },

         factor = {

           kidids[sp$index == 1] <- 1
           kidids[sp$index == 2] <- 2

         }
         )

  # if all the observations belong to the same node, no split is done
  if (all(kidids == 1) | all(kidids == 2))
    return(partynode(id = id))

  # setup all daugther nodes
  kids <-
    vector(mode = "list", length = max(kidids, na.rm = TRUE))

  print(c('length is',length(kids)))

  for (kidid in 1:length(kids)) {
    # select observations for current node
    w <- case.weights
    w[kidids != kidid] <- 0

    # get next node id
    if (kidid > 1) {
      myid <- max(nodeids(kids[[kidid - 1]]))
    } else{
      myid <- id
    }


    # start recursion on this daugther node
    kids[[kidid]] <-
      growtree(
        id = as.integer(myid + 1),
        response,
        covariates,
        w,
        minbucket,
        alpha,
        R,
        rnd.sel,
        rnd.splt ,
        n.var = n.var
      )
  }


  # return nodes
  return(partynode(
    id = as.integer(id),
    split = sp,
    kids = kids,
    info = list(p.value = min(info_split(sp)$p.value, na.rm = TRUE))
  ))

}

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


# mytree ------------------------------------------------------------------

# implement a conditional tree -- switches between different kind of covariates

mytree <- function(response,
                   covariates,
                   case.weights = NULL,
                   minbucket = 1,
                   alpha = 0.05,
                   R = 1000,
                   rnd.sel = T,
                   rnd.splt = TRUE,
                   nb = 5) {

  if(!is.list(covariates)) stop("Argument 'covariates' must be provided as a list")

  # number of covariates
  n.var = length(covariates)

  # if the case weights are not provided, they are all initialized as 1
  if (is.null(case.weights))
    case.weights <- rep(1L, as.numeric(length(response)))

  # new list of covariates (initialization)
  newcovariates=list()

  # trasformations based on the variables' nature
  newcovariates = lapply(covariates, function(j){
    if(class(j) == 'fdata'){

      foo <- fda.usc::min.basis(j, numbasis = nb)
      fd3 <- fda.usc::fdata2fd(foo$fdata.est,
                               type.basis = "bspline",
                               nbasis = foo$numbasis.opt)
      foo$coef <- t(fd3$coefs)
      return(foo)

    } else
      if(class(j) == 'list' &
         all(sapply(j, class) == 'igraph')){

        shell <- graph.to.shellness.distr.df(j)
        return(shell)

      } else {

        return(j)

      }
  }
  )

  nodes <- growtree(id = 1L,
                    response = response,
                    covariates = covariates,
                    case.weights = case.weights,
                    minbucket = minbucket,
                    alpha = alpha,
                    R = R,
                    rnd.sel = rnd.sel,
                    rnd.splt = rnd.splt,
                    n.var = n.var)


  # compute terminal node number for each observation
  # m.data <- lapply(1:n.var, function(j){
  #   if(class(covariates[[j]]) == 'fdata'){
  #
  #     foo <- newcovariates[[j]]$coef
  #     colnames(foo) <- paste(names(covariates)[j], colnames(newcovariates[[j]]$coef),sep = ".")
  #     return(foo)
  #
  #   } else {
  #
  #     foo <- newcovariates[[j]]
  #     colnames(foo) <- paste(names(covariates)[j], colnames(newcovariates[[j]]), sep = ".")
  #     return(foo)
  #
  #   }
  # }
  # )
  ### m.data probably not needed, thinking about working directly on newcovariates

  fitted.obs <- fitted_node(nodes,
                            data = newcovariates[[1]]$fdataobj$data,
                            vmatch = 1:length(newcovariates[[1]]$fdataobj$data),
                            obs = 1:nrow(newcovariates[[1]]$fdataobj$data))

  # # old fitted.obs (original data)
  # fitted.obs <- fitted_node(nodes,
  #                           data = covariates[[1]]$data,
  #                           vmatch = 1:length(covariates[[1]]$data),
  #                           obs = 1:nrow(covariates[[1]]$data))

  prova.dati = data.frame(y = response)
  for(i in 1:n.var){
    prova.dati[i+1] = I(covariates[i])
  }


  # return rich constparty object
  data1 = cbind('response'=as.data.frame(response), as.data.frame(newcovariates[[1]]$fdataobj$data))
  ret <- party(nodes, data = as.data.frame(newcovariates[[1]]$fdataobj$data),
               fitted = data.frame("(fitted)" = fitted.obs,
                                   "(response)" = as.data.frame(response),
                                   "(case.weights)" = case.weights,
                                   check.names = FALSE),
               terms = terms(response ~ ., data = data1))


  # ?? as.constparty(ret)

  return(ret)
}

