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
split.opt <- function(y, x, split.type = "coeff"){

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
                 x1 = x$coef
                 bselect <- 1:dim(x1)[2]
                 p1 <- c()
                 p1 <- sapply(bselect, function(i) mytestREG(x1[, i], y, R = R))
                 colnames(p1) <- colnames(x1)
                 if (length(which(p1[2,] == min(p1[2,], na.rm = T))) > 1) {
                   bselect <- which.max(p1[1,])
                 } else{
                   bselect <- which.min(p1[2,])
                 }

                 sel.coeff = x1[,bselect]
                 s  <- sort(sel.coeff)
                 comb = sapply(s, function(j) sel.coeff<j)
                 xp.value <- apply(comb, 2, function(q) mytestREG(x = q, y = y)[2])
                 splitindex <- s[which.min((xp.value))]
             }
             else if(split.type == "cluster"){
               cl.fdata = kmeans.fd(x, ncl=2, draw = FALSE, par.ini=list(method="exact"))
               splitindex <- cl.fdata
             }

           })
  return(splitindex)
}


# distances ---------------------------------------------------------------
# CHECK
compute.dissimilarity <- function(x,  lp = 2, case.weights, p = 2, dimension = 1){


    switch(class(x),
           logical    = dist((x[case.weights])),
           factor     = daisy(as.data.frame(x[case.weights,])),
           numeric    = dist(x[case.weights]),  # TBC: controlla se si possono accorpare condizioni sullo switch
           integer    = dist(x[case.weights]),
           data.frame = dist(x[case.weights]),
           matrix     = dist(x[case.weights]),
           fdata      = metric.lp(x[case.weights], lp=lp),
           list       = {
             if(attributes(x[[1]])$names == "diagram"){
               d1 = x[case.weights]
               k.fun = function(i, j) TDA::wasserstein(d1[[i]], d1[[j]], p=p, dimension = dimension)
               k.fun = Vectorize(k.fun)
               d.idx = seq_along(d1)
               outer(d.idx,d.idx, k.fun)

             }

           })

}



# Test --------------------------------------------------------------------
# CHECK
mytestREG <- function(x, y, R = 1000, dist.types = c("default", "default"), lp = c(2,2), case.weights) {

  d1 = compute.dissimilarity(x, dist.type = dist.types[1], lp = lp[1], case.weights = case.weights)
  d2 = compute.dissimilarity(y, dist.type = dist.types[2], lp = lp[2], case.weights = case.weights)

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
                      dist.types = rep("default",2),
                      lp = rep(2,2)) {

  if(!is.list(covariates)) stop("Argument 'covariates' must be provided as a list")

  n.var = length(covariates) #TBC: switch for type of covariates?

  # selects one covariate and performs an independece test
  p = sapply(covariates, function(sel.cov) mytestREG(x = sel.cov,
                                                     y = response,
                                                     R = R,
                                                     dist.types = dist.types,
                                                     lp = lp,
                                                     case.weights = case.weights))

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

   splitindex = split.opt(y = response, x = x, split.type = "coeff")

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
            return(list(
              sp = partysplit(
                varid = as.integer(bselect),
                index = splitindex,
                info = list(p.value = 1 - (1 - p[2,]) ^
                              sum(!is.na(p[2,])))),
                varselect = xselect
            ))
            }

   )

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
  n.var = length(covariates) #TBC: switch for type of covariates?

  # if the case weights are not provided, they are all initialized as 1
  if (is.null(case.weights))
    case.weights <- rep(1L, length(response))

  # new list of covariates (initialization)
  newcovariates=list()

  # trasformations based on the variables' nature
  ### covariates[[j]]
  newcovariates = lapply(covariates, function(j){
    switch(class(j),
           logical    = return(covariates[[j]]),
           factor     = return(covariates[[j]]),
           numeric    = return(covariates[[j]]),
           integer    = return(covariates[[j]]),
           data.frame = return(covariates[[j]]),
           matrix     = return(covariates[[j]]),
           fdata      = {
             foo <- fda.usc::min.basis(covariates[[j]], numbasis = nb)
             fd3 <- fda.usc::fdata2fd(foo$fdata.est,
                                      type.basis = "bspline",
                                      nbasis = foo$numbasis.opt)
             foo$coef <- t(fd3$coefs)
             return(foo)
           },
           list       = {
             if(all(sapply(covariates[[j]], class) == 'igraph')){
               shell <- graph.to.shellness.distr.df(covariates[[j]])
               return(shell)
             }
           }
    )
  }
  )

nodes <- growtree(id = 1L,
                  response = response,
                  covariates = newcovariates,
                  case.weights = case.weights,
                  minbucket = minbucket,
                  alpha = alpha,
                  R = R,
                  rnd.sel = rnd.sel,
                  rnd.splt = rnd.splt,
                  n.var = n.var)

  # compute terminal node number for each observation
  m.data <- c()

  newcovariates = lapply(covariates, function(j){
    switch(class(j),
           logical    = return(newcovariates[[j]]),
           factor     = return(newcovariates[[j]]),
           numeric    = return(newcovariates[[j]]),
           integer    = return(newcovariates[[j]]),
           data.frame = return(newcovariates[[j]]),
           matrix     = return(newcovariates[[j]]),
           fdata      = {
             foo <- newcovariates[[j]]$coef
             colnames(foo) <- paste(names(covariates)[j], colnames(newcovariates[[j]]$coef),sep = ".")
             return(foo)
           },
           list       = {
             if(all(sapply(covariates[[j]], class) == 'igraph')){
               foo <- newcovariates[[j]]
               colnames(foo) <- paste(names(covariates)[j], colnames(newcovariates[[j]]), sep = ".")
               return(foo)
             }
           }
    )
  }
  )

    # fill in m.data or initialize it
    if (!is.null(m.data)) {
      m.data <- cbind(m.data, foo)
    }
    else {
      m.data <- foo
    }
  }
  fitted <- fitted_node(nodes, data = data.frame(m.data))

  # return rich constparty object
  data1 = cbind(response, m.data)
  ret <- party(nodes, data = data.frame(m.data),
    fitted = data.frame("(fitted)" = fitted,
                        "(response)" = response,
                        "(case.weights)" = case.weights,
                        check.names = FALSE),
    terms = terms(response ~ ., data = data1))
  as.constparty(ret)
  }

  return(ret)
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
                     n.var) {
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

  # no split found, stop here
  if (is.null(res))
    return(partynode(id = id))

  # separately saving res_splt outputs
  sp <- res_splt$sp
  varselect <- res_splt$varselect


  #kidids_split suppressed since data may be arbirarily complex
  ### isn't it possible to use kidids_split anyway?
  kidids <- c()

  switch(class(covariates[[varselect]]),
         fdata = {
           if(split.type == "coeff"){

             # observations before the split point are assigned to node 1
               kidids[which(covariates[[varselect]]$coef[, sp$varid] <= sp$breaks)] <- 1
            #  observations before the split point are assigned to node 2
              kidids[which(covariates[[varselect]]$coef[, sp$varid] > sp$breaks)] <- 2

           # number of observations assigned to node 1
           sum1 <-
             length(which(covariates[[varselect]]$coef[which(case.weights == 1), sp$varid] <= sp$breaks))
           # number of observations assigned to node 2
           sum2 <-
             length(which(covariates[[varselect]]$coef[which(case.weights == 1), sp$varid] > sp$breaks))

           # vector containing the optimal number of basis for each covariate
           nb = sapply(1:n.var, function(i) covariates[[i]]$numbasis.opt)
           ###partiva da 0! dove viene utilizzato?

           # shift the varid of the tree based on the quantity of the
           # previous features/basis
           # Ex: if variable 3 is selected for splitting (variable 1
           # is the response, it's ignored), then shift varid by the
           # number of basis of variable 2 (if it's functional) or the
           # maximum k_core found in the graphs (if it's a graph)
           total_features <- lapply(covariates[2:n.var],
                                    function(v){switch(class(v),
                                                       logical    = 1,
                                                       factor     = 1,
                                                       numeric    = 1,
                                                       integer    = 1,
                                                       data.frame = ncol(covariates[[v]]),
                                                       matrix     = ncol(covariates[[v]]),
                                                       fdata      = covariates[[v]]$numbasis.opt
                                    )
                                    }
           )

           step <- sum(total_features[2:n.var[which(2:n.var < varselect)]], na.rm = T)
           sp$varid = sp$varid + as.integer(step)

           } else if (split.type == "cluster"){

             kidids[sp$index == 1] <- 1
             kidids[sp$index == 2] <- 2


         }
           },
         numeric = {
           kidids[(which(covariates[[varselect]][, sp$varid] <= sp$breaks))] <- 1
           kidids[(which(covariates[[varselect]][, sp$varid] > sp$breaks))] <- 2

           sum1 <-
             length(which(covariates[[varselect]][, sp$varid][which(case.weights == 1)] <= sp$breaks))
           sum2 <-
             length(which(covariates[[varselect]][, sp$varid][which(case.weights == 1)] > sp$breaks))
         },

         integer = {
           kidids[(which(covariates[[varselect]][, sp$varid] <= sp$breaks))] <- 1
           kidids[(which(covariates[[varselect]][, sp$varid] > sp$breaks))] <- 2

           sum1 <-
             length(which(covariates[[varselect]][, sp$varid][which(case.weights == 1)] <= sp$breaks))
           sum2 <-
             length(which(covariates[[varselect]][, sp$varid][which(case.weights == 1)] > sp$breaks))
         },

         factor = {

           kidids[sp$index == 1] <- 1
           kidids[sp$index == 2] <- 2

         }
         )

  # if all the observations belong to the same node, no split is done
  if (all(kidids == 1) | all(kidids == 2))
    return(partynode(id = id))

  if ((sum1 == 0 |
       sum2 == 0)) {
    ###differenza col precedente if?
    return(partynode(id = id))
  }

  # setup all daugther nodes
  kids <-
    vector(mode = "list", length = max(kidids, na.rm = TRUE))

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
