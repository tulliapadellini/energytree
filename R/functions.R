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
split.opt <- function(y, x, split.type = "coeff", rnd = T){

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

           })
  return(splitindex)
}




# distances ---------------------------------------------------------------
# CHECK
compute.dissimilarity <- function(x, dist.type = "default", lp = 2, case.weights){

  if(dist.type == "default"){
    switch(class(x),
           logical    = dist((x[case.weights])),
           factor     = daisy(as.data.frame(x[case.weights,])),
           numeric    = dist(x[case.weights]),  # TBC: controlla se si possono accorpare condizioni sullo switch
           integer    = dist(x[case.weights]),
           data.frame = dist(x[case.weights]),
           matrix     = dist(x[case.weights]),
           fdata      = metric.lp(x[case.weights], lp=lp))
  }
  else{
    if(!is.list(x)) stop("if distance function is arbitrary, argument x must be provided as a list")
    dist.type(x[case.weights]) # TBC: mettilo in modo che restituisca matrice di distanze
  }
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

   splitindex = split.opt(y = y, x = x, split.type = "coeff")



  # return split as partysplit object
  if (is.numeric(x)) {
    return(partysplit(
      varid = as.integer(xselect),
      breaks = splitindex,
      info = list(p.value = 1 - (1 - p) ^ sum(!is.na(p)))
    ))
  }
  if (is.data.frame(x)) {
    temp = list (
      sp = partysplit(
        varid = as.integer(cselect),
        breaks = splitindex,
        info = list(p.value = 1 - (1 - p) ^ sum(!is.na(p)))
      ),
      varselect = xselect
    )
    return(temp)
  }
  if (is.factor(x)) {
    return(partysplit(
      varid = as.integer(xselect),
      index = splitindex,
      info = list(p.value = 1 - (1 - p) ^ sum(!is.na(p)))
    ))
  }
  if (is.list(x)) {
    if (is.fdata(x$fdata.est)) {
      return(list(
        sp = partysplit(
          varid = as.integer(bselect),
          breaks = splitindex,
          info = list(p.value = 1 - (1 - p[2,]) ^
                        sum(!is.na(p[2,])))
        ),
        varselect = xselect
      ))
    }
  }
}






#implement a conditional tree -- switches between different kind of covariates

mytree <- function(Y,  # nome della variabile risposta
                   data, # lista contentente risposta e covariate
                   weights = NULL,
                   minbucket = 1,
                   alpha = 0.05,
                   R = 1000,
                   rnd.sel = T,
                   rnd.splt = TRUE,
                   nb = 5) {

  # change with checks if response is in the right format
  response <- data[[which(names(data) == Y)]]

  if (is.null(weights))
    weights <- rep(1L, length(response))

  # length data
  n.var <- which(names(data) != Y)

  # change with switch?
  if (class(data) == "list") {
    datanew <- list("response" = response)
    for (j in n.var) {

      if (class(data[[j]]) == "fdata") {
        foo <- min.basis(data[[j]], numbasis = nb)
        fd3 <-
          fdata2fd(foo$fdata.est,
                   type.basis = "bspline",
                   nbasis = foo$numbasis.opt)
        foo$coef <- t(fd3$coefs)
        datanew[[j]] <- foo
      }
      else if (class(data[[j]]) == "list" &
               class(data[[j]][[1]]) == "igraph") {
        datanew[[j]] <- graph.to.shellness.distr.df(data[[j]])
      }

      # if persistence diagram
      else if (class(data[[j]]) == "data.frame") {
        datanew[[j]] = data[[j]]
      }
    }
    names(datanew)[-1] <- names(data)[-1]
  }
  else if (class(data) == "data.frame") {
    datanew = data
    colnames(datanew)[colnames(datanew) == "Y"] <- "response"
  }

  nodes <-
    growtree(
      id = 1L,
      response = datanew$response,
      data = datanew,
      weights,
      minbucket = minbucket,
      alpha = alpha,
      R = R,
      rnd.sel = rnd.sel,
      rnd.splt = rnd.splt,
      n.var = n.var
    )

  # compute terminal node number for each observation
  response <- response
  response <- data.frame(response)
  y = response
  m.data <- c()

  for (j in n.var) {
    if (class(data[[j]]) == "fdata") {
      foo <- datanew[[j]]$coef
      colnames(foo) <-
        paste(names(data)[j], colnames(datanew[[j]]$coef),
              sep = ".")
    }
    else if (class(data[[j]]) == "data.frame" |
             (class(data[[j]]) == "list" &
              class(data[[j]][[1]]) == "igraph")) {
      foo <- datanew[[j]]
      colnames(foo) <-
        paste(names(data)[j], colnames(datanew[[j]]), sep = ".")
    }

    # fill in m.data or initialize it
    if (!is.null(m.data)) {
      m.data <- cbind(m.data, foo)
    }
    else {
      m.data <- foo
    }
  }

  data1 = cbind(response, m.data)
  m.data = m.data
  data1 = data1

  fitted <- fitted_node(nodes, data = data.frame(m.data))
  formula = response ~ .

  # return rich constparty object
  ret <- party(
    nodes,
    data = data.frame(m.data),
    fitted = data.frame(
      "(fitted)" = fitted,
      "(response)" = data1$response,
      "(weights)" = weights,
      check.names = FALSE
    ),
    terms = terms(formula, data = data1)
  )

  as.constparty(ret)
}
