
# Inizialization ----------------------------------------------------------

# Loading the libraries
library(fda.usc)
library(energy)
library(entropy)
library(partykit)
source("functions.R")

# Loading the dataset
load("sim_data.RData")

# Errors
ACC_etree <- c()
ACC_ftree <- c()


# Response and covariates lists construction ------------------------------

# Restriction on the observations' number (for computational reasons)
obs <- sample(1:798, 150)

# Response
resp <- lapply(data, function(x) x$cls[obs])[[1]]
#remark: we take the first element since the dataset contains 100 simulations

# Classification with one functional predictor
cov.list <- lapply(data, function(x) fdata(x[obs,2:129]))[1]

# Classification with many functional predictors
cov.list <- list(lapply(data, function(x) fdata(x[obs,2:129]))[[1]],
                 lapply(data, function(x) fdata(x[obs,2:129]))[[2]])

# Classification with two functional predictors and a numeric one
#the numeric is one of the basis of another simulation of the same dataset
foo <- fda.usc::min.basis(lapply(data, function(x) fdata(x[obs,2:129]))[[2]], numbasis = 15)
fd3 <- fda.usc::fdata2fd(foo$fdata.est,
                         type.basis = "bspline",
                         nbasis = foo$numbasis.opt)
foo$coef <- t(fd3$coefs)
cov.list <- list(lapply(data, function(x) fdata(x[obs,2:129]))[[1]],
                 lapply(data, function(x) fdata(x[obs,2:129]))[[2]],
                 foo$coef[,8])


# Model fitting -----------------------------------------------------------

# Number of basis
n.bas <- 15

### CLASSIFICATION ENERGY TREE ###
etree_fit <- etree(response = resp,
                   covariates = cov.list,
                   case.weights = NULL,
                   minbucket = 5,
                   alpha = 0.05,
                   R = 1000,
                   split.type = 'coeff',
                   coef.split.type = 'test',
                   nb = n.bas)
plot(etree_fit)

### FUNCTIONAL CLASSIFICATION TREE ###

# One functional covariate
list_all <- list("df" = as.data.frame(resp), "x" = cov.list[[1]])
fvb <- create.fdata.basis(cov.list[[1]],l=1:15)
ct_fit <- classif.tree(resp ~ x, data = list_all, basis.b = list(x =fvb))

# # più covariate funzionali
# list_all <- list("df" = as.data.frame(cbind(resp, cov.list[[3]])), "x" = list(cov.list[[1]], cov.list[[2]]))
# fvb1 <- create.fdata.basis(cov.list[[1]],l=1:15)
# a2 <- classif.tree(resp ~ x, data = list_all, basis.b = list(x1 = fvb1, x2 = fvb2))
# # errore: può essere che si possa utilizzare al massimo una covariata funzionale?
#
# # una coavariata funzionale e una numerica
# list_all <- list("df" = as.data.frame(cbind(resp, cov.list[[3]])), "x" = cov.list[[1]])
# fvb <- create.fdata.basis(cov.list[[1]],l=1:15)
# a2 <- classif.tree(resp ~ list_all$df$V2 + x, data = list_all, basis.b = list(x =fvb))
# predict(a2)
# t <- table(predict(a2))
# acc_m1[i] <- sum(diag(t))/(length(y))
# # errore: forse non funziona neanche questo caso?

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


### FUNCTIONAL CLASSIFICATION TREE ###
t <- table(predict(ct_fit))
ACC_ftree <- sum(diag(t))/(length(y))


# Storing the results
save(ACC_etree, ACC_ftree, file = "results.RData")



