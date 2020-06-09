
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


# Response and covariates lists construction ------------------------------

# Restriction on the observations' number (for computational reasons)
obs <- sample(1:798, 150)

# Response
resp <- lapply(data, function(x) x$cls[obs])[[1]]
#remark: we take the first element since the dataset contains 100 simulations

### Classification with a functional predictors and a numeric one ###
#the numeric is one of the basis of another simulation of the same dataset
foo <- fda.usc::optim.basis(lapply(data, function(x) fdata(x[obs,2:129]))[[2]], numbasis = 15)
fd3 <- fda.usc::fdata2fd(foo$fdata.est,
                         type.basis = "bspline",
                         nbasis = foo$numbasis.opt)
foo$coef <- t(fd3$coefs)
cov.list <- list(lapply(data, function(x) fdata(x[obs,2:129]))[[1]],
                 foo$coef[,8])

### Classification with three covariates: functional, graph & numeric ###

# Generation of the graphs with a different connection probability for each class
graph.list <- lapply(resp,
                     function(c){
                       if (c == 'Bel'){
                         sample_gnp(100, 0.10)
                       } else if (c == 'Cyl'){
                         sample_gnp(100, 0.115)
                       } else if (c == 'Fun'){
                         sample_gnp(100, 0.13)
                       }
                     })

# Covariates for the full mixed model
cov.list <- list(lapply(data, function(x) fdata(x[obs,2:129]))[[1]], graph.list, foo$coef[,8])


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
                   coef.split.type = 'test')
plot(etree_fit)


# Prediction --------------------------------------------------------------

### ETREE CLASSIFICATION PREDICTION ###

# Prediction
y_pred <- predict(etree_fit)

# Prediction with newdata
graph.list2 <- lapply(resp,
                     function(c){
                       if (c == 'Bel'){
                         sample_gnp(100, 0.10)
                       } else if (c == 'Cyl'){
                         sample_gnp(100, 0.125)
                       } else if (c == 'Fun'){
                         sample_gnp(100, 0.15)
                       }
                     })
new.cov.list <- list(lapply(data, function(x) fdata(x[obs,2:129]))[[2]], graph.list2, foo$coef[,8])
y_pred2 <- predict(etree_fit, newdata = new.cov.list)

# Error
y <- resp
t <- table(y_pred, y)
ACC_etree <- sum(diag(t))/(length(y))


