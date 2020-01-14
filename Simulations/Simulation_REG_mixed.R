
# Inizialization ----------------------------------------------------------

# Loading the libraries
library(fda.usc)
library(roahd)
library(energy)
library(entropy)
library(partykit)
library(igraph)
library(NetworkDistance)
source("functions.R")

# Loading the dataset
load("Regression_Sim_Dataset.RData")

# Errors
MEP_etree <- c()


# Response and covariates lists construction ------------------------------

# Response
resp <- lapply(data, function(x) x$Y)[[1]]
#remark: we take the first element since the dataset contains 100 simulations

### Regression with a numeric predictor and a functional one (two alternatives) ###

#1) numerical variable: response+noise
cov.list <- list(lapply(data, function(x) x$V1)[[1]], resp+rnorm(200, 0, 10))
#ok

#2) numerical variable: one of the basis of another simulation of the same dataset
foo <- fda.usc::optim.basis(lapply(data, function(x) x$V1)[[2]], numbasis = 15)
fd3 <- fda.usc::fdata2fd(foo$fdata.est,
                         type.basis = "bspline",
                         nbasis = foo$numbasis.opt)
foo$coef <- t(fd3$coefs)
cov.list <- list(lapply(data, function(x) x$V1)[[1]], foo$coef[,15])
#ok

### Regression with three covariates: functional, graph & numeric ###

# Normalization of resp to make it the connection probability of each graph
norm01 <- function(x){(x-min(x))/(max(x)-min(x))}
conn.prob <- norm01(resp)

# Avoiding too low connection probabilities
conn.prob[which(conn.prob < 0.05)] <- conn.prob[which(conn.prob < 0.05)] + 0.1

# Generation of the graphs
graph.list <- lapply(conn.prob, function(p){sample_gnp(100, p)})

# Covariates for the full mixed model
cov.list <- list(lapply(data, function(x) x$V1)[[1]], graph.list, foo$coef[,15])


# Model fitting -----------------------------------------------------------

# Number of basis
n.bas <- 15

### REGRESSION ENERGY TREE ###
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


# Prediction --------------------------------------------------------------

### ETREE PREDICTION ###

# Prediction
y_pred <- predict(etree_fit)

# Prediction with newdata
graph.list2 <- lapply(conn.prob, function(p){sample_gnp(100, p)})
new.cov.list <- list(lapply(data, function(x) x$V1)[[2]], graph.list2, foo$coef[,15])
y_pred2 <- predict(etree_fit, newdata = new.cov.list)


# Error
y <- resp
MEP_etree <- (sum((y-y_pred)^2)/length(y))/(var(y))
