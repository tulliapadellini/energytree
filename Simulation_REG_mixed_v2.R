
# Load packages, functions and data -------------------------------------------

# Packages
library(fda.usc)
library(energy)
library(entropy)
library(partykit)
library(cluster)
library(igraph)
library(NetworkDistance)
library(ggparty)
library(checkmate)

# Functions
source("functions_v2.R")
source("node_v2.R")
source("split_v2.R")
source("party_v2.R")
source("plot_v2.R")
source("ggparty_v1.R")
source("get_plot_data_v1.R")

# Data (following CATS by Serban and Wasserman, 2005)
load("Regression_Sim_Dataset.RData")


# Response and covariates -----------------------------------------------------

# Response
resp <- lapply(data, function(x) x$Y)[[1]]
#remark: we take the first element since the dataset contains 100 simulations

### Numerical covariate ###

# One of the basis of another simulation of the same dataset
foo <- fda.usc::optim.basis(lapply(data, function(x) x$V1)[[1]], numbasis = 15)
fd3 <- fda.usc::fdata2fd(foo$fdata.est,
                         type.basis = "bspline",
                         nbasis = foo$numbasis.opt)
foo$coef <- t(fd3$coefs)
num.cov <- apply(foo$coef, 1, mean)

### Graph covariate ###

# Normalization of resp+noise to have the connection probability of each graph
new_min <- 0.2
new_max <- 0.8
norm_gen <- function(x){((x-min(x))/(max(x)-min(x)))*(new_max-new_min)+new_min}
conn.prob <- norm_gen(resp + rnorm(length(resp), 0, 7))

# Generation of the graphs
graph.list <- lapply(conn.prob, function(p){sample_gnp(100, p)})

### Functional covariate ###

# Functional data for the first simulation
fun.list <- lapply(data, function(x) x$V1)[[1]]

### Covariates list ###

# Covariates list: functional, graph, numeric
cov.list <- list(fun.list, graph.list, num.cov)


# Model fitting -----------------------------------------------------------

# Number of basis
n.bas <- 15

### Regression Energy Tree ###
etree_fit <- etree(response = resp,
                   covariates = cov.list,
                   case.weights = NULL,
                   minbucket = 5,
                   alpha = 0.8,
                   R = 1000,
                   split.type = 'cluster',
                   coef.split.type = 'test',
                   p.adjust.method = 'fdr',
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
