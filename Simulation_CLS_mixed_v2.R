
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
load("sim_data.RData")


# Response and covariates -----------------------------------------------------

# Restriction on the observations' number (for computational reasons)
obs <- sample(1:798, 150)

# Response
resp <- lapply(data, function(x) x$cls[obs])[[1]]
#remark: we take the first element since the dataset contains 100 simulations

### Numerical covariate ###

# One of the basis of another simulation of the same dataset
foo <- fda.usc::optim.basis(lapply(data, function(x) fdata(x[obs,2:129]))[[2]], numbasis = 15)
fd3 <- fda.usc::fdata2fd(foo$fdata.est,
                         type.basis = "bspline",
                         nbasis = foo$numbasis.opt)
foo$coef <- t(fd3$coefs)
num.cov <- apply(foo$coef, 1, mean)

### Graph covariate ###

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

### Functional covariate ###

# Functional data for the first simulation
fun.list <- lapply(data, function(x) fdata(x[obs,2:129]))[[1]]

### Covariates list ###

# Covariates list: functional, graph, numeric
cov.list <- list(fun.list, graph.list, num.cov)


# Model fitting -----------------------------------------------------------

# Number of basis
n.bas <- 15

### Classification Energy Tree ###
set.seed(123)
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


