
# Load packages and functions ---------------------------------------------

library(fda.usc)
library(roahd)
library(energy)
library(entropy)
library(partykit)
library(cluster)
library(igraph)
library(NetworkDistance)
source("functions_v2.R")
source("node_v2.R")
source("split_v2.R")
source("party_v2.R")
source("plot_v2.R")


# Response and covariates lists construction ------------------------------

# Graph covariate
graph.list <- list()
n <- 5 #number of graphs for each class
for(i in 1:n){
  graph.list[[i]] <- sample_gnp(100,0.1)    #type1
  graph.list[[n+i]] <- sample_gnp(100,0.3)  #type2
  graph.list[[2*n+i]] <- sample_gnp(100,0.3)  #type2
}

# Functional covariate
m1 <- matrix(rnorm(200), nrow = 10)
m2 <- matrix(rnorm(100, mean = 2), nrow = 5)
fdata.list <- fdata(rbind(m1, m2))

# Covariates list
cov.list <- list(graph.list, fdata.list)

# Response
resp <- as.factor(c(rep('less_dense', n), rep('more_dense', n), rep('functional', n)))


# Model fitting -----------------------------------------------------------

### CLASSIFICATION ENERGY TREE ###
etree_fit <- etree(response = resp,
                   covariates = cov.list,
                   case.weights = NULL,
                   minbucket = 5,
                   alpha = 0.05,
                   R = 1000,
                   split.type = 'cluster',
                   coef.split.type = 'test')
plot(etree_fit)



# Prediction --------------------------------------------------------------

### ETREE PREDICTION ###

# Prediction
y_pred <- predict(etree_fit)

# Prediction with newdata
graph.list3 <- list()
for(i in 1:n){
  graph.list3[[i]] <- sample_gnp(100,0.1)    #type1
  graph.list3[[n+i]] <- sample_gnp(100,0.12)  #type2
}
graph.list4 <- list()
for(i in 1:n){
  graph.list4[[i]] <- sample_gnp(100,0.1)    #type1
  graph.list4[[n+i]] <- sample_gnp(100,0.12)  #type2
}
new.cov.list <- list(graph.list3, graph.list4)
y_pred2 <- predict(etree_fit, newdata = new.cov.list)

# Error
y <- resp
MEP_etree <- (sum((y-y_pred)^2)/length(y))/(var(y))
