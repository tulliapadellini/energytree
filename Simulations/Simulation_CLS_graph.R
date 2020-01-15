
# Inizialization ----------------------------------------------------------

# Loading the libraries
library(energy)
library(partykit)
library(cluster)
library(igraph)
library(NetworkDistance)
source("functions.R")

# Error(s)
ACC_etree <- c()


# Response and covariates lists construction ------------------------------

# Graph simulations
graph.list <- list()
n <- 5 #number of graphs of type 1 & number of graphs of type 2
for(i in 1:n){
  graph.list[[i]] <- sample_gnp(100,0.1)    #type1
  graph.list[[n+i]] <- sample_gnp(100,0.12)  #type2
}

# Response
resp <- as.factor(c(rep('less_dense', n), rep('more_dense', n)))

# Only one covariate
cov.list <- list(graph.list)

# Two covariates
graph.list2 <- list()
for(i in 1:n){
  graph.list2[[i]] <- sample_gnp(100,0.1)    #type1
  graph.list2[[n+i]] <- sample_gnp(100,0.12)  #type2
}
cov.list <- list(graph.list, graph.list2)


# Model fitting -----------------------------------------------------------

### REGRESSION ENERGY TREE ###
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
