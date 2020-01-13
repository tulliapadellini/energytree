
# Inizialization ----------------------------------------------------------

# Loading the libraries
library(energy)
library(partykit)
library(cluster)
library(igraph)
library(NetworkDistance)
#source("functions.R")

# Error(s)
MEP_etree <- c()


# Response and covariates lists construction ------------------------------

# Graph simulations
graph.list <- list()
n <- 5 #number of graphs of type 1 & number of graphs of type 2
for(i in 1:n){
  graph.list[[i]] <- sample_gnp(100,0.1)    #type1
  graph.list[[n+i]] <- sample_gnp(100,0.12)  #type2
}

# Response
resp <- sapply(graph.list, ecount) #number of edges in each graph

# Only one covariate
cov.list <- list(graph.list)

# Two covariates
graph.list2 <- lapply(resp, function(n.edges){sample_gnm(100, n.edges)})
# so that we have different graphs, with the same number of edges (i.e. resp) for i=1,...,n
cov.list <- list(graph.list, graph.list2)


# Model fitting -----------------------------------------------------------

### REGRESSION ENERGY TREE ###
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
graph.list3 <- lapply(resp, function(n.edges){sample_gnm(100, n.edges)})
graph.list4 <- lapply(resp, function(n.edges){sample_gnm(100, n.edges)})
new.cov.list <- list(graph.list3, graph.list4)
y_pred <- predict(etree_fit, newdata = new.cov.list)


# Error
y <- resp
MEP_etree <- (sum((y-y_pred)^2)/length(y))/(var(y))
