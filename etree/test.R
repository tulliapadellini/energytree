library(etree)

data("sim_reg")

resp <- lapply(sim_reg, function(x) x$Y)[[1]]

# Two functional predictors
cov.list <- list(lapply(sim_reg, function(x) x$V1)[[1]], lapply(sim_reg, function(x) x$V1)[[2]])


etree_fit <- etree(response = resp,
                   covariates = cov.list,
                   case.weights = NULL,
                   minbucket = 5,
                   alpha = 0.05,
                   R = 1000,
                   split.type = 'cluster',
                   coef.split.type = 'test')

y_pred = predict(etree_fit)
y_pred2 = predict(etree_fit, newdata = cov.list, split.type = "cluster", nb = 5)
y_pred3 = predict(etree_fit, newdata = cov.list, split.type = "cluster", nb = 5)
y_pred - y_pred2
y_pred3 - y_pred2

mean((resp-y_pred)^2)
mean((resp-y_pred2)^2)




# graphs ------------------------------------------------------------------

graph.list <- list()
n <- 5 #number of graphs of type 1 & number of graphs of type 2
for(i in 1:n){
  graph.list[[i]] <- sample_gnp(100,0.1)    #type1
  graph.list[[n+i]] <- sample_gnp(100,0.12)  #type2
}

# Response

# Only one covariate
cov.list <- list(graph.list)


resp <- sapply(graph.list, ecount) #number of edges in each graph

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

y_pred = predict(etree_fit)
y_pred2 = predict(etree_fit, newdata = cov.list, split.type = "cluster", nb = 5)
y_pred - y_pred2

mean((resp-y_pred)^2)
mean((resp-y_pred2)^2)



# diagrams ----------------------------------------------------------------


diag.list <- list()
n <- 5 #number of graphs of type 1 & number of graphs of type 2
for(i in 1:n){
  kk = TDA::circleUnif(50)    #type1
  diag.list[[i]] <- TDA::ripsDiag(kk, 1,1)
  kk2 = rbind(TDA::circleUnif(50), TDA::circleUnif(50)+2)
  diag.list[[n+i]] <- TDA::ripsDiag(kk2, 1,1)  #type2
}

resp <- rep(c(1,2), c(5,5))

# Two functional predictors
cov.list <- list(diag.list)

etree_fit <- etree(response = resp,
                   covariates = cov.list,
                   case.weights = NULL,
                   minbucket = 5,
                   alpha = 0.05,
                   R = 1000,
                   split.type = 'cluster',
                   coef.split.type = 'test')

y_pred = predict(etree_fit)
y_pred2 = predict(etree_fit, newdata = cov.list, split.type = "cluster", nb = 5)
y_pred3 = predict(etree_fit, newdata = cov.list, split.type = "cluster", nb = 5)
y_pred - y_pred2
y_pred3 - y_pred2

mean((resp-y_pred)^2)
mean((resp-y_pred2)^2)


