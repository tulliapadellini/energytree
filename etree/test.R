library(etree)

data("sim_reg")

resp <- lapply(sim_reg, function(x) x$Y)[[1]]

# Two functional predictors
cov.list <- list(lapply(sim_reg, function(x) x$V1)[[1]], lapply(sim_reg, function(x) x$V1)[[2]])
cov.list_small = list(lapply(sim_reg, function(x) x$V1)[[1]][1:10, ], lapply(sim_reg, function(x) x$V1)[[2]][1:10,])

plot(cov.list_small[[1]])

etree_fit <- etree(response = resp,
                   covariates = cov.list,
                   case.weights = NULL,
                   minbucket = 5,
                   alpha = 0.05,
                   R = 1000,
                   split.type = 'cluster',
                   coef.split.type = 'test')

y_pred = predict(etree_fit)
y_pred2 = predict(etree_fit, newdata = cov.list_small, split.type = "cluster", nb = 5)
y_pred3 = predict(etree_fit, newdata = cov.list, split.type = "cluster", nb = 5)
y_pred - y_pred3
y_pred3 - y_pred2

mean((resp-y_pred)^2)
mean((resp-y_pred3)^2)




# Graph -------------------------------------------------------------------

graph.list <- list()
n <- 10 #number of graphs of type 1 & number of graphs of type 2
for(i in 1:n){
  graph.list[[i]] <- sample_gnp(100,0.1)    #type1
  graph.list[[n+i]] <- sample_gnp(100,0.12)  #type2
}

# Response
resp <- sapply(graph.list, ecount) #number of edges in each graph

cov.list <- list(graph.list)


etree_fit <- etree(response = resp,
                   covariates = cov.list,
                   case.weights = NULL,
                   minbucket = 5,
                   alpha = 0.05,
                   R = 1000,
                   split.type = 'cluster',
                   coef.split.type = 'test')

plot(etree_fit)


y_pred = predict(etree_fit)
y_pred3 = predict(etree_fit, newdata = cov.list, split.type = "cluster", nb = 5)
y_pred - y_pred3
y_pred3 - y_pred2



