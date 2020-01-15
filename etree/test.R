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
<<<<<<< HEAD

mean((resp-y_pred)^2)
mean((resp-y_pred2)^2)
=======

mean((resp-y_pred)^2)
mean((resp-y_pred2)^2)


new.cov.list <- list(lapply(sim_reg, function(x) x$V1)[[51]], lapply(sim_reg, function(x) x$V1)[[52]])
y_pred4 = predict(etree_fit, newdata = new.cov.list, split.type = "cluster", nb = 5)

>>>>>>> master
