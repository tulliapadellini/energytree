
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
ACC_ftree <- c()


# Response and covariates lists construction ------------------------------

# Restriction on the observations' number (for computational reasons)
obs <- sample(1:798, 150)

# Response
resp <- lapply(data, function(x) x$cls[obs])[[1]]
#remark: we take the first element since the dataset contains 100 simulations

# Only one functional predictor
cov.list <- lapply(data, function(x) fdata(x[obs,2:129]))[1]

# Two functional predictors
cov.list <- list(lapply(data, function(x) fdata(x[obs,2:129]))[[1]],
                 lapply(data, function(x) fdata(x[obs,2:129]))[[2]])


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
                   split.type = 'cluster',
                   coef.split.type = 'test')
plot(etree_fit)

### FUNCTIONAL CLASSIFICATION TREE ###

# One functional covariate
list_all <- list("df" = as.data.frame(resp), "x" = cov.list[[1]])
fvb <- create.fdata.basis(cov.list[[1]],l=1:15)
ct_fit <- classif.tree(resp ~ x, data = list_all, basis.b = list(x =fvb))

# # più covariate funzionali
# list_all <- list("df" = as.data.frame(cbind(resp, cov.list[[3]])), "x" = list(cov.list[[1]], cov.list[[2]]))
# fvb1 <- create.fdata.basis(cov.list[[1]],l=1:15)
# a2 <- classif.tree(resp ~ x, data = list_all, basis.b = list(x1 = fvb1, x2 = fvb2))
# # errore: può essere che si possa utilizzare al massimo una covariata funzionale?
#
# # una coavariata funzionale e una numerica
# list_all <- list("df" = as.data.frame(cbind(resp, cov.list[[3]])), "x" = cov.list[[1]])
# fvb <- create.fdata.basis(cov.list[[1]],l=1:15)
# a2 <- classif.tree(resp ~ list_all$df$V2 + x, data = list_all, basis.b = list(x =fvb))
# predict(a2)
# t <- table(predict(a2))
# acc_m1[i] <- sum(diag(t))/(length(y))
# # errore: forse non funziona neanche questo caso?

# Prediction --------------------------------------------------------------

### ETREE CLASSIFICATION PREDICTION ###

# Prediction
y_pred <- predict(etree_fit)

# Prediction with newdata
new.cov.list <- list(lapply(data, function(x) fdata(x[obs,2:129]))[[3]],
                     lapply(data, function(x) fdata(x[obs,2:129]))[[4]])
y_pred2 <- predict(etree_fit, newdata = new.cov.list)


# Error
y <- resp
t <- table(y_pred, y)
ACC_etree <- sum(diag(t))/(length(y))


### FUNCTIONAL CLASSIFICATION TREE PREDICTION ###
t <- table(predict(ct_fit))
ACC_ftree <- sum(diag(t))/(length(y))


# Storing the results
save(ACC_etree, ACC_ftree, file = "results.RData")



