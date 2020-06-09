

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
source("NKI_data_import.R")


# Import data -------------------------------------------------------------

nki <- generate_dataset(data_folder = 'NKI_Rockland/',
                        y_filename = 'NKI_clinical_information.txt',
                        y_column = 'WASI_FULL_4',
                        output_filename = 'NKIdata.RData',
                        output_folder = ".",
                        ext_save = FALSE)



# Dataset construction ----------------------------------------------------

# Response
resp <- nki$y

# Covariates list
cov.list <- list(lapply(nki$structural, function(g) igraph::graph_from_adjacency_matrix(g, weighted = T)),
                 lapply(nki$functional, function(g) igraph::graph_from_adjacency_matrix(g, weighted = T)))


# Energy Tree fit ---------------------------------------------------------

# Fit
set.seed(2948)
etree_fit <- etree(response = resp,
                   covariates = cov.list,
                   case.weights = NULL,
                   minbucket = 5,
                   alpha = 0.5,
                   R = 1000,
                   split.type = 'cluster',
                   coef.split.type = 'test')

# Plot
plot(etree_fit)

# Fitted values
y_fitted <- predict(etree_fit)

# Mean Error Prediction
(MEP_etree <- (sum((resp-y_fitted)^2)/length(resp))/(var(resp)))

# Root Mean Square Error
(MEP_etree <- sqrt(sum((resp-y_fitted)^2)/length(resp)))

# Mean Square Percentage Error
(MEP_etree <- sum(((resp-y_fitted)/resp)^2)/length(resp))


# Prediction --------------------------------------------------------------

# Predicted values
y_pred <- predict(etree_fit, newdata = cov.list)



##### NEW #####

# Mixed covariates -----------------------------------------------------------

# Response
resp <- nki$y

# Graph covariates
nki_str <- lapply(nki$structural,
                  function(g) igraph::graph_from_adjacency_matrix(g, weighted = T))
nki_fun <- lapply(nki$functional,
                  function(g) igraph::graph_from_adjacency_matrix(g, weighted = T))

# Import clinical information
NKI_clinical <- read.csv("NKI_clinical_information.txt")

# Select individuals who also have structural and functional matrices
sel_id <- names(nki$structural)
sel_idx <- (NKI_clinical$Subject %in% sel_id)
#remark: already done for str, fun and wasi coming from nki

# Covariates list
cov.list <- list(SCM = nki_str,
                 FCM = nki_fun,
                 Gender = NKI_clinical$Gender[sel_idx],
                 Age = NKI_clinical$Age[sel_idx])

# Fit
set.seed(2948)
etree_fit <- etree(response = resp,
                   covariates = cov.list,
                   case.weights = NULL,
                   minbucket = 5,
                   alpha = 0.5,
                   R = 1000,
                   split.type = 'cluster',
                   coef.split.type = 'test')

# Plot
plot(etree_fit)

# Fitted values
y_fitted <- predict(etree_fit)

# Mean Error Prediction
(MEP_etree <- (sum((resp-y_fitted)^2)/length(resp))/(var(resp)))

# Root Mean Square Error
(RMSE_etree <- sqrt(sum((resp-y_fitted)^2)/length(resp)))

# Mean Square Percentage Error
(MSPE_etree <- sum(((resp-y_fitted)/resp)^2)/length(resp))

# Predicted values
y_pred <- predict(etree_fit, newdata = cov.list)


# CV evaluation ---------------------------------------------------------------

# Loading caret (Classification And REgression Training) package
library(caret)

# Folds
set.seed(12345)
folds <- caret::createFolds(resp, k = 10, list = TRUE, returnTrain = FALSE)

# Cross-validation: fit without fold and predict on fold (x10)
e_cv <- lapply(folds,
               function(f){
                 etree_fit <- etree(response = resp[-f],
                                    covariates = lapply(cov.list,
                                                        function(c){
                                                          c[-f]
                                                        }),
                                    case.weights = NULL,
                                    minbucket = 5,
                                    alpha = 0.5,
                                    R = 1000,
                                    split.type = 'cluster',
                                    coef.split.type = 'test')
                 y_pred <- predict(etree_fit, newdata = lapply(cov.list,
                                                               function(c){
                                                                 c[f]
                                                               }))
                 list(e_fit = etree_fit,
                      e_pred = y_pred,
                      e_resp = resp[f],
                      e_summ = defaultSummary(data.frame(obs = resp[f],
                                                         pred = y_pred)))
               })

# Summary performance measures
(e_cv_summ <- sapply(e_cv, function(e) e$e_summ))
