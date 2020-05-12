

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
