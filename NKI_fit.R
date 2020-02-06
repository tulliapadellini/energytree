

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
cov.list <- list(lapply(nki$structural, function(g) graph_from_adjacency_matrix(g, weighted = T)),
                 lapply(nki$functional, function(g) graph_from_adjacency_matrix(g, weighted = T)))


# Energy Tree fit ---------------------------------------------------------

# Fit
etree_fit <- etree(response = resp,
                   covariates = cov.list,
                   case.weights = NULL,
                   minbucket = 5,
                   alpha = 0.6,
                   R = 1000,
                   split.type = 'cluster',
                   coef.split.type = 'test')

# Plot
plot(etree_fit)

# Fitted values
y_fitted <- predict(etree_fit)

# Error
(MEP_etree <- (sum((resp-y_fitted)^2)/length(resp))/(var(resp)))


# Prediction --------------------------------------------------------------

# Predicted values
y_pred <- predict(etree_fit, newdata = cov.list)
