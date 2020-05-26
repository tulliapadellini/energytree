
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


# Initialization and dataset construction -------------------------------------

# Error
MEP_etree <- c()

# Data import
nki <- read.csv("~/PhD/NIK-ULTRA/energytree/NKI_clinical_information.txt")

# Reduced dataset
nki_red <- nki[c('Subject', 'Gender', 'Group', 'Hand', 'Operator', 'WASI_FULL_4')]

# Exclude NAs
nki2 <- nki_red[!is.na(nki_red$WASI_FULL_4),]

# Response
resp <- nki2$WASI_FULL_4

# Covariates
cov.list <- list(gender = nki2$Gender)


# More than one factor covariates --------------------------------------------

# Covariates
cov.list <- list(gender = nki2$Gender,
                 group = nki2$Group,
                 hand = nki2$Hand,
                 operator = nki2$Operator)
#in addition to having more than one covariate, also the number of levels is
#generalized: gender has 2, group has 11, hand has 4, operator has 6


# Model fitting -----------------------------------------------------------

### REGRESSION ENERGY TREE ###
etree_fit <- etree(response = resp,
                   covariates = cov.list,
                   case.weights = NULL,
                   minbucket = 5,
                   alpha = 0.5,
                   R = 1000,
                   split.type = 'coeff',
                   coef.split.type = 'test')
plot(etree_fit)


# Prediction --------------------------------------------------------------

### ETREE PREDICTION ###

# Prediction
y_pred <- predict(etree_fit)

# Prediction with newdata (no meaningful covariates, just to see if it works)
new.cov.list <- list(gender = sample(nki2$Gender, replace = T),
                     group = sample(nki2$Group, replace = T),
                     hand = sample(nki2$Hand, replace = TRUE),
                     operator = sample(nki2$Operator, replace = TRUE))
y_pred2 <- predict(etree_fit, newdata = new.cov.list)


# Error
y <- resp
MEP_etree <- (sum((y-y_pred)^2)/length(y))/(var(y))
