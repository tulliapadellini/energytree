
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
resp <- nki2$Gender

# Covariates
cov.list <- list(gender = nki2$Gender)


# More than one factor covariates --------------------------------------------

# Covariates
cov.list <- list(#gender = nki2$Gender,
                 group = nki2$Group,
                 hand = nki2$Hand,
                 operator = nki2$Operator)
#in addition to having more than one covariate, also the number of levels is
#generalized: gender has 2, group has 11, hand has 4, operator has 6


# Real dataset: Titanic ------------------------------------------------------
## same as that used by constparty vignette ##

# Import dataset
data("Titanic", package = "datasets")
ttnc <- as.data.frame(Titanic)
ttnc <- ttnc[rep(1:nrow(ttnc), ttnc$Freq), 1:4]
names(ttnc)[2] <- "Gender"

# Reduce dataset (for computational reasons)
set.seed(12345)
nobs <- 100
newobs <- sample(1:(dim(ttnc)[1]), size = nobs, replace = FALSE)
ttnc_red <- ttnc[newobs,]

# Response
resp <- ttnc_red$Survived

# Covariates
cov.list <- list(Class = ttnc_red$Class,
                 Gender = ttnc_red$Gender,
                 Age = ttnc_red$Age)


# Model fitting -----------------------------------------------------------

### CLASSIFICATION ENERGY TREE ###
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

# Prediction with newdata (no meaningful covariates, just to see if it works)
new.cov.list <- list(gender = sample(nki2$Gender, replace = T),
                     group = sample(nki2$Group, replace = T),
                     hand = sample(nki2$Hand, replace = TRUE),
                     operator = sample(nki2$Operator, replace = TRUE))
y_pred2 <- predict(etree_fit, newdata = new.cov.list)


# Error
y <- resp
t <- table(y_pred, y)
ACC_etree <- sum(diag(t))/(length(y))
