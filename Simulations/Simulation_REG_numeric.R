
# Inizialization ----------------------------------------------------------

# Loading the libraries
library(energy)
library(entropy)
library(partykit)
library(AppliedPredictiveModeling)
library(gdata)

# Loading the dataset
data("abalone")

# Errors
MEP_etree <- c()
MEP_lm <- c()


# Response and covariates lists construction ------------------------------

# Less observations for the moment
data <- abalone[1:100,]

# Response
resp <- data$Rings

# Covariates
x1 <- data[,3]
x2 <- data[,4]
x3 <- data[,5]
x4 <- data[,6]
x5 <- data[,7]
x6 <- data[,8]
cov.list <- list(x1, x2, x3, x4, x5, x6)


# Model fitting -----------------------------------------------------------

### REGRESSION ENERGY TREE ###
etree_fit <- etree(response = resp,
                   covariates = cov.list,
                   case.weights = NULL,
                   minbucket = 1,
                   alpha = 0.05,
                   R = 1000)
plot(etree_fit)

### LINEAR REGRESSION ###
lm_fit <- lm (resp ~ do.call(cbind, cov.list))


# Prediction --------------------------------------------------------------

### ETREE PREDICTION ###

# Prediction
y_pred <- predict(etree_fit)

# Error
y <- resp
MEP_etree <- (sum((y-y_pred)^2)/length(y))/(var(y))


### LINEAR MODEL PREDICTION ###
y_pred <- predict(lm_fit)
MEP_lm <- (sum((y-y_pred)^2)/length(y))/(var(y))

