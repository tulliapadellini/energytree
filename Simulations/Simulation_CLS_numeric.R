
# Inizialization ----------------------------------------------------------

# Loading the libraries
library(fda.usc)
library(energy)
library(entropy)
library(partykit)
library(mlr)
library(fastDummies)
source("functions.R")

# Loading the dataset
data('iris')

# Errors
ACC_etree <- c()
ACC_mlc <- c()


# Response and covariates lists construction ------------------------------

# Response
resp <- iris$Species

# Covariates
cov.list <- list('Sepal.Length' = iris$Sepal.Length, 'Sepal.Width' = iris$Sepal.Width, 'Petal.Length' = iris$Petal.Length, 'Petal.Width' = iris$Petal.Width)


# Model fitting -----------------------------------------------------------

### CLASSIFICATION ENERGY TREE ###
etree_fit <- etree(response = resp,
                   covariates = cov.list,
                   case.weights = NULL,
                   minbucket = 1,
                   alpha = 0.05,
                   R = 1000)
plot(etree_fit)

### MULTILABEL CLASSIFICATION VIA CLASSIFIER CHAINS (MLR PACKAGE) ###
# Transforming response into logical dummy variables (required by makeMultilabelTask)
resp_dummy <- dummy_cols(resp)[,-1]
resp_dummy <- sapply(resp_dummy, as.logical)
colnames(resp_dummy) <- c('setosa', 'versicolor', 'virginica')
# All data (coariates and dummies for the response) together
all_data <- cbind(as.data.frame(do.call(cbind, cov.list)), resp_dummy)
# Multilabel infrastucture (i.e. setting data and target)
iris_task <- makeMultilabelTask(data = all_data, target = colnames(resp_dummy))
# Base Learner
binary.learner = makeLearner("classif.rpart")
# Multilabel learner (wrapper of repetitions of the base one)
mcc = makeMultilabelClassifierChainsWrapper(binary.learner)
# Train and test sets
n = getTaskSize(iris_task)
train_set = sample(1:n, round(n*0.8))
test_set = (1:n)[!(1:n) %in% train_set]
# Train the multilabel learner
iris_train = train(mcc, iris_task, subset = train_set)
# Prediction
iris_pred = predict(iris_train, task = iris_task, subset = test_set)
# Accuracy
performance(iris_pred, measures = list(multilabel.acc))

# Prediction --------------------------------------------------------------

### ETREE CLASSIFICATION PREDICTION ###

# Prediction
y_pred <- predict(etree_fit)

# Error
y <- resp
t <- table(y_pred, y)
ACC_etree <- sum(diag(t))/(length(y))


### MULTILABEL CLASSIFICATION VIA CLASSIFIER CHAINS PREDICTION ###
ACC_mlc <- performance(iris_pred, measures = list(multilabel.acc))


# Storing the results
save(ACC_etree, ACC_mlc, file = "results.RData")



