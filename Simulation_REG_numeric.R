
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
                   R = 1000,
                   rnd.sel = T,
                   rnd.splt = TRUE,
                   nb = 15)
plot(etree_fit)

### LINEAR REGRESSION ###
lm_fit <- lm (resp ~ do.call(cbind, cov.list))


# Prediction --------------------------------------------------------------

### MYREG PREDICTION ###
new.cov.list = lapply(cov.list, function(j){
  if(class(j) == 'fdata'){

    foo <- fda.usc::min.basis(j, numbasis = n.bas)
    fd3 <- fda.usc::fdata2fd(foo$fdata.est,
                             type.basis = "bspline",
                             nbasis = foo$numbasis.opt)
    foo$coef <- t(fd3$coefs)
    return(foo$coef)

  } else if(class(j) == 'list' &
            all(sapply(j, class) == 'igraph')){

      shell <- graph.to.shellness.distr.df(j)
      return(shell)

    } else {

      return(j)

    }
}
)

# New covariates dataframe
new.cov.df <- as.data.frame(do.call(cbind, new.cov.list))
names(new.cov.df) <- 1:ncol(new.cov.df)

# Prediction
y_pred <- predict(etree_fit, newdata = new.cov.df)

# Error
y <- resp
MEP_etree <- (sum((y-y_pred)^2)/length(y))/(var(y))


### LINEAR MODEL PREDICTION ###
y_pred <- predict(lm_fit)
MEP_lm <- (sum((y-y_pred)^2)/length(y))/(var(y))

