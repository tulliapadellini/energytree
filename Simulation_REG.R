
# Inizialization ----------------------------------------------------------

# Loading the libraries
library(fda.usc)
library(roahd)
library(energy)
library(entropy)
library(partykit)

# Loading the dataset
load("Regression_Sim_Dataset.RData")

# Errors
MEP_etree <- c()
MEP_funbasis <- c()
MEP_funpc <- c()


# Response and covariates lists construction ------------------------------

# Response
resp <- lapply(data, function(x) x$Y)[[1]]
#remark: we take the first element since the dataset contains 100 simulations

# Regression with one functional predictor
cov.list <- lapply(data, function(x) x$V1)[1]

# Regression with many functional predictors
cov.list <- list(lapply(data, function(x) x$V1)[[1]], lapply(data, function(x) x$V1)[[2]])

# Regression with a numeric predictor and a functional one (three alternatives)

#1) numerical variable: random
cov.list <- list(lapply(data, function(x) x$V1)[[1]], rnorm(200, 10, 1))
#ok, ma non splitta (quasi mai) sulla numerica

#2) numerical variable: response+noise
cov.list <- list(lapply(data, function(x) x$V1)[[1]], resp+rnorm(200, 0, 10))
#ok

#3) numerical variable: one of the basis of another simulation of the same dataset
foo <- fda.usc::min.basis(lapply(data, function(x) x$V1)[[2]], numbasis = 15)
fd3 <- fda.usc::fdata2fd(foo$fdata.est,
                         type.basis = "bspline",
                         nbasis = foo$numbasis.opt)
foo$coef <- t(fd3$coefs)
cov.list <- list(lapply(data, function(x) x$V1)[[1]], lapply(data, function(x) x$V1)[[2]], foo$coef[,15])
#ok


# Model fitting -----------------------------------------------------------

# Number of basis
n.bas <- 15

### REGRESSION ENERGY TREE ###
etree_fit <- etree(response = resp,
                   covariates = cov.list,
                   case.weights = NULL,
                   minbucket = 5,
                   alpha = 0.05,
                   R = 1000,
                   split.type = 'coeff',
                   coef.split.type = 'test',
                   nb = n.bas)
plot(etree_fit)


### FUNCTIONAL LINEAR MODEL WITH BASIS ###

# Covariate
x <- cov.list[[1]]   #remark: here the covariate must be functional

# Response
y <- resp

# Time points of the functional covariate
tt <- x[["argvals"]]

# B-spline basis creation
basis1 <- create.bspline.basis(rangeval = range(tt), nbasis = n.bas)

# Functional regression with scalar response using basis representation
funbasis_fit <- fregre.basis(x, y, basis.x = basis1)


### FUNCTIONAL LINEAR MODEL WITH PRINCIPAL COMPONENTS ###

# Functional regression with scalar response using Principal Components Analysis
funpc_fit <- fregre.pc(x, y, kmax = 7)



# Prediction --------------------------------------------------------------

### MYREG PREDICTION ###

# New covariates
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


### FLMB PREDICTION ###
y_pred1 <- predict.fregre.fd(funbasis_fit, new.fdataobj = cov.list[[1]])
MEP_funbasis <- (sum((y-y_pred1)^2)/length(y))/(var(y))


### FLMPC PREDICTION ###
y_pred2 <- predict.fregre.fd(funpc_fit, cov.list[[1]])
MEP_funpc <- (sum((y-y_pred2)^2)/length(y))/(var(y))


# Storing the results
save(MEP_etree, MEP_funbasis, MEP_funpc, file = "results.RData")

