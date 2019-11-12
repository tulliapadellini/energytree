
# Inizialization ----------------------------------------------------------

# Loading the libraries
library(energy)
library(partykit)
library(cluster)
library(igraph)
library(NetworkDistance)
source("functions.R")

# Error(s)
MEP_etree <- c()


# Response and covariates lists construction ------------------------------

# Graph simulations
graph.list <- list()
n <- 5 #number of graphs of type 1 & number of graphs of type 2
for(i in 1:n){
  graph.list[[i]] <- sample_gnp(100,0.1)    #type1
  graph.list[[n+i]] <- sample_gnp(100,0.15)  #type2
}

# Response
resp <- sapply(graph.list, ecount) #number of edges in each graph

# Covariate
cov.list <- list(graph.list)

# Model fitting -----------------------------------------------------------

### REGRESSION ENERGY TREE ###
etree_fit <- etree(response = resp,
                   covariates = cov.list,
                   case.weights = NULL,
                   minbucket = 5,
                   alpha = 0.05,
                   R = 1000,
                   split.type = 'cluster',
                   coef.split.type = 'test')
plot(etree_fit)



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


# Storing the results
save(MEP_etree, file = "results.RData")

