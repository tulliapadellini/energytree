### S1 - Simulations based on Hothorn et al. (2006) ###

# Load packages and functions -------------------------------------------------

# Packages
library(fda.usc)
library(energy)
library(entropy)
library(partykit)
library(cluster)
library(igraph)
library(NetworkDistance)
library(ggparty)
library(checkmate)

# Functions
source("functions_v2.R")
source("node_v2.R")
source("split_v2.R")
source("party_v2.R")
source("plot_v2.R")
source("ggparty_v1.R")
source("get_plot_data_v1.R")


# Data simulation -------------------------------------------------------------

# Number of observations
n_obs <- 100

# Number of simulations
n_sim <- 10000

# Mean of the normal distribution by which obs of the second class are modelled
mu = 0

# Covariates simulations
ctree_sim <- list()
for(i in 1:n_sim){
  set.seed(i)
  # Uniform between 0 and 1
  x1 <- runif(n_obs, min = 0, max = 1)
  # Uniform between 0 and 1
  x2 <- runif(n_obs, min = 0, max = 1)
  # Uniform between 0 and 1
  x3 <- runif(n_obs, min = 0, max = 1)
  ## Uniform between 0 and 1, with 25% of the values being missing at random
  #x4 <- replace(runif(n_obs, min = 0, max = 1),
  #              sample(1:n_obs, n_obs/4, replace = FALSE),
  #              NA)
  # Uniform between 0 and 1, with values rounded to one digit
  x5 <- round(runif(n_obs, min = 0, max = 1),
              1)
  # Covariate that defines the two classes
  x6 <- c(rep(0, n_obs/2),
          rep(1, n_obs/2))
  # Covariates list for this simulation
  ctree_sim[[i]] <- list(x1 = x1, x2 = x2, x3 = x3, x5 = x5, x6 = x6)
}

# Response simulations
ctree_sim_resp <- list()
for(i in 1:n_sim){
  set.seed(i)
  ctree_sim_resp[[i]] <- c(rnorm(n_obs/2, mean = 0, sd = 1),
                           rnorm(n_obs/2, mean = mu, sd = 1))
}


# Model fitting ---------------------------------------------------------------

### Classification Energy Tree ###
etree_sim <- mapply(function(covs, resp){
  set.seed(123)
  etree(response = resp,
        covariates = covs,
        case.weights = NULL,
        minbucket = 1,
        alpha = 0.05,
        R = 1000,
        split.type = 'cluster',
        coef.split.type = 'test',
        p.adjust.method = 'fdr')
},
ctree_sim,
ctree_sim_resp)



# Variable selection probabilities (under independence) -----------------------

# First variable selection
source('etree_sel.R') #starts to fit an etree, but only return the first var's id
etree_first_var <- mapply(function(covs, resp){
  set.seed(123)
  etree_sel(response = resp,
            covariates = covs,
            case.weights = NULL,
            minbucket = 1,
            alpha = 1,
            R = 1000,
            split.type = 'cluster',
            coef.split.type = 'test',
            p.adjust.method = 'fdr')
},
ctree_sim,
ctree_sim_resp,
SIMPLIFY = TRUE)
