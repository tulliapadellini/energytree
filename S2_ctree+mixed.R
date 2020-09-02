### S2 - Extension of simulations from Hothorn et al. (2006) ###
###      to the mixed and structured case ###

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
etree_sim_cov <- list()
for(i in 1:n_sim){
  set.seed(i)
  # Numeric: uniform between 0 and 1
  x1 <- runif(n_obs, min = 0, max = 1)
  # Numeric: uniform between 0 and 1, with values rounded to one digit
  x2 <- round(runif(n_obs, min = 0, max = 1),
              1)
  # Nominal: defines the two classes
  x3 <- factor(c(rep(0, n_obs/2),
                 rep(1, n_obs/2)))
  # Nominal: defines two subclasses for each class
  #to be replaced with a covariate with NAs?
  x4 <- factor(c(rep(1, n_obs/4),
                 rep(2, n_obs/4),
                 rep(3, n_obs/4),
                 rep(4, n_obs/4)))
  # Graph: Erdos-Renyi (1959) model with same connection prob
  x5 <- lapply(1:n_obs, function(j) igraph::sample_gnp(100, 0.2))
  # Graph: Small-World model (Watts-Strogatz, 1998) with same rewiring prob
  x6 <- list(lapply(1:n_obs, function(j) igraph::sample_smallworld(1, 100, 5, 0.1)))
  # Graph: Preferential Attachment model (Barabasi-Albert, 1999)
  x7 <- list(lapply(1:n_obs, function(j) igraph::sample_pa(100)))
  # Functional: Ornstein-Uhlenbeck random process
  x8 <- fda.usc::rproc2fdata(n_obs, seq(0, 1, len = 20), sigma = 'OrnsteinUhlenbeck')
  # Functional: Wiener random process
  x9 <- fda.usc::rproc2fdata(n_obs, seq(0, 1, len = 20), sigma = 'brownian')
  # Covariates list for this simulation
  etree_sim_cov[[i]] <- list(x1 = x1, x2 = x2,
                         x3 = x3, x4 = x4,
                         x5 = x5, x6 = x6, x7 = x7,
                         x8 = x8, x9 = x9)
  print(i)
}

# Response simulations
etree_sim_resp <- list()
for(i in 1:n_sim){
  set.seed(i)
  etree_sim_resp[[i]] <- c(rnorm(n_obs/2, mean = 0, sd = 1),
                           rnorm(n_obs/2, mean = mu, sd = 1))
}


# Variable selection probabilities (under independence) -----------------------

# First variable selection
source('etree_sel.R') #starts to fit an etree, but only return the first var's id
library(pbapply)
etree_first_var <- pbmapply(function(covs, resp){
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
etree_sim_cov,
etree_sim_resp,
SIMPLIFY = TRUE)

# Frequency plot
barplot(table(etree_first_var))

# Proportions
prop.table(table(etree_first_var))
