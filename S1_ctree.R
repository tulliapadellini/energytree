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
library(pbapply)

# Functions
source("functions_v2.R")
source("node_v2.R")
source("split_v2.R")
source("party_v2.R")
source("plot_v2.R")
source("ggparty_v1.R")
source("get_plot_data_v1.R")
source('etree_sel.R') #starts to fit an etree, but only return the first var's id


# Data simulation -------------------------------------------------------------

# Number of observations
n_obs <- 100

# Number of simulations
n_sim <- 10000

# Mean of the normal distribution by which obs of the second class are modelled
mu = 0

# Covariates simulations
ctree_sim_cov <- list()
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
  ctree_sim_cov[[i]] <- list(x1 = x1, x2 = x2, x3 = x3, x5 = x5, x6 = x6)
}

# Response simulations
ctree_sim_resp <- list()
for(i in 1:n_sim){
  set.seed(i)
  ctree_sim_resp[[i]] <- c(rnorm(n_obs/2, mean = 0, sd = 1),
                           rnorm(n_obs/2, mean = mu, sd = 1))
}


# Variable selection probabilities (under independence) -----------------------

# First-variable selection
ctree_first_var <- pbmapply(resp = ctree_sim_resp,
                            covs = ctree_sim_cov,
                            function(resp, covs){
                              set.seed(123)
                              etree_sel(response = resp,
                                        covariates = covs,
                                        case.weights = NULL,
                                        minbucket = 10,
                                        alpha = 1,
                                        R = 1000,
                                        split.type = 'cluster',
                                        coef.split.type = 'test',
                                        p.adjust.method = 'fdr')
                            },
                            SIMPLIFY = TRUE)

# Frequency plot
barplot(table(ctree_first_var))

# Proportions
prop.table(table(ctree_first_var))
xtable(t(prop.table(table(ctree_first_var))*100))
xtable(prop.table(table(ctree_first_var)), digits = 4)


# Power and conditional probability (growing association) ---------------------

# Grid values for mu
mu_grid = seq(0, 1, 0.05)

# etree fit: n_sim simulations for each value of mu
ctree_mu_grid <- lapply(mu_grid,
                        function(mu) {
                          set.seed(123)
                          #generate response based on mu
                          ctree_sim_resp2 <- list()
                          for(i in 1:n_sim){
                            set.seed(i)
                            ctree_sim_resp2[[i]] <- c(rnorm(n_obs/2, mean = 0,
                                                            sd = 1),
                                                      rnorm(n_obs/2, mean = mu,
                                                            sd = 1))
                          }
                          #first-variable selection
                          thismu_varsel <- pbmapply(resp = ctree_sim_resp2,
                                                    covs = ctree_sim_cov,
                                                    function(resp, covs){
                                                      #first-variable selection
                                                      etree_sel(response = resp,
                                                                covariates = covs,
                                                                case.weights = NULL,
                                                                minbucket = 10,
                                                                alpha = 0.05,
                                                                R = 1000,
                                                                split.type = 'coeff',
                                                                coef.split.type = 'test',
                                                                p.adjust.method = 'fdr')
                                                    },
                                                    SIMPLIFY = TRUE)
                          return(thismu_varsel)

                        })


### Power ###

# Power calculation
ctree_power <- sapply(ctree_mu_grid,
                      function(thismu_varsel) {
                        #proportion of not-null values
                        pow = 1 - sum(sapply(thismu_varsel, is.null)) / n_sim
                        return(pow)
                      })
names(ctree_power) <- as.character(mu_grid)

# Plot
plot(y = ctree_power, x = names(ctree_power), ylab = 'Simulated Power',
     xlab = expression(mu), main = 'Energy Trees', type = 'l')
abline(h = 0.05, lty = 2)


### Conditional probability ###

# Index of x3, that is the discriminating covariate
dis_var_idx <- as.integer(which(names(ctree_sim_cov[[1]]) == 'x6'))

# Conditional probability calculation
ctree_cp <- sapply(ctree_mu_grid,
                   function(thismu_varsel) {
                     #variables' id (only when a variable is selected)
                     var_ifsel = unlist(thismu_varsel[!sapply(thismu_varsel,
                                                              is.null)])
                     #number of times a variable is selected
                     n_var_ifsel = length(var_ifsel)
                     #proportion of times the correct variable is selected
                     cp = sum(var_ifsel == dis_var_idx) / n_var_ifsel
                     return(cp)
                   })
names(ctree_cp) <- as.character(mu_grid)

# Plot
plot(y = ctree_cp, x = names(ctree_cp), ylab = 'Conditional Probability of Correct Split', xlab = expression(mu), main = 'Energy Trees', type = 'l')


### Side-by-side plot ###
old_par = par()
par(mfrow = c(1, 2), pty = 's')
#first plot: power
plot(y = ctree_power, x = names(ctree_power), ylab = 'Simulated Power',
     xlab = expression(mu), type = 'l')
abline(h = 0.05, lty = 2)
#second plot: conditional probability
plot(y = ctree_cp, x = names(ctree_cp), ylab = 'Conditional Probability of Correct Split', xlab = expression(mu), type = 'l')
par(old_par)

