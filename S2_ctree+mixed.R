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
n_sim <- 1000

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
  #x6 <- list(lapply(1:n_obs, function(j) igraph::sample_smallworld(1, 100, 5, 0.1)))
  # Graph: Preferential Attachment model (Barabasi-Albert, 1999)
  #x7 <- list(lapply(1:n_obs, function(j) igraph::sample_pa(100)))
  # Functional: Ornstein-Uhlenbeck random process
  x8 <- fda.usc::rproc2fdata(n_obs, seq(0, 1, len = 20), sigma = 'OrnsteinUhlenbeck')
  # Functional: Wiener random process
  #x9 <- fda.usc::rproc2fdata(n_obs, seq(0, 1, len = 20), sigma = 'brownian')
  # Covariates list for this simulation
  etree_sim_cov[[i]] <- list(x1 = x1, x2 = x2,
                         x3 = x3, x4 = x4,
                         x5 = x5, #x6 = x6, x7 = x7,
                         x8 = x8)#, x9 = x9)
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

# First-variable selection
etree_first_var <- pbmapply(resp = etree_sim_resp,
                            covs = etree_sim_cov,
                            function(resp, covs){
                              set.seed(123)
                              etree_sel(response = resp,
                                        covariates = covs,
                                        case.weights = NULL,
                                        minbucket = 10,
                                        alpha = 0.05,
                                        R = 1000,
                                        split.type = 'cluster',
                                        coef.split.type = 'test',
                                        p.adjust.method = 'fdr')
                            },
                            SIMPLIFY = TRUE)

# Frequency plot
barplot(table(etree_first_var))

# Proportions
prop.table(table(etree_first_var))


# Power and conditional probability (growing association) ---------------------

# Covariates: keep the first three for numeric and only one for each structured
etree_sim_cov2 <- lapply(etree_sim_cov,
                         function(i_sim) i_sim[c('x1', 'x2', 'x3', 'x5', 'x8')])

# Grid values for mu
mu_grid = seq(0, 1, 0.05)

# etree fit: n_sim simulations for each value of mu
etree_mu_grid <- lapply(mu_grid,
                        function(mu) {
                          set.seed(123)
                          #generate response based on mu
                          etree_sim_resp2 <- list()
                          for(i in 1:n_sim){
                            set.seed(i)
                            etree_sim_resp2[[i]] <- c(rnorm(n_obs/2, mean = 0,
                                                            sd = 1),
                                                      rnorm(n_obs/2, mean = mu,
                                                            sd = 1))
                          }
                          #first-variable selection
                          thismu_varsel <- pbmapply(resp = etree_sim_resp2,
                                                    covs = etree_sim_cov2,
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
etree_power <- sapply(etree_mu_grid,
                      function(thismu_varsel) {
                        #proportion of not-null values
                        pow = 1 - sum(sapply(thismu_varsel, is.null)) / n_sim
                        return(pow)
                      })
names(etree_power) <- as.character(mu_grid)

# Plot
plot(y = etree_power, x = names(etree_power), ylab = 'Simulated Power',
     xlab = expression(mu), main = 'Energy Trees', type = 'l')
abline(h = 0.05, lty = 2)


### Conditional probability ###

# Index of x3, that is the discriminating covariate
dis_var_idx <- as.integer(which(names(etree_sim_cov2[[1]]) == 'x3'))

# Conditional probability calculation
etree_cp <- sapply(etree_mu_grid,
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
names(etree_cp) <- as.character(mu_grid)

# Plot
plot(y = etree_cp, x = names(etree_cp), ylab = 'Conditional Probability of Correct Split', xlab = expression(mu), main = 'Energy Trees', type = 'l')


### Side-by-side plot ###
old_par = par()
par(mfrow = c(1, 2), pty = 's')
#first plot: power
plot(y = etree_power, x = names(etree_power), ylab = 'Simulated Power',
     xlab = expression(mu), type = 'l')
abline(h = 0.05, lty = 2)
#second plot: conditional probability
plot(y = etree_cp, x = names(etree_cp), ylab = 'Conditional Probability of Correct Split', xlab = expression(mu), type = 'l')
par(old_par)
