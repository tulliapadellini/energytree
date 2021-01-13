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
source('etree_sim_fun.R')


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

# Save
save(ctree_first_var, file = 'sim_results/ctree_indep_sim.RData')

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

# Save
save(ctree_mu_grid, file = 'sim_results/ctree_powercp_sim.RData')

### Power ###

# Power calculation
ctree_power <- sapply(ctree_mu_grid,
                      function(thismu_varsel) {
                        #proportion of not-null values
                        pow = 1 - sum(sapply(thismu_varsel, is.null)) / n_sim
                        return(pow)
                      })
names(ctree_power) <- as.character(mu_grid)

# Power confidence intervals
c_int_power <- sapply(ctree_power, FUN = function(power){
  round(power + c(-1, 1) * z * sqrt(power * (1 - power) / n), 3)
})


### Conditional probability ###

# Index of x6, that is the discriminating covariate
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

# CP confidence intervals
c_int_cp <- sapply(ctree_cp, FUN = function(cp){
  round(cp + c(-1, 1) * z * sqrt(cp * (1 - cp) / n), 3)
})


### Plot ###

# Setup pdf
pdf('sim_results_plots/traditional.pdf', width = 8.4, height = 4)
par(mfrow = c(1, 2), mar = c(4, 4, 0.5, 2) + 0.1,  pty = 's')

# First plot: Power (with confidence intervals and 0.05 baseline)
plot(y = ctree_power, x = names(ctree_power), ylab = 'Power',
     xlab = expression(mu), type = 'l', lty = 2, ylim = c(0, 1),
     col = rgb(56, 180, 24, maxColorValue = 255))
polygon(c(mu_grid, rev(mu_grid)),
        c(c_int_power[1,], rev(c_int_power[2,])),
        col = rgb(0, 255, 6, alpha = 30, maxColorValue = 255),
        border = NA)
abline(h = 0.05, lty = 3, col = 'Red')

# Second plot: Conditional Probability (with conf. intervals and 0.2 baseline)
plot(y = ctree_cp, x = names(ctree_cp), ylab = 'Conditional Probability of Correct Split', xlab = expression(mu), type = 'l', lty = 2, ylim = c(0, 1),
     col = rgb(56, 180, 24, maxColorValue = 255))
polygon(c(mu_grid, rev(mu_grid)),
        c(c_int_cp[1,], rev(c_int_cp[2,])),
        col = rgb(0, 255, 6, alpha = 30, maxColorValue = 255),
        border = NA)
abline(h = 0.2, lty = 3, col = 'Red')

# Close pdf
dev.off()


