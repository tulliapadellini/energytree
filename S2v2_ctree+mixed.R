### S2 - Extension of simulations from Hothorn et al. (2006) ###
###      to the mixed and structured case ###

# Load functions
source('etree_sim_future.R')
future::plan(multisession)


# Initialization --------------------------------------------------------------

# Initialize number of observations and number of simulations.
# The response variable is the same for all the independence cases.
# All the cases in this file are by the default run with both 'coeff' and
# 'cluster' methods.

# Number of observations
n_obs <- 100

# Number of simulations
n_sim <- 10000

# Response variable for all the independence cases
indep_resp <- list()
set.seed(123)
for(i in 1:n_sim){
  indep_resp[[i]] <- c(rnorm(n_obs, mean = 0, sd = 1))
}

# Grid values for mu (needed for power+cp analysis)
mu_grid = seq(0, 1, 0.05)


# Section zero: one structured covariate at a time ----------------------------

# Having only one covariate, it only makes sense to perform a power analysis.

##### S2.01: Graphs #####

# Dataset simulation
onegraph <- list()
for(i in 1:n_sim){
  set.seed(i)
  # Erdos-Renyi (1959) model with different connection prob for the two classes
  x1 <- c(lapply(1:(n_obs/2), function(j) igraph::sample_gnp(100, 0.05)),
          #0.05 is the approx value for the threshold to have a connected graph
          lapply(1:(n_obs/2), function(j) igraph::sample_gnp(100, 0.95)))
  # Covariates list
  onegraph[[i]] <- list(x1 = x1)
}

# Power
onegraph_power_sim <- powercp_sim(covariates = onegraph,
                                  ass_cov_idx = 1,
                                  mu_grid = mu_grid,
                                  split.type = c('coeff', 'cluster'),
                                  minbucket = 10,
                                  alpha = 0.05)

# Save
save(onegraph_power_sim, file = 'sim_results/onegraph_power_sim.RData')


##### S2.02: Functional #####

# Dataset
onefunctional <- list()
for(i in 1:n_sim){
  set.seed(i)
  # Gaussian random process
  x1 <- c(fda.usc::rproc2fdata(n_obs/2, seq(0, 1, len = 100), sigma = 1),
          #second half have a mean equal to 3 so that they are well differentiated
          fda.usc::rproc2fdata(n_obs/2, seq(0, 1, len = 100), mu = rep(3, 100),
                               sigma = 1))
  # Covariates list
  onefunctional[[i]] <- list(x1 = x1)
}

# Power
onefunctional_power_sim <- powercp_sim(covariates = onefunctional,
                                       ass_cov_idx = 1,
                                       mu_grid = mu_grid,
                                       split.type = c('coeff', 'cluster'),
                                       minbucket = 10,
                                       alpha = 0.05)

# Save
save(onefunctional_power_sim, file = 'sim_results/onefunctional_power_sim.RData')



##### S2.03: Persistence #####

# Dataset
onepersistence <- list()
if(FALSE){
  diag_gen = function(i){
    set.seed(i)
    x_data1 = lapply(rep(100, n_obs/2), function(x) TDA::circleUnif(x))
    x_data2 = lapply(rep(50, n_obs/2), function(x) rbind(TDA::circleUnif(x), TDA::circleUnif(x)+2))

    x_data = c(x_data1, x_data2)

    x1 = lapply(x_data, TDA::ripsDiag, maxdimension = 1, maxscale = 3) # 6 seconds
    #only covariate to have obs divided into two classes
    return(list(x1=x1))
  }

# Covariates list
  persistence_cov <- pbmcapply::pbmclapply(1:n_sim, diag_gen, mc.cores = 20)
  saveRDS(persistence_cov, file = "sim/persistence_x1.rds")
}

persistence_cov = readRDS("sim/persistence_x1.rds")

# Power
onepersistence_power_sim <- powercp_sim(covariates = onepersistence,
                                        ass_cov_idx = 1,
                                        mu_grid = mu_grid,
                                        split.type = c('cluster'),
                                        minbucket = 10,
                                        alpha = 0.05)

# Save
save(onepersistence_power_sim, file = 'sim_results/onepersistence_power_sim.RData')



# First section: one-type structured covariates -------------------------------

# Here the focus is on a single type of structured covariate at a time.
# For each type, both independence and power+cp analyses are carried out.

##### S2.A: Graphs #####

# Dataset simulation
graph_cov <- list()
for(i in 1:n_sim){
  set.seed(i)
  # Erdos-Renyi (1959) model with different connection prob for the two classes
  x1 <- c(lapply(1:(n_obs/2), function(j) igraph::sample_gnp(100, 0.05)),
          #0.05 is the approx value for the threshold to have a connected graph
          lapply(1:(n_obs/2), function(j) igraph::sample_gnp(100, 0.95)))
  # Small-World model (Watts-Strogatz, 1998) with same rewiring prob
  x2 <- lapply(1:n_obs, function(j) igraph::sample_smallworld(1, 100, 5, 0.1))
  # Preferential Attachment model (Barabasi-Albert, 1999)
  x3 <- lapply(1:n_obs, function(j) igraph::sample_pa(100))
  # Forest Fire Network model (Leskovec, Kleinberg, Faloutsos, 2005)
  x4 <- lapply(1:n_obs, function(j) igraph::sample_forestfire(100, 0.2,
                                                              directed = FALSE))
  # Geometric random graphs on a torus
  x5 <- lapply(1:n_obs, function(j) igraph::sample_grg(100, 0.2, torus = TRUE))
  # Covariates list
  graph_cov[[i]] <- list(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5)
}

# Independence
graph_indep_sim <- indep_sim(covariates = graph_cov,
                             response = indep_resp,
                             split.type = c('coeff', 'cluster'))

# Power and conditional probability
graph_powercp_sim <- powercp_sim(covariates = graph_cov,
                                  ass_cov_idx = 1,
                                  mu_grid = mu_grid,
                                  split.type = c('coeff', 'cluster'),
                                  minbucket = 10,
                                  alpha = 0.05)

# Save
save(graph_indep_sim, file = 'sim_results/graph_indep_sim.RData')
save(graph_powercp_sim, file = 'sim_results/graph_powercp_sim.RData')


##### S2.B: Functional #####

# Dataset
functional_cov <- list()
for(i in 1:n_sim){
  set.seed(i)
  # Gaussian random process
  x1 <- c(fda.usc::rproc2fdata(n_obs/2, seq(0, 1, len = 100), sigma = 1),
          #second half have a mean equal to 3 so that they are well differentiated
          fda.usc::rproc2fdata(n_obs/2, seq(0, 1, len = 100), mu = rep(3, 100),
                               sigma = 1))
  # Wiener random process
  x2 <- fda.usc::rproc2fdata(n_obs, seq(0, 1, len = 100), sigma = 'brownian')
  # Ornstein-Uhlenbeck random process
  x3 <- fda.usc::rproc2fdata(n_obs, seq(0, 1, len = 100), sigma = 'OrnsteinUhlenbeck',
                             par.list = list(theta = 1)) #theta to diff from x2
  # Fractional brownian random process (H = 0.8)
  x4 <- fda.usc::rproc2fdata(n_obs, seq(0, 1, len = 100), sigma = 'fbrownian',
                             par.list = list(H = 0.8)) #H to diff from x2
  # Gaussian random process with exponential variogram
  x5 <- fda.usc::rproc2fdata(n_obs, seq(0, 1, len = 100), sigma = 'vexponential')
  # Covariates list
  functional_cov[[i]] <- list(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5)
}

# Independence
functional_indep_sim <- indep_sim(covariates = functional_cov,
                                  response = indep_resp,
                                  split.type = c('coeff', 'cluster'))

# Power and conditional probability
functional_powercp_sim <- powercp_sim(covariates = functional_cov,
                                      ass_cov_idx = 1,
                                      mu_grid = mu_grid,
                                      split.type = c('coeff', 'cluster'),
                                      minbucket = 10,
                                      alpha = 0.05)

# Save
save(functional_indep_sim, file = 'sim_results/functional_indep_sim.RData')
save(functional_powercp_sim, file = 'sim_results/functional_powercp_sim.RData')


##### S2.C: Persistence #####

# Dataset
persistence_cov <- readRDS(file = "sim/persistence_allcov.rds")


# Independence
persistence_indep_sim <- indep_sim(covariates = persistence_cov,
                                   response = indep_resp,
                                   split.type = c('coeff', 'cluster'))

# Power and conditional probability
persistence_powercp_sim <- powercp_sim(covariates = persistence_cov,
                                       ass_cov_idx = 1,
                                       mu_grid = mu_grid,
                                       split.type = c('coeff', 'cluster'),
                                       minbucket = 10,
                                       alpha = 0.05)

# Save
save(persistence_indep_sim, file = 'sim_results/persistence_indep_sim.RData')
save(persistence_powercp_sim, file = 'sim_results/persistence_powercp_sim.RData')



# Second section: mixed-type covariates ---------------------------------------

# Here the focus is on mixed-type covariates, one for each type.
# Independence analysis is carried out only once.
# Power+cp analysis is carried out once for each covariate type, and that covariate
# is the one increasingly associated with the response (by varying the response).
# Independence dataset have obs divided into two classes for each covariate;
# power+cp datasets have obs divided into two classes only for the associated cov.

##### S2.D: Independence analysis with mixed covariates #####

# Dataset
#division into two classes is not really necessary, it's just for parallelism
mixed_cov <- list()
for(i in 1:n_sim){
  set.seed(i)
  # Numeric: uniform between 0 and 1
  x1 <- c(runif(n_obs/2, min = 0, max = 0.5),
          runif(n_obs/2, min = 0.5, max = 1))
  # Nominal: 0s and 1s
  x2 <- factor(c(rep(0, n_obs/2),
                 rep(1, n_obs/2)))
  # Graph: Erdos-Renyi (1959) model
  x3 <- c(lapply(1:(n_obs/2), function(j) igraph::sample_gnp(100, 0.05)),
          lapply(1:(n_obs/2), function(j) igraph::sample_gnp(100, 0.95)))
  # Functional: Gaussian process
  x4 <- c(fda.usc::rproc2fdata(n_obs/2, seq(0, 1, len = 100), sigma = 1),
          fda.usc::rproc2fdata(n_obs/2, seq(0, 1, len = 100), mu = rep(3, 100),
                               sigma = 1))
  # Persistence: #covariate that had obs divided into two classes in S2.C
  x5 <- readRDS("sim/persistence_x1.rds")
  x5 <- x5[[i]]
  # Covariates list
  mixed_cov[[i]] <- list(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5)
}

# Independence
mixed_indep_sim <- indep_sim(covariates = mixed_cov,
                             response = indep_resp,
                             split.type = c('coeff', 'cluster'))

# Save
save(mixed_indep_sim, file = 'sim_results/mixed_indep_sim.RData')


##### S2.E: Response associated with a numeric covariate #####

# Dataset
assnum_cov <- list()
for(i in 1:n_sim){
  set.seed(i)
  # Numeric: uniform between 0 and 1
  x1 <- c(runif(n_obs/2, min = 0, max = 0.5),
          runif(n_obs/2, min = 0.5, max = 1))
  # Nominal: 0s and 1s
  x2 <- factor(rbinom(n_obs, 1, 0.5))
  # Graph: Erdos-Renyi (1959) model
  x3 <- lapply(1:n_obs, function(j) igraph::sample_gnp(100, 0.05))
  # Functional: Gaussian process
  x4 <- fda.usc::rproc2fdata(n_obs, seq(0, 1, len = 100), sigma = 1)
  # Persistence:
  x5 <- readRDS("sim/persistence_x1_wodistinction.rds")  #covariate that had obs divided into two classes in S2.C, but w/o division
  x5 <- x5[[i]]
  # Covariates list
  assnum_cov[[i]] <- list(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5)
}

# Power and conditional probability
assnum_powercp_sim <- powercp_sim(covariates = assnum_cov,
                                  ass_cov_idx = 1,
                                  mu_grid = mu_grid,
                                  split.type = c('coeff', 'cluster'),
                                  minbucket = 10,
                                  alpha = 0.05)

# Save
save(assnum_powercp_sim, file = 'sim_results/assnum_powercp_sim.RData')


##### S2.F: Response associated with a nominal covariate #####

# Dataset
assnom_cov <- list()
for(i in 1:n_sim){
  set.seed(i)
  # Numeric: uniform between 0 and 1
  x1 <- runif(n_obs, min = 0, max = 1)
  # Nominal: 0s and 1s
  x2 <- factor(c(rep(0, n_obs/2),
                 rep(1, n_obs/2)))
  # Graph: Erdos-Renyi (1959) model
  x3 <- lapply(1:n_obs, function(j) igraph::sample_gnp(100, 0.05))
  # Functional: Gaussian process
  x4 <- fda.usc::rproc2fdata(n_obs, seq(0, 1, len = 100), sigma = 1)
  # Persistence:
  x5 <- readRDS("sim/persistence_x1_wodistinction.rds")  #covariate that had obs divided into two classes in S2.C, but w/o division
  x5 <- x5[[i]]
  assnom_cov[[i]] <- list(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5)
}

# Power and conditional probability
assnom_powercp_sim <- powercp_sim(covariates = assnom_cov,
                                  ass_cov_idx = 2,
                                  mu_grid = mu_grid,
                                  split.type = c('coeff', 'cluster'),
                                  minbucket = 10,
                                  alpha = 0.05)

# Save
save(assnom_powercp_sim, file = 'sim_results/assnom_powercp_sim.RData')


##### S2.G: Response associated with a graph covariate #####

# Dataset
assgph_cov <- list()
for(i in 1:n_sim){
  set.seed(i)
  # Numeric: uniform between 0 and 1
  x1 <- runif(n_obs, min = 0, max = 1)
  # Nominal: 0s and 1s
  x2 <- factor(rbinom(n_obs, 1, 0.5))
  # Graph: Erdos-Renyi (1959) model
  x3 <- c(lapply(1:(n_obs/2), function(j) igraph::sample_gnp(100, 0.05)),
          lapply(1:(n_obs/2), function(j) igraph::sample_gnp(100, 0.95)))
  # Functional: Gaussian process
  x4 <- fda.usc::rproc2fdata(n_obs, seq(0, 1, len = 100), sigma = 1)
  # Persistence:
  x5 <- readRDS("sim/persistence_x1_wodistinction.rds")  #covariate that had obs divided into two classes in S2.C, but w/o division
  x5 <- x5[[i]]
  assgph_cov[[i]] <- list(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5)
}

# Power and conditional probability
assgph_powercp_sim <- powercp_sim(covariates = assgph_cov,
                                  ass_cov_idx = 3,
                                  mu_grid = mu_grid,
                                  split.type = c('coeff', 'cluster'),
                                  minbucket = 10,
                                  alpha = 0.05)

# Save
save(assgph_powercp_sim, file = 'sim_results/assgph_powercp_sim.RData')


##### S2.H: Response associated with a functional covariate #####

# Dataset
assfun_cov <- list()
for(i in 1:n_sim){
  set.seed(i)
  # Numeric: uniform between 0 and 1
  x1 <- runif(n_obs, min = 0, max = 1)
  # Nominal: 0s and 1s
  x2 <- factor(rbinom(n_obs, 1, 0.5))
  # Graph: Erdos-Renyi (1959) model
  x3 <- lapply(1:n_obs, function(j) igraph::sample_gnp(100, 0.05))
  # Functional: Gaussian process
  x4 <- c(fda.usc::rproc2fdata(n_obs/2, seq(0, 1, len = 100), sigma = 1),
          fda.usc::rproc2fdata(n_obs/2, seq(0, 1, len = 100), mu = rep(3, 100),
                               sigma = 1))
  # Persistence:
  x5 <- readRDS("sim/persistence_x1_wodistinction.rds")  #covariate that had obs divided into two classes in S2.C, but w/o division
  x5 <- x5[[i]]
  assfun_cov[[i]] <- list(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5)
}

# Power and conditional probability
assfun_powercp_sim <- powercp_sim(covariates = assfun_cov,
                                  ass_cov_idx = 4,
                                  mu_grid = mu_grid,
                                  split.type = c('coeff', 'cluster'),
                                  minbucket = 10,
                                  alpha = 0.05)

# Save
save(assfun_powercp_sim, file = 'sim_results/assfun_powercp_sim.RData')


##### S2.I: Response associated with a persistence covariate #####

# Dataset
assper_cov <- list()
for(i in 1:n_sim){
  set.seed(i)
  # Numeric: uniform between 0 and 1
  x1 <- runif(n_obs, min = 0, max = 1)
  # Nominal: 0s and 1s
  x2 <- factor(rbinom(n_obs, 1, 0.5))
  # Graph: Erdos-Renyi (1959) model
  x3 <- lapply(1:n_obs, function(j) igraph::sample_gnp(100, 0.05))
  # Functional: Gaussian process
  x4 <- fda.usc::rproc2fdata(n_obs, seq(0, 1, len = 100), sigma = 1)
  # Persistence:
  x5 <- readRDS("sim/persistence_x1.rds")  #covariate that had obs divided into two classes in S2.C, but w/o division
  x5 <- x5[[i]]
  # Covariates list
  assper_cov[[i]] <- list(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5)
}

# Power and conditional probability
assper_powercp_sim <- powercp_sim(covariates = assper_cov,
                                  ass_cov_idx = 5,
                                  mu_grid = mu_grid,
                                  split.type = c('coeff', 'cluster'),
                                  minbucket = 10,
                                  alpha = 0.05)

# Save
save(assper_powercp_sim, file = 'sim_results/assper_powercp_sim.RData')


