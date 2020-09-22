
# Load packages and functions -------------------------------------------------

# Packages
library(fda.usc)
library(energy)
library(entropy)
library(partykit)
library(cluster)
library(igraph)
library(NetworkDistance)
library(checkmate)
library(future.apply)
future::plan(multisession)

# Functions
source("functions_v2.R")
source("node_v2.R")
source("split_v2.R")
source("party_v2.R")
source("plot_v2.R")


# First-selected variable id --------------------------------------------------

# Function to return the id of the first-selected variable in an etree fit
etree_sel <- function(response,
                      covariates,
                      case.weights = NULL,
                      minbucket = 5,
                      alpha = 0.05,
                      R = 1000,
                      split.type = 'coeff',
                      coef.split.type = 'test',
                      p.adjust.method = 'fdr',
                      nb = 5) {

  # Check whether covariates is a list
  if(!is.list(covariates)) stop("Argument 'covariates' must be provided as a list")

  # Number of covariates
  n.var = length(covariates)

  # If the case weights are not provided, they are all initialized as 1
  if(is.null(case.weights))
    case.weights <- rep(1L, as.numeric(length(response)))

  # New list of covariates (needed here to build the df used by party)
  newcovariates = lapply(covariates, function(j){

    if(class(j) == 'fdata'){

      if(split.type == "coeff"){

        foo <- fda.usc::optim.basis(j, numbasis = nb)
        fd3 <- fda.usc::fdata2fd(foo$fdata.est,
                                 type.basis = "bspline",
                                 nbasis = foo$numbasis.opt)
        foo <- t(fd3$coefs)

      } else if(split.type == "cluster"){

        foo <- as.factor(1:length(response))

      }

      attr(foo, 'cov.type') <- 'fdata'
      return(foo)

    } else if(class(j) == 'list' & all(sapply(j, function(x) attributes(x)$names) == 'diagram')){

      foo <- as.factor(1:length(response))
      attr(foo, 'cov.type') <- 'diagram'
      return(foo)

    } else if(class(j) == 'list' &
              all(sapply(j, class) == 'igraph')){

      if(split.type == "coeff"){

        foo <- graph.shell(j)

      } else if(split.type == "cluster"){

        foo <- as.factor(1:length(response))

      }

      attr(foo, 'cov.type') <- 'graph'
      return(foo)

    } else {

      return(j)

    }
  })

  # Covariates name
  if(!is.null(names(covariates))){
    names(newcovariates) <- names(covariates)
    #control if any name is void, i.e. if it is ''
    no_name <- which(sapply(names(covariates), function(n) n == '', USE.NAMES = FALSE))
    names(newcovariates) <- replace(names(newcovariates),
                                    no_name,
                                    as.factor(1:length(no_name)))
  } else {
    warning('No names available for covariates. Numbers are used instead.')
    names(newcovariates) <- 1:length(newcovariates)
  }

  # Distances
  cov.distance <- lapply(covariates, compute.dissimilarity)

  # Large list with covariates, newcovariates and distances
  covariates.large = list('cov' = covariates, 'newcov' = newcovariates, 'dist' = cov.distance)

  # Number of original covariates
  n.cov = length(covariates.large$cov)

  # Independence test between the response and each covariate
  p = lapply(covariates.large$dist,
             function(cov.dist) {
               ct <- energy::dcor.test(cov.dist,
                                       compute.dissimilarity(response),
                                       R = R)
               if (!is.na(ct$statistic)) {
                 return(c(ct$statistic, ct$p.value))
               } else{
                 c(NA, NA)
               }
             }
  )

  # Create matrix with test stats and p-values; check if p-values are all NULL
  p = t(matrix(unlist(p), ncol = 2, byrow = T))
  rownames(p) <- c("statistic", "p-value")
  if (all(is.na(p[2,]))) return(NULL)

  # Multiple testing correction
  adj_p <- p.adjust(p[2,], method = p.adjust.method)

  # Stop criterion
  if (min(adj_p) > alpha) return(NULL)

  # Variable selection (based on original p-values)
  if (length(which(p[2,] == min(p[2,], na.rm = T))) > 1) {
    xselect <- which.max(p[1,])
    #in case of multiple minima, take that with the highest test statistic
  } else{
    xselect <- which.min(p[2,])
  }

  return(xselect)
}


# Independence analysis -------------------------------------------------------
#or the first-selected variable under independence

# Function to compute the first-selected variable
indep_sim <- function(covariates,
                      response,
                      split.type = c('coeff', 'cluster')){

  # First-selected variable under independence
  indep_first_var <- future_lapply(split.type,
                                   function(type){

                                     # etree fit for each couple of response and cov
                                     mapply(resp = response,
                                            covs = covariates,
                                            function(resp, covs){
                                              etree_sel(response = resp,
                                                        covariates = covs,
                                                        case.weights = NULL,
                                                        minbucket = 10,
                                                        alpha = 1,
                                                        R = 1000,
                                                        split.type = type,
                                                        coef.split.type = 'test',
                                                        p.adjust.method = 'fdr')
                                            },
                                            SIMPLIFY = TRUE)

                                   },
                                   future.seed = TRUE)

  # Rename results based on the split type
  names(indep_first_var) <- split.type

  # Return results
  return(indep_first_var)

}


# Power and conditional probability (growing association) ---------------------

# Function to compute power and cp for each mu in mu_grid
powercp_sim <- function(covariates, #response is built inside (based on mu)
                        ass_cov_idx, #index of the covariate associated with the response
                        mu_grid = seq(0, 1, 0.05),
                        split.type = c('coeff', 'cluster'),
                        minbucket = 10,
                        alpha = 0.05){

  powercp_mu <- future_lapply(mu_grid,
                              function(mu) {

                                # Generate response based on mu
                                resp_thismu <- list()
                                for(i in 1:n_sim){
                                  resp_thismu[[i]] <- c(rnorm(n_obs/2, mean = 0,
                                                              sd = 1),
                                                        rnorm(n_obs/2, mean = mu,
                                                              sd = 1))
                                }

                                # First-selected variable, power and cp for this mu
                                # (and both split.types)
                                powercp_thismu <- lapply(split.type, function(type){

                                  # First-sel variable, power and cp for this mu and type
                                  pcp_mu_type <- powercp_mu_type(covariates = covariates,
                                                                 response = resp_thismu,
                                                                 ass_cov_idx = ass_cov_idx,
                                                                 spl.type = type,
                                                                 minbucket = minbucket,
                                                                 alpha = alpha)

                                  # Return results
                                  return(pcp_mu_type)

                                })

                                # Rename results based on the split type
                                names(powercp_thismu) <- split.type

                                # Return results
                                return(powercp_thismu)

                              },
                              future.seed = TRUE)

  # Rename based on mu grid values
  names(powercp_mu) <- as.character(mu_grid)

  # Return results
  return(powercp_mu)

}

powercp_mu_type <- function(covariates,
                            response,
                            ass_cov_idx,
                            spl.type,
                            minbucket = 10,
                            alpha = 0.05){

  ## First-selected variable ##
  firstvar_mu_type <- mapply(resp = response,
                             covs = covariates,
                             function(resp, covs){
                               set.seed(123)
                               etree_sel(response = resp,
                                         covariates = covs,
                                         case.weights = NULL,
                                         minbucket = minbucket,
                                         alpha = alpha,
                                         R = 1000,
                                         split.type = spl.type,
                                         coef.split.type = 'test',
                                         p.adjust.method = 'fdr')
                             },
                             SIMPLIFY = TRUE)

  ## Power ##
  #Proportion of not-null values
  power_mu_type = 1 - sum(sapply(firstvar_mu_type, is.null)) / n_sim

  ## Conditional probability ##
  #Variables' id (only when a variable is selected)
  var_ifsel = unlist(firstvar_mu_type[!sapply(firstvar_mu_type, is.null)])
  #Number of times a variable is selected
  n_var_ifsel = length(var_ifsel)
  #Proportion of times the correct variable is selected
  cp_mu_type = sum(var_ifsel == ass_cov_idx) / n_var_ifsel

  # Return first variable's ids, power and cp for this mu
  return(list('first_var' = firstvar_mu_type,
              'power' = power_mu_type,
              'cp' = cp_mu_type))

}

