# Number of observations
n_obs <- 100

# Number of simulations
n_sim <- 1000



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
  persistence_cov <- pbmcapply::pbmclapply(1:n_sim, diag_gen, mc.cores = 30)
  saveRDS(persistence_cov, file = "sim/persistence_x1.rds")
}


if(FALSE){
  diag_gen2 = function(i){
    set.seed(i)
    x_data1 = lapply(rep(100, n_obs), function(x) TDA::circleUnif(x))

    x1 = lapply(x_data1, TDA::ripsDiag, maxdimension = 1, maxscale = 3) # 6 seconds
    #only covariate to have obs divided into two classes
    return(list(x1=x1))
  }

  # Covariates list
  persistence_cov2 <- pbmcapply::pbmclapply(1:n_sim, diag_gen2, mc.cores = 30)
  saveRDS(persistence_cov2, file = "sim/persistence_x1_wodistinction.rds")
}




if(FALSE){
  diag_gen3 = function(i){
    set.seed(i)
    x_data2 = lapply(rep(100, n_obs), function(x) tdaunif::sample_swiss_roll(x))
    x_data3 = lapply(rep(100, n_obs), function(x) tdaunif::sample_trefoil(x))
    x_data4 = lapply(rep(100, n_obs), function(x) tdaunif::sample_torus(x))
    x_data5 = lapply(rep(100, n_obs), function(x) tdaunif::sample_klein(x))


    x2 = lapply(x_data2, TDA::ripsDiag, maxdimension = 1, maxscale = 3)
    x3 = lapply(x_data3, TDA::ripsDiag, maxdimension = 1, maxscale = 3)
    x4 = lapply(x_data4, TDA::ripsDiag, maxdimension = 1, maxscale = 3)
    x5 = lapply(x_data5, TDA::ripsDiag, maxdimension = 1, maxscale = 3)

    x1 <- readRDS("sim/persistence_x1_small.rds")
    x1 <- x1[[i]]

    return(list(x1 = x1, x2 = x2, x3 = x3, x4 = x4))
  }

  # Covariates list
  persistence_cov3 <- pbmcapply::pbmclapply(1:n_sim, diag_gen3, mc.cores = 20)
  saveRDS(persistence_cov3, file = "sim/persistence_allcov.rds")
}

