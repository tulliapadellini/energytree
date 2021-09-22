eforest <- function(response,
                    covariates,
                    weights = NULL,
                    ntrees = 100,
                    ncores = 1L,
                    minbucket = 1,
                    alpha = 1,
                    R = 500,
                    split_type = 'cluster',
                    coeff_split_type = 'test',
                    p_adjust_method = 'fdr',
                    perf_metric = NULL,
                    random_covs = NULL) {

  # If perf_metric is NULL, set it to default choices
  if (is.null(perf_metric)) {
    if (is.factor(response)) perf_metric <- 'BAcc' else
      if (is.numeric(response)) perf_metric <- 'RMSPE'
  } else {
    # If perf_metric is given, check that it has an acceptable value
    if (is.factor(response)) {
      if (isFALSE(perf_metric %in% c('BAcc', 'WBAcc')))
        stop("The value provided for 'perf_metric' is not available for this
             task")
    } else if (is.numeric(response)) {
      if (isFALSE(perf_metric %in% c('MAPE', 'RMSPE', 'NRMSE', 'MAE',
                                     'MedianAE', 'MedianAPE', 'MSE', 'RAE',
                                     'RMSE', 'RMLSE')))
        stop("The value provided for 'perf_metric' is not available for this
             task")
    }
  }

  # Number of observations
  nobs <- length(response)

  # New list of covariates
  newcovariates <- .create_newcov(covariates = covariates,
                                  response = response,
                                  split_type = split_type)

  # Distances
  cov_distance <- lapply(covariates, dist_comp)

  # Large list with covariates, newcovariates and distances
  covariates_large <- list('cov' = covariates,
                           'newcov' = newcovariates,
                           'dist' = cov_distance)

  # Large list with response and the corresponding distances
  response_large <- list('response' = response,
                         'response_dist' = dist_comp(response))

  # Generate B bootstrap samples
  set.seed(12345)
  boot_idx <- lapply(1:ntrees,
                     function(b) sample.int(nobs, replace = TRUE))

  # Energy Tree fits for each bootstrap sample
  etree_boot_fits <- parallel::mclapply(boot_idx, function(b_i) {

    # Covariates and response for this bootstrap sample
    boot_cov_large <- lapply(covariates_large[1:2],
                             function(cl) lapply(cl, function(cov) {
                               if (class(cov) == 'data.frame') cov[b_i, ]
                               else cov[b_i]
                             }
                             ))
    # Re-index newcov only if using 'cluster'
    if (split_type == 'cluster') {
      boot_cov_large$newcov[[1]] <- factor(1:nobs)
      boot_cov_large$newcov[[2]] <- factor(1:nobs)
    }
    boot_cov_large$dist <- lapply(covariates_large[[3]],
                                  function(cov_dist) {
                                    boot_dist <- usedist::dist_subset(cov_dist,
                                                                      b_i)
                                    boot_dist <-
                                      usedist::dist_setNames(boot_dist, 1:nobs)
                                    return(as.matrix(boot_dist))
                                  })
    resp_dist <- usedist::dist_subset(response_large$response_dist, b_i)
    resp_dist <- as.matrix(usedist::dist_setNames(resp_dist, 1:nobs))
    boot_resp_large <- list('response' = response[b_i],
                            'response_dist' = resp_dist)

    # Energy Trees fit
    e_fit <- etree(response = boot_resp_large,
                   covariates = boot_cov_large,
                   weights = weights,
                   minbucket = minbucket,
                   alpha = alpha,
                   R = R,
                   split_type = split_type,
                   coeff_split_type = coeff_split_type,
                   p_adjust_method = p_adjust_method,
                   random_covs = random_covs)

    # Remove .Environment attribute from terms
    attr(e_fit$terms, '.Environment') <- NULL
    print('one round again')
    # Return
    return(e_fit)

  },
  mc.cores = ncores,
  mc.set.seed = TRUE)

  # List containing for each obs the predicted response with all the corresponding OOB trees
  oob_pred_resp <- lapply(1:nobs, function(i) {

    # All covariates for observation i (to be used as newdata in predict())
    covs_i <- lapply(covariates, function(cov) cov[i])

    # Is the observation Out Of Bag in each bootstrap sample?
    is.oob <- sapply(boot_idx, function(b) !(i %in% b))

    # Indices of bootstrap samples for which the observation is OOB
    is.oob_idx <- which(is.oob)

    # Predict response for this obs only with trees whose index is is.oob_idx
    oob_pred_resp <- sapply(is.oob_idx, function(o) {
      predict(etree_boot_fits[[o]],
              newdata = covs_i)
    })

    # Return predicted response
    return(oob_pred_resp)

  })


  # Predicted responses and OOB performance metric calculation
  if (is.factor(response)) {

    ## Classification ##

    # Predicted response: majority voting rule
    pred_resp <- factor(sapply(oob_pred_resp,
                               function(i) names(which.max(table(i)))
    ))

    # OOB performance metric (measured via BAcc or WBAcc)
    if (perf_metric == 'BAcc' || perf_metric == 'WBAcc') {

      # Balanced Accuracy for each class (each given by (TP/P + TN/N) / 2)
      bal_accs <- sapply(levels(response),
                         function(lev) {
                           # Sensitivity
                           true_pos <- sum(pred_resp == lev &
                                             response == lev)
                           pos <- sum(response == lev)
                           sens <- true_pos / pos
                           # Specificity
                           true_neg <- sum(pred_resp != lev &
                                             response != lev)
                           neg <- sum(response != lev)
                           spec <- true_neg / neg
                           # Balanced Accuracy (with the current class as pos)
                           bal_acc_lev <- (sens + spec) / 2
                           # Return Balanced Accuracy
                           return(bal_acc_lev)
                         })

      # OOB performance metric
      #Balanced Accuracy -> arithmetic mean
      #Weighted Balanced Accuracy -> weighted mean (with class sizes as weights)
      oob_perf_metric <- switch(perf_metric,
                                BAcc = mean(bal_accs),
                                WBAcc = sum(bal_accs * table(response)) /
                                  length(response))

    }

  } else if (is.numeric(response)) {

    ## Regression ##

    # Predicted response: average
    pred_resp <- sapply(oob_pred_resp, mean)

    # OOB performance metric: various choices (default is 'RMSPE')
    oob_perf_metric <-
      switch(perf_metric,
             MAPE = MLmetrics::MAPE(pred_resp, response),
             RMSPE = MLmetrics::RMSPE(pred_resp, response),
             NRMSE = MLmetrics::RMSE(pred_resp, response) / mean(response),
             MAE = MLmetrics::MAE(pred_resp, response),
             MedianAE = MLmetrics::MedianAE(pred_resp, response),
             MedianAPE = MLmetrics::MedianAPE(pred_resp, response),
             MSE = MLmetrics::MSE(pred_resp, response),
             RAE = MLmetrics::RAE(pred_resp, response),
             RMSE = MLmetrics::RMSE(pred_resp, response),
             RMLSE = MLmetrics::RMSLE(pred_resp, response))

  }

  # Return couple (etree_boot_fits, oob_perf_metric)
  return(list(ensemble =  etree_boot_fits, oob_perf_metric = oob_perf_metric))

}

predict_eforest <- function(eforest_obj, newdata = NULL) {

  # Individual predictions with base learners
  #(newdata check, split retrieval, newcov computation are all done in predict)
  ind_pred_resp <- sapply(eforest_obj$ensemble, function(tree) {
    predict(tree, newdata = newdata)
  })

  # Predict response, differenty based on the type of problem (CLS or REG)
  response <- eforest_obj$ensemble[[1]]$fitted$`(response)`
  if (is.factor(response)) {

    # Majority voting rule
    pred_resp <- apply(ind_pred_resp, 1, function(i) names(which.max(table(i))))
    pred_resp <- factor(pred_resp)

  } else if (is.numeric(response)) {

    # Average
    pred_resp <- apply(ind_pred_resp, 1, mean)

  }

  # Return predicted response
  return(pred_resp)

}

