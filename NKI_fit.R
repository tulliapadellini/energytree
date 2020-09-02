

# Load packages and functions ---------------------------------------------

library(fda.usc)
library(energy)
library(entropy)
library(partykit)
library(cluster)
library(igraph)
library(NetworkDistance)
library(ggparty)
library(checkmate)
source("functions_v2.R")
source("node_v2.R")
source("split_v2.R")
source("party_v2.R")
source("plot_v2.R")
source("ggparty_v1.R")
source("get_plot_data_v1.R")
source("NKI_data_import.R")


# Import data -------------------------------------------------------------

nki <- generate_dataset(data_folder = 'NKI_Rockland/',
                        y_filename = 'NKI_clinical_information.txt',
                        y_column = 'WASI_FULL_4',
                        output_filename = 'NKIdata.RData',
                        output_folder = ".",
                        ext_save = FALSE)



# Dataset construction ----------------------------------------------------

# Response
resp <- nki$y

# Covariates list
cov.list <- list(lapply(nki$structural, function(g) igraph::graph_from_adjacency_matrix(g, weighted = T)),
                 lapply(nki$functional, function(g) igraph::graph_from_adjacency_matrix(g, weighted = T)))


# Energy Tree fit ---------------------------------------------------------

# Fit
set.seed(2948)
etree_fit <- etree(response = resp,
                   covariates = cov.list,
                   case.weights = NULL,
                   minbucket = 5,
                   alpha = 0.5,
                   R = 1000,
                   split.type = 'cluster',
                   coef.split.type = 'test',
                   p.adjust.method = 'fdr')

# Plot
plot(etree_fit)

# Fitted values
y_fitted <- predict(etree_fit)

# Mean Error Prediction
(MEP_etree <- (sum((resp-y_fitted)^2)/length(resp))/(var(resp)))

# Root Mean Square Error
(MEP_etree <- sqrt(sum((resp-y_fitted)^2)/length(resp)))

# Mean Square Percentage Error
(MEP_etree <- sum(((resp-y_fitted)/resp)^2)/length(resp))


# Prediction --------------------------------------------------------------

# Predicted values
y_pred <- predict(etree_fit, newdata = cov.list)



##### NEW #####

# Mixed covariates -----------------------------------------------------------

# Response
resp <- nki$y

# Graph covariates
nki_str <- lapply(nki$structural,
                  function(g) igraph::graph_from_adjacency_matrix(g, weighted = T))
nki_fun <- lapply(nki$functional,
                  function(g) igraph::graph_from_adjacency_matrix(g, weighted = T))

# Import clinical information
NKI_clinical <- read.csv("NKI_clinical_information.txt",
                         stringsAsFactors = TRUE)

# Select individuals who also have structural and functional matrices
sel_id <- names(nki$structural)
sel_idx <- (NKI_clinical$Subject %in% sel_id)
#remark: already done for str, fun and wasi coming from nki

# Covariates list
cov.list <- list(SCM = nki_str,
                 FCM = nki_fun,
                 Gender = NKI_clinical$Gender[sel_idx],
                 Age = NKI_clinical$Age[sel_idx])

# Fit
set.seed(2948)
etree_fit <- etree(response = resp,
                   covariates = cov.list,
                   case.weights = NULL,
                   minbucket = 5,
                   alpha = 0.5,
                   R = 1000,
                   split.type = 'cluster',
                   coef.split.type = 'test')

# Plot
plot(etree_fit)

# Fitted values
y_fitted <- predict(etree_fit)

# Mean Error Prediction
(MEP_etree <- (sum((resp-y_fitted)^2)/length(resp))/(var(resp)))

# Root Mean Square Error
(RMSE_etree <- sqrt(sum((resp-y_fitted)^2)/length(resp)))

# Mean Square Percentage Error
(MSPE_etree <- sum(((resp-y_fitted)/resp)^2)/length(resp))

# Predicted values
y_pred <- predict(etree_fit, newdata = cov.list)


# ggplot ----------------------------------------------------------------------

## Basic ##
ggparty(etree_fit) +
  geom_edge() +
  geom_edge_label() +
  geom_node_splitvar() +
  geom_node_info()

## More advanced ##

# Function to plot asterisks instead of pvalues
asterisk_sign <- function(p_value) {
  if (p_value < 0.05) return(c("***"))
  if (p_value < 0.1) return(c("**"))
  if (p_value < 0.2) return(c("*"))
  else return("")
}

# Main call
gg <- ggparty(etree_fit,
              terminal_space = 0.25,
              add_vars = list(nodedata_resp =
                                function(data, node){
                                  list(data_party(node)$'(response)')
                                }))
gg +
  # Plot edges
  geom_edge(size = 1) +
  # Edges' labels
  geom_edge_label(colour = "gray48", size = 4) +
  #geom_edge_label(mapping = aes(label = nodesize), colour = "gray48", size = 4) +
  # Boxplots for terminal nodes
  geom_node_plot(gglist = list(geom_boxplot(aes(x = '', y = resp)),
                               # 'resp' comes from 'nodedata_resp'
                               labs(x = '', y = 'Response'),
                               theme_bw(base_size = 10)),
                 ids = "terminal",
                 shared_axis_labels = TRUE) +
  # Inner nodes' labels
  geom_node_label(aes(col = splitvar),
                  # label nodes with ID, split variable and pvalue
                  line_list = list(aes(label = paste("Node", id)),
                                   aes(label = splitvar),
                                   aes(label = asterisk_sign(p.value))),
                  # set graphical parameters for each line
                  line_gpar = list(list(size = 8, col = "black", fontface = "bold"),
                                   list(size = 12),
                                   list(size = 8)),
                  ids = "inner") +
  # Terminal nodes' labels
  geom_node_label(aes(label = paste0("Node ", id, ", N = ", nodesize)),
                  fontface = "bold",
                  ids = "terminal",
                  size = 3,
                  # 0.01 nudge_y to be above the node plot
                  nudge_y = 0.01,
                  # 0.005 nudge_x to center terminal nodes' labels
                  nudge_x = 0.005) +
  theme(legend.position = "none")


# CV evaluation ---------------------------------------------------------------

# Loading caret (Classification And REgression Training) package
library(caret)

# Folds
set.seed(12345)
folds <- caret::createFolds(resp, k = 10, list = TRUE, returnTrain = FALSE)

# Cross-validation: fit without fold and predict on fold (x10)
e_cv <- lapply(folds,
               function(f){
                 set.seed(123)
                 etree_fit <- etree(response = resp[-f],
                                    covariates = lapply(cov.list,
                                                        function(c){
                                                          c[-f]
                                                        }),
                                    case.weights = NULL,
                                    minbucket = 5,
                                    alpha = 0.5,
                                    R = 1000,
                                    split.type = 'cluster',
                                    coef.split.type = 'test')
                 y_pred <- predict(etree_fit, newdata = lapply(cov.list,
                                                               function(c){
                                                                 c[f]
                                                               }))
                 list(e_fit = etree_fit,
                      e_pred = y_pred,
                      e_resp = resp[f],
                      e_summ = defaultSummary(data.frame(obs = resp[f],
                                                         pred = y_pred)))
               })

# Summary performance measures
(e_cv_summ <- sapply(e_cv, function(e) e$e_summ))

# Other performance metrics
library(MLmetrics)
e_perf <- sapply(e_cv,
                function(fold){
                  pred <- mean(fold$e_resp)
                  true <- fold$e_resp
                  out <- rbind(MAPE = MLmetrics::MAPE(pred, true),
                               RMSPE = MLmetrics::RMSPE(pred, true),
                               NRMSE = MLmetrics::RMSE(pred, true)/mean(resp))
                })
apply(e_perf, 1, mean)
