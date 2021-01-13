
# Load packages, functions and data -------------------------------------------

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

# Data (following CATS by Serban and Wasserman, 2005)
load("sim_data.RData")


# Response and covariates -----------------------------------------------------

# Restriction on the observations' number (for computational reasons)
obs <- sample(1:798, 150)

# Response
resp <- lapply(data, function(x) x$cls[obs])[[1]]
#remark: we take the first element since the dataset contains 100 simulations

### Numerical covariate ###

# One of the basis of another simulation of the same dataset
foo <- fda.usc::optim.basis(lapply(data, function(x) fdata(x[obs,2:129]))[[2]], numbasis = 15)
fd3 <- fda.usc::fdata2fd(foo$fdata.est,
                         type.basis = "bspline",
                         nbasis = foo$numbasis.opt)
foo$coef <- t(fd3$coefs)
num.cov <- foo$coef[,3]
#num.cov <- apply(foo$coef, 1, mean)

### Graph covariate ###

# Generation of the graphs with a different connection probability for each class
graph.list <- lapply(resp,
                     function(c){
                       if (c == 'Bel'){
                         sample_gnp(100, 0.10)
                       } else if (c == 'Cyl'){
                         sample_gnp(100, 0.125)
                       } else if (c == 'Fun'){
                         sample_gnp(100, 0.15)
                       }
                     })

### Functional covariate ###

# Functional data for the first simulation
fun.list <- lapply(data, function(x) fdata(x[obs,2:129]))[[1]]

### Covariates list ###

# Covariates list: functional, graph, numeric
cov.list <- list(fun.list, graph.list, num.cov)


# Model fitting -----------------------------------------------------------

# Number of basis
n.bas <- 15

### Classification Energy Tree ###
set.seed(123)
etree_fit <- etree(response = resp,
                   covariates = cov.list,
                   case.weights = NULL,
                   minbucket = 5,
                   alpha = 0.05,
                   R = 500,
                   split.type = 'coeff',
                   coef.split.type = 'test',
                   p.adjust.method = 'fdr',
                   nb = n.bas)
plot(etree_fit)


# Prediction --------------------------------------------------------------

### ETREE CLASSIFICATION PREDICTION ###

# Prediction
y_pred <- predict(etree_fit)

# Prediction with newdata
graph.list2 <- lapply(resp,
                      function(c){
                        if (c == 'Bel'){
                          sample_gnp(100, 0.10)
                        } else if (c == 'Cyl'){
                          sample_gnp(100, 0.125)
                        } else if (c == 'Fun'){
                          sample_gnp(100, 0.15)
                        }
                      })
new.cov.list <- list(lapply(data, function(x) fdata(x[obs,2:129]))[[2]], graph.list2, foo$coef[,8])
y_pred2 <- predict(etree_fit, newdata = new.cov.list)

# Error
y <- resp
t <- table(y_pred, y)
ACC_etree <- sum(diag(t))/(length(y))


# ggparty ----------------------------------------------------------------------

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

# Change name to the variables so that they are numbered (instead of b.spl)
gg$data$splitvar <- c(2L, 1L, 1L, NA, NA, NA, 3L, NA, 1L, NA, NA)

# Tedious round of the split points
gg$data$breaks_label <- lapply(gg$data$breaks_label,
                               function(s){
                                 #take first part of breaks_label
                                 paste0(sub('\\* .*', '', s),
                                        #add central part
                                        '* ',
                                        #add rounded split point (3 digits)
                                        round(as.numeric(sub('.*\\* ', '', s)), 3))
                               })



#for(idx in c(3,4,7,8,10)){
#  gg$data$breaks_label[[idx]] <- paste('n =', gg$data$nodesize[idx])
#}


gg +
  # Plot edges
  geom_edge(size = 1) +
  # Edges' labels
  geom_edge_label(colour = "gray48", size = 4) +
  #geom_edge_label(mapping = aes(label = nodesize), colour = "gray48", size = 4) +
  # Boxplots for terminal nodes
  geom_node_plot(gglist = list(geom_bar(aes(x = '', fill = resp),
                                        position = position_fill()),
                               # 'resp' comes from 'nodedata_resp'
                               labs(x = '', y = 'Response', fill = ''),
                               theme_bw(base_size = 12), scale_fill_brewer(palette='YlGn', type = 'div')),
                 ids = "terminal",
                 shared_axis_labels = TRUE) +
  # Inner nodes' labels
  geom_node_label(aes(col = factor(splitvar)),
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
                  nudge_x = 0.015) +
  theme(legend.position = "none")
