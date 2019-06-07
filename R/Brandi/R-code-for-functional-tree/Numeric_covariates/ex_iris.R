library(energy)
library(entropy)
library(partykit)
set.seed(1234)

data=iris
source("functions11.R")

myS  <-  mytree     (group = "Species", data = data, weights = NULL, 
                     minbucket = 1, 
                     alpha = 0.05, R = 1000, 
                     rnd.sel = T, rnd.splt = TRUE, nb=15)

