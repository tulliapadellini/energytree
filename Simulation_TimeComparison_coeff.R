
# Inizialization ----------------------------------------------------------

# Loading the libraries
library(fda.usc)
library(roahd)
library(energy)
library(entropy)
library(partykit)
library(igraph)
library(NetworkDistance)
library(microbenchmark)

# Loading the dataset
load("Regression_Sim_Dataset.RData")



# Datasets ----------------------------------------------------------------

norm01 <- function(x){(x-min(x))/(max(x)-min(x))}
n.bas <- 15

### DATASET 1 (20 observations) ###

# Response
resp1 <- lapply(data, function(x) x$Y[1:20])[[1]]

# Building the numerical variable
foo <- fda.usc::min.basis(lapply(data, function(x) x$V1)[[2]], numbasis = 15)
fd3 <- fda.usc::fdata2fd(foo$fdata.est,
                         type.basis = "bspline",
                         nbasis = foo$numbasis.opt)
foo$coef <- t(fd3$coefs)

# Building the graph variable
conn.prob <- norm01(resp1)
conn.prob[which(conn.prob < 0.05)] <- conn.prob[which(conn.prob < 0.05)] + 0.1
set.seed(12345)
graph.list <- lapply(conn.prob, function(p){sample_gnp(100, p)})

# Covariates for the full mixed model
cov.list1 <- list(lapply(data, function(x) x$V1[1:20])[[1]], graph.list, foo$coef[1:20,15])


### DATASET 2 (100 observations) ###

# Response
resp2 <- lapply(data, function(x) x$Y[1:100])[[1]]

# Building the numerical variable
foo <- fda.usc::min.basis(lapply(data, function(x) x$V1)[[2]], numbasis = 15)
fd3 <- fda.usc::fdata2fd(foo$fdata.est,
                         type.basis = "bspline",
                         nbasis = foo$numbasis.opt)
foo$coef <- t(fd3$coefs)

# Building the graph variable
conn.prob <- norm01(resp2)
conn.prob[which(conn.prob < 0.05)] <- conn.prob[which(conn.prob < 0.05)] + 0.1
set.seed(12345)
graph.list <- lapply(conn.prob, function(p){sample_gnp(100, p)})

# Covariates for the full mixed model
cov.list2 <- list(lapply(data, function(x) x$V1[1:100])[[1]], graph.list, foo$coef[1:100,15])


### DATASET 3 (200 observations) ###

# Response
resp3 <- lapply(data, function(x) x$Y)[[1]]

# Building the numerical variable
foo <- fda.usc::min.basis(lapply(data, function(x) x$V1)[[2]], numbasis = 15)
fd3 <- fda.usc::fdata2fd(foo$fdata.est,
                         type.basis = "bspline",
                         nbasis = foo$numbasis.opt)
foo$coef <- t(fd3$coefs)

# Building the graph variable
conn.prob <- norm01(resp3)
conn.prob[which(conn.prob < 0.05)] <- conn.prob[which(conn.prob < 0.05)] + 0.1
set.seed(12345)
graph.list <- lapply(conn.prob, function(p){sample_gnp(100, p)})

# Covariates for the full mixed model
cov.list3 <- list(lapply(data, function(x) x$V1)[[1]], graph.list, foo$coef[,15])



# functions ---------------------------------------------------------------
## distances & expansions inside
source("functions.R")

mb_co_fun_1 <- microbenchmark(etree(response = resp1,
                                 covariates = cov.list1,
                                 case.weights = NULL,
                                 minbucket = 5,
                                 alpha = 0.05,
                                 R = 1000,
                                 split.type = 'coeff',
                                 coef.split.type = 'test',
                                 nb = n.bas),
                           times = 10L)

mb_co_fun_2 <- microbenchmark(etree(response = resp2,
                                 covariates = cov.list2,
                                 case.weights = NULL,
                                 minbucket = 5,
                                 alpha = 0.05,
                                 R = 1000,
                                 split.type = 'coeff',
                                 coef.split.type = 'test',
                                 nb = n.bas),
                           times = 10L)

mb_co_fun_3 <- microbenchmark(etree(response = resp3,
                                 covariates = cov.list3,
                                 case.weights = NULL,
                                 minbucket = 5,
                                 alpha = 0.05,
                                 R = 1000,
                                 split.type = 'coeff',
                                 coef.split.type = 'test',
                                 nb = n.bas),
                           times = 10L)


# fun_test1 ---------------------------------------------------------------
## distances inside, expansions outside
source("fun_test1.R")

mb_co_ft1_1 <- microbenchmark(etree(response = resp1,
                                 covariates = cov.list1,
                                 case.weights = NULL,
                                 minbucket = 5,
                                 alpha = 0.05,
                                 R = 1000,
                                 split.type = 'coeff',
                                 coef.split.type = 'test',
                                 nb = n.bas),
                           times = 10L)

mb_co_ft1_2 <- microbenchmark(etree(response = resp2,
                                 covariates = cov.list2,
                                 case.weights = NULL,
                                 minbucket = 5,
                                 alpha = 0.05,
                                 R = 1000,
                                 split.type = 'coeff',
                                 coef.split.type = 'test',
                                 nb = n.bas),
                           times = 10L)

mb_co_ft1_3 <- microbenchmark(etree(response = resp3,
                                 covariates = cov.list3,
                                 case.weights = NULL,
                                 minbucket = 5,
                                 alpha = 0.05,
                                 R = 1000,
                                 split.type = 'coeff',
                                 coef.split.type = 'test',
                                 nb = n.bas),
                           times = 10L)


# fun_test2 ---------------------------------------------------------------
## variable selection distances & expansions outside
source("fun_test2.R")

mb_co_ft2_1 <- microbenchmark(etree(response = resp1,
                                 covariates = cov.list1,
                                 case.weights = NULL,
                                 minbucket = 5,
                                 alpha = 0.05,
                                 R = 1000,
                                 split.type = 'coeff',
                                 coef.split.type = 'test',
                                 nb = n.bas),
                           times = 10L)

mb_co_ft2_2 <- microbenchmark(etree(response = resp2,
                                 covariates = cov.list2,
                                 case.weights = NULL,
                                 minbucket = 5,
                                 alpha = 0.05,
                                 R = 1000,
                                 split.type = 'coeff',
                                 coef.split.type = 'test',
                                 nb = n.bas),
                           times = 10L)

mb_co_ft2_3 <- microbenchmark(etree(response = resp3,
                                 covariates = cov.list3,
                                 case.weights = NULL,
                                 minbucket = 5,
                                 alpha = 0.05,
                                 R = 1000,
                                 split.type = 'coeff',
                                 coef.split.type = 'test',
                                 nb = n.bas),
                           times = 10L)


# fun_test3 ---------------------------------------------------------------
## (variable selection + clustering) distances & expansions outside
source("fun_test3.R")

mb_co_ft3_1 <- microbenchmark(etree(response = resp1,
                                 covariates = cov.list1,
                                 case.weights = NULL,
                                 minbucket = 5,
                                 alpha = 0.05,
                                 R = 1000,
                                 split.type = 'coeff',
                                 coef.split.type = 'test',
                                 nb = n.bas),
                           times = 10L)

mb_co_ft3_2 <- microbenchmark(etree(response = resp2,
                                 covariates = cov.list2,
                                 case.weights = NULL,
                                 minbucket = 5,
                                 alpha = 0.05,
                                 R = 1000,
                                 split.type = 'coeff',
                                 coef.split.type = 'test',
                                 nb = n.bas),
                           times = 10L)

mb_co_ft3_3 <- microbenchmark(etree(response = resp3,
                                 covariates = cov.list3,
                                 case.weights = NULL,
                                 minbucket = 5,
                                 alpha = 0.05,
                                 R = 1000,
                                 split.type = 'coeff',
                                 coef.split.type = 'test',
                                 nb = n.bas),
                           times = 10L)


# fun_test4 ---------------------------------------------------------------
## (variable selection + clustering + expansions) distances & expansions outside
source("fun_test4.R")

mb_co_ft4_1 <- microbenchmark(etree(response = resp1,
                                 covariates = cov.list1,
                                 case.weights = NULL,
                                 minbucket = 5,
                                 alpha = 0.05,
                                 R = 1000,
                                 split.type = 'coeff',
                                 coef.split.type = 'test',
                                 nb = n.bas),
                           times = 10L)

mb_co_ft4_2 <- microbenchmark(etree(response = resp2,
                                 covariates = cov.list2,
                                 case.weights = NULL,
                                 minbucket = 5,
                                 alpha = 0.05,
                                 R = 1000,
                                 split.type = 'coeff',
                                 coef.split.type = 'test',
                                 nb = n.bas),
                           times = 10L)

mb_co_ft4_3 <- microbenchmark(etree(response = resp3,
                                 covariates = cov.list3,
                                 case.weights = NULL,
                                 minbucket = 5,
                                 alpha = 0.05,
                                 R = 1000,
                                 split.type = 'coeff',
                                 coef.split.type = 'test',
                                 nb = n.bas),
                           times = 10L)



# Comparison --------------------------------------------------------------

et_co <- rbind(mb_co_fun_1$time,
               mb_co_fun_2$time,
               mb_co_fun_3$time,
               mb_co_ft1_1$time,
               mb_co_ft1_2$time,
               mb_co_ft1_3$time,
               mb_co_ft2_1$time,
               mb_co_ft2_2$time,
               mb_co_ft2_3$time,
               mb_co_ft3_1$time,
               mb_co_ft3_2$time,
               mb_co_ft3_3$time,
               mb_co_ft4_1$time,
               mb_co_ft4_2$time,
               mb_co_ft4_3$time)

et_co <- et_co/1000000000

summary(et_co)

save(mb_co_fun_1,
     mb_co_fun_2,
     mb_co_fun_3,
     mb_co_ft1_1,
     mb_co_ft1_2,
     mb_co_ft1_3,
     mb_co_ft2_1,
     mb_co_ft2_2,
     mb_co_ft2_3,
     mb_co_ft3_1,
     mb_co_ft3_2,
     mb_co_ft3_3,
     mb_co_ft4_1,
     mb_co_ft4_2,
     mb_co_ft4_3,
     file = 'TimeComparison_coeff.RData')

### Legend ###
# co: coeff
# fun: distances & expansions inside
# ft1: distances inside, expansions outside
# ft2: variable selection distances & expansions outside
# ft3: (variable selection + clustering) distances & expansions outside
# ft4: (variable selection + clustering + expansions) distances & expansions outside
# _1: small dataset (20 observations)
# _2: medium dataset (100 observations)
# _3: large dataset (200 observations)
