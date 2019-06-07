#Load function to create Saito dataset
source("create_saito.R")

#Define the number of observations for each class
n <- 266
#data <- shape.data(n.samp=c(n,n,n))

#Load functions for mytree and packages
source("functions10.R")
library(fda.usc)
library(energy)
library(entropy)
library(partykit)

#data <- list()
#Number of dataset simulations
M=100

#for(i in 1:M){
#  data[[i]] <- shape.data(n.samp=c(n,n,n))
#}
#save(data, file="sim_data.RData")

#Load the simulated data
load("sim_data.RData")


#Define the accuracy slot
acc_my <- c()
acc_m1 <- c()


#start the simulation study

for(i in 1:M){


id.train <- sample(1:(3*n), size=0.7*3*n)
train <- data[[i]][id.train,]
test <- data[[i]][-id.train,]


nb <- 15
data_tot <- list(class=train$cls,V1=fdata(train[,2:129]))


####    MYTREE
myS  <-  mytree     (group = "class", data = data_tot, weights = NULL, 
                     minbucket = 20, 
                     alpha = 0.05, R = 1000, 
                     rnd.sel = T, rnd.splt = TRUE, nb=nb)

#### CLASSIF.TREE
class <- data_tot$class
dataf <- data.frame(class)
V1 <- data_tot$V1
dat <- list("df" = dataf, "x" = V1)
v1b <- create.fdata.basis(V1,l=1:15)
a2 <- classif.tree(class~x, data = dat,
                   basis.b = list(V1 = v1b))





###  CLASSIFY
test_class <- test$cls
f.test <- fdata(test[,2:129])
foo <- min.basis(f.test, numbasis = nb)
fd3 <- fdata2fd(foo$fdata.est, type.basis = "bspline", nbasis = foo$numbasis.opt)
m.coef <- data.frame(t(fd3$coefs))
for(j in 1: dim(m.coef)[2]){
  names(m.coef)[j] <- paste("V1",names(m.coef)[j],sep = ".")
}

newdat <- list("x" = fdata(test[2:129]))

t1 <- table(predict(myS, newdata = m.coef), test_class)
acc_my[i] <-  sum(diag(t1))/(nrow(test))
t2 <- table(predict.classif(a2, newdat,type = "class"), test_class)
acc_m1[i] <- sum(diag(t2))/(nrow(test))

}



