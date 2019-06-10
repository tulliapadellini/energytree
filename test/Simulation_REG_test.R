##########################
#Load the libraries

library(fda.usc)
library(roahd)
library(energy)
library(entropy)
library(partykit)

##Load the classification tree function
source("functionsREG.R")


################
#Build the  functional object ----- define parameters

size <- c(50,50,50)

P <- 1e2
n.var <- 1


############
#Load the function that create the multivariate functional object

source("gen_data_REG.R")
#Generate data
data=list()

M=100
y_mean=c(10,30,50)
y_sd=c(2,3,2.5)


for(i in 1:M){
 alpha <- matrix( round(runif(3,0.1,1),2),
                  nrow=3, ncol=1)
 beta <- matrix( round(runif(3,0.1,1),2),
                 nrow=3, ncol=1)
 b <- matrix( sample(c(1,1.5,2,2.5,3,3.5,4),size=3, replace=T),
              nrow=3, ncol=1)
 a <- matrix( sample(c(0,1,2,3),size=3, replace=T),
              nrow=3, ncol=1)
 data[[i]]=gen_data_reg(y_mean, y_sd ,size=size,P=P,n.var=n.var,alpha=alpha,beta=beta,
               a=a,b=b)


}

#save(data, file = "Regression_Sim_Dataset.RData")
load("Regression_Sim_Dataset.RData")

MEPmy <- c()

MEP1 <- c()

for(i in 1:M){
  print(i)

  #select train e test set
  id.train <- sample(1:(3*size[1]), size=0.7*3*size[1])
  list.train <- lapply(data[[i]], function(x) {x[id.train]})
  list.test <- lapply(data[[i]], function(x) {x[-id.train]})

  nb <- 15

  #REGRESSION FUNCTIONAL TREE
  myREG <- mytree(list.train$Y, list.train$V1, case.weights = NULL,
                  minbucket = 10,
                  alpha = 0.05, R = 1000,
                  rnd.sel = T, rnd.splt = TRUE, nb=nb)
  plot(myREG)
}

