##########################
#Load the libraries

library(fda.usc)
library(roahd)
library(energy)
library(entropy)
library(partykit)

##Load the classification tree function
source("functions10.R")

################
#Build the multivariate functional object ----- define parameters

class <- c("A","B","C")
size <- c(50,50,50)

P <- 1e2
n.var <- 3


############
#Load the function that create the multivariate functional object

source("gen_data_3class.R")
#Generate data
#data=list()

M=100

#for(i in 1:M){
#  alpha <- matrix( round(runif(9,0.1,1),2),
#                   nrow=3, ncol=3)
#  beta <- matrix( round(runif(9,0.1,1),2),
#                  nrow=3, ncol=3)
#  b <- matrix( sample(c(1,1.5,2,2.5,3,3.5,4),size=9, replace=T),
#               nrow=3, ncol=3)
#  a <- matrix( sample(c(0,1,2,3),size=9, replace=T),
#               nrow=3, ncol=3)
#  data[[i]]=gen_data_3(class=class,size=size,P=P,n.var=n.var,alpha=alpha,beta=beta,
#                a=a,b=b)
#  
#    
#}

#save(data, file = "multivariate_dataset.RData")
load("multivariate_dataset.RData")

accMy <- c()

accM1 <- c()

for(i in 1:M){
  print(i)
  #select train e test set
  id.train <- sample(1:(3*size[1]), size=0.7*3*size[1])
  list.train <- lapply(data[[i]], function(x) {x[id.train]})
  list.test <- lapply(data[[i]], function(x) {x[-id.train]})
  
  nb <- 15
  
  
  
  
  ##MYTREE
  myt  <-  mytree     (group="class", data=list.train, weights = NULL, 
                       minbucket = 1, 
                       alpha = 0.05, R = 1000, 
                       rnd.sel = T, rnd.splt = TRUE, nb=nb)
  
  ##CLASSIF.TREE
  class=list.train$class
  dataf=data.frame(class)
  V1=list.train$V1
  V2=list.train$V2
  V3=list.train$V3
  dat=list("df"=dataf,"x"=V1,"y"=V2,"z"=V3)
  a2<-classif.tree(class~x+y+z,data=dat)
  
  
  ##CLASSIFICATION
  test_class <- list.test$class
  ###Transform test data into functional object and take coefficient for prediction
  
  ##V1
  foo <- min.basis(list.test$V1, numbasis = nb)
  fd3 <- fdata2fd(foo$fdata.est, type.basis = "bspline", nbasis = foo$numbasis.opt)
  m.coef1 <- data.frame(t(fd3$coefs))
  for(j in 1: dim(m.coef1)[2]){
    names(m.coef1)[j] <- paste("V1",names(m.coef1)[j],sep = ".")
  }
  
  ##V2
  foo <- min.basis(list.test$V2, numbasis = nb)
  fd3 <- fdata2fd(foo$fdata.est, type.basis = "bspline", nbasis = foo$numbasis.opt)
  m.coef2 <- data.frame(t(fd3$coefs))
  for(j in 1: dim(m.coef2)[2]){
    names(m.coef2)[j] <- paste("V2",names(m.coef2)[j],sep = ".")
  }
  
  ##V3
  foo <- min.basis(list.test$V3, numbasis = nb)
  fd3 <- fdata2fd(foo$fdata.est, type.basis = "bspline", nbasis = foo$numbasis.opt)
  m.coef3 <- data.frame(t(fd3$coefs))
  for(j in 1: dim(m.coef3)[2]){
    names(m.coef3)[j] <- paste("V3",names(m.coef3)[j],sep = ".")
  }
  
  
  ##Create the matrix of all coefficients
  m.coef <- cbind(m.coef1,m.coef2,m.coef3)
  
  
  
  #newdata for classif.tree
  newdat <- list("x"=list.test$V1, "y"=list.test$V2, "z"=list.test$V3)
  
  
  t1 <- table(predict(myt, newdata = m.coef), test_class)
  accMy[i] <- sum(diag(t1))/(length(test_class))
  t2 <- table(predict.classif(a2, newdat,type = "class"), test_class)
  accM1[i] <- sum(diag(t2))/(length(test_class))
  
  
}

###F-GLM
acc_m2 <- c()
for(i in 1:M){
  print(i)
  #select train e test set
  id.train <- sample(1:(3*size[1]), size=0.7*3*size[1])
  list.train <- lapply(data[[i]], function(x) {x[id.train]})
  list.test <- lapply(data[[i]], function(x) {x[-id.train]})
  
  nb <- 15
  class=list.train$class
  class=class[order(class)]
  dataf=data.frame(class)
  V1=list.train$V1[order(class)]
  V2=list.train$V2[order(class)]
  V3=list.train$V3[order(class)]
  dat=list("df"=dataf,"x"=V1,"y"=V2,"z"=V3)
  a3<-classif.glm(class~x+y+z,data=dat)
  #newdata for classif.glm
  newdat <- list("x"=list.test$V1, "y"=list.test$V2, "z"=list.test$V3)
  p1 = round( 
    predict.classif(a3, newdat, type = "probs")$prob.group,
    2) ##solo NA!!!
  c1=factor(colnames(p1)[apply(p1, 1, which.max)])
  t3 <- table(c1, list.test$class)
  acc_m2[i] <- sum(diag(t3))/(length(list.test$class))
  
}
