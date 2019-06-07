#Set the seed for simulation
set.seed(1234)

##########################
#Load the libraries

library(fda.usc)
library(roahd)
library(energy)
library(entropy)
library(partykit)
library(cvTools)
################
#Build the multivariate functional object ---<- define parameters

class  c("A","B","C")
size <- c(50,50,50)
alpha <- matrix( c(0.2, 0.5, 0.25, 0.7, 1, 0.1, 0.7, 0.45, 0.6),
                nrow=3, ncol=3)
beta <- matrix( c(0.3, 0.1, 0.6, 0.2, 0.6, 1, 0.1, 0.5, 0.9),
               nrow=3, ncol=3)
b <- matrix( c(2, 3, 1.5, 3, 1.5, 3.5, 2, 4, 1),
            nrow=3, ncol=3)
a <- matrix( c(0, 0, 1, 1, 1, 1, 2, 3, 2),
            nrow=3, ncol=3)
P <- 1e2
n.var <- 3

############
#Load the function that create the multivariate functional object

source("gen_data_3class.R")

#Generate data
d2=gen_data_3(class=class,size=size,P=P,n.var=n.var,alpha=alpha,beta=beta,
            a=a,b=b)

##########
##Load the classification tree function
source("functions10.R")


#############
#Generate folds for 10 fold cross validation

k <- 10 #the number of folds
set.seed(1234)
folds <- cvFolds(dim(d2$V1)[1], K=k)
acc1 <- c()
acc2 <- c()

##########
#Start the classification/validation cicle


for( i in 1:k){
  
  #Define test and train set
  
  list.train <- lapply(d2, function(x) {x[folds$subsets[folds$which != i]]})
  list.test <- lapply(d2, function(x) {x[folds$subsets[folds$which == i]]})
  
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
  acc1[i] <- sum(diag(t1))/(length(test_class))
  t2 <- table(predict.classif(a2, newdat,type = "class"), test_class)
  acc2[i] <- sum(diag(t2))/(length(test_class))
  


}


pdf(file = "Accuracy-PlotMultivariate.pdf")
plot(acc1, type = "l", ylim = c(0,1), ylab = "Accuracy",
     xlab = "Fold")
lines(acc2, col = "red")
legend(2, 0.5, c("MyTree","Classif.tree"), cex=0.8, fill = c("black","red"), border = "white",
       bty = "n")
dev.off()
