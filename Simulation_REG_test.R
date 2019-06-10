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
  #FUNCTIONAL LINEAR MODEL
  x <- list.train$V1
  y <- list.train$Y
  tt=x[["argvals"]]

  dataf=as.data.frame(y)

  nbasis.x=nb
  nbasis.b=nb
  basis1=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)
  basis2=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.b)

  f=y~x
  basis.x=list("x"=basis1)
  basis.b=list("x"=basis2)
  ldata=list("df"=dataf,"x"=x)
  res=fregre.lm(f,ldata,basis.x=basis.x,basis.b=basis.b)

  ###################
  ##PREDICTION MYREG##

  test_y <- list.test$Y
  f.test <- list.test$V1
  foo <- min.basis(f.test, numbasis = nb)
  fd3 <- fdata2fd(foo$fdata.est, type.basis = "bspline", nbasis = foo$numbasis.opt)
  m.coef <- data.frame(t(fd3$coefs))
  for(j in 1: dim(m.coef)[2]){
    names(m.coef)[j] <- paste("V1",names(m.coef)[j],sep = ".")
  }
  y_pred=predict(myREG, newdata = m.coef)
  MEPmy[i] <- (sum((test_y-y_pred)^2)/length(test_y))/(var(test_y))

  #################
  ##PREDICTION FLM
  newldata=list( "x"=list.test$V1)
  y_pred1 <- predict.fregre.lm(res, newx=newldata)
  MEP1[i] <- (sum((test_y-y_pred1)^2)/length(test_y))/(var(test_y))


}

save(MEPmy,MEP1, file="results.RData")
