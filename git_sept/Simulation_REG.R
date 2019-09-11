##########################
#Load the libraries

library(fda.usc)
library(roahd)
library(energy)
library(entropy)
library(partykit)

#Load the classification tree functions
#source("functions_test.R")


################
#Build the  functional object ----- define parameters

size <- c(50,50,50,50)

P <- 50
n.var <- 1


############
#Load the function that create the multivariate functional object

source("gen_data_REG_CATS.R")
#Generate data
data=list()

M <- 1



for(i in 1:M){
  alpha <- matrix( round(runif(4,0.1,1),2),
                   nrow=4, ncol=1)
  beta <- matrix( round(runif(4,0.1,1),2),
                  nrow=4, ncol=1)
  data[[i]]=gen_data_reg(size=size,P=P,n.var=n.var,alpha=alpha,beta=beta,
                 sigma=1)


}

#save(data, file = "Regression_Sim_Dataset.RData")
load("Regression_Sim_Dataset.RData")

MEPmy <- c()

MEP_b <- c()

MEP_pc <- c()


for(i in 1:M){
  print(i)

  #select train e test set
  #id.train <- sample(1:(3*size[1]), size=0.7*3*size[1])
  #list.train <- lapply(data[[i]], function(x) {x[id.train]})
  #list.test <- lapply(data[[i]], function(x) {x[-id.train]})

  nb <- 15

  resp.list <- lapply(data, function(x) x$Y)
  cov.list <- lapply(data, function(x) x$V1)

  #REGRESSION FUNCTIONAL TREE
  myREG <- mytree(response = resp.list[[i]],
                  covariates = cov.list[i],
                  case.weights = NULL,
                  minbucket = 1,
                  alpha = 0.05,
                  R = 1000,
                  rnd.sel = T,
                  rnd.splt = TRUE,
                  nb = nb)
  plot(myREG)

  #FUNCTIONAL LINEAR MODEL BASIS
  x <- cov.list[[i]]
  y <- resp.list[[i]]
  tt=x[["argvals"]]

  dataf=as.data.frame(y)

  nbasis.x=15

  basis1=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)


  basis.x=list("x"=basis1)


  res=fregre.basis(x,y,basis.x=basis1)



  #########FUNCTIONAL LINEAR MODEL --- PC


  res2 <- fregre.pc(x, y, kmax = 7)




  ###################
  ##PREDICTION MYREG##


  foo <- min.basis(cov.list[[i]], numbasis = nb)
  fd3 <- fdata2fd(foo$fdata.est, type.basis = "bspline", nbasis = foo$numbasis.opt)
  m.coef <- data.frame(t(fd3$coefs))
  for(j in 1: dim(m.coef)[2]){
    names(m.coef)[j] <- j
  }
  y_pred=predict(myREG, newdata = m.coef)
  MEPmy[i] <- (sum((y-y_pred)^2)/length(y))/(var(y))

  #################
  ##PREDICTION FLM

  y_pred1 <- predict.fregre.fd(res, new.fdataobj=data[[i]]$V1)
  MEP_b[i] <- (sum((y-y_pred1)^2)/length(y))/(var(y))


  ###PREDICTION PC
  y_pred2 <- predict.fregre.fd(res2,data[[i]]$V1)
  MEP_pc[i] <- (sum((y-y_pred2)^2)/length(y))/(var(y))



}

save(MEPmy,MEP_b,MEP_pc, file="results.RData")
