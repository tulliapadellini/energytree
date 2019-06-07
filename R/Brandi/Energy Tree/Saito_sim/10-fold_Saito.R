source("create_saito.R")

n <- 266
data <- shape.data(n.samp=c(n,n,n))
source("functions10.R")


library(cvTools)
library(fda.usc)
library(energy)
library(entropy)
library(partykit)


k <- 10 #the number of folds
set.seed(1234)
folds <- cvFolds(dim(data)[1], K=k)
acc1 <- c()
acc2 <- c()


for( i in 1:k){
  
  train  <- data[folds$subsets[folds$which != i],]
  test   <- data[folds$subsets[folds$which == i],]
  
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
  acc1[i] <- sum(diag(t1))/(nrow(test))
  t2 <- table(predict.classif(a2, newdat,type = "class"), test_class)
  acc2[i] <- sum(diag(t2))/(nrow(test))
  
}

pdf(file = "Accuracy-PlotSaito.pdf")
plot(acc1, type = "l", ylim = c(0,1), ylab = "Accuracy",
     xlab = "Fold")
lines(acc2, col = "red")
legend(2, 0.5, c("MyTree","Classif.tree"), cex=0.8, fill = c("black","red"), border = "white",
       bty = "n")
dev.off()
