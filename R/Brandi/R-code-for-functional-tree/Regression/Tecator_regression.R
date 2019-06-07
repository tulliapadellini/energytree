library(fda.usc)
library(energy)
library(partykit)
data("tecator")


Y=tecator$y$Fat
func.data=tecator$absorp.fdata
data=list(Y=Y, X=func.data)

source("functionsREG.R")


myREG <- mytree(Y="Y", data=data,minbucket = 20)
plot(myREG)


test_y <- data$Y
f.test <- data$X
foo <- min.basis(f.test, numbasis = nb)
fd3 <- fdata2fd(foo$fdata.est, type.basis = "bspline", nbasis = foo$numbasis.opt)
m.coef <- data.frame(t(fd3$coefs))
for(j in 1: dim(m.coef)[2]){
  names(m.coef)[j] <- paste("X",names(m.coef)[j],sep = ".")
}
y_pred=predict(myREG, newdata = m.coef)
MEP <- (sum((test_y-y_pred)^2)/length(test_y))/(var(test_y))

