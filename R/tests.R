### Test each function

# compute.dissimilarity ---------------------------------------------------


xx =c(1:10)
yy =c(2:11)

library(fda.usc)
data(phoneme)
mlearn<-phoneme$learn[1:4,1:150]

compute.dissimilarity(mlearn, dist.type = "default", lp = 4, case.weights = c(3,2))



foo = function(x)


# mytestREG  --------------------------------------------------------------


mytestREG(mlearn, rbind(xx,xx,xx+2,xx))
compute.dissimilarity(rbind(xx,xx,xx+2,xx), dist.type = "default", lp = 4)




fda.usc::metric.dist(data.frame(yy))
dist(yy)



# findsplit ---------------------------------------------------------------

xxx = list(mlearn, 1:4, matrix(runif(16), 4,4))
yyy = runif(4)

findsplit(response = 1:4, covariates= xxx, case.weights=1:4,
                      alpha=1,
                      R=1000,
                      rnd.sel=1,
                      rnd.splt=1,
                      dist.types = rep("default",2),
                      lp = rep(2,2))


foo.x = as.factor(letters[sample(1:5, 10, replace =T)])

foo.lev = levels(foo.x)
length(foo.lev)

comb = do.call("c", lapply(1:4, function(x) combn(foo.lev, x, simplify = F)))


combn(foo.lev, 1, simplify = F)
combn(foo.lev, 2, simplify = F)


library(fda.usc)
data(phoneme)
mlearn<-phoneme$learn[c(1:50,101:150,201:250),]

mlearn[out.fd1$cluster==1]
out.fd1=kmeans.fd(mlearn,ncl=3,draw=TRUE)


