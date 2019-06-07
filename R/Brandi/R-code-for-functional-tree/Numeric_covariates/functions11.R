
#### ver. "2018-01-24 12:20:14 CET"

# perform 2 sample test of x
mytest <- function(x, y, R = 1000){
  library(cluster)
  library(fda.usc)
  if( is.factor(x) ){
    d1 = daisy(as.data.frame(x))
  }
  if( is.numeric(x) ){
    d1 = daisy(as.data.frame(x))
  }
  if( is.fdata(x) ){
    d1 = metric.lp(x)
  }
  ct <- eqdist.etest(d1, summary(y), distance=T , method = "discoB", R = R)
  if ( !is.na(ct$statistic) ){
    ct$p.value
  } else{
    NA
  }
}

mytest2 <- function(x, y, R = 1000){
  library(cluster)
  library(fda.usc)
  if( is.factor(x) ){
    d1 = daisy(as.data.frame(x))
  }
  if( is.numeric(x) ){
    d1 = daisy(as.data.frame(x))
  }
  if( is.fdata(x) ){
    d1 = metric.lp(x)
  }
  ct <- eqdist.etest(d1, summary(y), distance=T , method = "discoB", R = R)
  if ( !is.na(ct$statistic) ){
    return(c(ct$statistic, ct$p.value))
  } else{
    NA
  }
}


##test for split nominal variables

##perform chi-squared test of yvs.x
mychisqtest<-function(x,y){
  x <- factor(x)
  if(length(levels(x))<2)return(NA)
  ct <- suppressWarnings(chisq.test(table(y,x), correct = FALSE))
  pchisq(ct$statistic, ct$parameter, log = TRUE, lower.tail = FALSE)
}

####

# See sample()'s surprise -- example in help file
resample <- function(x, ...) x[sample.int(length(x), ...)]

####

mytree <- function(group, data, weights = NULL, 
                   minbucket = 1, 
                   alpha = 0.05, R = 1000, 
                   rnd.sel = T, rnd.splt = TRUE, nb=5) {
  
  # name of the response variable
  #response <- all.vars(formula)[1]
  response <- data[[which(names(data)==group)]]
  index  <- order(response)
  
  #names(response)=rownames(data$V2$data)
  # data without missing values, response comes last
  #data <- data[ complete.cases(data), c(all.vars(terms(formula, data = data))[-1], response) ]
  
  if ( is.null(weights) ) weights <- rep(1L,length(response))
  
  
  n.var <- which(names(data)!=group)
  # weights are case weights, i.e.,integers
  #stopifnot( length(weights) == nrow(data) &
  #             max( abs(weights - floor(weights)) ) < .Machine$double.eps )
  
  
  if(class(data)=="list"){
  datanew <- list("class"=response[index])
  for(j in n.var){
    if(class(data[[j]])=="fdata"){
      foo <- min.basis(data[[j]][index], numbasis = nb)
    fd3 <- fdata2fd(foo$fdata.est, type.basis = "bspline", nbasis = foo$numbasis.opt)
    foo$coef <- t(fd3$coefs)
    datanew[[j]] <- foo
    }
  }
  names(datanew) <- names(data)
  }
  
  if(class(data)=="data.frame"){
    datanew <- list("class"=response[index])
    for(j in n.var){
      datanew[[j+1]] <- data[index,j]
    names(datanew)[j+1] <- names(data)[j]
    }
  }
  
  
  
  
  
  # growtree
  nodes <- growtree( id = 1L, response=datanew$class, data = datanew, weights, 
                     minbucket = minbucket, 
                     alpha = alpha, R = R, 
                     rnd.sel = rnd.sel, rnd.splt = rnd.splt,
                     n.var = n.var)
  
  # compute terminal node number for each observation
  response <- response[index]
  response=data.frame(response)
  y=response
  #m.coef <- c()
  #for(j in n.var){
  #  foo <- datanew[[j]]$coef
  #  colnames(foo) <- paste(names(data)[j], colnames(datanew[[j]]$coef), sep = ".")
  #  m.coef <- cbind(m.coef,foo)
  #}
  #data1=cbind(response,m.coef)
  #m.coef=m.coef[order(index),]
  #data1=data1[order(index),]
  
  fitted <- fitted_node(nodes, data = data.frame(datanew))
  formula=datanew$class~.
  
  # return rich constparty object
  ret <- party(nodes, data = data.frame(datanew), 
               fitted = data.frame("(fitted)" = fitted, 
                                   "(response)" = datanew$class,
                                   "(weights)" = weights,
                                   check.names = FALSE),
               terms = terms(formula,data=data.frame(datanew))
  )
  as.constparty(ret)
}

####

growtree <- function(id=1L, response, data, weights, 
                     minbucket, 
                     alpha, R, 
                     rnd.sel, rnd.splt,n.var) {
  
  
  # for less than <minbucket> observations stop here (for ctree() is 7 in ?ctree_control)
  if (sum(weights) <= minbucket) return( partynode(id = id) )
  # Preventive stop :: if the best split would create small children, stop here
  yt <- table(as.factor( rep(response, weights) ) )
  #if (any(yt <= minbucket)) return( partynode(id = id) )
  
  if(length(which(yt==0))+1==length(levels(response))) return( partynode(id = id) )
  # find best split
  res <- findsplit( response, data, weights, 
                    alpha = alpha, R = R, 
                    rnd.sel = rnd.sel, rnd.splt = rnd.splt,
                    n.var = n.var )
  if(class(res)!="partysplit"){
  sp <- res$sp
  varselect <- res$varselect
  }else{
    sp <- res
    varselect <- res$varid
  }
  # no split found, stop here
  if (is.null(sp)) return( partynode(id = id) )
  
 
  
  kidids <- c()
  if(class(data[[varselect]])=="fdata"){
  kidids[which(data[[varselect]]$coef[,sp$varid]<=sp$breaks)] <- 1
  kidids[which(data[[varselect]]$coef[,sp$varid]>sp$breaks)] <- 2
  
  if(all(kidids==1) | all(kidids==2)) return( partynode(id = id) )
  sum1 <- length(which(data[[varselect]]$coef[which(weights==1),
                                              sp$varid]<=sp$breaks))
  sum2 <- length(which(data[[varselect]]$coef[which(weights==1),
                                              sp$varid]>sp$breaks))
  if(sum1 == 0 | sum2 == 0) return( partynode(id = id) )
  nb=0
  for(i in n.var){
    k=data[[i]]$numbasis.opt
    nb=c(nb,k)
  }
  
  if(varselect!=min(n.var)){
    step=sum(nb[1:(max(n.var[which(n.var<varselect)]))])
    sp$varid=sp$varid+as.integer(step)
  }
  }else{
    kidids[which(data[[varselect]]<sp$breaks)] <- 1
    kidids[which(data[[varselect]]>=sp$breaks)] <- 2
    
    if(all(kidids==1) | all(kidids==2)) return( partynode(id = id) )
    sum1 <- length(which(data[[varselect]][which(weights==1)]<sp$breaks))
    sum2 <- length(which(data[[varselect]][which(weights==1)]>=sp$breaks))
    if(sum1 == 0 | sum2 == 0) return( partynode(id = id) )
  }
  
  
  # setup all daugther nodes
  kids <- vector(mode = "list", length = max(kidids, na.rm = TRUE) )
  
  for ( kidid in 1:length(kids) ){
    # select observations for current node
    w <- weights
    w[kidids != kidid] <- 0
    # get next node id
    if(kidid > 1){
      myid <- max( nodeids( kids[[kidid - 1]] ) )
    } else{
      myid <- id
    }
    # Start recursion on this daugther node
    kids[[kidid]] <- growtree( id = as.integer(myid + 1), response, data, w, 
                               minbucket, 
                               alpha, R, 
                               rnd.sel, rnd.splt , n.var = n.var)
  }
  
  # return nodes
  if(class(res)!="partysplit"){
  return( partynode( id = as.integer(id), split = sp, kids = kids,
                     info = list( p.value = min( info_split(sp)$p.value, na.rm = TRUE ) )
  ) )
  }else{
    return( partynode( id = as.integer(id), split = res, kids = kids,
                       info = list( p.value = min( info_split(res)$p.value, na.rm = TRUE ) )
    ) )
  }
}

####

findsplit <- function( response, data, weights, 
                       alpha, R, 
                       rnd.sel, rnd.splt, n.var ) {
  
  
  # extract response values from data / as.factor retains all levels
  y <- response[which(weights==1)]
  # something to split?
  #if (!all(summary(y) > 0)) return(NULL)
  
  #data1    <- data[which(weights==1),]
  #data     <- data1[ order(rep(group, weights)), ]      # stack vals for each class one under the other
  #val      <- val[which(weights==1)]
  #xselect  <- which( names(data) != response )
  #p        <- sapply( xselect, function(i) mytest( data[[ i ]], y, R = R) )
  #names(p) <- names(data)[ xselect ]
  #p         <- mytest(x=data,y=y,R=R)
  nv=n.var+1
  p <- matrix(NA, nrow = length(n.var), ncol=2)
  colnames(p) <- c("statistic","p-value")
  
  for( i in nv){
    if(class(data[[i]])=="fdata"){
  p[i-1,] <- mytest2(x=data[[i]]$fdata.est[which(weights==1)],y=y,R=R)
    }else{
      p[i-1,] <- mytest2(x=data[[i]][which(weights==1)], y=y, R=R)
    }
  }
  
  # Bonferroni-adjusted p-value small enough?
  if ( all(is.na(p[,2])) ) return(NULL)
  
  minp <- min(p[,2], na.rm = TRUE)
  minp <- 1 - (1 - minp)^sum( !is.na(p[,2]) )
  if ( minp > alpha ) return(NULL)
  
  # for selected variable, search for split minimizing p-value
  #if (!rnd.sel) {
  #  xselect <- which.min(p)  # not randomize
  #} else {
  #  xselect <- resample(which(p == min(p, na.rm = T)), 1)   # randomize
  #}
  xselect <- which(names(data)!="class")
  if(length(which(p[,2] == min(p[,2], na.rm = T)))>1){
    xselect <- which.max(p[,1])+1
  }else{
    xselect <- which.min(p[,2])+1
  }
  #print( c( round(minp,5), round(min(p),5), names(p)[xselect] ) )
  x <-  data[[xselect]]
  if(is.list(x)){
    if(is.fdata(x$fdata.est)) x1=x$coef[which(weights==1),]
  } 
  
  # split into two groups minimizing entropy
  if(is.factor(x)){
    ##setup all possible splits in two kid nodes
    lev<-levels(x[drop = TRUE])
    if(length(lev) == 2){
      splitpoint <- lev[1]
    }else{
      comb <- do.call("c", lapply(1:(length(lev)-2),
                                  function(x) combn(lev, x, simplify = FALSE)))
      xlogp <- sapply(comb, function(q) mychisqtest(x%in%q, y))
      splitpoint <- comb[[which.min(xlogp)]]
    }
    
    ##split into two groups (setting groups that do not occur to NA)
    splitindex <- !(levels(data[[xselect]])%in%splitpoint)
    splitindex[!(levels(data[[xselect]])%in%lev)] <- NA_integer_
    splitindex <- splitindex-min(splitindex, na.rm=TRUE)+1L
  }
  if(is.numeric(x)){
    splitindex <- s.opt(response, x[which(weights==1)], rnd.splt)
  }
  if(is.list(x)){
    if(is.fdata(x$fdata.est)){
      bselect <- 1:dim(x1)[2]
      p1 <- c()
      p1        <- sapply( bselect, function(i) mytest2( x1[,i], y, R = R) )
      colnames(p1) <- colnames(x1)
      if(length(which(p1[2,] == min(p1[2,], na.rm = T))) > 1){
        bselect <- which.max(p1[1,])
      }else{
        bselect <- which.min(p1[2,])
      }
      
      
      splitindex <- s.opt(y=y, X=x1[, bselect], rnd.splt)
    }
  }
  
  # return split as partysplit object
  if(is.numeric(x)){
    return( partysplit( varid = as.integer(xselect),
                        breaks = splitindex,
                        info = list( p.value = 1 - (1 - p[xselect-1,2])^sum( !is.na(p[,2]) ) )
    ) )
  }
  if(is.factor(x)){
    return(partysplit(varid = as.integer(xselect),
                      index = splitindex,
                      info = list( p.value = 1 - (1 - p)^sum( !is.na(p) ) )
    ) )
  }
  
  
  if(is.list(x)){
    if(is.fdata(x$fdata.est)){
      return(list(sp = partysplit(varid = as.integer(bselect),
                                  breaks = splitindex,
                                  info=list( p.value = 1 - (1 - p[2,])^sum( !is.na(p[2,]) ) )
      ),varselect = xselect  ))
    }
  }
}

####

s.opt <- function(y, X, rnd = T){
  # find the split minimizing entropy
  s  <- sort(X)
  obj <- c()
  for(i in 1:length(s)){
    data1  <- y[ which(X <  s[i]) ]
    data2  <- y[ which(X >= s[i]) ]
    freqs1 <- table(data1)/length(data1)
    freqs2 <- table(data2)/length(data2)
    e1 <- entropy.empirical( freqs1, unit = "log2")
    e2 <- entropy.empirical( freqs2, unit = "log2")
    obj[i] = ( length(data1)*e1 + length(data2)*e2 )/length(y)
  }
  #if (!rnd){
  #  return( s[ which.min(obj) ] )
  #} else {
  #  return( s[ resample( which(obj == min(obj, na.rm = T)), 1) ] )
  #}
  return( s[ which.min(obj) ] )
}


