shape.data <- function(t.step = 128, n.samp = c(10, 10, 10),
                       eta.par = c(0,1), eps.par = c(0,1),
                       my.seed = NULL,
                       scrumble = FALSE,
                       plotta = TRUE){
  # Setup   \
  x   <-  1:t.step
  n.c <- n.samp[1]
  n.b <- n.samp[2]
  n.f <- n.samp[3]
  if (!is.null(my.seed)) set.seed(my.seed)
  
  # Sample from cylinder
  a <- sample(16:32, n.c, replace = TRUE)
  a <- matrix(a, nrow = n.c, ncol = t.step)
  delta <- sample(32:96, n.c, replace = TRUE)
  delta <- matrix(delta, nrow = n.c, ncol = t.step)  
  b     <- a + delta
  eta   <- rnorm(n.c, mean = eta.par[1], sd = eta.par[2])
  eta   <- matrix(eta, nrow = n.c, ncol = t.step)
  eps   <- matrix(rnorm(t.step*n.c, mean = eps.par[1], sd = eps.par[2]), 
                  nrow = n.c, ncol = t.step)
  # Build
  cyl <- matrix(0,nrow=n.c,ncol=t.step)
  for(i in 1:n.c){
    cyl[i,] <- (6 + eta[i,])*(x >= a[i,])*(x <= b[i,]) + eps[i,]
  }
  
  # Sample from bell
  a <- sample(16:32, n.b, replace = TRUE)
  a <- matrix(a, nrow = n.b, ncol = t.step)
  delta <- sample(32:96, n.b, replace = TRUE)
  delta <- matrix(delta, nrow = n.b, ncol = t.step)  
  b     <- a + delta
  eta   <- rnorm(n.b, mean = eta.par[1], sd = eta.par[2])
  eta   <- matrix(eta, nrow = n.b, ncol = t.step)
  eps   <- matrix(rnorm(t.step*n.b, mean = eps.par[1], sd = eps.par[2]), 
                  nrow = n.b, ncol = t.step)
  # Build
  bell <- matrix(0, nrow =n.b, ncol=t.step)
  for(i in 1:n.b){
    bell[i,] <- (6 + eta[i,])*(x >= a[i,])*(x <= b[i,])*(x - a[i,])/delta[i,] + eps[i,]
  }
  
  # Sample from funnel
  a <- sample(16:32, n.f, replace = TRUE)
  a <- matrix(a, nrow = n.f, ncol = t.step)
  delta <- sample(32:96, n.f, replace = TRUE)
  delta <- matrix(delta, nrow = n.f, ncol = t.step)  
  b     <- a + delta
  eta   <- rnorm(n.f, mean = eta.par[1], sd = eta.par[2])
  eta   <- matrix(eta, nrow = n.f, ncol = t.step)
  eps   <- matrix(rnorm(t.step*n.f, mean = eps.par[1], sd = eps.par[2]), 
                  nrow = n.f, ncol = t.step)
  # Build
  fun <- matrix(0, nrow =n.f, ncol=t.step)
  for(i in 1:n.f){
    fun[i,] <- (6 + eta[i,])*(x >= a[i,])*(x <= b[i,])*(b[i,] - x)/delta[i,] + eps[i,]
  }
  # Collect
  cls  <- c(rep("Cyl", n.c), rep("Bel", n.b), rep("Fun", n.f))
  out  <- rbind(cyl, bell, fun)
  cls  <- factor(cls)
  cls  <- data.frame(cls)
  out  <- cbind(cls, out)
  mean.cyl <- apply(cyl,2,mean)
  mean.bell <- apply(bell,2,mean)
  mean.fun <- apply(fun,2,mean)
  # Scrumble
  if (scrumble){
    idx <- sample(1:nrow(out))
    out <- out[idx,]
  }
  
  # Plot
  if (plotta){
    camp <- sample(1:t.step, 1)
    par(mfrow = c(2,2))
    matplot(t(cyl), lty = 1, type = "l", 
            main = "Cylinder", xlab = "", ylab = "",
            sub = paste("num. samples:", n.c),
            col="yellow")
    lines(cyl[camp,], lty=1)
    matplot(t(bell), lty = 1, type = "l", 
            main = "Bell", xlab = "", ylab = "",
            sub = paste("num. samples:", n.b),
            col="red")
    lines(bell[camp,], lty=1)
    matplot(t(fun), lty = 1, type = "l", 
            main = "Funnel", xlab = "", ylab = "",
            sub = paste("num. samples:", n.f),
            col="blue")
    lines(fun[camp,], lty=1)
    matplot(mean.cyl, lty = 1, type="l",
            main = "Means", xlab = "", ylab = "",
            sub = paste("num. samples:", n.f),col="yellow")
    lines(mean.bell,col="red")
    lines(mean.fun,col="blue")
  }
  # Output
  return(out)
}
