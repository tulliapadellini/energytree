gen_data_reg  <- function( size, P, n.var, alpha, beta,
                      a, b, sigma=0.1){
  
  
  if(dim(alpha)[2]!=n.var | dim(beta)[2]!=n.var) stop("Dimensions Error")
  if(dim(alpha)[1]!=length(size) | dim(beta)[1]!=length(size)) stop("Dimensions Error")
  if(dim(a)[2]!=n.var | dim(b)[2]!=n.var) stop("Dimensions Error")
  if(dim(a)[1]!=length(size) | dim(b)[1]!=length(size)) stop("Dimensions Error")
  
  nclass <- length(class)

  data <- list()
  for( i in 1:n.var){
    
    
    grid <- seq(0, 1, length.out = P)
    Cov <- exp_cov_function(grid, alpha = alpha[1,i], beta = beta[1,i])
    Data <- generate_gauss_fdata(size[1], 
                                 centerline = cos(a[1,i]+b[1,i]*pi*grid), Cov = Cov)
    fD1 <- fData( grid, Data )
    
    Cov <- exp_cov_function(grid, alpha=alpha[2,i], beta=beta[2,i])
    Data <- generate_gauss_fdata(size[2], 
                                 centerline = sin(a[2,i]+b[2,i]*pi*grid), Cov = Cov)
    fD2 <- fData( grid, Data )
    
    Cov <- exp_cov_function(grid, alpha=alpha[3,i], beta=beta[3,i])
    Data <- generate_gauss_fdata(size[3], 
                                 centerline = cos(a[3,i]+b[3,i]*pi*grid)*sin(a[3,i]+b[3,i]*pi*grid), Cov = Cov)
    fD3 <- fData( grid, Data )
    
    eps1   <- matrix(rnorm(P*size[1], mean = 0, sd = sigma), 
                     nrow = size[1], ncol = P)
    eps2   <- matrix(rnorm(P*size[2], mean = 0, sd = sigma), 
                     nrow = size[2], ncol = P)
    eps3   <- matrix(rnorm(P*size[3], mean = 0, sd = sigma), 
                     nrow = size[3], ncol = P)
    
    k1 <- fD1$values + eps1
    k2 <- fD2$values + eps2
    k3 <- fD3$values + eps3
    k <- rbind(k1,k2,k3)
    
    
    
    data[[paste("V",i,sep="")]] <- fdata(k)
    
    
    data[["Y"]] <- rowMeans(k)
  }
  
  return(data)
  
  
}

