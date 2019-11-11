gen_data_reg  <- function( size, P, n.var, alpha, beta,  sigma=0.1){
  
  
  nclass <- length(class)
  
  data <- list()
  for( i in 1:n.var){
    
    
    grid <- seq(0, 1, length.out = P)
    Cov <- exp_cov_function(grid, alpha = alpha[1,i], beta = beta[1,i])
    Data <- generate_gauss_fdata(size[1], 
                                 centerline = ((2-5*grid)/2)*((5*grid-2)/2)^2+ sin(5*pi*grid/2), Cov = Cov)
    fD1 <- fData( grid, Data )
    
    Cov <- exp_cov_function(grid, alpha=alpha[2,i], beta=beta[2,i])
    Data <- generate_gauss_fdata(size[2], 
                                 centerline = -((2-5*grid)/2)*((5*grid-2)/2)^2+ sin(5*pi*grid/2), Cov = Cov)
    fD2 <- fData( grid, Data )
    
    Cov <- exp_cov_function(grid, alpha=alpha[3,i], beta=beta[3,i])
    Data <- generate_gauss_fdata(size[3], 
                                 centerline = cos(2*pi*grid), Cov = Cov)
    fD3 <- fData( grid, Data )
    
    Cov <- exp_cov_function(grid, alpha=alpha[4,i], beta=beta[4,i])
    Data <- generate_gauss_fdata(size[4], 
                                 centerline = -cos(2*pi*grid), Cov = Cov)
    fD4 <- fData( grid, Data )
    
    eps1   <- matrix(rnorm(P*size[1], mean = 0, sd = sigma), 
                     nrow = size[1], ncol = P)
    eps2   <- matrix(rnorm(P*size[2], mean = 0, sd = sigma), 
                     nrow = size[2], ncol = P)
    eps3   <- matrix(rnorm(P*size[3], mean = 0, sd = sigma), 
                     nrow = size[3], ncol = P)
    eps4   <- matrix(rnorm(P*size[4], mean = 0, sd = sigma), 
                     nrow = size[4], ncol = P)
    
    k1 <- fD1$values + eps1
    k2 <- fD2$values + eps2
    k3 <- fD3$values + eps3
    k4 <- fD4$values + eps4
    k <- rbind(k1,k2,k3,k4)
    
    
    
    data[[paste("V",i,sep="")]] <- fdata(k)
    
    
    data[["Y"]] <- c(rnorm(50,-10,3),rnorm(50,10,3),rnorm(50,15,5),rnorm(50,-15,5))
  }
  
  return(data)
  
  
}
