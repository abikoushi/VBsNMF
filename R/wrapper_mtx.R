rand_pick <- function(N1, b_size){
  rind <-sample.int(N1)
  if(N1%%b_size==0){
    sb <- floor(N1/b_size)
    subind <- vector("list",sb)
    for(i in 1:sb){
      subind[[i]] <- rind[1:b_size+b_size*(i-1)]  
    }
  }else{
    sb <- floor(N1/b_size) 
    subind <- vector("list",sb+1)
    for(i in 1:sb){
      subind[[i]] <- rind[1:b_size+b_size*(i-1)]  
    }
    subind[[sb+1]] <- rind[(b_size*sb):length(rind)] 
  }
  return(subind)
}


SVBNMF <- function(file_path, rank,
                   subiter = 5,
                   n_epochs = 100,
                   b_size = 10000,
                   delay=1, forgetting=0.9,
                   prior_shape=1, prior_rate=1,
                   shuffle = TRUE,
                   learning_rate = lr_default,
                   dataloader = dataloader_mtx){
  #get matrix size
  dims = size_mtx(file_path) #get matrix size
  rowsize <- dims[1]
  colsize <- dims[2]
  N1 <- dims[3]
  b_size <- min(b_size, N1)
  a <- prior_shape
  b <- prior_rate
  alpha_z = matrix(a, rowsize, rank)
  beta_z = matrix(b, 1, rank)
  alpha_w = matrix(a, colsize, rank)
  beta_w = matrix(b, 1, rank)
  subind <- rand_pick(N1,b_size)
  lp <- numeric(n_epochs)
  pb <- txtProgressBar(0, n_epochs, style = 3)
  for(ep in 1:n_epochs){
    rho <- learning_rate(ep, delay=delay, forgetting=forgetting)
    for(k in 1:length(subind)){
      Y <- dataloader_mtx(file_path, subind[[k]])
      out_t <- doVB_pois_s(Y[,3], Y[,1]-1L, Y[,2]-1L,
                           L=rank, iter = subiter,
                           a, b,
                           N1,
                           alpha_z, beta_z,
                           alpha_w, beta_w)
      Ns <- length(subind[[k]])
      rho2 <- rho*(Ns/N1)
      alpha_z <- (1-rho2)*alpha_z + rho2*out_t$shape_row
      beta_z <- (1-rho2)*beta_z + rho2*out_t$rate_row
      alpha_w <- (1-rho2)*alpha_w + rho2*out_t$shape_col
      beta_w <- (1-rho2)*beta_w + rho2*out_t$rate_col
      lp[ep] <- lp[ep] + out_t$logprob[subiter]
    }
    if(shuffle){
      subind <- rand_pick(N1,b_size)
    }
    setTxtProgressBar(pb,ep)
  }
  out <- list(shape_row = alpha_z, rate_row=beta_z,
              shape_col = alpha_w, rate_col=beta_w,
              logprob=lp, family="poisson")
  return(out)
}
