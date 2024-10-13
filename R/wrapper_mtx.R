learning_rate <- function(t, delay=1, forgetting=0.9){
  (t+delay)^(-forgetting)
}

svb_nmf_pois_mtx <- function(file_path, rank,
                        b_size = 100,
                        subiter = 10,
                        n_epochs = 10,
                        delay=1, forgetting=0.9,
                        prior_shape=1, prior_rate=1){
  #get matrix size
  con <- file(file_path, open = "r") #Open for reading in text mode
  rowsize <- scan(con, what=integer(), comment.char = "%", nmax=1,  quiet = TRUE)
  colsize <- scan(con, what=integer(), nmax=1, quiet = TRUE)
  N1 <- scan(con, what=integer(), nmax=1, quiet = TRUE)
  close(con)
  
  a <- prior_shape
  b <- prior_rate
  alpha_z = matrix(a, rowsize, rank)
  beta_z = matrix(b, 1, rank)
  alpha_w = matrix(a, colsize, rank)
  beta_w = matrix(b, 1, rank)
  rind <-sample.int(N1)
  sb <- floor(N1/b_size)
  m_sb <-b_size*sb
  subind <- vector("list",sb+1)
  for(i in 1:sb){
    subind[[i]] <- rind[1:b_size+b_size*(i-1)]  
  }
  subind[[sb+1]] <- rind[(b_size*sb):length(rind)]
  lp <- numeric(n_epochs)
  pb <- txtProgressBar(0, n_epochs, style = 3)
  for(ep in 1:n_epochs){
    rho <- learning_rate(ep, delay=delay, forgetting=forgetting)
    for(k in 1:length(subind)){
      Y <- dataloader_mtx(file_path, subind[[k]])
      #print(Y)
      out_t <- doVB_pois_s(Y[,3], Y[,1]-1L, Y[,2]-1L,
                            L=rank, iter = subiter,
                                    a,b,
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
    setTxtProgressBar(pb,ep)
  }
  out <- list(shape_row = alpha_z, rate_row=beta_z,
              shape_col = alpha_w, rate_col=beta_w,
              logprob=lp)
  class(out) <- "nmf_pois_posterior"
  return(out)
}
