vb_nmf_pois <- function(Y, rank, iter=100, 
         prior_shape=0.5, prior_rate=1e-5){
  if(!any(class(Y)=="dgTMatrix")){
    Y <- as(Y, "TsparseMatrix")    
  }
  alpha_z = matrix(prior_shape, Y@Dim[1], rank)
  beta_z = matrix(prior_rate, 1, rank)
  alpha_w = matrix(prior_shape, Y@Dim[2], rank)
  beta_w = matrix(prior_rate, 1, rank)
  out <- doVB_pois(y = Y@x, rowi = Y@i, coli = Y@j,
                   L = rank, iter=iter,
                   a = prior_shape, b = prior_rate,
                   alpha_z, beta_z, alpha_w, beta_w)
  #colnames(out$logprob) <- c("loglik1", "loglik2", "KLD")
  #out$elbo <- rowSums(out$logprob)
  class(out) <- "nmf_pois_posterior"
  return(out)
}

learning_rate <- function(t, delay=1, forgetting=0.9){
  (t+delay)^(-forgetting)
}

#write file loader

svb_nmf_pois<- function(Y, rank,
                        b_size = 100,
                        subiter = 10,
                        n_epochs = 10,
                        delay=1, forgetting=0.9,
                        prior_shape=0.5, prior_rate=0.1){
  if(!any(class(Y)=="dgTMatrix")){
    Y <- as(Y, "TsparseMatrix")    
  }
  a <- prior_shape
  b <- prior_rate
  alpha_z = matrix(a, Y@Dim[1], rank)
  beta_z = matrix(b, 1, rank)
  alpha_w = matrix(a, Y@Dim[2], rank)
  beta_w = matrix(b, 1, rank)
  
  N1 <- length(Y@x)
  rind <-sample.int(N1)

  sb <- floor(N1/b_size)
  m_sb <-b_size*sb
  subind <- vector("list",sb+1)
  for(i in 1:sb){
    subind[[i]] <- rind[1:b_size+b_size*(i-1)]  
  }
  subind[[sb+1]] <- rind[(b_size*sb):length(rind)]
  
  lp <- numeric(n_epochs)
  pb <- txtProgressBar(1, n_epochs, style = 3)
  for(ep in 1:n_epochs){
    rho <- learning_rate(ep, delay = delay, forgetting = forgetting)
    for(k in 1:length(subind)){
      out_t <- VBsNMF:::doVB_pois_s(Y@x[subind[[k]]],Y@i[subind[[k]]],Y@j[subind[[k]]],
                                    L=2, iter = subiter,
                                    a,b,
                                    N1,
                                    alpha_z, beta_z,
                                    alpha_w, beta_w)
      Ns <- length(subind[[k]])
      alpha_z <- (1-rho)*alpha_z + rho*out_t$shape_row
      beta_z <- (1-rho)*beta_z + rho*out_t$rate_row
      alpha_w <- (1-rho)*alpha_w + rho*out_t$shape_col
      beta_w <- (1-rho)*beta_w + rho*out_t$rate_col
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

basemean <- function(obj){
  stopifnot(any(class(obj)=="nmf_pois_posterior"))
  with(obj, sweep(shape_row,1,c(rate_row),"/"))
} 

coefmean <- function(obj){
  stopifnot(any(class(obj)=="nmf_pois_posterior"))
  t(with(obj, sweep(shape_col, 1, c(rate_col),"/")))
}
