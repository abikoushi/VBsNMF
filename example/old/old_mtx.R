rand_pick <- function(N1, b_size){
  rind <-sample.int(N1)
  rem <- N1%%b_size
  if(rem==0){
    sb <- floor(N1/b_size)
    subind <- vector("list",sb)
    for(i in 1:sb){
      subind[[i]] <- rind[1:b_size+b_size*(i-1)]  
    }
  }else{
    sb <- N1%/%b_size
    subind <- vector("list", sb+1)
    for(i in 1:sb){
      subind[[i]] <- rind[1:b_size+b_size*(i-1)]  
    }
    r_rem  <-  rind[1:rem]
    subind[[sb+1]] <- c(rind[(b_size*sb):length(rind)], r_rem)
  }
  return(subind)
}

fac_int <- function(x, uid){
  as.integer(factor(x, levels = uid))-1L
}

setinit_default <- function(prior_shape, prior_rate, rowsize, colsize, rank){
  paramlist <- list(alpha_z = matrix(prior_shape, rowsize, rank),
                    beta_z = matrix(prior_rate, 1, rank),
                    alpha_w = matrix(prior_shape, colsize, rank),
                    beta_w = matrix(prior_rate, 1, rank),
                    L = rank)
  return(paramlist)
}



draw_param_default <- function(paramlist, uid_r, uid_c){
  paramlist2 <- paramlist
  paramlist2$alpha_z <- paramlist$alpha_z[uid_r,]
  paramlist2$alpha_w <- paramlist$alpha_w[uid_c,]
  return(paramlist2)
}


updatefun_default <- function(paramlist, paramlist2, upparmalist, uid_r, uid_c, rho, rho2){
  paramlist$alpha_z[uid_r,] <- rho2*paramlist$alpha_z[uid_r,] + rho*upparmalist$shape_row
  paramlist$beta_z <- rho2*paramlist$beta_z + rho*upparmalist$rate_row
  paramlist$alpha_w[uid_c,] <- rho2*paramlist$alpha_w[uid_c, ] + rho*upparmalist$shape_col
  paramlist$beta_w <- rho2*paramlist$beta_w + rho*upparmalist$rate_col
  return(paramlist)
}

# SVBNMF <- function(file_path, rank,
#                    subiter = 5,
#                    n_epochs = 100,
#                    b_size = 10000,
#                    delay=1, forgetting=0.9,
#                    prior_shape=1, prior_rate=1,
#                    shuffle = TRUE,
#                    learning_rate = lr_default,
#                    dataloader = dataloader_mtx,
#                    setinit = setinit_default,
#                    updatefun = updatefun_default,
#                    draw_param = draw_param_default){
#   #get matrix size
#   dims = size_mtx(file_path) #get matrix size
#   rowsize <- dims[1]
#   colsize <- dims[2]
#   N1 <- dims[3]
#   b_size <- min(b_size, N1)
#   paramlist <- setinit(prior_shape, prior_rate, rowsize, colsize, rank)
#   subind <- rand_pick(N1, b_size)
#   lp <- numeric(n_epochs)
#   pb <- txtProgressBar(0, n_epochs, style = 3)
#   for(ep in 1:n_epochs){
#     rho <- learning_rate(ep, delay=delay, forgetting=forgetting)
#     for(k in 1:length(subind)){
#       uid_r = unique(Y[,1])
#       uid_c = unique(Y[,2])
#       paramlist2 <- draw_param(paramlist, uid_r, uid_c)
#       out_t <- doVB_pois_s(y = Y[,3], 
#                            rowi = fac_int(Y[,1], uid_r),
#                            coli = fac_int(Y[,2], uid_c),
#                            L = rank, iter = subiter,
#                            a = prior_shape, b = prior_rate,
#                            N1 = N1,
#                            paramlist2$alpha_z, paramlist2$beta_z,
#                            paramlist2$alpha_w, paramlist2$beta_w)
#       Ns <- length(subind[[k]])
#       rho2 <- 1-rho
#       paramlist <- updatefun(paramlist, paramlist2, out_t, uid_r, uid_c, rho, rho2)
#       lp[ep] <- lp[ep] + out_t$logprob[subiter]
#     }
#     if(shuffle){
#       subind <- rand_pick(N1, b_size)
#     }
#     setTxtProgressBar(pb,ep)
#   }
#   out <- list(shape_row = paramlist$alpha_z,
#               rate_row=paramlist$beta_z,
#               shape_col = paramlist$alpha_w,
#               rate_col=paramlist$beta_w,
#               logprob=lp, family="poisson")
#   return(out)
# }
# 
