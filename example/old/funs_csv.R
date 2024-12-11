setinit_csv <- function(prior_shape, prior_rate, rowsize, colsize, rank){
  pdir <- "VBparams"
  if(dir.exists(pdir)){
    system(paste("rm -r", pdir))
  }
  dir.create(pdir)
  
  A <- rep(prior_shape, rank)
  for(i in 1:rowsize){
    cat(A, file = paste(pdir,"shape_row.csv", sep="/"), sep = ", ", append = TRUE)
    cat("\n", file = paste(pdir,"shape_row.csv", sep="/"), append = TRUE)    
  }
  
  for(i in 1:colsize){
    cat(A, file = paste(pdir,"shape_col.csv", sep="/"), sep = ", ", append = TRUE)
    cat("\n", file = paste(pdir,"shape_col.csv", sep="/"), append = TRUE)
  }
  
  Bz = matrix(prior_rate, 1, rank)
  
  cat(Bz, file = paste(pdir,"rate_row.csv", sep="/"), sep = ", ", append = TRUE)
  cat("\n", file = paste(pdir,"rate_row.csv", sep="/"), append = TRUE)
  
  Bw = matrix(prior_rate, 1, rank)  
  cat(Bw, file = paste(pdir,"rate_col.csv", sep="/"), sep = ", ", append = TRUE)
  cat("\n", file = paste(pdir,"rate_col.csv", sep="/"), append = TRUE)
  
  paramlist <-  list(alpha_z = NULL,
                     beta_z = Bz,
                     alpha_w = NULL,
                     beta_w = Bw,
                     L = rank)
  return(paramlist)
}

query_awk <- function(vpar, pdir, filename, uid){
  chl <- paste0("cat ", pdir, "/", filename)
  for(i in 1:nrow(vpar)){
    chr <- paste0("\"", paste(vpar[i,], collapse = ", "), "\"")
    awk <- paste0(" | awk 'NR==", uid[i], "{sub('/.*/',", chr, ")}1'")
    chl <- paste0(chl, awk)
  }
  chl <- paste0(chl, " > ", pdir, "/tmp.csv")
  mvit <- paste0("mv ", pdir, "/tmp.csv " , pdir, "/", filename)
  system(chl)
  system(mvit)
}

updatefun_csv <- function(paramlist, paramlist2, upparmalist, uid_r, uid_c, rho, rho2){
  paramlist2$alpha_z <- rho2*paramlist2$alpha_z + rho*upparmalist$shape_row
  query_awk(paramlist2$alpha_z, "VBparams", "shape_row.csv", uid_r)
  
  paramlist2$alpha_w <- rho2*paramlist2$alpha_w + rho*upparmalist$shape_col
  query_awk(paramlist2$alpha_w, "VBparams", "shape_col.csv", uid_c)
  
  paramlist2$beta_z <- rho2*paramlist2$beta_z + rho*upparmalist$rate_row
  query_awk(paramlist2$beta_z, "VBparams", "rate_row.csv", 1)
  
  paramlist2$beta_w <- rho2*paramlist2$beta_w + rho*upparmalist$rate_col
  query_awk(paramlist2$beta_w, "VBparams", "rate_col.csv", 1)
  return(paramlist2)
}


load_alpha_csv <- function(file_path, bag, L){
  out <- matrix(0, length(bag), L)
  bag <- sort(bag)
  con <- file(file_path, open = "r") #Open for reading in text mode
  newL <- scan1_csv(con, skip = bag[1]-1L)
  bag <- diff(bag)
  out[1,] <- newL
  for(i in 1:length(bag)){
    newL <- scan1_csv(con, skip = bag[i]-1L)
    out[i+1L,] <- newL
  }
  close(con)
  return(out)
}

draw_param_csv <- function(paramlist, uid_r, uid_c){
  path <- "VBparams/shape_row.csv"
  paramlist$alpha_z <- load_alpha_csv(path, bag=uid_r, paramlist$L)
  path <- "VBparams/shape_col.csv"
  paramlist$alpha_w <- load_alpha_csv(path, bag=uid_c,  paramlist$L)
  return(paramlist)
}


