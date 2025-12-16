size_mtx <- function(file_path){
  con <- file(file_path, open = "r") #Open for reading in text mode
  #get matrix size
  rowsize <- scan(con, what=integer(), comment.char = "%", nmax=1, quiet = TRUE)
  colsize <- scan(con, what=integer(), nmax=1, quiet = TRUE)
  len <- scan(con, what=integer(), nmax=1, quiet = TRUE)
  close(con)
  c(row=rowsize, column=colsize, nonzero=len)
}

SVBNMF_mtx <- function(file_path, rank,
                       n_epochs,
                       b_size,
                       lr_param,
                       lr_type="exponential",
                       subiter = 1L,
                       prior_shape=1, prior_rate=1,
                       display_progress=TRUE){
  size <- size_mtx(file_path)
  out <- doVB_pois_s_mtx(file_path, L=rank, iter = n_epochs, 
                         subiter = subiter,
                         a =  prior_shape, b =  prior_rate, N1 = size[3],
                         Nr = size[1], Nc = size[2],
                         ns = b_size, 
                         lr_param,
                         lr_type,
                         display_progress=display_progress)
  #class(out) <- "nmf_pois_posterior"
  return(out)
}

