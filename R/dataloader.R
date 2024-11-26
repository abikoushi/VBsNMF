lr_default <- function(t, delay=1, forgetting=0.9){
  (t+delay)^(-forgetting)
}

scan1_mtx <- function(con, skip = 0){
  base::scan(con, nmax = 1, quiet=TRUE, 
             what = list(i=integer(), j=integer(), v=numeric()),
             skip = skip, comment.char = "%")
}

size_mtx <- function(file_path){
  con <- file(file_path, open = "r") #Open for reading in text mode
  #get matrix size
  rowsize <- scan(con, what=integer(), comment.char = "%", nmax=1,  quiet = TRUE)
  colsize <- scan(con, what=integer(), nmax=1, quiet = TRUE)
  len <- scan(con, what=integer(), nmax=1, quiet = TRUE)
  close(con)
  c(row=rowsize, column=colsize, nonzero=len)
}

dataloader_mtx <- function(file_path, bag){
  dims = size_mtx(file_path) #get matrix size
  rowsize <- dims[1]
  colsize <- dims[2]
  len <- dims[3]
  bag <- sort(bag)
  con <- file(file_path, open = "r") #Open for reading in text mode
  newL <- scan1_mtx(con, skip = bag[1]+1L) #initialize
  out <- matrix(0, length(bag), 3)
  out[1,] <- unlist(newL)
  bag <- diff(bag)
  for(i in 1:length(bag)){
    newL <- scan1_mtx(con, skip = bag[i] - 1L)
    out[i+1,] <- unlist(newL)
  }
  close(con)
  return(out)
}
