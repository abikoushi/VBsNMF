scan1_mtx <- function(con){
  base::scan(con, nmax = 1, quiet=TRUE, what=list(i=integer(), j=integer(), v=numeric()))
}

dataloader_mtx <- function(file_path, bag){
  con <- file(file_path, open = "r") #Open for reading in text mode
  #get matrix size
  rowsize <- scan(con, what=integer(), comment.char = "%", nmax=1,  quiet = TRUE)
  colsize <- scan(con, what=integer(), nmax=1, quiet = TRUE)
  len <- scan(con, what=integer(), nmax=1, quiet = TRUE)
  newL <- scan1_mtx(con) #initialize
  out <- matrix(0, length(bag), 3)
  k <- 1
  bag <- sort(bag)
  for(i in 1:len){
    if(bag[k]==i){
      out[k,] <- unlist(newL)
      k <- k+1
      if(k>length(bag)){
        break
      }
    }
    newL <-scan1_mtx(con) #next one
  }
  close(con)
  return(out)
}
