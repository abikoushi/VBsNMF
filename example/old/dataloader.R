lr_default <- function(t, delay=1, forgetting=0.9){
  (t+delay)^(-forgetting)
}

scan1_mtx <- function(con, skip = 0){
  base::scan(con, nmax = 1, quiet=TRUE, 
             what = list(i=integer(), j=integer(), v=numeric()),
             skip = skip, comment.char = "%")
}

scan1_csv <- function(con, skip = 0, nlines=1){
  base::scan(con, quiet=TRUE, 
             what = numeric(),
             nlines = nlines,
             skip = skip, comment.char = "%", sep = ",")
}

dataloader_mtx <- function(file_path, bag){
  dims = size_mtx(file_path) #get matrix size
  rowsize <- dims[1]
  colsize <- dims[2]
  len <- dims[3]
  bag <- sort(bag)
  con <- file(file_path, open = "r") #Open for reading in text mode
  newL <- scan1_mtx(con, skip = bag[1]+1L) #initialize
  
  rowi = integer(length(bag))
  coli = integer(length(bag))
  yvec = numeric(length(bag))
  # out <- matrix(0, length(bag), 3)
  #out[1,] <- unlist(newL)
  rowi = newL$i
  coli = newL$j
  yvec = newL$v
  bag <- diff(bag)
  for(i in 1:length(bag)){
    newL <- scan1_mtx(con, skip = bag[i] - 1L)
    #out[i+1,] <- unlist(newL)
    rowi[i] = newL$i
    coli[i] = newL$j
    yvec[i] = newL$v
  }
  close(con)
  return(list(rowi, coli, yvec))
}


