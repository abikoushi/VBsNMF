scan1_mtx <- function(con, skip = 0){
  base::scan(con, nmax = 1, quiet=TRUE, 
             what = list(i=integer(), j=integer(), v=numeric()),
             skip = skip)
}

dataloader_mtx <- function(file_path, bag){
  con <- file(file_path, open = "r") #Open for reading in text mode
  #get matrix size
  rowsize <- scan(con, what=integer(), comment.char = "%", nmax=1,  quiet = TRUE)
  colsize <- scan(con, what=integer(), nmax=1, quiet = TRUE)
  len <- scan(con, what=integer(), nmax=1, quiet = TRUE)
  bag <- sort(bag)
  newL <- scan1_mtx(con, skip = bag[1] - 1L) #initialize
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

# dataloader_mtx <- function(file_path, bag){
#   con <- file(file_path, open = "r") #Open for reading in text mode
#   #get matrix size
#   rowsize <- scan(con, what=integer(), comment.char = "%", nmax=1,  quiet = TRUE)
#   colsize <- scan(con, what=integer(), nmax=1, quiet = TRUE)
#   len <- scan(con, what=integer(), nmax=1, quiet = TRUE)
#   newL <- scan1_mtx(con) #initialize
#   out <- matrix(0, length(bag), 3)
#   k <- 1
#   bag <- sort(bag)
#   for(i in 1:len){
#     if(bag[k]==i){
#       out[k,] <- unlist(newL)
#       k <- k+1
#       if(k>length(bag)){
#         break
#       }
#     }
#     newL <-scan1_mtx(con) #next one
#   }
#   close(con)
#   return(out)
# }
# 
