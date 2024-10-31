vb_nmf_pois <- function(Y, rank, iter=100,
                        prior_shape=1, prior_rate=1){
  if(!any(class(Y)=="dgTMatrix")){
    Y <- as(Y, "TsparseMatrix")
  }
  wch <- which(is.na(Y@x))
  if(length(wch)>0){
    ni <- rep(Y@Dim[2],Y@Dim[1])
    nj <- rep(Y@Dim[1],Y@Dim[2])
    for(k in 1:length(wch)){
      ni[Y@i[wch[k]]] <- ni[Y@i[wch[k]]]-1
      nj[Y@j[wch[k]]] <- nj[Y@j[wch[k]]]-1  
    }
    wrow <- ni/Y@Dim[1]
    wcol <- nj/Y@Dim[2]
    out <- doVB_pois_na(y = Y@x[-wch], rowi = Y@i[-wch], coli = Y@j[-wch],
                        Y@Dim[1], Y@Dim[2],
                        L = rank, iter=iter, 
                        a = prior_shape, b = prior_rate,
                        wrow = wrow, wcol = wcol)
                        #alpha_z, beta_z, alpha_w, beta_w
  }else{
    out <- doVB_pois(y = Y@x, rowi = Y@i, coli = Y@j,
                     Y@Dim[1], Y@Dim[2],
                     L = rank, iter=iter,
                     a = prior_shape, b = prior_rate)
                     #alpha_z, beta_z, alpha_w, beta_w  
  }
  #colnames(out$logprob) <- c("loglik1", "loglik2", "KLD")
  #out$elbo <- rowSums(out$logprob)
  class(out) <- "nmf_pois_posterior"
  return(out)
}

basemean <- function(obj){
  stopifnot(any(class(obj)=="nmf_pois_posterior"))
  with(obj, sweep(shape_row, 2, c(rate_row),"/"))
} 

coefmean <- function(obj){
  stopifnot(any(class(obj)=="nmf_pois_posterior"))
  t(with(obj, sweep(shape_col, 2, c(rate_col),"/")))
}

###

em_nmf_pois <- function(Y, rank, iter=100, prior_shape=1, prior_rate=1){
  if(!any(class(Y)=="dgTMatrix")){
    Y <- as(Y, "TsparseMatrix")
  }
  out <- doEM_pois(y = Y@x, rowi = Y@i, coli = Y@j,
                     Y@Dim[1], Y@Dim[2],
                     L = rank, iter=iter,
                     a = prior_shape, b = prior_rate)

  return(out)
}
