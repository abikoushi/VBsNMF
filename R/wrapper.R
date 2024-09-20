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
                   N = Y@Dim[1], M = Y@Dim[2],
                   L = rank, iter=iter,
                   a = prior_shape, b = prior_rate,
                   alpha_z, beta_z, alpha_w, beta_w)
  colnames(out$logprob) <- c("loglik1", "loglik2", "KLD")
  return(out)
}
