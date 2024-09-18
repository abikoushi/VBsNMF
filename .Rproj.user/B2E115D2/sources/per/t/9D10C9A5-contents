vb_nmf_pois <- function(Y, L, iter=100, 
         prior_shape=0.5, prior_rate=1e-5,
         rowinit=NULL, colinit=NULL){
  if(!any(class(Y)=="dgTMatrix")){
    Y <- as(dat$Y, "TsparseMatrix")    
  }
  if(is.null(rowinit)){
    rowinit = matrix(rexp(Y@Dim[1]*L), Y@Dim[1], L)    
  }
  if(is.null(colinit)){
    colinit = matrix(rexp(Y@Dim[2]*L), Y@Dim[2], L)
  }
  out <- doVB_pois(y = Y@x, rowi = Y@i, coli = Y@j,
                   N = Y@Dim[1], M = Y@Dim[2],
                   L = L, iter=iter,
                   a = prior_shape, b = prior_rate,
                   Z = rowinit, W = colinit)
  colnames(out$logprob) <- c("loglik1", "loglik2", "KLD")
  return(out)
}
