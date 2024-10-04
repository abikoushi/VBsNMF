Y <- as(log(as.matrix(iris[,-5])), "TsparseMatrix")

out <- VBsNMF:::doVB_norm(Y@x, rowi=Y@i, coli = Y@j,
                          Nr = Y@Dim[1], Nc=Y@Dim[2], 
                          L=2,iter=5,prior_prec=1,
                          a=1,b=1)
plot(out$lp)
plot(as.matrix(iris[,-5]),out$mean_z%*%t(out$mean_w))
