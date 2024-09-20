
Y <- as(dat$Y, "dgTMatrix")
a <- 1
b <- 1
rank <- 2
alpha_z = matrix(a, Y@Dim[1], rank)
beta_z = matrix(b, 1, rank)
alpha_w = matrix(a, Y@Dim[2], rank)
beta_w = matrix(b, 1, rank)
for(k in 1:length(Y@x)){
  out_t <- VBsNMF:::doVB_pois(Y@x[k],Y@i[k],Y@j[k],
                              100,100,
                              2,1,
                              1,1,
                              alpha_z, beta_z,
                              alpha_w, beta_w)
  alpha_z[Y@i[k],] = out_t$shape_row[Y@i[k],]
  beta_z = out_t$rate_row
  alpha_w[Y@j[k],] = out_t$shape_col[Y@j[k],]
  beta_w = out_t$rate_col
}

Zhat <- sweep(alpha_z,1,c(beta_z),"/")
What <- sweep(alpha_w,1,c(beta_w),"/")
plot(Zhat%*%t(What), dat$Y, col=rgb(0,0,0,0.1))

