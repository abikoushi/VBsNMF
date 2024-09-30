library(Matrix)
library(VBsNMF)
library(NMF)

dens <- seq(0.1, 0.9, by=0.1)
Y <- rsparsematrix(100, 100, 0.9,
                   rand.x= function(n){rpois(n,3)+1},
                   repr = "T")
Y2 <- as.matrix(Y)

meth <- nmfAlgorithm()

t1 <- system.time({out = vb_nmf_pois(Y, rank=3, iter=100)})
t2 <- system.time({out2 = nmf(Y2, rank=3, maxIter=100, eps=0,
                              method = "lee")})
sqrt(mean((Y2-basis(out2)%*%coef(out2))^2))
sqrt(mean((Y2-basemean(out)%*%coefmean(out))^2))
mean(dpois(Y2,basis(out2)%*%coef(out2), log = TRUE))
mean(dpois(Y2,basemean(out)%*%coefmean(out), log = TRUE))

