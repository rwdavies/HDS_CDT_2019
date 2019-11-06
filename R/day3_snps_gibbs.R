setwd("/Users/robert davies/Google Drive/Presentations and conferences/2019_11_03 CDT Lectures/Day3 - Gibbs sampling/")
source("/Users/robert davies/proj/HDS_CDT_2019/R/day3_snp_gibbs_functions.R")
source("/Users/robert davies/proj/HDS_CDT_2019/R/day3_stephens.R")
main <- ""
library("xtable")



##
## simple simulate
##
T <- 6
K <- 2
theta <- t(rbind(
    c(0.5,0.5,0.5,0.5,0.5,0.5),
    c(0.001,0.999,0.001,0.999,0.001,0.999)
))
N <- 50
true_Z <- sample(c(1, 2), N, replace = TRUE)
H <- array(NA, c(T, N))
for(i in 1:N) {
    H[, i] <- as.integer(runif(T) < theta[, true_Z[i]])
}
## sample H
out <- run_gibbs(H, K = K, nits = 50)
image(H)
Zstore <- out$Zstore
lls <- out$lls
plot_groups(Zstore, lls, true_Z, pngname = "eleph1.png")

## keep going, do 1000
out <- run_gibbs(H, K = K, nits = 1000)
## averaged state
w <- seq(10, 1000, 20)
estimated_Z <- colMeans(out$Zstore[w, ])
post_Z_is_1 <- colSums(out$Zstore[w, ] == 1) / length(w)
xtable(table(true_Z, post_Z_is_1))

1



##
## still simple, more
##
T <- 100
K <- 2
theta <- array(runif(T * K), c(T, K))
N <- 50
true_Z <- sample(c(1, 2), N, replace = TRUE)
H <- array(NA, c(T, N))
for(i in 1:N) {
    H[, i] <- as.integer(runif(T) < theta[, true_Z[i]])
}
## sample H
out <- run_gibbs(H, K = K, nits = 50)
image(H)
Zstore <- out$Zstore
lls <- out$lls
plot_groups(Zstore, lls, true_Z, pngname = "eleph2.png")

out <- run_gibbs(H, K = K, nits = 1000)

estimated_Z <- colMeans(out$Zstore[seq(1, 1000, 20), ])
xtable(table(true_Z, estimated_Z))
post_Z_is_1 <- colSums(out$Zstore[w, ] == 1) / length(w)
xtable(table(true_Z, post_Z_is_1))


1
## 











##
## humans
##
## sample 100 SNPs per chromosome, one per hundred kbp
hAll <- get_data_per_chr(20) ## h_{t, i}
hAll[is.na(hAll)] <- 0 ## arrrrgh




## can look at correlation here, pretty independent
round(cor(t(h[1:20, ])) ** 2, 2)


##
## choose easily discernable subsets, 50 of each
##
H <- hAll[, c(
    grep("CEU", colnames(hAll))[1:50],
    grep("JPT", colnames(hAll))[1:50],
    grep("YRI", colnames(hAll))[1:50]
)]
true_z <- c(rep(1, 50), rep(2, 50), rep(3, 50))

for(f in 1:2) {
    for(i in 1:5) {
        for(K in c(3, 6)) {
            set.seed(91 * i * (f - 1))
            if (f == 1) { Hl <- H }
            if (f == 2) { Hl <- H[1:10, ]}
            out <- run_gibbs(H, K = K)
            Zstore <- out$Zstore
            lls <- out$lls
            plot_groups(Zstore, lls, true_z, paste0("easy", K, ".group1", i, ".", f, ".png"))
        }
    }
}

## kind of sucks! do a better init, using hacky method
Z_init <- tiny_pca(H, K)
out <- run_gibbs(H, Z_init = Z_init, K = 3)
Zstore <- out$Zstore
lls <- out$lls
plot_groups(Zstore, lls, true_z, paste0("easy3.group1.pca.png"))






##
## heck, try a bunch of different pops
##
cols <- sapply(c("GBR", "ASW", "YRI"), function(code) {
    grep(code, colnames(hAll))[1:50]
})

H <- hAll[seq(1, 562, 5), c(cols)]
true_z <- c(rep(1, 50), rep(2, 50), rep(3, 50))
out <- run_gibbs(H, K = 2, nits = 20)
Zstore <- out$Zstore
lls <- out$lls
plot_groups(Zstore, lls, true_z, pngname = "asw.png")












##
## closer example
##
H <- hAll[, c(
    grep("CEU", colnames(hAll))[1:100],
    grep("TSI", colnames(hAll))[1:100]
)]
K <- 2
true_z <- c(rep(1, 100), rep(2, 100))
Z_init <- tiny_pca(H, K)
## Z_init <- true_z
out <- run_gibbs(H, Z_init = Z_init, K = K, nits = 50)
Zstore <- out$Zstore
lls <- out$lls
res <- gibbs(t(H), niter = 50) ## matthew stephens version

png("weird.png", res = 150, units = "in", height = 4, width = 8)
par(mfrow = c(1, 2))
plot_groups(Zstore, lls, true_z)
plot_groups(res$z, array(1, length(lls)), true_z)
dev.off()

## OK, so indeed, doesn't seem to work as well as expected!
## weird property, when difficult to distinguish, collapses!
## write up where I am now
1 



##
## closer example
##
H <- hAll[, c(
    grep("YRI", colnames(hAll))[1:50]
)
]
Z_init <- tiny_pca(H, K = 3)
true_z <- rep(1, 50)
## Z_init <- true_z
out <- run_gibbs(H, Z_init = Z_init, K = 3, nits = 50)
Zstore <- out$Zstore
lls <- out$lls
plot_groups(Zstore, lls, true_z) #, paste0("easy3.group2.pca.png"))





##
## scratch
##
##
## simple example 
##
out <- run_gibbs(H, Z_init = Z_init, K = K, nits = 50)
Zstore <- out$Zstore
lls <- out$lls
res <- gibbs(t(H), niter = 50) ## matthew stephens version




