source("/Users/robert davies/proj/HDS_CDT_2019/R/day3_phasing_functions.R")
setwd("/Users/robert davies/Google Drive/Presentations and conferences/2019_11_03 CDT Lectures/Day3 - Gibbs sampling/")
library("xtable")

## 
## "prove" the two versions make sense
## by comparing exhaustive augmented version against proposal
##
nReads <- 5
nSNPs <- 5
seed <- 10
eps <- 0.1
out <- sim_reads(nReads, nSNPs, seed, eps)
s <- out$s
u <- out$u

y <- t(sapply(0:(2 ** nReads - 1), function(iRead) {
    z_proposed <- as.integer(intToBits(iRead)[1:nReads])
    ## p_z_given_o(z_proposed)
    c(p_z_given_o_augmented(z_proposed, s, u), p_z_given_o(z_proposed, s, u))
}))

## check here, stop if they are not the same
xtable(t(head(y)), digits = 5)
stopifnot(abs(y[, 1] - y[, 2]) < 1e-8)
stopifnot(abs(sum(y[, 1]) - 1) < 1e-8)









## 
## small gibbs with slow use of the above
##
## 
nReads <- 10
nSNPs <- 5
seed <- 12
eps <- 0.1
out <- sim_reads(nReads, nSNPs, seed, eps)
s <- out$s
u <- out$u
true_s <- out$true_s
true_theta <- out$true_theta

## 
outX <- slow_gibbs(s, u, nReads, nits = 10000)
out <- outX[["out"]]
lls <- outX[["lls"]]
results <- outX[["results"]]
time <- outX[["time"]]; time


## 
png("phaser_cali.png", height = 4, width = 8, units = "in", res = 300)
par(mfrow = c(1, 2))
par(oma = c(0, 0, 0, 0))
par(mar = c(5, 4, 1, 1))
jitter <- runif(nits) / 5
plot(log(lls) + jitter, xlab = "Full iteration", ylab = "-log(ll) + jitter", pch = 16, col = rgb(red = 0, green = 0, blue = 1, alpha = 0.1))
par(mar = c(5, 4, 1, 1))
plot(out[, "prob"], out[, "counts"] / nits, xlab = "Exact probability", ylab = "Observed counts")
abline(a = 0, b = 1)
dev.off()

## some top examples
out[sort(sample(1:nrow(out), 10, prob = out[, "prob"] / sum(out[, "prob"]))), ]


## can exhaustively get all probabilities
y <- t(sapply(0:(2 ** nReads - 1), function(iRead) {
    z_proposed <- as.integer(intToBits(iRead)[1:nReads])
    c(z_proposed, p_z_given_o(z_proposed, s, u))
}))
f <- function(a, b) {
    sum(y[y[, 1] == a & y[, 2] == b, nReads + 1])
}

## truth
xtable(rbind(c(f(0, 0), f(0, 1)), c(f(1, 0), f(1, 1))))
## empirical
xtable(table(results[, 1], results[, 2]) / nits, digits = 4)












##
## more efficient Gibbs
##
nReads <- 10
nSNPs <- 5
seed <- 12
eps <- 0.1
out <- sim_reads(nReads, nSNPs, seed, eps)
s <- out$s
u <- out$u

outX <- fast_gibbs(s, u, nReads, nits = 10000)
out <- outX[["out"]]
lls <- outX[["lls"]]
results <- outX[["results"]]
time <- outX[["time"]]; time




##
## larger Gibbs
##
nReads <- 100
nSNPs <- 20
seed <- 12
eps <- 0.1
out <- sim_reads(nReads, nSNPs, seed, eps)
s <- out$s
u <- out$u

## 
outX <- fast_gibbs(s, u, nReads, nits = 10000)
out <- outX[["out"]]
lls <- outX[["lls"]]
results <- outX[["results"]]
time <- outX[["time"]]; time

## look at probabilities, specifically first SNP
out <- count.rows(results[, u == 1])
x <- array(0, c(nrow(out), nReads))
x[, 1:(ncol(out) - 1)] <- out[, -1]
a <- cbind(prob = apply(x, 1, p_z_given_o, s = s, u = u), out)
a[, "prob"] <- a[, "prob"] / sum(a[, "prob"])

## 
png("phaser_cali2.png", height = 4, width = 8, units = "in", res = 300)
par(mfrow = c(1, 2))
par(oma = c(0, 0, 0, 0))
par(mar = c(5, 4, 2, 1))
jitter <- runif(nits) / 5
plot(lls + jitter, xlab = "Full iteration", ylab = "-log(ll) + jitter", pch = 16, col = rgb(red = 0, green = 0, blue = 1, alpha = 0.1))
par(mar = c(5, 4, 2, 1))
plot(a[, "prob"], out[, "counts"] / nits, xlab = "Exact probability", ylab = "Observed counts", main = "First SNP only")
abline(a = 0, b = 1)
dev.off()



##
## block problems
##
png("phaser_block.png", height = 4, width = 8, units = "in", res = 300)
par(mfrow = c(1, 2))
par(oma = c(0, 0, 0, 0))
nReads <- 10
nSNPs <- 1
nits <- 10000
for(nReads in c(10, 100)) {
    s <- c(rep(0, nReads / 2), rep(1, nReads / 2))
    u <- rep(1, nReads)
    outX <- fast_gibbs(s, u, nReads, nits)
    results <- outX[["results"]]
    ##
    par(mar = c(2, 2, 2, 2))
    plot_rect_z(nits, nReads, results, main = paste0(nReads,  " reads"))
}
dev.off()

