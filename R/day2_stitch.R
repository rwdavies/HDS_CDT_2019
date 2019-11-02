setwd("/Users/robert davies/Google Drive/Presentations and conferences/2019_11_03 CDT Lectures/Day2 - HMM")
source("/Users/robert davies/Google Drive/Presentations and conferences/2019_11_03 CDT Lectures/CDT_code/read_functions.R")
library("STITCH")


##
## make reference panel (at random!)
##
set.seed(10110)
nSNPs <- 1000
read_length <- 100
chr_length <- 1000000
L <- sort(sample(chr_length, nSNPs, replace = FALSE))
K <- 8
nGen <- 10 ## how many generations of mixing
## make truth parameters
theta <- array(rbeta(n = nSNPs * K, shape1 = 0.2, shape2 = 0.2), c(K, nSNPs)) ## column - major on SNPs
eps <- 0.01 ## bound
theta[theta < eps] <- eps
theta[theta > (1 - eps)] <- (1 - eps)
## partly random
pi <- rbeta(n = 8, shape1 = 1, shape2 = 1)
pi <- pi / sum(pi)
## make alpha based on this plus some noise
alpha <- array(NA, c(K, nSNPs - 1))
for(k in 1:K) {
    alpha[k, ] <- pi[k] + 0.1 * runif(nSNPs - 1)
}
alpha <- t(t(alpha) / colSums(alpha)) ## normalize in a slightly ridiculous way
##
average_recomb_rate <- 10 ## cM / Mb, fairly high
rho <- diff(L) / 1e6 * (average_recomb_rate  / 100) ## recomb rate
## prod(exp(-rho * nGen)) ## yup


##
## checks
## 
stopifnot(min(alpha) > 0)
stopifnot(max(alpha) < 1)
stopifnot(min(theta) > 0)
stopifnot(max(theta) < 1)
stopifnot(sum(abs(colSums(alpha) - 1) > 1e-8) == 0) ## check all 0
stopifnot(abs(sum(pi) - 1) < 1e-8)
stopifnot(sum(rho < 0) == 0)


##
## sample some hidden states and then plot
##
set.seed(101)




##
## generative model - sample underlying truth for sample
##
set.seed(912)
q_true <- sample_hidden_state(
    rho = rho,
    nGen = nGen,
    pi = pi,
    alpha = alpha
)
## choose underlying haplotypes based on this
## (note - here assumption a little weaker, should sample haplotypes now arguably)
h <- theta[cbind(q_true, 1:nSNPs)]
## sample a true H once
sampled_truth_h <- as.integer(runif(nSNPs) < h)
gamma <- array(0, c(K, nSNPs))
gamma[cbind(q_true, 1:nSNPs)] <- 1
plot_posterior_probabilities(gamma, png_name = paste0("stitch_posterior_truth.png"), main = "Truth", xlab = "SNP index")

for(depth in c(0, 0.01, 0.1, 0.5)) {

    ##
    ## sample reads
    ##
    out <- simulate_sampleReads(
        nSNPs = nSNPs,
        read_length = read_length,
        depth = depth,
        chr_length = chr_length,
        L = L,
        seed = 410,
        truth_haplotypes = matrix(h, ncol = 1),
        bq_lambda = 30
    )
    sampleReads <- out[["sampleReads"]]
    nReads <- out[["nReads"]]
    
    ##
    ## run forward-backwards
    ##
    out <- run_forward_backwards(
        sampleReads = sampleReads,
        nReads = nReads,
        pi = pi,
        rho = rho,
        nGen = nGen,
        alpha = alpha,
        theta = theta
    )
    ## dosage <- out[["dosage"]]
    gamma <- out[["gamma"]]
    
    plot_posterior_probabilities(gamma, png_name = paste0("stitch_posterior_depth", depth, ".png"), main = paste0("Depth = ", depth), xlab = "SNP index")

}




##
## run EM (on theta only, since have posteriors)
##
set.seed(10110)
nSNPs <- 10
read_length <- 100
chr_length <- 1000
L <- sort(sample(chr_length, nSNPs, replace = FALSE))
K <- 4
nGen <- 10 ## how many generations of mixing
## make truth parameters
theta <- array(rbeta(n = nSNPs * K, shape1 = 0.2, shape2 = 0.2), c(K, nSNPs)) ## column - major on SNPs
eps <- 0.01 ## bound
theta[theta < eps] <- eps
theta[theta > (1 - eps)] <- (1 - eps)
## partly random
pi <- rbeta(n = K, shape1 = 1, shape2 = 1)
pi <- pi / sum(pi)
## make alpha based on this plus some noise
alpha <- array(NA, c(K, nSNPs - 1))
for(k in 1:K) {
    alpha[k, ] <- pi[k] + 0.1 * runif(nSNPs - 1)
}
alpha <- t(t(alpha) / colSums(alpha)) ## normalize in a slightly ridiculous way
##
average_recomb_rate <- 10 ## cM / Mb, fairly high
rho <- diff(L) / 1e6 * (average_recomb_rate  / 100) ## recomb rate
truth_theta <- theta


##
## sample truth
##
set.seed(919)
nSamples <- 200
depth <- 1
data <- lapply(1:nSamples, function(iSample) {
    q_true <- sample_hidden_state(
        rho = rho,
        nGen = nGen,
        pi = pi,
        alpha = alpha
    )
    out <- simulate_sampleReads(
        nSNPs = nSNPs,
        read_length = read_length,
        depth = depth,
        chr_length = chr_length,
        L = L,
        seed = 410 * iSample,
        truth_haplotypes = matrix(theta[cbind(q_true, 1:nSNPs)], ncol = 1),
        bq_lambda = 30
    )
    sampleReads <- out[["sampleReads"]]
    nReads <- out[["nReads"]]
    return(list(sampleReads = sampleReads, q = q_true))
})
## make truth dosages from this
truth_dosages <- array(0, c(nSamples, nSNPs))
for(iSample in 1:nSamples) {
    truth_dosages[iSample, ] <- theta[cbind(data[[iSample]]$q, 1:nSNPs)]
}


nits <- 20
lls <- array(0, nits)
## initialize new theta
theta <- array(runif(nSNPs * K), c(K, nSNPs)) ## column - major on SNPs
r2s <- array(0, nits)
for(it in 1:nits) {
    ##
    ## combined E-step, M-step
    ##
    gamma_numer <- array(0, c(K, nSNPs)) ## for M-step
    gamma_denom <- array(0, c(K, nSNPs))
    all_dosages <- array(0, c(nSamples, nSNPs))
    ##
    ## 
    for(iSample in 1:nSamples) {
        ##
        ## E-step bit
        ##
        sampleReads <- data[[iSample]]$sampleReads
        nReads <- length(sampleReads)
        out <- run_forward_backwards(
            sampleReads = sampleReads,
            nReads = nReads,
            pi = pi,
            rho = rho,
            nGen = nGen,
            alpha = alpha,
            theta = theta
        )
        dosage <- out[["dosage"]]
        gamma <- out[["gamma"]]
        lls[it] <- lls[it] + out[["ll"]]
        all_dosages[iSample, ] <- out[["dosage"]]
        ##
        ## M-step bit
        ##
        out <- update_theta(
            K = K,
            nSNPs = nSNPs,
            gamma = gamma,
            sampleReads = sampleReads,
            nReads = nReads,
            gamma_numer = gamma_numer,
            gamma_denom = gamma_denom
        )
        gamma_numer <- out[["gamma_numer"]]
        gamma_denom <- out[["gamma_denom"]]
    }
    ## now make transformation
    theta <- gamma_numer / gamma_denom
    ## check accuracy
    per_snp_r2 <- sapply(1:nSNPs, function(iSNP) {
        cor(all_dosages[, iSNP], truth_dosages[, iSNP]) ** 2
    })    
    r2s[it] <- mean(per_snp_r2)
    print(paste0("it = ", it, ", ll = ", round(lls[it], 2), ", r2 = ", round(r2s[it], 2)))
}

## argh, bug
cbind(it = 1:nits, ll = lls, r2 = r2s)
png("em_result.png", height = 4, width = 8, units = "in", res = 100)
par(mfrow = c(1, 3))
plot(-lls, xlab = "Iteration", ylab = "-log(ll)", type = "l")
plot(r2s, xlab = "Iteration", ylab = "r2", type = "l", ylim = c(0, 1))
plot(per_snp_r2, xlab = "SNP", ylab = "per-snp r2", ylim = c(0, 1))
dev.off()

round(truth_theta, 2)
round(theta, 2)
