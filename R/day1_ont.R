setwd("/Users/robert davies/Google Drive/Presentations and conferences/2019_11_03 CDT Lectures/Day1 - EM")
source("/Users/robert davies/Google Drive/Presentations and conferences/2019_11_03 CDT Lectures/CDT_code/read_functions.R")



##
## first part, simulation to look at allele balance
##
set.seed(9009)
L <- 1e6
het <- 0.001
nSNPs <- 2 * L * het
error <- 0.127
## conditional on site being het,
is_snp <- sample(c(FALSE, TRUE), L, prob = c(1 - het, het), replace = TRUE)
a1 <- rep(1, L); a1[is_snp] <- sample(c(1, 2), sum(is_snp), replace = TRUE)
a2 <- rep(1, L); a2[is_snp] <- sample(c(1, 2), sum(is_snp), replace = TRUE)
## sample reads conditional on that
depth <- rpois(L, 60)
a1_depth <- sapply(1:L, function(i) {
    rbinom(1, size = depth[i], prob = 0.5)
})
a2_depth <- depth - a1_depth

f <- function(a1, a1_depth) {
    t(sapply(1:L, function(iSNP) {
        a <- a1[iSNP]
        d <- a1_depth[iSNP]
        if (a == 1) {
            e <- error
        } else {
            e <- 1 - error
        }
        nref <- rbinom(n = 1, size = d, prob = e)
        nalt <- d - nref
        return(c(nref, nalt))
    }))
}

a1_alleles <- f(a1, a1_depth)
a2_alleles <- f(a2, a2_depth)

png("ont_depth.png", height = 4, width = 12, units = "in",res = 150)
par(mfrow = c(1, 3))
## all sites
x <- (a1_alleles[, 1] + a2_alleles[, 1]) /  rowSums(a1_alleles + a2_alleles)
breaks <- seq(0, 1, length.out = 21)
ylim <- c(0, max(table(cut(x, breaks))))
hist(x[(a1 + a2) == 2], breaks = breaks, col = rgb(1,0,0,0.5), ylim = ylim)
hist(x[(a1 + a2) == 3], breaks = breaks, col = rgb(0,1,0,0.5), ylim = ylim, add = TRUE)
hist(x[(a1 + a2) == 4], breaks = breaks, col = rgb(0,0,1,0.5), ylim = ylim, add = TRUE)
legend(
    "topright",
    c("Hom ref", "Het", "Hom alt"),
    col = c(rgb(1,0,0,0.5), rgb(0, 1, 0,0.5), rgb(0, 0, 1, 0.5)),
    lwd = 4
)
## true positive sites
x <- (a1_alleles[, 1] + a2_alleles[, 1]) /  rowSums(a1_alleles + a2_alleles)
breaks <- seq(0, 1, length.out = 21)
ylim <- c(0, max(table(cut(x[(a1 + a2) >= 3], breaks))))
hist(x[(a1 + a2) == 2], breaks = breaks, col = rgb(1,0,0,0.5), ylim = ylim)
hist(x[(a1 + a2) == 3], breaks = breaks, col = rgb(0,1,0,0.5), ylim = ylim, add = TRUE)
hist(x[(a1 + a2) == 4], breaks = breaks, col = rgb(0,0,1,0.5), ylim = ylim, add = TRUE)
legend(
    "topright",
    c("Hom ref", "Het", "Hom alt"),
    col = c(rgb(1,0,0,0.5), rgb(0, 1, 0,0.5), rgb(0, 0, 1, 0.5)),
    lwd = 4
)
## with (latent!) haplotypes
col <- array(NA, L)
col[(a1 + a2) == 2] <- rgb(1,0,0,0.5)
col[(a1 + a2) == 3] <- rgb(0,1,0,0.5)
col[(a1 + a2) == 4] <- rgb(0,0,1,0.5)
plot(a1_alleles[, 1] / rowSums(a1_alleles), a2_alleles[, 1] / rowSums(a2_alleles), col = col, xlab = "A1 ref proportion", ylab = "A2 ref proportion")
## 
dev.off()













##
## implement phasing EM in R
##
out <- simulate_sampleReads(
    nSNPs = 100,
    read_length = 10000,
    depth = 40,
    chr_length = 100000,
    seed = 400
)
sampleReads <- out[["sampleReads"]]
nReads <- out[["nReads"]]






##
## define mixed E-step and M-step
##
perform_mixed_E_and_M_step <- function(theta, sampleReads) {
    K <- 2
    log_likelihood <- 0
    nSNPs_output <- 1
    N_r_output <- 1
    calculate_updates <- TRUE
    calculate_read_probabilities <- TRUE
    N_r <- length(sampleReads)
    if (calculate_updates) {
        nSNPs_output <- nSNPs
    }
    if (calculate_read_probabilities) {
        N_r_output <- N_r
    }
    t_break <- 0
    t_reading_offset <- 0
    t_writing_offset <- 0
    eHapsUpdate_numer <- array(0, c(nSNPs, K))
    eHapsUpdate_denom <- array(0, c(nSNPs, K))
    p_reads <- array(0, c(1, N_r))        
    p_reads_given_hap_k <- array(0, c(N_r, 2))
    p_h_given_O <- array(0, c(N_r, 2))
    eHapsCurrent <- theta
    ##
    ## mixed e-step and m-step!
    ##
    for(r in 1:N_r) {
        ##
        ## E-step
        ##
        sampleRead <- sampleReads[[r]]
        J_r <- sampleRead$Jr - 1 ## 0-based
        bq <- matrix(sampleRead$b, ncol = 1)
        s <- sampleRead$s
        bq[s == 0] <- -bq[s == 0]
        u <- sampleRead$u - 1
        pr <- STITCH::convertScaledBQtoProbs(bq)
        flip_k <- array(0, J_r + 1)
        if (t_break > 0) {
            flip <- (u + 1) > t_break
            flip_k[flip] <- 1            
        }
        p_read_given_hap_k_components <- array(0, c(J_r + 1, K))
        for(k in 1:K) {
            use_k <- k
            for(j in 1:(J_r + 1)) {
                u_j <- u[j] + 1                
                u_j_reading <- u[j] + 1 - t_reading_offset
                if (flip_k[j] == 1)
                    use_k <- 3 - k
                p_read_given_hap_k_components[j, k] <-
                    eHapsCurrent[u_j_reading, use_k] * pr[j, 2] +
                    (1 - eHapsCurrent[u_j_reading, use_k]) * pr[j, 1]
            }
        }
        p_read_given_hap_k <- apply(p_read_given_hap_k_components, 2, prod)
        p_read <- 0.5 * sum(p_read_given_hap_k)
        log_likelihood <- log_likelihood + -sum(log(p_read))
        if (calculate_read_probabilities) {
            p_reads[1, r] <- p_read
            p_reads_given_hap_k[r, ] <- p_read_given_hap_k
        }
        ##
        ## M-step
        ##
        for(j in 1:(J_r + 1)) {
            u_j_reading <- u[j] + 1 - t_reading_offset
            u_j_writing <- u[j] + 1 - t_writing_offset                    
            for(k in 1:K) {
                use_k <- k
                if (flip_k[j] == 1)
                    use_k <- 3 - k
                p_r_no_j <- p_read_given_hap_k[k] / p_read_given_hap_k_components[j,k]
                p_h_and_g_1 <- p_r_no_j * 0.5 * pr[j, 2] *
                    eHapsCurrent[u_j_reading, use_k] / p_read
                p_h_and_g_2 <- p_r_no_j * 0.5 * pr[j, 1] *
                (1 - eHapsCurrent[u_j_reading, use_k]) / p_read
                if (calculate_read_probabilities)
                    p_h_given_O[r, k] <- p_h_and_g_1 + p_h_and_g_2;
                eHapsUpdate_numer[u_j_writing, k] <-
                    eHapsUpdate_numer[u_j_writing, k] + p_h_and_g_1
                eHapsUpdate_denom[u_j_writing, k] <-
                    eHapsUpdate_denom[u_j_writing, k] +
                    p_h_and_g_1 + p_h_and_g_2
            }
        }
    }
    new_theta <- eHapsUpdate_numer / eHapsUpdate_denom
    ll <- sum(log(p_reads)) ## from before
    return(
        list(
            new_theta = new_theta,
            p_reads_given_hap_k = p_reads_given_hap_k,
            ll = ll
        )
    )
}
calculate_pse <- function(test, truth) {
    which_sites <- rowSums(truth == 0 | truth == 1) == 2 & rowSums(truth) == 
        1
    truth <- truth[which_sites, ]
    test <- test[which_sites, ]
    if (nrow(test) == 0) 
        return(NA)
    test[, 1] <- as.integer(round(test[, 1]))
    test[, 2] <- as.integer(round(test[, 2]))
    choose_at_random <- which(rowSums(test) != 1)
    if (length(choose_at_random) > 0) {
        test[choose_at_random, ] <- 0
        r <- sample(c(1, 2), length(choose_at_random), replace = TRUE)
        test[cbind(choose_at_random, r)] <- 1
    }
    if (test[1, 1] != truth[1, 1]) {
        test <- test[, c(2, 1)]
    }
    n_bad <- sum(diff(abs(test[, 1] - truth[, 1])) != 0)
    return(c(n_bad/(nrow(test) - 1)))
}


results <- lapply(1:10, function(init) {
    ##
    ## do 10 starts
    ## 
    theta <- matrix(runif(nSNPs * 2), ncol = 2)
    nits <- 10
    pses <- lls <- array(0, nits)
    for(it in 1:10) {
        ##
        ## mixed E and M step
        ##
        out <- perform_mixed_E_and_M_step(theta, sampleReads)
        theta <- out$new_theta
        lls[it] <- out$ll
        pses[it] <- calculate_pse(test = theta, truth = th)
    }
    return(
        list(
            pses = pses,
            lls = lls,
            p_reads_given_hap_k = out$p_reads_given_hap_k,
            theta = theta
        )
    )
})

## 
png("ont_sim1.png", height = 4, width = 8, units = "in",res = 150)
lls <- results[[1]]$lls
pses <- results[[1]]$pses
par(mfrow = c(1, 2))
plot(lls, xlab = "Iteration", ylab = "ll", type = "l")
plot(pses, xlab = "Iteration", ylab = "Phase switch error", type = "l")
dev.off()

##
png("ont_sim2.png", height = 4, width = 8, units = "in",res = 150)
lls <- sapply(results, function(x) x$lls)
pses <- sapply(results, function(x) x$pses)
par(mfrow = c(1, 2))
plot(lls[, 1], xlab = "Iteration", ylab = "ll", type = "l", ylim = range(lls))
for(i in 2:ncol(lls)) {
    lines(lls[, i])
}
plot(pses[, 1], xlab = "Iteration", ylab = "Phase switch error", type = "l", ylim = range(pses))
for(i in 2:ncol(pses)) {
    lines(pses[, i])
}
dev.off()

## 
png("ont_sim3.png", height = 4, width = 4, units = "in",res = 150)
plot(theta[, 1], theta[, 2], xlab = "P(Allele 1 is alt)", ylab = "P(Allele 2 is alt)")
dev.off()


##
## ugh, not quite working
##

## also, p_read vs truth?
p_reads_given_hap_k <- results[[1]]$p_reads_given_hap_k
p_reads_given_hap_k <- p_reads_given_hap_k / rowSums(p_reads_given_hap_k)
theta <- results[[1]]$theta
test <- theta
truth <- th
which_sites <- rowSums(truth == 0 | truth == 1) == 2 & rowSums(truth) == 1
truth <- truth[which_sites, ]
test <- test[which_sites, ]
test[, 1] <- as.integer(round(test[, 1]))
test[, 2] <- as.integer(round(test[, 2]))
if (test[1, 1] != truth[1, 1]) {
    test <- test[, c(2, 1)]
}
flip_at <- which(diff(abs(test[, 1] - truth[, 1])) > 0)
flip <- 1:which(which_sites)[flip_at]
truth_h <- sapply(sampleReads, function(x) x$h)
truth_h[flip] <- 3 - truth_h[flip]

## right, what if align theta
hist(p_reads_given_hap_k[truth_h == 1, 1])
hist(p_reads_given_hap_k[truth_h == 2, 2])

plot(p_reads_given_hap_k[truth_h == 1, 1])
plot(p_reads_given_hap_k[truth_h == 2, 1])

sapply(results, function(x) x$pse)


test = theta
truth = th

