simulate_sampleReads <- function(
    nSNPs = 100,
    read_length = 10000,
    depth = 40,
    chr_length = 100000,
    seed = 400,
    truth_haplotypes = NA,
    L = NA,
    bq_lambda = 8
) {
    if (depth == 0) {
        return(list(
            sampleReads = list(),
            nReads = 0
        ))
    }
    set.seed(seed)
    if (is.na(L[1])) {
        L <- sort(sample(chr_length, nSNPs))
    }
    nReads <- chr_length * depth / read_length
    if (is.na(truth_haplotypes[1])) {
        ## truth haplotypes
        th <- cbind(
            sample(c(0, 1), nSNPs, replace = TRUE),
            sample(c(0, 1), nSNPs, replace = TRUE)
        )
    } else {
        th <- truth_haplotypes
    }
    ## 
    ## sample truth reads
    ## 
    sampleReads <- lapply(1:nReads, function(iRead) {
        ## 
        middle <- sample(chr_length, 1)
        read_start <- floor(middle - read_length / 2)
        read_end <- floor(middle + read_length / 2)
        ##
        u <- which((read_start <= L) & (L <= read_end))
        if (length(u) == 0) {
            return(NULL)
        }
        ##
        h <- sample(ncol(truth_haplotypes), 1) ## which haplotype
        ## usually 0/1, can be 0-1
        ## underlying genotype
        Jr <- length(u)
        g <- sapply(1:Jr, function(j) sample(c(0, 1), size = 1, prob = c(1 - th[u[j], h], th[u[j], h])))
        bq <- rpois(Jr, lambda = bq_lambda)
        s <- sapply(1:Jr, function(j) {
            eps <- 10 ** (-bq[j] / 10)
            p <- rep(eps / 3, 4)
            p[g[j] + 1] <- 1 - eps
            s <- sample(0:3, 1, prob = p)
            return(s)
        })
        ##
        keep <- s <= 1
        u <- u[keep] ## 1-based
        s <- s[keep]
        g <- g[keep]
        bq <- bq[keep]
        Jr <- length(bq)
        cr <- u[getCentral(L[u])]
        return(
            list(
                u = u,
                s = s,
                g = g,
                bq = bq,
                Jr = Jr,
                h = h,
                cr = cr
            )
        )
    })
    sampleReads <- sampleReads[!sapply(sampleReads, is.null)]
    sampleReads <- sampleReads[!sapply(sampleReads, function(x) length(x[["cr"]]) == 0)]
    sampleReads <- sampleReads[order(unlist(sapply(sampleReads, function(x) x[["cr"]])))]
    nReads <- length(sampleReads)
    return(
        list(
            sampleReads = sampleReads,
            nReads = nReads
        )
    )
}
    

getCentral <- function(LL) {
    if(length(LL)==1) return(1)
    if(length(LL)==2) return(sample(2,1))
    if(length(LL)>=3) return(which.min(sapply(1:length(LL),function(i) sum(abs(LL[i]-LL[-i])))))
}



make_eMatGrid_t <- function(K, nReads, sampleReads, theta) {
    eMatGrid_t <- array(1, c(K, nSNPs))
    if (nReads == 0) {
        return(eMatGrid_t)
    }
    for(iRead in 1:nReads) {
        u <- sampleReads[[iRead]]$u
        cr <- sampleReads[[iRead]]$cr
        bq <- matrix(sampleReads[[iRead]]$bq, ncol = 1)
        s <- sampleReads[[iRead]]$s ## s = 0 means ref
        Jr <- sampleReads[[iRead]]$Jr
        ## convert
        bq[s == 0] <- -bq[s == 0]
        pr <- STITCH::convertScaledBQtoProbs(bq)
        for(j in 1:Jr) {
            for(k in 1:K) {
                eMatGrid_t[k, cr] <- eMatGrid_t[k, cr] *
                    (theta[k, u[j]] * pr[j, 2] + (1 - theta[k, u[j]]) *  pr[j, 1])
            }
        }
    }
    eMatGrid_t <- apply(eMatGrid_t, 2, function(x) {return(x * (1 / max(x)))})
    eMatGrid_t[eMatGrid_t < 1e-10] <- 1e-10
    return(eMatGrid_t)
}

forward_haploid <- function(
    eMatGrid_t,
    pi,
    theta,
    alpha,
    rho,
    nGen
) {
    ## re-name here
    alphaMat_t <- alpha
    x <- exp(-rho * nGen)
    transMatRate_t_H <- rbind(x, 1 - x)
    pi <- pi
    ##
    nSNPs <- ncol(theta)
    T <- nSNPs
    K <- nrow(theta)
    alphaHat_t <- array(0, c(K, T))
    c <- array(0, nSNPs)
    ## 
    for(k1 in 0:(K - 1)) {
        alphaHat_t[k1 + 1,0 + 1] <- pi[k1 + 1] * eMatGrid_t[k1 + 1, 0 + 1]
    }
    c[1] <- 1 / sum(alphaHat_t[, 1])
    alphaHat_t[, 1] <- alphaHat_t[, 1] * c[1]
    ## 
    ## 
    alphaConst <- 0
    for(t_0_based in 1:(T - 1)) {
        alphaConst <-
            transMatRate_t_H[1 + 1, t_0_based - 1 + 1] *
            sum(alphaHat_t[, t_0_based - 1 + 1])
        ## 
        alphaHat_t[, t_0_based + 1] <-
            eMatGrid_t[, t_0_based + 1] * (
                transMatRate_t_H[0 + 1, t_0_based - 1 + 1] *
                alphaHat_t[, t_0_based - 1 + 1] + 
                alphaConst *
                alphaMat_t[, t_0_based - 1 + 1]
            )
        ## 
        c[t_0_based + 1] <- 1 / sum(alphaHat_t[, t_0_based + 1])
        alphaHat_t[, t_0_based + 1] <- alphaHat_t[, t_0_based + 1] * c[t_0_based + 1]
    }
    ## check
    ## alphaHat2_t <- array(0, c(K, T))
    ## c2 <- array(0, nSNPs)
    ## f <- function(alpha) {
    ##     alphaMatCurrent_tc <- alpha
    ##     dim(alphaMatCurrent_tc) <- c(dim(alphaMatCurrent_tc), 1)
    ##     return(alphaMatCurrent_tc)
    ## }
    ## pi2 <- array(pi, c(length(pi), 1))
    ## Rcpp_run_forward_haploid(
    ##     alphaHat2_t,
    ##     c2,
    ##     eMatGrid_t,
    ##     f(alpha),
    ##     f(transMatRate_t_H),
    ##     pi2,
    ##     s = 0
    ## )
    return(
        list(
            alphaHat_t = alphaHat_t,
            c = c
        )
    )
}


backward_haploid <- function(
    eMatGrid_t,
    pi,
    theta,
    alpha,
    rho,
    nGen,
    c
) {
    ## re-name here
    alphaMat_t <- alpha
    x <- exp(-rho * nGen)
    transMatRate_t_H <- rbind(x, 1 - x)
    pi <- pi
    ##
    nSNPs <- ncol(theta)    
    T <- nSNPs
    K <- nrow(theta)
    betaHat_t <- array(0, c(K, T))
    betaHat_t[, T] <- c[T]
    ## 
    for(t_0_based in (T - 2):0) {
        t <- t_0_based + 1
        e_times_b <- eMatGrid_t[, t + 1] * betaHat_t[, t + 1]
        x <- transMatRate_t_H[1 + 1, t] * sum(alphaMat_t[, t] * e_times_b)
        betaHat_t[, t] <- c[t] * (x + transMatRate_t_H[0 + 1, t] * e_times_b)
    }
    ## check
    ## betaHat2_t <- array(0, c(K, T))
    ## betaHat2_t[, T] <- c[T]
    ## f <- function(alpha) {
    ##     alphaMatCurrent_tc <- alpha
    ##     dim(alphaMatCurrent_tc) <- c(dim(alphaMatCurrent_tc), 1)
    ##     return(alphaMatCurrent_tc)
    ## }
    ## pi2 <- array(pi, c(length(pi), 1))
    ## Rcpp_run_backward_haploid(
    ##     betaHat2_t,
    ##     c,
    ##     eMatGrid_t,
    ##     f(alpha),
    ##     f(transMatRate_t_H),
    ##     s = 0
    ## )
    return(list(betaHat_t = betaHat_t))
}


sample_hidden_state <- function(
    rho,
    nGen,
    pi,
    alpha
) {
    ## chose recombination points
    nRecombs <- rpois(n = 1, lambda = sum(rho * nGen)) ## 0-based
    ## where recomb takes place
    where_recomb <- which.max(runif(1) < c(cumsum(rho * nGen), 1))
    q <- array(NA, nSNPs)
    where_recomb <- sort(unique(c(1, where_recomb, nSNPs + 1)))
    for(i in 1:(length(where_recomb) - 1)) {
        if (i == 1) {
            ql <- sample(1:K, size = 1, prob = pi)
        } else {
            ql <- sample(1:K, size = 1, prob = alpha[, where_recomb[i] - 1])
        }
        q[where_recomb[i]:(where_recomb[i + 1] - 1)] <- ql
    }
    return(q)
}



run_forward_backwards <- function(
    sampleReads,
    nReads,
    pi,
    rho,
    nGen,
    alpha,
    theta
) {
    ## convert sampleReads into something easier to work with
    eMatGrid_t <- make_eMatGrid_t(K, nReads, sampleReads, theta)
    stopifnot(max(eMatGrid_t) <= 1)
    stopifnot(min(eMatGrid_t) > 0)
    ## 
    out <- forward_haploid(
        eMatGrid_t = eMatGrid_t,
        pi = pi,
        theta = theta,
        alpha = alpha,
        rho = rho,
        nGen = nGen
    )
    alphaHat_t <- out[["alphaHat_t"]]
    c <- out[["c"]]
    out <- backward_haploid(
        eMatGrid_t = eMatGrid_t,
        pi = pi,
        theta = theta,
        alpha = alpha,
        rho = rho,
        nGen = nGen,
        c = c
    )
    betaHat_t <- out[["betaHat_t"]]
    ## make gamma
    gamma <- alphaHat_t * betaHat_t
    for(t in 1:nSNPs) {
        gamma[, t] <- gamma[, t] / c[t]
    }
    stopifnot(sum(abs(colSums(gamma) - 1) > 1e-8) == 0) ## check all 0
    ## make dosage
    dosage <- colSums(theta * gamma)
    return(
        list(
            gamma = gamma,
            ll = sum(log(c)),
            dosage = dosage
        )
    )
}


## 
update_theta <- function(K, nSNPs, gamma, sampleReads, nReads, gamma_numer, gamma_denom) {
    ## 
    ## for every read, go over gamma
    ## remove influence of read on gamma and sum appropriately
    for(iRead in 1:length(sampleReads)) {
        u <- sampleReads[[iRead]]$u
        cr <- sampleReads[[iRead]]$cr
        bq <- matrix(sampleReads[[iRead]]$bq, ncol = 1)
        s <- sampleReads[[iRead]]$s ## s = 0 means ref
        Jr <- sampleReads[[iRead]]$Jr
        ## convert
        bq[s == 0] <- -bq[s == 0]
        pr <- STITCH::convertScaledBQtoProbs(bq)
        for(j in 1:Jr) {
            t <- u[j]
            for(k in 1:K) {
                left <- theta[k, u[j]] * pr[j, 2]
                right <- (1 - theta[k, u[j]]) *  pr[j, 1]
                gamma_numer[k, t] <- gamma_numer[k, t] + gamma[k, t] * left / (left + right)
                gamma_denom[k, t] <- gamma_denom[k, t] + gamma[k, t]
            }
        }
    }
    return(
        list(
            gamma_numer = gamma_numer,
            gamma_denom = gamma_denom
        )
    )
}


plot_multiple_hidden_states <- function() {


    set.seed(9919)
    qs <- sapply(1:10, function(i) {
        q_true <- sample_hidden_state(
            rho = rho,
            nGen = nGen,
            pi = pi,
            alpha = alpha
        )
    })

    ##
    ## first bit, hidden states
    ##
    png_name <- "mult_hidden.png"
    main <- NULL
    xlab <- "SNP"
    png(png_name, height = 3, width = 8, units = "in", res = 300)
    par(oma = c(0, 0, 0, 0))
    par(mar = c(5, 4, 1, 1))
    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    xlim <- c(1 - 0.5, nSNPs + 0.5)
    ylim <- c(1, ncol(qs) + 1)
    plot(x = 0, y = 0, xlim = xlim, ylim = ylim, axes = FALSE, xlab = xlab, ylab = "Hidden state", col = "white", main = main)
    axis(1)
    axis(2)
    for(iSample in 1:ncol(qs)) {
        gamma <- array(0, c(K, nSNPs))
        gamma[cbind(qs[, iSample], 1:nSNPs)] <- 1
        goodies <- which(rowSums(gamma) > 0)
        for(t in 1:nSNPs) {
            ybottom <- iSample
            ytop <- iSample
            for(k in 1:K) {
                ytop <- ytop + gamma[k, t]
                if (k %in% goodies) {
                    rect(xleft = t - 0.5, xright = t + 0.5, ybottom = ybottom, ytop = ytop, col = cbPalette[k], border = NA)
                }
                ybottom <- ytop
            }
        }
    }
    legend("topright", paste0("K", 1:K), col = cbPalette, lwd = 4, bg = "white")
    dev.off()
    
    ##
    ## actual SNP
    ##
    png_name <- "mult_haps.png"
    main <- NULL
    xlab <- "SNP"
    png(png_name, height = 3, width = 8, units = "in", res = 300)
    par(oma = c(0, 0, 0, 0))
    par(mar = c(5, 4, 1, 1))
    xlim <- c(1 - 0.5, nSNPs + 0.5)
    ylim <- c(1, ncol(qs) + 1)
    plot(x = 0, y = 0, xlim = xlim, ylim = ylim, axes = FALSE, xlab = xlab, ylab = "Haplotype", col = "white", main = main)
    axis(1)
    axis(2)
    for(iSample in 1:ncol(qs)) {
        h <- theta[cbind(qs[, iSample], 1:nSNPs)]
        for(t in 1:nSNPs) {
            ybottom <- iSample
            ymedium <- iSample + h[t]
            ytop <- iSample + 1
            rect(xleft = t - 0.5, xright = t + 0.5, ybottom = ybottom, ytop = ymedium, col = cbPalette[2], border = NA)
            rect(xleft = t - 0.5, xright = t + 0.5, ybottom = ymedium, ytop = ytop, col = cbPalette[3], border = NA)            
        }
    }
    dev.off()

    ##
    ## sample reads from this, visualize! alt vs ref base
    ##
    super_out <- lapply(1:ncol(qs), function(iSample) {
        h <- theta[cbind(qs[, iSample], 1:nSNPs)]
        out <- simulate_sampleReads(
            nSNPs = nSNPs,
            read_length = read_length,
            depth = 0.1,
            chr_length = chr_length,
            L = L,
            seed = 411 * iSample,
            truth_haplotypes = matrix(h, ncol = 1),
            bq_lambda = 30
        )
        return(out)
    })

    png_name <- "mult_reads.png"
    main <- NULL
    xlab <- "SNP"
    png(png_name, height = 6, width = 8, units = "in", res = 300)
    par(oma = c(0, 0, 0, 0))
    par(mar = c(5, 4, 1, 1))
    xlim <- c(1 - 0.5, nSNPs + 0.5)
    frac <- 0.2
    xlim <- round(c(1 - 0.5 + nSNPs * (1 - frac) / 2, (1 + frac) / 2 * nSNPs + 0.5)    )
    ylim <- c(1, ncol(qs) + 1)
    plot(x = 0, y = 0, xlim = xlim, ylim = ylim, axes = FALSE, xlab = xlab, ylab = "Haplotype", col = "white", main = main)
    axis(1)
    axis(2)
    for(iSample in 1:ncol(qs)) {
        sampleReads <- super_out[[iSample]]$sampleReads
        for(iRead in 1:length(sampleReads)) {
            sampleRead <- sampleReads[[iRead]]
            ## plot!
            u <- sampleRead$u
            s <- sampleRead$s
            Jr <- sampleRead$Jr
            for(j in 1:Jr) {
                rect(xleft = u, xright = u + 1, ybottom = iSample, ytop = iSample + 1, col = cbPalette[s + 2], border = NA)
            }
        }
    }
    dev.off()
    
    
    
}



plot_posterior_probabilities <- function(gamma, png_name, main = NULL, xlab = "Time") {
    png(png_name, height = 4, width = 8, units = "in", res = 300)
    par(oma = c(0, 0, 0, 0))
    par(mar = c(5, 4, 1, 1))
    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    xlim <- c(1 - 0.5, nSNPs + 0.5)
    ylim <- c(0, 1)
    plot(x = 0, y = 0, xlim = xlim, ylim = ylim, axes = FALSE, xlab = xlab, ylab = "Posterior probability", col = "white", main = main)
    axis(1)
    axis(2)
    for(t in 1:nSNPs) {
        ybottom <- 0
        ytop <- 0
        for(k in 1:K) {
            ytop <- ytop + gamma[k, t]
            rect(xleft = t - 0.5, xright = t + 0.5, ybottom = ybottom, ytop = ytop, col = cbPalette[k], border = NA)
            ybottom <- ytop
        }
        stopifnot(abs(ytop - 1) < 1e-8)
    }
    legend("topleft", paste0("K", 1:K), col = cbPalette, lwd = 4, bg = "white")
    dev.off()
}


