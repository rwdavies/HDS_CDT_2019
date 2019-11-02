ll <- function(theta, x) {
    pi1 <- theta[1]
    pi2 <- 1 - pi1
    mu1 <- theta[2]
    mu2 <- theta[3]
    sd1 <- theta[4]
    sd2 <- theta[5]
    l <- log(
        pi1 * dnorm(x = x, mean = mu1, sd = sd1) + 
        pi2 * dnorm(x = x, mean = mu2, sd = sd2)
    )
    return(-sum(l))
}
dll <- function(theta, x) {
    pi1 <- theta[1]
    pi2 <- 1 - pi1
    mu1 <- theta[2]
    mu2 <- theta[3]
    sd1 <- theta[4]
    sd2 <- theta[5]
    d1 <- dnorm(x = x, mean = mu1, sd = sd1)
    d2 <- dnorm(x = x, mean = mu2, sd = sd2)
    bottom <- pi1 * d1 + pi2 * d2
    dll_pi1 <- sum((1 / bottom) * (d1 - d2))
    ## dll_pi2 <- sum((1 / bottom) * d2)
    dll_mu1 <- sum((1 / bottom) * (pi1 * d1 * (x - mu1) / (2 * sd1 ** 2)))
    dll_mu2 <- sum((1 / bottom) * (pi2 * d2 * (x - mu2) / (2 * sd2 ** 2)))    
    ##
    dll_sigma1 <- sum((1 / bottom) * (d1 * ( (-1 / sd1) + ((x - mu1) ** 2 / sd1 ** 3))))
    dll_sigma2 <- sum((1 / bottom) * (d2 * ( (-1 / sd2) + ((x - mu2) ** 2 / sd2 ** 3))))    
    return(-c(dll_pi1, dll_mu1, dll_mu2, dll_sigma1, dll_sigma2))
}
## f <- function(pi1, x) {
##     pi2 <- 1 - pi1
##     mu1 <- Mtheta[2]
##     mu2 <- Mtheta[3]
##     sd1 <- Mtheta[4]
##     sd2 <- Mtheta[5]
##     d1 <- dnorm(x = x, mean = mu1, sd = sd1)
##     d2 <- dnorm(x = x, mean = mu2, sd = sd2)
##     bottom <- pi1 * d1 + pi2 * d2
##     r <- sum((1 / bottom) * (d1 - d2))
##     return(r)
## }
## optimize(function(pi1, x) { abs(f(pi1, x))}, interval = c(0, 1), x = x)
## a <- seq(0.2, 0.8, 0.001)
## plot(a, sapply(a, f, x = x))
## abline(v = 0.5)
## abline(h = 0)


perform_constrOptim <- function() {
    ##
    k <- 2
    p <- 5
    ui <- matrix(0, nrow = k, ncol = p)
    ## ui[1, 1:2] <- c(1, 1)
    ## ui[2, 1:2] <- c(-1, -1)
    ui[1, 1] <- 1
    ui[2, 1] <- -1
    ## ui[5, 2] <- 1
    ## ui[6, 2] <- -1
    ##
    ci <- matrix(0, nrow = k, ncol = 1)
    ci[1] <- 0
    ci[2] <- -1
    ## eps <- 0
    ## ci[1] <- 1 - eps
    ## ci[2] <- -1 - eps
    ## ci[3] <- 0
    ## ci[4] <- -1
    ## ci[5] <- 0
    ## ci[6] <- -1
    results <- optim(par = theta, fn = ll, gr = dll, x = x, control = list(trace = 100, factr = 1e7), lower = c(0.1, 50, 50, 1, 1), upper = c(0.9, 80, 80, 10, 10), method = "L-BFGS-B")
    return(results)
}


manual_gd <- function() {
    ## simple? start, do direction, line fit, done
    start_theta <- theta
    n_total_its <- 1000
    results <- array(0, c(n_total_its, 7))
    results[1, ] <- c(1, ll(start_theta, x), start_theta)
    total_it <- 2
    for(total_it in 2:n_total_its) {
        ## derivative
        start_theta <- results[total_it - 1, 3:7]
        der <- -dll(start_theta, x = x)
        ## line search
        eps <- 0.01
        is_done <- FALSE
        ##
        it <- 0
        lower <- c(0.1, 50, 50, 1, 1)
        upper <- c(0.9, 100, 100, 10, 10)
        d <- (upper - start_theta) / der
        min_mult <- - 0.2 * min(abs(d))
        d <- (lower - start_theta) / der
        max_mult <- 0.2 * min(abs(d))
        r <- array(0, c(20, 7))
        ## do equally spaced to start
        for(mult in seq(min_mult, max_mult, length.out = 11)) {
            it <- it + 1
            proposed_theta <- start_theta + der * mult
            proposed_ll <- ll(proposed_theta, x = x)
            r[it, ] <- c(mult, proposed_ll, proposed_theta)
        }
        ## refine
        while(it < 20) {
            ## chose new one to try
            it <- it + 1
            mult <- mean(r[order(r[1:(it - 1), 2])[1:2], 1])
            proposed_theta <- start_theta + der * mult
            proposed_ll <- ll(proposed_theta, x = x)        
            r[it, ] <- c(mult, proposed_ll, proposed_theta)
        }
        ## reset
        results[total_it, ] <- c(total_it, r[nrow(r), 2], r[nrow(r), 3:7])
        if ((total_it %% 100) == 0) {
            print(results[total_it, ])
        }
    }
    return(
        list(
            par = results[nrow(results), 3:7],
            results = results
        )
    )
}


EM.iter <- function(theta, x) {
    ## 
    pi1 <- theta[1]
    pi2 <- 1 - pi1
    mu1 <- theta[2]
    mu2 <- theta[3]
    sd1 <- theta[4]
    sd2 <- theta[5]
    ## E-step: compute E_{Z|X,w0}[I(Z_i = k)]
    d1 <- dnorm(x = x, mean = mu1, sd = sd1)
    d2 <- dnorm(x = x, mean = mu2, sd = sd2)
    z_ik <- cbind(pi1 * d1, pi2 * d2)
    z_ik <- z_ik / rowSums(z_ik) ## yes, rowSums is right!
    ## M-step - calculate updates
    Nk <- colSums(z_ik)
    N <- length(x)
    Npi1 <- Nk[1] / N
    Npi2 <- Nk[2] / N
    Nmu1 <- (1 / Nk[1]) * sum(z_ik[, 1] * x)
    Nmu2 <- (1 / Nk[2]) * sum(z_ik[, 2] * x)
    Nsd1 <- sqrt((1 / Nk[1]) * sum(z_ik[, 1] * (x - mu1) ** 2))
    Nsd2 <- sqrt((1 / Nk[2]) * sum(z_ik[, 2] * (x - mu2) ** 2))
    ## 
    Ntheta <- c(
        Npi1,
        Nmu1,
        Nmu2,
        Nsd1,
        Nsd2
    )
    return(Ntheta)
}

em <- function(theta_init) {

    ## have some initial theta
    theta <- theta_init
    ## store log-likehoods for each iteration
    maxits <- 1000
    results <- array(0, c(maxits, 7))
    results[1, ] <- c(1, ll(theta, x), theta)
    delta.ll <- 1
    it <- 1
    while ((1e-5 < delta.ll) & (it < maxits)) {
        it <- it + 1
        theta <- EM.iter(theta, x)
        results[it, ] <- c(it, ll(theta, x), theta)
        delta.ll <- results[it - 1, 2] - results[it, 2]
    }
    results <- results[1:it, ]
    colnames(results) <- c("it", "ll", "pi1", "mu1", "mu2", "sd1", "sd2")
    xtable(head(results))

    ## at the start, z's?
    png("heightEM1.png", height = 4, width = 8, units = "in", res = 100)
    par(mfrow = c(1, 2))
    plot(-results[, 2], xlab = "Iteration", ylab = "Negative log likelihoodd", type = "l")
    plot(log10(-diff(results[, 2])), xlab = "Iteration", ylab = "log10 Difference in likelihood", type = "l")
    dev.off()

    ## compare
## height here!
    final_theta <- list(par = results[nrow(results), 3:7])
png("heightEM2.png", height = 4, width = 8, units = "in", res = 100)
par(mfrow = c(1, 2))
xlim <- range(c(male_height, female_height))
breaks <- seq(floor(min(x)), ceiling(max(x)))
ylim <- c(0, 100)
h1 <- hist(c(male_height, female_height), col=rgb(1,1,1,0.25), xlim = xlim, main = "Histogram of height", xlab = "Height", breaks = breaks, ylim = ylim)
h2 <- hist(male_height, col=rgb(1,0,0,0.5), xlim = xlim, breaks = breaks, main = "Histogram of height", xlab = "Height", ylim = ylim)
f <- function(mu, sd, N, col = "red", lty = 1) {
    a <- seq(breaks[1], tail(breaks, 1), 0.05)
    y <- dnorm(a, mean = mu, sd = sd) * N
    ## m <- N / max(h1$counts) 
    lines(a, y, col = col, lty = lty, lwd = 3)
}
f(Tmu1, Tsd1, N_males, col = "red", lty = 1)
f(mu1, sd1, N_males, col = "red", lty = 2)
f(final_theta$par[2], final_theta$par[4], N_males, col = "red", lty = 3)
## 
hist(female_height, col=rgb(0,0,1,0.5), add=TRUE, breaks = breaks, ylim = ylim)
f(Tmu2, Tsd2, N_males, col = "blue", lty = 1)
f(mu2, sd2, N_males, col = "blue", lty = 2)
f(final_theta$par[3], final_theta$par[5], N_males, col = "blue", lty = 3)
legend("topright", c("Male", "Female"), col = c(rgb(1,0,0,0.5), rgb(0, 0, 1,0.5)), lwd = 4)
legend("topleft", c("True", "Initial", "Final"), col = "black", lwd = 4, lty = c(1, 2, 3))
dev.off()
    
    
}
