sim_reads <- function(nReads = 5, nSNPs = 5, seed = 10, eps = 0.1) {
    set.seed(seed)
    ##nReads <- 5
    ## eps <- 0.1 ## error rate
    ## nSNPs <- 5
    ## G is entirely het, so don't need hap2
    x <- sample(c(0, 1), nSNPs, replace = TRUE)
    true_theta <- rbind(x, 1 - x) ## hap1, hap2
    ## sample reads
    true_z <- sample(c(0, 1), nReads, replace = TRUE)
    u <- sample(1:nSNPs, nReads, replace = TRUE)
    true_z <- true_z[order(u)]
    u <- u[order(u)]
    ## observed sequencing reads
    true_s <- true_theta[cbind(true_z + 1, u)]
    s <- sapply(1:nReads, function(iRead) {
        b <- true_s[iRead]
        s <- sample(c(b, 1 - b), 1, prob = c(1 - eps, eps)) ## observed
        return(s)
    })
    return(
        list(
            s = s,
            u = u,
            true_z = true_z,
            true_theta = true_theta,
            true_s = true_s
        )
    )
}




## augmented probability - sum across all possible options
p_z_given_o_augmented <- function(z_proposed, s, u) {
    ## sum across all options for theta
    y <- sapply(0:(2 ** nSNPs - 1), function(iSNP) {
        x <- as.integer(intToBits(iSNP)[1:nSNPs])
        proposed_theta <- rbind(x, 1 - x)
        suggested_s <- proposed_theta[cbind(z_proposed + 1, u)]
        agree <- suggested_s == s
        (1 - eps) ** sum(agree) * (eps) ** sum(!agree)
    })
    sum(y) / 2 ** nSNPs
}

## alt - how dissimilar from each other?
p_z_given_o <- function(z_proposed, s, u) {
    ## 
    d <- (sapply(1:nSNPs, function(iSNP) {
        w <- u == iSNP
        o1 <- sum(z_proposed[w] == s[w])
        o2 <- sum(z_proposed[w] == (1 - s[w]))
        return(c(sum(o1), sum(o2)))
    }))
    ##
    c1 <- d[1, ]
    c2 <- d[2, ]
    y1a <- (1 - eps) ** (c1) * eps ** c2
    y1b <- (eps) ** (c1) * (1 - eps) ** c2
    prod(((y1a + y1b) / 2)[colSums(d) > 0])
}


## most likely option
## https://stat.ethz.ch/pipermail/r-help/2008-January/151489.html
count.rows <- function(x)   {
    order.x <- do.call(order,as.data.frame(x))
    equal.to.previous <-
        rowSums(x[tail(order.x,-1),] != x[head(order.x,-1),])==0
    tf.runs <- rle(equal.to.previous)
    counts <- c(1,
                unlist(mapply( function(x,y) if (y) x+1 else (rep(1,x)),
                              tf.runs$length, tf.runs$value )))
    counts <- counts[ c(diff(counts) <= 0, TRUE ) ]
    unique.rows <- which( c(TRUE, !equal.to.previous ) )
    cbind( counts, x[order.x[ unique.rows ], ,drop=F] )
}




slow_gibbs <- function(s, u, nReads, nits = 10000, block_resample_its = NULL) {
    ## do a lot of iterations, so can compare
    results <- array(0, c(nits, nReads))
    z_proposed <- sample(c(0, 1), nReads, replace = TRUE)
    lls <- array(0, nits)
    time <- system.time({
    for(it in 1:nits) {
        for(iRead in 1:nReads) {
            z <- z_proposed
            ## calculate two probabilities
            z[iRead] <- 0
            p1 <- p_z_given_o(z, s, u)
            z[iRead] <- 1
            p2 <- p_z_given_o(z, s, u)        
            ##
            z_proposed[iRead] <- sample(c(0, 1), 1, prob = c(p1, p2) / (p1 + p2))
        }
        if (it %in% block_resample_its) {
            for(uu in unique(u)) {
                z <- z_proposed
                w <- which(u == uu)
                ## calculate two probabilities
                p1 <- p_z_given_o(z, s, u)
                z[w] <- 1 - z[w]
                p2 <- p_z_given_o(z, s, u)        
                ##
                a <- sample(c("same", "switch"), 1, prob = c(p1, p2) / (p1 + p2))
                if (a == "switch") {
                    z_proposed[w] <- z[w]
                }
            }
        }
        lls[it] <- p_z_given_o(z_proposed, s, u)
        results[it, ] <- z_proposed
    }
    })
    out <- count.rows(results)
    out <- cbind(prob = apply(out[, -1], 1, p_z_given_o, s = s, u = u), out)
    return(
        list(
            lls = lls,
            results = results,
            out = out,
            time = time
        )
    )
}



fast_gibbs <- function(s, u, nReads, nits = 10000) {
    ##
    ## do a lot of iterations, so can compare
    ##
    results <- array(0, c(nits, nReads))
    z_proposed <- sample(c(0, 1), nReads, replace = TRUE)
    ##
    ## initialize
    ##
    time <- system.time({    
    d <- (sapply(1:nSNPs, function(iSNP) {
        w <- u == iSNP
        o1 <- sum(z_proposed[w] == s[w])
        o2 <- sum(z_proposed[w] == (1 - s[w]))
        return(c(sum(o1), sum(o2)))
    }))
    ##
    c1 <- d[1, ]
    c2 <- d[2, ]
    y1a <- (1 - eps) ** (c1) * eps ** c2
    y1b <- (eps) ** (c1) * (1 - eps) ** c2
    y <- log((y1a + y1b) / 2)
    ll <- sum(y)
    lls <- array(0, nits)
    ##
    for(it in 1:nits) {
        for(iRead in 1:nReads) {
            ##
            y1aa <- y1a[u[iRead]]
            y1bb <- y1b[u[iRead]]
            if (z_proposed[iRead] == s[iRead]) {
                y1aa <- y1aa * eps / (1 - eps)
                y1bb <- y1bb * (1 - eps) / eps
            } else {
                y1aa <- y1aa * (1 - eps) / eps
                y1bb <- y1bb * eps / (1 - eps)
            }
            ll_alt <- ll - y[u[iRead]] + log((y1aa + y1bb) / 2)
            if (z_proposed[iRead] == 0) {
                p1 <- exp(ll)            
                p2 <- exp(ll_alt)
            } else {
                p1 <- exp(ll_alt)
                p2 <- exp(ll)
            }
            ## 
            ## compare here as well
            if (it == 1) {
                z <- z_proposed
                z[iRead] <- 0
                p1ori <- p_z_given_o(z, s, u)
                z[iRead] <- 1
                p2ori <- p_z_given_o(z, s, u)
                m <- c(p1, p1ori, p2, p2ori); m * 1000
                stopifnot(abs((p1 - p1ori) / p1) < 1e-8)
                stopifnot(abs((p2 - p2ori) / p2) < 1e-8)
            }
            ##
            new_z <- sample(c(0, 1), 1, prob = c(p1, p2) / (p1 + p2))
            if (new_z != z_proposed[iRead]) {
                ll <- ll_alt
                y[u[iRead]] <- log((y1aa + y1bb) / 2)
                y1a[u[iRead]] <- y1aa
                y1b[u[iRead]] <- y1bb
            }
            z_proposed[iRead] <- new_z
        }
        ## 
        results[it, ] <- z_proposed
        lls[it] <- ll
        ##
    }
    })
    out <- count.rows(results)
    out <- cbind(prob = apply(out[, -1], 1, p_z_given_o, s = s, u = u), out)
    return(
        list(
            lls = lls,
            results = results,
            out = out,
            time = time
        )
    )
}


plot_rect_z <- function(nits, nReads, results, main = "", u = NULL) {
    ylim <- c(0, nits + 1)
    xlim <- c(0, nReads + 1)
    plot(x = 0, y = 0, xlim = xlim, ylim = ylim, xlab = "", ylab = "", axes = FALSE, col = "white", main = main)
    ybottom <- 1:(nits) - 0.5
    ytop <- ybottom + 0.5
    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    for(icol in 1:nReads) {
        col <- cbPalette[results[, icol] + 2] 
        rect(xleft = icol - 0.5, icol + 0.5, ybottom = ybottom, ytop = ytop, col = col, border = NA)
    }
    if (!is.null(u)) {
        text(x = 0:nReads, y = -0.5, labels = c("u", u))
    }
}


auto_corr_results <- function(results) {
    w <- sort(sample(1:(nits - 100), 100))
    m <- sapply(w, function(iSNP) {
        sapply(0:20, function(j) {
            cor(results[iSNP, ], results[iSNP + j, ]) ** 2
        })
    })
    return(rowMeans(m))
}
