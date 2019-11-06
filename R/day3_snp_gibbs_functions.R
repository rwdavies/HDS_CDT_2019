get_url <- function(chr) return(paste0("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr", chr, ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"))


get_data_per_chr <- function(chr) {
    message("get whole chr file")
    setwd(tempdir())
    system(paste0("curl -O ", get_url(chr)))
    system(paste0("curl -O ", get_url(chr), ".tbi"))
    ##
    message("subset")
    ## 
    Rfile <- tempfile()
    m <- cbind(chr, seq(1 + 1, 60e6, 10000), seq(1 + 100, 60e6, 10000))
    write.table(
        m,
        file = Rfile,
        row.names = FALSE,
        sep = "\t",
        quote = FALSE,
        col.names = FALSE
    )
    ##
    message("extract")
    ##
    outfile <- tempfile()    
    system(paste0("bcftools view -R ", Rfile, " ", basename(get_url(chr)), " | grep -v '^##' > ", outfile))
    ##
    message("load")
    ##
    data <- data.table::fread(outfile, data.table = FALSE)
    ## keep higher freq, at least 1 percent global freq
    m <- t(sapply(strsplit(data[, "INFO"], ";"), function(x) x[grep("AF", x)]))
    data <- data[as.numeric(sapply(strsplit(m[, 1], "="), function(x) x[[2]])) > 0.01, ]
    ## subset distance
    x <- floor(data[, "POS"] / 100000)
    data <- data[match(unique(x), x), ]
    ## make haplotypes from this
    h <- array(0, c(nrow(data), (ncol(data) - 9) * 2))
    for(i_col in 1:(ncol(data) - 9)) {
        h[, 2 * i_col - 1] <- as.integer(substr(data[, 9 + i_col ], 1, 1))
        h[, 2 * i_col - 0] <- as.integer(substr(data[, 9 + i_col ], 3, 3)        )
    }
    ## further subset, one variant every 50 kbp (nearly independent)
    ## return!
    message("get panel")
    ## 
    panelfile <- tempfile()
    panel_url <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
    system(paste0("curl -O ", panel_url))
    panel <- read.table(basename(panel_url), header = TRUE)
    ## 
    message("relabel")
    ## 
    ## relabel columns of h with population and number
    panel_pops <- as.character(panel[match(colnames(data), panel[, 1]), "pop"][-c(1:9)])
    x <- rep(panel_pops, each = 2)
    for(p in unique(x)) {
        x[x == p] <- paste0(p, 1:sum(x == p))
    }
    ## ok!
    colnames(h) <- x
    return(h)
}

p_theta_given_pi_H_Z <- function(pi, H, Z) {
    ## rbeta function
    K <- length(pi)
    T <- nrow(H)    
    theta <- array(NA, c(T, K))
    for(k in 1:K) {
        n_group <- sum(Z == k)
        if (n_group == 0) {
            num_of_ones <- rep(0, T)
        } else {
            num_of_ones <- rowSums(H[, Z == k, drop = FALSE])
        }
        theta[, k] <- rbeta(T, 1 + num_of_ones, 1 + n_group - num_of_ones)
    }
    return(theta)    
}
p_Z_given_pi_theta_H <- function(pi, theta, H) {
    N <- ncol(H)
    K <- length(pi)    
    Z <- array(NA, N)
    for(i in 1:N) {
        h <- H[, i]
        y0 <- colSums(log(1 - theta[h == 0, , drop = FALSE]))
        y1 <- colSums(log(theta[h == 1, , drop = FALSE]))
        y <- y0 + y1
        y <- y - max(y)
        b <- pi * exp(y)
        p <- b / sum(b)
        Z[i] <- sample(1:K, 1, prob = p)
    }
    return(Z)
}
## complete data likelihood
get_ll <- function(pi, theta, H, Z) {
    N <- ncol(H)
    ll <- sum(sapply(1:N, function(i) {
        h <- H[, i]
        y0 <- sum(log(1 - theta[h == 0, Z[i]]))
        y1 <- sum(log(theta[h == 1, Z[i]]))
        y2 <- log(pi[Z[i]])
        ## Z bit
        return(y0 + y1 + y2)
    }))
    ## Z
    return(ll)
}


## p_pi_given_H_Z_theta <- function(H, Z, theta) {
##     ## hmm?
##     pi <- array(0, K)
## }


run_gibbs <- function(H, Z_init = NULL, K = 3, nits = 10) {
    N <- ncol(H)
    T <- nrow(H)
    ## random probabilities
    pi <- runif(K)
    pi <- rep(1 / K, K)
    pi <- pi / sum(pi)
    theta <- array(runif(T * K), c(T, K))
    if (is.null(Z_init)) {
        Z <- sapply(1:N, function(i) sample(1:K, 1, prob = pi))
    } else {
        Z <- Z_init
    }
    lls <- array(NA, nits)
    Zstore <- array(NA, c(nits, N))
    lls[1] <- get_ll(pi, theta, H, Z)
    Zstore[1, ] <- Z
    ## 
    for(it in 2:nits) {
        theta <- p_theta_given_pi_H_Z(pi, H, Z)
        Z <- p_Z_given_pi_theta_H(pi, theta, H) 
        lls[it] <- get_ll(pi, theta, H, Z)
        Zstore[it, ] <- Z
    }
    return(
        list(
            lls = lls,
            Zstore = Zstore,
            Z = Z,
            theta = theta,
            pi = pi
        )
    )
}


plot_groups <- function(Zstore, lls, true_z, pngname = NULL, remove_first_lls = TRUE) {
    if (!is.null(pngname)) {
        png(pngname, height = 4, width = 8, units = "in", res = 150)
        par(oma = c(0, 0, 0, 0))
    }
    results <- rbind(true_z, NA, Zstore)
    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    par(mar = c(5, 4, 1, 1))
    ylim <- c(0, nrow(results) + 3)
    mult <- 1.4
    xlim <- c(0, mult * (ncol(results) + 1))
    plot(x = 0, y = 0, xlim = xlim, ylim = ylim, ylab = "Iteration", xlab = "", axes = FALSE, col = "white", main = main)
    ## plot each row
    ytop <- nrow(results):1 + 0.5
    ybottom <- ytop - 1
    for(icol in 1:ncol(results)) {
        col <- cbPalette[results[, icol] + 2] 
        rect(xleft = icol - 0.5, xright = icol + 0.5, ybottom = ybottom, ytop = ytop, col = col, border = NA)
    }
    ## add in abline
    for(irow in 1:nrow(results)) {    
        rect(xleft = 0.5, xright = ncol(results) + 0.5, ybottom = irow - 0.5, ytop = irow + 0.5, col = NA, border = "red")
    }
    ##text(x = 1:ncol(results), y = ylim[1] - 1, labels = 1:ncol(results), cex = 1.5)
    ##axis(2)
    ## add likelihoods
    ## place onto
    ## re-scale lls too!
    if (remove_first_lls) {
        lls <- lls[-1]
        y <- (nrow(results) - 3):1
    } else {
        y <- (nrow(results) - 2):1
    }
    lls <- lls - max(lls)
    yy <- lls
    yy <- yy - min(yy)
    yy <- xlim[2] / mult * 1.05 + yy * 1 / max(yy) * (xlim[2] / mult * (mult - 1.05))
    ## par(new = FALSE)
    points(x = yy, y = y, lwd = 2, type = "l", col = "black")
    ##
    labels <- round(seq(min(lls), max(lls), length.out = 2))
    at <- seq(xlim[2] / mult * 1.05 , xlim[2], length.out = 2)
    axis(side = 1, at = at, labels = labels)
    if (!is.null(pngname)) {
        dev.off()
    }
}


tiny_pca <- function(H, K) {
    H2 <- H
    for(irow in 1:nrow(H2)) {
        H2[irow, ] <- H2[irow, ] - mean(H2[irow, ])
    }
    out <- eigen(t(H2) %*% H2)
    v <- out$vectors[, 1:2]
    a <- kmeans(v, centers = K)
    return(a$cluster)
}
