setwd("/Users/robert davies/Google Drive/Presentations and conferences/2019_11_03 CDT Lectures/Day2 - HMM")
library("xtable")
library("parallel")


set.seed(757) ## R 3.6.1
##
## simple markov example
##
pi <- c(0.6, 0.4)
A <- rbind(
    c(0.8, 0.2),
    c(0.6, 0.4)
)


sample_markov <- function(T = 10) {
    K <- length(pi)
    x <- array(0, T)
    x[1] <- sample(1:K, 1, prob = pi)
    for(t in 2:T) {
        x[t] <- sample(1:K, 1, prob = A[, x[t - 1]])
    }
    return(x)
}


set.seed(92)
xtable(rbind(
    sample_markov(),
    sample_markov(),
    sample_markov()
), digits = 0)






##
## simple hidden markov example
##
pi <- c(0.6, 0.4)
A <- rbind(
    c(0.8, 0.2),
    c(0.6, 0.4)
)
## three observations
B <- rbind(
    c(0.6, 0.3, 0.1),
    c(0.15, 0.35, 0.5)
)
## 
T <- 10
K <- length(pi)


sample_hidden_markov <- function(T = 10) {
    K <- length(pi)
    M <- ncol(B)
    q <- o <- array(0, T)
    q[1] <- sample(1:K, 1, prob = pi)
    o[1] <- sample(1:M, 1, prob = B[q[1], ])
    for(t in 2:T) {
        q[t] <- sample(1:K, 1, prob = A[, q[t - 1]])
        o[t] <- sample(1:M, 1, prob = B[q[t], ])
    }
    return(rbind(q, o))
}

set.seed(191)
xtable(sample_hidden_markov(), digits = 0)
xtable(sample_hidden_markov(), digits = 0)
third <- sample_hidden_markov()
xtable(third, digits = 0)

## calculate likelihood
o <- third["o", ]
##
ls <- array(0, 2 ** T)
for(i in 0:((2 ** T) - 1)) {
    ## 
    q <- as.integer(intToBits(i)[1:T]) + 1
    l <- pi[q[1]]
    for(t in 2:T) {
        l <- l * A[q[t - 1], q[t]]
    }
    for(t in 1:T) {
        l <- l * B[q[t], o[t]]
    }
    ls[i + 1] <- l
}
sum(ls)


make_forward <- function(o) {
    alpha <- array(0, c(K, T))
    for(k in 1:K) {
        alpha[k, 1] <- pi[k] * B[k, o[1]]
    }
    for(t in 2:T) {
        for(k_to in 1:K) {
            for(k_from in 1:K) {
                alpha[k_to, t] <- alpha[k_to, t] +
                    (alpha[k_from, t - 1] * A[k_from, k_to])
            }
            alpha[k_to, t] <- alpha[k_to, t] * B[k_to, o[t]]
        }
    }
    return(alpha)
}


## it works!
alpha <- make_forward(o)
sum(alpha[, T])


make_backward <- function(o) {
    beta <- array(0, c(K, T))
    beta[, T] <- 1
    for(t in (T - 1):1) {
        for(k_to in 1:K) {
            for(k_from in 1:K) {
                beta[k_from, t] <- beta[k_from, t] +
                    A[k_from, k_to] * B[k_to, o[t + 1]] * beta[k_to, t + 1]
            }
        }
    }
    return(beta)
}

beta <- make_backward(o)
## Yup, adds up
colSums(alpha * beta)
## yup
p_o_given_lambda <- sum(beta[, 1] * pi * B[, o[1]])
gamma <- alpha * beta / p_o_given_lambda

## visualize
png("hmm_posterior.png", height = 4, width = 8, units = "in", res = 300)
par(oma = c(0, 0, 0, 0))
par(mar = c(5, 4, 1, 1))
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
xlim <- c(1 - 0.5, T + 0.5)
ylim <- c(0, 1)
plot(x = 0, y = 0, xlim = xlim, ylim = ylim, axes = FALSE, xlab = "Time", ylab = "Posterior probability", col = "white")
axis(1)
axis(2)
for(t in 1:T) {
    ybottom <- 0
    ytop <- 0
    for(k in 1:K) {
        ytop <- ytop + gamma[k, t]
        rect(xleft = t - 0.5, xright = t + 0.5, ybottom = ybottom, ytop = ytop, col = cbPalette[k + 1], border = NA)
        ybottom <- ytop
    }
}
legend("topright", c("sunny", "rainy"), col = cbPalette[2:3], lwd = 2)
dev.off()


sum(alpha[, T])
xtable(log10(alpha))
xtable(log10(beta))
## check range 
.Machine$double.xmin
.Machine$double.xmax





##
## viterbi
## 
viterbi <- function(o) {
    delta <- array(0, c(K, T))
    phi <- array(0, c(K, T))
    for(k in 1:K) {
        delta[k, 1] <- pi[k] * B[k, o[1]]
    }
    phi[, 1] <- 0
    for(t in 2:T) {
        for(k in 1:K) {
            delta[k, t] <- max(delta[, t - 1] * A[, k]) * B[k, o[t]]
            phi[k, t] <- which.max(delta[, t - 1]* A[, k])
        }
    }
    ##
    q_star <- array(0, T)
    q_star[T] <- which.max(delta[, T])
    for(t in (T - 1):1) {
        q_star[t] <- phi[q_star[t + 1], t + 1]
    }
    return(q_star)
}

## 
q_star <- viterbi(o)
## confirm using exhaustive approach
p_o_q_given_l <- array(0, 2 ** T)
p_q_given_l <- array(0, 2 ** T)
for(i in 0:((2 ** T) - 1)) {
    ## 
    q <- as.integer(intToBits(i)[1:T]) + 1
    l <- pi[q[1]]
    for(t in 2:T) {
        l <- l * A[q[t - 1], q[t]]
    }
    p_q_given_l[i + 1] <- l
    for(t in 1:T) {
        l <- l * B[q[t], o[t]]
    }
    p_o_q_given_l[i + 1] <- l
}
stopifnot(abs(sum(    p_q_given_l) - 1) < 1e-10)
stopifnot(abs(sum(p_o_q_given_l) - p_o_given_lambda) < 1e-10)
x <- p_o_q_given_l * p_q_given_l
x <- x / sum(x)
i_wrong <- which.max(x / sum(x)) ## for argmax_q P(Q | O, \lambda)
wrong_q_star <- as.integer(intToBits(i_wrong- 1)[1:T]) + 1
i <- which.max(p_o_q_given_l) ## for argmax_q P(Q, O | \lambda)
q_star_exhaustive <- as.integer(intToBits(i - 1)[1:T]) + 1
##

xtable(rbind(o, q_star, wrong_q_star), digits = 0)
m <- cbind(p_o_q_given_l, p_q_given_l, p_q_given_o_l = x)[c(i, i_wrong), ]
xtable(m, digits = 3)
xtable(m, digits = -10)
## add probabilities






##
## scaled variables
##
## work with scaled alpha, beta
make_forward_scaled <- function(o, pi, A, B) {
    alphaHat <- array(0, c(K, T))
    c <- array(0, T)
    for(k in 1:K) {
        alphaHat[k, 1] <- pi[k] * B[k, o[1]]
    }
    c[1] <- 1 / sum(alphaHat[, 1])
    alphaHat[, 1] <- alphaHat[, 1] * c[1]
    for(t in 2:T) {
        for(k_to in 1:K) {
            for(k_from in 1:K) {
                alphaHat[k_to, t] <- alphaHat[k_to, t] +
                    (alphaHat[k_from, t - 1] * A[k_from, k_to])
            }
            alphaHat[k_to, t] <- alphaHat[k_to, t] * B[k_to, o[t]]
        }
        c[t] <- 1 / sum(alphaHat[, t])
        alphaHat[, t] <- alphaHat[, t] * c[t]
    }
    return(list(alphaHat = alphaHat, c = c))
}

out <- make_forward_scaled(o, pi, A, B)
alphaHat <- out$alphaHat
c <- out$c
xtable(rbind(alphaHat, c), digits = 2)

## same 
exp(sum(-log(c)))
p_o_given_lambda


make_backward_scaled <- function(o, c, pi, A, B) {
    betaHat <- array(0, c(K, T))
    betaHat[, T] <- 1
    betaHat[, T] <- c[T] * betaHat[, T]
    for(t in (T - 1):1) {
        for(k_to in 1:K) {
            for(k_from in 1:K) {
                betaHat[k_from, t] <- betaHat[k_from, t] +
                    A[k_from, k_to] * B[k_to, o[t + 1]] * betaHat[k_to, t + 1]
            }
        }
        betaHat[, t] <- c[t] * betaHat[, t]        
    }
    return(betaHat)
}

betaHat <- make_backward_scaled(o, c, pi, A, B)
xtable(rbind(alphaHat, c, betaHat), digits = 2)

gamma <- alphaHat * betaHat
for(t in 1:T) {
    gamma[, t] <- gamma[, t] / c[t]
}

## all 1
stopifnot(sum(colSums(gamma) != 1) > 0)








##
## EM here
##
Estep <- function(o, epi, eA, eB) {
    ## 
    out <- make_forward_scaled(o, epi, eA, eB)
    alphaHat <- out[["alphaHat"]]
    c <- out[["c"]]
    betaHat <- make_backward_scaled(o, c, epi, eA, eB)
    gamma <- alphaHat * betaHat
    for(t in 1:T) {
        gamma[, t] <- gamma[, t] / c[t]
    }
    ## note - sometimes (like here), more pragmatic to make xi on the fly
    xi <- array(0, c(T, K, K))
    for(t in 1:(T - 1)) {
        for(i in 1:K) {
            for(j in 1:K) {
                xi[t, i, j] <- alphaHat[i, t] * eA[i, j] * eB[j, o[t + 1]] * betaHat[j, t + 1]
            }
        }
    }
    ## note - make xi 
    ## p_o_given_l <- exp(sum(-log(c)))
    return(list(alphaHat = alphaHat, betaHat = betaHat, c = c, gamma = gamma, xi = xi))
}
Mstep <- function(o, alphaHat, betaHat, c, gamma, xi, epi, eA, eB) {
    new_epi <- gamma[, 1]
    ## new A
    new_A_numer <- new_A_denom <- array(0, dim(A))    
    for(t in 1:(T - 1)) {
        ## sum(xi) should be 1,
        ## rowSums(xi) should equal gamma[, 1]
        for(i in 1:K) {
            for(j in 1:K) {
                new_A_numer[i, j] <- new_A_numer[i, j] + xi[t, i, j]
                new_A_denom[i, ] <- new_A_denom[i, ] + xi[t, i, j]
            }
            ## new_A_denom[i, ] <- new_A_denom[i, ] + sum(xi[t, i, ])
        }
    }
    new_A <- new_A_numer / new_A_denom
    ##
    ## new B
    ##
    new_B_numer  <- array(0, dim(B))
    new_B_denom <- array(0, dim(B))
    for(t in 1:T) {
        for(k in 1:K) {
            j <- o[t]
            new_B_numer[k, j] <- new_B_numer[k, j] + gamma[k, t]
            new_B_denom[k, ] <- new_B_denom[k, ] + gamma[k, t]
        }
    }
    return(
        list(
            new_epi = new_epi,
            new_A_numer = new_A_numer,
            new_A_denom = new_A_denom,            
            new_B_numer = new_B_numer,
            new_B_denom = new_B_denom
        )
    )
}


##
## do first one just using o
## 
set.seed(91)
o <- third["o", ]
T <- length(o)
eA <- array(runif(K * K), c(K, K))
for(k in 1:K) {
    eA[k, ] <- eA[k, ] / sum(eA[k, ])
}
eB <- array(runif(prod(dim(B))), dim(B))
for(k in 1:K) {
    eB[k, ] <- eB[k, ] / sum(eB[k, ])
}
epi <- runif(K)
epi <- epi / sum(K)
nits <- 30
ll <- array(0, nits)
for(it in 1:nits) {
    ##
    ## E-step calculate posterior
    ##
    out <- Estep(o, epi, eA, eB)
    alphaHat <- out[["alphaHat"]]
    betaHat <- out[["betaHat"]]
    c <- out[["c"]]
    gamma <- out[["gamma"]]
    xi <- out[["xi"]]
    ll[it] <- sum(-log(c))
    ##
    ## M-step! 
    ##
    out <- Mstep(o, alphaHat, betaHat, c, gamma, xi, epi, eA, eB)
    ## rename here
    epi <- out[["new_epi"]]
    eA <- out[["new_A_numer"]] / out[["new_A_denom"]]
    eB <- out[["new_B_numer"]] / out[["new_B_denom"]]    
}

png("hmm_em.png", height = 4, width = 4, units = "in", res = 300)
par(oma = c(0, 0, 0, 0))
par(mar = c(5, 4, 1, 1))
plot(ll, xlab = "Iteration", ylab = "Log likelihood")
dev.off()

## build manually
epi
eA
eB
xtable(eA, digits = 3)
xtable(eB, digits = -1)




##
## EM using multiple, much longer sequences. also more clearly obvious A, B. slow, poorly written
##
pi <- c(0.4, 0.6)
A[1, ] <- c(0.9, 0.1)
A[2, ] <- c(0.2, 0.8)
B[1, ] <- c(0.6, 0.3, 0.05)
B[2, ] <- c(0.2, 0.1, 0.7)
T <- 100
N <- 1000
out <- parallel::mclapply(c(411, 919, 110, 303, 309, 110, 303, 556), mc.cores = 4, function(seed) {
    set.seed(seed)
    O <- array(0, c(N, T))
    for(n in 1:N) {
        O[n, ] <- sample_hidden_markov(T = T)["o", ]
    }
    ## initialize, partly randomly
    eA <- array(runif(K * K), c(K, K))
    for(k in 1:K) {
        eA[k, ] <- eA[k, ] / sum(eA[k, ])
    }
    eB <- array(runif(prod(dim(B))), dim(B))
    for(k in 1:K) {
        eB[k, ] <- eB[k, ] / sum(eB[k, ])
    }
    epi <- runif(K)
    epi <- epi / sum(epi)
    nits <- 30
    ll <- array(0, nits)
    params <- list()
    params <- append(params, list(list(epi = epi, eA = eA, eB = eB)))
    ##
    for(it in 1:nits) {
        new_epi <- array(0, length(epi))
        new_eA_numer <- new_eA_denom <- array(0, dim(eA))
        new_eB_numer <- new_eB_denom <- array(0, dim(eB))    
        for(n in 1:N) {
            ##
            ## E-step calculate posterior
            ##
            o <- O[n, ]
            out <- Estep(o, epi, eA, eB)
            alphaHat <- out[["alphaHat"]]
            betaHat <- out[["betaHat"]]
            c <- out[["c"]]
            gamma <- out[["gamma"]]
            xi <- out[["xi"]]
            ll[it] <- ll[it] + sum(-log(c))
            ##
            ## M-step! 
            ##
            out <- Mstep(o, alphaHat, betaHat, c, gamma, xi, epi, eA, eB)
            ## rename here
            new_epi <- new_epi + out[["new_epi"]]
            new_eA_numer <- new_eA_numer + out[["new_A_numer"]]
            new_eA_denom <- new_eA_denom + out[["new_A_denom"]]
            new_eB_numer <- new_eB_numer + out[["new_B_numer"]]
            new_eB_denom <- new_eB_denom + out[["new_B_denom"]]        
        }
        ## check again
        ## for(i in 1:4) {
        ##     print(sum(sapply(1:N, function(n) {
        ##         if (i == 1) {out <- Estep(O[n, ], epi, eA, eB)}
        ##     if (i == 2) {out <- Estep(O[n, ], (new_epi / N), eA, eB)} 
        ##     if (i == 3) {out <- Estep(O[n, ], (new_epi / N), (new_eA_numer / new_eA_denom), eB)}
        ##     if (i == 4) {out <- Estep(O[n, ], (new_epi / N), (new_eA_numer / new_eA_denom), (new_eB_numer / new_eB_denom))        }
        ##     c <- out[["c"]]
        ##     return(sum(-log(c)))
        ##     })))
        ## }
        ## switch over here
        epi <- new_epi / N
        eA <- new_eA_numer / new_eA_denom
        eB <- new_eB_numer / new_eB_denom
        params <- append(params, list(list(epi = epi, eA = eA, eB = eB))    )
        stopifnot(abs(sum(epi) - 1) < 1e-8)
        stopifnot(sum(abs(rowSums(eA) - 1)) < 1e-8)
        stopifnot(sum(abs(rowSums(eB) - 1)) < 1e-8)    
    }
    return(
        list(
            epi = epi,
            eA = eA,
            eB = eB,
            ll = ll,
            params = params
        )
    )
})

## plot the different likelihoods
## surprinsingly large difference (can recover)
ll <- sapply(out, function(x) x$ll)
ll[nrow(ll), ]
i_best <- which.max(ll[nrow(ll), ])
ylim <- range(ll)
xlim <- range(1:nrow(ll))
png("mult_likelihoods.png", height = 4, width = 8, units = "in", res = 300)
par(mfrow = c(1, 2))
## likelihoods straight up
plot(x = 0, y = 0, col = "white", xlab = "Iteration", ylab = "-ll", axes = FALSE, xlim = xlim, ylim = ylim, main = "ll")
for(i in 1:ncol(ll)) {
    lines(ll[, i], col = cbPalette[i])
}
axis(1)
axis(2)
## likelihoods change
ylim <- range(apply(ll[-1, ], 2, diff))
plot(x = 0, y = 0, col = "white", xlab = "Iteration", ylab = "-ll", axes = FALSE, xlim = xlim, ylim = ylim, main = "Change in ll\n(excluding initialization)")
for(i in 1:ncol(ll)) {
    lines(diff(ll[-1, i]), col = cbPalette[i])
}
axis(1)
axis(2)
dev.off()


##
## for the best result, add to table manually
##
pi
xtable(A, digits = 3)
xtable(B, digits = -1)
j <- length(out[[i_best]]$params)
out[[i_best]]$params[[j]]$epi
xtable(out[[i_best]]$params[[j]]$eA, digits = 3)
xtable(out[[i_best]]$params[[j]]$eB, digits = -1)


##
##
##
plot(ll, xlab = "Iteration", ylab = "Log likelihood")
plot(log10(diff(ll)), xlab = "Iteration", ylab = "Log likelihood")

png("hmm_em_multiple.png", height = 4, width = 4, units = "in", res = 300)
plot(ll, xlab = "Iteration", ylab = "Log likelihood")
dev.off()
