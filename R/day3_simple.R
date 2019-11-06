sample_x_given_y <- function(y) {
    if (y == 1) {
        p <- 0.4 / (0.4 + 0.2)
    } else if (y == 2) {
        p <- 0.1 / (0.1 + 0.3)
    }
    return(sample(c(1, 2), 1, prob = c(p, 1 - p)))
}
sample_y_given_x <- function(x) {
    if (x == 1) {
        p <- 0.4 / (0.4 + 0.1)
    } else if (x == 2) {
        p <- 0.2 / (0.2 + 0.3)
    }
    return(sample(c(1, 2), 1, prob = c(p, 1 - p)))
}

set.seed(9919)
nits <- 10000
xs <- array(NA, nits)
ys <- array(NA, nits)
xs[1] <- sample(c(1, 2), 1)
ys[1] <- sample_y_given_x(xs[1])
for(it in 2:nits) {
    xs[it] <- sample_x_given_y(ys[it - 1])
    ys[it] <- sample_y_given_x(xs[it])    
}
table(xs) / nits
table(ys) / nits

## example
m <- rbind(head(xs), head(ys))
rownames(m) <- c("x", "y")
xtable(m, digits = 0)

xtable(table(xs, ys) / nits, digits = 3)

library("xtable")
