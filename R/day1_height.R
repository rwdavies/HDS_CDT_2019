setwd("/Users/robert davies/Google Drive/Presentations and conferences/2019_11_03 CDT Lectures/Day1 - EM")
source("/Users/robert davies/Google Drive/Presentations and conferences/2019_11_03 CDT Lectures/CDT_code/day1_functions.R")


## 
## unrelated to height - genotype probabilities
##
set.seed(9)
f <- function(p) {
    c(p ** 2, 2 * p * (1 - p), (1 - p) ** 2)
}
a <- sample(c(0, 1, 2), 600, prob = f(0.2), replace = TRUE)
b <- sample(c(0, 1, 2), 400, prob = f(0.7), replace = TRUE)
table(a)
table(b)
table(c(a, b))






##
## true values
## 
Tpi1 <- 0.5
Tpi2 <- 1 - Tpi1
Tmu1 <- 5 * 12 + 12
Tmu2 <- 5 * 12 + 4
Tsd1 <- 3
Tsd2 <- 3




##
## make plot, get numbers
##
for(i in 1:2) {
    set.seed(10)
    N <- 1000
    N_males <- rbinom(n = 1, size = N, prob = Tpi1)
    N_females <- N - N_males
    if (i == 1) {
        TX <- 5 * 12 + 10
        png_name <- "height_hist.png"
    } else {
        ## the numbers above
        TX <- Tmu1
        png_name <- "height_hist2.png"        
    }
    male_height <- rnorm(n = N_males, mean = TX, sd = Tsd1)
    female_height <- rnorm(n = N_females, mean = Tmu2, sd = Tsd2)
    x <- c(male_height, female_height)
    png(png_name, height = 4, width = 8, units = "in", res = 100)
    xlim <- range(c(male_height, female_height))
    ylim <- c(0, max(table(cut(x, breaks))))
    breaks <- seq(50, 85, 1)
    par(mfrow = c(1, 2))
    hist(c(male_height, female_height), col=rgb(1,1,1,0.25), xlim = xlim, main = "Histogram of height", xlab = "Height", breaks = breaks, ylim = ylim)
    hist(male_height, col=rgb(1,0,0,0.5), xlim = xlim, breaks = breaks, main = "Histogram of height", xlab = "Height", ylim = ylim)
    hist(female_height, col=rgb(0,0,1,0.5), add=TRUE, breaks = breaks, ylim = ylim)
    legend("topright", c("Male", "Female"), col = c(rgb(1,0,0,0.5), rgb(0, 0, 1,0.5)), lwd = 4)
    dev.off()
}


##
## initial parameters
## 
pi1 <- 0.6
pi2 <- 1 - pi1
mu1 <- 70
mu2 <- 60
sd1 <- 5
sd2 <- 5

## 
theta <- c(pi1, mu1, mu2, sd1, sd2)
Ttheta <- c(Tpi1, Tmu1, Tmu2, Tsd1, Tsd2)
Mtheta <- c(N_males / N, mean(male_height), mean(female_height), sd(male_height), sd(female_height))

pi1 <- 0.6
pi2 <- 1 - pi1
mu1 <- 70
mu2 <- 60
sd1 <- 5
sd2 <- 5

##
## fit model using constrained optimization, or manual approach
##
results_constrOptim <- perform_constrOptim()
results_manual <- manual_gd()

rbind(
    theta,
    Ttheta,
    Mtheta,
    results_constrOptim$par
)

xtable(rbind(
    c(ll(theta, x), theta),
    c(ll(Ttheta, x),    Ttheta),
    c(ll(Mtheta, x),    Mtheta),
    c(ll(results_constrOptim$par, x),    results_constrOptim$par)
))






ll(theta, x)
ll(Ttheta, x)
ll(Mtheta, x)
ll(results_constrOptim$par, x)

ll(results_manual$par, x)

dll(theta, x)
dll(Ttheta, x)
dll(Mtheta, x)
dll(results_constrOptim$par, x)
dll(results_manual$par, x)

## height here!
results <- results_constrOptim
png("height_hist3.png", height = 4, width = 8, units = "in", res = 100)
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
f(results$par[2], results$par[4], N_males, col = "red", lty = 3)
## 
hist(female_height, col=rgb(0,0,1,0.5), add=TRUE, breaks = breaks, ylim = ylim)
f(Tmu2, Tsd2, N_males, col = "blue", lty = 1)
f(mu2, sd2, N_males, col = "blue", lty = 2)
f(results$par[3], results$par[5], N_males, col = "blue", lty = 3)
legend("topright", c("Male", "Female"), col = c(rgb(1,0,0,0.5), rgb(0, 0, 1,0.5)), lwd = 4)
legend("topleft", c("True", "Initial", "Final"), col = "black", lwd = 4, lty = c(1, 2, 3))
dev.off()

png("height_hist4.png", height = 4, width = 8, units = "in", res = 100)
par(mfrow = c(1, 2))
##
## pi
## 
m <- seq(0.3, 0.7, 0.01)
y <- sapply(m, function(a) {
    b <- results$par
    b[1] <- a
    -ll(b, x) / log(10)
})
plot(m, y, xlab = "pi value", ylab = "Likelihood")
abline(h = max(y) - qchisq(p = 0.975, df = 5), col = "red")
##
## first normal 
##
m <- seq(70, 75, length.out = 100)
y <- sapply(m, function(a) {
    b <- results$par
    b[2] <- a
    -ll(b, x) / log(10)
})
plot(m, y, xlab = "pi value", ylab = "Likelihood")
abline(h = max(y) - qchisq(p = 0.975, df = 5), col = "red")
dev.off()
