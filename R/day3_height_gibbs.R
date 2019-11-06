setwd("/Users/robert davies/Google Drive/Presentations and conferences/2019_11_03 CDT Lectures/Day1 - EM")
source("/Users/robert davies/proj/HDS_CDT_2019/R/day3_height_functions.R")


##
## true values
## 
Tpi1 <- 0.5
Tpi2 <- 1 - Tpi1
Tmu1 <- 5 * 12 + 12
Tmu2 <- 5 * 12 + 4
Tsd1 <- 3
Tsd2 <- 3

set.seed(10)
N <- 1000
N_males <- rbinom(n = 1, size = N, prob = Tpi1)
N_females <- N - N_males
TX <- Tmu1
male_height <- rnorm(n = N_males, mean = TX, sd = Tsd1)
female_height <- rnorm(n = N_females, mean = Tmu2, sd = Tsd2)
x <- c(male_height, female_height)

## x <- rmix(n = 1000, pi = c(pi1, pi2), mu = c(mu1, mu2), s = c(sd1, sd2))


##
## initial parameters
## 
pi1 <- 0.6
pi2 <- 1 - pi1
mu1 <- 70
mu2 <- 60
## sd1 <- 5
## sd2 <- 5



res <- gibbs(x, 2)

x,k,niter =1000,muprior = list(mean=0,prec=0.1)
plot(res$mu[,1],ylim=c(-4,4),type="l")
lines(res$mu[,2],col=2)

## gibbs!

