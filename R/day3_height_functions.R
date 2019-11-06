##
## https://stephens999.github.io/fiveMinuteStats/gibbs2.html
##

## set.seed(33)

# generate from mixture of normals
#' @param n number of samples
#' @param pi mixture proportions
#' @param mu mixture means
#' @param s mixture standard deviations
rmix = function(n,pi,mu,s){
  z = sample(1:length(pi),prob=pi,size=n,replace=TRUE)
  x = rnorm(n,mu[z],s[z])
  return(x)
}

normalize <- function(x){return(x/sum(x))}
  
#' @param x an n vector of data
#' @param pi a k vector
#' @param mu a k vector
sample_z = function(x,pi,mu){
    dmat = outer(mu,x,"-") # k by n matrix, d_kj =(mu_k - x_j)
    p.z.given.x = as.vector(pi) * dnorm(dmat,0,1) 
    p.z.given.x = apply(p.z.given.x,2,normalize) # normalize columns
    z = rep(0, length(x))
    for(i in 1:length(z)){
        z[i] = sample(1:length(pi), size=1,prob=p.z.given.x[,i],replace=TRUE)
    }
    return(z)
}


#' @param z an n vector of cluster allocations (1...k)
#' @param k the number of clusters
sample_pi = function(z,k){
    counts = colSums(outer(z,1:k,FUN="=="))
    pi = gtools::rdirichlet(1,counts+1)
    return(pi)
}

#' @param x an n vector of data
#' @param z an n vector of cluster allocations
#' @param k the number o clusters
#' @param prior.mean the prior mean for mu
#' @param prior.prec the prior precision for mu
sample_mu = function(x, z, k, prior){
    df = data.frame(x=x,z=z)
    mu = rep(0,k)
    for(i in 1:k){
        sample.size = sum(z==i)
        sample.mean = ifelse(sample.size==0,0,mean(x[z==i]))
        
        post.prec = sample.size+prior$prec
        post.mean = (prior$mean * prior$prec + sample.mean * sample.size)/post.prec
        mu[i] = rnorm(1,post.mean,sqrt(1/post.prec))
    }
    return(mu)
}

gibbs = function(x,k,niter =1000,muprior = list(mean=0,prec=0.1)){
    pi <- rep(1/k,k) # initialize
    ## mu <- rnorm(k,0,10)
    mu <- c(70, 60)
    z <- sample_z(x,pi,mu)
    res <- list(mu=matrix(nrow=niter, ncol=k), pi = matrix(nrow=niter,ncol=k), z = matrix(nrow=niter, ncol=length(x)))
    res$mu[1,]=mu
    res$pi[1,]=pi
    res$z[1,]=z 
    for(i in 2:niter){
        pi = sample_pi(z,k)
        mu = sample_mu(x,z,k,muprior)
        z = sample_z(x,pi,mu)
        res$mu[i,] = mu
        res$pi[i,] = pi
        res$z[i,] = z
    }
    return(res)
}
