## from 
## https://github.com/stephens999/fiveMinuteStats

# generate from mixture of normals
#' @param n number of samples
#' @param P a 2 by R matrix of allele frequencies
r_simplemix = function(n,P){
  R = ncol(P)
  z = sample(1:2,prob=c(0.5,0.5),size=n,replace=TRUE) #simulate z as 1 or 2
  x = matrix(nrow = n, ncol=R)
  for(i in 1:n){
    x[i,] = rbinom(R,rep(1,R),P[z[i],])
  }
  return(list(x=x,z=z))
}


#' @param x an R vector of data
#' @param P a K by R matrix of allele frequencies
#' @return the log-likelihood for each of the K populations
log_pr_x_given_P = function(x,P){
  tP = t(P) #transpose P so tP is R by K
  return(colSums(x*log(tP)+(1-x)*log(1-tP)))
}

normalize = function(x){return(x/sum(x))} #used in sample_z below

#' @param x an n by R matrix of data
#' @param P a K by R matrix of allele frequencies
#' @return an n vector of group memberships
sample_z = function(x,P){
  K = nrow(P)
  loglik_matrix = apply(x, 1, log_pr_x_given_P, P=P)
  lik_matrix = exp(loglik_matrix) 
  p.z.given.x = apply(lik_matrix,2,normalize) # normalize columns
  z = rep(0, nrow(x))
  for(i in 1:length(z)){
    z[i] = sample(1:K, size=1,prob=p.z.given.x[,i],replace=TRUE)
  }
  return(z)
}


#' @param x an n by R matrix of data
#' @param z an n vector of cluster allocations
#' @return a 2 by R matrix of allele frequencies
sample_P = function(x, z){
  R = ncol(x)
  P = matrix(ncol=R,nrow=2)
  for(i in 1:2){
    sample_size = sum(z==i)
    if(sample_size==0){
      number_of_ones=rep(0,R) 
    } else {
      number_of_ones = colSums(x[z==i,])
    }
    P[i,] = rbeta(R,1+number_of_ones,1+sample_size-number_of_ones) 
  }
  return(P)
}

gibbs = function(x,niter = 100){
  z = sample(1:2,nrow(x),replace=TRUE)
  res = list(z = matrix(nrow=niter, ncol=nrow(x)))
  res$z[1,]=z 
  
  for(i in 2:niter){
    P = sample_P(x,z)
    z = sample_z(x,P)
    res$z[i,] = z
  }
  return(res)
}

if (1 == 0) {
    
    set.seed(33)
    P <- rbind(
        c(0.5,0.5,0.5,0.5,0.5,0.5),
        c(0.001,0.999,0.001,0.999,0.001,0.999)
    )
    ## so exceptionally distributed
    sim <- r_simplemix(n=50,P)
    x <- sim$x ## columns are SNPs
    
    res <- gibbs(x, 100)
    
    table(res$z[1,],sim$z)
    table(res$z[nrow(res$z),],sim$z) ## better
    image(t(res$z))
    
    res <- gibbs(x, 100)
}
