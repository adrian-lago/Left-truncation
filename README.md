#### left truncation ####
# Lynden-Bell estimator
# sample = (truncation, t.to event)
lb.estimator <- function(sample){
  
  dist.fail <- unique(sample[,2])
  di <- ni <- prob <- numeric( length(dist.fail) )
  j <- 1
  
  for (t in sort(dist.fail)) {
    di[j] <- sum( sample[,2] == t )
    ni[j] <- sum( sample[,1] <= t & t <= sample[,2]  )  
    j <- j+1
  }
  
  cum.prod <- c(1,cumprod(1-di/ni))
  
  for (i in 1:length(dist.fail) ){
    prob[i] <- cum.prod[i]-cum.prod[i+1]
  }
  
  
  
  lista <- list('fail.time'= sort( dist.fail ) , 'estimation'= cumprod(1-di/ni),
                'prob'=prob)
  return(lista)
}

lb.estimator.mod <- function(sample){
  
  dist.fail <- unique(sample[,2])
  di <- ni <- prob <- numeric( length(dist.fail) )
  j <- 1
  
  for (t in sort(dist.fail)) {
    di[j] <- sum( sample[,2] == t )
    ni[j] <- sum( sample[,1] <= t & t <= sample[,2]  ) + 1 # to avoid holes
    j <- j+1
  }
  
  # ni[length(dist.fail)] <- 1
  
  cum.prod <- c(1,cumprod(1-di/ni))
  
  for (i in 1:length(dist.fail) ){
    prob[i] <- cum.prod[i]-cum.prod[i+1]
  }
  
  
  
  lista <- list('fail.time'= sort( dist.fail ) , 'estimation'= cumprod(1-di/ni),
                'prob'=prob)
  return(lista)
}


# sample = ( truncation, failure times )

ks.stat <- function(sample1,sample2,rep=F){ # rep = T o F para representar o no
  
  est1 <- lb.estimator.mod(sample1); nodos1 <- est1$fail.time; widehat.s1 <- est1$estimation
  n <- nrow(sample1) 
  
  est2 <- lb.estimator.mod(sample2); nodos2 <- est2$fail.time; widehat.s2 <- est2$estimation
  m <- nrow(sample2)
  
  nodos  <- c(0,unique(sort(c(nodos1,nodos2),decreasing=F))) # total de tiempos de fallo
  
  eval1 <- eval2 <- numeric(length(nodos))
  
  for (i in 1:length(nodos)){
    
    if (nodos[i] < nodos1[1]){
      eval1[i] <- 1
    }
    
    if (nodos[i] < nodos2[1]){
      eval2[i] <- 1
    }
    
    if (nodos[i] >= nodos1[length(nodos1)]){
      eval1[i] <- 0
    }
    
    if (nodos[i] >= nodos2[length(nodos2)]){
      eval2[i] <- 0
    }
    
    for ( j in 1:(length(nodos1) - 1) ){
      
      if(nodos1[j] <= nodos[i] & nodos[i] < nodos1[j+1]){
        eval1[i] <- widehat.s1[j]
      } 
      
    }
    
    for ( j in 1:(length(nodos2)-1) ){
      
      if(nodos2[j] <= nodos[i] & nodos[i] < nodos2[j+1]){
        eval2[i] <- widehat.s2[j]
      } 
    }
  }
  
  statistic <- max(abs(eval1-eval2)) # statistic
  tt1 <- tt2 <- tt12 <- 0
  
  # agujeros
  #for (i in 1:(length(nodos1)-1)){if(widehat.s1[i]==0){tt1 <- min(which(widehat.s1==0))}}
  #for (i in 1:(length(nodos2)-1)){if(widehat.s2[i]==0){tt2 <- min(which(widehat.s2==0))}}
  #if (tt1 != 0 & tt2 != 0){tt12 <- min(tt1,tt2)}
  
  # representacion si rep=T
  if (rep==T){
    plot(nodos,eval1,type='s',lwd=2,ylim=c(-0.1,1.1),col='firebrick3',xlab='t',ylab='Estimated survival functions',main='')
    points(nodos,eval2,type='s',lwd=2,col='forestgreen' )
  }
  lista <- list('statistic'=statistic,'tt1'=tt1,'tt2'=tt2,'tt12'=tt12)
  return(lista)
}

sim.exp <- function(n,prob=0.5,lambda=1,umin=0,umax=1,eps=0){
  x <- rexp(n=n,rate=lambda) + eps
  u <- runif(n=n,min=umin,max=umax)
  z <- rbinom(n=n,size=1,prob=prob)
  for (i in 1:n){
    while(u[i] > x[i]){
      u[i] <- runif(1,umin,umax)
      x[i] <- rexp(1,rate=lambda) + eps
    }
  }
  # ord <- sort(x,index.return=T)
  return( data.frame('truncation'=u,'times'=x,'group'=z) )
}

sim.gamma <- function(n,prob=0.5,shape=1,rate=1,umin=0,umax=1){
  u <- runif(min=0,max=umax,n=n)
  x <- rgamma(n=n,shape=shape,rate=rate)
  for (i in 1:n){
    while(u[i] > x[i]){
      u[i] <- runif(min=0,max=umax,n=1)
      x[i] <- rgamma(n=1,shape=shape,rate=rate)
    }
  }
  group <- rbinom(n,1,prob)
  lista <- data.frame( 'truncation'=u , 'times'=x , 'group'= group)
  return(lista)
}

sim.unif <- function(n,prob=0.5,min.x=0,max.x=1,min.u=0,max.u=1){
  x <- runif(n,min.x,max.x)
  u <- runif(n,min.u,max.u)
  z <- rbinom(n,1,prob) # group
  for (i in 1:n){
    while (u[i] > x[i]){
      u[i] <- runif(1,min.u,max.u)
      x[i] <- runif(1,min.x,max.x)
    }
  }
  return(data.frame('truncation'=u,'times'=x,'group'=z))
}

sim.wei <- function(n,prob=0.5,shape=1,scale=1,min.u=0,max.u=1){
  x <- rweibull(n=n,shape=shape,scale=scale)
  u <- runif(n=n,min=min.u,max=max.u)
  for (i in 1:n){
    while(u[i] > x[i]){
      u[i] <- runif(1,min.u,max.u)
      x[i] <- rweibull(1,shape,scale)
    }
  }
  group <- rbinom(n,1,prob)
  return(  data.frame('truncation'=u,'times'=x,'group'=group)  )
}

sim.beta <- function(n,prob=0.5,shape1=1,shape2=1,umin=0,umax=1){
  u <- runif(n,umin,umax)
  x <- rbeta(n,shape1,shape2)
  for (i in 1:n){
    while(u[i] > x[i]){
      u[i] <- runif(n=1,umin,umax)
      x[i] <- rbeta(n=1,shape1,shape2)
    }
  }
  group <- rbinom(n,1,prob)
  return(  data.frame('truncation'=u,'times'=x,'group'=group)  )
}

sim.norm <- function(n,prob=0.5,mu=0,sigma=1,umin=0,umax=1){
  u <- runif(n,umin,umax)
  x <- rnorm(n,mu,sigma)
  
  for (i in 1:n){
    while(u[i] > x[i]){
      u[i] <- runif(1,umin,umax)
      x[i] <- rnorm(1,mu,sigma)
    }
  }
  group <- rbinom(n,1,prob)
  return(  data.frame('truncation'=u,'times'=x,'group'=group)  )
}

sim.mix <- function(n,prob=0.5,lambda=1,min.x=1,max.x=2,p=0.5,umin=0,umax=2){
  
  u <- runif(n,umin,umax)
  ind <- rbinom(n,1,p)
  x <- ind*rexp(n,lambda)+(1-ind)*runif(n,min.x,max.x)
  
  for(i in 1:n){
    while(u[i] > x[i]){
      u[i] <- runif(1,umin,umax)
      ind <- rbinom(1,1,p)
      x[i] <- ind*rexp(1,lambda)+(1-ind)*runif(1,min.x,max.x)
    }
  }
  group <- rbinom(n,1,prob)
  return( data.frame('truncation'=u,'times'=x,'group'=group) ) 
}

norm.exp <- function(n,prob=0.5,mean=0,sd=1,lambda=1){
  x <- rnorm(n,mean,sd)
  u <- rexp(n,lambda)
  
  for (i in 1:n){
    while(u[i] > x[i]){
      x[i] <- rnorm(1,mean,sd)
      u[i] <- rexp(1,lambda)
    }
  }
  return( data.frame('truncation'=u,'times'=x,'group'=rbinom(n,1,prob)) )
}

wei.norm <- function(n,prob=0.5,shape=1,scale=1,mean=0,sd=1){
  u <- rnorm(n,mean,sd)
  x <- rweibull(n,shape,scale)
  
  for (i in 1:n){
    while(u[i] > x[i]){
      u[i] <- rnorm(1,mean,sd)
      x[i] <- rweibull(1,shape,scale)
    }
  }
  return( data.frame('truncation'=u,'times'=x,'group'=rbinom(n,1,prob)) )
}

wei.exp <- function(n,prob=0.5,shape=1,scale=1,lambda=1){
  x <- rweibull(n,shape,scale)
  u <- rexp(n,lambda)
  
  for (i in 1:n){
    while(u[i] > x[i]){
      x[i] <- rweibull(1,shape,scale)
      u[i] <- rexp(1,lambda)
    }
  }
  return( data.frame('truncation'=u,'times'=x,'group'=rbinom(n,1,prob)) )
}

sim.norm2 <- function(n,prob=0.5,mu=0,sigma=1,mu.u,sigma.u){
  u <- rnorm(n,mu.u,sigma.u)
  x <- rnorm(n,mu,sigma)
  
  for (i in 1:n){
    while(u[i] > x[i]){
      u[i] <- rnorm(1,mu.u,sigma.u)
      x[i] <- rnorm(1,mu,sigma)
    }
  }
  return(  data.frame('truncation'=u,'times'=x,'group'=rbinom(n,1,prob))  )
}

exp.exp <- function(n,prob=0.5,lambda.x=1,lambda.u=1){
  u <- rexp(n,lambda.u)
  x <- rexp(n,lambda.x)
  
  for (i in 1:n){
    while(u[i] > x[i]){
      u[i] <- rexp(1,lambda.u)
      x[i] <- rexp(1,lambda.x)
    }
  }
  return(  data.frame('truncation'=u,'times'=x,'group'=rbinom(n,1,prob))  )
}


#### p-value for Kolmogorov-Smirnov ####
# sample = (truncation, t.to event)
ks.pv <- function(s0.0,s1.0,B=200,rep=F){
  
  n <- nrow(s0.0); m <- nrow(s1.0)

  if (rep == T){
    stat0 <- ks.stat(s0.0,s1.0,rep=T)$statistic
  } else {
    stat0 <- ks.stat(s0.0,s1.0,rep=F)$statistic
  }
  
  est0 <- lb.estimator.mod(s0.0)
  est1 <- lb.estimator.mod(s1.0)
  
  ttevent <- c(est0$fail.time,est1$fail.time)
  prob.ev <- c(n*est0$prob/(n+m), m*est1$prob/(n+m))
  
  est.trunc0 <- lb.estimator( -s0.0[,2:1] )
  est.trunc1 <- lb.estimator( -s1.0[,2:1] )
  
  prob0.trunc <- est.trunc0$prob
  prob1.trunc <- est.trunc1$prob
  
  trunc0 <- -est.trunc0$fail.time
  trunc1 <- -est.trunc1$fail.time
  
  ## bootstrap
  stat <- numeric(B)
  for (b in 1:B){ 
    sample0 <- data.frame( sample(trunc0,size=n,replace=T,prob=prob0.trunc),
                           sample(ttevent,size=n,replace=T,prob=prob.ev))
    
    sample1 <- data.frame( sample(trunc1,size=m,replace=T,prob=prob1.trunc),
                           sample(ttevent,size=m,replace=T,prob=prob.ev))
    
    for (i in 1:n){
      while(sample0[i,1] > sample0[i,2]){
        sample0[i,1] <- sample(trunc0,size=1,prob=prob0.trunc)
        sample0[i,2] <- sample(ttevent,size=1,prob=prob.ev)
      }
    }
    
    for (i in 1:m){
      while(sample1[i,1] > sample1[i,2]){
        sample1[i,1] <- sample(trunc1,size=1,prob=prob1.trunc)
        sample1[i,2] <- sample(ttevent,size=1,prob=prob.ev)
      }
    }

    stat[b] <- ks.stat(sample0,sample1,rep=F)$statistic
    
  } # fin del bootstrap 
  
  # p-valor 
  pv <- mean(stat >= stat0)
  
  return( list( 'statistic'= stat0,'p.value'= pv, 'boot.stat'= stat) )
}


#### log-rank ####
# muestra = (truncation,t. to event)
log.rank <- function(p.sample,correction=F){
  
  times <- unique(p.sample[,2])
  groups <- sort(unique(p.sample[,3]))
  K <- length(groups)
  D <- length(times)
  d <- n <- matrix(numeric(K*D),nrow=D)
  
  for ( i in 1:K ){
    ind <- groups[i]
    for (j in 1:D){
      d[j,i] <- sum( p.sample[p.sample[,3]==ind,2] == times[j] )
      n[j,i] <- sum( p.sample[p.sample[,3]==ind,1] <= times[j] & times[j] <= p.sample[p.sample[,3]==ind,2] )
    }
  }
  
  d.pooled <- rowSums(d) # # failures in the pooled sample
  n.pooled <- rowSums(n) # # individuals at risk in the pooled sample
  
  Z <- numeric(K-1)
  varcov <- matrix(numeric((K-1)^2),ncol=K-1)
  
  for (i in 1:(K-1)){
    Z[i] <- sum( d[,i] - d.pooled*n[,i]/n.pooled )
    diag(varcov)[i] <- sum( n[,i]*d.pooled*( 1- n[,i]/n.pooled )/n.pooled  )
    if (i < (K-1)){
      for ( j in (i+1):(K-1) ){
        varcov[i,j] <- - sum( d.pooled*n[,i]*n[,j]/n.pooled^2  )
        varcov[j,i] <- varcov[i,j]
      }
      
    }
  }
  
  # corrector <- 1
  # if (correction==T){
  #   corrector <- (ni-di)/(ni-1) # va a valer 1 siempre que no haya empates
  #   corrector[which(is.na(corrector))] <- 1
  # }
  
  statistic <- Z %*% solve(varcov) %*% Z
  
  pvalue <- 1-pchisq(q=statistic,df=(K-1))
  return(list('statistic'=statistic,'p.value'=pvalue))
  
}
#### test for independence (from Tsai(1990)) ####
# sample = (truncation times, event times)
# H0: U y X independientes
indep <- function(sample){
  n <- nrow(sample)
  r <- s <- numeric(n)
  
  for (i in 1:n){
    risk.i <- which(  sample[,1] <= sample[i,2] & sample[i,2] <= sample[,2]  )
    s[i] <- sum( sign(sample[risk.i,1] - sample[i,1])  )
    r[i] <- length( risk.i  )
  }
  
  stat <- sum(s)/( sqrt( sum(r^2-1)/3 )  )
  
  pv <- 2*(1-pnorm(abs(stat)))
  
  return(list('statistic'=stat,'p.value'=pv))
}

#### linear transformation of the truncation times (Chiou et al. (2019))####
# data = (left truncation, times to event, ... )
trans.lin <- function(data,plots=F,cor=F){
  if (cor==T){
    a <- -cor(data[,1],data[,2])*sd(data[,1])/sd(data[,2])
  } else {
    library(tranSurv)
    a <- trSurvfit(data[,1],data[,2],plots=plots)$byP$par
  }
  data[,1] <- (data[,1] + a* data[,2]) / ( 1 + a)
  return(data)
}

#### left truncation and right censorship ####
ltrc.estimator <- function(sample){
  
  dist.fail <- unique(sample[,2])
  di <- ni <- prob <- numeric( length(dist.fail) )
  j <- 1
  
  for (t in sort(dist.fail)) {
    di[j] <- sum( sample[,2] == t  & sample[,3] == 1)
    ni[j] <- sum( sample[,1] <= t & t <= sample[,2]  )
    j <- j+1
  }
  
  cum.prod <- c(1,cumprod( (1-di/ni)^(sample[,3]) ))
  
  for (i in 1:length(dist.fail) ){
    prob[i] <- cum.prod[i]-cum.prod[i+1]
  }
  
  lista <- list('fail.time'= sort( dist.fail ) , 'estimation'= cumprod(1-di/ni),
                'prob'=prob)
  return(lista)
}

simulacion.ltrc <- function(n,prob=0.5,lambda=1,beta=0,umin=0,umax=1,c.min=1,c.max=2.5){
  lt.sample <- simulacion.x(n,prob=0.5,lambda=1,beta=0,umin=0,umax=1)
  censoring.times <- runif(n,c.min,c.max)
  sample <- data.frame( 'truncation'=lt.sample[,1], 'observed time'=min(lt.sample[,2],censoring.times), 
                        'indicator'=ifelse(lt.sample[,2] <= censoring.times,1,0),'group'=lt.sample[,3] )
  return(sample)
}
#### statistic based on density estimators ####
depa <- function(x){ 0.75*(1-x^2)*(abs(x) <= 1)}

fun_integrate <- function(x, t.evento0, probabilities0, t.evento, probabilities,h.sub,h) { 
  ( sum(probabilities0*dnorm(x,t.evento0,h.sub)) - sum(probabilities*dnorm(x,t.evento,h)) )^2
}

library(DTDA)

stat.int <- function(sample,H){
  
  ind <- unique(sample[,3])
  K <- length(ind)
  aux <- numeric(K)
  
  est <- lb.estimator(sample)
  n <- length(est$fail.time)
  
  for (j in 1:K){
    sub.s <- sample[sample[,3]==ind[j],]
    nsub <- nrow(sub.s)
    est0 <- lb.estimator(sub.s)
    aux[j] <- nsub/n * integrate(Vectorize(fun_integrate,vectorize.args='x'), lower = -Inf, upper = Inf, t.evento0 = est0$fail.time, 
                  probabilities0 =  est0$prob, t.evento = est$fail.time, probabilities =  est$prob, H[j], H[K+1],rel.tol=1e-6,subdivisions=2000)$value
  }
  return('statistic'= sum(aux))
  
}

stat.int.mod <- function(sample,H){
  
  ind <- unique(sample[,3])
  K <- length(ind)
  aux <- numeric(K)
  
  est <- lb.estimator.mod(sample)
  n <- length(est$fail.time)
  
  for (j in 1:K){
    sub.s <- sample[sample[,3]==ind[j],]
    nsub <- nrow(sub.s)
    est0 <- lb.estimator.mod(sub.s)
    aux[j] <- nsub/n * integrate(Vectorize(fun_integrate,vectorize.args='x'), lower = -Inf, upper = Inf, t.evento0 = est0$fail.time, 
                                 probabilities0 =  est0$prob, t.evento = est$fail.time, probabilities =  est$prob, H[j], H[K+1],
                                 rel.tol=1e-6,subdivisions=2000)$value
  }
  return('statistic'= sum(aux))
  
}


stat.sum <- function(sample,H){
  
  groups <- unique( sort(sample[,3]) ) # indicators of groups
  K <- length(groups) # number of groups
  
  # if ( length(H) != (K+1) ){ cat('Bandwidth(s) missing') }
  
  est <- lb.estimator(sample) # Lynden-Bell estimator of the pooled sample
  aux <- numeric(K) # vector of the sum statistic in each group
  
  for (j in 1:K){
    ind <- groups[j]
    est.sub <- lb.estimator( sample[sample[,3]==ind,] )
    aux[j] <- length(which(sample[,3]==ind)) * 
      sum( est$prob * ( colSums( dnorm( outer(est.sub$fail.time,est$fail.time,'-')/  H[j]   )* est.sub$prob  )/H[j] -
                          colSums( dnorm( outer(est$fail.time,  est$fail.time,'-')  / H[K+1] ) *  est$prob    )/H[K+1] )^2 )
  }
  return( sum(aux)/nrow(sample) )
}

fun_integrate3 <- function(x, t.evento1, probabilities1, t.evento2, probabilities2,t.evento3, probabilities3,h1,h2,h3,p1,p2,p3) { 
  ( (1-p1)*sum(probabilities1*dnorm(x,t.evento1,h1)) - ( p2*sum(probabilities2*dnorm(x,t.evento2,h2))+p3*sum(probabilities3*dnorm(x,t.evento3,h3)) ) )^2
}

stat.int3 <- function(sample,H,weights){
  
  ind <- unique(sample[,3])
  K <- length(ind) # es 3
  # aux <- numeric(K)
  N <- numeric(K)
  
  # para 3 muestras
  
  #M1
  sub.s <- sample[sample[,3]==ind[1],]
  est <- lb.estimator.mod(sub.s)
  omega1 <- est$prob
  x1 <- est$fail.time
  N[1] <- length(x1)
  
  #M2
  sub.s <- sample[sample[,3]==ind[2],]
  est <- lb.estimator.mod(sub.s)
  omega2 <- est$prob
  x2 <- est$fail.time
  N[2] <- length(x2)
  
  #M3
  sub.s <- sample[sample[,3]==ind[3],]
  est <- lb.estimator.mod(sub.s)
  omega3 <- est$prob
  x3 <- est$fail.time
  N[3] <- length(x3)
  
  n <- sum(N)
  
  # p <- numeric(3)
  # if (weights==NULL){
  #   p[1] <- n1/n
  #   p[2] <- n2/n
  #   p[3] <- n3/n
  # } else{
  #   p <- weights
  # }
  
  p <- weights
  
  h1 <- H[1]; h2 <- H[2]; h3 <- H[3]
  p1 <- p[1]; p2 <- p[2]; p3 <- p[3]
  
  # para primera muestra (j=1)
  
  statistic <- N[1]/n * integrate(Vectorize(fun_integrate3,vectorize.args='x'), lower = -Inf, upper = Inf, x1, omega1, 
                                  x2, omega2,x3, omega3,h1,h2,h3,p1,p2,p3,subdivisions=2000,rel.tol=1e-6)$value+
    N[2]/n * integrate(Vectorize(fun_integrate3,vectorize.args='x'), lower = -Inf, upper = Inf, x2, omega2, 
                       x1, omega1,x3, omega3,h2,h1,h3,p2,p1,p3,subdivisions=2000,rel.tol=1e-6)$value+
    N[3]/n * integrate(Vectorize(fun_integrate3,vectorize.args='x'), lower = -Inf, upper = Inf, x3, omega3, 
                       x2, omega2,x1, omega1,h3,h2,h1,p3,p2,p1,subdivisions=2000,rel.tol=1e-6)$value
  
  return('statistic'= statistic)
  
}

fun_integrate2 <- function(x, t.evento1, probabilities1, t.evento2, probabilities2,h1,h2,p1,p2) { 
  ( (1-p1)*sum(probabilities1*dnorm(x,t.evento1,h1)) - ( p2*sum(probabilities2*dnorm(x,t.evento2,h2)) ) )^2
}

stat.int2 <- function(sample,H,weights){
  
  ind <- unique(sample[,3])
  K <- length(ind) # es 3
  # aux <- numeric(K)
  N <- numeric(K)
  
  # para 3 muestras
  
  #M1
  sub.s <- sample[sample[,3]==ind[1],]
  est <- lb.estimator.mod(sub.s)
  omega1 <- est$prob
  x1 <- est$fail.time
  N[1] <- length(x1)
  
  #M2
  sub.s <- sample[sample[,3]==ind[2],]
  est <- lb.estimator.mod(sub.s)
  omega2 <- est$prob
  x2 <- est$fail.time
  N[2] <- length(x2)
  
  n <- sum(N)
  
  p <- weights
  
  h1 <- H[1]; h2 <- H[2]
  p1 <- p[1]; p2 <- p[2]
  
  # para primera muestra (j=1)
  
  statistic <- N[1]/n * integrate(Vectorize(fun_integrate2,vectorize.args='x'), lower = -Inf, upper = Inf, x1, omega1, 
                                  x2, omega2,h1,h2,p1,p2,rel.tol=1e-6,subdivisions=2000)$value+
    N[2]/n * integrate(Vectorize(fun_integrate2,vectorize.args='x'), lower = -Inf, upper = Inf, x2, omega2, 
                       x1, omega1,h2,h1,p2,p1,rel.tol=1e-6,subdivisions=2000)$value
  
  return('statistic'= statistic)
  
}

representacion <- function(sample,h,t){
  est <- lb.estimator(sample)
  M <- matrix(numeric(length(t)*length(est$fail.time)),
              ncol=length(t),nrow=length(est$fail.time))
  for (i in 1:length(t)){
    M[,i] <- dnorm( t[i],est$fail.time,h )
  }
  return(colSums( M*est$prob ))
  
}


#### estimation ####
est.moments <- function(sample){
	est <- lb.estimator(sample)
	mu <- sum( est$fail.time * est$prob  )
	sigma2 <- sum( est$prob * (est$fail.time-mu)^2 ) 
	return( list('mean'=mu,'sigma2'=sigma2) )
}

gamma.hat <- function(sample){
  
  est <- lb.estimator.mod(-sample[,2:1])
  failures <- rev(-est$fail.time) # ttevent U
  estimation <- rev(est$estimation) # estimation U
  
  t <- unique( sample[,2] ) ## ttevent X
  aux <- numeric( length(t) ) # vector for the 1/G(Xi)
  
  for (i in 1:length(aux) ){
    
    # big values
    if ( t[i] >= max(failures) ){aux[i] <- 1} 
    
    for ( j in 1:(length(failures)-1) ){
      if ( (failures[j] <= t[i] & t[i] < failures[j+1]) == T){
        aux[i] <- 1/estimation[j]
      }
    }
    
  }
  
  est.int <- sum( lb.estimator.mod(sample)$prob*aux )
  
  return( data.frame('gamma'=1/mean(aux),'est.int'=est.int) )
  
}

hrot <- function(sample){
  c <- gamma.hat(sample)
  n <- nrow(sample)
  hnr <- (0.8*c$gamma*c$est.int)^(0.2)*sqrt(est.moments(sample)$sigma2)*n^(-0.2)
  return(hnr)
}
