  ### R functions ###
  
  rm(list=ls())
  
  # path <- ("C:/Users/T16/Desktop/Code") 
  path <- ("R/Replicability Assessment") #document directory, please check it before you run this code!
  dyn.load(paste(path,"/ConcordanceX64.dll",sep="")) 
  
  #modified from package mclust
  unmap <- function(classification){
    n <- length(classification)
    u <- sort(unique(classification))
    labs <- as.character(u)
    k <- length(u)
    z <- matrix(0, n, k)
    for (j in 1:k) z[classification == u[j], j] <- 1
    dimnames(z) <- list(NULL, labs)
    return(z)
  }
  
  em.normal.partial.concordant <- function(data, class, tol=0.000001, restriction=0, constrain=0, iteration=1000){
    n <- as.integer(dim(data)[1])
    g <- as.integer(nlevels(as.factor(class)))
    
    yl.outer <- function(k, zx, zy){
      return( c(zx[k,] %o% zy[k,]) )
    }
    yl.diag <- function(k, z){
      return( c(diag(z[k,])) )
    }
    
    zx <- unmap(class[,1])
    zy <- unmap(class[,2])
    zxy <- sapply(1:dim(zx)[1], yl.outer, zx, zy)
    
    pi <- double(g*g)
    mu <- double(g)
    sigma <- double(g)
    nu <- double(g)
    tau <- double(g)
    loglik <- double(1)
    convergence <- integer(1)
    if(restriction>0){
      if(length(constrain)==g){
        results <- .C("em_normal_partial_concordant", as.double(data[,1]), as.double(data[,2]), as.double(t(zxy)), n, pi, mu, sigma, nu, tau, g, loglik, as.double(tol), as.integer(restriction), as.integer(constrain), as.integer(iteration), convergence)
      }else{
        print("Error with constrain!")
        return(0)
      }
    }else{
      results <- .C("em_normal_partial_concordant", as.double(data[,1]), as.double(data[,2]), as.double(t(zxy)), n, pi, mu, sigma, nu, tau, g, loglik, as.double(tol), as.integer(restriction), as.integer(constrain), as.integer(iteration), convergence)
    }
    print(paste("convergence within", results[[15]], "run?", results[[16]], sep=" "))
    return(list(model="PCD", convergence=results[[16]], pi=t(array(results[[5]],dim=c(g,g))), mu_sigma=rbind(results[[6]], results[[7]]), nu_tau=rbind(results[[8]], results[[9]]), loglik=results[[11]], class=apply(array(results[[3]], dim=c(n,g*g)),1,order,decreasing=T)[1,], z=array(results[[3]], dim=c(n,g*g))))
  }

  #Example 1: Simulation data  
  Pd1 <- c(rnorm(n = 200,mean = 3,sd = 1),rnorm(n = 200,mean = -3,sd = 1),rnorm(n = 600,mean = 0,sd = 1))
  #200*N(3,1)+200*N(-3,1)+600*N(0,1)
  Pd2 <- c(rnorm(n = 200,mean = 3,sd = 1),rnorm(n = 200,mean = -3,sd = 1),rnorm(n = 600,mean = 0,sd = 1))
  #200*N(3,1)+200*N(-3,1)+600*N(0,1)
  pd <- cbind(Pd1,Pd2)
  P_IR <- c()
  R_IR <- c()

  class <- cbind(rep(0, dim(pd)[1]), rep(0, dim(pd)[1]))
  class[pd[,1]>0+1.96,1] <- 2
  class[pd[,1]<0-1.96,1] <- 1
  class[pd[,2]>0+1.96,2] <- 2
  class[pd[,2]<0-1.96,2] <- 1
  
  pd_aal2 <- em.normal.partial.concordant(pd, class, tol=0.001, restriction=1, constrain=c(0,-1,1), iteration=5000)
  P <- pd_aal2[["pi"]] 
  P_IR<- (P[1,2]+P[1,3]+P[2,1]+P[2,3]+P[3,1]+P[3,2])/(1-P[1,1]) 
  
  R <- pd_aal2[["z"]] 
  R_IR <- R[,1]+R[,2]+R[,3]+R[,4]+R[,6]+R[,7]+R[,8]
  
  
  #Example 2: Paired z-scores for Proteomics-based association study

  load(file.path("data/Z_Example.rda"))
  P_IR <- c()
  R_IR <- c()
  
  class <- cbind(rep(0, dim(pd)[1]), rep(0, dim(pd)[1]))
  class[pd[,1]>0+1.96,1] <- 2
  class[pd[,1]<0-1.96,1] <- 1
  class[pd[,2]>0+1.96,2] <- 2
  class[pd[,2]<0-1.96,2] <- 1
  
  pd_aal2 <- em.normal.partial.concordant(pd, class, tol=0.001, restriction=1, constrain=c(0,-1,1), iteration=5000)
  P <- pd_aal2[["pi"]] 
  P_IR<- (P[1,2]+P[1,3]+P[2,1]+P[2,3]+P[3,1]+P[3,2])/(1-P[1,1]) 
  
  R <- pd_aal2[["z"]] 
  R_IR <- R[,1]+R[,2]+R[,3]+R[,4]+R[,6]+R[,7]+R[,8]
