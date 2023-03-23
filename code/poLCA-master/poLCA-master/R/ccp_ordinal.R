class.cond.prob <- function(params, K.j, J, R){
  
  alpha <- params[1:sum(K.j-1)] # a vector of \sum(K_j-1)
  tau <- params[-c(1:sum(K.j-1))] # a vector of J*(R-1) 
  ind.alpha<- rep(1:J, K.j-1)
  ind.tau <- rep(1:J,each=R-1)
  probs.ccp <- list()
  for(j in 1:ncol(K.j)){
    alpha.j <- alpha[ind.alpha==j]
    tau.j <- tau[ind.tau==j]
    K_j=K.j[j]
    ccp <- matrix(NA, nrow=R, ncol = K_j)
    
    for(r in 1:R){
      if(r==R){
        alpha.sum <- rev(cumsum(rev(alpha.j)))
      }else{
        alpha.sum <- rev(cumsum(rev(alpha.j))) + rev(1:(K_j-1))*tau.j[r]  
      }
      pi.jr <- exp(c(alpha.sum,0))
      pi.jr <- pi.jr/sum(pi.jr)
      ccp[r, ] <- pi.jr
    }  
    probs.ccp[[j]] <- ccp
  }
  return(probs.ccp)
}


grad.act <- function(params, rgivy,w, y, K.j, J , R, N){
  probs.test <- class.cond.prob(params, K.j, J, R)
  grad.act <- rep(NA, sum(K.j-1)+J*(R-1))
  count=1
  for(j in 1:J){
    for(m in 1:(K.j[j]-1)){
      grad.alpha=0
      for(r in 1:R){
          for(l in 1:m){
            grad.alpha=grad.alpha+ sum(w*rgivy[,r]*((y[,j]==l)-probs.test[[j]][r,l]))
          }
      }
      grad.act[count] <- grad.alpha
      count=count+1
    }
  }
  for(j in 1:J){
    for(r in 1:(R-1)){
      grad.tau=0
        for(k in 1:K.j[j]){
          grad.tau=grad.tau+ sum(w*rgivy[,r]*(K.j[j]-k)*((y[,j]==k)-probs.test[[j]][r,k]))
      }
      grad.act[count] <- grad.tau
      count=count+1
    }
  }
  return(grad.act)
}



jackknife.se <- function(formula,data,w,nclass, maxiter, 
                         tol, lt,nrep=nrep,seed=seed){
  
  strata_cluster <- table(data$strata,data$cluster)
  # Remove entries with only one cluster in a strata
  inds <- sum(unlist(apply(strata_cluster, 1, function(x){if(sum(x==0)!=1){sum(x!=0)}})))
  
  mod.probs <- vector('list', inds)
  mod.P <- vector('list', inds)
  mod.probs.unvectorize <- vector('list', inds)
  c <- 1
  jk.w <- rep(NA, inds)
  for(i in 1:nrow(strata_cluster)){
  #for(i in 1:5){
    k_i <- sum(strata_cluster[i,]!=0)
    for(j in 1:ncol(strata_cluster)){
      if(k_i==1 | strata_cluster[i,j]==0 ){
        next
        print("skipped")
      }else{
        if(c %% 1 ==0) {print(paste("PSU Number:" ,c))}
        obs_removed_ind <- which(data$strata==as.numeric(rownames(strata_cluster))[i] & 
                                   data$cluster == as.numeric(colnames(strata_cluster))[j])
        obs_reweighted <- which(data$strata==as.numeric(rownames(strata_cluster))[i] & 
                                  data$cluster != as.numeric(colnames(strata_cluster))[j])
        jk.w[c] <- (1+length(obs_removed_ind)/length(obs_reweighted)) # 1+sum(w[obs_removed_ind])/sum(w[obs_reweighted])
        w_rw <- w
        w_rw[obs_reweighted] <- w[obs_reweighted]*jk.w[c] #*jk.w[c]/(jk.w[c]-1)
        
        data_sc <- data[-obs_removed_ind, ]
        data_sc <- data.frame(data_sc)
        w_sc <- w_rw[-obs_removed_ind]
        #w_sc <- w_sc/sum(w_sc)*length(w_sc)
        ff <- as.formula(formula)
        environment(ff) <- environment()
        mod <- poLCA.ordinal(formula=ff, data = data_sc,weights = w_sc, nclass = nclass, 
                             nrep = nrep, calc.se = FALSE,
                             maxiter = maxiter, tol = tol,lt=lt,
                             verbose = FALSE,seed = seed)
        mod.probs.unvectorize[[c]] <- mod$probs
        mod.probs[[c]] <- poLCA.vectorize(mod$probs)$vecprobs
        mod.P[[c]] <- mod$P
        
        print(head(mod.probs[[c]]))
        print(mod.P[[c]])
        
        c <- c+1
      }
    }
  }
  # calculate estimator and se from above objects for probs and P
  mod.P.mat <- do.call(rbind, mod.P)
  mod.probs.mat <- do.call(rbind, mod.probs)
  
  jk.P.mean <- colMeans(mod.P.mat)
  jk.probs.mean <- colMeans(mod.probs.mat)
  jk.P.sqe <- t(apply(mod.P.mat, 1, function(x){(x-jk.P.mean)^2}))
  jk.probs.sqe <- t(apply(mod.probs.mat, 1, function(x){(x-jk.probs.mean)^2}))
  jk.P.se <- sqrt(matrix(1/jk.w,nrow=1) %*% jk.P.sqe)
  jk.probs.se <- sqrt(matrix(1/jk.w,nrow=1) %*% jk.probs.sqe)
  
  vpse <- list("vecprobs"=as.vector(jk.probs.se),
               "numChoices"=poLCA.vectorize(mod$probs)$numChoices,
               "classes"=poLCA.vectorize(mod$probs)$classes)
  jk.probs.se <- poLCA.unvectorize(vpse)
  
  vpm <- list("vecprobs"=as.vector(jk.probs.mean),
               "numChoices"=poLCA.vectorize(mod$probs)$numChoices,
               "classes"=poLCA.vectorize(mod$probs)$classes)
  jk.probs.mean <- poLCA.unvectorize(vpm)
  
  res <- list("P.jk"=jk.P.mean,
              "probs.jk"=jk.probs.mean,
              "P"=jk.P.se,
              "probs"=jk.probs.se)
  return(res)
}
