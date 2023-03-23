library(MASS)

poLCA.ordinal <-
function(formula,data,weights,nclass=2,maxiter=1000,graphs=FALSE,tol=1e-10,
                na.rm=TRUE,probs.start=NULL,nrep=1,verbose=TRUE,
         calc.se=FALSE,lt=0.001,ltbeta=0.1,seed=1234) {
    set.seed(seed)
    starttime <- Sys.time()
    
    if(calc.se & !(sum(c("strata","cluster") %in% colnames(data))==2)){
      stop("Jackknife se requires strata and cluster info")
    }
    # from lm
    mfall <- match.call(expand.dots = FALSE)
    m <- match(
      c(
        "formula", "data", "subset", "weights", "na.action",
        "offset"
      ),
      names(mfall),
      0L
    )
    mframe <- mfall[c(1L, m)]
    mframe$drop.unused.levels <- TRUE
    mframe[[1L]] <- quote(stats::model.frame)
    mframe <- eval(mframe, parent.frame())
    #mframe <- model.frame(formula,data,na.action=NULL,weights=weights)
    mf <- model.response(mframe)
    if (any(mf<1,na.rm=TRUE) | any(round(mf) != mf,na.rm=TRUE)) {
        cat("\n ALERT: some manifest variables contain values that are not
    positive integers. For poLCA to run, please recode categorical
    outcome variables to increment from 1 to the maximum number of
    outcome categories for each variable. \n\n")
        ret <- NULL
    } else {
    data <- data[rowSums(is.na(model.matrix(formula,mframe)))==0,]
    if (na.rm) {
      # from lm
      mfall <- match.call(expand.dots = FALSE)
      m <- match(
        c(
          "formula", "data", "subset", "weights", "na.action",
          "offset"
        ),
        names(mfall),
        0L
      )
      mframe <- mfall[c(1L, m)]
      mframe$drop.unused.levels <- TRUE
      mframe[[1L]] <- quote(stats::model.frame)
      mframe <- eval(mframe, parent.frame())
    
      y <- model.response(mframe)
      w <- model.weights(mframe)
      if(is.null(w)){
        w <- rep(1,nrow(y))
      }
    } else {
      # from lm
      mfall <- match.call(expand.dots = FALSE)
      m <- match(
        c(
          "formula", "data", "subset", "weights", "na.action",
          "offset"
        ),
        names(mfall),
        0L
      )
      mframe <- mfall[c(1L, m)]
      mframe$drop.unused.levels <- TRUE
      mframe[[1L]] <- quote(stats::model.frame)
      mframe <- eval(mframe, parent.frame())
      
      #mframe <- model.frame(formula,data,na.action=NULL)
      y <- model.response(mframe)
      w <- model.weights(mframe)
      y[is.na(y)] <- 0
      if(is.null(w)){
        w <- rep(1,nrow(y))
        }
        
    }
    if (any(sapply(lapply(as.data.frame(y),table),length)==1)) {
        y <- y[,!(sapply(apply(y,2,table),length)==1)]
        cat("\n ALERT: at least one manifest variable contained only one
    outcome category, and has been removed from the analysis. \n\n")
    }
    x <- model.matrix(formula,mframe)
    N <- nrow(y)
    J <- ncol(y)
    K.j <- t(matrix(apply(y,2,max)))
    R <- nclass
    S <- ncol(x)
    if (S>1) { calc.se <- TRUE }
    eflag <- FALSE
    probs.start.ok <- TRUE
    ret <- list()
    if (R==1) {
        ret$probs <- list()
        for (j in 1:J) {
            ret$probs[[j]] <- matrix(NA,nrow=1,ncol=K.j[j])
            for (k in 1:K.j[j]) { ret$probs[[j]][k] <- sum((y[,j]==k)*w)/sum((y[,j]>0)*w) }
        }
        ret$probs.start <- ret$probs
        ret$P <- 1
        ret$posterior <- ret$predclass <- prior <- matrix(1,nrow=N,ncol=1)
        ret$llik <- sum(w*log(poLCA.ylik.C(poLCA.vectorize(ret$probs),y)))
        if (calc.se) {
            se <- jackknife.se(formula,data,weights, nclass, maxiter=maxiter, tol=tol, lt=lt,nrep=nrep)#poLCA.se(y,x,ret$probs,prior,ret$posterior)
            ret$probs.se <- se$probs           # standard errors of class-conditional response probabilities
            ret$P.se <- se$P                   # standard errors of class population shares
        } else {
            ret$probs.se <- NA
            ret$P.se <- NA
        }
        ret$numiter <- 1
        ret$probs.start.ok <- TRUE
        ret$coeff <- NA
        ret$coeff.se <- NA
        ret$coeff.V <- NA
        ret$eflag <- FALSE
        if (S>1) {
            cat("\n ALERT: covariates not allowed when nclass=1; will be ignored. \n \n")
            S <- 1
        }
    } else {
        if (!is.null(probs.start)) { # error checking on user-inputted probs.start
            if ((length(probs.start) != J) | (!is.list(probs.start))) {
                probs.start.ok <- FALSE
            } else {
                if (sum(sapply(probs.start,dim)[1,]==R) != J) probs.start.ok <- FALSE
                if (sum(sapply(probs.start,dim)[2,]==K.j) != J) probs.start.ok <- FALSE
                if (sum(round(sapply(probs.start,rowSums),4)==1) != (R*J)) probs.start.ok <- FALSE
            }
        }
        ret$llik <- -Inf
        ret$attempts <- NULL
        for (repl in 1:nrep) { 
          # automatically reestimate the model multiple times to locate the global max llik
            error <- TRUE; firstrun <- TRUE
            probs <- probs.init <- probs.start
            while (error) { # error trap
                error <- FALSE
                b <- rep(0,S*(R-1))
                prior <- poLCA.updatePrior(b,x,R)
                if ((!probs.start.ok) | (is.null(probs.start)) | (!firstrun) | (repl>1)) { # only use the specified probs.start in the first nrep
                    probs <- list()
                    # here
                    params <- c(rnorm(sum(K.j-1), 0, 2),rnorm(J*(R-1),0,1))
                    probs <- class.cond.prob(params, K.j=K.j, J=J, R=R)

                    probs.init <- probs
                }
                vp <- poLCA.vectorize(probs)
                iter <- 1
                llik <- matrix(NA,nrow=maxiter,ncol=1)
                llik[iter] <- -Inf
                dll <- Inf
                while ((iter <= maxiter) & (dll > tol) & (!error)) {
                    iter <- iter+1
                    rgivy <- poLCA.postClass.C(prior,vp,y) # calculate posterior
                    
                    # update probs
                    params <- params + lt*grad.act(params,rgivy,w=w, y, K.j=K.j, J=J , R=R , N=N)
                    probs.new <- class.cond.prob(params, K.j=K.j, J=J , R=R)
                    vp$vecprobs <- poLCA.vectorize(probs.new)$vecprobs
                    if (S>1) {
                        dd <- poLCA.dLL2dBeta.C(rgivy,prior,x)
                        b <- b + ltbeta*ginv(-dd$hess) %*% dd$grad     # update betas
                        prior <- poLCA.updatePrior(b,x,R)       # update prior
                    } else {
                        prior <- matrix(w %*% rgivy/sum(w),nrow=N,ncol=R,byrow=TRUE)
                    }
                    llik[iter] <- sum(w*log(rowSums(prior*poLCA.ylik.C(vp,y))))
                    dll <- llik[iter]-llik[iter-1]
                    # if(dll<tol){
                    #   print("gradient optimized")
                    # }
                    if (is.na(dll)) {
                        error <- TRUE
                    } else if ((S>1) & (dll < -1e-7)) {
                        error <- TRUE
                    }
                }
                # if (!error) { 
                #     if (calc.se) {
                #         se <- jackknife.se(formula=formula,data=data,w=w, nclass=nclass, maxiter=maxiter, tol=tol, lt=lt,nrep=nrep)#poLCA.se(y,x,poLCA.unvectorize(vp),prior,rgivy)
                #     } else {
                #         se <- list(probs=NA,P=NA,probs.jk=NA,P.jk=NA)
                #     }
                # } else {
                #     eflag <- TRUE
                # }
                if(error){eflag <- TRUE}
                firstrun <- FALSE
            } 
            # finish estimating model without triggering error
            ret$attempts <- c(ret$attempts,llik[iter])
            if (llik[iter] > ret$llik) {
                ret$posterior <- rgivy           # NxR matrix of posterior class membership probabilities
                ret$P <- colSums(diag(w)%*% ret$posterior)/sum(w) # estimated class population shares
                class_order <- order(ret$P)
                ret$P <- ret$P[class_order] # reorder P
                ret$posterior <- ret$posterior[ ,class_order] # reorder posterior
                ret$llik <- llik[iter]      # maximum value of the log-likelihood
                ret$probs.start <- lapply(probs.init,function(x) x[class_order, ])  #reorder probs.init    # starting values of class-conditional response probabilities
                ret$probs <- poLCA.unvectorize(vp) # estimated class-conditional response probabilities
                ret$probs <- lapply(ret$probs,function(x) x[class_order, ]) # reorder probs
                ret$predclass <- apply(ret$posterior,1,which.max)   # Nx1 vector of predicted class memberships, by modal assignment
                ret$numiter <- iter-1              # number of iterations until reaching convergence
                ret$probs.start.ok <- probs.start.ok # if starting probs specified, logical indicating proper entry format
                # ret$params <- params # this
                if (S>1) {
                    b <- matrix(b,nrow=S)
                    rownames(b) <- colnames(x)
                    rownames(se$b) <- colnames(x)
                    ret$coeff <- b                 # coefficient estimates (when estimated)
                    ret$coeff.se <- se$b           # standard errors of coefficient estimates (when estimated)
                    ret$coeff.V <- se$var.b        # covariance matrix of coefficient estimates (when estimated)
                } else {
                    ret$coeff <- NA
                    ret$coeff.se <- NA
                    ret$coeff.V <- NA
                }
                ret$eflag <- eflag                 # error flag, true if estimation algorithm ever needed to restart with new initial values
            }
            if (nrep>1 & verbose) { cat("Model ",repl,": llik = ",llik[iter]," ... best llik = ",ret$llik,"\n",sep=""); flush.console() }
        } # end replication loop
        
        if (calc.se) {
            se <- jackknife.se(formula=formula,data=data,w=w, 
                               nclass=nclass, maxiter=maxiter, tol=tol, 
                               lt=lt,nrep=nrep,seed=seed)#poLCA.se(y,x,poLCA.unvectorize(vp),prior,rgivy)
        } else {
            se <- list(probs=NA,P=NA,probs.jk=NA,P.jk=NA)
        }
        
        ret$probs.se <- se$probs           # standard errors of class-conditional response probabilities
        ret$P.se <- se$P                   # standard errors of class population shares
        ret$P.jk <- se$P.jk
        ret$probs.jk <- se$probs.jk
        
    }
    names(ret$probs) <- colnames(y)
    if (calc.se) { 
      names(ret$probs.se) <- colnames(y)
      names(ret$probs.jk) <- colnames(y)
    }
    ret$npar <- (sum(K.j-1))+ J*(R-1) + (R-1)                  # number of degrees of freedom used by the model (number of estimated parameters)
    if (S>1) { ret$npar <- ret$npar + (S*(R-1)) - (R-1) }
    ret$aic <- (-2 * ret$llik) + (2 * ret$npar)         # Akaike Information Criterion
    ret$bic <- (-2 * ret$llik) + (log(N) * ret$npar)    # Schwarz-Bayesian Information Criterion
    ret$Nobs <- sum(rowSums(y==0)==0)                   # number of fully observed cases (if na.rm=F)
    # if (all(rowSums(y==0)>0)) { # if no rows are fully observed
    #     ret$Chisq <- NA
    #     ret$Gsq <- NA
    #     ret$predcell <- NA
    # } else {
    #     compy <- poLCA.compress(y[(rowSums(y==0)==0),])
    #     datacell <- compy$datamat
    #     rownames(datacell) <- NULL
    #     freq <- compy$freq
    #     if (!na.rm) {
    #         fit <- matrix(ret$Nobs * (poLCA.ylik.C(poLCA.vectorize(ret$probs),datacell) %*% ret$P))
    #         ret$Chisq <- sum((freq-fit)^2/fit) + (ret$Nobs-sum(fit)) # Pearson Chi-square goodness of fit statistic for fitted vs. observed multiway tables
    #     } else {
    #         fit <- matrix(N * (poLCA.ylik.C(poLCA.vectorize(ret$probs),datacell) %*% ret$P))
    #         ret$Chisq <- sum((freq-fit)^2/fit) + (N-sum(fit))
    #     }
    #     ret$predcell <- data.frame(datacell,observed=freq,expected=round(fit,3)) # Table that gives observed vs. predicted cell counts
    #     ret$Gsq <- 2 * sum(freq*log(freq/fit))  # Likelihood ratio/deviance statistic
    # }
    y[y==0] <- NA
    ret$y <- data.frame(y)             # outcome variables
    ret$x <- data.frame(x)             # covariates, if specified
    ret$w <- w
    for (j in 1:J) {
        rownames(ret$probs[[j]]) <- paste("class ",1:R,": ",sep="")
        if(calc.se) {
          rownames(ret$probs.se[[j]]) <- paste("class ",1:R,": ",sep="")
          rownames(ret$probs.jk[[j]]) <- paste("class ",1:R,": ",sep="")
          }
        if (is.factor(data[,match(colnames(y),colnames(data))[j]])) {
            lev <- levels(data[,match(colnames(y),colnames(data))[j]])
            colnames(ret$probs[[j]]) <- lev
            ret$y[,j] <- factor(ret$y[,j],labels=lev)
        } else {
            colnames(ret$probs[[j]]) <- paste("Pr(",1:ncol(ret$probs[[j]]),")",sep="")
            if(calc.se) {
              colnames(ret$probs.se[[j]]) <-   paste("Pr(",1:ncol(ret$probs.se[[j]]),")",sep="")
              colnames(ret$probs.jk[[j]]) <-   paste("Pr(",1:ncol(ret$probs.jk[[j]]),")",sep="")
              }
        }
    }
    ret$N <- N                         # number of observations
    ret$maxiter <- maxiter             # maximum number of iterations specified by user
    ret$resid.df <- min(ret$N,(prod(K.j)-1))-ret$npar # number of residual degrees of freedom
    class(ret) <- "poLCA"
    if (graphs) plot.poLCA(ret)
    if (verbose) print.poLCA(ret)
    ret$time <- Sys.time()-starttime   # how long it took to run the model
    }
    ret$call <- match.call()
    return(ret)
}
