###### non-penalized DREGAR Algorithm #####
dregar6 <- function(data , da , ar , mselection=4 , type = 'alasso', 
                    normalize = FALSE , 
                    iteration = 15, intercept=FALSE){
  penalized = TRUE
  if(!penalized) {stop('This algorithm does not work when there is no penalization in the model! Use dregar2() function instead. ')}
  x = data[,-c(1,2),drop=FALSE];  y = data[,2]; n=length(y); BIC=10^63
  p=da;q=ar
  if(normalize){
    y = y-mean(y)
    x = scale(x)
  }
  final.phi=final.theta=final.beta=0
  #CP=1; AIC=2; GCV=3; BIC=4; type='enet', 'alasso'
  if(any(c(p,q)<1)) stop ('p and q must be greater than 1!')
  lambda = tau = gamma = 0
  counter = 2
  OPR = c(1:3,rep(4:6,iteration))
  msgps2 <- function(X,y,intercept=FALSE,penalty,criteria,penalized){
    y = as.vector(y)
    X = as.matrix(X)
    if(!penalized){
      if(intercept!=0){
        o= lm(y~X)
      }else{
        o= lm(y~X+0)  
      }
      co=coef(o)
      tun=NA
    }else{
      o = msgps::msgps(X = X,y = y,intercept = intercept,penalty = penalty)
      co=coef(o)[,criteria]
      tun=o[[criteria+3]]$tuning
    }
    return(list(model=o,coef=co,tun=tun))
  }
  cbind2 = function(a,b,add=TRUE){
    if(add){
      r=cbind(a,b)
    }else{
      r=b
    }
    return(r)
  }
  cut<-function(x,intercept){
    if(intercept){
      a=x[1]
      b=x[-1]
    }else{
      a=0
      b=x
    }
    return(list(x=b,inc=a))
  }
  for (operate in OPR){
    if (operate == 1){
      # STEP 1 #############################
      Hp  = makematrix(y = y,p = max(p),rname = 'phi')[,p,drop=FALSE]
      os1 = msgps2(X = Hp,y = y,intercept = intercept,penalty = type,criteria = mselection,penalized = penalized)
      gamma[counter] = os1$tun
      
    }else if(operate == 2){
      # STEP 2 #############################
      es2 = y - cbind2(1,Hp,intercept) %*% os1$coef
      os2 = msgps2(X = x,y = es2,intercept = intercept,penalty = type,criteria = mselection,penalized = penalized)
      lambda[counter] = os2$tun
      
    }else if(operate==3){
      # STEP 3 #############################
      es3 = es2 - cbind2(1,x,intercept) %*% os2$coef
      Hq  = makematrix(es3,p = max(q),rname = 'theta')[,q,drop=FALSE]
      os3 = msgps2(X = Hq,y = es3,intercept = intercept,penalty = type,criteria = mselection,penalized = penalized)
      tau[counter] = os3$tun
      
    }else if (operate==4){
      # STEP 4 #############################
      es4 = y - cbind2(1,Hq,intercept) %*% os3$coef
      os4 = msgps2(X = x,y = es4,intercept = intercept,penalty = type,criteria = mselection,penalized = penalized)
      lambda[counter] = os4$tun
      
    }else if (operate == 5){
      # STEP 5 #############################
      es5 = es4 - cbind2(1,x,intercept)  %*% os4$coef
      os5 = msgps2(X = Hp,y = es5,intercept = intercept,penalty = type,criteria = mselection,penalized = penalized)
      gamma[counter] = os5$tun
      
    }else if(operate==6){
      # STEP 6 #############################
      es6 = y-cbind2(1,Hp,intercept) %*% os5$coef - cbind2(1,x,intercept)  %*% os4$coef
      Hq  = makematrix(es6,p = max(q),rname = 'theta')[,q,drop=FALSE]
      os6 = msgps2(X = Hq,y = es6,intercept = intercept,penalty = type,criteria = mselection,penalized = penalized)
      tau[counter] = os6$tun
      
    }
    
    if(counter>6){
      est.y = cbind2(1,x,intercept) %*% os4$coef + cbind2(1,Hp,intercept) %*% os5$coef + cbind2(1,Hq,intercept) %*% os6$coef
      e = y-est.y
      cf= c(os4$coef,os5$coef,os6$coef)
      BIC.temp = -2 * sum(dnorm(e,mean(e),sd(e),log = TRUE)) +  log(n) * sum(cf[cf!=0])
      if(BIC.temp <= BIC){
        final.phi=os5$coef
        final.theta=os6$coef
        final.beta=os4$coef
        BIC = BIC.temp
      }
    }
    
    cat('\r Iteration ',round((counter-3)/3),' /',iteration)
    counter = counter + 1
  }
  
  lambda = na.omit(lambda)
  gamma  = na.omit(gamma)
  tau    = na.omit(tau)
  
  return(list(phi=cut(os5$coef,intercept)$x,
              theta=cut(os6$coef,intercept)$x,
              beta=cut(os4$coef,intercept)$x,
              bic.phi=cut(final.phi,intercept)$x,
              bic.theta=cut(final.theta,intercept)$x,
              bic.beta=cut(final.beta,intercept)$x,
              lambda=tail(lambda,1),gamma=tail(gamma,1),tau=tail(tau,1),
              trace.l=as.vector(lambda),trace.g=as.vector(gamma),trace.t=as.vector(tau),
              mod.beta=os4$model, mod.phi=os5$model, mod.theta=os6$model,
              x=x,y=y,ntheta=q,nphi=p,
              intercept= cut(os4$coef,intercept)$inc+cut(os5$coef,intercept)$inc+cut(os6$coef,intercept)$inc,
              #normalize = normalize,
              # pen = penalized,
              # rep = iteration,
              sel = c('CP', 'AICc', 'GCV', 'BIC')[mselection] ,
              est.y=est.y,
              res=e
  )
  )
}
