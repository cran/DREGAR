findlastZero <- function  ( x ,force = FALSE) { 
  lx = length  ( x ) 
  out.p = lx
  for  ( i in 1:lx ) { 
    if  ( round( x [ lx - i + 1], 6) == 0 ) { 
      out.p = out.p - 1
    } else { 
      if(force == FALSE){
        break
      }
    } 
  } 
  return  ( out.p ) 
} 
fillWithZero <- function ( x, l ) { 
  nx = length ( x )
  if ( l > nx ) { 
    if ( nx != l ) { 
      otp = c ( x, rep ( 0, l-nx ) )
    }else{ 
      otp = x
    } 
  }else{ 
    otp = x    
  } 
  return ( otp )
} 
countTonumber<-function(x){
  if(x <= 0){
    return(0)
  }else{
    return(1:x)
  }
}

splitname<-function(btas,x,ar,ma){
  btas = as.matrix(btas,nrow=1)
  p=ar
  q=ma
  if(p>0 & q>0){  
    ar.v=btas[1:p]
    ma.v=btas[(p+1):(p+q)]
    if  ( length ( colnames ( x )  )  != 0 ) {
      rownames ( btas ) <-c ( paste ( 'Hp.', 1:ar, sep = '' ) , paste ( 'Hq.', 1:ma, sep = '' ) , colnames ( x )  ) 
    }else{
      rownames ( btas ) <-c ( paste ( 'Hp.', 1:ar, sep = '' ) , paste ( 'Hq.', 1:ma, sep = '' ) , paste ( 'Var.', 1:ncol ( x ) , sep = '' )  ) 
    }
    btas=btas[-c(1:(p+q)),]
  }else if(p>0 & q<1){
    ar.v=btas[1:p]
    ma.v=0
    if  ( length ( colnames ( x )  )  != 0 ) {
      rownames ( btas ) <-c ( paste ( 'Hp.', 1:ar, sep = '' ) , 'H.q' , colnames ( x )  ) 
    }else{
      rownames ( btas ) <-c ( paste ( 'Hp.', 1:ar, sep = '' ) , 'H.q' , paste ( 'Var.', 1:ncol ( x ) , sep = '' )  ) 
    }  
    btas=btas[-c(1:(p+1)),]
  }else if(p<1 & q>0){
    ar.v=0
    ma.v=btas[2:(q+1)]
    if  ( length ( colnames ( x )  )  != 0 ) {
      rownames ( btas ) <-c ( 'Hp' , paste ( 'Hq.', 1:ma, sep = '' ) , colnames ( x )  ) 
    }else{
      rownames ( btas ) <-c ( 'Hp' , paste ( 'Hq.', 1:ma, sep = '' ) , paste ( 'Var.', 1:ncol ( x ) , sep = '' )  ) 
    }
    btas=btas[-c(1:(q+1)),]
  }else{
    ar.v=0
    ma.v=0
    if  ( length ( colnames ( x )  )  != 0 ) {
      rownames ( btas ) <-c ( 'Hp','Hq', colnames ( x )  ) 
    }else{
      rownames ( btas ) <-c ( 'Hp','Hq', paste ( 'Var.', 1:ncol ( x ) , sep = '' )  ) 
    }
    btas=btas[-c(1:2),]
  }
  return(list(ar.v=ar.v,ma.v=ma.v,btas=btas))
}


vectorOutliers<-function(x,percent=FALSE){
  if(length(x)==0 || is.null(x)==TRUE){
    stop ('X must be a vector')
  }
  if(percent==TRUE){
    n=length(x)
  }
  o=c()
  qnt <- quantile(x, probs=c(.25, .75), na.rm = TRUE)  
  H <- 1.5 * IQR(x, na.rm = TRUE)
  otl=which(x < (qnt[1] - H)) 
  otg=which(x > (qnt[2] + H)) 
  o=c(otl,otg)
  if (percent==FALSE){n=1}
  lo=length(o)/n
  if(length(o)>=0){
    return(lo)
  }else{
    return(0)
  }
}



# ---- PRODUCE Hp and Hq ----#
`makematrix`<-function (y,p=3,beforestart.mean=NA,rname='H'){
  if (p>0 && is.null(p)==FALSE && length(p)>0){
    ny=length(y)
    H=matrix(ncol=p,nrow=ny)
    progress=1
    for (i in 1:ny){
      if(progress <= p){
        for (j in progress:1){
          H[i,progress-j]=y[j]  #+1
        }
      }else {
        col=1
        
        for (j in progress:(progress - p +1)){ 
          H[i,col]=y[j-1] #-1
          col=col+1
        }
      }
      progress=progress+1
    }
    if(is.na(beforestart.mean)){
      f=function(x){
        x<-as.numeric(as.character(x)) #first convert each column into numeric if it is from factor
        x[is.na(x)] = mean(x, na.rm=TRUE,trim = .1) #convert the item with NA to median value from the column
        x #display the column
      }
      H=apply(H,2,f)
    }else{
      H[is.na(H)]= beforestart.mean
    }
    colnames(H)=paste(rname,1:p,sep='.')
    return(H)
  }else{
    return (0)
  }
}

