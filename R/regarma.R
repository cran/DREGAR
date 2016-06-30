`regarma`<-function ( x, y, ar = 0, ma = 0, 
                      normalize = FALSE   ,
                      debug = FALSE, rep = 50 , 
                      pen=0, criteria=3) {
  if (pen) requireNamespace('msgps')
  BIC=c() ; pen = pen; criteria = criteria ; rep=max(rep,1)#CP AICc GCV BIC
  if ( debug == TRUE ) {Rprof ( tf <- "nplregarma.log",  memory.profiling = TRUE ) } 
  #------ Initializing ------  #  
  x = as.matrix ( x ) ; temX = x ; temY = y ; r = ncol ( x )  ; ly = length ( y ) ; 
  ar.value = ma.value = c () ;  intcpt = 0 ;  lar=length(ar) ; lma=length(ma) ; mar=max(ar) ; mma=max(ma)
  if(lar==1 && ar==0) {lar=0} ; if(lma==1 && ma==0){lma=0}
  coef.matrix = matrix(NA, nrow = (max(1,lar)+max(1,lma)+max(r,1)+2) ,ncol=rep)
  if ( normalize == TRUE ) {x = scale ( x ) ; y =  y - mean(y) }
  if ( max (mar , mma) >= ly ) {stop ( 'Error : p + q > T!' ,call. = 0) }
  #--------- REGARMA STEP 1 ----------
  if ( mar > 0 ) {
    Hp = makematrix ( y, mar, NA , 'Hp' ) [ , ar]
    if(pen){
      s1.regar = msgps::msgps(X = cbind ( Hp , x  ) , y =  as.vector ( y  ), 
                              penalty = 'alasso' , intercept = 1 )
      s1.coeff = coef(s1.regar)[,criteria]
    }else{
      s1.regar = lm ( as.vector ( y  ) ~ cbind ( Hp  , x  )   ) 
      s1.coeff = coef(s1.regar)
    }
    s1.haty  = cbind ( 1, Hp  , x  )  %*% s1.coeff #
    s1.error = y  - s1.haty  # residuals
  }else{ 
    if(pen){
      s1.regar = msgps::msgps(X = x , y =  as.vector ( y ), 
                              penalty = 'alasso' , intercept = 1)
      s1.coeff = coef(s1.regar)[,criteria]
    }else{
      s1.regar = lm ( as.vector ( y ) ~ cbind ( x )   ) 
      s1.coeff = coef(s1.regar)
    }
    s1.haty  = cbind ( 1, x )  %*% s1.coeff  
    s1.error = y - s1.haty # residuals
  }
  #--------- REGARMA STEP 2 -----#
  for ( j in 1:rep ) {    
    if ( mma>0 ) {
      if ( mar>0 ) {
        Hq = makematrix ( s1.error, mma, 0 , 'Hq' ) [ , ma]
        if(pen){
          s2.regarma = msgps::msgps(X = cbind ( Hp , Hq , x  ) ,y =  as.vector ( y ), 
                                    penalty = 'alasso' , intercept = 1 , STEP.max  =  10^5 )
          s2.regarma.coeffs = coef(s2.regarma)[,criteria]
        }else{
          s2.regarma = lm ( as.vector ( y ) ~ cbind ( Hp  , Hq , x )   ) 
          s2.regarma.coeffs = coef(s2.regarma)
        }
        s2.y.hat = cbind ( 1 , Hp , Hq , x )  %*% s2.regarma.coeffs # + s2.regarma$intercept
        s2.error = y - s2.y.hat
        intcpt   = as.vector ( s2.regarma.coeffs )  [ 1           ] 
        ar.value = as.vector ( s2.regarma.coeffs )  [ 2 : (lar+1) ] 
        ma.value = as.vector ( s2.regarma.coeffs )  [  (  lar  + 2 ) : (  lar   + lma +1)  ] 
        btas =     as.matrix ( s2.regarma.coeffs ) 
        btas = btas [ - ( 1: ( lar + lma + 1)  ) ,  ] 
      }else{
        Hq = makematrix ( s1.error, mma, 0 , 'Hq' ) [ , ma]
        if(pen){
          s2.regarma = msgps::msgps(X = cbind ( Hq , x  ), y = as.vector ( y ), 
                                    penalty = 'alasso' , intercept = 1 )
          s2.regarma.coeffs = coef(s2.regarma)[,criteria]
        }else{
          s2.regarma = lm ( as.vector ( y  ) ~ cbind ( Hq , x   )   )
          s2.regarma.coeffs = coef(s2.regarma)
        }
        s2.y.hat = cbind ( 1, Hq , x  )  %*% s2.regarma.coeffs 
        s2.error = y  - s2.y.hat      
        ar.value = 0
        intcpt   = as.vector ( s2.regarma.coeffs )[1]
        ma.value = as.vector ( s2.regarma.coeffs )  [ 2:(lma+1) ]  
        btas = as.matrix ( s2.regarma.coeffs  ) 
        btas = btas [ - ( 1: (lma+1)  ) ,  ]  #REMOVE AR/MA Coefficients
      }
    }else{
      if  ( mar>0 ) {  
        s2.regarma = s1.regar
        s2.regarma.coeffs = s1.coeff
        s2.y.hat = s1.haty 
        s2.error = s1.error
        intcpt   = as.vector ( s2.regarma.coeffs )[1]
        ar.value = as.vector ( s2.regarma.coeffs )  [ 2:(lar+1) ] 
        ma.value = 0
        btas = as.matrix ( s2.regarma.coeffs   ) 
        btas = btas [ - ( 1: (lar+1)  ) ,  ]  #REMOVE AR/MA Coefficients
      }else{
        s2.regarma = s1.regar
        s2.regarma.coeffs = s1.coeff
        s2.y.hat = s1.haty
        s2.error = s1.error
        ar.value = 0; ma.value = 0
        intcpt =as.matrix ( s2.regarma.coeffs  ) [1]
        btas = as.matrix ( s2.regarma.coeffs  ) [-1]
      }
    }
    new.coefs=c( 0, j,  fillWithZero ( ar.value, lar ) , fillWithZero ( ma.value, lma ), btas )
    coef.matrix[,j] =  new.coefs 
    L=sum(dnorm(s2.error, mean(s2.error,trim = .05), sd(s2.error), log=T))
    k=sum(s2.regarma.coeffs!=0)
    BIC[j]           = -2*L + k *  log(ly)
    if ( mar>0 && rep>1 && mma>0) {
      s1.haty  = cbind ( 1, Hp , x  )  %*% c ( intcpt , ar.value, btas ) 
      s1.error = y  - s1.haty 
      stop.criteria.value = sum(abs(  new.coefs[-c(1,2)] - coef.matrix[ -c(1,2), max(j - 1,1) ] ) )
      cat('\r ',j,'/' , rep, ' - ERR : ',stop.criteria.value,'')
    }else{
      break
    }
  } 
  b = which(BIC==min(BIC))[1];#b=rep
  #   if(length(BIC)>1){
  #     plot(BIC,type='b',xlim=c(0,rep+1),xlab='Iteration'); abline(v=b,col=2);text(x = b,y=max(BIC),labels=paste('--',(b),'--'),bg='white',font = 2,col=2)
  #   }
  btas = as.vector(coef.matrix[-c(1:2),b])
  vars = splitname(btas = btas , x , lar , lma)
  btas = vars$btas ;   ar.value = vars$ar.v ;   ma.value = vars$ma.v
  if ( debug == TRUE ) {Rprof ( NULL ) ; cat ( '\n ------- DEBUG DATA ------- \n' ) ; print ( summaryRprof ( tf )  ) }
  return ( list
           ( 
           obj = s2.regarma,
           y = temY,  x = temX, 
           nphi = ar,  ntheta = ma, 
           phi   = as.matrix (ar.value  ), 
           intercept= intcpt             ,
           theta = as.matrix ( ma.value ), 
           beta  = as.matrix( btas      ),
           sel = c('CP', 'AICc', 'GCV', 'BIC')[criteria] ,
           res   = as.matrix (s2.error  ), 
           est.y = as.matrix( s2.y.hat) 
           ) 
  )
}
