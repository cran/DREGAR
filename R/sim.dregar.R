generateAR <- function ( n = 1 , l = -1 , u = 1, 
                         min.distance = .Machine$double.eps, 
                         sort.coeff = FALSE ) { 
  if (n > 0) {
    cat('\r please wait ... ')
    if (((u - l)/n) <= min.distance) {
      stop("Data are not possible to be generated. please decrease min.distance!")
    }
    l=rep(l,n)
    u=rep(u,n)
    k=(1:n)^(sqrt(n)/2)
    repeat {
      minroots = 0
      if (sort.coeff == FALSE) {
        ar = runif(n, l/k, u/k)
      }
      else {
        ar = sort(runif(n, l/k, u/k), decreasing = TRUE)
      }
      minroots <- min(abs(polyroot(z = c(1, -ar))))
      if (n > 1) {
        if (abs(minroots) > 1 & min(diff(abs(sort(ar)))) > 
            min.distance) {
          cat('\r                    \n')
          break
        }
      }
      else {
        if (abs(minroots) > 1) {
          break
        }
      }
    }
    return(ar)
  }
  else {
    return(0)
  }
} 


sim.dregar <- function ( n = 500, beta = 1 , ind  = FALSE, 
                         phi = .3 , theta = .5 , var = 1, 
                         n.z.coeffs = 0, shuffle = FALSE, 
                         plot = FALSE ) {    
  if ( n.z.coeffs > 0 ) { 
    beta = c ( beta, rep ( 0, n.z.coeffs ) ) 
    if ( shuffle == TRUE ) { 
      beta = sample ( beta, length ( beta ), replace  =  FALSE ) 
    } 
  } 
  p = length ( phi ) 
  q = length ( theta ) 
  lbeta = length ( beta ) 
  
  if ( ind  == TRUE ) { 
    x = rnorm ( n*length ( beta ),0, sd = 1 ) #sqrt( abs ( beta ) ) ) 
    x = matrix ( x, ncol = lbeta, nrow = n ) 
    x = rbind ( matrix ( 0, ncol = lbeta, nrow = p + q ), x ) 
    colnames ( x ) = paste ( 'X', 1:length ( beta ), sep = '.' ) 
  } else { 
    x = matrix ( 0, ncol = lbeta, nrow = n ) 
    for ( i in 1:lbeta ) { 
      x [, i ] = arima.sim ( n = n, list ( ar = sign(runif(1,-1,1))*runif(0,1)) , sd = 1 ) #sqrt( abs ( beta ) ) ) 
    } 
    x = rbind ( matrix ( 0, ncol = lbeta, nrow = p + q ), x ) 
    colnames ( x ) = paste ( 'X', 1:length ( beta ), sep = '.' ) 
  } 
  y = rep ( 0, p + q ) 
  t0 = p + q
  distrurbance = c () 
  for ( t in 1:n ) { 
    rs = 0
    rs = x [ t0 + t, ]  %*% as.matrix ( beta ) 
    
    js = 0
    js = phi %*% y [ ( t + t0 - 1 ) : ( t + t0 - p ) ] 
    
    
    ks = 0
    ks = theta  %*% ( as.matrix ( y [ ( t + t0 - 1 ) : ( t + t0 - q ) ] )  - x [ ( t + t0 - 1 ) : ( t + t0 - q ), ]  %*% as.matrix ( beta ) ) 
    
    s2k = 0
    for ( l in 1:q ) { 
      for ( i in 1:p ) { 
        s2k = s2k + theta [ l ] *phi [ i ] *y [ t + t0 - l - i ] 
      } 
    } 
    ks = ks - s2k
    distrurbance [ t ] = rnorm ( 1, 0, sqrt ( var ) ) 
    y [ t + t0 ] = rs + js + ks + distrurbance [ t ] 
  }  
  
  y = y [  -  ( 1:t0 ) ] 
  x = x [  -  ( 1:t0 ), ] 
  if ( plot == TRUE ) { 
    plot ( spline ( y ), type = 'l', ylab = 'y', xlab = 'Time', main = paste ( 'DA ( ', toString ( round ( phi, 2 ) ), ' ) \n AR ( ', toString ( round ( theta, 2 ) ), sep = ' ', ' ) ' ) ) 
  } 
  otpMatrix = cbind ( 1:length ( y ), y, x ) 
  colnames ( otpMatrix ) = c ( 'T', 'Y', paste ( 'X.', 1:lbeta, sep = '' ) ) 
  
  return ( list ( rawdata = otpMatrix, y = y, x = x, beta = beta, phi = phi, theta = theta, error = distrurbance ) ) 
} 
