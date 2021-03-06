\name{dregar2}
\alias{dregar2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Estimating (just) adaptive-DREGAR coefficients using an iterative 2-step procedure
}
\description{
 Estimating coefficients for penalized/non-penalized dynamic regression in the presence of autocorrelated residuals using an iterative 2-step procedure.
}
\usage{
dregar2(data, da = 0, ar = 0, mselection = 4, 
       normalize = FALSE, penalized = TRUE, 
       iteration = 15)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
  Data matrix of order (time, response, covariates)
}
  \item{da}{
%%     ~~Describe \code{da} here~~
   A vector of lags. Autoregressive orders for response. For example 1:p for all lags from 1 to p
}
  \item{ar}{
%%     ~~Describe \code{ar} here~~
  A vector of lags. Autoregressive orders for residuals. For example 1:q for all lags from 1 to q
}
  \item{mselection}{
%%     ~~Describe \code{mselection} here~~
  Model selection criteria. Choosing among 1 (CP), 2 (AIC), 3 (GCV) and 4 (BIC)
}
  \item{normalize}{
%%     ~~Describe \code{normalize} here~~
  Logical flag. Setting to TRUE to normalize data prior to analysis
}
  \item{penalized}{
%%     ~~Describe \code{penalized} here~~
  Logical flag. Setting to TRUE to estimate coefficients through penalized likelihood. Otherwise the algorithm applies iterative OLS.
}
  \item{iteration}{
%%     ~~Describe \code{iteration} here~~
The number of iterations
}
}

\author{
  %%  ~~who you are~~
    Hamed Haselimashhadi <hamedhaseli@gmail.com>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

  \seealso{
    %% ~~objects to See Also as \code{\link{help}}, ~~~
    \code{\link{dregar6}} ,
    \code{\link{generateAR}} ,
    \code{\link{sim.dregar}}
  }
  
\examples{
par(mfrow=c(2,2))
  m=sim.dregar(n=500 ,  beta=1:4, phi=generateAR(2), theta=.3, 
               n.z.coeffs=3 , plot=TRUE) # generating data
  r=dregar2(data = m$rawdata,da = 1:3,ar = 1:2,mselection = 4,
            penalized = 1 )# estimating parameters using2-step adaptive-DREGAR
  round(cbind(
    true      = c(phi=c(m$phi,0),theta=c(m$theta,0),beta=m$beta),
    estimates = c(phi=r$phi,theta=r$theta,beta=r$beta)
  )
  ,3
  )
  plot(r$obj)
  acf(r$res, main='Residual ACF')
  pacf(r$res,main='Residual PACF')
}
