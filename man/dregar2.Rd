\name{dregar6}
\alias{dregar6}
\title{
Estimating adaptive/non-adaptive DREGAR coefficients using an iterative 6-step procedure
}
\description{
 Estimating coefficients for penalized dynamic regression in the presence of autocorrelated residuals using an iterative 6-step procedure.
}
\usage{
dregar6(data , da, ar, mselection = 4, type = "alasso", 
      normalize = FALSE, iteration = 15,  intercept=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
  Data matrix of order (time, response, covariates)
}
  \item{da}{
      A vector of lags. Autoregressive orders for the response. For example 1:p for all lags from 1 to p
}
  \item{ar}{
  A vector of lags. Autoregressive orders for residuals. For example 1:q for all lags from 1 to q
}
  \item{mselection}{
  Model selection criteria. Choosing among 1 (CP), 2 (AIC), 3 (GCV) and 4 (BIC)
}
  \item{type}{
  Type of penalty. Choosing between 'enet' and 'alasso' for DREGAR and adaptive-DREGAR penalties.
}
  \item{normalize}{
  Logical flag. Setting to TRUE to normalise data prior to analysis
}
  \item{iteration}{
  The number of iterations
}
  \item{intercept}{
  Logical flag. Setting to TRUE to have intercept in the model.
}

}

\author{
  Hamed Haselimashhadi <hamedhaseli@gmail.com>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    %% ~~objects to See Also as \code{\link{help}}, ~~~
    \code{\link{dregar2}},
    \code{\link{generateAR}} ,
    \code{\link{sim.dregar}}
  }
  
\examples{
 par(mfrow=c(2,2))
  m=sim.dregar(n=500 ,  beta=1:4, phi=generateAR(2), theta=.1, 
               n.z.coeffs=3 , plot=TRUE) # generating data
  r=dregar6(data=m$rawdata, da = 1:3,
      ar = 1:2,mselection = 4,
      type='alasso')# estimating parameters using (non-apdative) DREGAR 
  round(cbind(
    true      = c(phi=c(m$phi,0),theta=c(m$theta,0),beta=m$beta),
    estimates = c(phi=r$phi,theta=r$theta,beta=r$beta)
  )
  ,3
  )
  plot(r$mod.phi,main='phi')
  plot(r$mod.theta,main='theta')
  plot(r$mod.beta,main='beta')
}
