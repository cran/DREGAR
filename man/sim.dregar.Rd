\name{sim.dregar}
\alias{sim.dregar}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
  %%  ~~function to do ... ~~
   Simulating data from DREGAR model
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Simulating a mean zero Gaussian lagged response regression in the presence of autocorrelated residuals
}
\usage{
  sim.dregar(n = 500    , beta = 1, ind = FALSE  , 
             phi = .3, theta = .5, var = 1 , n.z.coeffs=0,
             shuffle = FALSE     , plot = FALSE  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    %%     ~~Describe \code{n} here~~
      The number of data points to be simulated
  }
  \item{beta}{
    %%     ~~Describe \code{beta} here~~
      Regression coefficients
  }
  \item{ind}{
    %%     ~~Describe \code{ind} here~~
   Logical flag. If TRUE then observations are assumed to be independent. Otherwise they are generated from random AR(1) processes. In both cases, variables are assumed to be mutually independent and follow Gaussian distribution.
  }
  \item{phi}{
    %%     ~~Describe \code{phi} here~~
      Dynamic coefficient(s)
  }
  \item{theta}{
    %%     ~~Describe \code{theta} here~~
      Residuals coefficient(s)
  }
  \item{var}{
    %%     ~~Describe \code{var} here~~
      Variance of the error term
  }
  \item{n.z.coeffs}{
    Number of zero coefficients if needed
  }
  \item{shuffle}{
    Logical flag. If TRUE shuffle coefficients. Otherwise data are grouped corresponded to non-zero and zero coefficients.
  }
  \item{plot}{
    %%     ~~Describe \code{plot} here~~
      Logical flag. Plot response
  }
}


\author{
  %%  ~~who you are~~
    Hamed Haselimashhadi <hamedhaseli@gmail.com>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~
  
  \seealso{
    %% ~~objects to See Also as \code{\link{help}}, ~~~
    %%\code{\link{regarma}}, \code{\link{iregarma}}
    \code{\link{dregar}} , 
    \code{\link{generateAR}}
  }
  
\examples{
  simdata=sim.dregar(n = 100 , beta = 1,
    ind = TRUE , phi = .40 , theta = -.25,
    var = 1 , plot = TRUE)
  str(simdata)
}

