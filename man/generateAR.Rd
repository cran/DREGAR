\name{generateAR}
\alias{generateAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
    Generating stationary autoregressive coefficients
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Generating stationary autoregressive coefficients 
}
\usage{
  generateAR(n = 1, l = -1, u = 1    , 
             min.distance = .Machine$double.eps , 
             sort.coeff = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    %%     ~~Describe \code{n} here~~
      The number of coefficients
  }
  \item{l}{
    %%     ~~Describe \code{l} here~~
      Lower bound for coefficients
  }
  \item{u}{
    %%     ~~Describe \code{u} here~~
      Upper bound for coefficients
  }
  \item{min.distance}{
    %%     ~~Describe \code{min.distance} here~~
      Minimum distance amongst all pair of coefficients
  }
  \item{sort.coeff}{
      Logical flag. If TRUE, then resulting coefficients are sorted descending
  }
}

\author{
  %%  ~~who you are~~
    Hamed Haselimashhadi <hamedhaseli@gmail.com>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
  
  \seealso{
    %% ~~objects to See Also as \code{\link{help}}, ~~~
    \code{\link{dregar}} , 
    \code{\link{sim.dregar}}
  }
\examples{
  generateAR(3 , -1 , 1 , .01 )
}
