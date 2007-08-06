\name{plot}

\alias{plot.varcheck}

\alias{plot.varest}

\alias{plot.varfevd}

\alias{plot.varirf}

\alias{plot.varprd}

\alias{plot.varstabil}

\alias{plot.vec2var}

\encoding{latin1}

\title{Plot methods for objects in vars}

\description{
  Plot method for objects with class attribute \code{varest},
  \code{vec2var}, \code{varcheck}, \code{varfevd}, \code{varirf},
  \code{varprd}, \code{varstabil}. 
}

\usage{
\method{plot}{varcheck}(x, ...)
\method{plot}{varest}(x, ...)
\method{plot}{vec2var}(x, ...)
\method{plot}{varfevd}(x, ...)
\method{plot}{varirf}(x, ...)
\method{plot}{varprd}(x, ...)
\method{plot}{varstabil}(x, ...)
}

\arguments{
  \item{x}{An object of one of the above classes.}
  \item{\dots}{Currently not used.}
}

\author{Bernhard Pfaff}

\references{
  Hamilton, J. (1994), \emph{Time Series Analysis}, Princeton
  University Press, Princeton.

  Lütkepohl, H. (2006), \emph{New Introduction to Multiple Time Series
    Analysis}, Springer, New York.

  Zeileis, A., F. Leisch, K. Hornik and C. Kleiber (2002), strucchange:
   An R Package for Testing for Structural Change in Linear Regression
   Models, \emph{Journal of Statistical Software}, \bold{7(2)}: 1-38,
   \url{http://www.jstatsoft.org/v07/i02/}  
}

\seealso{
  \code{\link{VAR}}, \code{\link{vec2var}}, \code{\link{fevd}},
  \code{\link{irf}}, \code{\link{predict}}, \code{\link{fanchart}},
  \code{\link{stability}}, \code{\link{arch}}, \code{\link{normality}},
  \code{\link{serial}}  
}

\examples{
\dontrun{
data(Canada)
var.2c <- VAR(Canada, p = 2, type = "const")
plot(var.2c)
## Diagnostic Testing
## ARCH test
archtest <- arch(var.2c)
plot(archtest)
## Normality test
normalitytest <- normality(var.2c)
plot(normalitytest)
## serial correlation test
serialtest <- serial(var.2c)
plot(serialtest)
## FEVD
var.2c.fevd <- fevd(var.2c, n.ahead = 5)
plot(var.2c.fevd)
## IRF
var.2c.irf <- irf(var.2c, impulse = "e",
response = c("prod", "rw", "U"), boot = FALSE)
plot(var.2c.irf)
## Prediction
var.2c.prd <- predict(var.2c, n.ahead = 8, ci = 0.95)
plot(var.2c.prd)
## Stability
var.2c.stabil <- stability(var.2c, type = "Rec-CUSUM")
plot(var.2c.stabil)
}
}

\keyword{regression}
