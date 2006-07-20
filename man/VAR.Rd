\name{VAR}

\alias{VAR}

\alias{print.varest}

\title{Estimation of a VAR(p)}

\encoding{latin1}

\description{
  Estimation of a VAR by utilising OLS per equation.
}

\usage{
VAR(y, p = 1, type = c("const", "trend", "both", "none"))
}

\arguments{
  \item{y}{Data item containing the endogenous variables}
  \item{p}{Integer for the lag order (default is p=1).}
  \item{type}{Type of deterministic regressors to include.}
}

\details{
  Estimates a VAR by OLS per equation. The model is of the following
  form:

  \deqn{
    \bold{y}_t = CD_t + A_1 \bold{y}_{t-1} + \ldots + A_p \bold{y}_{t-p} + \bold{u}_t
  }

  where \eqn{\bold{y}_t} is a \eqn{K \times 1} vector of endogenous
  variables and \eqn{u_t} assigns a spherical disturbance term of the
  same dimension. The coefficient matrices \eqn{A_1, \ldots, A_p} are of
  dimension \eqn{K \times K}. In addition, either a constant and/or a
  trend can be included as deterministic regressors (term \eqn{CD_T}, by
  setting the \code{type} argument to the corresponding value. The
  default is \code{const}. 
}

\value{
  A list with class attribute \sQuote{\code{varest}} holding the
  following elements:\cr

  \item{varresult}{list of \sQuote{\code{lm}} objects.}
  \item{resid}{A matrix object of the residuals.}
  \item{datamat}{The data matrix of the endogenous and explanatory variables.}
  \item{y}{The data matrix of the endogenous variables}
  \item{type}{A character, specifying the deterministic regressors.}
  \item{p}{An integer specifying the lag order.}
  \item{K}{An integer specifying the dimension of the VAR.}
  \item{obs}{An integer specifying the number of used observations.}
  \item{totobs}{An integer specifying the total number of observations.}
  \item{restrictions}{Either \code{NULL} or a matrix object containing
    the zero restrictions of the VAR(p).}
  \item{call}{The \code{call} to \command{VAR()}.}
}

\author{Bernhard Pfaff}

\references{
  Hamilton, J. (1994), \emph{Time Series Analysis}, Princeton
  University Press, Princeton.

  Lütkepohl, H. (2006), \emph{New Introduction to Multiple Time Series
    Analysis}, Springer, New York.
}

\seealso{
  \code{\link{summary.varest}}, \code{\link{plot.varest}},
  \code{\link{predict.varest}}  
}

\examples{
data(Canada)
VAR(Canada, p = 2, type = "none")
VAR(Canada, p = 2, type = "const")
VAR(Canada, p = 2, type = "trend")
VAR(Canada, p = 2, type = "both")
}

\keyword{regression}
\concept{VAR}
\concept{Vector autoregressive model}
\concept{VAR}

