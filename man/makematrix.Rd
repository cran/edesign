\name{makematrix}
\alias{makematrix}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Covariance matrix}
\description{
  This function generates a symmetic positive definite (\eqn{n,n})
  covariance matrix using the algorithm in Ko et al. 
}
\usage{
makematrix(n)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{n}{number of rows and columns of the matrix}
}

\value{
  (\eqn{n,n}) covariance matrix.
}
\references{ 
Ko, Lee, Queyranne, An exact algorithm for maximum entropy sampling, Operations Research 43 (1995), 684-691. 
Gebhardt, C.: Bayessche Methoden in der geostatistischen Versuchsplanung. PhD Thesis, Univ. Klagenfurt, Austria, 2003
O.P. Baume, A. Gebhardt, C. Gebhardt, G.B.M. Heuvelink and J. Pilz: Network optimization algorithms and scenarios in the context of automatic mapping. Computers & Geosciences 37 (2011) 3, 289-294 
}
\author{C. Gebhardt}

\examples{
x.cov<-makematrix(5)
}

%\seealso{ \code{\link[kdesign]{kdesign.cov}} }

\keyword{ spatial }%-- one or more ...
\keyword{ design }

