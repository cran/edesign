\name{interchange}
\alias{interchange}
\title{ Interchange algorithm for maximum entropy sampling }
\usage{
interchange(A, nf, ns, S.start, etol=0, mattest = TRUE)
}
\arguments{
 \item{A}{ Spatial covariance matrix \eqn{A}. }
 \item{nf}{ Number of stations are forced into every feasible solution. }
 \item{ns}{ Number of stations have to be added to the existing network. }
 \item{S.start}{ Vector that gives the \eqn{ns} indices contained in
   the initial solution of the dimension \eqn{dim(A)[1]-nf} that should
   to be improved.  }
  \item{etol}{ Tolerance for checking positve definiteness (default 0)}
 \item{mattest}{ Toggles testing matrix \code{A} for symmetry and
   positive definiteness (default \code{T})}
}
\description{
  Starting point is a network \eqn{A[F]} with \eqn{nf} points.
  Now one has to select
\eqn{ns} points of a set of candidate sites to augment the existing
network. The aim of maximum entropy sampling is to select a feasible
D-optimal design that maximizes the logarithm of the determinant of all
principal submatrices of \eqn{A} arising by this expansion.

The interchange algorithm improves a feasible initial solution directly 
given or obtained by the greedy or dual greedy algorithm for maximum
entropy sampling.

It is also possible to improve the initial solution for the
construction of a completely new network, that means \eqn{nf=0}, but in
this case the interchange algorithm fails for \eqn{ns=1}.
}
  \details{
 \eqn{A[F]} denotes the principal submatrix of \eqn{A} having rows and columns
   indexed by \eqn{1..nf}. 
}
\value{
  A object of class \code{\link{monet}} containing the following
  elements:
  \item{S.start}{Vector containing the indices of the added sites in the
    initial solution or 0 for the other sites.}
  \item{S}{Vector containing the indices of the added sites in the solution
     or 0 for the other sites.}
  \item{det}{Determinant of the principal submatrix indexed by the
    initial solution.}
}
\references{ 
Ko, Lee, Queyranne, An exact algorithm for maximum entropy sampling,
Operations Research 43 (1995), 684-691.

Gebhardt, C.: Bayessche Methoden in der geostatistischen
Versuchsplanung. PhD Thesis, Univ. Klagenfurt, Austria, 2003

O.P. Baume, A. Gebhardt, C. Gebhardt, G.B.M. Heuvelink and J. Pilz:
Network optimization algorithms and scenarios in the context of
automatic mapping. Computers & Geosciences 37 (2011) 3, 289-294
}
\author{ C. Gebhardt }

\seealso{ \code{\link{greedy}}, \code{\link{dualgreedy}}}

\examples{
x <- c(0.97900601,0.82658702,0.53105628,0.91420190,0.35304969,
       0.14768239,0.58000004,0.60690101,0.36289026,0.82022147,
       0.95290664,0.07928365,0.04833764,0.55631735,0.06427738,
       0.31216689,0.43851418,0.34433556,0.77699357,0.84097327)
y <- c(0.36545512,0.72144122,0.95688671,0.25422154,0.48199229,
       0.43874199,0.90166634,0.60898628,0.82634713,0.29670695,
       0.86879093,0.45277452,0.09386800,0.04788365,0.20557817,
       0.61149264,0.94643855,0.78219937,0.53946353,0.70946842)
A <- outer(x, x, "-")^2 + outer(y, y, "-")^2
A <- (2 - A)/10
diag(A) <- 0
diag(A) <- 1/20 + apply(A, 2, sum)

S.c<-c(0,7,0,9,0,11,0,13,14,0,0,0,0,0,0)
interchange(A,5,5,S.c)
interchange(A,5,5,greedy(A,5,5)$S)
interchange(A,5,5,dualgreedy(A,5,5)$S)
}
\keyword{ spatial }%-- one or more ...
\keyword{ design }