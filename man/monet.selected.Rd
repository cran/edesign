\name{monet.selected}
\alias{monet.selected}
\alias{monet.solution}
\title{
  extract the selected station indices from a \code{\link{monet}} object
}
\description{
  This function can be used to get an index vector containing the
  additionally selected (\code{monet.selected}) rows/columns of a covariance matrix used in
  \code{"greedy"}, \code{"dual greedy"}, \code{"interchange"} or
  \code{"maxentropy"} (or the complete network, \code{monet.solution}).
}
\usage{
monet.selected(x)
monet.solution(x)
}
\arguments{
  \item{x}{
    object of class \code{\link{monet}}
  }
}
\value{
  integer vector
}
\author{
C. Gebhardt
}
\seealso{
  \code{\link{monet}}
}
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

MN <- greedy(A,5,5)
# internal representaion in monet object:
MN$S
# reusable index, only added stations:
monet.selected(MN)
# ... complete network:
monet.solution(MN)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }
