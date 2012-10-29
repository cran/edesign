\name{monet}
\title{Monet (Monitoring Network) object}
\author{C. Gebhardt}
\alias{monet}
\arguments{
  depending on the specific \code{method} more fields may exist, see
  apropriate help pages. All objects at least contain:
  \item{A}{covariance matrix of all stations}
  \item{na}{ Number of all stations. }
  \item{nf}{ Number of stations that are forced into every feasible solution. }
  \item{ns}{ Number of stations that have to be added to the existing
    network. }
  \item{ne}{ Number of eligible stations. }
  \item{method}{ One of \code{"greedy"}, \code{"dual greedy"},
    \code{"interchange"} (heuristics) or \code{"maxentropy"} (exact
    solution). 
  }
  \item{S}{selected subset of the eligible points, use
    \code{\link{monet.selected}} to convert into reusable index vector.}
  \item{det}{determinant of selected covariance matrix}
}
\description{
  An \code{monet} object is created with \code{\link{greedy}},
  \code{\link{dualgreedy}}, \code{\link{interchange}} or
  \code{\link{maxentropy}}.

  The object only contains indices of the rows/columns of the given
  covariance matrix \eqn{A}.

  Indices of fixed stations start with 1, indices of fixed and eligible
  stations are not mixed, i.e. the last \eqn{ne} stations are the
  eligible ones, the first \eqn{nf} are the fixed ones.

  The solution contains the indices counted absolute within all stations.
}
\seealso{\code{\link{print.monet}}
}
\keyword{spatial}
\keyword{ design }
