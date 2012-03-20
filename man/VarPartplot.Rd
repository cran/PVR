\name{VarPartplot}
\alias{VarPartplot}

\title{
Plot function for objects of class PVR to show variation partition components.
}
\description{
Plot a, b, c and d components of traits variation partition between phylogeny and environment (as defined by Legendre and Legendre, 1998).
}
\usage{
VarPartplot(x, ...)
}

\arguments{
  \item{x}{
An object of class PVR with a non NULL slot VarPart.
}
  \item{\dots}{
Parameters passed to the plot function.
}
}

\references{
Diniz-Filho, J.A.F., Sant'Ana, C.E.R. and Bini, L.M. (1998). An eigenvector method for estimating phylogenetic inertia. Evolution 52:1247-1262.

Legendre, P. and Legendre, L. (1998). Numerical ecology, 2nd Englished. Elsevier.

Diniz-filho, J.A.F., Rangel, T.F., Santos, T. and Bini, L.M. (2012). Exploring patterns of interespecific variation in quantitative traits using sequential phylogenetic eigenvector regressions. Evolution, 66(4):1079-1090.

Diniz-filho, J.A.F., Bini, L.M., Rangel, T.F., Morales-Castilla, I., Olalla-Tarraga, M.A., Rodriguez, M.A. and Hawkins, B.A. (2012). On the selection of phylogenetic eigenvectors for ecological analyses. Ecography, 35:239-249.
}
\author{
Santos, T; Diniz-Filho, J.A.F.; Rangel, T.F.; Bini, L.M.
}

\seealso{
 \code{\link{PSR}},  \code{\link{PVRdecomp}},  \code{\link{PVR}},  \code{\link{PSRplot}}
}
\examples{
library(ape)
tree <- rcoal(10)
#Decomposing phylogenetic distance matrix derived from tree into a set of orthogonal vectors
x <- PVRdecomp(tree)
trait <- runif(10)
envvar <- runif(10)
y <- PVR(x, trait = trait, envVar = envvar, method = "moran")
VarPartplot(y)
}

\keyword{ Phylogenetic signal }
\keyword{ PSR }
