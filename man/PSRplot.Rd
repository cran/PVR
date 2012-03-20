\name{PSRplot}
\alias{PSRplot}

\title{
Plot function for objects of class PSR 
}
\description{
Plot PSR curve along with null and neutral expectancies curves.
}
\usage{
PSRplot(x, info = c("area", "null", "neutral", "both"), ...)
}

\arguments{
  \item{x}{
An object of class PSR
}
  \item{info}{
Type of information to be ploted. It can be "area" (to plot only the PSR curve), "null" (to plot PSR curve and null PSR curve), "neutral" (to plot PSR curve and neutral PSR curve) and "both" (to plot PSR curve and both null and neutral PSR curve).
}
  \item{\dots}{
Parameters passed to the plot function.
}
}

\references{
Diniz-Filho, J.A.F., Sant'Ana, C.E.R. and Bini, L.M. (1998). An eigenvector method for estimating phylogenetic inertia. Evolution 52:1247-1262.

Diniz-filho, J.A.F., Rangel, T.F., Santos, T. and Bini, L.M. (2012). Exploring patterns of interespecific variation in quantitative traits using sequential phylogenetic eigenvector regressions. Evolution, 66(4):1079-1090.

Diniz-filho, J.A.F., Bini, L.M., Rangel, T.F., Morales-Castilla, I., Olalla-Tarraga, M.A., Rodriguez, M.A. and Hawkins, B.A. (2012). On the selection of phylogenetic eigenvectors for ecological analyses. Ecography, 35:239-249.
}
\author{
Santos, T; Diniz-Filho, J.A.F.; Rangel, T.F.; Bini, L.M.
}

\seealso{
 \code{\link{PSR}},  \code{\link{PVRdecomp}},  \code{\link{PVR}},  \code{\link{VarPartplot}}
}
\examples{
library(splancs)
#Creating a 10 tips ultrametric random phylogeny
library(ape)
tree <- rcoal(10)
#Decomposing phylogenetic distance matrix derived from tree into a set of orthogonal vectors
x <- PVRdecomp(tree)
trait <- runif(10)
res <- PSR(x, trait = trait, null.model = TRUE, neutral.model = TRUE, times = 10)
PSRplot(res, info = "both")
}

\keyword{ PSR area }
\keyword{ PVR }
