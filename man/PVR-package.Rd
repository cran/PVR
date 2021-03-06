\name{PVR-package}
\alias{PVR-package}
\docType{package}
\title{
Computes \bold{P}hylogenetic eigen\bold{V}ectors \bold{R}egression and \bold{P}hylogentic \bold{S}ignal-\bold{R}epresentation curve (with null and neutral expectations).
}
\description{
Computes PVR and PSR curve along with some plot utilitties.
}
\details{
\tabular{ll}{
Package: \tab PVR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2018-05-29\cr
License: \tab GPL >= 2\cr
}

}
\author{
Santos, T; Diniz-Filho, J.A.F.; Rangel, T.F.; Bini, L.M.

Maintainer: Thiago Santos <thiago.santos@ufvjm.edu.br>
}
\references{
Diniz-Filho, J.A.F., Sant'Ana, C.E.R. and Bini, L.M. (1998). An eigenvector method for estimating phylogenetic inertia. Evolution 52:1247-1262.

Diniz-filho, J.A.F., Rangel, T.F., Santos, T. and Bini, L.M. (2012). Exploring patterns of interespecific variation in quantitative traits using sequential phylogenetic eigenvector regressions. Evolution, 66(4):1079-1090.

Diniz-filho, J.A.F., Bini, L.M., Rangel, T.F., Morales-Castilla, I., Olalla-Tarraga, M.A., Rodriguez, M.A. and Hawkins, B.A. (2012). On the selection of phylogenetic eigenvectors for ecological analyses. Ecography, 35:239-249.
}

\keyword{ package }
\seealso{

\code{\link{PVR-class}}, \code{\link{PSR-class}}, \code{\link{PVR}}, \code{\link{PSR}},  \code{\link{PVRdecomp}},  \code{\link{PSRplot}},  \code{\link{VarPartplot}}
}
\examples{
library(splancs)
#Creating a 10 tips ultrametric random phylogeny
library(ape)
tree <- rcoal(10)
#Decomposing phylogenetic distance matrix derived from tree into a set of orthogonal vectors
x <- PVRdecomp(tree)
trait <- runif(10)
y <- PSR(x, trait)
y
}
