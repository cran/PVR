\name{PVR}
\alias{PVR}

\title{
Phylogenetic eigenvectors regression.
}
\description{
The phylogenetic eigenvector regression (PVR) starts by performing an eigendecomposition of a pairwise double-centered phylogenetic distance matrix between species. The eigenvectors (representing the traits under analysis) estimated values express phylogenetic trends in data and residuals express independent evolution of each species.
}
\usage{
PVR(x, phy = NULL, trait = NULL, envVar = NULL, method = "moran", weights = NULL, scaled = FALSE, significance = TRUE, alternative = "two.sided", sig.treshold = 0.05, MI.treshold = 0.05, psr.treshold = 0.01, accvalue.treshold = 0.9, ...)
}

\arguments{
  \item{x}{
An object of class PVR (created by the PVRdecomp function) or class PSR (requiered by the "PSR" method).
}
  \item{phy}{
An object of class phylo that contains an ultrametric phylogeny.
}
  \item{trait}{
A vector, data frame or matrix that contains traits sets (for data frames and matrices, each column must represent a trait set).
}
  \item{envVar}{
A vector, data frame or matrix that contains environmental variables. Used to estimates the variation of a trait set that is explained by phylogeny and by environment.
}
  \item{method}{
Character string. A name for the eigenvectors selection method. It can be "moran", "stepwise", "psr" or "sequential".
}
  \item{weights}{
Weighting matrix based on phylogentic distances used in the "moran" method. If no weights matrix is provided, weights will be set to max(D) - Dij, where D is the phylogenetic distance matrix. 
}
  \item{scaled}{
Logical. Should the phylogenetic distances be scaled into the range of 0 to 1. Default is FALSE.
}
  \item{significance}{
Logical. Should the eigenvectors selected by the "moran" method be selected by the significance of residuals autocorrelation. If FALSE the eigenvectors will be selected by Moran's I values.
}
  \item{alternative}{
The alternative hypothesis used to compute the significance level when selecting eigenvectors by the "moran" method.
}
  \item{MI.treshold}{
Minimum residuals Moran's I value used to select eigenvectors when significance is FALSE.
}
  \item{sig.treshold}{
The significance treshold used to select eigenvectors by the "moran" method.
}
  \item{psr.treshold}{
The minimum acumulate R2 gain treshold used to select eigenvectors by the "PSR" method.
}
  \item{accvalue.treshold}{
Relative accumulated eigenvalue treshold use to select the eigenvectors by the "sequential" method.
}
  \item{\dots}{
Parameters passed to the stepwise regression used in the "AIC" method
}
}

\value{
A PVR class object.
}
\references{
Diniz-Filho, J.A.F., Sant'Ana, C.E.R. and Bini, L.M. (1998). An eigenvector method for estimating phylogenetic inertia. Evolution, 52:1247-1262.

Legendre, P. and Legendre, L. (1998). Numerical ecology, 2nd Englished. Elsevier.

Desdevises, Y., Legendre, P., Azouzi, L. and Morand, S. (2003). Quantifying phylogenetic structured environmental variation, Evolution, 57(11):2647-2652

Diniz-filho, J.A.F., Rangel, T.F., Santos, T. and Bini, L.M. (2012). Exploring patterns of interespecific variation in quantitative traits using sequential phylogenetic eigenvector regressions. Evolution, 66(4):1079-1090.

Diniz-filho, J.A.F., Bini, L.M., Rangel, T.F., Morales-Castilla, I., Olalla-Tarraga, M.A., Rodriguez, M.A. and Hawkins, B.A. (2012). On the selection of phylogenetic eigenvectors for ecological analyses. Ecography, 35:239-249.
}
\author{
Santos, T; Diniz-Filho, J.A.F.; Rangel, T.F.; Bini, L.M.
}

\seealso{
 \code{\link{PSR}},  \code{\link{PVRdecomp}},  \code{\link{PSRplot}},  \code{\link{VarPartplot}}
}
\examples{
library(ape)
tree <- rcoal(10)
#Decomposing phylogenetic distance matrix derived from tree into a set of orthogonal vectors
x <- PVRdecomp(tree)
trait <- runif(10)
y <- PVR(x, trait = trait, method = "moran")
str(y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Phylogenetic signal }
\keyword{ Phylogenetic models }% __ONLY ONE__ keyword per line
