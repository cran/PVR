\name{PSR}
\alias{PSR}

\title{
Creates a phylogenetic signal-representation (PSR) curve for a given phylogeny and trait set and computes its area.
}
\description{
Phylogenetic signal-representation (PSR) curve, built upon phylogenetic eigenvector regression (PVR). Sequential PVR models are fitted after successively increasing the number of eigenvectors and plotting their R2 against the accumulated eigenvalues. The PSR area, expressing deviations from Brownian motion, is strongly correlated with Blomberg's K-statistics, so nonlinear PSR curves reveal if traits are evolving at a slower or higher rate than expected. The PSR area is also correlated with phylogenetic half-life under an OU process, so both methods describe the relationship between interspecific variation and time since divergence among species shape.
}
\usage{
PSR(x, trait = NULL, null.model = FALSE, Brownian.model = FALSE, times = 1000)
}

\arguments{
  \item{x}{
A object of class PVR (created using PVRdecomp function).
}
  \item{trait}{
A vector, data frame or matrix that contains contiuous or binary traits sets (for data frames and matrices, each column must represent a trait set). For now, only the first colunm will be used.
}

  \item{null.model}{
Logical. Should the function estimates the null expectancy for the PSR area. Default is FALSE.
}
  \item{Brownian.model}{
Logical. Should the function estimates the Brownian (neutral) expectancy for the PSR area. Default is FALSE.
}
  \item{times}{
Number of iterations used in null.model and Brownian.model.
}
}

\value{
PSR function returns a object of (S4) class PSR that extends the (S4) class PVR, by adding: a PSR slot that contains the area (and p if null.model = TRUE) for each trait set, the cumulative eigenvalues and R2 that determine the curve; a slot with the null and neutral expectancies (if any).

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
 \code{\link{PVR}},  \code{\link{PVRdecomp}},  \code{\link{PSRplot}},  \code{\link{VarPartplot}}
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
plot(y)
}

\keyword{phylogenetic signal}
\keyword{phylogeny}
