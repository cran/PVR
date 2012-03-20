\name{PVRdecomp}
\alias{PVRdecomp}

\title{
Phylogenetic distances matrix (eigen)decomposition.
}
\description{
The PVRdecomp function decompose phylogenetic distance matrices (computed based on phylogenies) into a set of orthongonal eigenvectors
}
\usage{
PVRdecomp(phy, type = "newick", dist = NULL, scale = FALSE, ...)
}

\arguments{
  \item{phy}{
An object of class phylo that contains a ultrametric phylogeny with branch lengths or a character string that represents a phylogeny file (must have full path) of type given by type.
}
  \item{type}{
A character string with the phylogeny file format. It can be "newick" or "nexus". Used only if phy is a character string, ignored otherwise.
}
  \item{dist}{
A phylogenetic distance matrix. Not used.
}
  \item{scale}{
Logical. Should the phylogenetic distances be scaled into the range of 0 to 1. Default is FALSE.
}
  \item{\dots}{
Not used.
}
}

\value{
The PVRdecomp funtion returns a object of (S4) class PVR with an Eigen slot that contains eigenvalues and eigenvectors computed after the matrix decomposition.

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
 \code{\link{PVR}},  \code{\link{PSR}}
}
\examples{
library(ape)
tree <- rcoal(10)
#Decomposing phylogenetic distance matrix derived from tree into a set of orthogonal vectors
x <- PVRdecomp(tree, scale = TRUE)
str(x)
}

\keyword{ phylogenetic signal }
\keyword{ phylogeny }
