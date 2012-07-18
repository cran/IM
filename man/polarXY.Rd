\name{polarXY}
\alias{polarXY}
\alias{polarXY-methods}
\alias{polarXY,matrix,numeric-method}
\title{
Calculate polar coordinates of elements in a matrix
}
\description{
Computes the radius and and angle of each element's location in a matrix with respect to the centroid. 
}
\usage{
polarXY(I, centroid)
}

\arguments{
  \item{I}{
	\code{"Matrix"} representation of an image, or any matrix.
}
  \item{centroid}{
	two numbers, the x coordinate of the centroid and the y coordinate of the centroid.
	The function \code{\link{calcCentroid}} may be used to estimate the image centroid.
}
}
\value{
\code{"list"}
\item{radius}{The radius of the polar coordinates of each element in the input matrix. }
\item{theta}{The angle of the polar coordinates of each element in the input matrix. }
}
\author{
Allison Irvine, Tan Dang
}
\seealso{
\code{\link{polarTransform}},
\code{\link{calcCentroid}}
}
\examples{
data(lena);
I=rowSums(img,dims=2);
center= calcCentroid(I);
result= polarXY(I, center);
}
\keyword{utilities}
