\name{calcMaxRadius}
\alias{calcMaxRadius}
\alias{calcMaxRadius-methods}
\alias{calcMaxRadius,matrix-method}
\title{
Compute maximum radius.
}
\description{
Compute maximum distance of all pixels from the centroid of the image. 
Alternatively, it is the radius of the the smallest disk that covers the entire image with its origin at the centroid. 
The radius is determined by:

\deqn{ r = \frac{m_{00}}{2} \sqrt{\frac{N}{M}+ \frac{M}{N}}}

where

\deqn{ m_{pq} = \sum_{1}^{N} \sum_{1}^{M} (x^p)(y^q)f(x,y) }


}
\usage{
calcMaxRadius(I, center)
}
\arguments{
\item{I}{
A 2-dimensional image matrix. 
}
\item{center}{
The x and y coordinates of the image centroid, input as \code{c(x0, y0)}. 
If missing, the centroid will be calculated automatically.
}
}

\value{
\item{maxRadius}{The radius of the the smallest disk that covers the entire image with its origin at the centroid.}
}
\author{
Tan Dang
}
\note{
The value calculated by this function is used to obtain the polar coordinates of the pixels in the image.
The coordinates of the pixels are scaled with this value so that they will be between -1 and 1 before 
the polar coordinates are calculated. This is necessary because complex moments are only defined on the unit disk. 
}
\seealso{
\code{\linkS4class{CmplxIm}},\code{\link{calcCentroid}},\code{\link{polarXY}}
}
\examples{
data(mandril)
I=rowSums(img,dims=2);

maxRadius= calcMaxRadius(I);
maxRadius= calcMaxRadius(I, calcCentroid(I));
}

\keyword{utilities}
