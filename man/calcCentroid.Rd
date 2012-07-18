\name{calcCentroid}
\alias{calcCentroid}
\alias{calcCentroid-methods}
\alias{calcCentroid,matrix-method}

\title{
Calculate the centroid of a grayscale image
}

\description{
Calculate the centroid of a grayscale image using geometric moments. The centroid is based on pixel intensity.
The centroid, located at pixel coordinates (x0,y0) is calculated using the following equation:

\deqn{ x0 = \frac{m_{10}}{m_{00}} }
\deqn{ y0 = \frac{m_{10}}{m_{00}} }

where

\deqn{ m_{pq} = \sum_{1}^{N} \sum_{1}^{M} (x^p)(y^q)f(x,y) }

N is the horizontal dimensionality of the image and M is the vertical dimensionality of the image.
\deqn{ f(x,y) } is the pixel intensity of the pixel at location (x,y) in the image.
}
\usage{
calcCentroid(I)
}
\arguments{
  \item{I}{An image represented as a 2-dimensional grayscale matrix.}
}

\value{
numeric: \code{(x0, y0)}
	\item{x0}{
		x-coordinate of centroid
	}
	\item{y0}{
		y-coordinate of centroid
	}
}
\author{
Allison Irvine, Tan Dang
}

\references{
Flusser, J., Suk, T., Zitova, B. Moments and Moment Invariants in Pattern Recognition. 2009. John Wiley & Sons Ltd. ISBN: 978-0-470-69987-4
}

\seealso{
\code{\linkS4class{Image}},\code{\linkS4class{CmplxIm}},\code{\linkS4class{OrthIm}},\code{\link{polarXY}}
}
\examples{
data(pirate)
I=rowSums(img,dims=2);
center= calcCentroid(I);
}

\keyword{multivariate}
\keyword{utilities}