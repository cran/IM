\name{polarTransform}
\alias{polarTransform}
\alias{polarTransform-methods}
\alias{polarTransform,matrix,numeric,numeric,numeric-method}
\alias{polarTransform,matrix,numeric-method}
\title{
	Represent an image by plotting radius against angle
}
\description{
Perform a polar transform of an image by calculating the polar coordinates of the pixels in an image and rearranging the pixels
into a plot of angle (theta) against radius. The polar coordinates are calculated with respect to the centroid, determined by 
\code{\link{calcCentroid}} or the center of the image space. All transformed images will be shifted to principal axis determined by:

\deqn{\theta=\frac{1}{2}arctan(\frac{2\mu_{00}}{\mu_{10}-\mu_{01}})}

where
\deqn{\mu_{pq}= \sum_{1}^{N} \sum_{1}^{M}(x-x0)^p(y-x0)^q f(x,y)}
with (x0, y0) are centroid computed by \code{\link{calcCentroid}}

}
\usage{
polarTransform(I, resolution, scale, center)
}

\arguments{
  \item{I}{
	Image represented as a matrix or numeric array
	}
  \item{resolution}{
	An integer which determines the number of angle values between 0 and 2pi to be represented in the transform.
	}
  \item{scale}{
	Integer. Default is 0, no scaling in the transformation. If set to a positive integer, the maximum radius 
	used in the transformation will be set to \code{scale}, and the transformed image will have dimensions 
	\code{scale, (scale*resolution)}
	}
  \item{center}{
	The x and y coordinates of the centroid to be used in the transform. If this argument is not provided,
	the center of the image space will be the centroid. 
	}

}
\details{
This transformation is performed by finding the Euclidean coordinates of points which are contained within the 
smallest circle with origin at the centroid or center of image which can be contained within the image boundaries. 
The pixels whose coordinates most closely match these Euclidean coordinates are plotted at the corresponding polar coordinates.
At each radius, the perimeter is divided into \code{radius * resolution} equally spaced point. Then, the Euclidean coordinates 
of the points will be used in the transformation. This results in a upper-triangular matrix.
This transform was created to approximate rotational invariance for orthogonal moments.
\code{revPolar} is the inverse transform of \code{polarTransform}.
}
\value{
	\item{PI}{The transformed image.}
	\item{pAxis}{The principal axis of the image.}
	\item{resolution}{The resolution used for the transformation.}
	\item{scale}{The scale used for the transformation.}
	\item{center}{The centroid used for the transformation.}
}

\author{
Allison Irvine, Tan Dang
}

\seealso{
	\code{\link{revPolar}}, \code{\link{polarXY}}, \code{\link{calcCentroid}}
}
\examples{
data(circles);
I=rowSums(img,dims=2)
R=polarTransform(I, 6, 100);
displayImg(list(I,R[[1]]));
}

\keyword{utilities}
