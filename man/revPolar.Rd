\name{revPolar}
\alias{revPolar}
\alias{revPolar-methods}
\alias{revPolar,numeric,list-method}
\title{
	Inverse transform of the \code{polarTransform} method
}
\description{
Performs an inverse polar transform of an image by calculating the Euclidean coordinates of the rearranged pixels in an image
transformed by the \code{polarTransform} method. See the documentation for \link{polarTransform} for details on the 
calculation of the polar image coordinates.
}
\usage{
	revPolar(d,params)
}

\arguments{
	\item{d}{
		The dimensions of the original image in Euclidean space. 
	}
	\item{params}{
		A list containing the resulting transformed image and parameters from the \link{polarTransform} method. 
		This list contains the following items:
		PI: The transformed image returned by the function \code{polarTransform}.
		pAxis: The principal axis used to shift the polar coordinates of the pixels in the transformed image.
		resolution: An integer which determines the number of angle values between 0 and 2pi which was used to create the transform.
		center: The x and y coordinates of the centroid used in the transform.  
		scale: The scaling factor used to transform the image.
	}
}

\details{
This function is only for recovering the original image from an image transformed using the \link{polarTransform} method. 
The input arguments to \code{revPolar}, \code{PI,pAxis,resolution,center,scale} are returned by the \code{polarTransform}
method in a list. Note that the dimensions of the original image, the first argument \code{d}, must also be specified.
}

\value{
	\item{IR}{The inverse transformed image; an approximation of the original image from the polar transformed image.}
}

\author{
Allison Irvine, Tan Dang
}

\seealso{
	\code{\link{polarTransform}}, \code{\link{polarXY}}, \code{\link{calcCentroid}}
}
\examples{
#perform a polar transform on the original image
data(circles);
I=rowSums(img,dims=2)
R=polarTransform(I, 20, 100);
\dontrun{displayImg(R[[1]]);}

#now reverse the transform
IR = revPolar(dim(I),R);
\dontrun{displayImg(IR);}

}

\keyword{utilities}
