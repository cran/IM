\name{OrthIm-class}
\Rdversion{1.1}
\docType{class}
\alias{OrthIm-class}

\alias{Moments<-,OrthIm-method}
\alias{Invariant<-,OrthIm-method}
\alias{Reconstruct<-,OrthIm-method}
\alias{initialize,OrthIm-method}
\alias{momentType<-,OrthIm-method}
\alias{plotPoly,OrthIm,numeric-method}
\alias{setOrder<-,OrthIm-method}
\alias{setParams<-,OrthIm-method}
\alias{transform<-,OrthIm-method}

\title{Class for handling continuous or discrete orthogonal moments analysis of an image}
\description{
This class contains an image and all information about the image required for computing a moments analysis
of the image. The object can be passed to functions to compute moments and reconstruct the image from the moments.
The types of moments that can be calculated using this object are: chebyshev, krawtchouk, hahn, continuous chebyshev,
legendre, and gegenbauer.
}
\section{Creating Objects}{
Objects can be created by calls of the form \code{new("OrthIm", img, filename)}.
First the constructor checks for an input image (img). This can be a multidimensional numeric array or a matrix. 
The image is converted to grayscale by summing over all colors and stored in slot "I". 
If an input image is not provided, a filename (filename) can be passed to the constructor.
The constructor will then read in the image file and convert it to grayscale. 
If a filename is given, it must have the extension ".jpg", ".png", or ".bmp".
}
\section{Slots}{
  \describe{
    \item{\code{polynomials}:}{Object of class \code{"list"} Contains the polynomials used to calculate moments. First list element contains the polynomials over the x coordinates, the second contains polynomials over the y coordinates. The number of rows is the pixel coordinates, the number of columns is the order of the polynomials. }
    \item{\code{params}:}{Object of class \code{"list"} Contains the parameters used to calculate the moments. First list element is the parameters for the moments over x coordinates, second is for the y coordinates. }
    \item{\code{order}:}{Object of class \code{"numeric"} Contains the order used to calculate moments. Two values may be given to calculate different orders in the x and y directions.  }
    \item{\code{I}:}{Object of class \code{"matrix"} Contains the original image, in grayscale. }
    \item{\code{dimensions}:}{Object of class \code{"numeric"} Contains the dimensions of the image. }
    \item{\code{centroid}:}{Object of class \code{"numeric"} Contains the x and y coordinates of the centroid (center point) of the image, used for calculating continuous orthogonal moments. }
    \item{\code{filename}:}{Object of class \code{"character"} Contains the filename of the image. }
    \item{\code{imType}:}{Object of class \code{"character"} Contains the original image file type. }
    \item{\code{momentType}:}{Object of class \code{"character"} Contains the types of moments being calculated. Up to two different types may be selected, the first being the type of moments calculated in the x direction, the second being for the y direction. If two different types are chosen, they must both be either discrete or continuous. Choices of discrete moments are: "cheby","krawt","hahn". Choices of continuous moments are: "chebycont","gegen","legendre". }
    \item{\code{moments}:}{Object of class \code{"matrix"} Contains the moments calculated from the image. Dimensions are order in y direction by order in x direction. }
    \item{\code{reconstruction}:}{Object of class \code{"matrix"} Contains the image reconstructed from the moments. }
    \item{\code{error}:}{Object of class \code{"character"} Contains the last error message produced from calling a class method. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{Image}"}, directly.
}
\section{Methods}{
  \describe{
    \item{Moments<-}{\code{signature(obj = "OrthIm")}: Calculate the moments of the image. Before calling this method the momentType, params, and order must be set. Usage: Moments(obj) <- 0 }
    \item{Invariant<-}{\code{signature(obj = "OrthIm")}: Calculate the invariants of the image by normalization before computing moments. \code{\link{polarTransform}} will be called with \code{resolution=value[1],scale=value[2]}. Before calling this method the momentType, params, and order must be set. Usage: Invariant(obj) <- c(6,100) }
    \item{Reconstruct<-}{\code{signature(obj = "OrthIm")}: Reconstruct the image from the moments, using moments up to a specified order. If no order value is given, the image will be reconstructed from all available moments. Useage: Reconstruct(obj) <- order }
    \item{initialize}{\code{signature(.Object = "OrthIm")}: Constructor. A call to new("OrthIm",img,filename) calls this function and creates an OrthIm object. The object slots I, centroid and dimensions are set. (And filename and imType if a filename was provided) }
    \item{momentType<-}{\code{signature(obj = "OrthIm")}: Set the moment type for the moments to be calculated. Takes one character string or an array of 2 character strings, for example c("type1","type2") }
    \item{plotPoly}{\code{signature(obj = "OrthIm", order = "numeric")}: Plot the polynomials used to calculate the moments in this object. The order, can be 1 number or several numbers, and the polynomials of those orders will be plotted for all x coordinates.}
    \item{setOrder<-}{\code{signature(obj = "OrthIm")}: Set the order up to which to calculate moments in the x and y directions. }
    \item{setParams<-}{\code{signature(obj = "OrthIm")}: Set the parameters to be used for calculating moments. }
    \item{transform<-}{\code{signature(obj = "OrthIm")}: Perform a transform of the image which creates an image of the pixels rearranged by polar coordinates. The rows are radius from the centroid and columns are angle. The method calls \code{\link{polarTransform}} with \code{padded=FALSE, scale=0} }
	 }
}
\author{
Allison Irvine
}

\seealso{
\code{\link{OrthMoments}},
\code{\link{OrthReconstruct}},
\code{\link{calcCentroid}},
\code{\link{plotMoment}},
\code{\link{displayImg}},
\code{\link{checkOrder}},
\code{\link{polarTransform}},
\code{\linkS4class{Image}},
\code{\linkS4class{CmplxIm}}
}
\examples{
data(mandril);
obj <- new("OrthIm",img = img,filename = "");
momentType(obj) = c("chebycont","legend");
setOrder(obj) = c(150,150);
Moments(obj) = NULL;
Reconstruct(obj) = c(125,125);
\dontrun{
	displayImg(list(obj@I,obj@moments,obj@reconstruction));
	plotPoly(obj,order=c(100,150));
}
}
\keyword{classes}