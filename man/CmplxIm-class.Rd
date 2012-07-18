\name{CmplxIm-class}
\Rdversion{1.1}
\docType{class}
\alias{CmplxIm-class}

\alias{Moments<-,CmplxIm-method}
\alias{Invariant<-,CmplxIm-method}
\alias{Reconstruct<-,CmplxIm-method}
\alias{initialize,CmplxIm-method}
\alias{momentType<-,CmplxIm-method}
\alias{setImage<-,CmplxIm-method}
\alias{setOrder<-,CmplxIm-method}
\alias{setParams<-,CmplxIm-method}
\alias{setPolar<-,CmplxIm-method}

\title{CmplxIm Class}
\description{
This class contains an image and all information about the image required for computing the moments 
of the image. The object has methods for computing moments, moment invariants, and reconstructing the image from moments.\cr
The types of moments that can be calculated using this object are: generalized Pseudo-Zernike, Fourier Mellin, Fourier Chebyshev, and Radial Harmonic Fourier.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CmplxIm", img, filename)}.
First the constructor checks for an input argument containing an image (img). This can be a multidimensional numeric array or a matrix. 
The image is converted to grayscale by summing over all colors and stored in slot "I". 
If an input image is not provided, a filename (filename) can be passed to the constructor.
The constructor will then read in the image file and convert it to grayscale. 
If a filename is given, it must have the extension ".jpg", ".png", or ".bmp".
If "img" and "filename" are not provided, and empty object will be created.
}
\section{Slots}{
  \describe{
    \item{\code{radius}:}{Object of class \code{"matrix"} Contains the radius of all of the pixels with respect to the centroid of the image. }
    \item{\code{theta}:}{Object of class \code{"matrix"} Contains the angle of all of the pixels with respect to the centroid. }
    \item{\code{constants}:}{Object of class \code{"matrix"} Contains constants used in the calculation of generalized pseudo-zernike moments and image reconstruction. }
    \item{\code{params}:}{Object of class \code{"numeric"} Contains the parameters used to calculate generalized pseudo-zernike moments. }
    \item{\code{order}:}{Object of class \code{"numeric"} Contains the order and repetition used to calculate moments. Two values may be given, which specify order and repetition, except in the case of generalized Pseudo-Zernike moments. In this case only one value should be provided, specifying order. Then repetition is constrained by the order. }
    \item{\code{I}:}{Object of class \code{"matrix"} Contains the original image, in grayscale. }
    \item{\code{dimensions}:}{Object of class \code{"numeric"} Contains the dimensions of the image. }
    \item{\code{centroid}:}{Object of class \code{"numeric"} Contains the x and y coordinates of the centroid (center point) of the image, 
    used for calculating continuous orthogonal moments. }
    \item{\code{filename}:}{Object of class \code{"character"} Contains the filename of the image. }
    \item{\code{imType}:}{Object of class \code{"character"} Contains the original image file type. }
    \item{\code{momentType}:}{Object of class \code{"character"} Contains the type of moments being calculated. 
    Choices are: "gpzm" (Generalized Pseudo-Zernike),"fm" (Fourier Mellin),"fc" (Fourier Chebyshev),"fr" (Harmonic Radial Fourier). }
    \item{\code{moments}:}{Object of class \code{"matrix"} Contains the moments calculated from the image. Dimensions are order by repetition. }
    \item{\code{invariant}:}{Object of class \code{"matrix"} Contains the invariants computed for the image. Dimensions are order in y direction by order in x direction. }
    \item{\code{reconstruction}:}{Object of class \code{"matrix"} Contains the image reconstructed from the moments. }
    \item{\code{error}:}{Object of class \code{"character"} containing a list of error messages produced from calling class methods. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{Image}"}, directly.
}
\section{Methods}{
  \describe{
    \item{Moments<-}{\code{signature(obj = "CmplxIm")}: Calculate the moments of the image. Before calling this method the momentType, params, and order must be set. Usage: Moments(obj) <- 0 }
    \item{Invariant<-}{\code{signature(obj = "CmplxIm")}: Calculate moment invariants of the image by taking the absolute values of the moments and normalizing them by the moment with order and repetition 0. Before calling this method the momentType, params, and order must be set. Usage: Invariant(obj) <- 0 }
    \item{Reconstruct<-}{\code{signature(obj = "CmplxIm")}: Reconstruct the image from the moments, using moments up to a specified order. If no order value is given, the image will be reconstructed from all available moments. Useage: Reconstruct(obj) <- order }
    \item{initialize}{\code{signature(.Object = "CmplxIm")}: Constructor. A call to new("CmplxIm",img,filename) calls this function and creates an CmplxIm object. The object slots I, centroid, radius, theta, and dimensions are set. (And filename and imType if a filename was provided) }
    \item{momentType<-}{\code{signature(obj = "CmplxIm")}: Set the moment type for the moments to be calculated. }
    \item{setImage<-}{\code{signature(obj = "CmplxIm")}: Set "I", the image in the object. Input can be a matrix or multidimensional numeric array. When this is called, the object slots centroid, radius, theta, and dimensions are set. }
    \item{setOrder<-}{\code{signature(obj = "CmplxIm")}: Set the order up to which to calculate moments in the x and y directions. }
    \item{setParams<-}{\code{signature(obj = "CmplxIm")}: Set the parameters to be used for calculating moments.}
    \item{setPolar<-}{\code{signature(obj = "CmplxIm")}: Calculate the polar coordinates of the pixels with respect to the centroid. The object slots "radius" and "theta" are set. }
	 }
}

\details{
	For generalized Pseudo-Zernike ("gpzm") moments, a single parameter is required, \deqn{a}. A larger value of \deqn{a} decreases the range of the polynomial values used for moment computation. 
}

\author{
Allison Irvine
}

\seealso{
\code{\link{calcCentroid}}, 
\code{\link{polarXY}}, 
\code{\link{plotMoment}}, 
\code{\link{displayImg}}, 
\code{\link{checkOrder}}, 
\code{\link{polarTransform}}, 
\code{\linkS4class{Image}}, 
\code{\linkS4class{OrthIm}}
}
\examples{



#load image
data(earth)
#initialize object
obj<- new("CmplxIm", img=img)
\dontrun{displayImg(obj@I)}
#set the moment type to generalized Pseudo-Zernike
momentType(obj)<- "gpzm"
#set the order
setOrder(obj)<- 25
#set the parameter  
setParams(obj)<- 1

\dontrun{
	#calculate moments of the image
	Moments(obj)<- NULL
	#calculate rotation invariants
	Invariant(obj) =NULL;
	#reconstruct the image from moments
	Reconstruct(obj)<- c(20,20)
}
\dontrun{plotMoment(obj);}
\dontrun{displayImg(obj@reconstruction);}
}
\keyword{classes}