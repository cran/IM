\name{MultiIm-class}
\Rdversion{1.1}
\docType{class}
\alias{MultiIm-class}

\alias{initialize,MultiIm-method}
\alias{checkSize<-,MultiIm-method}
\alias{momentType<-,MultiIm-method}
\alias{setOrder<-,MultiIm-method}
\alias{setParams<-,MultiIm-method}
\alias{Moments<-,MultiIm-method}
\alias{Invariant<-,MultiIm-method}
\alias{plotPoly,MultiIm,Numbers-method}
\alias{Reconstruct<-,MultiIm-method}
\alias{setPoly<-,MultiIm-method}
\alias{transform<-,MultiIm-method}

\title{MultiIm Class}
\description{
The class contains a list of images and all information required for computing the moments of the images. 
The images must have the same dimensions and centroid. The object has methods for computing moments and 
reconstructing the image from moments. Types of moments that can be calculated for this class are: 
generalized pseudo-zernike, fourier mellin, fourier chebyshev, radial fourier, krawtchouk, hahn, 
continuous chebyshev, legendre, and gegenbauer.
}
\section{Creating objects}{
Objects can be created by calls of the form \code{new("MultiIm", images)} where \code{"images"} is list of matrices. 
The dimensions of the images will be checked for consistency when the constructor is called. Moment type, order, and parameters can 
be set using the class methods.
}
\section{Slots}{
  \describe{
    \item{\code{imageList}:}{Object of class \code{"list"} Contains a list of image matrices of the same dimension.}
    \item{\code{dimension}:}{Object of class \code{"numeric"} Contains the dimensions of the images. }
    \item{\code{polynomials}:}{Object of class \code{"list"} Contains the polynomial values used for computing moments.}
    \item{\code{storePoly}:}{Object of class \code{"logi"} Default is FALSE. Set to \code{"TRUE"} to store polynomial values. For gpzm, fm, fr, and fc moment types, if object@storePoly== TRUE, polynomial values will be stored in the slot "polynomials". 
	Polynomial values for other moment types will always be stored. Set directly by object@storePoly<- TRUE}
    \item{\code{momentType}:}{Object of class \code{"character"} Contains the moment type.}
    \item{\code{order}:}{Object of class \code{"numeric"} Contains the maximum order up to which to calculate moments. }
    \item{\code{Params}:}{Object of class \code{"list"} Contains the parameters used to calculate the moments. First list element is the parameters for the moments over x coordinates, second is for the y coordinates. }
    \item{\code{moments}:}{Object of class \code{"list"} Contains a list of the moments calculated for each image. }
    \item{\code{invariant}:}{Object of class \code{"list"} Contains a list of the invariants calculated for each image. }
    \item{\code{reconstruction}:}{Object of class \code{"list"} Contains the list of reconstructed images.}
    \item{\code{error}:}{Object of class \code{"character"} Contains the last error message produced from calling a class method. }
  }
}
\section{Methods}{
  \describe{
    \item{momentType<-}{\code{signature(obj = "MultiIm")}: Set the moment type for the moments to be calculated. 
    Takes one character string or an array of 2 character strings, for example c("type1","type2"). 
    Usage: momentType(obj)<- c("krawt", "hahn")  or momentType(obj)<-"gpzm"}
 \item{setOrder<-}{\code{signature(obj = "MultiIm")}: Set the orders up to which to calculate moments in the x and y 
 directions in the case of real orthogonal moments, or specify order and repetition for complex moment types. Usage: setOrder(obj)<- c(50,50)}
    \item{setParams<-}{\code{signature(obj = "MultiIm")}: Set the parameters to be used for calculating moments. 
    If two types of moments are being used, this should be a list of the parameters to be used for each moment type 
    being used. Usage: setParams(obj)<- list(c(0.5), c(1,1))}
    \item{transform<-}{\code{signature(obj = "MultiIm")}: Perform a transform of the images which creates an image 
    of the pixels rearranged by polar coordinates. The rows are radius from the centroid and columns are angle. See \code{\link{polarTransform}}.}
    \item{Moments<-}{\code{signature(obj = "MultiIm")}: Calculate the moments of a list of images. Before calling 
    this method the momentType, params, and order must be set. Usage: Moments(obj)<- 0. }
    \item{Invariant<-}{\code{signature(obj = "MultiIm")}: Calculate the moment invariants w.r.t. rotation of a list 
    of images. Before calling this method the momentType, params, and order must be set. The method assumes that all 
    images have the same centroid location. For complex moments, the magnitude of the moments are invariant to rotation. 
    Usage: Invariant(obj)<- NULL. \cr
    For real orthogonal moments (continuous or discrete), the method normalizes all images prior to computing moments 
    using the same procedure as \code{\link{Invariant<-,OrthIm-method}}. Usage: Invariant(obj)<- c(5,120) }
    \item{plotPoly}{\code{signature(obj = "MultiIm", order = "list")}: Plot the polynomials used to calculate moments 
    of a specified order or array of orders. Usage: plotPoly(obj, list(3:6, 1:4)).}
    \item{Reconstruct<-}{\code{signature(obj = "MultiIm")}: Reconstruct the images from moments. 
    Any order lower than or equal to the maximum order calculated can be used for reconstruction. Usage: Reconstruct(obj)<- c(50,50) }
	 }
}

\author{
Tan Dang
}

\seealso{
\code{\link{plotPoly}}
\code{\link{displayImg}}
\code{\linkS4class{Image}}, 
\code{\linkS4class{OrthIm}}
\code{\linkS4class{CmplxIm}}, 
}
\examples{

# Load sample data
data(bacteria);

#create MultiIm object, calculate generalized Pseudo-Zernike moments, invariants, and reconstruction for all images
obj <- new("MultiIm", img);
momentType(obj)= c("gpzm");
setOrder(obj)= 20;
setParams(obj)= 1;	# <--- GPZM has 1 parameter
obj@storePoly= TRUE;

\dontrun{
	Moments(obj)= NULL;
	Invariant(obj)= NULL;
	Reconstruct(obj)= NULL;
	#display a subset of polynomials stored in the object
	plotPoly(obj, order=rbind(c(3,2), c(4,2), c(5,3), c(6,3)));
	displayImg(obj@reconstruction[1:5]);
}

\dontrun{
	#create MultiIm object, calculate discrete Chebyshev moments, invariants, and reconstruction for all images
	obj <- new("MultiIm", img);
	momentType(obj)= "cheby"
	setOrder(obj)= 100;
	Moments(obj)= 0;
	Reconstruct(obj)= 75;
	Invariant(obj)= c(7,100); # resolution is 7, scale is 100
	displayImg(obj@reconstruction[1:4]);
}

\dontrun{
	#create MultiIm object, calculate continuous Chebyshev - Legendre moments, invariants, and reconstruction for all images
	obj <- new("MultiIm", img);
	momentType(obj)= c("chebycont", "legend");
	setOrder(obj)= c(100, 110);	
	transform(obj)= 8; 	# perform polar transform on images
	Moments(obj)= 0;
	Reconstruct(obj)= NULL;

	# order of pairs (20,10), (21, 11),...(28,18) will be plotted
	plotPoly(obj, order=cbind(20:28, 10:18));
	displayImg(obj@reconstruction[1:6]);
}
}
\keyword{classes}
