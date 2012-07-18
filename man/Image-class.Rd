\name{Image-class}
\Rdversion{1.1}
\docType{class}
\alias{Image-class}
\alias{initialize,Image-method}
\alias{plotMoment,Image-method}
\alias{setCentroid<-,Image-method}
\alias{setImage<-,Image-method}

\title{Basic class to handle image moment analysis}
\description{
This class contains general information about an image, and is inherited by the \code{"CmplxIm"} and \code{"OrthIm"} classes.
It contains basic utilities and object slots which are needed for both continuous complex moment analysis and
discrete/continuous real orthogonal moment analysis. An object of this class can be cast to either a CmplxIm type object
or a OrthIm type object. 
}

\section{Creating objects}{
Objects will be created by calls of the form \code{new("Image", img, filename)}. When the constructor is called and a filename is given, 
the input file will be loaded and the image will be converted to grayscale and stored in slot "I". Only files of type "jpeg", "bmp" or "png" may be imported this way. If the argument "img", an image in matrix or numeric form, is provided
when the constructor is called, the image will be converted to grayscale and stored in slot "I" of the object. When this slot value is set, the centroid and dimensions of the image will be calculated. 
The constructor for this class will be called when objects of the inheriting classes \code{"CmplxIm"} and \code{"OrthIm"} are created. See \code{\linkS4class{OrthIm}}, \code{\linkS4class{CmplxIm}}.
}

\section{Slots}{
  \describe{
    \item{\code{I}:}{Object of class \code{"matrix"} representing the image.}
    \item{\code{dimensions}:}{Object of class \code{"numeric"} containing dimension of the image.}
    \item{\code{centroid}:}{Object of class \code{"numeric"} containing centroid of the image. See \code{\link{calcCentroid}}.}
    \item{\code{filename}:}{Object of class \code{"character"} containing name of the image}
    \item{\code{imType}:}{Object of class \code{"character"} containing image extension. }
    \item{\code{momentType}:}{Object of class \code{"character"} containing moment type. Either \code{"gpzm", "fm", "fr", "fc", "krawt", "hahn", "cheby", "gegen", "legend", "chebycont"}. See \code{\linkS4class{OrthIm}}, \code{\linkS4class{CmplxIm}}.}
    \item{\code{moments}:}{Object of class \code{"matrix"} containing computed moments from the image.}
    \item{\code{invariant}:}{Object of class \code{"matrix"} containing computed invariants from the image.}
    \item{\code{reconstruction}:}{Object of class \code{"matrix"} representing the reconstructed image from the computed moment.}
    \item{\code{error}:}{Object of class \code{"character"} containing a list of error messages produced from calling class methods.}
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "Image")}: Construct the \code{"Image"} object. To be used internally when constructing \code{"OrthIm"} and \code{"CmplxIm"} object only. }
    \item{plotMoment}{\code{signature(obj = "Image")}: Plot the calculated moments. Usage: plotMoment(obj).}
    \item{setCentroid<-}{\code{signature(obj = "Image")}: Set centroid of image to an array of 2 values (the x and y pixel coordinates of the centroid). If a value is not supplied, the centroid will be calculated using geometric moments. Usage: setCentroid(obj)<- c(value1,value2). }
    \item{setImage<-}{\code{signature(obj = "Image")}: Set the image slot of the object to an image in matrix or numeric form. Usage: setImage(obj)<- img.}
	 }
}
\note{
This class does not have functions available for calculating moments or reconstruction from moments. 
In order to calculate moments you must use an object of type "OrthIm" or "CmplxIm".
}
\author{
Allison Irvine
}
\seealso{
\code{\linkS4class{OrthIm}}, 
\code{\linkS4class{CmplxIm}}, 
\code{\link{calcCentroid}}, 
\code{\link{plotMoment}}
}
\examples{
# create an object of class "Image"
data(lena)
obj = new("Image",img=img)
#display image
\dontrun{displayImg(obj@I)}
#convert the object to the class "OrthIm"
obj=as(obj,"OrthIm")
}
\keyword{classes}
