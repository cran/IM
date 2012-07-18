\name{displayImg}
\alias{displayImg}
\alias{displayImg-methods}
\alias{displayImg,Numbers-method}
\alias{displayImg,list-method}
\title{
Display an image in grayscale
}
\description{
Display an image in grayscale. If the image is in color, it will still be displayed in grayscale. All images displayed will be processed with histogram equalization.
}
\usage{
displayImg(img)
}
\arguments{
  \item{img}{ A matrix or numeric array representation of the image. If the image is not in grayscale, it will still be displayed in grayscale. If using a CmplxIm, Image, or OrthIm object, 
  the original image would be represented as \code{obj@I} and the reconstructed image would be \code{obj@reconstruction}. 
  For MultiIm object, function takes in list of images. In this case the input argument would be \code{obj@imageList} to access the original images 
  or \code{obj@reconstruction}to access the reconstructed images. If the input is a list of images, they will all be plotted in one window.
}
}


\author{
Allison Irvine, Tan Dang
}


\seealso{
\code{\linkS4class{Image}},\code{\linkS4class{CmplxIm}},\code{\linkS4class{OrthIm}}, \code{\linkS4class{MultiIm}}
}
\examples{
#### display the original image in a "OrthIm" object.
data(pirate);

#create OrthIm object from image
obj <- new("OrthIm",img = img,filename = "");
momentType(obj) = "cheby";
setOrder(obj) = c(200,200);

#display original image
displayImg(obj@I);

\dontrun{
	#analyze image
	Moments(obj) = 0;
	Reconstruct(obj) = NULL;
	
	#display reconstructed image
	displayImg(obj@reconstruction)
}


\dontrun{
	#### display a list of images from a "MultiIm" object.
	data(characters)
	#take a small subset
	img = img[1:5];
	#analyze images
	obj= momentObj(I=img, type = "cheby", order = c(dim(img[[1]])[2],dim(img[[1]])[1]) );
	Reconstruct(obj)= NULL;
	
	#display original images
	displayImg(obj@imageList);
	
	#display reconstructed images
	displayImg(obj@reconstruction);
}
}
\keyword{misc}
