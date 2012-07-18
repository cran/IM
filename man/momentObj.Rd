\name{momentObj}
\alias{momentObj}
\alias{momentObj-methods}
\alias{momentObj,list,character,Numbers-method}
\alias{momentObj,Numbers,character,Numbers-method}
\title{
Calculate moments of image 
}
\description{
Computes continuous complex or discrete/continuous orthogonal moments for an image or list of images.}
\usage{
momentObj(I, type, order, paramsX, paramsY)
}

\arguments{
  \item{I}{
	Either an image in matrix representation or a list of images.
}
  \item{type}{
	One string or an array of two strings indicate the moment types to be calculated. 
	Two types can be given only if they are both either continuous or discrete orthogonal moment types (not complex types).
	In this case, the moments in the x direction will be of first type and the moments in the y direction will be of the second type.
}
  \item{order}{
	In the case of real orthogonal moments, one or two numbers may be provided to indicate the orders up to which to calculate moments in the x and y directions of the image.
	If one value is provided, the order will be the same in the x and y directions. In the case of complex moments (except for generalized Pseudo-Zernike), 
	two values should be provided to specify the order and repetition. If generalized Pseudo-Zernike moments are being calculated, repetition is constrained by order,
	so only one value should be specified to indicate order.
}
  \item{paramsX}{
  	Parameters required for calculation of moments in the x direction. These are only required for certain types of moments:
  	generalized Pseudo-Zernike, Hahn, Krawtchouk, and Gegenbauer. If the same moment types and parameters are to be used to calculate moments in 
  	both the x and y directions, the parameters may be specified only by this argument.
}
  \item{paramsY}{
	Parameters required for calculation of moments in the y direction. These are only required for certain types of moments:
  	generalized Pseudo-Zernike, Hahn, Krawtchouk, and Gegenbauer. 
}
}

\value{
If only one image is provided to the function in the first argument "I", an object of type \code{"CmplxIm"} or \code{"OrthIm"} will be returned, depending on the moment type specified. \cr
If the argument "type", specifying the moment type to be calculated, is "gpzm","fm","fc",or "fr", an object of type \code{"CmplxIm"} will be returned. \cr
If the argument "type", specifying the moment type to be calculated, is "cheby","chebycont","hahn","krawtchouk","gegen", or "legend", an object of type \code{"OrthIm"} will be returned.\cr 
If a list of images is provided to the function in the first argument "I", an object of type \code{"MultiIm"} will be returned.
}

\details{
	For generalized Pseudo-Zernike ("gpzm") moments, a single parameter is required, \deqn{a}. A larger value of \deqn{a} decreases the range of the polynomial values used for moment computation. \cr
	For Gegenbauer moments, a single parameter is required, \deqn{a} with the constraints \deqn{a > 1} and \deqn{a ~= 0}.A larger value of \deqn{a} increases the values of the polynomials used for 
	moment calculation. If \deqn{a = 1}, the moments are equivalent to continuous Chebyshev moments.\cr
	For Hahn moments, two parameters are required, \deqn{a} and \deqn{c} with the constraints \deqn{a > -\frac{1}{2}} and \deqn{a > abs(c)-1}.
	\deqn{a} specifies where the polynomials are centered and the difference between \deqn{a} and \deqn{c} is positively correlated to the range of the polynomial values used to calculate moments.\cr
	For Krawtchouk moments, a single parameter is required, \deqn{a} with the constraint \deqn{0 < a < 1}.\deqn{a} specifies where the polynomials are centered, relative to the image centroid.\cr
	See references for further details on these parameters.\cr
	See the documentation for \code{\linkS4class{OrthIm}},\code{\linkS4class{CmplxIm}},and \code{\linkS4class{MultiIm}} for details about image reconstruction.
}

\author{
Allison Irvine, Tan Dang
}
\seealso{
\code{\linkS4class{OrthIm}},\code{\linkS4class{CmplxIm}},\code{\linkS4class{MultiIm}}
}

\examples{
#compute chebyshev moments of an image 
data(lena)
\dontrun{displayImg(obj@I)}
obj = momentObj(I=img,type="cheby",order=c(500,500));
\dontrun{plotMoment(obj)}
#reconstruct the image from the moments
Reconstruct(obj) = c(200,200);
\dontrun{displayImg(obj@reconstruction)}


\dontrun{
	#calculate bivariate Legendre/Gegenbauer moments up to orders 50 and 100 with parameter 1 for Gegenbauer moments.
	data(lena);
	displayImg(obj@I);
	obj= momentObj(img,c("legend", "gegen"),c(200, 300),NULL,2);
	plotMoment(obj);
	#reconstruct the image from all the moments
	Reconstruct(obj) = NULL;
	displayImg(obj@reconstruction);
}


\dontrun{
	#compute Radial Harmonic-Fourier moment invariants of a list of images, using up to order 10 and repetition 10
	data(bacteria);
	obj= momentObj(I=img[1:10],type="fr",order=c(10,10));  
	Invariant(obj) = NULL;
}
}


\keyword{misc}