\name{plotMoment}
\alias{plotMoment}
\alias{plotMoment-methods}
\alias{plotMoment,CmplxIm-method}
\alias{plotMoment,OrthIm-method}
\alias{plotMoment,MultiIm-method}
\title{
Display a heat map of image moments
}
\description{
Displays image moments as a heat map. A inverse hyperboic sine transform is used to scale the values.
If the moments are complex numbers, the absolute value is plotted.
}
\usage{
plotMoment(obj)
}
\arguments{
  \item{obj}{
	A \code{Image}, \code{OrthIm}, \code{CmplxIm}, \code{MultiIm} class object
}
}
\author{
Allison Irvine
}
\seealso{
\code{\linkS4class{Image}}, \code{\linkS4class{CmplxIm}}, \code{\linkS4class{OrthIm}}, \code{\linkS4class{MultiIm}}, \code{\link{momentObj}}
}
\examples{
data(livingroom)
Obj= momentObj(img, type="gpzm", order=20, 0);
plotMoment(Obj);

\dontrun{
	images= list(img,img,img,img);
	Obj= momentObj(images, type="krawt", order=100, 0.5);
	plotMoment(Obj);
}
}
\keyword{misc}
\keyword{hplot}