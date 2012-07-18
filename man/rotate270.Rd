\name{rotate270}
\alias{rotate270}
\alias{rotate270-methods}
\alias{rotate270,Numbers-method}

\title{
Rotate an image (or matrix) 270 degrees.
}
\description{
Rotate an image (or matrix) 270 degrees.
}

\usage{
rotate270(img)
}

\arguments{
	\item{img}{ A matrix or numeric array representation of an image. 
	}
}


\value{
	An image of equal dimensionality to the input, rotated 270 degrees.
}

\author{
	Allison Irvine
}


\seealso{
	\code{\link{displayImg}}
}

\examples{
	data(pirate);
	#perform histogram equalization
	img = histeq(img);
	#convert to grayscale
	img = rowSums(img, dims=2)/3;
	#rotate image 270 degrees so it appears upright
	img = rotate270(img);
	#set colors
	levels = seq(0,1,.0000001);
	g = gray(levels);
	image(img,col=g)
}
\keyword{misc}
