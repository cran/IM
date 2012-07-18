\name{histeq}
\alias{histeq}
\alias{histeq-methods}
\alias{histeq,Numbers-method}

\title{
	Perform histogram equalization on an image.
}
\description{
	Perform histogram equalization on an image.
}
\usage{
	histeq(I)
}
\arguments{
	\item{I}{ A matrix or numeric array representation of an image. 
	}
}

\value{
	An image of equal dimensions to the input, with histogram equalization applied to all color channels.
}


\author{
	 Tan Dang
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
	\dontrun{image(img,col=g)}
}
\keyword{misc}
