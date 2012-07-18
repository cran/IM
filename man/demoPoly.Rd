\name{demoPoly}
\alias{demoPoly}
\alias{demoPoly-methods}
\alias{demoPoly,Numbers,character,numeric-method}
\title{
Calculate and display polynomials used to calculate image moments
}
\description{
Displays a plot of polynomials of a given order or plots of multiple orders. The polynomials are plotted
over a range of pixel coordinates. Discrete and continuous orthogonal polynomials are plotted for a given 
order over a range of a single dimension, and complex polynomials are plotted over 2 dimensions. 
}
\usage{
demoPoly(order, type, N, params)
}
\arguments{
  \item{order}{
The order of polynomials to be displayed. A single order or array of orders must be specified. 
}
  \item{type}{
The type of polynomial which will be displayed.
}
  \item{N}{
The maximum dimension up to which polynomials will be calculated. The polynomials displayed will be plotted from 0 to N along the x-axis. This argument represents 
the dimensionality of the image for which moments would be calculated.  
}
  \item{params}{
The parameters of the polynomials to be calculated. The parameters required depend on the polynomial type. 
}
}
\details{
For complex polynomials, (generalized Pseudo-Zernike,Fourier Chebyshev,Radial Harmonic Fourier, and Fourier Mellin), both order and repetition must be specified for each plot, 
for example as: \code{(order, repetition)}. To plot multiple orders and repetition for all complex polynomials, the input argument 
would be in the form of a 2-dimensional matrix where the first column is order and the second column is the corresponding repetition.\cr 
For real polynomials, (Chebyshev, Dual Hahn, Krawtchouk, Gegenbauer, Legendre, continuous Chebyshev), 
each polynomials will be plotted over a single dimension of pixel coordinates for a single set of parameters and a specified range of orders.
All of the polynomials in this case will be plotted in the same graph.
}

\value{
polynomials: the displayed polynomials will be returned as a list, with each element in the list being a polynomial
of a specific order. This value is only returned for real polynomials, (Chebyshev, Dual Hahn, Krawtchouk, Gegenbauer, Legendre, continuous Chebyshev).
}
\author{
Allison Irvine, Tan Dang
}
\seealso{
\code{\linkS4class{OrthIm}}, \code{\linkS4class{CmplxIm}}, \code{\linkS4class{MultiIm}} 
}
\examples{
#display chebyshev polynomials of orders 20 to 25 for dimensionality 256
demoPoly(order = 20:25,type = "cheby",N = 256)

\dontrun{
	#display generalized Pseudo-Zernike polynomials of order and repetition (3,2), (4,2), (5,3), (6,3) 
	#for pixel dimensions of 300x300, with paramater a=0
	demoPoly(rbind(c(3,2), c(4,2), c(5,3), c(6,3)), "gpzm", 300, 0);
}

\dontrun{
	#display Radial Fourier polynomials of order and repetition (3,2), (4,2), (5,3), (6,3)  
	#for pixel dimensions of 300x300,
	demoPoly(rbind(c(3,2), c(4,2), c(5,3), c(6,3)), "fr", 300);
}

}
\keyword{utilities}
\keyword{misc}
