\name{checkOrder}
\alias{checkOrder}
\alias{checkOrder-methods}
\alias{checkOrder,matrix-method}
\title{
Find the maximum order of finite moments
}
\description{
Function to find the maximum order (and repetition for complex moments) of moments such that all moment values up to that order (in both the x and y directions) are finite and not Nan values. 
The maximum order will be at most the maximum order represented in the matrix of moment values. 
This function looks for the largest size rectangular block in the moment matrix such that no values are infinity or NaN.
}

\usage{
checkOrder(moments)
}
\arguments{
  \item{moments}{
Matrix containing image moments from order 0 up to some order \code{p} in the x-direction and \code{q} in the y-direction, 
or in the case of complex moments, up to order \code{p} and repetition \code{q}.
}
}
\value{
\item{orderP}{The maximum order along the x-axis of the moments matrix, which is also the maximum number of valid columns minus 1.}
\item{orderQ}{The maximum order along the y-axis of the moments matrix, which is also the maximum number of valid rows minus 1.}
}

\author{
Allison Irvine
}
\note{
The maximum order of \code{p} and \code{q} returned will be the same in the current implementation, representing the upper bounds of a square block within the moments matrix.
}

\seealso{
\code{\linkS4class{OrthIm}}
}
\examples{
data(lena);

obj= momentObj(img, "cheby", c(300,200));

maxOrder= checkOrder(obj@moments);
}
\keyword{misc}