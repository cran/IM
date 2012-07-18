\name{IM-package}
\alias{IM-package}
\alias{IM}
\docType{package}
\title{
Image Moment (IM) - Package
}
\description{
Compute Generalized-Pseudo Zernike, Fourier-Chebyshev, Fourier-Mellin, Radial Harmonic Fourier, Dual Hahn, Discrete Chebyshev, Krawtchouk, Gegenbauer, Legendre, Chebyshev moment type of single and multiple images.\cr \cr
Legend: gpzm, fc, fm, fr, hahn, cheby, krawt, gegen, legend, chebycont.
}
\details{
\tabular{ll}{
Package: \tab IM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-07-18\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
Depends: \tab R (>= 2.9.0),png (>= 0.1-4),jpeg (>= 0.1-2),bmp (>= 0.1)\cr
}
Class \code{linkS4class{Image}} parent class for CmplxIm and OrthIm classes\cr
Class \code{\linkS4class{CmplxIm}} handles moment analysis of single image with type : Generalized Pseudo-Zernike, Fourier-Chebyshev, Fourier-Mellin, Radial Harmonic Fourier.\cr
Class \code{\linkS4class{OrthIm}} handles moment analysis of single image with type: Dual Hahn, Discrete Chebyshev, Krawtchouk, Gegenbauer, Legendre, Chebyshev.\cr
Class \code{\linkS4class{MultiIm}} handles moments analysis of list of matrices with type: Generalized Pseudo-Zernike, Fourier-Chebyshev, Fourier-Mellin, Radial Harmonic Fourier, Dual Hahn, Discrete Chebyshev, Krawtchouk, Gegenbauer, Legendre, Chebyshev.\cr \cr
Function \code{\link{momentObj}} provides easy construction of object of different moment type and parameter.\cr\cr
Parameter requirements:
\tabular{ll}{
Generalized Pseudo-Zernike: \tab a is natural number. When a=0, the polynomials are equivalent to Pseudo-Zernike \cr
Dual Hahn: \tab a, c are integers, with a >= 0 and a > abs(c)-1 \cr
Krawtchouk: \tab a is a real number between 0 and 1 \cr
Gegenbauer: \tab a is real number and a > 1/2 \cr
}
See references for more details.
}
\author{
Allison Irvine, Tan Dang.
}
\references{
  Hosny, K.M.(2011). Image representation using accurate orthogonal Gegenbauer moments. \emph{Pattern Recognition Letters}, 32(6), 795--804.\cr

  Mukundan, R., Ong, S.H., and Lee, P.A.(2001). Image analysis by Tchebichef moments. \emph{Image Processing, IEEE Transactions on}, 10(9), 1357--1364.\cr

  Ping, Z., Wu, R., and Sheng, Y. (2002). Image description with Chebyshev-Fourier moments. \emph{Journal of the Optical Society of America A}, 19(9), 1748--1754.\cr

  Ren, H., Liu, A., Zou, J., Bai, D., and Ping, Z.(2007). Character Reconstruction with Radial-Harmonic-Fourier Moments. \emph{Fuzzy Systems and Knowledge Discovery}, 2007. FSKD 2007. \cr

  Sheng, Y., and Shen, L. (2007). Orthogonal Fourier-Mellin moments for invariant pattern recognition. \emph{Journal of the Optical Society of America A}, 11(6), 1748--1757.\cr

  Xia, T., Zhu, H., Shu, H., Haigron, P., and Luo, L. (2007). Image description with generalized pseudo-Zernike moments. \emph{Journal of the Optical Society of America A}, 24(1), 50--59.\cr

  Venkataramana, A., and Raj, P.A.(2007). Image Watermarking Using Krawtchouk Moments. \emph{Computing: Theory and Applications, 2007. ICCTA '07. International Conference on}, 5-7 March 2007, 676--680.\cr

  Zhu, H.(2012). Image representation using separable two-dimensional continuous and discrete orthogonal moments. \emph{Pattern Recognition}, 45(4), 1540--1558.\cr

  Zhu, H., Shu, H., Zhou, J., Luo, L. and Coatrieux, J.L.(2007) Image Analysis by Discrete Orthogonal Dual Hahn Moments. \emph{Pattern Recognition Letters}, 28(13), 1866--1704.\cr
}
\keyword{ package }