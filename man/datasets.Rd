\name{datasets}
\docType{data}
\alias{lena}
\alias{pirate}
\alias{mandril}
\alias{livingroom}
\alias{earth}
\alias{earth_200}
\alias{bacteria}
\alias{characters}
\alias{img}
\alias{labels}
\title{Image Sample Dataset}
\description{
       Sample images provided with this package are lena, pirate, mandril, livingroom, earth, and earth_200. 
       Two sets of images used to demonstrate classification are provided. One dataset, "bacteria", contains 
       light scatter images of 4 different types of bacteria. There are 5 images from each class. When the data is loaded, 
       the variable "img" contains the images, and "labels" contains the type of bacteria of each image in the image list.
       The "characters" dataset contains images of 8 different handwritten Tamil characters. There are 40 images for each class.
}
\usage{
   data(lena)
   data(pirate)
   data(mandril)
   data(livingroom)
   data(earth)
   data(earth_200)
   data(bacteria)
   data(characters)
}

\format{RBG image}
\source{http://www.imageprocessingplace.com/root_files_V3/image_databases.htm\cr
		 http://en.wikipedia.org/wiki/Earth\cr
		 Bartek Rajwa, Bindley Bioscience Center, Discovery Park, Purdue University\cr
		 http://www.hpl.hp.com/india/research/penhw-interfaces-1linguistics.html}
\references{
       Image Processing Place. Retrieved on 07/18/2012 from http://www.imageprocessingplace.com/root_files_V3/image_databases.htm\cr
       HP Labs India. Retrived on 08/13/2012 from http://www.hpl.hp.com/india/research/penhw-interfaces-1linguistics.html	
}
\keyword{datasets}