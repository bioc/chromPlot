\name{mm10_cytoBandIdeo}
\alias{mm10_cytoBandIdeo}
\docType{data}
\title{
cytoBandIdeo
}
\description{
 Describes the positions of cytogenetic bands with a chromosome of mouse.
}
\usage{data(mm10_cytoBandIdeo)}
\format{
  A data frame with 448 observations on the following 5 variables.
  \describe{
    \item{\code{Chrom}}{a character vector}
    \item{\code{Start}}{a numeric vector}
    \item{\code{End}}{a numeric vector}
    \item{\code{Name}}{a character vector}
    \item{\code{gieStain}}{a character vector}
  }
}
\details{
This file describes the cytogenetics positions on  chromosomes of mouse.
 specifically  it has 5 columns Chrom,  Start, End, Name y Group.
 Chrom  referenced to the chromosomes, the Start and End columns indicated
 the start and end positions, while that Name indicated the name of
 cytogenetics Bands, finally, Group column contains informations associated 
 to the  Giemsa stain results.
}
\value{
data.frame that contain the positions of cytogenetic bands with a chromosome
 of mouse.
}
\examples{
data(mm10_cytoBandIdeo)
## maybe str(mm10_cytoBandIdeo) ; plot(mm10_cytoBandIdeo) ...
}
\keyword{datasets}
