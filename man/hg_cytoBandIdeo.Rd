\name{hg_cytoBandIdeo}
\alias{hg_cytoBandIdeo}
\docType{data}
\title{
cytoBandIdeo human
}
\description{
Describes the positions of cytogenetic bands with a chromosome of human.
}
\usage{data(hg_cytoBandIdeo)}
\format{
  A data frame with 862 observations on the following 5 variables.
  \describe{
    \item{\code{Chrom}}{a character vector}
    \item{\code{Start}}{a numeric vector}
    \item{\code{End}}{a numeric vector}
    \item{\code{Name}}{a character vector}
    \item{\code{gieStain}}{a character vector}
  }
}
\details{
This file describes the cytogenetics positions on  chromosomes of human.
 specifically  it has 5 columns Chrom,  Start, End, Name y Group.
 Chrom  referenced to the chromosomes, the Start and End columns indicated
 the start and end positions, while that Name indicated the name of cytogenetics
 Bands, finally, Group column contains informations associated to the  Giemsa
 stain results.
}
\value{
data.frame that contain the positions of cytogenetic bands with a chromosome 
 of human.
}
\examples{
data(hg_cytoBandIdeo)
## maybe str(hg_cytoBandIdeo) ; plot(hg_cytoBandIdeo) ...
}
\keyword{datasets}
