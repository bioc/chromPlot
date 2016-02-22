\name{hg_gap}
\alias{hg_gap}
\docType{data}
\title{
Human Gap
}
\description{
This track depicts gaps in the assembly of the human genome.
}
\usage{data(hg_gap)}
\format{
  A data frame with 457 observations on the following 4 variables.
  \describe{
    \item{\code{Chrom}}{a character vector}
    \item{\code{Start}}{a numeric vector}
    \item{\code{End}}{a numeric vector}
    \item{\code{Name}}{a character vector}
  }
}
\details{
 This track depicts gaps in the assembly of human genome.
 Gaps are represented as black boxes in this track. This assembly contains
 the following principal types of gaps: 
 (In this context, a contig is a set of overlapping sequence reads.)
 Clone - gaps between clones (114 gaps).
 Contig - gaps between map contigs, various sizes (104 gaps).
 Telomere - 42 gaps for telomeres (100,000 Ns)
 Centromere - 20 gaps for centromeres (size: 2,890,000 Ns)
 Short_arm - 21 gaps for the short arm (10,000 Ns) at base positions
 100, 001-110,000 of each chromosome.
 other - sequence of Ns in the assembly that were not marked as gaps in the
 AGP assembly definition file, various sizes (384 gaps).
 Fragment - a single gap of 31 bases in chrX_GL456233_random. 
}
\value{
data.frame that contain the gaps in the assembly of the human genome.
}
\examples{
data(hg_gap)
## maybe str(hg_gap) ; plot(hg_gap) ...
}
\keyword{datasets}
