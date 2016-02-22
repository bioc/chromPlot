\name{mm10_gap}
\alias{mm10_gap}
\docType{data}
\title{
Gaps
}
\description{
This track depicts gaps in the assembly of the mouse genome.
}
\usage{data(mm10_gap)}
\format{
  A data frame with 686 observations on the following 4 variables.
  \describe{
    \item{\code{Chrom}}{a character vector}
    \item{\code{Start}}{a numeric vector}
    \item{\code{End}}{a numeric vector}
    \item{\code{Name}}{a character vector}
  }
}
\details{
This track depicts gaps in the assembly (Dec. 2011, Genome Reference Consortium
 Mouse Build 38 (GCA_000001635.2)) of the mouse genome.
 Gaps are represented as black boxes in this track. This assembly contains the 
 following principal types of gaps: (In this context, a contig is a set of
  overlapping sequence reads.)
 Clone - gaps between clones (114 gaps).
 Contig - gaps between map contigs, various sizes (104 gaps).
 Telomere - 42 gaps for telomeres (100,000 Ns)
 Centromere - 20 gaps for centromeres (size: 2,890,000 Ns)
 Short_arm - 21 gaps for the short arm (10,000 Ns) at base positions
  100,001-110,000 of each chromosome.
 other - sequence of Ns in the assembly that were not marked as gaps in the AGP 
 assembly definition file, various sizes (384 gaps).
 Fragment - a single gap of 31 bases in chrX_GL456233_random. 
}
\value{
data.frame that contain the gaps in the assembly of the mouse genome.
}
\examples{
data(mm10_gap)
## maybe str(mm10_gap) ; plot(mm10_gap) ...
}
\keyword{datasets}
