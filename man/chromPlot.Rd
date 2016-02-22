\name{chromPlot}
\alias{chromPlot}
\title{Global visualization tool of genomic data}
\description{
Package designed to visualize genomic data along the chromosomes,
where the vertical chromosomes are sorted by number, with sex chromosomes
at the end.\cr
}
\usage{
chromPlot(annot1, annot2, annot3, annot4, stat, stat2, 
scale.title="Counts", statTyp="p", scex=1, spty=20, statCol, statCol2,
statName="Statistic", statName2="Statistic2", bands, bandsDesc, gaps,
gapsDesc, segment, segmentDesc, segment2=NULL, segment2Desc=NULL, chr,
bin=1e6, yAxis=TRUE, figCols=NULL, colBand="lightgray", colAnnot1="brown", 
colAnnot2="gold", colAnnot3="darkgreen", colAnnot4="blue", colSegments=c("darkgreen",
"orange", "blue",  "darkslategray2", "cyan", "blueviolet", "goldenrod3", "darkseagreen4",
"red", "green", "salmon", "darkolivegreen", "maroon", "purple"), 
colSegments2=colSegments[-1L], colStat="blue", colStat2="orange", title=NULL,
plotRndchr=FALSE, maxSegs=200, noHist=FALSE, segLwd=3, sortSegs=TRUE,
chrSide=c(-1, -1, -1, -1, 1, 1, -1, 1), cex=0.75, legChrom, org=NULL, strand=NULL,
stack=TRUE, statThreshold=NULL, statThreshold2=NULL, statSumm="none")
}
\arguments{
  \item{annot1}{Genome annotations}
  \item{annot2}{Genome annotations, subset of annot1}
  \item{annot3}{Genome annotations, subset of annot2}
  \item{annot4}{Genome annotations, subset of annot3}
  \item{stat}{Genome annotations associated to quantitative values}
  \item{stat2}{Second track of genome annotations associated to
quantitative values}
  \item{statCol}{Name column in stat with the values to plot}
  \item{statCol2}{Name column in stat2 with the values to plot}
  \item{statTyp}{Type of plot for stat ("l", "p", NULL)}
  \item{statName}{Description for stat (default="Statistic")}
  \item{statName2}{Description for stat2 (default="Statistic")}
  \item{bands}{Genome annotations to be plotted on chromosomal body
(e.g G bands) }
  \item{bandsDesc}{Description for bands}
  \item{gaps}{Chromosome alignment gaps (only centromers and telomers
used)}
  \item{gapsDesc}{Description for gaps}
  \item{segment}{Genomic segments. Can contain a 'Group' column with
categories}
  \item{segmentDesc}{Description for segment}
  \item{segment2}{second track of genomic segments. Can contain a
'Group' column with categories}
  \item{segment2Desc}{Description for segment2}
  \item{chr}{Vector of chromosome names to plotted (optional)}
  \item{bin}{Bin size for histograms in base pairs}
  \item{yAxis}{Should I draw the y-axis (logical)}
  \item{figCols}{Maximum number of chromosomes in a row}
  \item{colBand}{Color for chromosome bands}
  \item{colAnnot1}{Color for histograms for annot1}
  \item{colAnnot2}{Color for histograms for annot2}
  \item{colAnnot3}{Color for histograms for annot3}
  \item{colAnnot4}{Color for histograms for annot4}
  \item{colSegments}{Color for chromosome segment (ignored if segment
are grouped (see details)}
  \item{colSegments2}{Color for chromosome segment2 (ignored if
segment2 are grouped (See details)}
  \item{colStat}{Color for stat}
  \item{colStat2}{Color for stat2}
  \item{title}{Plot title}
  \item{plotRndchr}{Include random scaffolds}
  \item{maxSegs}{Maximum number of segments. If the segment or
segment2 tracks contain more segments than this value, a histogram of
segments is drawn instead}
  \item{noHist}{If TRUE, segments are never drawn as histograms, even
they are more than maxSegs or if the largest segment is smaller than
the bin size.}
  \item{segLwd}{Line width for segments}
  \item{sortSegs}{Sort overlapping segments by size}
  \item{chrSide}{Chromosome side where to draw annot1, annot2, annot3,
annot4, segments, segments2, stat and stat2, respectively. 1=right,
-1=left}
  \item{cex}{Cex for plot (see ?par for details)}
  \item{legChrom}{Legend chromosome (character string). Place legend
after this chromosome}
  \item{scale.title}{Title for histograms scales}
  \item{scex}{Cex for stat track}
  \item{spty}{A character specifying the type of plot region to be
used in stat}
  \item{org}{Organism name, e.g. mmusculus, hsapiens }
  \item{strand}{Strand "+" or "-" for local view using GenomeGraphs}
  \item{stack}{Stack overlapping segments in segment and segment2 in
clusters}
  \item{statThreshold}{Only plot segments in stat with values above
this threshold}
  \item{statThreshold2}{Only plot segments in stat2 with values
above this threshold}
 \item{statSumm}{Type of statistical function for apply to the data ("mean", "median","sum", "none"),
if the value is 'none', chromPlot will not apply some statistical function.}
}
\details{
chromPlot package creates an idiogram with all chromosomes including
the sex chromosomes. The package is able to plot genomic data on both
sides of chromosome as histograms or vertical segments. Histograms represent
the number of genomic elements in each bin of size bin. The parameters annot1, annot2,
annot3, annot4, segment, segment2, stat, stat2, band,
gaps should be data.frames with at leas these columns: 'Chrom', 'Start', 'End'. The
gaps and bands arguments are used to plot the chromosomal ideogram.The argument
band should also have a 'Group' column with categories for classifying
each annotation element. Arguments stat and stat2 should have a statCol and stat2Col
column respectively with continuoos values.

If plotted on the same chromosomal side, tracks will be plotted on top of each other, in
the order they are in the function's syntax. This can be used for plotting stacked barplots
if, for instance, annot1, annot2, annot3, and annot4 are supersets of ecah other.
This, however, is not enforced nor checked. An alternative way to create a stacked histogram is
providing a single track with Group category. The user can modify the side tracks are
plotted on by modifying chrSide.

The segment and segment2 tracks are plotted as vertical bars by default. However, the
their elements exceed in number given to maxSegs or if the maximum segment size is smaller
than bin, they are plotted as histograms. This behaviour can be modified by setting
noHist = TRUE.

For more details and usage examples see the vignette.}
\value{
Karyotype diagram in device.
}
\author{Ricardo Verdugo and Karen Orostica}
\examples{
data(hg_cytoBandIdeo)
data(hg_gap)
chromPlot(bands=hg_cytoBandIdeo, gaps=hg_gap)


}
\keyword{Histograms}

