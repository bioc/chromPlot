chromPlot <- function(annot1, annot2, annot3, annot4, stat, stat2, 
scale.title="Counts", statTyp="p", scex=1, spty=20, statCol, statCol2,
statName="Statistic", statName2="Statistic2", bands, bandsDesc, gaps,
gapsDesc, segment, segmentDesc, segment2=NULL, segment2Desc=NULL, chr,
bin=1e6, yAxis=TRUE, figCols=NULL, colBand="lightgray", colAnnot1="brown", 
colAnnot2="gold", colAnnot3="darkgreen", colAnnot4="blue",
colSegments=c("darkgreen","orange", "blue",  "darkslategray2", "cyan", 
"blueviolet", "goldenrod3", "darkseagreen4","red", "green", "salmon",
"darkolivegreen", "maroon", "purple"), colSegments2=colSegments[-1L],
colStat="blue", colStat2="orange", title=NULL, plotRndchr=FALSE, 
maxSegs=200, noHist=FALSE, segLwd=3, sortSegs=TRUE, 
chrSide=c(-1, -1, -1, -1, 1, 1, -1, 1), cex=0.75, legChrom, org=NULL,
strand=NULL, stack=TRUE, statThreshold=NULL, statThreshold2=NULL,
statSumm="none") {
## chromPlot          
##
## Chromosome summary plot
## ---------- ------- ----
##
## Description:
## -----------
##   Karyotype diagram with genomic elements in the complete genome.
##
## Usage:
## -----
##   chromPlot(annot1, annot2, annot3, annot4, stat, statTyp, scex = 1,
## statCol, bands, bandsDesc, gaps, gapsDesc, segment, segmentDesc, chr,
 ##bin = 1e6, yAxis = T, figCols = 11, colBand = "##9e9f93", 
##colAnnot1 = "brown4", colAnnot2 = "gold1", colAnnot3 = "darkgreen", 
##colAnnot4 = "blue", colSegments = "darkblue", colStat = "blue",
##title = "Chromosome Plot", plotRndchr=FALSE, legChrom)
##
## Arguments:
## ---------
## annot1          genome annotations
## annot2          genome annotations, subset of annot1
## annot3          genome annotations, subset of annot2
## annot4          genome annotations, subset of annot3
## stat            genome annotations associated to quantitative values
## stat2           second track of genome annotations associated to quantitative
##                 values
## statCol         name column in stat with the values to plot
## statCol2        name column in stat2 with the values to plot
## statTyp         type of plot for stat ("p", "l", "hst", "seg", NULL)
## statName        description for stat (default="Statistic")
## statName2       description for stat2 (default="Statistic")
## bands           genome annotations to be plotted on chromosomal body 
##                 (e.g G bands) 
## bandsDesc       description for bands
## gaps            chromosome alignment gaps (only centromers and telomers used)
## gapsDesc        description for gaps
## segment         genomic segments. Can contain a 'Group'column with categories
## segmentDesc     description for segment
## segment2        second track of genomic segments. Can contain a 'Group'
##                 column with categories
## segment2Desc    description for segment2
## chr             vector of chromosome names to plotted (optional)
## bin             bin size for histograms in base pairs
## yAxis           should I draw the y-axis (logical)
## figCols         maximum number of chromosomes in a row
## colBand         color for chromosome bands
## colAnnot1       color for histograms for annot1
## colAnnot2       color for histograms for annot2
## colAnnot3       color for histograms for annot3
## colAnnot4       color for histograms for annot4
## colSegments     color for chromosome segment (ignored if segment are
##                 grouped (see details)
## colSegments2    Color for chromosome segment2 (ignored if segment2 are
##                 grouped (see details)
## colStat         color for stat
## colStat2        color for stat2
## title           plot title
## plotRndchr      include random scaffolds
## maxSegs         maximum number of segments. If the segment or segment2 
##                 tracks contain more segments than this value, a histogram of 
##                 segments is drawn instead
## noHist          If TRUE, segments are never drawn as histograms, even they
##                 are more than maxSegs or if the largest segment is smaller
##                 than the bin size.
## segLwd          Line width for segments
## sortSegs        Sort overlapping segments by size
## chrSide         Chromosome side where to draw annot1,annot2, annot3,
##                 annot4, Segments, segments2, stat, stat2, 
##                 respectively. 1=right, -1=left.
## cex             cex for plot (see ?par for details)
## legChrom        legend chromosome (character string). Place legend after this
##                 chromosome
## scex            cex for stat track
## spty            a character specifying the type of plot region to be used in 
##                 stat
## org             organism name, e.g. mmusculus, hsapiens 
## strand          Strand "+" or "-" for local view using GenomeGraphs
## stack          stack overlapping segments in segment and segment2 in clusters
## statThreshold   only plot segments in stat with values above this threshold
## statThreshold2  only plot segments in stat2 with values above this threshold
## statSumm        Type of statistical function for apply to the data ("mean",
##                 "median","sum", "none"), if the value is 'none', chromPlot
##                  will not apply some statistical function.   
##
## Details:
## -------
## chromPlot package creates an idiogram with all chromosomes including
## the sex chromosomes. The package is able to plot genomic data on both
## sides of chromosome as histograms or vertical segments. Histograms represent
## the number of genomic elements in each bin of size 'bin'.The parameters 
##'annot1','annot2', 'annot3', 'annot4', 'segment', 'segment2'', 'stat', 
##'stat2', 'band', 'gaps' should be data.frames with at leas these columns:
##'Chrom', 'Start', 'End'. The 'gaps' and 'bands' arguments are used to plot the
## chromosomal ideogram. Arguments and 'band' should also have a 'Group' column
## with categories for classifying each annotation element. Arguments 'stat' and
## 'stat2' should have a statCol and stat2Col column respectively with
## continuous values.
##
## If plotted on the same chromosomal side, tracks will be plotted on top of
## each other, in the order they are in the function's syntax. This can be used
## for plotting stacked barplots if, for instance, annot1, annot2, annot3, and 
## annot4 are supersets of each other. This, however, is not enforced or 
## checked. An alternative way to create stacked histogram is 
## providing a single track to segment or segment2 with a Group categorical
## variable. If a Group2 column is present, it is used for doble-categorization
## (only supported for plotting as segments, not histograms).
## The user can modify the side tracks are plotted on by modifying chrSide.
##
## The 'segment' and 'segment2' tracks are plotted as vertical bars by default.
## However, the their elements exceed in number given to maxSegs or if the
## maximum segment size is smaller than bin, they are plotted as histograms. 
##This behaviour can be modified by setting noHist = TRUE.
##
## For more details and usage examples see the vignette.

## Value:
## -----
##   A text string with a figure caption (invisible if not assigned to an 
##object).  
 
  ## Declaration of constants
    imgscl           <- cex
    upper            <- 0.1
    marAxis          <- 0.1
    margin           <- 0.1
    segColors        <- NULL
    minGeneCount     <- 0
    maxGeneCount     <- 1 
    plot_range       <- c(margin, maxGeneCount)
    flagBands        <- 0
    pos              <- 0
    Chrom            <- "Chrom"
    Start            <- "Start"
    End              <- "End"
    PointPosition    <- "Mid"
    Position         <- c(PointPosition, Start)
    segmentCols      <- c(Chrom, Start, End)
    statCols         <- c(Chrom, Start)
    CategoryName     <- "Group"
    Category2Name    <- "Group2"
    Name             <- "Name"
    ID               <- "ID"
    Colors           <- "Colors"
    aligncolors      <- colors()[c(8, 12, 26, 31, 33, 68, 125, 90, 258, 640,
                                598, 652, 562, 259, 619, 153, 79, 200, 526, 128,
                                104, 116, 382, 383, 635, 646)]
    colChrom         <- "lightgray"
                                  
    ## check missing parameters
    if(missing(annot1)) 
    annot1      <- NULL  
    if(missing(annot2))
    annot2   <- NULL    
    if(missing(annot3))
    annot3     <- NULL
    if(missing(annot4))
    annot4     <- NULL
    if(missing(stat))
    stat      <- NULL     
    if(missing(statCol)) {
        if(!is.null(stat)) {
            stop("You must provide statCol if stat is not null.")
        } else {
            statCol    <- NULL
        }
    }
    if(missing(stat2))
    stat2     <- NULL     
    if(missing(statCol2)) {
        if(!is.null(stat2)) {
            stop("You must provide statCol2 if stat2 is not null.")
        } else {
            statCol2   <- NULL
        }
    }
    if(missing(bands))
    bands      <- NULL  
    if(missing(bandsDesc) ){
        if(!is.null(bands)){
            bandsDesc   <- deparse(substitute(bands))
        }
        else  bandsDesc <- NULL
    }  
    if(missing(chr))
    chr        <- NULL 
    if(missing(gaps))
    gaps       <- NULL   
    if(missing(gapsDesc) ){
        if(!is.null(gaps)){
            gapsDesc    <- deparse(substitute(gaps))
        }
        else  gapsDesc  <- NULL
    }
    if(missing(segment)) 
    segment         <- NULL                
    if(missing(segmentDesc)) {
        segmentDesc     <- "Segments"
    } 
    if(missing(segment2Desc)) {
        segment2Desc    <- "Segments 2"
    }
    if(missing(legChrom))
    legChrom        <- NULL
    if(missing(org))
    org             <- NULL
    if(missing(statThreshold)) 
    statThreshold  <- NULL
    if(missing(statThreshold2)) 
    statThreshold2 <- NULL
    if(is.null(bandsDesc))
    bandsDesc<-"you forgot to provide a description of the chromosome bands file" 
    
    # Define objects
    chr.length       <- NULL
    colScaleStat     <- NULL
    colScaleStat2    <- NULL
    centroms         <- NULL
    dataBarplot2     <- NULL
    annot1Hist       <- NULL
    annot2Hist       <- NULL
    annot3Hist       <- NULL
    annot4Hist       <- NULL
    chrbands         <- NULL
    chrsegs          <- NULL
    chrsegs2         <- NULL
    trackStat        <- NULL
    trackStat2       <- NULL
    realMaxStat      <- NULL 
    realMaxStat2     <- NULL
    annot1_plot_range <- NULL
        
    # Sanity checks
    if(is.null(annot1) & is.null(bands) & is.null(gaps) & is.null(org)) {
        stop("Enter at least one of the following parameteres for plotting:",
        " bands, gaps, org, annot1.") 
    }
    if(is.null(org)) {
        if(is.null(annot1) & !is.null(annot2)) {
            stop("An object must be provided to annot1 if annot2 is provided.")
        }
        if(is.null(annot2) & !is.null(annot3)) {
            stop("An object must be provided to annot2 if annot3 is provided.")
        }
        if(is.null(annot3) & !is.null(annot4)) {
            stop("An object must be provided to annot3 if annot4 is provided.")
        }
    }
    if(!is.null(stat) & is.null(statCol)) {
        stop("The column name with statistics values in the stat object must",
        " be specified in the statCol argument.")      
    }
    if(!is.null(stat2) & is.null(statCol2)) {
        stop("The column name with statistics values in the stat2 object must", 
        " be specified in the statCol2 argument.")      
    }

   # checkStructure function searches Chrom, Start and End column from data, and 
   # it removes NA and random scaffolds chromosomes.   
    gaps            <- checkStructure(gaps, rmrnd=!plotRndchr,
                                        columnNames=segmentCols)
    bands            <- checkStructure(bands, rmrnd=!plotRndchr,
                                        columnNames=segmentCols,
                                        optional.columns=c(CategoryName, Colors)) 
    annot1            <- checkStructure(annot1, rmrnd=!plotRndchr,
                                        columnNames=segmentCols)
    annot2            <- checkStructure(annot2, rmrnd=!plotRndchr, 
                                        columnNames=segmentCols)
    annot3            <- checkStructure(annot3, rmrnd=!plotRndchr,
                                        columnNames=segmentCols)
    annot4            <- checkStructure(annot4, rmrnd=!plotRndchr, 
                                        columnNames=segmentCols)
    stat            <- checkStructure(stat, rmrnd=!plotRndchr, 
                                        columnNames=statCols,
                                        optional.columns=c(ID, Colors))
    stat2            <- checkStructure(stat2, rmrnd=!plotRndchr, 
                                        columnNames=statCols,
                                        optional.columns=c(ID, Colors))
    segment        <- checkStructure(segment, rmrnd=!plotRndchr, 
                                        columnNames=segmentCols,
                                        optional.columns=c(CategoryName,
                                        Category2Name, Colors))
    segment2        <- checkStructure(segment2, rmrnd=!plotRndchr,
                                        columnNames=segmentCols,
                                        optional.columns=c(CategoryName, 
                                        Category2Name, Colors))
                                           
    if(!is.null(gaps)){
        centroms  <- subset(gaps, Name%in%"centromere")
    }
           
    if(!is.null(org)){
        if(!is.null(annot1)){
            stop("Choose org or annot arguments to plot annotation data")
        } else{
            print(org)
            annot1 <- getAnnotBiomaRt(org, annotBiomaRt=TRUE)
        }
    } 
        
    if(!is.null(bands)) {  
        chrbands <- build.track(bands, bycol=Chrom, poscol=Position,
                                start.col=Start, end.col=End,
                                category.col=CategoryName, no.histogram=TRUE)    
    } 
    
    # Set chromosome lengths
    if(!is.null(gaps)) {
        chr.length <- as.list(tapply(gaps[ ,End], gaps[ ,Chrom], max)) 
    } else if(!is.null(bands)) {
        chr.length <- as.list(tapply(bands[ ,End], bands[ ,Chrom], max))
    } else if(!is.null(annot1)) {  
        chr.length  <- as.list(tapply(annot1[, End], annot1[ ,Chrom], max))  
    } else {
        stop("I cannot find any track to set the chromosomes lengths.")
    }
    
    # Set subset of chromosomes
    chr.length     <- sort.chrom(chr.length)
    if(!is.null(chr)) {
        if( any ( !chr %in% names(chr.length))) 
            stop("Some chromsomes in 'chr' cannot be drawn.") 
         
        chr.length <- chr.length[names(chr.length) %in% chr]
    }  
    nchrom       <- length(chr.length)

    # Se numer of figure rows
    if(is.null(figCols )) {
        figCols      <- ceiling(nchrom/2)
    }
    maxlen       <- max (unlist(chr.length)) 
    lastChrom    <- names (chr.length)[nchrom]
    smallestChr  <- which.min(unlist(chr.length))
    figRows      <- ceiling (nchrom/figCols)

    # Set chromosome to place legends
    if(is.null(legChrom)) {
        legChrom  <- names (chr.length)[max(smallestChr-1L, 1L)]
    } else if(!is.na(legChrom)){
        if(!legChrom %in% names (chr.length))
           stop("The chromosome provided in legChrom is not present:", legChrom) 
    }
    leg2Chrom <- names(chr.length)[smallestChr]
    legChromBand <-  names(chr.length)[max(smallestChr-2L, 1L)]

    # Create data tracks     
    annot1Hist  <- build.track(annot1, bycol=Chrom, poscol=Position,
                                start.col=Start, end.col=End, binsize=bin,
                                no.histogram = FALSE, max.segments=maxSegs,
                                plot.range=plot_range, chroms=chr)
    annot1_plot_range  <- attr(annot1Hist, "count_range")

    
    annot2Hist <- build.track(annot2, bycol=Chrom, poscol=Position, 
                                start.col=Start,
                                end.col=End, binsize=bin, max.segments=maxSegs,
                                no.histogram=FALSE, plot.range=plot_range,
                                orig.range=annot1_plot_range, chroms=chr)

    annot3Hist <- build.track(annot3, bycol=Chrom, poscol=Position, 
                                start.col=Start,
                                end.col=End, binsize=bin, max.segments=maxSegs,
                                no.histogram=FALSE, plot.range=plot_range,
                                orig.range=annot1_plot_range, chroms=chr)
               
    annot4Hist <- build.track(annot4, bycol=Chrom, poscol=Position, 
                                start.col=Start,
                                end.col=End, binsize=bin, max.segments=maxSegs,
                                no.histogram=FALSE, plot.range=plot_range,
                                orig.range=annot1_plot_range, chroms=chr)
                                      

    chrsegs <- build.track(segment, bycol=Chrom, poscol=Position, 
                            start.col=Start, end.col=End, binsize=bin, 
                            max.segments=maxSegs, category.col=CategoryName,
                            category.col2=Category2Name, no.histogram=noHist,
                            plot.range=plot_range, chroms=chr) 
      
    chrsegs2 <- build.track(segment2, bycol=Chrom, poscol=Position, 
                            start.col=Start, end.col=End, binsize=bin, 
                            max.segments=maxSegs, category.col=CategoryName,
                            category.col2=Category2Name, no.histogram=noHist,
                            plot.range=plot_range, chroms=chr)       
    
    if(!is.null (stat)) {
        stat         <- stat[stat$Chrom %in% names(chr.length),]
        trackStat    <- build.stat(stat=stat, statcol=statCol, bycol=Chrom, 
                                    poscol=Position, start.col=Start,
                                    end.col=End,id.col=ID, 
                                    maxgenecount=maxGeneCount, stat.typ=statTyp, 
                                    scex=scex, spty=spty, binsize=bin, 
                                    margin=margin, statthreshold=statThreshold,
                                    col.stat=colStat, statSumm=statSumm, 
                                    colors.col=Colors, chroms=chr) 
        maxstat       <- attr(trackStat, "maxstat")
        realMinStat   <- attr(trackStat, "Realminstat")
        realMaxStat   <- attr(trackStat, "Realmaxstat")
        minstat       <- attr(trackStat, "minstat")
        colScaleStat  <- attr(trackStat, "colScaleStat")
    }

    if(!is.null (stat2)) {
        stat2         <- stat2[stat2$Chrom %in% names (chr.length),]
        trackStat2    <- build.stat(stat2, statcol=statCol2, bycol=Chrom, 
                                    poscol=Position, start.col=Start, 
                                    end.col=End, id.col=ID,
                                    maxgenecount=maxGeneCount, stat.typ=statTyp, 
                                    scex=scex, spty=spty, binsize=bin,
                                    margin=margin, statthreshold=statThreshold,
                                    col.stat=colStat2, statSumm= statSumm, 
                                    colors.col=Colors, chroms=chr)
        maxstat2       <- attr(trackStat2, "maxstat")
        realMinStat2   <- attr(trackStat2, "Realminstat")
        realMaxStat2   <- attr(trackStat2, "Realmaxstat")
        minstat2       <- attr(trackStat2, "minstat")
        colScaleStat2  <- attr(trackStat2, "colScaleStat")

    }
    
    # Start plotting
    par(xpd=FALSE, mar=c(0, 0, 0, 0), cex=imgscl, mfrow=c(figRows, figCols),
    oma=c(4, 6*imgscl*.8, 4, 2), xpd=NA)     
    # add 'upper' margin on top and bottom
    plot.window (ylim=c( maxlen*(1+upper), -maxlen * upper),
    xlim=c(-maxGeneCount, maxGeneCount))

    xylims   <- par("usr")
    chromwd  <- usr2pt ( margin, axisunit="x", lpi=96 )  
    segLwd   <- segLwd*((11+5)/(figCols+5))
    trckMar  <- chrplotMar ( margin, chrSide, tracks=list(annot1Hist, 
    annot2Hist, annot3Hist, annot4Hist, chrsegs, chrsegs2, trackStat, 
    trackStat2), segLwd, lpi=96, stacked=stack)
    
    tickLoc     <- floor (seq(maxGeneCount, margin, length=3)) 
    tickLab     <- -floor(tickLoc-margin) 
    i           <- 0
    gnotPlotted <- TRUE
    lnotPlotted <- TRUE
    lnotPlotted <- TRUE
    snotPlotted <- TRUE
    lanotPlotted<- TRUE

    for(chrom in names(chr.length)){
        cat ("Chrom", chrom, ":", chr.length[[chrom]], "bp\n")
        plot (c(pos, pos), c(pos, chr.length[[chrom]]), lwd=chromwd, type="l",
        axes=FALSE, xlab="", ylab="", ylim=c( maxlen*(1+upper),-maxlen * upper), 
        xlim=c( -maxGeneCount, maxGeneCount), col=colChrom)
      
        # Print Chrom name
        text(0, -maxlen*.1, chrom, cex=cex*1.5)
        ## plot bands
        plot.track(chrbands, chrom, pos, side=1, lwd=chromwd, sort_segs=FALSE,
                    col_segs=aligncolors, ylims=c(0, chr.length[[chrom]]),
                    overlap.ok=TRUE, stack=FALSE, legchr=legChromBand, cex=cex, 
                    title="Chromosome banding")

        ## plot annot1
        plot.track(annot1Hist, chrom, margin, side=chrSide[1], 
                    lwd=segLwd, sort_segs=sortSegs, col_segs=colAnnot1, 
                    ylims=c(0, chr.length[[chrom]]), bin=bin,
                    trckmar=trckMar[1])
          
        ## plot annot2
        plot.track (annot2Hist, chrom, trckMar[2], side=chrSide[2],
                    lwd=segLwd, sort_segs=sortSegs, col_segs=colAnnot2,
                    ylims=c(0, chr.length[[chrom]]), bin=bin)
        
        ## plot annot3
        plot.track(annot3Hist, chrom, trckMar[3], side=chrSide[3], 
                   lwd=segLwd, sort_segs=sortSegs, col_segs=colAnnot3, bin=bin)

        ## plot annot4
        plot.track(annot4Hist, chrom, trckMar[4], side=chrSide[4],
                    lwd=segLwd, sort_segs=sortSegs, col_segs=colAnnot4, 
                    ylims=c(0, chr.length[[chrom]]), bin=bin)
          
        ## plot stats
        plot.track(trackStat, chrom, trckMar[7], side=chrSide[7],
                    lwd=segLwd, sort_segs=sortSegs, col_segs=colStat,
                    cex=imgscl)

        ##plot stat2
        plot.track(trackStat2, chrom, trckMar[8], side=chrSide[8],
                    lwd=segLwd,  sort_segs=sortSegs, col_segs=colStat2,
                    cex=imgscl)

        ## segments
        plot.track (chrsegs, chrom, trckMar[5], side=chrSide[5], 
                    lwd=segLwd, sort_segs=sortSegs, col_segs=colSegments, 
                    ylims=c(0, chr.length[[chrom]]), stack=stack, cex=cex,
                    legchr=legChrom, title=segmentDesc)

        ## segments2
        plot.track(chrsegs2, chrom, trckMar[6], side=chrSide[6],
                   lwd=segLwd, sort_segs=sortSegs, col_segs=colSegments2,
                   ylims=c(0, chr.length[[chrom]]), stack=stack,
                   legchr=leg2Chrom, cex=cex, title=segment2Desc)
        
        #plot centroms
        if(!is.null(centroms)){ 
            if(chrom %in% centroms$Chrom) {
                points(pos, centroms[centroms[,Chrom] %in% chrom, Start],
                pch=19,bg="black", cex=chromwd*.3 )              
            } else {
                warning("Chromosome ", chrom, " missing from gaps file.")
            }
        }
  
        ## Y-axis
        if((i%%figCols) == 0 & yAxis) {
            tickLoc <- axTicks (2)  
            tickLab <- sprintf("%.0f", tickLoc/1e6)
            axis(2, at=tickLoc, labels=tickLab, las=2, cex.axis=cex)
            mtext("Mb", 2, line=2, las=2, cex=cex)
        }
        #plot annot count scale 
        if(!is.null(annot1)) {
            if(gnotPlotted & (i == smallestChr-2 )) {
                draw.scale(y=0.80*chr.length[[chrom]] + 0.20 * xylims[3],
                minval=margin, maxval=maxGeneCount, minlab=annot1_plot_range[1], 
                maxlab=annot1_plot_range[2], lwd=4, col=colAnnot1, cex=cex, 
                title=scale.title)
                
                gnotPlotted <- FALSE
            }
        }
        ## Stat Scale
        if(snotPlotted & (i == smallestChr-1) & !is.null(stat))
            if(!is.null(attr(trackStat, "is_int")) & attr(trackStat, "sdstat")>0) {    
                draw.scale(y= 0.80 * chr.length[[chrom]] + 0.20 * xylims[3],
                lwd=4, col=colScaleStat, minval=margin, maxval=maxGeneCount, 
                minlab=realMinStat, maxlab=realMaxStat, title=statName, cex=cex,
                as.is=TRUE)
              
                snotPlotted <- FALSE
        }
        ## Stat Scale2
        if(!is.null(stat2) ) snotPlotted <- TRUE
        if(snotPlotted & (chrom == lastChrom) & !is.null(stat2))
            if(attr(trackStat2, "sdstat")>0) {
                draw.scale(y= 0.65 * chr.length[[chrom]] + 0.35 * xylims[3],
                lwd=4, col=colScaleStat, minval=margin, maxval=maxGeneCount, 
                minlab=realMinStat2, maxlab=realMaxStat2, title=statName2, 
                cex=cex, as.is=TRUE) 
                
                snotPlotted <- FALSE
        }
        if(!is.null(segment) & (i == smallestChr - 1 | chrom == lastChrom)) {
            # Histogram, plot scale
            if(!attr(chrsegs, "is_int")){
                draw.scale(y = 0.80 * chr.length[[chrom]] + 0.20 * xylims[3],
                minval=margin, maxval=maxGeneCount,
                minlab=attr(chrsegs, "count_range")[1],
                maxlab=attr(chrsegs, "count_range")[2],
                lwd=segLwd, col=colSegments, title=segmentDesc,  cex=cex)
            }  
        }  
        if(!is.null(segment2) & (i == smallestChr - 1 | chrom == lastChrom)) {
            # Histogram, plot scale
            if(!attr(chrsegs2, "is_int")){
                draw.scale(y = 0.60 * chr.length[[chrom]] + 0.40 * xylims[3],
                minval=margin, maxval=maxGeneCount,
                minlab=attr(chrsegs2, "count_range")[1],
                maxlab=attr(chrsegs2, "count_range")[2],
                lwd=segLwd, col=colSegments2, title=segment2Desc,  cex=cex)
            }    
        }  

    i = i + 1
    } #End loop
    
    # Plot a title
    if(!is.null(title))
        mtext(title, cex=cex*.9, outer=TRUE, line=1)

} ## chromPlot

