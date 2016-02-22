## Track management functions
## ----- ---------- ---------

build.track <- function(datatab, bycol="Chrom", poscol=c("Mid", "Start"),
start.col="Start", end.col="End", category.col=NULL, binsize=1e6,
no.histogram=FALSE, max.segments=200, check.unit=FALSE, 
category.col2=NULL, plot.range=NULL, orig.range=NULL, chroms=NULL) {
    if(is.null(datatab))
        return(NULL)
    if(nrow(datatab)==0L)
        return(NULL)
        
    # Filter chromosomes
    if(!is.null(chroms)) {
        if(any(!chroms %in% datatab[,bycol]))
            warning("Chromosomes not in the dataset: ", paste(setdiff(chroms, 
                    datatab[,bycol]), collapse=", "))
        datatab <- datatab[datatab[,bycol] %in% chroms,]
    }
        
    groupby <- NULL
    Group  <- match(category.col, colnames(datatab))
    Group  <- Group[!is.na(Group)][1]#

    # Find position column
    poscoli <- match(poscol, colnames(datatab))
    if(sum(is.na(poscoli))==length(poscol)) {
        stop("'", paste(poscol[is.na(poscoli)], collapse=", "),
        "' is not a column of datatab.")
    } else {
       poscol  <- poscol[!is.na(poscoli)][1L]
    }
    
    # Remove rows with missing values for postions columns
    datatab <- datatab[!is.na(datatab[,bycol]) & !is.na(datatab[,poscol]),]
    if(check.unit)
        datatab   <- checkunit(datatab)
    
    maxseglen <- max(datatab[,end.col] - datatab[,start.col] - 1)

    
    if(!is.null(category.col)) {
        if(category.col %in% colnames(datatab)) {
            col.idx <- match(category.col, colnames(datatab))
            colnames(datatab)[col.idx] <- "Category"
        }
    }

    if(!is.null(category.col2)) {    
        if(category.col2 %in% colnames(datatab)) {    
            col.idx2 <- match(category.col2, colnames(datatab))
            colnames(datatab)[col.idx2] <- "Category2"
        }
    }

    if("gieStain" %in% names(datatab)) {
        datatab     <- cytoBands(datatab)
    } 
    
    if((maxseglen < binsize | nrow(datatab) >= max.segments) & !no.histogram){
    
        # build track histogram
        if("Category" %in% colnames(datatab)) {
            groupby=split(datatab$Category, datatab[,bycol])
        }
        track_ls <- list()
        for(chr in sort(unique(datatab[,bycol]))) {
            rows <- datatab[,bycol] %in% chr
            track_ls[[chr]] <- split.bins(datatab[rows,poscol], step=binsize, 
                                            fun=length, group.by=groupby[[chr]])     
        }
        track_ls <- track_ls[!sapply(track_ls, is.null)]
        track_ls <- sort.chrom(track_ls)
        # Adjust values to preset scale if provided
        if(!is.null(plot.range)){
            if(length(plot.range)!=2L) {
                stop("plot.range must ve a vector with 2 values, min and max.")
            }
            if(is.null(orig.range)) {
                if(is.null(dim(track_ls[[1]]$counts))){
                    cmin <- min(sapply(track_ls, function(x) min(x$counts)))
                    cmax <- max(sapply(track_ls, function(x) max(x$counts)))
                } else {
                    cmin <- min(sapply(track_ls,
                                function(x) min(rowSums(x$counts))))
                    cmax <- max(sapply(track_ls, 
                                function(x) max(rowSums(x$counts)))) 
                }
            } else {
                cmin <- orig.range[1]
                cmax <- orig.range[2]
            }
            # forces min count to 0 for histograms
            track_ls <- lapply(track_ls, 
                            function(x){x$counts=scale.vec(x$counts, 0, cmax,
                            plot.range[1], plot.range[2]);x}) 
            attr(track_ls, "count_range") <- c(cmin, cmax)
        }
        attr(track_ls, "is_int") <- FALSE
    }else{
        track_ls <- split(datatab, datatab[,bycol])
        track_ls <- sort.chrom(track_ls)
        attr(track_ls, "is_int") <- TRUE        
    }
         
    return(track_ls)
} 

build.stat <- function(stat=NULL, statcol=NULL, maxgenecount=5, stat.typ="p", 
scex=1, spty=20, binsize=1e6, margin=margin, statthreshold=NULL,
col.stat="blue", col.stat2="gray57",  statSumm="none", bycol="Chrom",
poscol=c("Mid", "Start"), start.col="Start", end.col="End", id.col="ID",
colors.col="Colors", chroms=NULL) {
    # Sanity checks
    if(!statcol %in% colnames(stat) & !id.col %in% colnames(stat)) {
        stop(paste(statcol, " not a column of stat object."))
    } 
    if(id.col %in% colnames(stat) & statSumm != "none") {
        warning("the statSumm argument should be 'none' if an 'ID'",
        " columns is present. statSum will be ignored for this track.")  
    }
    # Filter chromosomes
    if(!is.null(chroms)) {
        if(any(!chroms %in% stat[,bycol]))
            warning("Chromosomes not in the dataset: ", paste(setdiff(chroms, 
                    stat[,bycol]), collapse=", "))
        stat <- stat[stat[,bycol] %in% chroms,]
    }
    # Find position column
    poscoli <- match(poscol, colnames(stat))
    if(sum(is.na(poscoli))==length(poscol)) {
        stop("'", paste(poscol[is.na(poscoli)], collapse=", "),
             "' is not a column of stat.")
    } else {
        poscol  <- poscol[!is.na(poscoli)][1L]
    }
        
    # Declare variables
    maxstat     <- NULL
    minstat     <- NULL
    Realmaxstat <-NULL
    # Always order stat by position
    stat      <- stat[order(stat[,bycol], stat[,poscol], decreasing=FALSE ) , ] 
        
    if(id.col %in% colnames(stat)) { 
        stat <- stat[!is.na(stat[,poscol]), ]
        stat <- stat[!is.na(stat[,bycol]), ]
        stat[, id.col]<-as.character(stat[,id.col])
        NaID   <- is.na(stat[,id.col])
        NULLID <- is.null(stat[,id.col]) 
        stat[NaID, id.col] <- "" 
        
        if(!statcol %in% colnames(stat)) {
            stat <- cbind(stat, margin)
            colnames(stat)[ncol(stat)] <- statcol
        }
        
        if( colors.col %in% colnames(stat)) {
            stat <-stat[,c(bycol, poscol, end.col, id.col, statcol, colors.col)]
        } else {
            stat[,colors.col] <-"black"
            stat <- stat[,c(bycol, poscol,end.col, id.col,statcol, colors.col)]
        }   
    } 

    if(statSumm== "none" | (id.col %in% colnames(stat))){
        stat_vals<-stat
        colnames(stat_vals)[match(c(poscol,statcol),
                                  colnames(stat_vals))]<-c("mids", "statSumm")
        trackStat     <- split (stat_vals, as.factor(stat_vals[,bycol]))
        maxstat        <- max(stat[,statcol])
        minstat        <- min(stat[,statcol])
        Realminstat    <- minstat           
        Realmaxstat    <- maxstat           
    } else{
        stat_vals  <- data.frame(as.numeric(stat[,poscol]),
                                as.numeric(stat[,statcol]))
        chromstat  <- split(stat_vals, as.factor(stat[,bycol]))   
        trackStat  <- lapply(chromstat, split.bins, step=binsize, fun=statSumm)
        trackStat  <- trackStat[!sapply(trackStat, is.null)]
        maxstat   <- max(unlist(sapply(trackStat, function(x) max(x$statSumm))))
        minstat   <- min(unlist(sapply(trackStat, function(x) min(x$statSumm))))
        Realminstat <- minstat         
        Realmaxstat <- maxstat
    } #else

    if(!is.null(statthreshold)){
        trackStat    <- lapply(trackStat, function(x){
                        x$Colors=ifelse(x$statSumm >= statthreshold,
                                        col.stat, col.stat2);x} )
    } else {
        trackStat    <- lapply(trackStat, function(x){
                        x$Colors=col.stat; x} )
    }    
    sdstat         <- sd(unlist(sapply(trackStat, function(x) x$statSumm)))

    if(sdstat ==0) { # Value was a constant
        trackStat <- lapply(trackStat, function(x){x$statSumm=margin;x})
    } else { # Value varied, must scale
#Use this code to increase observable differences by removing the minimum value.  
         trackStat  <- lapply(trackStat, 
                          function(x){x$statSumm=(x$statSumm-minstat)/
                                      (maxstat-minstat)*
                                      (maxgenecount-margin)+margin;x})
    }
                                
    maxstat  <- max(unlist(sapply(trackStat, function(x) max(x$statSumm))))
    minstat  <- min(unlist(sapply(trackStat, function(x) min(x$statSumm))))

    attr(trackStat, "colScaleStat")  <- col.stat
    attr(trackStat, "is_int" )       <- FALSE
    attr(trackStat, "stats.typ" )    <- stat.typ
    attr(trackStat, "scex" )         <- scex #
    attr(trackStat, "spty" )         <- spty #   
    attr(trackStat, "maxstat" )      <- maxstat #
    attr(trackStat, "Realminstat" )  <- Realminstat #
    attr(trackStat, "Realmaxstat" )  <- Realmaxstat #
    attr(trackStat, "minstat" )      <- minstat
    attr(trackStat, "sdstat")        <- sdstat
     
    return(trackStat)
}



## Plotting functions
## -------- ---------

plot.track <- function(track_ls, which.track, margin, side=1, lwd=3,
     sort_segs=TRUE, col_segs=NULL, ylims=NULL, stack=FALSE, form=NULL, bin=1e6,
     cex=NULL, trckmar=0, legchr=NULL, title=NULL, overlap.ok=FALSE) {
    
    if(is.null(track_ls))
        return(NULL)
    
    track   <- track_ls[[which.track]]
    
    if("ID" %in% colnames(track_ls[[which.track]]) ){ 
        if(any(track$ID != "")) {
            linelength <- strwidth("xx", units="user", cex=cex)
            labtrack <- track[track$ID != "",]
            unstacked <- unstacked(x= 1, y=labtrack$mids, trackId= labtrack, 
            cex=cex)
            if("Colors" %in% colnames(track)){
                 lab_colors <- as.character(track$Colors)
            } else {
                lab_colors <- col_segs
            }
            segments (side*labtrack$statSumm, labtrack$mids,
                        side*(labtrack$statSumm+linelength),
                        unstacked$y, col="black", lwd=0.7, lend=1)
            text(side*(labtrack$statSumm+linelength), unstacked$y, 
                unstacked$ID, col=unstacked$Colors, adj=0.5-(side*0.5),
                lwd=lwd, cex=cex)
        }
    }

    typ     <- attr(track_ls, "stats.typ")
    scex    <- attr(track_ls, "scex")
    spty    <- attr(track_ls, "spty")
    flag    <- 0 

    if(attr(track_ls, "is_int")) {
        if("Colors" %in% colnames(track_ls[[1]])){
             col_segs   <- sort(unique(unlist(sapply(track_ls, "[", "Colors"))))
        }
        if("Category" %in% colnames(track_ls[[1]])) {
            all_categories1 <- unlist(sapply(track_ls, "[", "Category"))
            slevs1 <- sort(unique(all_categories1))
            nlevs1 <- length(slevs1)
            if(length(col_segs)<nlevs1) {
                 stop("Not enough colors provided. There are ", nlevs1, 
                      " categories of elements and ", length(col_segs),
                      " colors.")
             }
        }
        if("Category2" %in% colnames(track_ls[[1]])){
            seg_form        <- c(15, 16, 17, 18, 25, 20:30)
            all_categories2 <- unlist(sapply(track_ls, "[", "Category2"))
            slevs2 <- sort(unique(all_categories2))
            nlevs2 <- length(slevs2)
            if(length(seg_form)<nlevs2) {
                stop("Not enough forms provided. There are ", nlevs2,
                    " categories of elements and ", length(seg_form),
                    " forms.")
            }
        }
    }
    
    if(!is.null(track) ) {
        if(attr(track_ls, "is_int")) {
            if("Colors" %in% colnames(track)){
                seg_colors <- as.character(track$Colors)            }
            if("Category" %in% colnames(track)) {
                seg_colors <- col_segs[match(track$Category, slevs1)]
            } else if(!"Colors" %in% colnames(track)) {
                seg_colors <- rep(col_segs[1L], nrow(track))     
            }
            if("Category2" %in% colnames(track)){
                form <- seg_form[match(track$Category2, slevs2)]
            }
        }
        if(!is.null(typ)) {   
            if("Colors" %in% names(track)) {
                seg_colors <- track$Colors
            } else {
                seg_colors <- rep(col_segs[1L], length(track$mids))     
            }
            flag=1
            switch(typ,
            l=lines(side * (track$statSumm), track$mids,
                    col=seg_colors, type="l"),
            p=lines(side * (track$statSumm), track$mids, 
                    col=seg_colors, type="p",
            n=NULL,
            pch=spty, cex=scex),
            hst=attr(track_ls, "is_int")<-FALSE,
            seg=attr(track_ls, "is_int")<-FALSE)     
        }

        if(attr(track_ls, "is_int")) {
            if(overlap.ok) {
                segments(margin*side, track[,"Start"],margin*side,track[,"End"],
                         col=seg_colors, lwd=lwd, lend=1)
            } else {
                if(!stack){
                    plot.lodclust(as.matrix(track[,c("Start", "End")]),
                                  col=seg_colors, lwd=lwd, ymax=margin*side,
                                  stack_ints=FALSE, form=form)                  
                } else {
                    plot.segments(track, seg_colors, startfrom=side*margin, 
                                  lwd=lwd, sort.segs=sort_segs, form=form)
                }
            }
        } else{
            cat("Small segments, I will plot histograms for this track:\n")
            if(flag==0){
                if(trckmar>margin) {
                    track$counts <- scale.vec(track$counts, margin, 1,trckmar,1)
                    margin <- trckmar
                }
                if(!is.null(dim(track$counts))) {
                     plot.stackedHist(track, margin, ylims=NULL, col=col_segs, 
                                      stepsize=bin, side=side)
                } else {
                     chromhist(side*margin, side*(track$counts),
                     track$mids, ylims=NULL, col=col_segs, stepsize=bin)#bin
                }
            }   
        }
    }    

    # Plot a legend for is_int with Category
    if(!is.null(legchr)) if(!is.na(legchr)) {
        if("Category"%in%colnames(track_ls[[1]])&attr(track_ls, "is_int") & 
            !is.null(legchr)){
            if(which.track == legchr) {
                lcols  <- col_segs[1:nlevs1]
                legend("bottom",legend=slevs1,ncol=max(1,floor(sqrt(nlevs1))-1), 
                       lwd=3, col=lcols, lty=1, xpd=NA,
                cex=cex, bg="white", title=title)
            }
        }

        if("Category2"%in%colnames(track_ls[[1]])&attr(track_ls, "is_int") &
            !is.null(legchr)){
            if(which.track == legchr) {
                lform  <- seg_form[1:nlevs2]
                legend("right",
                legend=slevs2, ncol=max(1,floor(sqrt(nlevs2))-1),
                lwd=3, col="black", lty=1, xpd=NA, pch=lform,
                cex=cex, bg="white", title=paste(title, "2nd Category"))
            }
        }
    }
} # plot.track


plot.stackedHist <- function(track=NULL, margin=NULL, ylims=NULL, col=NULL, 
stepsize=1e6, side=1) {
    data <- t(track[[2]]) # assuming that first elements is always mids

    for(i in 1:nrow(data)){
        if(i==1){
            x1 <- margin
            x2 <- as.numeric(data[i,])    
        } else{
            x1 <- x2
            x2 <- x1 + as.numeric(data[i,]) - margin
        }
        chromhist(side*x1, side*x2, track$mids, col=col[i], stepsize=stepsize,
                ylims=ylims)
    }
} # plot.stackedHist

plot.int <- function(intervals, intclusters, col, lwd=3, ymax, stack_ints=TRUE, 
form) {
    intervals <- as.matrix(intervals)
    if(ncol(intervals)==1)
        intervals <- t(intervals)
    for(i in 1:length(intclusters)){
        incluster=intclusters[[i]]
        plot.lodclust(intervals[incluster,], col=col[incluster], lwd, ymax,
        stack_ints=stack_ints, form=form[incluster])
    }
} 
chromhist <- function(xleft, xright, y, ylims=NULL, col="red", stepsize=1) {
    ytop <- y - stepsize/2
    ybot <- y + stepsize/2 
    if(!is.null(ylims)) {
        ytop <- unlist(sapply(ytop, max, ylims[1]))
        ybot <- unlist(sapply(ybot, min, ylims[2]))
    }
    rect(xleft, ybot, xright, ytop, col=col, border=NA)
}


draw.scale <-function(y, minval, maxval, lwd, col, minlab=minlab, maxlab=maxval,
title="Counts", cex, ...) { 
    lht      <- sign(diff(par("usr")[3:4])) * strheight(title)

    tickpos  <- pretty_ticks(minval, maxval, 2, ...)
    ticklab  <- pretty_ticks(minlab, maxlab, 2, ...)
    
    ticklab <- format.bignum(ticklab, unit="", sep="")
    tickpos <- tickpos - tickpos[length(tickpos)] *0.5 # move it 50% to the left
    
    lines(c(tickpos[1], tickpos[length(tickpos)]), rep(y+lht*2,2), col=col,
    lwd=lwd, lend=1)
    
    segments(tickpos, y+lht*1, tickpos, y+lht*3, col="black", lwd=1)
    text(tickpos[1] + (tickpos[length(tickpos)]-tickpos[1])*.7, y+lht*0,
    title, adj=.5, font=1, cex=cex) 
    text(tickpos, y+lht*4, ticklab, cex=cex) # Bottom side of scale
} 

plot.segments <- function(segobj, colvec, startfrom, lwd=3, sort.segs=TRUE,
form=NULL) {
    intsmat <- segobj[, c("Start", "End")]
    intslen <- intsmat[,2]-intsmat[,1]
    clusts  <- find.intclusters(intsmat, sort.by.size=sort.segs)
    plot.int(intsmat, clusts, colvec, lwd=lwd, startfrom, stack_ints=TRUE,
            form=form)
}

plot.lodclust<-function(intmat, col, lwd=3, ymax, stack_ints=TRUE, form=NULL){
    intmat  <- matrix(intmat, ncol=2)
    lwdusr  <- pt2usr(lwd,  axisunit="x", lpi=96)
    onestep <- lwdusr * 1.2 * sign(ymax)
    y       <- ymax
    idx     <- 1L:nrow(intmat)
    assigned_levs <- idx
    for(i in 1L:nrow(intmat)) {
        if(i > 1L) {
            if(stack_ints) {
                used_lev<-assigned_levs[sapply(idx,function(x) overlap(intmat[i,],
                intmat[x,]) & x!=i)]
                clear <- min(setdiff(idx, used_lev))
                assigned_levs[i] <- clear
            } else {
                clear <- i
            }   
            y <- ymax + onestep * (clear-1L)
        }
    
        plot.lodint(intmat[i,], col=col[i], lwd=lwd, y=y,form=form[i])
    }
}

plot.lodint<-function(lodint, lty=2, col=1, lwd=3, y=NULL, form=1) {
    lowint<-numeric()
    highlim<-numeric()
    if(length(dim(lodint))>1) {
        lowlim<-lodint[1,2]
        highlim<-lodint[nrow(lodint),2]
    }else{
        lowlim<-lodint[1]
        highlim<-lodint[2]
    }
    if(!is.null(form)) {
       line_col <- "darkgray"
    } else {
       line_col <- col
    }
     
    lines(rep(y, 2), c(lowlim, highlim), col=line_col, lwd=lwd, lend=1)
    if(!is.null(form)){
        points(y, mean(c(lowlim, highlim)), col=col, lwd=lwd, lend=1,
        pch=form, cex=1)
    }
     
}

## Auxiliary functions
## --------- ---------

split.bins <-function(dataframe, column=NULL, step, fun=length, group.by=NULL){
## splits data by bins of 'step' size and apply 'fun' to each cell resulting 
## from the combination of levels bin and the group.by column in dataframe
## (if provided)
##
## dataframe    object to split n bins
## column       column to use as variable to create bins. Can 
##                be a integer or a string with colname
## step         bin size. values > 1 are used as absolute
##                size or a fraction of the range in column

    if(is.null(fun))
        return(dataframe)
        
#     if(is.character(fun))
#         fun <- get(fun)
        
    if(is.null(column))
        x<-dataframe
    else
        x<-dataframe[,column] 
        
    if(is.data.frame(x) | is.matrix(x)) {
        if(ncol(x)>1) {
            y=x[,2]
            x=x[,1]
        } else {
            y=x
        }
    }else{
        y=x
    } 
    y<-y[!is.na(y)]
    x<-x[!is.na(x)] 
    if(identical(length(x), 0L))
        return(NULL) 
    span<-diff(range(x, na.rm=TRUE))   
    if(step<=1)
        step<-span*step  
    bin<-ceiling(x/step)  
    mids<-round((bin*step+(bin-1)*step)/2)

    if(is.null(group.by)) {
        byfactor <- list(mids)
    } else {
        byfactor <- list(mids, group.by)
    }

    out<-tapply(y, byfactor, fun)
    if(!is.null(group.by) & identical(fun, length)) {
        out[is.na(out)] <- 0
    }


    outname<-deparse(substitute(fun))
    if(grepl('length', outname)) outname<-'counts'
    if(is.null(group.by))
        out <- as.numeric(out)
        
    output<-list(mids=sort(unique(mids)), out)
    names(output)[length(output)]=outname
    
    return(output)
} # split.bins

chrplotMar <- function(mar, sides, tracks, lwd, lpi=96, stacked=FALSE) {
## Calculate appropriate margin for each track according to what other tracks
## are plotted on the same chromosomal side. If a histogram and segments are
## plotted on the same side, the histogram start where segments end.
    if(length(sides)!=length(tracks))
        stop("Arguments sides and tracks must be of same length.")
    if(any(!sides %in% c(-1,1)))
        stop("values of sides must be in c(-1, 1)")
        
    out      <- rep(mar, length(sides))
    os       <- pt2usr(lwd,  axisunit="x", lpi=lpi)
    are_drawn<- !sapply(tracks, is.null)
    are_ints <- sapply(tracks, attr, "is_int")
    are_ints <- unlist(sapply(are_ints, function(x) if(is.null(x)) FALSE else x))
    off_set  <- which(are_drawn & !are_ints)
    side_height <- list(`-1`=NULL, `1`=NULL)
    

    ## If plotting int segments and stacking, calculate max cluster size to
    ## adjust margin
    track_height <- rep(1, length(tracks))
    if(sum(are_ints)>0 & stacked) {
      track_height[are_ints] <- sapply(tracks[are_ints], max.cluster.height)
    }
    # Max height per chrom side
    side_height_tmp <- tapply(track_height, sides, max) 
    side_height[names(side_height_tmp)] <- side_height_tmp
    
    # Add the hight of one line each segment in the same side as a hist
    if(length(off_set)>0L) {
        for (i in 1:length(off_set)) {
            if(any(sides[-i] %in% sides[i] & are_ints[-i]))
                out[i] <- mar + os * side_height[[sides[i]*.5+1.5]]
        }
    }
    return(out)
}

format.bignum <- function(number, unit="b", sep=" ") {
    if(length(number)>1)
        return(unlist(sapply(number, format.bignum, unit, sep)))
    
    gigas <- number/1e9
    megas <- number/1e6
    kilos <- number/1e3
    gig_u <- paste(sep="", "G", unit)  
    meg_u <- paste(sep="", "M", unit)
    kil_u <- paste(sep="", "K", unit) 
    if(gigas[length(gigas)] == 0)
        return(0)
    if(gigas[length(gigas)] == as.integer(gigas[length(gigas)]))
        return(paste(gigas, gig_u, sep=sep))
    else if(megas[length(megas)] == as.integer(megas[length(megas)]))
        return(paste(megas, meg_u, sep=sep))
    else if(kilos[length(kilos)] == as.integer(kilos[length(kilos)]))
        return(paste(kilos, kil_u, sep=sep))
    else
        return(paste(format(number, big.mark=",", trim=TRUE), unit, sep=sep))
}

max.cluster.size <- function(segslist) {
    ints_by_chrom   <- lapply(segslist, function(x) x[,c("Start", "End")])
    clusts_by_chrom <- lapply(ints_by_chrom, find.intclusters,
                                sort.by.size=FALSE)
    output <- max(sapply(clusts_by_chrom[[1]], length))
    return(output)
}

max.cluster.height <- function(segslist) {
    ints_by_chrom   <- lapply(segslist, function(x) x[,c("Start", "End")])
    clusts_by_chrom <- lapply(ints_by_chrom, find.intclusters, 
                                sort.by.size=FALSE)

    ints_by_clust   <- list()
    for (chr in 1:length(ints_by_chrom)) {
        ints_by_clust[[chr]] <- lapply(clusts_by_chrom[[chr]], 
                                        function(x) ints_by_chrom[[chr]][x,])
    }

    output <- max(unlist(sapply(ints_by_clust,
                                function(x) sapply(x, maxclustsize))))
    return(output)
}
## this is just an approximation to the number of levels needed for plotting
## a clusters
maxclustsize <- function(ints) {
   tot_size <- sum(apply(ints, 1, diff))
   
    min_possible_size <- max(ints[,2])-min(ints[,1])
    max_possible_size <- mean(ints[,2])-mean(ints[,1])
    estimate1 <- tot_size/min_possible_size
    estimate2 <- tot_size/max_possible_size
    estimated_size <- ceiling((estimate1+estimate2)/2)
    return(estimated_size)
}
scale.track <- function(track, stat.col="statSumm", mintarget, maxtarget) {
## Converts scale of values in a track to a new range of values
## (mintarget,maxtarget)
    maxstat     <- max(unlist(sapply(track, function(x) max(x$statSumm))))
    minstat     <- min(unlist(sapply(track, function(x) min(x$statSumm))))
  
    for(i in length(track)) {
        track[[i]] <- scale.vec(track[[i]][,stat.col], minstat, maxstat,
                                mintarget, maxtarget)
    }

    return(track)
}

scale.vec <- function(vec, minstat, maxstat, mintarget, maxtarget) {
    (vec-minstat)/(maxstat-minstat)*(maxtarget-mintarget)+mintarget
}

checkunit <- function(featable, column=c("Start", "End"), minlen=10e6) {
    if(max(featable[,column], na.rm=TRUE)<minlen) { # it is not Mb
            featable[,column] <- featable[,column] * 1e6
    }
return(featable)
}

sort.chrom <- function(chroms) {
    if(is.list(chroms))
        return(sort.chrom.list(chroms))
    chroms <- as.character(chroms)
    dchroms <-sub("chr", "", chroms)  
    arenum <- grep("^[[:blank:]]*[[:digit:]]+[[:blank:]]*$", dchroms)
    aresex <- setdiff((1:length(chroms)), arenum)
    outchroms <- chroms[arenum][order(as.numeric(dchroms[arenum]))]
    outchroms <- c(outchroms, sort(chroms[aresex]))
    return(outchroms)
}

sort.chrom.list <- function(chromList) {
    chroms <- names(chromList)
    return(chromList[sort.chrom(chroms)])
}

pretty.choose <- function(number, num.options) {
    if(missing(num.options))
        num.options <- c(0, 5, 10, 20, 50, seq(100, 1000, by=100),
        seq(2e3, 1e4, by=1e3), seq(2e4, 1e5, by=1e4), seq(2e5, 1e6, by=1e5),
        seq(1e-01, 1, by=1e-01))
    if(length(number)>1)
        return(unique(unlist(sapply(number, pretty.choose, num.options))))
    output <- number
    num.options <- sort(num.options, decreasing=TRUE)
    for(i in 1:length(num.options)) {
        if(number > num.options[i] | i==length(num.options)) {
            output <- num.options[i]
        break
        } else if(number > num.options[i+1]) {
            if((num.options[i] - number) <= (number - num.options[i+1]))
                output <- num.options[i]
            else
                output <- num.options[i+1]
            break
        } else 
    next
    }
    return(output)
}

pretty_ticks <- function(minval, maxval, howmany=3, as.is=FALSE, digits=2, ...){
    
    outticks <- signif(seq(minval, maxval, length=howmany), digits=digits)
    
    if(!as.is) {
        outticks <- pretty.choose(outticks)
    }

    return(outticks)
}

pt2usr <- function(lwd,  axisunit="x", lpi=96) {
    if(!axisunit %in% c("x", "y"))
        stop("Argument 'axisunit' must be in c('x', 'y').")
    pusr   <- par("usr")
    pinch  <- par("pin")
    lwdin  <- lwd * 1/lpi
    if(axisunit=="x")
        usr.in <- (pusr[2]-pusr[1])/pinch[1]
    else
        usr.in <- (pusr[4]-pusr[3])/pinch[1]
        lwdusr <- lwdin * usr.in
    return(lwdusr)
}

usr2pt <- function(lwdusr, axisunit="x", lpi=96) {
# lwdusr  line width in user coordinates units
# lpi     lines per inch (device dependent)
    if(!axisunit %in% c("x", "y"))
        stop("Argument 'axisunit' must be in c('x', 'y').") 
    pusr   <- par("usr")
    pinch  <- par("pin")  
    if(axisunit == "x")
        usr.in <- (pusr[2]-pusr[1])/pinch[1]
    else
        usr.in <- (pusr[4]-pusr[3])/pinch[1]
    lwdin <- lwdusr / usr.in
    lwdpt <- lwdin * lpi
    return(lwdpt)
}
find.intclusters <- function(intmat, sort.by.size=TRUE){
    indices  <- 1:nrow(intmat)
    clusters <- list()
    i=0
    while(length(indices)>0L) {
        i=i+1
    
        lrgst    <- which.max(intmat[indices,2]-intmat[indices,1])
        members  <- indices[lrgst]
        secindx <- setdiff(indices, members)
        if(length(secindx)>0L) {
            member_mem <- 0
            while(!all(members %in% member_mem)) {
                member_mem <- members
                for(j in secindx) {
                    if(overlap(intmat[members,], intmat[j,])) {
                        members <- unique(append(members, j))
                        secindx <- setdiff(indices, members)            
                    }
                }
            }
        }
        if(sort.by.size)
            members <- members[order(intmat[members,2]-intmat[members,1],
            decreasing=TRUE)]
            clusters[[i]] <- members
            indices=setdiff(indices, members)
    }
    cat("  ", length(clusters), "segment clusters.\n")
    return(clusters)
}
overlap<-function(int1, int2) {
    if(length(int1)==0L | length(int2)==0L)
        return(FALSE)
    if(!is.null(dim(int1)))
        int1 <- c(min(int1[,1]), max(int1[,2]))
    if(!is.null(dim(int2)))
        int2 <- c(min(int2[,1]), max(int2[,2]))    
    int1 <- as.numeric(sort(int1))
    int2 <- as.numeric(sort(int2))
    len1 <- diff(int1) 
    len2 <- diff(int2)
    total_len <- len1 + len2  
    span <- max(int1, int2) - min(int1, int2)
    if(length(total_len)!=0L & length(span)!=0L)
        if(total_len >= span) {
            return(TRUE)
        }
    return(FALSE)
}

checkStructure <- function(data = NULL, rmrnd=FALSE, 
                           columnNames=c("Chrom", "Start", "End"),
                           optional.columns=NULL) {

    if(is.null(data)) return(NULL)

    # GenomicRanges
    if("GRanges" %in% class(data)) {
        data <- data.frame(seqnames(data), start(data), end(data), mcols(data))
        colnames(data)[1:3] <- columnNames
    }

    if(is.character(data)) {
        # I will assume that this is a file name
        if(!file.exists(data) &  !grepl("http", data))
            stop("File ", data, "does not exist.")
            
        Lines<-read.delim(data, sep="\n", header=FALSE, fill=FALSE, 
                            comment.char="#")
        if(grepl("[A|C|G|T]{1,}", Lines[2,]) ){
            cat ("Reading AXT file", "\n")
            data <- transformAlign(Lines)    
        } else {
            cat("Reading file ", data, " assuming bed format...", "\n")
            data <- read.delim(data, sep="\t", header=TRUE, fill=FALSE,
                                comment.char="#")
        }
    } 
        
    name <- deparse(substitute(data))

    if(!is.data.frame(data)){
        # Everything failed....
        stop(name, " could not be read.")
    }

    datacols <- colnames(data)
    
    if(!all(columnNames %in% datacols)){ 
        stop(name, " is missing some column(s): ", 
             paste(setdiff(columnNames, datacols), collapse=", "))
    }
    areNA <- apply(is.na(data[,columnNames]), 1, any)
    if(sum(areNA)>0L) {
        warning(sum(areNA), " rows have missing values for columns ",
                paste(columnNames, collapse=", "), 
                "and will be removed.")
        data <- data[!areNA,]
    }
    data[,columnNames[1]] <- as.character(sub("^chr", "",
                                          as.character(data[,columnNames[1]])))
   
    if(rmrnd) {
        #Remove bad chrom
        bad_chroms <-grepl("_",  data[,columnNames[1]])
        if( length ( bad_chroms )>0L )
            data <- data[!bad_chroms,]
    }
    
    if(!is.null(optional.columns)) {
        for(i in 1: length(optional.columns)) {
            if(optional.columns[i] %in% datacols)
                data[,optional.columns[i]] <- as.character(data[,optional.columns[i]])            
        }
    }

    return(data)
}


selectchr <- function(chrom, data) {
    if( any (!chrom %in% names(data))) 
        stop("Some chromsomes in 'chr' cannot be drwan.") 
        data <- data[names (data) %in% chrom]
    return (data)
}
transformAlign <- function(filename = NULL){
    y      <- nrow(filename) # size file
    row    <- seq(1, y, by=3) # vector with the number of rows that
    ## contains chromosomal position
    data <- filename[row, ]
    data<-as.character(data)
    data <- strsplit(data, " ", fixed=TRUE)
    matrix <- t(sapply(data, unlist)) #list to matrix
    aux    <-matrix[, c(2:5)]
    dataBed <- as.data.frame(aux) 
    colnames(dataBed) <- c("Chrom", "Start", "End", "Group")
    dataBed$Start     <-as.numeric(as.character(dataBed$Start))
    dataBed$End       <-as.numeric(as.character(dataBed$End))
    dataBed$Group     <-as.character(dataBed$Group) 
    dataBed$Chrom     <-sub("^chr([^:digit:]|[M|X|Y|].+)", "\\1", dataBed$Chrom)
       
    return(dataBed)
}
getAnnotBiomaRt <-function(org=NULL, annotBiomaRt=NULL ) {  
    if(!is.null(org) ){ 
        if(annotBiomaRt){

            dataset <-paste(org, "_gene_ensembl", sep="")
            ensembl <- useMart("ENSEMBL_MART_ENSEMBL",dataset=dataset,
                        host="www.ensembl.org")
            annot<-getBM(attributes=c("chromosome_name","start_position", 
            "end_position","ensembl_gene_id", "strand"),
            filters="chromosome_name",
            values=list("1","2","3", "4", "5","6", "7", "8", "9", "10","11",
            "12","13","14","15","16","17", "18","19", "20", "21", "22", "X",
            "Y"), mart = ensembl)
            
            colnames(annot)[1]<-"Chrom"
            colnames(annot)[2]<-"Start"
            colnames(annot)[3]<-"End"
            colnames(annot)[4]<-"Name"
            colnames(annot)[5]<-"Strand"    
            return(annot)
        } 
    }
}


#For plotting ID from tables
unstacked<-function(x=NULL, y=NULL, trackId=NULL, cex=NULL){

    rest      <- NULL
    trackId$x <- x
    trackId$y <- y
 
    for(i in 1:nrow(trackId)){
        rest <- 0  
        if(i<=nrow(trackId)-1){
            rest<-trackId[i+1,"mids"] - trackId[i,"mids"]
            offset <- trackId[i,"mids"] + abs(strheight(trackId[i,"ID"], cex=cex)*1.05)- trackId[i+1,"mids"]
            if(offset>0){
              trackId[i+1, "y"]<- trackId[i+1, "y"]+offset    
            } 
        }
     
    }
    return(trackId)
}
cytoBands<- function(data=NULL){
    if(!"Colors" %in% names(data)) {
        if("gieStain" %in% names(data)) {
            data        <- data[!is.na(data[ ,"gieStain"]), ]
            type.bands <- match(data$gieStain, c("acen", "gneg", "gpos",
            "gvar", "stalk", "gpos25", "gpos50","gpos66", "gpos75", "gpos100"))
            bandcol<- c("#CCCCCC", "grey96", "#000000", "grey96",  "#CCCCCC", 
                    "grey90", "grey70","#666666", "grey40", "grey0")[type.bands]
            data$Colors<-bandcol
        } else{
                stop("The parameter bands does not contain a color column or",
                " Giemsa stain results (UCSC format).")
          }
    } else{
        if("Group" %in% names(data)) {
            data<-NULL
        }
    }
    return(data)
}

