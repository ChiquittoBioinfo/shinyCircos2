## setwd("absolute path of a directory containing the input data files")
options(warn=-1)
library(circlize)
library(RColorBrewer)
library(GenomicRanges)
library(data.table)
library(RLumShiny)
library(grDevices)

plotcircos <- function(x, color, height, plotTypes, units, rotation, gap.width, labeltextchr, poslabelschr, heightlabelschr, marginlabelschr, data.CN){
  circos.par("start.degree"=90-rotation, "gap.degree"=gap.width, cell.padding=c(0,0,0,0), track.margin=c(0,0))
  circos.genomicInitialize.new(x,plotType=plotTypes,unit=units)
  if(!is.null(data.CN) && ncol(data.CN)==4 && labeltextchr==1 && poslabelschr=="outer"){
    circos.genomicLabels(data.CN, labels.column=4, connection_height=heightlabelschr, track.margin=c(0.01,marginlabelschr), side="outside")
  }		
  circos.genomicTrackPlotRegion(ylim = c(0, 1),bg.col = color, bg.border = NA, track.height = height)	
  if(!is.null(data.CN) && ncol(data.CN)==4 && labeltextchr==1 && poslabelschr=="inner"){
    circos.genomicLabels(data.CN, labels.column=4, connection_height=heightlabelschr, track.margin=c(0.01,marginlabelschr), side="inside")
  }		
}

plotcircos.notrack <- function(x, plotTypes, units, rotation, gap.width, data.CN, labeltextchr, poslabelschr, heightlabelschr, marginlabelschr){
  circos.par("start.degree"=90-rotation, "gap.degree"=gap.width, cell.padding=c(0,0,0,0), track.margin=c(0,0))
  circos.genomicInitialize.new(x,plotType=plotTypes,unit=units)
  if(!is.null(data.CN) && ncol(data.CN)==4 && labeltextchr==1 && poslabelschr=="inner"){
    circos.genomicLabels(data.CN, labels.column = 4, connection_height = heightlabelschr, track.margin = c(0.01,marginlabelschr), side = "inside")
  }
}

plotcircos.font <- function(x, color, height, plotTypes, units, rotation, gap.width, cexLabel, labeltextchr, poslabelschr, heightlabelschr, marginlabelschr, data.CN){
  circos.par("start.degree"=90-rotation, "gap.degree"=gap.width, cell.padding=c(0,0,0,0), track.margin=c(0,0))
  circos.genomicInitialize.new.font(x, plotType=plotTypes, unit=units, cexlabel=cexLabel)
  if(!is.null(data.CN) && ncol(data.CN)==4 && labeltextchr==1 && poslabelschr=="outer"){
    circos.genomicLabels(data.CN, labels.column = 4, connection_height = heightlabelschr, track.margin = c(0.01,marginlabelschr), side = "outside")
  }	
  circos.genomicTrackPlotRegion(ylim = c(0, 1),bg.col = color, bg.border = NA, track.height = height)
  if(!is.null(data.CN) && ncol(data.CN)==4 && labeltextchr==1 && poslabelschr=="inner"){
    circos.genomicLabels(data.CN, labels.column = 4, connection_height = heightlabelschr, track.margin = c(0.01,marginlabelschr), side = "inside")
  }	
}

plotcircos.notrack.font <- function(x, plotTypes, units, rotation, gap.width, cexLabel, data.CN, labeltextchr, poslabelschr, heightlabelschr, marginlabelschr){  
  circos.par("start.degree"=90-rotation, "gap.degree"=gap.width, cell.padding=c(0,0,0,0), track.margin=c(0,0))
  circos.genomicInitialize.new.font(x, plotType=plotTypes, unit=units, cexlabel=cexLabel)
  if(!is.null(data.CN) && ncol(data.CN)==4 && labeltextchr==1 && poslabelschr=="inner"){
    circos.genomicLabels(data.CN, labels.column = 4, connection_height = heightlabelschr, track.margin = c(0.01,marginlabelschr), side = "inside")
  }	
}

plotcircos.cyto <- function(x, height, plotTypes, units, rotation, gap.width, labeltextchr, poslabelschr, heightlabelschr, marginlabelschr, data.CN){ 
  circos.par("start.degree"=90-rotation, "gap.degree"=gap.width, cell.padding=c(0,0,0,0), track.margin=c(0,0))
  circos.genomicInitialize.new(x, plotType = plotTypes, unit=units)
  if(!is.null(data.CN) && ncol(data.CN)==4 && labeltextchr==1 && poslabelschr=="outer"){
    circos.genomicLabels(data.CN, labels.column = 4, connection_height = heightlabelschr, track.margin = c(0.01,marginlabelschr), side = "outside")
  }	
  circos.genomicTrackPlotRegion(x, ylim = c(0, 1), bg.border = NA, 
                                track.height = height, panel.fun = function(region, value, ...){
                                  col = cytoband.col(value[[2]])
                                  circos.genomicRect(region, value, ybottom = 0, 
                                                     ytop = 1, col = col, border = NA, ...)
                                  xlim = get.cell.meta.data("xlim")
                                  circos.rect(xlim[1], 0, xlim[2], 1, border = "black")
                                }, cell.padding = c(0, 0, 0, 0))  
  if(!is.null(data.CN) && ncol(data.CN)==4 && labeltextchr==1 && poslabelschr=="inner"){
    circos.genomicLabels(data.CN, labels.column = 4, connection_height = heightlabelschr, track.margin = c(0.01,marginlabelschr), side = "inside")
  }	
}

plotcircos.cyto.font <- function(x, height, plotTypes, units, rotation, gap.width, cexLabel, labeltextchr, poslabelschr, heightlabelschr, marginlabelschr, data.CN){
  circos.par("start.degree"=90-rotation, "gap.degree"=gap.width, cell.padding=c(0,0,0,0), track.margin=c(0,0))
  circos.genomicInitialize.new.font(x, plotType=plotTypes, unit=units, cexlabel=cexLabel)  
  if(!is.null(data.CN) && ncol(data.CN)==4 && labeltextchr==1 && poslabelschr=="outer"){
    circos.genomicLabels(data.CN, labels.column = 4, connection_height = heightlabelschr, track.margin = c(0.01,marginlabelschr), side = "outside")
  }	
  circos.genomicTrackPlotRegion(x, ylim = c(0, 1), bg.border = NA, 
                                track.height = height, panel.fun = function(region, value, ...){
                                  col = cytoband.col(value[[2]])
                                  circos.genomicRect(region, value, ybottom = 0, 
                                                     ytop = 1, col = col, border = NA, ...)
                                  xlim = get.cell.meta.data("xlim")
                                  circos.rect(xlim[1], 0, xlim[2], 1, border = "black")
                                }, cell.padding = c(0, 0, 0, 0))
  if(!is.null(data.CN) && ncol(data.CN)==4 && labeltextchr==1 && poslabelschr=="inner"){
    circos.genomicLabels(data.CN, labels.column = 4, connection_height = heightlabelschr, track.margin = c(0.01,marginlabelschr), side = "inside")
  }	
}

circos.genomicInitialize.new <- 
  function (data, sector.names = NULL, major.by = NULL, unit = "", plotType, tickLabelsStartFromZero = TRUE, track.height = 0.05, 
            ...) 
  {
    if(is.factor(data[[1]])){
      fa = levels(data[[1]])
    }
    else {
      fa = unique(data[[1]])
    }
    if(!is.null(sector.names)){
      if(length(sector.names) != length(fa)){
        stop("length of `sector.names` and length of sectors differ.")
      }
    }
    else {
      sector.names = fa
    }
    names(sector.names) = fa
    x1 = tapply(data[[2]], data[[1]], min)[fa]
    x2 = tapply(data[[3]], data[[1]], max)[fa]
    op = circos.par("cell.padding")
    ow = circos.par("points.overflow.warning")
    circos.par(cell.padding = c(0, 0, 0, 0), points.overflow.warning = FALSE)
    circos.initialize(factor(fa, levels = fa), xlim = cbind(x1, 
                                                            x2), ...)
    if(any(plotType %in% c("axis", "labels"))){
      circos.genomicTrackPlotRegion(data, ylim = c(0, 1), bg.border = NA, 
                                    track.height = track.height, panel.fun = function(region, 
                                                                                      value, ...){
                                      sector.index = get.cell.meta.data("sector.index")
                                      xlim = get.cell.meta.data("xlim")
                                      if(tickLabelsStartFromZero){
                                        offset = xlim[1]
                                        if(is.null(major.by)){
                                          xlim = get.cell.meta.data("xlim")
                                          major.by = .default.major.by()
                                        }
                                        major.at = seq(xlim[1], xlim[2], by = major.by)
                                        major.at = c(major.at, major.at[length(major.at)] + 
                                                       major.by)
                                        if(major.by > 1e+06){
                                          major.tick.labels = paste((major.at - offset)/1e+06, 
                                                                    "MB", sep = "")
                                        }
                                        else if(major.by > 1000){
                                          major.tick.labels = paste((major.at - offset)/1000, 
                                                                    "KB", sep = "")
                                        }
                                        else {
                                          major.tick.labels = paste((major.at - offset), 
                                                                    "bp", sep = "")
                                        }
                                      }
                                      else {
                                        if(is.null(major.by)){
                                          xlim = get.cell.meta.data("xlim")
                                          major.by = .default.major.by()
                                        }
                                        major.at = seq(floor(xlim[1]/major.by) * major.by, 
                                                       xlim[2], by = major.by)
                                        major.at = c(major.at, major.at[length(major.at)] + 
                                                       major.by)
                                        if(major.by > 1e+06){
                                          major.tick.labels = paste(major.at/1e+06, 
                                                                    "MB", sep = "")
                                        }
                                        else if(major.by > 1000){
                                          major.tick.labels = paste(major.at/1000, 
                                                                    "KB", sep = "")
                                        }
                                        else {
                                          major.tick.labels = paste(major.at, "bp", 
                                                                    sep = "")
                                        }
                                      }
                                      
                                      if(unit==""){ major.tick.labels <- gsub("[mkbp]","",major.tick.labels,ignore.case = T)}
                                      
                                      if(all(c("axis", "labels") %in% plotType)){
                                        circos.axis(h = 0, major.at = major.at, labels = major.tick.labels, 
                                                    labels.cex = 0.49 * par("cex"), labels.facing = "clockwise", 
                                                    major.tick.percentage = 0.2)
                                        circos.text(mean(xlim), 1.2, labels = sector.names[sector.index], 
                                                    cex = par("cex")-0.1, adj = c(0.5, -0.1*par("cex")*6-(par("cex")-1)*3), niceFacing = TRUE)
                                      }
                                      else if("labels" %in% plotType){
                                        circos.text(mean(xlim), 0, labels = sector.names[sector.index], 
                                                    cex = par("cex")-0.1, adj = c(0.5, -0.1*par("cex")*6-(par("cex")-1)*3), niceFacing = TRUE)
                                      }
                                      else if("axis" %in% plotType){
                                        circos.axis(h = 0, major.at = major.at, labels = major.tick.labels, 
                                                    labels.cex = 0.49 * par("cex"), labels.facing = "clockwise", 
                                                    major.tick.percentage = 0.2)
                                      }
                                    })
    }
    circos.par(cell.padding = op, points.overflow.warning = ow)
    return(invisible(NULL))
  }
  
circos.genomicInitialize.new.font <- 
  function (data, sector.names = NULL, major.by = NULL, unit = "", plotType, tickLabelsStartFromZero = TRUE, track.height = 0.05, cexlabel, 
            ...) 
  {
    if(is.factor(data[[1]])){
      fa = levels(data[[1]])
    }
    else {
      fa = unique(data[[1]])
    }
    if(!is.null(sector.names)){
      if(length(sector.names) != length(fa)){
        stop("length of `sector.names` and length of sectors differ.")
      }
    }
    else {
      sector.names = fa
    }
    names(sector.names) = fa
    x1 = tapply(data[[2]], data[[1]], min)[fa]
    x2 = tapply(data[[3]], data[[1]], max)[fa]
    op = circos.par("cell.padding")
    ow = circos.par("points.overflow.warning")
    circos.par(cell.padding = c(0, 0, 0, 0), points.overflow.warning = FALSE)
    circos.initialize(factor(fa, levels = fa), xlim = cbind(x1, 
                                                            x2), ...)
    if(any(plotType %in% c("axis", "labels"))){
      circos.genomicTrackPlotRegion(data, ylim = c(0, 1), bg.border = NA, 
                                    track.height = track.height, panel.fun = function(region, 
                                                                                      value, ...){
                                      sector.index = get.cell.meta.data("sector.index")
                                      xlim = get.cell.meta.data("xlim")
                                      if(tickLabelsStartFromZero){
                                        offset = xlim[1]
                                        if(is.null(major.by)){
                                          xlim = get.cell.meta.data("xlim")
                                          major.by = .default.major.by()
                                        }
                                        major.at = seq(xlim[1], xlim[2], by = major.by)
                                        major.at = c(major.at, major.at[length(major.at)] + 
                                                       major.by)
                                        if(major.by > 1e+06){
                                          major.tick.labels = paste((major.at - offset)/1e+06, 
                                                                    "MB", sep = "")
                                        }
                                        else if(major.by > 1000){
                                          major.tick.labels = paste((major.at - offset)/1000, 
                                                                    "KB", sep = "")
                                        }
                                        else {
                                          major.tick.labels = paste((major.at - offset), 
                                                                    "bp", sep = "")
                                        }
                                      }
                                      else {
                                        if(is.null(major.by)){
                                          xlim = get.cell.meta.data("xlim")
                                          major.by = .default.major.by()
                                        }
                                        major.at = seq(floor(xlim[1]/major.by) * major.by, 
                                                       xlim[2], by = major.by)
                                        major.at = c(major.at, major.at[length(major.at)] + 
                                                       major.by)
                                        if(major.by > 1e+06){
                                          major.tick.labels = paste(major.at/1e+06, 
                                                                    "MB", sep = "")
                                        }
                                        else if(major.by > 1000){
                                          major.tick.labels = paste(major.at/1000, 
                                                                    "KB", sep = "")
                                        }
                                        else {
                                          major.tick.labels = paste(major.at, "bp", 
                                                                    sep = "")
                                        }
                                      }
                                      
                                      if(unit==""){ major.tick.labels <- gsub("[mkbp]","",major.tick.labels,ignore.case = T)}
									  
                                      if(all(c("axis", "labels") %in% plotType)){
                                        circos.axis(h = 0, major.at = major.at, labels = major.tick.labels, 
                                                    labels.cex = 0.49 * cexlabel, labels.facing = "clockwise", 
                                                    major.tick.percentage = 0.2)
                                        circos.text(mean(xlim), 1.2, labels = sector.names[sector.index], 
                                                    cex = cexlabel, adj = c(0.5, -0.1*cexlabel*6-(cexlabel-1)*3), niceFacing = TRUE)
                                      }
                                      else if("labels" %in% plotType){
                                        circos.text(mean(xlim), 0, labels = sector.names[sector.index], 
                                                    cex = cexlabel, adj = c(0.5, -0.1*cexlabel*6-(cexlabel-1)*3), niceFacing = TRUE)
                                      }
                                      else if("axis" %in% plotType){
                                        circos.axis(h = 0, major.at = major.at, labels = major.tick.labels, 
                                                    labels.cex = 0.49 * cexlabel, labels.facing = "clockwise", 
                                                    major.tick.percentage = 0.2)
                                      }
                                    })
    }
    circos.par(cell.padding = op, points.overflow.warning = ow)
    return(invisible(NULL))
  }
  
.default.major.by = function(sector.index = get.cell.meta.data("sector.index"),
	track.index = get.cell.meta.data("track.index")){
	d = circos.par("major.by.degree")
	cell.start.degre = get.cell.meta.data("cell.start.degree", sector.index, track.index)
	tm = reverse.circlize(c(cell.start.degre, cell.start.degre-d), rep(get.cell.meta.data("cell.bottom.radius", sector.index = sector.index, track.index = track.index), 2))
	major.by = abs(tm[1, 1] - tm[2, 1])
	digits = as.numeric(gsub("^.*e([+-]\\d+)$", "\\1", sprintf("%e", major.by)))
	major.by = round(major.by, digits = -1*digits)
	return(major.by)
}

get_most_inside_radius = function() {
	tracks = get.all.track.index()
	if(length(tracks) == 0) {
	   1
	}else{
	   n = length(tracks)
	   get.cell.meta.data("cell.bottom.radius", track.index = tracks[n]) - get.cell.meta.data("track.margin", track.index = tracks[n])[1] - circos.par("track.margin")[2]
	}
}

data.C.name <- "chr.csv"
data.C <- data.frame(fread(data.C.name),stringsAsFactors=F)
data.C[,2] <- as.numeric(data.C[,2])
data.C[,3] <- as.numeric(data.C[,3])
data.T.file <- c("chart1_mRNA-DensityMap.csv","chart1_ncRNA-DensityMap.csv","chart1_TE-DensityMap.csv")
data.T <- lapply(1:length(data.T.file),function(x){
		  if(!is.null(data.T.file[x])){
		  data.frame(fread(data.T.file[x]),stringsAsFactors=F)
		  }
		  })
data.CN <- NULL
data.N.file <- c("","","","","","","","","","")
uploadtrack <- c(2,2,2,1,1,1,1,1,1,1)
data.N <- lapply(1:10,function(x){
			 if(uploadtrack[x] == 2 && nchar(data.N.file[x])>0){	  
		     data.frame(fread(data.N.file[x]),stringsAsFactors=F)
			 }
			 })
trackindx <- c(1,2,3)
data.N <- data.N[trackindx]
data.L <- NULL
for(i in 1:length(data.T.file)){
  assign(paste("hltdata",i,sep=""),"")
}
hltdata1 <- ""
hltdata2 <- ""
hltdata3 <- ""
hltregion.List <- list()
if(!is.null(data.T)){
			for(k in 1:length(data.T)){
			data.TT <- data.T[[k]]
			hltregion.List[[k]] <- ""
if(nchar(get(paste("hltdata",k,sep="")))>0){
tmp <- matrix(strsplit(get(paste("hltdata",k,sep="")), "\n")[[1]])
            myColnames <- c("chr","start","end","color")
            data <- matrix(0, length(tmp), length(myColnames))
            colnames(data) <- myColnames
            for(p in 1:length(tmp)){
                 myRow <- strsplit(tmp[p], ",")[[1]]
                  if(length(myRow)==4){                                        
                    data[p,] <- myRow
                  }
               }
            data <- data.frame(data,stringsAsFactors=F)
            data$start <- as.numeric(data$start)
            data$end <- as.numeric(data$end)
			query <- GRanges(seqnames = data$chr,ranges=IRanges(start=data$start,end=data$end),seqinfo=NULL)
            subj <- GRanges(seqnames = data.TT[,1],ranges=IRanges(start=data.TT[,2],end=data.TT[,3]),seqinfo=NULL) 
            indx <- findOverlaps(query,subj)
            indx <- data.frame(indx,stringsAsFactors=F)
			indx$queryHits <- as.numeric(indx$queryHits)
			indx$subjectHits <- as.numeric(indx$subjectHits)
            hltregion <- data.TT[indx$subjectHits,]
			hltregion$color <- data$color[indx[,1]]
			hltregion$id <- paste(hltregion[,1],hltregion[,2],hltregion[,3],sep="")
			hltregion.List[[k]] <- hltregion
			}
			}
			}

pdf("shinyCircos.pdf", width=750/72, height=750/72)
## svg("shinyCircos.svg", width=750/72, height=750/72)
fontSize <- 1
par(oma=c(0,0,0,0), mar=c(9,0.5,1,9.5), xpd=TRUE, cex=fontSize-0.05)
trackChr <- "track"
plotTypes <- "labels"
plotTypes <- "axis"
unitChr <- "unit"
rotation <- 0.5
gap.width <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
labeltextchr <- 2
poslabelschr <- "inner"
heightlabelschr <- 0.06
marginlabelschr <- 0.01
colorChr <- c("grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey")
heightChr <- 0.05
plotcircos(data.C, height=heightChr, color=colorChr, plotTypes=plotTypes, units=unitChr, rotation=rotation, gap.width=gap.width, labeltextchr=labeltextchr, poslabelschr=poslabelschr, heightlabelschr=heightlabelschr, marginlabelschr=marginlabelschr, data.CN=data.CN)
takindx <- 1
typeTrack <- c("bar","bar","bar")
i <- 1
data.TT <- data.T[[i]]
	tktype <- typeTrack[i]
	data.TT[,2] <- as.numeric(data.TT[,2])
	data.TT[,3] <- as.numeric(data.TT[,3])
	data.NN <- data.N[[i]]
	data.TT$num <- 1:nrow(data.TT)
data.TTC <- NULL
coltypeTrack <- 2
tkcolor <- c("black")
data.TT$num <- NULL
tkbgcol <- c("grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95")
tkmargin <- 0.01
tkheight <- 0.06
tklinecoord <- c(0.25,0.75)
tklinecolor <- c("grey","grey")
hmapcols <- c("blue","white","red")
tkborder <- ""
innergap <- 0.5
tkbordercol <- NA
tkbardir <- 1
tkrectcol <- 1
selrectcol <- 1
rectcols <- c("#EDEDFD","#6969F5","#00008B")
tktransparency <- 1
tkcolor <- c("#000000FF")
data.TTT <- data.T[[i]]
	data.TTT$id <- paste(data.TTT[,1],data.TTT[,2],data.TTT[,3],sep="")
	data.TTT$num <- 1:nrow(data.TTT)
transparencyHlt <- c(1,1,1)
lkmargin <- 0
tkborder <- NA
columns <- c(4)
data.TT[,ncol(data.TT)] <- as.numeric(data.TT[,ncol(data.TT)])
			circos.genomicTrackPlotRegion(data.TT, track.height = tkheight, track.margin = c(lkmargin,tkmargin), bg.col = tkbgcol, bg.border = tkborder, panel.fun = function(region,value,...){
			        if(!("color" %in% colnames(data.T[[i]])) && !("cols" %in% colnames(data.TTC))){
			        if(length(columns)==1 && tkbardir==1){
					if(nchar(tklinecolor[1])!=0){               
					  xlim <- get.cell.meta.data("xlim")
                      ylim <- get.cell.meta.data("ylim")
					  for(k in 1:length(tklinecoord)){
                      y1 <- as.numeric(quantile(ylim,probs=tklinecoord[k]))
                      circos.lines(x=xlim,y=c(y1,y1), col=tklinecolor[k], lwd=0.1)
					  }
				    }
					if(coltypeTrack==1){
                       circos.genomicRect(region, value, numeric.column=columns-3, ytop.column = 1, ybottom = min(data.TT[,4]), col=tkcolor, border = NA, ...)
					}
                    }else if(length(columns)==1 && tkbardir==2){
					if(nchar(tklinecolor[1])!=0){               
					  xlim <- get.cell.meta.data("xlim")
                      ylim <- get.cell.meta.data("ylim")
					  for(k in 1:length(tklinecoord)){
                      y1 <- as.numeric(quantile(ylim,probs=tklinecoord[k]))
                      circos.lines(x=xlim,y=c(y1,y1), col=tklinecolor[k], lwd=0.1)
					  }
				    }
				    tkbarvalue <- as.numeric(tkbarvalue)
					indx <- value[,1] > tkbarvalue
                    if(length(value[indx,])!=0 && length(value[!indx,])!=0){
                         circos.genomicRect(region[indx,], value[indx,], ytop.column = 1, ybottom = tkbarvalue, col=tkbarcol1, border = NA, ...)
                         circos.genomicRect(region[!indx,], value[!indx,], ytop.column = 1, ybottom =  tkbarvalue, col=tkbarcol2, border = NA, ...)
                    }else if(length(value[indx,])!=0 && length(value[!indx,])==0){
                         circos.genomicRect(region[indx,], value[indx,], ytop.column = 1, ybottom = tkbarvalue, col=tkbarcol1, border = NA, ...)
                    }else if(length(value[indx,])==0 && length(value[!indx,])!=0){
                         circos.genomicRect(region[!indx,], value[!indx,], ytop.column = 1, ybottom =  tkbarvalue, col=tkbarcol2, border = NA, ...)
                    }
					}else if(length(columns)==2 && tkbardir==1){
					if(nchar(tklinecolor[1])!=0){               
					   xlim <- get.cell.meta.data("xlim")
                       ylim <- get.cell.meta.data("ylim")
					   for(k in 1:length(tklinecoord)){
                       y1 <- as.numeric(quantile(ylim,probs=tklinecoord[k]))
                       circos.lines(x=xlim,y=c(y1,y1), col=tklinecolor[k], lwd=0.1)
					   }
				    }
				  if(coltypeTrack!=3){					
					    tkcolor <- c(tkcolor,rep("grey",length(columns)))
					    tkcolor <- tkcolor[1:length(columns)]				
						circos.genomicRect(region, value, ytop.column = 1, ybottom = min(c(data.TT[,4],data.TT[,5])), col=tkcolor[1], border = NA, ...)
                        circos.genomicRect(region, value, ytop.column = 2, ybottom = min(c(data.TT[,4],data.TT[,5])), col=tkcolor[2], border = NA, ...)
                  }						
					}else if(length(columns)==2 && tkbardir==2){
					if(nchar(tklinecolor[1])!=0){               
					  xlim <- get.cell.meta.data("xlim")
                      ylim <- get.cell.meta.data("ylim")
					  for(k in 1:length(tklinecoord)){
                      y1 <- as.numeric(quantile(ylim,probs=tklinecoord[k]))
                      circos.lines(x=xlim,y=c(y1,y1), col=tklinecolor[k], lwd=0.1)
					  }
				    }
						circos.genomicRect(region, value, ytop.column = 1, ybottom = min(c(data.TT[,4],data.TT[,5])), col=tkbarcol1, border = NA, ...)
                        circos.genomicRect(region, value, ytop.column = 2, ybottom = min(c(data.TT[,4],data.TT[,5])), col=tkbarcol2, border = NA, ...)						
					}
					}
			        if(length(columns)==1 && tkbardir==1 && coltypeTrack==2){
					if(nchar(tklinecolor[1])!=0){               
					  xlim <- get.cell.meta.data("xlim")
                      ylim <- get.cell.meta.data("ylim")
					  for(k in 1:length(tklinecoord)){
                      y1 <- as.numeric(quantile(ylim,probs=tklinecoord[k]))
                      circos.lines(x=xlim,y=c(y1,y1), col=tklinecolor[k], lwd=0.1)
					  }
				    }
                      circos.genomicRect(region, value, numeric.column=columns-3, ytop.column = 1, ybottom = min(data.TT[,4]), col=tkcolor[1], border = NA, ...)
                    }	
					})
poslabels <- c("outer","outer","outer")
if(poslabels[i]=="inner"){
			    takindx <- takindx+3
			}else{
			    takindx <- takindx+1
			}
i <- 2
data.TT <- data.T[[i]]
	tktype <- typeTrack[i]
	data.TT[,2] <- as.numeric(data.TT[,2])
	data.TT[,3] <- as.numeric(data.TT[,3])
	data.NN <- data.N[[i]]
	data.TT$num <- 1:nrow(data.TT)
data.TTC <- NULL
coltypeTrack <- 2
tkcolor <- c("red")
data.TT$num <- NULL
tkbgcol <- c("grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95")
tkmargin <- 0.01
tkheight <- 0.06
tklinecoord <- c(0.25,0.75)
tklinecolor <- c("grey","grey")
hmapcols <- c("blue","white","red")
tkborder <- ""
innergap <- 0.5
tkbordercol <- NA
tkbardir <- 1
tkrectcol <- 1
selrectcol <- 1
rectcols <- c("#EDEDFD","#6969F5","#00008B")
tktransparency <- 1
tkcolor <- c("#FF0000FF")
data.TTT <- data.T[[i]]
	data.TTT$id <- paste(data.TTT[,1],data.TTT[,2],data.TTT[,3],sep="")
	data.TTT$num <- 1:nrow(data.TTT)
transparencyHlt <- c(1,1,1)
lkmargin <- 0
tkborder <- NA
columns <- c(4)
data.TT[,ncol(data.TT)] <- as.numeric(data.TT[,ncol(data.TT)])
			circos.genomicTrackPlotRegion(data.TT, track.height = tkheight, track.margin = c(lkmargin,tkmargin), bg.col = tkbgcol, bg.border = tkborder, panel.fun = function(region,value,...){
			        if(!("color" %in% colnames(data.T[[i]])) && !("cols" %in% colnames(data.TTC))){
			        if(length(columns)==1 && tkbardir==1){
					if(nchar(tklinecolor[1])!=0){               
					  xlim <- get.cell.meta.data("xlim")
                      ylim <- get.cell.meta.data("ylim")
					  for(k in 1:length(tklinecoord)){
                      y1 <- as.numeric(quantile(ylim,probs=tklinecoord[k]))
                      circos.lines(x=xlim,y=c(y1,y1), col=tklinecolor[k], lwd=0.1)
					  }
				    }
					if(coltypeTrack==1){
                       circos.genomicRect(region, value, numeric.column=columns-3, ytop.column = 1, ybottom = min(data.TT[,4]), col=tkcolor, border = NA, ...)
					}
                    }else if(length(columns)==1 && tkbardir==2){
					if(nchar(tklinecolor[1])!=0){               
					  xlim <- get.cell.meta.data("xlim")
                      ylim <- get.cell.meta.data("ylim")
					  for(k in 1:length(tklinecoord)){
                      y1 <- as.numeric(quantile(ylim,probs=tklinecoord[k]))
                      circos.lines(x=xlim,y=c(y1,y1), col=tklinecolor[k], lwd=0.1)
					  }
				    }
				    tkbarvalue <- as.numeric(tkbarvalue)
					indx <- value[,1] > tkbarvalue
                    if(length(value[indx,])!=0 && length(value[!indx,])!=0){
                         circos.genomicRect(region[indx,], value[indx,], ytop.column = 1, ybottom = tkbarvalue, col=tkbarcol1, border = NA, ...)
                         circos.genomicRect(region[!indx,], value[!indx,], ytop.column = 1, ybottom =  tkbarvalue, col=tkbarcol2, border = NA, ...)
                    }else if(length(value[indx,])!=0 && length(value[!indx,])==0){
                         circos.genomicRect(region[indx,], value[indx,], ytop.column = 1, ybottom = tkbarvalue, col=tkbarcol1, border = NA, ...)
                    }else if(length(value[indx,])==0 && length(value[!indx,])!=0){
                         circos.genomicRect(region[!indx,], value[!indx,], ytop.column = 1, ybottom =  tkbarvalue, col=tkbarcol2, border = NA, ...)
                    }
					}else if(length(columns)==2 && tkbardir==1){
					if(nchar(tklinecolor[1])!=0){               
					   xlim <- get.cell.meta.data("xlim")
                       ylim <- get.cell.meta.data("ylim")
					   for(k in 1:length(tklinecoord)){
                       y1 <- as.numeric(quantile(ylim,probs=tklinecoord[k]))
                       circos.lines(x=xlim,y=c(y1,y1), col=tklinecolor[k], lwd=0.1)
					   }
				    }
				  if(coltypeTrack!=3){					
					    tkcolor <- c(tkcolor,rep("grey",length(columns)))
					    tkcolor <- tkcolor[1:length(columns)]				
						circos.genomicRect(region, value, ytop.column = 1, ybottom = min(c(data.TT[,4],data.TT[,5])), col=tkcolor[1], border = NA, ...)
                        circos.genomicRect(region, value, ytop.column = 2, ybottom = min(c(data.TT[,4],data.TT[,5])), col=tkcolor[2], border = NA, ...)
                  }						
					}else if(length(columns)==2 && tkbardir==2){
					if(nchar(tklinecolor[1])!=0){               
					  xlim <- get.cell.meta.data("xlim")
                      ylim <- get.cell.meta.data("ylim")
					  for(k in 1:length(tklinecoord)){
                      y1 <- as.numeric(quantile(ylim,probs=tklinecoord[k]))
                      circos.lines(x=xlim,y=c(y1,y1), col=tklinecolor[k], lwd=0.1)
					  }
				    }
						circos.genomicRect(region, value, ytop.column = 1, ybottom = min(c(data.TT[,4],data.TT[,5])), col=tkbarcol1, border = NA, ...)
                        circos.genomicRect(region, value, ytop.column = 2, ybottom = min(c(data.TT[,4],data.TT[,5])), col=tkbarcol2, border = NA, ...)						
					}
					}
			        if(length(columns)==1 && tkbardir==1 && coltypeTrack==2){
					if(nchar(tklinecolor[1])!=0){               
					  xlim <- get.cell.meta.data("xlim")
                      ylim <- get.cell.meta.data("ylim")
					  for(k in 1:length(tklinecoord)){
                      y1 <- as.numeric(quantile(ylim,probs=tklinecoord[k]))
                      circos.lines(x=xlim,y=c(y1,y1), col=tklinecolor[k], lwd=0.1)
					  }
				    }
                      circos.genomicRect(region, value, numeric.column=columns-3, ytop.column = 1, ybottom = min(data.TT[,4]), col=tkcolor[1], border = NA, ...)
                    }	
					})
poslabels <- c("outer","outer","outer")
if(poslabels[i]=="inner"){
			    takindx <- takindx+3
			}else{
			    takindx <- takindx+1
			}
i <- 3
data.TT <- data.T[[i]]
	tktype <- typeTrack[i]
	data.TT[,2] <- as.numeric(data.TT[,2])
	data.TT[,3] <- as.numeric(data.TT[,3])
	data.NN <- data.N[[i]]
	data.TT$num <- 1:nrow(data.TT)
data.TTC <- NULL
coltypeTrack <- 2
tkcolor <- c("blue")
data.TT$num <- NULL
tkbgcol <- c("grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95","grey95")
tkmargin <- 0.01
tkheight <- 0.06
tklinecoord <- c(0.25,0.75)
tklinecolor <- c("grey","grey")
hmapcols <- c("blue","white","red")
tkborder <- ""
innergap <- 0.5
tkbordercol <- NA
tkbardir <- 1
tkrectcol <- 1
selrectcol <- 1
rectcols <- c("#EDEDFD","#6969F5","#00008B")
tktransparency <- 1
tkcolor <- c("#0000FFFF")
data.TTT <- data.T[[i]]
	data.TTT$id <- paste(data.TTT[,1],data.TTT[,2],data.TTT[,3],sep="")
	data.TTT$num <- 1:nrow(data.TTT)
transparencyHlt <- c(1,1,1)
lkmargin <- 0.01
tkborder <- NA
columns <- c(4)
data.TT[,ncol(data.TT)] <- as.numeric(data.TT[,ncol(data.TT)])
			circos.genomicTrackPlotRegion(data.TT, track.height = tkheight, track.margin = c(lkmargin,tkmargin), bg.col = tkbgcol, bg.border = tkborder, panel.fun = function(region,value,...){
			        if(!("color" %in% colnames(data.T[[i]])) && !("cols" %in% colnames(data.TTC))){
			        if(length(columns)==1 && tkbardir==1){
					if(nchar(tklinecolor[1])!=0){               
					  xlim <- get.cell.meta.data("xlim")
                      ylim <- get.cell.meta.data("ylim")
					  for(k in 1:length(tklinecoord)){
                      y1 <- as.numeric(quantile(ylim,probs=tklinecoord[k]))
                      circos.lines(x=xlim,y=c(y1,y1), col=tklinecolor[k], lwd=0.1)
					  }
				    }
					if(coltypeTrack==1){
                       circos.genomicRect(region, value, numeric.column=columns-3, ytop.column = 1, ybottom = min(data.TT[,4]), col=tkcolor, border = NA, ...)
					}
                    }else if(length(columns)==1 && tkbardir==2){
					if(nchar(tklinecolor[1])!=0){               
					  xlim <- get.cell.meta.data("xlim")
                      ylim <- get.cell.meta.data("ylim")
					  for(k in 1:length(tklinecoord)){
                      y1 <- as.numeric(quantile(ylim,probs=tklinecoord[k]))
                      circos.lines(x=xlim,y=c(y1,y1), col=tklinecolor[k], lwd=0.1)
					  }
				    }
				    tkbarvalue <- as.numeric(tkbarvalue)
					indx <- value[,1] > tkbarvalue
                    if(length(value[indx,])!=0 && length(value[!indx,])!=0){
                         circos.genomicRect(region[indx,], value[indx,], ytop.column = 1, ybottom = tkbarvalue, col=tkbarcol1, border = NA, ...)
                         circos.genomicRect(region[!indx,], value[!indx,], ytop.column = 1, ybottom =  tkbarvalue, col=tkbarcol2, border = NA, ...)
                    }else if(length(value[indx,])!=0 && length(value[!indx,])==0){
                         circos.genomicRect(region[indx,], value[indx,], ytop.column = 1, ybottom = tkbarvalue, col=tkbarcol1, border = NA, ...)
                    }else if(length(value[indx,])==0 && length(value[!indx,])!=0){
                         circos.genomicRect(region[!indx,], value[!indx,], ytop.column = 1, ybottom =  tkbarvalue, col=tkbarcol2, border = NA, ...)
                    }
					}else if(length(columns)==2 && tkbardir==1){
					if(nchar(tklinecolor[1])!=0){               
					   xlim <- get.cell.meta.data("xlim")
                       ylim <- get.cell.meta.data("ylim")
					   for(k in 1:length(tklinecoord)){
                       y1 <- as.numeric(quantile(ylim,probs=tklinecoord[k]))
                       circos.lines(x=xlim,y=c(y1,y1), col=tklinecolor[k], lwd=0.1)
					   }
				    }
				  if(coltypeTrack!=3){					
					    tkcolor <- c(tkcolor,rep("grey",length(columns)))
					    tkcolor <- tkcolor[1:length(columns)]				
						circos.genomicRect(region, value, ytop.column = 1, ybottom = min(c(data.TT[,4],data.TT[,5])), col=tkcolor[1], border = NA, ...)
                        circos.genomicRect(region, value, ytop.column = 2, ybottom = min(c(data.TT[,4],data.TT[,5])), col=tkcolor[2], border = NA, ...)
                  }						
					}else if(length(columns)==2 && tkbardir==2){
					if(nchar(tklinecolor[1])!=0){               
					  xlim <- get.cell.meta.data("xlim")
                      ylim <- get.cell.meta.data("ylim")
					  for(k in 1:length(tklinecoord)){
                      y1 <- as.numeric(quantile(ylim,probs=tklinecoord[k]))
                      circos.lines(x=xlim,y=c(y1,y1), col=tklinecolor[k], lwd=0.1)
					  }
				    }
						circos.genomicRect(region, value, ytop.column = 1, ybottom = min(c(data.TT[,4],data.TT[,5])), col=tkbarcol1, border = NA, ...)
                        circos.genomicRect(region, value, ytop.column = 2, ybottom = min(c(data.TT[,4],data.TT[,5])), col=tkbarcol2, border = NA, ...)						
					}
					}
			        if(length(columns)==1 && tkbardir==1 && coltypeTrack==2){
					if(nchar(tklinecolor[1])!=0){               
					  xlim <- get.cell.meta.data("xlim")
                      ylim <- get.cell.meta.data("ylim")
					  for(k in 1:length(tklinecoord)){
                      y1 <- as.numeric(quantile(ylim,probs=tklinecoord[k]))
                      circos.lines(x=xlim,y=c(y1,y1), col=tklinecolor[k], lwd=0.1)
					  }
				    }
                      circos.genomicRect(region, value, numeric.column=columns-3, ytop.column = 1, ybottom = min(data.TT[,4]), col=tkcolor[1], border = NA, ...)
                    }	
					})
poslabels <- c("outer","outer","outer")
if(poslabels[i]=="inner"){
			    takindx <- takindx+3
			}else{
			    takindx <- takindx+1
			}
n <- 4
xleft <- 1.2
xright <- 1.22
ybottom <- -0.12286
ytop <- 0.11714
len <- 0.24
gap <- 0.075
legendtext <- c("Chromosome","Coding RNA","ncRNA","TE")
for(i in 1:n){
			   assign(paste("n",i,sep=""),legendtext[i])
			}
rect(xleft, ybottom, xright, ytop, col = "black")
            polygon(x=c(xleft-0.01,(xleft+xright)/2,xright+0.01), y=c(ybottom,ybottom-0.02,ybottom), col="black")
            text(x=xleft-0.08, y=ybottom, labels="inner", cex=0.95)
text(x=xleft-0.08, y=ytop-0.02, labels="outer", cex=0.95)
for(i in 1:n){
                    text(x=xright+0.028, y=ytop-gap*(i-1)-0.025, labels=get(paste("n",i,sep="")), cex=1, adj=c(0,0))
                }
dev.off()
circos.clear()
