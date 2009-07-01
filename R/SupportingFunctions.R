CNSeg <- function(segList){
    return(new("CNSeg", segList = segList))
}

RS <- function(rs, by){
    return(new("RS", rs = rs, by = by))
}

seg2RS <- function(segData, by = c("region", "gene", "pair"), 
    input = TRUE, XY = FALSE, geneMap){
  by <- match.arg(by)
  if(by == "gene" & missing(geneMap)){
      stop("Need geneMap when by == \"gene\"")
  }
  if(!XY){
    seg <- segList(segData)
    segList(segData) <- seg[which(!seg[, "chrom"] 
        %in% c("X", "Y", "x", "y")),]
  }
  rs <- switch(by,
		   region = getCommonSegValues(segList(segData), 
                   drop = TRUE),
               gene = collapseSegList(segList(segData), geneMap),
               pair = getPairwise(segList(segData)))
  if(input & by != "pair"){
     rs <- convertRS(rs, sampleStart = ifelse(by == "region", 4, 6)) 
  }
  return(RS(rs, by))
}

## segList - "output" of CBS
## drop - drop reduced segments with all NAs
## by - column name for chromosome
getCommonSegValues <- function(segList, by = "chrom", 
    drop = FALSE, verbose = TRUE){
  splited <- split.data.frame(segList, factor(segList[, by]))
  if(verbose){
  	cat("Processing samples ...")
  }
  segsByChroms <- lapply(splited, findOverlapingSegs)
  matchedSegs <- lapply(segsByChroms, matchSegValues, segList = segList)
  
  segWithValues <- do.call("rbind", args = matchedSegs)
  if(drop){
    toKeep <- apply(segWithValues[, 4:ncol(segWithValues)], 1,
                  FUN = function(x) !all(is.na(x)))
    segWithValues <- segWithValues[toKeep, ]
  }
  if(verbose){
  	cat(" Done\n")
  }
  return(sortByChromNLoc(segWithValues, by1 = "chrom", by2 = "start"))
}

findOverlapingSegs <- function(segListByChrom){
  breakPoints <- sort(unique(c(as.numeric(
                     as.vector(segListByChrom[, "loc.start"])),
                     as.numeric(as.vector(segListByChrom[, "loc.end"])))))
  newSegs <- cbind(c(breakPoints[1], breakPoints[2:(length(breakPoints) - 1)]
                     + 1), breakPoints[2:length(breakPoints)])
  newSegs <- cbind(as.vector(segListByChrom[1, "chrom"]), newSegs)
  colnames(newSegs) <- c("chrom", "start", "end")
  
  return(newSegs)
}


# segsByChrom - ouput of the function findOverLapingSegs
matchSegValues <- function(segsByChrom, segList){
  
  getSegValues <- function(segData){
    
    mappedSegs[as.numeric(mappedSegs[, "start"]) >=
               as.numeric(segData["loc.start"]) &
               as.numeric(mappedSegs[, "end"]) <=
               as.numeric(segData["loc.end"]),
               which(as.character(colnames(mappedSegs)) ==
                     as.character(segData["ID"]))
               ] <<- as.numeric(segData["seg.mean"])
  }
  
  mappedSegs <- cbind(segsByChrom,
                      matrix(NA, ncol = length(unique(segList[, "ID"])),
                             nrow = nrow(segsByChrom)))
  colnames(mappedSegs) <- c(colnames(segsByChrom),
                            as.vector(unique(segList[, "ID"])))
  # only keeps the segements within a given chromosome
  tempSegList <- segList[gsub(" ", "", segList[, "chrom"]) == gsub(" ", "", segsByChrom[1, "chrom"]), ]
  junk <- apply(tempSegList, 1, getSegValues)  
  
  return(mappedSegs)
}


# Sort mapped in order by chromsome and then by location
# A data frame or matrix to be sorted by columns defined by by1 and by2
sortByChromNLoc <- function(sortMe, by1 = "Ch", by2 = "Pos"){
  splited <- split.data.frame(sortMe, factor(sortMe[, by1]))
  sorted <- NULL
  for(chrom in c(1:22, "X", "Y")){
    if(!is.null(splited[[chrom]])){
      sorted <- rbind(sorted,
              splited[[chrom]][order(as.numeric(splited[[chrom]][, by2])), ])
    }
  }

  return(sorted)
}


convertRS <- function(rs, sampleStart = 4){
 
  splited <- split.data.frame(rs, factor(rs[, "chrom"]))
  splited <- lapply(splited, convertNAByChrom, sampleStart = sampleStart)

  return(do.call("rbind", args = splited))  
}


convertNAByChrom <- function(rsByChrom, sampleStart = 4, 
                             loc = 2){
  rsByChrom <- rsByChrom[order(as.numeric(rsByChrom[, loc])), ] 
  converted <- apply(rsByChrom[, sampleStart:ncol(rsByChrom)], 2, convertNA)
  
  return(cbind(rsByChrom[, 1:(sampleStart - 1)], converted))
}

convertNA <- function(convertMe){
 
  while(length(toConvert <- which(is.na(convertMe))) > 0){
    for(index in toConvert){
      if(index == 1){
        convertMe[index] <- convertMe[index + 1]
        next
      }
      if(index == length(convertMe)){
        convertMe[index] <- convertMe[index - 1]
        next
      }
      temp <- convertMe[c(index + 1, index - 1)]
      temp <- temp[!is.na(temp)]
      if(length(temp) == 0){
        next
      }else{
        convertMe[index] <- min(as.numeric(temp))
      }
    }
  }
  
  return(convertMe)
}


collapseSegList <- function(segList, geneMap){
  splited <- split.data.frame(segList, factor(segList[, "ID"]))
  template <- as.matrix(geneMap)
  rownames(template) <- gsub(" ", "", template[, "geneid"])
  cat("\nProcessing samples .")
  for(sp in names(splited)){
    cat(".")
    template <- cbind(template, rep(NA, nrow(template)))
    colnames(template) <- c(colnames(template)[-length(colnames(template))], 
                            sp)
    filled <- getSegMeanByGene(splited[[sp]], geneMap)
    template[as.character(gsub(" ", "", filled[, "geneid"])), sp] <- 
                            as.numeric(filled[, "seg.mean"])
  }
  cat(" Done\n")
  return(template[apply(template[, 7:ncol(template)], 1, 
                  FUN = function(x) !all(is.na(x))), ])
}


# Maps chromosomal regions to genes
# segList - output of CBS 
# geneInfo - a matrix or data frame with the following columns:
#            chrom, star, end, geneName, entrezID.
#
getSegMeanByGene <- function (sampleSeg, geneInfo) {

  mapMCA <- function(x){
    found <- geneInfo[as.character(gsub(" ", "", geneInfo[, "chrom"])) ==
                     gsub(" ", "", as.character(x["chrom"])) &
                     as.numeric(geneInfo[, "end"]) > 
                     as.numeric(x["loc.start"]) & 
                     as.numeric(geneInfo[, "start"]) < 
                     as.numeric(x["loc.end"]), , drop = FALSE]
    if(nrow(found) == 0){
      #return(matrix(c(id = NA, x[c("chrom", "loc.start", "loc.end")], 
      #    geneName = NA, entrezID = NA, x["seg.mean"]), nrow = 1, 
      #     dimnames = list(NULL, c("chrom", "start", "end", "geneName",
      #                     "entrezID", "seg.mean"))))
      return(NA)
      
    }
    
    return(cbind(found, seg.mean = as.numeric(as.vector(x["seg.mean"]))))
  }
  sampleSeg[, "loc.start"] <- gsub(" ", "", sampleSeg[, "loc.start"])
  sampleSeg[, "loc.end"] <- gsub(" ", "", sampleSeg[, "loc.end"])
  mapped <- do.call("rbind", args = apply(sampleSeg, 1, mapMCA))
  
  return(collapseDupGenes(mapped[which(!is.na(mapped[, "geneid"])), ]))
}


collapseDupGenes <- function(mapped){
  if(any(duplicated(mapped[, "geneid"]))){
      dups <- mapped[mapped[, "geneid"] %in% 
          mapped[which(duplicated(mapped[, "geneid"])), 
          "geneid"], c("geneid", "seg.mean")]
      dups <- split.data.frame(dups, factor(dups[, "geneid"])) 
  
      dups <- sapply(dups, FUN = function(x){
            x <- as.matrix(x)
            temp <- as.numeric(x[which(abs(as.numeric(x[, "seg.mean"])) 
                == max(abs(as.numeric(x[, "seg.mean"])), na.rm = TRUE)), "seg.mean"])
            return(temp[1])
          })
      mapped <- as.matrix(mapped[!duplicated(mapped[, "geneid"]), ])
      rownames(mapped) <- gsub(" ", "", mapped[, "geneid"])
      mapped[as.character(names(dups)), "seg.mean"] <- as.numeric(dups)
  }
  return(mapped)
}


getPairwise <- function(segList){
    rsList <- list()
    uniqueSamples <- unique(segList[, "ID"])
    cat("\nProcessing samples .")
    for(i in 1:(length(uniqueSamples) -1)){ 
        cat(".")
        for(j in (i + 1):length(uniqueSamples)){ 
            rs <- getCommonSegValues(segList[which(segList[, "ID"] 
            %in% c(uniqueSamples[i], uniqueSamples[j])), ], 
            drop = TRUE, verbose = FALSE)
            rs <- rs[!(is.na(rs[, uniqueSamples[i]]) | 
               is.na(rs[, uniqueSamples[j]])), ] 
            rsList[[paste(i, j, sep = "")]] <- rs[, -c(1:3)]
        }
    }
    cat(" done\n")
    return(rsList)
}


getPairDist <- function(pairList, method = "euclidean", diag = FALSE,
    upper = FALSE, p = 2){

    pair2Dist <- function(pair){
      tempDis <- dist(t(pair), 
        method = METHODS[method], diag = diag, upper = upper, p = 2)
      d[colnames(pair)[2], colnames(pair)[1]]  <<- tempDis
      d[colnames(pair)[1], colnames(pair)[2]] <<- tempDis
      return(invisible())
    }

    if (!is.na(pmatch(method, "euclidian"))) 
        method <- "euclidean"
    METHODS <- c("euclidean", "maximum", "manhattan", "canberra", 
        "binary", "minkowski")
    method <- pmatch(method, METHODS)
    if (is.na(method)) 
        stop("invalid distance method")
    if (method == -1) 
        stop("ambiguous distance method")
    sNames <- unique(unlist(lapply(pairList, 
        FUN = function(x) return(colnames(x)))))
    d <- matrix(data = 1, nrow = length(sNames), ncol = length(sNames), 
        dimnames = list(row = sNames, col = sNames))
    sapply(pairList, pair2Dist)

    return(as.dist(d))
}


getPairCor <- function(pairList, use = "everything", 
    method = c("pearson", "kendall", "spearman")){

    pair2Cor <- function(pair){
      tempCor <- cor(pair[, 1], pair[, 2], use = use, method = method)
      co[colnames(pair)[2], colnames(pair)[1]]  <<- tempCor
      co[colnames(pair)[1], colnames(pair)[2]] <<- tempCor
      return(invisible())
    }

    sNames <- unique(unlist(lapply(pairList, 
        FUN = function(x) return(colnames(x)))))
    co <- matrix(data = 1, nrow = length(sNames), ncol = length(sNames), 
        dimnames = list(row = sNames, col = sNames))
    junk <- sapply(pairList, pair2Cor)

    return(co)
}


filterPair <- function(rsObj, flist){
  tempPair <- lapply(rs(rsObj), FUN = function(x) 
      x[genefilter:::genefilter(x, flist), ])
  return(RS(rs = tempPair, by = segBy(rsObj)))  
}

filterRS <- function(rsObj, flist){
    if(segBy(rsObj) == "region"){
        drop <- 1:3
    }else{
        drop <- 1:5
    }
    tempRS <- rs(rsObj)

    return(RS(rs = rs(rsObj)[genefilter(tempRS[, -drop], flist), ], 
       by = segBy(rsObj)))
}

filterByMad <- function(rsObj, cutoff = 0.80){
  if(segBy(rsObj) == "pair"){
    stop("The function is not applicable to paired RS")
  }
  if(segBy(rsObj) == "region"){
        drop <- 1:3
  }else{
        drop <- 1:5
  }
  tempRS <- rs(rsObj)
  mads <- apply(apply(tempRS[, -drop], 2, as.numeric), 1, mad, na.rm = TRUE)
  return(RS(rs = tempRS[mads >= quantile(mads, cutoff), ], 
               by = segBy(rsObj)))
}


diffBy <- function(marging = 0.1){
    function(x){
        abs(as.numeric(x[1]) - as.numeric(x[2])) >= marging
    }
}







procRefGene <- function(refGene, refLink){
  # Process UCSC gene file
  refgene <- read.delim(refGene, sep = "\t", header = FALSE, 
    as.is = TRUE)[, c(2, 3, 5, 6, 13)]
  colnames(refgene) <- c("refSeq", "chrom", "start", "end", "geneName")
  refgene[, "chrom"] <- gsub("chr", "", refgene[, "chrom"])
  # drop random and hap entries
  refgene <- refgene[-grep("_random|_hap", refgene[, "chrom"]), ]
  # Process the link file
  reflink <- read.delim(refLink, sep = "\t", header = FALSE, 
    as.is = TRUE)[, c(3, 7)]
  colnames(reflink) <- c("refSeq", "entrezID") 
  # Merge the two by refSeq ids
  merged <- merge(refgene, reflink, by = "refSeq", all.x = TRUE)[, -1]
  
  return(sortByChromNLoc(merged, by1 = "chrom", by2 = "start"))
}


# Sort mapped in order by chromsome and then by location
sortByChromNLoc <- function(sortMe, by1 = "Ch", by2 = "Pos"){
  splited <- split.data.frame(sortMe, factor(sortMe[, by1]))
  sorted <- NULL
  for(chrom in c(1:22, "X", "Y")){
    if(!is.null(splited[[chrom]])){
      sorted <- rbind(sorted,
                      splited[[chrom]][
                         order(as.numeric(splited[[chrom]][, by2])), ])
    }
  }

  return(sorted)
}



# rois - output by function findROI
plotSegByChrom <- function(rois){
  drawPolygon <- function(x, isGain){
    if(all(x == "") | all(x == 0)){
      return(invisible())
    }
    if(isGain){
      color <- "red"
    }else{
      color = "green"
    }
    if(as.numeric(x["seg.mean"]) != 0){
      poly <- region2Polygon(x)
      polygon(as.numeric(poly[, "x"]), 
              as.numeric(poly[, "y"]), col = color, border = color)
    }
   
    return(invisible())
  }
  
  chromLength <- getChromLength(rois)
  par(mfrow = c(ceiling(length(chromLength)/2), 2),
      mar = c(0.1, 4, 0.1, 1), oma = c(0.1, 0.1, 0.1, 0.5))
 
  for(chrom in names(chromLength)){
    temp <- segList[gsub(" ", "", as.character(segList[, "chrom"])) == gsub(" ", "", as.character(chrom)), 
             , drop = FALSE]
    plot.new()
    plot.window( xlim = c(0, chromLength[chrom]),
         ylim = range(as.numeric(as.vector(temp[, "seg.mean"])), 
                na.rm = TRUE))
    apply(temp[as.numeric(as.vector(temp[, "seg.mean"])) >= 0, , 
                drop = FALSE], 1, drawPolygon, isGain = TRUE)
    apply(temp[as.numeric(as.vector(temp[, "seg.mean"])) < 0, , 
                drop = FALSE], 1, drawPolygon, isGain = FALSE)
    axis(2, lwd = 1.5, tck = -0.05, labels = FALSE)
    mtext(chrom, 4, cex = 0.8, adj = 1)
    abline(h = 0)
  }

  return(invisible())
}

getChromLength <- function(segList, by = "chrom"){
  splited <- split.data.frame(segList, factor(segList[, by]))
  chroms <- c(1:22, "X", "Y")
  temp <- sapply(splited, FUN = function(x) max(as.numeric(as.vector(x[,
                           "loc.end"]))))

  return(temp[chroms[which(chroms %in% names(temp))]])
}

# Convert start and end location data to a matrix so that polygons can be
# drawn
region2Polygon <- function(startNend, isGain = TRUE){
  return(cbind(x = c(startNend["loc.start"], startNend["loc.start"],
                   startNend["loc.end"], startNend["loc.end"]),
                 y = c(0, startNend["seg.mean"], startNend["seg.mean"], 0)))
 
}




