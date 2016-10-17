CNSeg <- function(segList, id = "ID", chromosome = "chrom", start = "loc.start", 
    end = "loc.end", segMean = "seg.mean"){
    return(new("CNSeg", segList = segList, chromosome = chromosome, id = id,
        start = start, end = end, segMean = segMean))
}

RS <- function(rs, by){
    return(new("RS", rs = rs, by = by))
}
    
seg2RS <- function(segData, by = c("region", "gene", "pair"), 
    imput = TRUE, XY = FALSE, geneMap, what = c("mean", "median", "max", "min"), 
    mapChrom = "chrom", mapStart = "start", mapEnd = "end"){
  by <- match.arg(by)
  what <- match.arg(what)
  if(by == "gene" & missing(geneMap)){
      stop("Need geneMap when by == \"gene\"")
  }
  if(!XY){
    seg <- segList(segData)
    segList(segData) <- seg[which(!seg[, chromosome(segData)] 
        %in% c("X", "Y", "x", "y")),]
  }
  rs <- switch(by,
	       region = getCommonSegValues(segList(segData), 
                   drop = TRUE, what = what, segID = id(segData), segChrom = chromosome(segData), 
                   segStart = start(segData), segEnd = end(segData),  
                   segMean = segMean(segData)),
               gene = getReducedSeg(segList(segData), geneMap, what = what, 
                   segID = id(segData), segChrom = chromosome(segData), 
                   segStart = start(segData), segEnd = end(segData),  
                   segMean = segMean(segData), mapChrom = mapChrom, mapStart = mapStart, 
                   mapEnd = mapEnd),
               pair = getPairwise(segData, imput = imput, XY = XY, what = what) 
  )
  
  if(imput & by != "pair"){
     rs <- convertRS(rs, sampleStart = ifelse(by == "region", 4, 6))
  }
  if(by != "pair"){
      return(RS(as.data.frame(rs), by))
  }
  return(RS(rs, by))
}

## segList - "output" of CBS
## drop - drop reduced segments with all NAs
## by - column name for chromosome
getCommonSegValues <- function(seglist,  
    drop = FALSE, what = c("mean", "median", "max", "min"), segID = "ID",  segChrom = "chrom",
    segStart = "start", segEnd = "end", segMean = "seg.mean", verbose = TRUE){
 what <- match.arg(what)
 
  if(verbose){
  	cat("Processing samples ...")
  }
  segsByChroms <- do.call("rbind", args = 
      lapply(split.data.frame(seglist, factor(seglist[, segChrom])), findOverlapingSegs, 
          segChrom = segChrom, segStart = segStart, segEnd = segEnd))
  segWithValues <- getReducedSeg(seglist, segsByChroms, what = what, segID = segID, 
      segChrom = segChrom, segStart = segStart, segEnd = segEnd, segMean = segMean, 
      mapChrom = "chrom", mapStart = "start", mapEnd = "end")
  
  if(drop){
    toKeep <- apply(segWithValues[, 4:ncol(segWithValues)], 1,
                  FUN = function(x) !all(is.na(x)))
    segWithValues <- segWithValues[toKeep, ]
  }
 
  sorted <- sortByChromNLoc(segWithValues, by1 = "chrom", by2 = "start")
  if(verbose){
  	cat(" Done\n")
  }
  
  return(sorted)
}

findOverlapingSegs <- function(segListByChrom, segChrom = "chrom", segStart = "loc.start", 
    segEnd = "Loc.end"){
  breakPoints <- sort(unique(c(as.numeric(
                     as.vector(segListByChrom[, segStart])),
                     as.numeric(as.vector(segListByChrom[, segEnd])))))
  newSegs <- cbind(c(breakPoints[1], breakPoints[2:(length(breakPoints) - 1)]
                     + 1), breakPoints[2:length(breakPoints)])
  newSegs <- cbind(as.vector(segListByChrom[1, segChrom]), newSegs)
  colnames(newSegs) <- c("chrom", "start", "end")
  
  return(newSegs)
}

sortByChromNLoc <- function(sortMe, by1 = "Ch", by2 = "Pos"){
  splited <- split.data.frame(sortMe, factor(sortMe[, by1]))
  sorted <- NULL
  for(chrom in c(1:25, "X", "x", "Y", "y")){
    if(!is.null(splited[[chrom]])){
      sorted <- rbind(sorted,
              splited[[chrom]][order(as.numeric(splited[[chrom]][, by2])), ])
    }
  }

  return(sorted)
}



# seglist - segment list with at least 4 columns for chromosome, start, end, and segment mean data
# map - matrix with at least 3 columns for chromosome, start, and end data
# what - the method used to calculate segment value in case there are multiple segment values
#             for a given region defined in the map
getReducedSeg <- function(seglist, map, what = c("mean", "median", "max", "min"), segID = "ID", 
    segChrom = "chrom", segStart = "loc.start", segEnd = "loc.end",  segMean = "seg.mean", 
    mapChrom = "chrom",  mapStart = "start", mapEnd = "end"){
    
    what <- match.arg(what)
    getGeneSegMean <- function(segData){
        segged <- rep(0, nrow(map))
        return( .C("getratios", as.character(map[, mapChrom]), as.double(map[, mapStart]), 
            as.double(map[, mapEnd]), as.integer(nrow(map)), as.character(segData[, segChrom]),
            as.double(segData[, segStart]), as.double(segData[, segEnd]), as.integer(nrow(segData)),
            as.double(segData[, segMean]), as.character(what), as.double(segged), 
            PACKAGE = "CNTools")[[11]])      
    }
    
    splited <- split.data.frame(seglist, factor(seglist[, segID]))
    
    return(cbind(map, do.call("cbind", args = lapply(splited, getGeneSegMean)))) 
}


getPairwise <- function(segData, imput, XY, what){    
    segs <- segList(segData)
    rsList <- list()
    uniqueSamples <- unique(segs[, id(segData)])
    cat("\nProcessing samples .")
    for(i in 1:(length(uniqueSamples) -1)){ 
        cat(".")
        for(j in (i + 1):length(uniqueSamples)){ 
            rs <- CNTools:::getCommonSegValues(segs[which(segs[, id(segData)] 
                %in% c(uniqueSamples[i], uniqueSamples[j])), ],  drop = TRUE, 
                verbose = FALSE, what = what, segID = id(segData),
                segChrom = chromosome(segData), segStart = start(segData), segEnd = end(segData), 
                segMean = segMean(segData))
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
  convertMe[which(convertMe == 0)] <- NA
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

