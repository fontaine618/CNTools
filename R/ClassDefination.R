# A class that contains the segmentation
# 
# Copyright 2009, Jianhua Zhang, all rights reserved
#

setClass("CNSeg", representation(segList = "data.frame", 
                                 id = "character",
                                 chromosome = "character",
                                 start = "character",
                                 end = "character",
                                 segMean = "character"))

setClass("RS", representation(rs = "ANY", 
                              by = "character"))

setGeneric("getRS", function(object, by = c("region", "gene", "pair"), imput = TRUE, XY = FALSE, 
    geneMap, what = c("mean", "max", "mini", "median"),  mapChrom = "chrom", 
    mapStart = "start", mapEnd = "end")
           standardGeneric("getRS"))
setMethod("getRS", "CNSeg", 
          function(object, by = c("region", "gene", "pair"), 
            imput = TRUE, XY = FALSE, geneMap, what = c("mean", "max", "mini", "median"), 
            mapChrom = "chrom",  mapStart = "start", mapEnd = "end")
              seg2RS(object, by, imput, XY, geneMap, what = what, mapChrom = mapChrom, 
              mapStart = mapStart, mapEnd = mapEnd))

setGeneric("segList", function(object)
           standardGeneric("segList"))
setMethod("segList", "CNSeg",
          function(object) object@segList)
          
setGeneric("chromosome", function(object)
           standardGeneric("chromosome"))
setMethod("chromosome", "CNSeg",
          function(object) object@chromosome) 
          
setGeneric("start", function(object)
	             standardGeneric("start"))
	  setMethod("start", "CNSeg",
          function(object) object@start)
          
setGeneric("end", function(object)
           standardGeneric("end"))
setMethod("end", "CNSeg",
          function(object) object@end)

setGeneric("segMean", function(object)
           standardGeneric("segMean"))
setMethod("segMean", "CNSeg",
          function(object) object@segMean)
          
setGeneric("id", function(object)
           standardGeneric("id"))
setMethod("id", "CNSeg",
          function(object) object@id)
 
setMethod("show", "CNSeg",
          function(object) {
            cat("Object of CNSeg\n")
            cat(paste("Number of samples:", 
                 length(unique(segList(object)[, id(object)])), 
                 "\n"), sep = " ")
            cat(paste("\tRow = ", nrow(segList(object)), 
                "\n", sep = ""))
            cat(paste("\tColumn = ", ncol(segList(object)), 
                "\n", sep = ""))
            if(nrow(segList(object)) <= 5){
              rows <- nrow(segList(object))
            }else{
              rows <- 5
            }
            print(segList(object)[1:rows, ])
            if(nrow(segList(object) > 5)){
              cat(".........................\n")
            }
          })

setGeneric("segList<-", function(object, value)
           standardGeneric("segList<-"))
setReplaceMethod("segList", "CNSeg", function(object, value){
  object@segList <- value; object})

setGeneric("rs", function(object)
           standardGeneric("rs"))
setMethod("rs", "RS", 
          function(object) object@rs)

setGeneric("rs<-", function(object, value)
           standardGeneric("rs<-"))
setReplaceMethod("rs", "RS", 
          function(object, value){ object@rs <- value; object})

setGeneric("segBy", function(object)
           standardGeneric("segBy"))
setMethod("segBy", "RS", 
          function(object) object@by)

setGeneric("madFilter", function(object, cutoff = 0.8)
           standardGeneric("madFilter"))
setMethod("madFilter", "RS", 
          function(object, cutoff = 0.8) filterByMad(object, cutoff))

#setGeneric("genefilter", function(expr, flist)
#           standardGeneric("genefilter"))
setMethod("genefilter", "RS",
          function(expr, flist){
            if(segBy(expr) == "pair"){
              return(filterPair(expr, flist))
            }else{
              return(filterRS(expr, flist))
            }
          })

setGeneric("getDist", function(x, method = "euclidean", diag = FALSE, 
           upper = FALSE, p = 2)
           standardGeneric("getDist")) 
setMethod("getDist", "RS",
          function(x, method = "euclidean", diag = FALSE, 
              upper = FALSE, p = 2){
            if(segBy(x) == "pair"){
              return(getPairDist(rs(x)))
            }else{
              if(segBy(x) == "region"){
                drop <- 1:3
              }else{
                drop <- 1:5
              }
              return(dist(t(rs(x)[, -drop]), method = method,
                  diag = diag, upper = upper, p = p))
            }})

setGeneric("getCor", function(x, y = NULL, use = "everything", 
           method = c("pearson", "kendall", "spearman"))
           standardGeneric("getCor"))
setMethod("getCor", "RS", 
          function(x, y = NULL, use = "everything", 
              method = c("pearson", "kendall", "spearman")){
              if(segBy(x) == "pair"){
                  return(getPairCor(rs(x), use = use, method = method))
              }else{
                  if(segBy(x) == "region"){
                      drop <- 1:3
                  }else{
                      drop <- 1:5
                  }
                  return(cor(rs(x)[, -drop], use = use, method = method))
              }
          })

setMethod("show", "RS",
          function(object) {
            cat("Object of RS\n")
            if(segBy(object) == "pair"){
              cat(paste("Number of pairs:", length(rs(object))))
            }else{
              if(segBy(object) == "region"){
                  feature <- "regions"
              }else{
                  feature <- "genes"
              }	       
              cat(paste("\tNumber of ", feature, " = ", nrow(rs(object)), 
                  "\n", sep = ""))
              cat(paste("\tNumber of samples = ", ncol(rs(object)) - 
                  ifelse(segBy(object) == "region", 3, 5), "\n", sep = ""))
            }
            cat(paste("\nSegment data converted based on", segBy(object)))
            cat("\n")
          })

