.onLoad <- function(lib, pkgname){
   require("methods", quietly = TRUE) || stop("methods package not found")
} 