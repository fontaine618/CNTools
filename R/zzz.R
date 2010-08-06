.onLoad <- function(lib, pkgname){
   require("methods", quietly = TRUE) || stop("methods package not found")
   require("tools", quietly = TRUE) || stop("methods package not found")
   library.dynam( "CNTools", package = "CNTools") 
}

.onUnload <- function( libpath ) {
   library.dynam.unload( "CNTools", libpath )
}

 