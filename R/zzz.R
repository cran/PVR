.onLoad <-function(libname, pkgname)
{
	ver <- read.dcf(file.path(libname, pkgname, "DESCRIPTION"), "Version")
	ver <- as.character(ver)	
	
	packageStartupMessage("PVR ", ver, " loaded\n")
}

setOldClass("phylo")

setClass("PVR", representation(
				Eigen = "list", phyDist = "matrix", phylo = "phylo", Selection = "list", PVR = "list", VarPart = "list")
)

setClass("PSR", representation(
				PSRarea = "data.frame", PSR = "data.frame", Expect.area.values = "list", nullPSR = "matrix", neutralPSR = "matrix"),
			contains = "PVR" 
)

setMethod("plot",
		signature(x = "PSR"),
		function (x, y, expected, ...) 
		{
			plot(slot(x, "PSR"), 
					xlab = "cumulative eigenvalues(%)", 
					ylab = "r squared", ylim = c(0,1), 
					xlim = c(0,1), bty = "L", expected = segments(0, 0, 1, 1), pch = 16, ...)
			
#			coords = data.frame(x = c(0,slot(x, "PSR")[,1]), y = c(0,slot(x, "PSR")[,2]))
#			polygon(coords)
		}
)

setMethod("show",
		signature(object = "PVR"),
		function (object) 
		{
			if(!is.null(object@PVR$R2)){
				
				cat("\n", "Selected vectors", "\n", sep = "")
				cat("Used selection method: ", object@Selection$Method, "\n", sep = "")
				cat("R2: ",  object@PVR$R2, "\n", sep = "")
				cat("Residuals: ",  "\n", object@PVR$Residuals, "\n", sep = "")	
			} else{
				
				cat("\n", "Unselected vectors", "\n", sep = "")
				print(object@Eigen)
			}					
		}
)

setMethod("show",
		signature(object = "PSR"),
		function (object) 
		{
			print(object@PSRarea)
		}
)
