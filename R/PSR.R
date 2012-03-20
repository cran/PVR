PSR <-
function(x, trait = NULL, null.model = FALSE, hypothesis = "two", neutral.model = FALSE, times = 1000){
####---------Computing psr curve
		
		pvr <- x@Eigen
		pD <- x@phyDist
		psr <- data.frame(r.squared = 0, eigenvalues = 0)

		#---------loading species traits
		if(is.character(trait)){
			
			trait <- read.table(trait)
		}
		
		trait <- as.matrix(trait)
		Ntraits <- ncol(trait)
		areas <- data.frame(Id = 0, PSR.area = 0, p = NA)
		nullareas <- data.frame(Id = 0, mean.area = 0, variance.area = NA)
		Naxis <- ncol(pvr$vectors)
		for(k in 1:Ntraits){
			
			pvr$vectors[ ,1:Naxis] <- t(t(pvr$vectors[ ,1:Naxis])/sqrt(pvr$values[1:Naxis]))
			relVal <- pvr$values/sum(pvr$values)
			for(i in 1:Naxis){
				
				psr[i, 1] <- summary(lm(trait[,k] ~ matrix(as.real(pvr$vectors[,1:i]), ncol = i)))$r.squared
				psr[i, 2] <- sum(relVal[1:i])			
			}
			
			#---------Computing psrarea
			coords = data.frame(x = c(0, psr[ ,2]), y = c(0, psr[ ,1]))
			while(round(coords[nrow(coords), 1],7) == round(1,7) & round(coords[nrow(coords), 2],7) == round(1,7)){
				
				coords <- coords[-nrow(coords),]
			}
			
			coords[(nrow(coords) + 1), ] <- 1
			psrarea <- .Area(as.matrix(coords))
			
			#---------Creating PSRarea null distribution
			if(null.model){
				
				nullDistribution <- matrix(0, nrow = Naxis, ncol = times)
				psrNull <- data.frame(r.squared = 0, eigenvalues = 0)
				nullPsrarea <- data.frame(a = 0)
				for(t in 1:times){
					
					#---------shuffling trait vector
					traitRand <- sample(trait[,k])
					
					#---------Computing psr random curve
					for(i in 1:Naxis){
						
						psrNull[i, 1] <- summary(lm(traitRand ~ matrix(as.real(pvr$vectors[,1:i]), ncol = i)))$r.squared
						psrNull[i, 2] <- sum(relVal[1:i])		
					}
					
					nullDistribution[ ,t] <- psrNull[ ,1]
					#---------Computing psrarea
					coordsRand = data.frame(x = c(0, psrNull[ ,2]), y = c(0, psrNull[ ,1]))
					while(round(coordsRand[nrow(coordsRand), 1],7) == round(1,7) & round(coordsRand[nrow(coordsRand), 2],7) == round(1,7)){
						
						coordsRand <- coordsRand[-nrow(coordsRand),]
					}
					coordsRand[(nrow(coordsRand) + 1), ] <- 1
					
					nullPsrarea[t,1] <- .Area(as.matrix(coordsRand))
				}
				
				#---------Computing p
				distribCum <- ecdf(nullPsrarea[,1])
				p2 <- 1 - distribCum(psrarea) 
			}
			areas[k, 1] <- paste("traits set ", k, sep = "")
			areas[k, 2] <- psrarea
			if(null.model){
				
				areas[k, 3] <- p2
				nullareas[k, 1] <- paste("traits.set.", k, sep = "")
				nullareas[k, 2] <- mean(nullPsrarea[ ,1])
				nullareas[k, 3] <- var(nullPsrarea[ ,1])
			}
			 
		}
			#---------Creating PSRarea neutral distribution
			if(neutral.model){
				
				neutralDistribution <- matrix(0, nrow = Naxis, ncol = times)
				psrNeutral <- data.frame(r.squared = 0, eigenvalues = 0)
				neutralPsrarea <- data.frame(a = 0)
				phy <- x@phylo
				for(t in 1:times){
					
					#---------shuffling trait vector
					traitRand <- rTraitCont(phy, mpdel = "BM")
					
					#---------Computing psr random curve
					for(i in 1:Naxis){
						
						psrNeutral[i, 1] <- summary(lm(traitRand ~ matrix(as.real(pvr$vectors[,1:i]), ncol = i)))$r.squared
						psrNeutral[i, 2] <- sum(relVal[1:i])		
					}
					
					neutralDistribution[ ,t] <- psrNeutral[ ,1]
					
					#---------Computing psrarea
					coordsNeutral = data.frame(x = c(0, psrNeutral[ ,2]), y = c(0, psrNeutral[ ,1]))
					while(round(coordsNeutral[nrow(coordsNeutral), 1],7) == round(1,7) & round(coordsNeutral[nrow(coordsNeutral), 2],7) == round(1,7)){
						
						coordsNeutral <- coordsNeutral[-nrow(coordsNeutral),]
					}
					
					coordsNeutral[(nrow(coordsNeutral) + 1), ] <- 1
					neutralPsrarea[t,1] <- .Area(as.matrix(coordsNeutral))
				}
			}
		results <- new("PSR", x)
		results@PSRarea <- areas
		results@PSR <- data.frame(Cumul.eigen.values = psr[,2], r.squared = psr[,1])

		Expect.area.values <- list(
				Null.model = data.frame(),
				meanNeutral = as.numeric(NA),
				varianceNeutral = as.numeric(NA)
				)
		
		if(null.model){
			
			Expect.area.values$Null.model <- nullareas 
			results@nullPSR <- nullDistribution
			
		}
		
		if(neutral.model){
			
			Expect.area.values$meanNeutral <- mean(neutralPsrarea[ ,1]) 
			Expect.area.values$varianceNeutral <- var(neutralPsrarea[ ,1])
			results@neutralPSR <- neutralDistribution
		}
		
		results@Expect.area.values <- Expect.area.values
		return(results)
}
