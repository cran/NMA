nmafunnel <- function(x, method="NH", legends="topright"){

	nmx <- nma(x,method=method)
	coefx <- nmx[[5]][,1]

	study <- x$study
	treat <- x$treat
	n <- x$n

	y <- x$y
	S <- x$S
	
	study <- as.numeric(factor(study))

    ####

	treat1 <- sort(unique(treat))
	
	N <- length(unique(study))
	p <- max(treat) - 1

	des <- n.arm <- nm <- numeric(N)
	Ti <- NULL

	for(i in 1:N){

		wi <- which(study==i)
		ti <- sort(treat[wi],decreasing=FALSE)
		Ti[[i]] <- ti

		di <- NULL
		for(j in 1:length(wi)){
			if(is.null(di)==FALSE) di <- paste0(di,"-",ti[j])
			if(is.null(di)) di <- paste0(di,ti[j])
		}

		des[i] <- di
		n.arm[i] <- length(wi)
		nm[i] <- sum(n[wi])
		
	}
	
	des0 <- sort(unique(des))
	
	w1 <- which(substr(des,1,2)=="1-")
	w0 <- which(substr(des0,1,2)=="1-")

	###

	L <- length(w0)
	n.des <- n1 <- rep(NA,times=L)
	
	for(i in w0){
	
		desi <- des0[i]
		wi <- which(des==desi)

		n.des[i] <- length(wi)
		n1[i] <- sum(nm[wi])

	}

	R1 <- data.frame(des0[w0],n.des,n1)
	colnames(R1) <- c("design","N","n")

	###

	y1 <- y[w1,]
	S1 <- S[w1,]

	R3 <- NULL

    for(i in 1:length(w1)){
      
      yi <- as.vector(y1[i,])
      wi <- which(is.na(yi)==FALSE)
      pl <- length(wi)
      
      Si <- vmat(S1[i,], p)
      
      yi <- yi[wi]
      Si <- pmat(Si,wi)

	  for(j in 1:length(wi)){

		ayi <- yi[j] - coefx[wi[j]]
		R3 <- rbind(R3,c(wi[j],yi[j],ayi,coefx[wi[j]],diag(Si)[j]))
	  
	  }

    }
    
	y2 <- R3[,3]
	s2 <- sqrt(R3[,5])
	
	plot(s2~y2, type="n", xlim=c(-max(abs(y2)),max(abs(y2))), ylim=c(max(s2),0), xlab="effect size centered at comparison-specific pooled effect", ylab="standard error of effect size", main="Comparison adjusted funnel plot")
	abline(v=0)
	lines(c(0,1.96*10^4),c(0,10^4),lty=2)
	lines(c(0,-1.96*10^4),c(0,10^4),lty=2)
	
	su3 <- sort(unique(R3[,1]))
	
	for(i in su3){
	
		wi <- which(R3[,1]==i)
		points(y2[wi],s2[wi],pch=19,col=i,cex=1)

	}
	
	su4 <- paste0("1 vs. ",(su3+1))
	
	legend(legends, su4, pch = 19, col=su3, bg = "transparent")

	message("Comparison adjusted funnel plot for the trials involving treatment 1 (as control)")

	R5 <- list("coding"=x$coding,"summary"=R1)

	return(R5)

}


