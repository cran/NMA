transitivity <- function(x, z, gcol="blue", yrange=NULL,digits=3){

	call <- match.call()

	covariate <- deparse(substitute(z))

	z <- x$Z[, deparse(substitute(z))]

	study <- x$study
	treat <- x$treat
	n <- x$n
	
	study <- as.numeric(factor(study))
	
    ####

	treat1 <- sort(unique(treat))
	
	N <- length(unique(study))
	p <- max(treat)
	q <- p - 1

	des <- n.arm <- zm <- nm <- numeric(N)
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
		zm[i] <- sum(z[wi]*n[wi])/sum(n[wi])
		nm[i] <- sum(n[wi])
	
	}
	
	des0 <- sort(unique(des))

	L <- length(des0)
	n.des <- z1 <- z2 <- z3 <- n1 <- rep(NA,times=L)
	
	for(i in 1:L){
	
		desi <- des0[i]
		wi <- which(des==desi)

		n.des[i] <- length(wi)
		z1[i] <- sum(zm[wi]*nm[wi])/sum(nm[wi])
		z2[i] <- min(zm[wi])
		z3[i] <- max(zm[wi])
		n1[i] <- sum(nm[wi])

	}

	z1 <- round(z1,digits)
	z2 <- round(z2,digits)
	z3 <- round(z3,digits)

	R1 <- data.frame(des0,n.des,n1,z1,z2,z3)
	colnames(R1) <- c("design","N","n","wt.mean","min","max")

	o4 <- order(R1[,4])
	R1 <- R1[o4,]

	if(is.null(yrange)) plot(1:L,R1$wt.mean,xaxt="n",xlab="design",ylab="study-level covariate",col=gcol,pch=15,cex=2,ylim=c(min(zm),max(zm)))
	if(is.null(yrange)==FALSE) plot(1:L,R1$wt.mean,xaxt="n",xlab="design",ylab="study-level covariate",col=gcol,pch=15,cex=2,ylim=yrange)
	axis(1,at=1:L,formatC(R1[,1]))
	
	for(k in 1:L){

		i <- o4[k]
		desi <- des0[i]
		wi <- which(des==desi)

		for(j in wi) points(k,zm[j],col="gray",pch=20)
	
	}

	#spf <- splinefun(1:L,R1$wt.mean)
	#x2 <- seq(1, L, by = 0.01)
	#y2 <- spf(x2)

	#lines(x2,y2,col=gcol,lty=2)

	legend("bottomright", c("Weighted Mean","Observations"), pch = c(15,20), col=c(gcol,"gray"), bg = "transparent")

	R4 <- list("coding"=x$coding,"covariate"=covariate,"summary"=R1,digits=digits,call=call)

	class(R4) <- "transitivity"
	return(R4)
	
}
	

print.transitivity <- function(x, digits = x$digits, ...) {

  cat("Call:\n")
  print(x$call,row.names=FALSE)
  cat("\n")
  
  cat("Coding:\n", sep = "")
  print(x$coding,row.names=FALSE)
  cat("\n")

  cat("covariate: ", sep = "")
  cat(x$covariate)
  cat("\n")
  cat("\n")

  cat("Summary: ", sep = "")
  cat("\n")
  A <- x[[3]]
  print(A,row.names=FALSE)
  cat("\n")
  
  invisible(x)
  
}
