netplot <- function(x,text=TRUE,col="black",bg="blue",base.lwd=1,base.cex=1){			# ネットワークプロットの作図

	xms <- x$measure

	if(xms=="OR"||xms=="RR"||xms=="RD"||xms=="HR"||xms=="SPD"){

	study <- x$study
	trt <- x$trt
	d <- x$d
	n <- x$n
	
	###
	
	A1 <- factor(trt)
	A2 <- levels(A1)

	l <- length(trt)
	m <- length(A2)
	
	z <- rep(NA,times=l)
	
	treat <- as.numeric(l)

	for(i in 1:m)	treat[A1==A2[i]] <- i

	N <- max(study)
	p <- max(treat)
	
	t1 <- seq(0, 2*pi, length=p+1)
	x1 <- cos(t1)
	y1 <- sin(t1)

	plot(x1, y1, type="n", axes=FALSE, xlab="", ylab="", xlim=c(-1.3,1.3), ylim=c(-1.3,1.3))

	###

	cmb <- combn(1:p,m=2)
	L <- dim(cmb)[2]
	lwd2 <- numeric(L)

	for(i in 1:L){
	
		k <- cmb[,i]
	
		for(j in unique(study)){
		
			wj <- which(study==j)
			ij <- intersect(k,treat[wj])
			
			if(length(ij)>=2){
			
				u1 <- which(treat[wj]==k[1])
				u2 <- which(treat[wj]==k[2])
				lwd2[i] <- lwd2[i] + n[wj][u1] + n[wj][u2]

			}

		}

	}
				
	sp3 <- approxfun(c(min(lwd2[lwd2>0]),max(lwd2)),c(1,10),method = "linear")

	lwd3 <- numeric(L)
	for(i in 1:L)	lwd3[i] <- base.lwd*sp3(lwd2[i])
		
	for(i in 1:L){
	
		if(is.na(lwd3[i])==FALSE){
		
			k <- cmb[,i]
			lines(x1[k],y1[k],lwd=lwd3[i])
		
		}
		
	}
	
	###
	
	np <- rep(NA,times=p)
	
	for(i in 1:p){

		wi <- which(treat==i)
		np[i] <- sum(n[wi])
	
	}
	
	sp2 <- approxfun(c(min(np),max(np)),c(3,8),method = "linear")

	cex2 <- numeric(p)
	for(i in 1:p)	cex2[i] <- base.cex*sp2(np[i])
	
	for(i in 1:p)	points(x1[i], y1[i], pch=21,col=col,bg=bg,cex=cex2[i],lwd=2)
	
	if(text==TRUE)	for(i in 1:p)	text(1.25*x1[i], 1.25*y1[i], A2[i])
	
	}

	if(xms=="MD"||xms=="SMD"){

	study <- x$study
	trt <- x$trt
	n <- x$n
	
	###
	
	A1 <- factor(trt)
	A2 <- levels(A1)

	l <- length(trt)
	m <- length(A2)
	
	z <- rep(NA,times=l)
	
	treat <- as.numeric(l)

	for(i in 1:m)	treat[A1==A2[i]] <- i

	N <- max(study)
	p <- max(treat)
	
	t1 <- seq(0, 2*pi, length=p+1)
	x1 <- cos(t1)
	y1 <- sin(t1)

	plot(x1, y1, type="n", axes=FALSE, xlab="", ylab="", xlim=c(-1.3,1.3), ylim=c(-1.3,1.3))

	###

	cmb <- combn(1:p,m=2)
	L <- dim(cmb)[2]
	lwd2 <- numeric(L)

	for(i in 1:L){
	
		k <- cmb[,i]
	
		for(j in unique(study)){
		
			wj <- which(study==j)
			ij <- intersect(k,treat[wj])
			
			if(length(ij)>=2){
			
				u1 <- which(treat[wj]==k[1])
				u2 <- which(treat[wj]==k[2])
				lwd2[i] <- lwd2[i] + n[wj][u1] + n[wj][u2]

			}

		}

	}
				
	sp3 <- approxfun(c(min(lwd2[lwd2>0]),max(lwd2)),c(1,10),method = "linear")

	lwd3 <- numeric(L)
	for(i in 1:L)	lwd3[i] <- base.lwd*sp3(lwd2[i])
		
	for(i in 1:L){
	
		if(is.na(lwd3[i])==FALSE){
		
			k <- cmb[,i]
			lines(x1[k],y1[k],lwd=lwd3[i])
		
		}
		
	}
	
	###
	
	np <- rep(NA,times=p)
	
	for(i in 1:p){

		wi <- which(treat==i)
		np[i] <- sum(n[wi])
	
	}
	
	sp2 <- approxfun(c(min(np),max(np)),c(3,8),method = "linear")

	cex2 <- numeric(p)
	for(i in 1:p)	cex2[i] <- base.cex*sp2(np[i])
	
	for(i in 1:p)	points(x1[i], y1[i], pch=21,col=col,bg=bg,cex=cex2[i],lwd=2)
	
	if(text==TRUE)	for(i in 1:p)	text(1.25*x1[i], 1.25*y1[i], A2[i])
	
	}
    
}
