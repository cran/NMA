nmaforest <- function(x,method="NH",col.plot="black",digits=3,ascending=TRUE){

	xms <- x$measure

	if((xms=="RD")||(xms=="MD")||(xms=="SMD")){
	
	nmx <- nma(x,method=method)
	
	treat <- x$coding[,2]

	x.logform <- FALSE

	A1 <- nmx[[5]]

	o1 <- order(A1[,1])	
	if(ascending==FALSE) o1 <- order(A1[,1],decreasing = TRUE)
	treat <- treat[-1]; treat <- treat[o1]
	
	A1 <- A1[o1,]
	
	p <- dim(A1)[1]
	
	A2 <- data.frame(coef=A1[,1],low=A1[,3],high=A1[,4]); A2 <- rbind(c(NA,NA,NA),A2)
	
	A2 <- B2 <- round(A2,digits)
	
	for(i in 2:dim(A2)[1]){
		for(j in 1:dim(A2)[2]){
		
			B2[i,j] <- rdc(A2[i,j],digits)

		}
	}
	
	B2 <- B2[-1,]
	
	if(is.null(treat))  A3 <- cbind(c("Treatment",o1),c("Estimates",B2[,1]),c("95%CI",paste0("(",B2[,2],", ",B2[,3],")")))
	if(is.null(treat)==FALSE){
		A3 <- cbind(c("Treatment",treat),c("Estimates",B2[,1]),c("95%CI",paste0("(",B2[,2],", ",B2[,3],")")))
	}
		
	S1 <- A1[,2]^-2

	sp1 <- approxfun(c(min(S1),max(S1)),c(.2,.3),method = "linear")
	S2 <- c(NA,sp1(S1))

	zeropoint <- 0
	
	}

	if((xms=="OR")||(xms=="RR")){

	nmx <- nma(x,method=method,eform=TRUE)
	
	treat <- x$coding[,2]

	x.logform <- TRUE

	A1 <- nmx[[5]]

	o1 <- order(A1[,1])	
	if(ascending==FALSE) o1 <- order(A1[,1],decreasing = TRUE)
	treat <- treat[-1]; treat <- treat[o1]
	
	A1 <- A1[o1,]
	
	p <- dim(A1)[1]
	
	A2 <- data.frame(coef=A1[,1],low=A1[,3],high=A1[,4]); A2 <- rbind(c(NA,NA,NA),A2)
	
	A2 <- B2 <- round(A2,digits)
	
	for(i in 2:dim(A2)[1]){
		for(j in 1:dim(A2)[2]){
		
			B2[i,j] <- rdc(A2[i,j],digits)

		}
	}
	
	B2 <- B2[-1,]
	
	if(is.null(treat))  A3 <- cbind(c("Treatment",o1),c("Estimates",B2[,1]),c("95%CI",paste0("(",B2[,2],", ",B2[,3],")")))
	if(is.null(treat)==FALSE){
		A3 <- cbind(c("Treatment",treat),c("Estimates",B2[,1]),c("95%CI",paste0("(",B2[,2],", ",B2[,3],")")))
	}
		
	S1 <- A1[,2]^-2

	sp1 <- approxfun(c(min(S1),max(S1)),c(.2,.3),method = "linear")
	S2 <- c(NA,sp1(S1))
	
	zeropoint <- 1

	}

	forestplot(
		labeltext=A3,
		A2,
		boxsize = S2,
		align=c("l","c","c"),
		zero = zeropoint,
		xlog = x.logform,
		col = fpColors(lines = col.plot, box = col.plot), 
		txt_gp = fpTxtGp(ticks=gpar(cex=1)))

	#R1 <- list(labeltext=A3,coef=A2,boxsize = S2)
	
	#return(R1)
	
}

