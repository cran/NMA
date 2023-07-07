nmaleague <- function(x, method="NH", eform=FALSE, digits=3, PI=FALSE, out.csv=NULL){	

	xms <- x$measure	

	if(xms=="OR"||xms=="RR"||xms=="RD"){

	study <- x$study
	treat <- x$treat
	d <- x$d
	n <- x$n

	treat1 <- sort(unique(treat))
	p <- length(treat1)
	
	data1 <- data.frame(study,treat,d,n)
	
	R <- matrix(rep(NA,times=p*p),p)

	diag(R) <- x$coding[,2]

	for(k in 1:p){

		T1 <- ttrt(treat, ref=treat1[k])	
		data1$treat1 <- T1$code

		edat <- setup(study=study,trt=treat1,d=d,n=n,measure=xms,ref=1,data=data1)
		
		if(PI==FALSE){

			R1 <- nma(edat, method=method, eform=eform)[[5]]
		
			Q1 <- rep(NA,times=p-1)
			for(j in 1:(p-1))  Q1[j] <- paste0(rdc(R1[j,1],digits)," (",rdc(R1[j,3],digits),", ",rdc(R1[j,4],digits),")")
		
			R[k,-k] <- Q1

		}
		
		if(PI==TRUE){

			R1 <- nma(edat, method=method, eform=eform)[[11]]
		
			Q1 <- rep(NA,times=p-1)
			for(j in 1:(p-1))  Q1[j] <- paste0("(",rdc(R1[j,1],digits),", ",rdc(R1[j,2],digits),")")
		
			R[k,-k] <- Q1

		}
		
	}
	
	if(is.null(out.csv)==FALSE)  write.csv(R, file=out.csv, row.names = FALSE)
	
	return(R)
	
	}

	if(xms=="MD"||xms=="SMD"){

	study <- x$study
	treat <- x$treat
	m <- x$m
	s <- x$s
	n <- x$n

	treat1 <- sort(unique(treat))
	p <- length(treat1)
	
	data1 <- data.frame(study,treat,m,s,n)
	
	R <- matrix(rep(NA,times=p*p),p)

	diag(R) <- x$coding[,2]

	for(k in 1:p){

		T1 <- ttrt(treat, ref=treat1[k])	
		data1$treat1 <- T1$code

		edat <- setup(study=study,trt=treat1,m=m,s=s,n=n,measure=xms,ref=1,data=data1)
		
		if(PI==FALSE){

			R1 <- nma(edat, method=method, eform=eform)[[5]]
		
			Q1 <- rep(NA,times=p-1)
			for(j in 1:(p-1))  Q1[j] <- paste0(rdc(R1[j,1],digits)," (",rdc(R1[j,3],digits),", ",rdc(R1[j,4],digits),")")
		
			R[k,-k] <- Q1

		}
		
		if(PI==TRUE){

			R1 <- nma(edat, method=method, eform=eform)[[11]]
		
			Q1 <- rep(NA,times=p-1)
			for(j in 1:(p-1))  Q1[j] <- paste0("(",rdc(R1[j,1],digits),", ",rdc(R1[j,2],digits),")")
		
			R[k,-k] <- Q1

		}
		
	}
	
	if(is.null(out.csv)==FALSE)  write.csv(R, file=out.csv, row.names = FALSE)
	
	return(R)
	
	}
	
}

