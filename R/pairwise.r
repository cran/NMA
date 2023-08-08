pairwise <- function(x,method="REML"){

	xms <- x$measure

	if(xms=="OR"||xms=="RR"||xms=="RD"){

	study <- x$study
	treat <- x$treat
	d <- x$d
	n <- x$n
	
	study <- as.numeric(factor(study))
	
	data1 <- data.frame(study,treat,d,n)

    ####

	treat1 <- sort(unique(treat))
	
	N <- length(unique(study))
	p <- max(treat)
	q <- p - 1

	des <- n.arm <- numeric(N)
	Ti <- NULL

	for(i in 1:N){

		wi <- which(study==i)
		ti <- sort(treat[wi],decreasing=FALSE)
		Ti[[i]] <- ti

		di <- ti[1]
		for(j in 2:length(wi)) di <- paste0(di,"-",ti[j])
		des[i] <- di
		n.arm[i] <- length(wi)
	
	}
	
	des0 <- sort(unique(des))

	Q1 <- rname <- NULL

	for(k in 1:(p-1)){
		
		for(h in (k+1):p){

			pair <- paste0(k,"-",h)
			i.pair <- str_detect(des, pattern=pair)
			w.pair <- which(i.pair)

			n.i <- sum(i.pair)
		
			if(sum(i.pair)>=1){
			
					d1 <- d2 <- n1 <- n2 <- NULL
				
					for(l in 1:n.i){
					
						w.l <- w.pair[l]
						w1 <- which((data1$study==w.l)&(data1$treat==k))
						w2 <- which((data1$study==w.l)&(data1$treat==h))

						d1 <- c(d1,data1$d[w1])
						d2 <- c(d2,data1$d[w2])
						n1 <- c(n1,data1$n[w1])
						n2 <- c(n2,data1$n[w2])
			
					}
					
					ai <- d1
					bi <- n1 - d1
					ci <- d2
					di <- n2 - d2
					
					ei <- escalc(xms,ai=ci,bi=di,ci=ai,di=bi)
					rmi <- rma(ei$yi,ei$vi,method=method)
					
					if(n.i>=3)  egi <- regtest(rmi, model="lm")
					
					R1 <- c(rmi$beta,rmi$ci.lb,rmi$ci.ub)
					if(xms=="OR"||xms=="RR")		R1 <- exp(R1)
					if(n.i<3)  R2 <- c(rmi$pval,rmi$tau2,.01*rmi$I2,rmi$H2,NA)
					if(n.i>=3)  R2 <- c(rmi$pval,rmi$tau2,.01*rmi$I2,rmi$H2,egi$pval)

					R3 <- c(n.i,R1,R2)
					
					Q1 <- rbind(Q1,R3)
					rname <- c(rname,paste0(k," vs. ",h))
					
			}
			
		}
		
	}
	
	if(is.null(rname)==FALSE){
	
		rownames(Q1) <- rname
		
		Q2 <- Q1[,1:5]
		colnames(Q2) <- c("N","estimate","95%CL","95%CU","P-value")
		
		Q3 <- Q1[,c(1,6:8)]
		colnames(Q3) <- c("N","tau^2","I^2","H^2")
		
		Q4 <- t(t(Q1[,c(1,9)]))
		colnames(Q4) <- c("N","P-value"	)
		
		R4 <- list("coding"=x$coding,"measure"=xms,"Summary effect measures"=Q2,"Heterogeneity measures"=Q3,"Egger test"=Q4)
	
		return(R4)
	
	}
	
	if(is.null(rname)){
	
		R9 <- "There are no corresponding pairs that can be analyzed on the network."
		return(R9)
		
	}
	
	}
	
	
	if(xms=="MD"||xms=="SMD"){

	study <- x$study
	treat <- x$treat
	m <- x$m
	s <- x$s
	n <- x$n
	
	study <- as.numeric(factor(study))
	
	data1 <- data.frame(study,treat,m,s,n)

    ####

	treat1 <- sort(unique(treat))
	
	N <- length(unique(study))
	p <- max(treat)
	q <- p - 1

	des <- n.arm <- numeric(N)
	Ti <- NULL

	for(i in 1:N){

		wi <- which(study==i)
		ti <- sort(treat[wi],decreasing=FALSE)
		Ti[[i]] <- ti

		di <- ti[1]
		for(j in 2:length(wi)) di <- paste0(di,"-",ti[j])
		des[i] <- di
		n.arm[i] <- length(wi)
	
	}
	
	des0 <- sort(unique(des))

	Q1 <- rname <- NULL

	for(k in 1:(p-1)){
		
		for(h in (k+1):p){

			pair <- paste0(k,"-",h)
			i.pair <- str_detect(des, pattern=pair)
			w.pair <- which(i.pair)

			n.i <- sum(i.pair)
		
			if(sum(i.pair)>=1){
			
					m1 <- m2 <- s1 <- s2 <- n1 <- n2 <- NULL
				
					for(l in 1:n.i){
					
						w.l <- w.pair[l]
						w1 <- which((data1$study==w.l)&(data1$treat==k))
						w2 <- which((data1$study==w.l)&(data1$treat==h))

						m1 <- c(m1,data1$m[w1])
						m2 <- c(m2,data1$m[w2])
						s1 <- c(s1,data1$s[w1])
						s2 <- c(s2,data1$s[w2])
						n1 <- c(n1,data1$n[w1])
						n2 <- c(n2,data1$n[w2])
			
					}
					
					ei <- escalc(xms,m1i=m2,sd1i=s2,n1i=n2,m2i=m1,sd2i=s1,n2i=n1)
					rmi <- rma(ei$yi,ei$vi,method=method)
					
					if(n.i>=3)  egi <- regtest(rmi, model="lm")
					
					R1 <- c(rmi$beta,rmi$ci.lb,rmi$ci.ub)
					if(n.i<3)  R2 <- c(rmi$pval,rmi$tau2,.01*rmi$I2,rmi$H2,NA)
					if(n.i>=3)  R2 <- c(rmi$pval,rmi$tau2,.01*rmi$I2,rmi$H2,egi$pval)

					R3 <- c(n.i,R1,R2)
					
					Q1 <- rbind(Q1,R3)
					rname <- c(rname,paste0(k," vs. ",h))
					
			}
			
		}
		
	}
	
	if(is.null(rname)==FALSE){
	
		rownames(Q1) <- rname
		
		Q2 <- Q1[,1:5]
		colnames(Q2) <- c("N","estimate","95%CL","95%CU","P-value")
		
		Q3 <- Q1[,c(1,6:8)]
		colnames(Q3) <- c("N","tau^2","I^2","H^2")
		
		Q4 <- t(t(Q1[,c(1,9)]))
		colnames(Q4) <- c("N","P-value"	)
		
		R4 <- list("coding"=x$coding,"measure"=xms,"Summary effect measures"=Q2,"Heterogeneity measures"=Q3,"Egger test"=Q4)
	
		return(R4)
	
	}
	
	if(is.null(rname)){
	
		R9 <- "There are no corresponding pairs that can be analyzed on the network."
		return(R9)
		
	}
	
	}
		
}

