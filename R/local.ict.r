local.ict <- function(x){

	xms <- x$measure

	if(xms=="OR"||xms=="RR"||xms=="RD"||xms=="HR"||xms=="SPD"){

	study <- x$study
	treat <- x$treat
	d <- x$d
	n <- x$n
	
	data1 <- data.frame(study,treat,d,n)

	treat1 <- sort(unique(treat))
	study1 <- as.numeric(factor(study))
	
	N <- length(levels(factor(study)))
	p <- max(treat)
	q <- p - 1

	des <- n.arm <- numeric(N)
	Ti <- NULL
	study.i <- rep(NA,times=length(study))

	for(i in 1:N){

		wi <- which(study1==i)
		ti <- sort(treat[wi],decreasing=FALSE)
		Ti[[i]] <- ti

		di <- ti[1]
		for(j in 2:length(wi)) di <- paste0(di,"-",ti[j])
		des[i] <- di
		n.arm[i] <- length(wi)
		study.i[wi] <- i
	
	}

	des0 <- sort(unique(des))

	rname <- R9 <- NULL

	for(k in 1:(p-2)){
		for(h in (k+1):(p-1)){
			for(l in (h+1):p){
			
				i.pair1 <- i.pair2 <- i.pair3 <- i.tri <- NULL
				
				for(i in 1:N){
				
					tri <- Ti[[i]]
					if( (sum(tri==k) + sum(tri==h))==2 ) i.pair1 <- c(i.pair1,i)
					if( (sum(tri==h) + sum(tri==l))==2 ) i.pair2 <- c(i.pair2,i)
					if( (sum(tri==l) + sum(tri==k))==2 ) i.pair3 <- c(i.pair3,i)
					if( (sum(tri==k) + sum(tri==h) + sum(tri==l))==3 ) i.tri <- c(i.tri,i)
				
				
				}
			
				w.pair1 <- setdiff(i.pair1,i.tri)
				w.pair2 <- setdiff(i.pair2,i.tri)
				w.pair3 <- setdiff(i.pair3,i.tri)
				w.tri <- i.tri
				
				cond1 <- (length(w.pair1)>=1)&&(length(w.pair2)>=1)&&(length(w.pair3)>=1)
				cond2 <- (length(w.pair1)>=1)&&(length(w.tri)>=1)
				cond3 <- (length(w.pair2)>=1)&&(length(w.tri)>=1)
				cond4 <- (length(w.pair3)>=1)&&(length(w.tri)>=1)
				
				if(cond1||cond2||cond3||cond4){
				
					w5 <- sort(unique(c(w.pair1,w.pair2,w.pair3,w.tri)))
					
					w6 <- NULL
					for(i in w5) w6 <- c(w6,which(study.i==i))

					w7 <- which( (treat==k)|(treat==h)|(treat==l) )
					
					w8 <- intersect(w6,w7)
					
					study8 <- as.numeric(factor(study[w8]))
					treat8 <- as.numeric(factor(treat[w8]))
					d8 <- d[w8]
					n8 <- n[w8]
					
					data8 <- data.frame(study8,treat8,d8,n8)

					edat8 <- setup(study=study8,trt=treat8,d=d8,n=n8,measure=xms,ref=1,data=data8)
					R8 <- global.ict(edat8)
					#R8 <- global.ict(study8, treat8, d8, n8,data=data8)

					R9 <- rbind(R9,c(R8[[3]],R8[[6]],R8[[8]],R8[[9]],R8[[10]]))
					rname <- c(rname,paste0(k,"-",h,"-",l))
				
				}
			}
		}
	}
	
	if(is.null(rname)){
	
		R9 <- "There are no closed loops on the network."
		return(R9)
		
	}

	if(is.null(rname)==FALSE){
	
		rname <- paste0("loop: ",rname)
		rownames(R9) <- rname
	
		colnames(R9) <- c("N","tau","X2-statistic","df","P-value")
	
		R10 <- list("coding"=x$coding,"reference"=x$reference,"loop inconsistency tests"=R9)
	
		return(R10)
		
	}

	}

	if(xms=="MD"||xms=="SMD"){

	study <- x$study
	treat <- x$treat
	m <- x$m
	s <- x$s
	n <- x$n
	
	data1 <- data.frame(study,treat,m,s,n)

	treat1 <- sort(unique(treat))
	study1 <- as.numeric(factor(study))
	
	N <- length(levels(factor(study)))
	p <- max(treat)
	q <- p - 1

	des <- n.arm <- numeric(N)
	Ti <- NULL
	study.i <- rep(NA,times=length(study))

	for(i in 1:N){

		wi <- which(study1==i)
		ti <- sort(treat[wi],decreasing=FALSE)
		Ti[[i]] <- ti

		di <- NULL
		for(j in 1:length(wi)) di <- paste0(di,ti[j])
		des[i] <- di
		n.arm[i] <- length(wi)
		study.i[wi] <- i
	
	}

	des0 <- sort(unique(des))

	rname <- R9 <- NULL

	for(k in 1:(p-2)){
		for(h in (k+1):(p-1)){
			for(l in (h+1):p){
			
				i.pair1 <- i.pair2 <- i.pair3 <- i.tri <- NULL
				
				for(i in 1:N){
				
					tri <- Ti[[i]]
					if( (sum(tri==k) + sum(tri==h))==2 ) i.pair1 <- c(i.pair1,i)
					if( (sum(tri==h) + sum(tri==l))==2 ) i.pair2 <- c(i.pair2,i)
					if( (sum(tri==l) + sum(tri==k))==2 ) i.pair3 <- c(i.pair3,i)
					if( (sum(tri==k) + sum(tri==h) + sum(tri==l))==3 ) i.tri <- c(i.tri,i)
				
				
				}
			
				w.pair1 <- setdiff(i.pair1,i.tri)
				w.pair2 <- setdiff(i.pair2,i.tri)
				w.pair3 <- setdiff(i.pair3,i.tri)
				w.tri <- i.tri
				
				cond1 <- (length(w.pair1)>=1)&&(length(w.pair2)>=1)&&(length(w.pair3)>=1)
				cond2 <- (length(w.pair1)>=1)&&(length(w.tri)>=1)
				cond3 <- (length(w.pair2)>=1)&&(length(w.tri)>=1)
				cond4 <- (length(w.pair3)>=1)&&(length(w.tri)>=1)
				
				if(cond1||cond2||cond3||cond4){
				
					w5 <- sort(unique(c(w.pair1,w.pair2,w.pair3,w.tri)))
					
					w6 <- NULL
					for(i in w5) w6 <- c(w6,which(study.i==i))

					w7 <- which( (treat==k)|(treat==h)|(treat==l) )
					
					w8 <- intersect(w6,w7)
					
					study8 <- as.numeric(factor(study[w8]))
					treat8 <- as.numeric(factor(treat[w8]))
					m8 <- m[w8]
					s8 <- s[w8]
					n8 <- n[w8]
					
					data8 <- data.frame(study8,treat8,m8,s8,n8)

					edat8 <- setup(study=study8,trt=treat8,m=m8,s=s8,n=n8,measure=xms,ref=1,data=data8)
					R8 <- global.ict(edat8)
					#R8 <- global.ict(study8, treat8, d8, n8,data=data8)

					R9 <- rbind(R9,c(R8[[3]],R8[[6]],R8[[8]],R8[[9]],R8[[10]]))
					rname <- c(rname,paste0(k,"-",h,"-",l))
				
				}
			}
		}
	}
	
	if(is.null(rname)){
	
		R9 <- "There are no closed loops on the network."
		return(R9)
		
	}

	if(is.null(rname)==FALSE){
	
		rname <- paste0("loop: ",rname)
		rownames(R9) <- rname
	
		colnames(R9) <- c("N","tau","X2-statistic","df","P-value")
		
		R10 <- list("coding"=x$coding,"reference"=x$reference,"loop inconsistency tests"=R9)
	
		return(R10)
		
	}

	}
				
}
