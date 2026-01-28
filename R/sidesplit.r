sidesplit <- function(x,digits=3){

	call <- match.call()

	xms <- x$measure

	MY <- max(x$y,na.rm=TRUE) - min(x$y,na.rm=TRUE)


	if(xms=="OR"||xms=="RR"||xms=="RD"||xms=="HR"||xms=="SPD"){

	study <- x$study
	treat <- x$treat
	d <- x$d
	n <- x$n
	
	study <- as.numeric(factor(study))
	
	data1 <- data.frame(study,treat,d,n)

####

vmat <- function(q2, p){
  
  i1 <- 1; i2 <- p
  
  Sg <- matrix(numeric(p*p),p,p)
  
  for(i in 1:p){
    
    Sg[i,(i:p)] <- q2[i1:i2]
    
    i1 <- i2 + 1 
    i2 <- i2 + p - i
    
  }
  
  Sg <- Sg + t(Sg); diag(Sg) <- diag(Sg)/2
  
  return(Sg)
  
}

cmat <- function(q2, p){
  
  Sg <- q2*diag(p)
  
  return(Sg)
  
}

pmat <- function(Si, wi){
  
  pl <- length(wi)
  
  R <- matrix(numeric(pl*pl),pl)
  
  for(i in 1:pl){
    for(j in 1:pl){
      
      R[i,j] <- Si[wi[i],wi[j]]
      
    }
  }
  
  return(R)
  
}

imat <- function(Si, wi, p){
  
  pl <- length(wi)
  
  R <- matrix(numeric(p*p),p)
  
  for(i in 1:pl){
    for(j in 1:pl){
      
      R[wi[i],wi[j]] <- Si[i,j]
      
    }
  }
  
  return(R)
  
}

ivec <- function(yi, wi, p){
  
  pl <- length(wi)
  
  R <- numeric(p)
  
  for(i in 1:pl) R[wi[i]] <- yi[i]
  
  return(R)
  
}

ivec2 <- function(yi, wi, p){
  
  pl <- length(wi)
  
  R <- rep(NA,times=p)
  
  for(i in 1:pl) R[wi[i]] <- yi[i]
  
  return(R)
  
}

gmat <- function(g1,g2,p){
  
  G <- diag(0, p) + g2
  diag(G) <- g1
  return(G)
  
}

QT <- function(x,x0){
  
  x1 <- sort(c(x,x0))
  w1 <- which(x1==as.numeric(x0))
  qt <- 1 - w1/(length(x)+1)
  return(qt)
  
}

####

REMLSS <- function(y,S,X1,maxitr=50){
  
  N <- dim(y)[1]
  p <- dim(y)[2]
  
  Q <- 0
  for(i in 1:p) Q <- Q + dim(X1[[i]])[1]
  
  Qp <- numeric(p)
  for(i in 1:p) Qp[i] <- dim(X1[[i]])[1]
  Qc2 <- cumsum(Qp)
  Qc1 <- c(1,(Qc2[1:p-1]+1))
  L1 <- Qc2 - Qc1 + 1

  mu <- rnorm(Q)	# initial values
  g1 <- 0.2
  g2 <- 0.1
  
  Qc0 <- c(mu,g1,g2)
  
  LL1 <- function(g){
    
    #G <- gmat(g,g2,p)
    G <- gmat(g,(g/2),p)
    
    ll1 <- 0
    
    for(i in 1:N){
      
      yi <- as.vector(y[i,])
      wi <- which(is.na(yi)==FALSE)
      pl <- length(wi)
      
      Si <- vmat(S[i,], p)
      
      yi <- yi[wi]
      Si <- pmat(Si,wi)
      Gi <- pmat(G,wi)

	  mui <- rep(NA,times=pl)

	  for(k in 1:pl){

	    j <- wi[k]

		qj <- Qc1[j]:Qc2[j]
		muj <- mu[qj]
		xj <- X1[[j]][,i]
		mui[k] <- muj%*%xj

	  }

      B1 <- (yi - mui)
      B2 <- ginv2(Gi + Si)
      
      A1 <- log(det(Gi + Si))
      A2 <- t(B1) %*% B2 %*% B1
      # A3 <- pl * log(2*pi)
      
      ll1 <- ll1 + A1 + A2 # + A3
      
    }
    
	A1 <- numeric(Q)
    A2 <- matrix(numeric(Q*Q),Q)
    
    for(i in 1:N){
      
      yi <- as.vector(y[i,])
      wi <- which(is.na(yi)==FALSE)
      pl <- length(wi)
      
      Si <- vmat(S[i,], p)
      
      yi <- yi[wi]
      Si <- pmat(Si,wi)
      Gi <- pmat(G,wi)
      
      Wi <- ginv2(Gi + Si)
 
	  Xj <- NULL

	  for(k in 1:pl){

		j <- wi[k]

		if(k==1){
			xj <- matrix(X1[[j]][,i])
			Xj <- xj
		}
		
		if(k>=2){
			xj <- matrix(X1[[j]][,i])
			dimj <- dim(Xj)
			qj <- length(xj)
			
			B1 <- matrix(numeric(dimj[2]*qj),qj)
			B2 <- matrix(numeric(dimj[1]))

			B3 <- rbind(Xj,B1)
			B4 <- rbind(B2,xj)
			
			Xj <- cbind(B3,B4)
		}

	  }
	  
	  a1 <- t(yi) %*% Wi %*% t(Xj)
	  a2 <- Xj %*% Wi %*% t(Xj)

	  L2 <- L1[wi]
	  
	  if(pl==1){

		j <- wi
		A2[Qc1[j]:Qc2[j],Qc1[j]:Qc2[j]] <- A2[Qc1[j]:Qc2[j],Qc1[j]:Qc2[j]] + a2
	  
	  }
	  	  
      if(pl>=2){

	   L3 <- cumsum(L2)			
       L4 <- c(1,(L3[1:(pl-1)]+1))	
  
		for(k in 1:pl){
			for(h in 1:pl){

				j1 <- wi[k]
				j2 <- wi[h]
		
				wj1 <- Qc1[j1]:Qc2[j1]
				wj2 <- Qc1[j2]:Qc2[j2]
		
				uj1 <- L4[k]:L3[k]
				uj2 <- L4[h]:L3[h]
		
				A2[wj1,wj2] <- A2[wj1,wj2] + a2[uj1,uj2]
			}
		}

	  }

	  for(k in 1:pl){

		j <- wi[k]

		Lk <- L1[j]
		
		A1[Qc1[j]:Qc2[j]] <- A1[Qc1[j]:Qc2[j]] + a1[1:Lk]

		dim2 <- length(a1)
		
		if(k!=pl){
			a1 <- a1[(Lk+1):dim2]
		}
		
	  }

    }
	
	ll1 <- ll1 + log(det(A2))
    
    return(ll1)
    
  }
 
  for(itr in 1:maxitr){
    
    A1 <- numeric(Q)
    A2 <- matrix(numeric(Q*Q),Q)
    
    G <- gmat(g1,g2,p)
    
    for(i in 1:N){
      
      yi <- as.vector(y[i,])
      wi <- which(is.na(yi)==FALSE)
      pl <- length(wi)
      
      Si <- vmat(S[i,], p)
      
      yi <- yi[wi]
      Si <- pmat(Si,wi)
      Gi <- pmat(G,wi)
      
      Wi <- ginv2(Gi + Si)
 
	  Xj <- NULL

	  for(k in 1:pl){

		j <- wi[k]

		if(k==1){
			xj <- matrix(X1[[j]][,i])
			Xj <- xj
		}
		
		if(k>=2){
			xj <- matrix(X1[[j]][,i])
			dimj <- dim(Xj)
			qj <- length(xj)
			
			B1 <- matrix(numeric(dimj[2]*qj),qj)
			B2 <- matrix(numeric(dimj[1]))

			B3 <- rbind(Xj,B1)
			B4 <- rbind(B2,xj)
			
			Xj <- cbind(B3,B4)
		}

	  }
	  
	  a1 <- t(yi) %*% Wi %*% t(Xj)
	  a2 <- Xj %*% Wi %*% t(Xj)

	  L2 <- L1[wi]
	  
	  if(pl==1){

		j <- wi
		A2[Qc1[j]:Qc2[j],Qc1[j]:Qc2[j]] <- A2[Qc1[j]:Qc2[j],Qc1[j]:Qc2[j]] + a2
	  
	  }
	  	  
      if(pl>=2){

	   L3 <- cumsum(L2)				
       L4 <- c(1,(L3[1:(pl-1)]+1))	
  
		for(k in 1:pl){
			for(h in 1:pl){

				j1 <- wi[k]
				j2 <- wi[h]
		
				wj1 <- Qc1[j1]:Qc2[j1]
				wj2 <- Qc1[j2]:Qc2[j2]
		
				uj1 <- L4[k]:L3[k]
				uj2 <- L4[h]:L3[h]
		
				A2[wj1,wj2] <- A2[wj1,wj2] + a2[uj1,uj2]
			}
		}

	  }

	  for(k in 1:pl){

		j <- wi[k]

		Lk <- L1[j]
		
		A1[Qc1[j]:Qc2[j]] <- A1[Qc1[j]:Qc2[j]] + a1[1:Lk]

		dim2 <- length(a1)
		
		if(k!=pl){
			a1 <- a1[(Lk+1):dim2]
		}
		
	  }

    }
    
	if(det(A2)>0){
	
		mu <- A1 %*% ginv2(A2)
		g1 <- optimize(LL1, lower = 0, upper = MY)$minimum
		g2 <- 0.5*g1
    
		V.mu <- ginv2(A2)
    
		Qc <- c(mu,g1,g2)
    
		rb <- abs(Qc - Qc0)/abs(Qc0); rb[is.nan(rb)] <- 0
		if(max(rb) < 10^-4) break
        
		Qc0 <- Qc
	
	}
	
	if(det(A2)<=0)	break
    
  }
  
  if(det(A2)>0){

	SE <- sqrt(diag(V.mu))
  
	R1 <- as.vector(mu)
	R2 <- as.vector(SE)
	R3 <- as.vector(mu - qnorm(.975)*SE)
	R4 <- as.vector(mu + qnorm(.975)*SE)
	P4 <- 2*(1-pnorm(abs(R1)/R2))

	R5 <- cbind(R1,R2,R3,R4,P4); colnames(R5) <- c("Coef.","SE","95%CL","95%CU","P-value")
  
	rname <- NULL
	for(i in 1:p) rname <- c(rname,rownames(X1[[i]]))

	kp <- which(rname=="i.direct")
	kq <- which(rname=="i.indirect")
  
	R6 <- R5[kp,]
	R7 <- R5[kq,]
  
	r1 <- R6[1] - R7[1]  
	r2 <- sqrt(R6[2]^2 + R7[2]^2)
	r3 <- r1 - qnorm(.975)*r2
	r4 <- r1 + qnorm(.975)*r2
	r5 <- 2*(1-pnorm(abs(r1)/r2))

	R8 <- data.frame(r1,r2,r3,r4,r5)
	colnames(R8) <- c("Coef.","SE","95%CL","95%CU","P-value")

	R9 <- list("Direct evidence"=R6,"Indirect evidence"=R7,"Difference"=R8)

	return(R9)
	
  }
  
  if(det(A2)<=0){

	R6 <- R7 <- R8 <- t(data.frame(c(NA,NA,NA,NA,NA)))
	colnames(R6) <- colnames(R7) <- colnames(R8) <- c("Coef.","SE","95%CL","95%CU","P-value")
	rownames(R6) <- rownames(R7) <- rownames(R8) <- NULL

  	R9 <- list("Direct evidence"=R6,"Indirect evidence"=R7,"Difference"=R8)
    
	return(R9)

  }
  
}

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

	R1 <- R2 <- R3 <- rname <- NULL

	for(k in 1:(p-1)){

		T2 <- ttrt(treat, ref=k)	
		treat2 <- T2$code
		
		data2 <- data.frame(study,treat2,d,n)

		edat <- setup(study=study,trt=treat2,d=d,n=n,measure=xms,ref=1,data=data2)
		#edat <- data.edit(study,treat2,d,n,data=data2)

		y1 <- edat$y
		S1 <- edat$S
		
		for(h in (k+1):p){

			w.pair <- NULL
				
			for(i in 1:N){
				
				tri <- Ti[[i]]
				if( (sum(tri==k) + sum(tri==h))==2 ) w.pair <- c(w.pair,i)
				
			}
			
			i.pair <- rep(FALSE,times=N)
			i.pair[w.pair] <- TRUE
	
			i.direct <- numeric(N); i.direct[w.pair] <- 1
			i.indirect <- numeric(N) + 1; i.indirect[w.pair] <- 0
			
			if(sum(i.pair)>=1){
			
				X1 <- NULL
				
				Q1 <- setdiff(1:q,(h-1))

				for(i in Q1)  X1[[i]] <- t(matrix(numeric(N) + 1))
				for(i in Q1)  rownames(X1[[i]]) <- paste0(i+1,": cons")
				
				X1[[h-1]] <- rbind(i.direct,i.indirect)
				
				wm <- which((i.pair)&(n.arm>=3))
				wd <- des[wm]
				uwd <- unique(wd)
				
				if(length(wm)>=1){
				
					for(l in 1:length(uwd)){
					
						des.l <- uwd[l]

						i.arm <- as.numeric(strsplit(des.l,"-")[[1]])

						wl <- setdiff(i.arm,c(k,h))
					
						for(a in wl){

							ww <- setdiff(which(str_detect(des, pattern=paste(a))),which(str_detect(des,pattern=des.l)))

							if(length(ww)>=1){

								i.adj <- numeric(N)
								i.adj[des==des.l] <- 1
								if(a>k) X1[[a-1]] <- rbind(X1[[a-1]] ,i.adj)
								if(a<k) X1[[a]] <- rbind(X1[[a]] ,i.adj)

							}

						}
									
					}				
				
				}
				
				R0 <- REMLSS(y1,S1,X1)

				R1 <- rbind(R1,R0[[1]])
				R2 <- rbind(R2,R0[[2]])
				R3 <- rbind(R3,R0[[3]])
				rname <- c(rname,paste0(k," vs. ",h))
				
				# message(paste0("The sidesplitting computation for the treatment pair ",pair," is completed."))

			}
			
		}
		
	}

	if(is.null(rname)==FALSE){
	
		# rname <- paste0("pair: ",rname)
		rownames(R1) <- rname
		rownames(R2) <- rname
		rownames(R3) <- rname

		if(sum(is.na(R3[,1]))>=1){
			R5 <- "For the NA components, the estimates were not calculable because the corresponding nodes were isolated or sufficient information did not exist."
			R4 <- list("coding"=x$coding,"reference"=x$reference,"Direct evidence"=R1,"Indirect evidence"=R2,"Difference"=R3,"Warning"=R5,digits=digits,call=call)
		}
		if(sum(is.na(R3[,1]))==0){
			R4 <- list("coding"=x$coding,"reference"=x$reference,"Direct evidence"=R1,"Indirect evidence"=R2,"Difference"=R3,"Warning"=NA,digits=digits,call=call)
		}

		class(R4) <- "sidesplit"  
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

vmat <- function(q2, p){
  
  i1 <- 1; i2 <- p
  
  Sg <- matrix(numeric(p*p),p,p)
  
  for(i in 1:p){
    
    Sg[i,(i:p)] <- q2[i1:i2]
    
    i1 <- i2 + 1 
    i2 <- i2 + p - i
    
  }
  
  Sg <- Sg + t(Sg); diag(Sg) <- diag(Sg)/2
  
  return(Sg)
  
}

cmat <- function(q2, p){
  
  Sg <- q2*diag(p)
  
  return(Sg)
  
}

pmat <- function(Si, wi){
  
  pl <- length(wi)
  
  R <- matrix(numeric(pl*pl),pl)
  
  for(i in 1:pl){
    for(j in 1:pl){
      
      R[i,j] <- Si[wi[i],wi[j]]
      
    }
  }
  
  return(R)
  
}

imat <- function(Si, wi, p){
  
  pl <- length(wi)
  
  R <- matrix(numeric(p*p),p)
  
  for(i in 1:pl){
    for(j in 1:pl){
      
      R[wi[i],wi[j]] <- Si[i,j]
      
    }
  }
  
  return(R)
  
}

ivec <- function(yi, wi, p){
  
  pl <- length(wi)
  
  R <- numeric(p)
  
  for(i in 1:pl) R[wi[i]] <- yi[i]
  
  return(R)
  
}

ivec2 <- function(yi, wi, p){
  
  pl <- length(wi)
  
  R <- rep(NA,times=p)
  
  for(i in 1:pl) R[wi[i]] <- yi[i]
  
  return(R)
  
}

gmat <- function(g1,g2,p){
  
  G <- diag(0, p) + g2
  diag(G) <- g1
  return(G)
  
}

QT <- function(x,x0){
  
  x1 <- sort(c(x,x0))
  w1 <- which(x1==as.numeric(x0))
  qt <- 1 - w1/(length(x)+1)
  return(qt)
  
}

####

REMLSS <- function(y,S,X1,maxitr=50){
  
  N <- dim(y)[1]
  p <- dim(y)[2]
  
  Q <- 0
  for(i in 1:p) Q <- Q + dim(X1[[i]])[1]
  
  Qp <- numeric(p)
  for(i in 1:p) Qp[i] <- dim(X1[[i]])[1]
  Qc2 <- cumsum(Qp)
  Qc1 <- c(1,(Qc2[1:p-1]+1))
  L1 <- Qc2 - Qc1 + 1

  mu <- rnorm(Q)	# initial values
  g1 <- 0.2
  g2 <- 0.1
  
  Qc0 <- c(mu,g1,g2)
  
  LL1 <- function(g){
    
    #G <- gmat(g,g2,p)
    G <- gmat(g,(g/2),p)
    
    ll1 <- 0
    
    for(i in 1:N){
      
      yi <- as.vector(y[i,])
      wi <- which(is.na(yi)==FALSE)
      pl <- length(wi)
      
      Si <- vmat(S[i,], p)
      
      yi <- yi[wi]
      Si <- pmat(Si,wi)
      Gi <- pmat(G,wi)

	  mui <- rep(NA,times=pl)

	  for(k in 1:pl){

	    j <- wi[k]

		qj <- Qc1[j]:Qc2[j]
		muj <- mu[qj]
		xj <- X1[[j]][,i]
		mui[k] <- muj%*%xj

	  }

      B1 <- (yi - mui)
      B2 <- ginv2(Gi + Si)
      
      A1 <- log(det(Gi + Si))
      A2 <- t(B1) %*% B2 %*% B1
      # A3 <- pl * log(2*pi)
      
      ll1 <- ll1 + A1 + A2 # + A3
      
    }
    
	A1 <- numeric(Q)
    A2 <- matrix(numeric(Q*Q),Q)
    
    for(i in 1:N){
      
      yi <- as.vector(y[i,])
      wi <- which(is.na(yi)==FALSE)
      pl <- length(wi)
      
      Si <- vmat(S[i,], p)
      
      yi <- yi[wi]
      Si <- pmat(Si,wi)
      Gi <- pmat(G,wi)
      
      Wi <- ginv2(Gi + Si)
 
	  Xj <- NULL

	  for(k in 1:pl){

		j <- wi[k]

		if(k==1){
			xj <- matrix(X1[[j]][,i])
			Xj <- xj
		}
		
		if(k>=2){
			xj <- matrix(X1[[j]][,i])
			dimj <- dim(Xj)
			qj <- length(xj)
			
			B1 <- matrix(numeric(dimj[2]*qj),qj)
			B2 <- matrix(numeric(dimj[1]))

			B3 <- rbind(Xj,B1)
			B4 <- rbind(B2,xj)
			
			Xj <- cbind(B3,B4)
		}

	  }
	  
	  a1 <- t(yi) %*% Wi %*% t(Xj)
	  a2 <- Xj %*% Wi %*% t(Xj)

	  L2 <- L1[wi]
	  
	  if(pl==1){

		j <- wi
		A2[Qc1[j]:Qc2[j],Qc1[j]:Qc2[j]] <- A2[Qc1[j]:Qc2[j],Qc1[j]:Qc2[j]] + a2
	  
	  }
	  	  
      if(pl>=2){

	   L3 <- cumsum(L2)	
       L4 <- c(1,(L3[1:(pl-1)]+1))	
  
		for(k in 1:pl){
			for(h in 1:pl){

				j1 <- wi[k]
				j2 <- wi[h]
		
				wj1 <- Qc1[j1]:Qc2[j1]
				wj2 <- Qc1[j2]:Qc2[j2]
		
				uj1 <- L4[k]:L3[k]
				uj2 <- L4[h]:L3[h]
		
				A2[wj1,wj2] <- A2[wj1,wj2] + a2[uj1,uj2]
			}
		}

	  }

	  for(k in 1:pl){

		j <- wi[k]

		Lk <- L1[j]
		
		A1[Qc1[j]:Qc2[j]] <- A1[Qc1[j]:Qc2[j]] + a1[1:Lk]

		dim2 <- length(a1)
		
		if(k!=pl){
			a1 <- a1[(Lk+1):dim2]
		}
		
	  }

    }
	
	ll1 <- ll1 + log(det(A2))
    
    return(ll1)
    
  }
 
  for(itr in 1:maxitr){
    
    A1 <- numeric(Q)
    A2 <- matrix(numeric(Q*Q),Q)
    
    G <- gmat(g1,g2,p)
    
    for(i in 1:N){
      
      yi <- as.vector(y[i,])
      wi <- which(is.na(yi)==FALSE)
      pl <- length(wi)
      
      Si <- vmat(S[i,], p)
      
      yi <- yi[wi]
      Si <- pmat(Si,wi)
      Gi <- pmat(G,wi)
      
      Wi <- ginv2(Gi + Si)
 
	  Xj <- NULL

	  for(k in 1:pl){

		j <- wi[k]

		if(k==1){
			xj <- matrix(X1[[j]][,i])
			Xj <- xj
		}
		
		if(k>=2){
			xj <- matrix(X1[[j]][,i])
			dimj <- dim(Xj)
			qj <- length(xj)
			
			B1 <- matrix(numeric(dimj[2]*qj),qj)
			B2 <- matrix(numeric(dimj[1]))

			B3 <- rbind(Xj,B1)
			B4 <- rbind(B2,xj)
			
			Xj <- cbind(B3,B4)
		}

	  }
	  
	  a1 <- t(yi) %*% Wi %*% t(Xj)
	  a2 <- Xj %*% Wi %*% t(Xj)

	  L2 <- L1[wi]
	  
	  if(pl==1){

		j <- wi
		A2[Qc1[j]:Qc2[j],Qc1[j]:Qc2[j]] <- A2[Qc1[j]:Qc2[j],Qc1[j]:Qc2[j]] + a2
	  
	  }
	  	  
      if(pl>=2){

	   L3 <- cumsum(L2)		
       L4 <- c(1,(L3[1:(pl-1)]+1))	
  
		for(k in 1:pl){
			for(h in 1:pl){

				j1 <- wi[k]
				j2 <- wi[h]
		
				wj1 <- Qc1[j1]:Qc2[j1]
				wj2 <- Qc1[j2]:Qc2[j2]
		
				uj1 <- L4[k]:L3[k]
				uj2 <- L4[h]:L3[h]
		
				A2[wj1,wj2] <- A2[wj1,wj2] + a2[uj1,uj2]
			}
		}

	  }

	  for(k in 1:pl){

		j <- wi[k]

		Lk <- L1[j]
		
		A1[Qc1[j]:Qc2[j]] <- A1[Qc1[j]:Qc2[j]] + a1[1:Lk]

		dim2 <- length(a1)
		
		if(k!=pl){
			a1 <- a1[(Lk+1):dim2]
		}
		
	  }

    }
    
    
	if(det(A2)>0){
	
		mu <- A1 %*% ginv2(A2)
		g1 <- optimize(LL1, lower = 0, upper = MY)$minimum
		g2 <- 0.5*g1
    
		V.mu <- ginv2(A2)
    
		Qc <- c(mu,g1,g2)
    
		rb <- abs(Qc - Qc0)/abs(Qc0); rb[is.nan(rb)] <- 0
		if(max(rb) < 10^-4) break
        
		Qc0 <- Qc
	
	}
	
	if(det(A2)<=0)	break
    
  }
  
  if(det(A2)>0){

	SE <- sqrt(diag(V.mu))
  
	R1 <- as.vector(mu)
	R2 <- as.vector(SE)
	R3 <- as.vector(mu - qnorm(.975)*SE)
	R4 <- as.vector(mu + qnorm(.975)*SE)
	P4 <- 2*(1-pnorm(abs(R1)/R2))

	R5 <- cbind(R1,R2,R3,R4,P4); colnames(R5) <- c("Coef.","SE","95%CL","95%CU","P-value")
  
	rname <- NULL
	for(i in 1:p) rname <- c(rname,rownames(X1[[i]]))

	kp <- which(rname=="i.direct")
	kq <- which(rname=="i.indirect")
  
	R6 <- R5[kp,]
	R7 <- R5[kq,]
  
	r1 <- R6[1] - R7[1]  
	r2 <- sqrt(R6[2]^2 + R7[2]^2)
	r3 <- r1 - qnorm(.975)*r2
	r4 <- r1 + qnorm(.975)*r2
	r5 <- 2*(1-pnorm(abs(r1)/r2))

	R8 <- data.frame(r1,r2,r3,r4,r5)
	colnames(R8) <- c("Coef.","SE","95%CL","95%CU","P-value")

	R9 <- list("Direct evidence"=R6,"Indirect evidence"=R7,"Difference"=R8)
    
	return(R9)
	
  }
  
  if(det(A2)<=0){

	R6 <- R7 <- R8 <- t(data.frame(c(NA,NA,NA,NA,NA)))
	colnames(R6) <- colnames(R7) <- colnames(R8) <- c("Coef.","SE","95%CL","95%CU","P-value")
	rownames(R6) <- rownames(R7) <- rownames(R8) <- NULL

  	R9 <- list("Direct evidence"=R6,"Indirect evidence"=R7,"Difference"=R8)
    
	return(R9)

  }
  
}

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

	R1 <- R2 <- R3 <- rname <- NULL

	for(k in 1:(p-1)){

		T2 <- ttrt(treat, ref=k)
		treat2 <- T2$code
		
		data2 <- data.frame(study,treat2,m,s,n)

		edat <- setup(study=study,trt=treat2,m=m,s=s,n=n,measure=xms,ref=1,data=data2)
		#edat <- data.edit(study,treat2,d,n,data=data2)

		y1 <- edat$y
		S1 <- edat$S
		
		for(h in (k+1):p){

			w.pair <- NULL
				
			for(i in 1:N){
				
				tri <- Ti[[i]]
				if( (sum(tri==k) + sum(tri==h))==2 ) w.pair <- c(w.pair,i)
				
			}
			
			i.pair <- rep(FALSE,times=N)
			i.pair[w.pair] <- TRUE
	
			i.direct <- numeric(N); i.direct[w.pair] <- 1
			i.indirect <- numeric(N) + 1; i.indirect[w.pair] <- 0
					
			if(sum(i.pair)>=1){
			
				X1 <- NULL
				
				Q1 <- setdiff(1:q,(h-1))

				for(i in Q1)  X1[[i]] <- t(matrix(numeric(N) + 1))
				for(i in Q1)  rownames(X1[[i]]) <- paste0(i+1,": cons")
				
				X1[[h-1]] <- rbind(i.direct,i.indirect)
				
				wm <- which((i.pair)&(n.arm>=3))
				wd <- des[wm]
				uwd <- unique(wd)
				
				if(length(wm)>=1){
				
					for(l in 1:length(uwd)){
					
						des.l <- uwd[l]

						i.arm <- as.numeric(strsplit(des.l,"-")[[1]])

						wl <- setdiff(i.arm,c(k,h))
					
						for(a in wl){

							ww <- setdiff(which(str_detect(des, pattern=paste(a))),which(str_detect(des,pattern=des.l)))

							if(length(ww)>=1){

								i.adj <- numeric(N)
								i.adj[des==des.l] <- 1
								if(a>k) X1[[a-1]] <- rbind(X1[[a-1]] ,i.adj)
								if(a<k) X1[[a]] <- rbind(X1[[a]] ,i.adj)

							}

						}
									
					}				
				
				}
				
				R0 <- REMLSS(y1,S1,X1)

				R1 <- rbind(R1,R0[[1]])
				R2 <- rbind(R2,R0[[2]])
				R3 <- rbind(R3,R0[[3]])
				rname <- c(rname,paste0(k," vs. ",h))
				
				# message(paste0("The sidesplitting computation for the treatment pair ",pair," is completed."))

			}
			
		}
		
	}

	if(is.null(rname)==FALSE){
	
		# rname <- paste0("pair: ",rname)
		rownames(R1) <- rname
		rownames(R2) <- rname
		rownames(R3) <- rname

		if(sum(is.na(R3[,1]))>=1){
			R5 <- "For the NA components, the estimates were not calculable because the corresponding nodes were isolated or sufficient information did not exist."
			R4 <- list("coding"=x$coding,"reference"=x$reference,"Direct evidence"=R1,"Indirect evidence"=R2,"Difference"=R3,"Warning"=R5,digits=digits,call=call)
		}
		if(sum(is.na(R3[,1]))==0){
			R4 <- list("coding"=x$coding,"reference"=x$reference,"Direct evidence"=R1,"Indirect evidence"=R2,"Difference"=R3,"Warning"=NA,digits=digits,call=call)
		}

		class(R4) <- "sidesplit"  
		return(R4)

	}
		
	if(is.null(rname)){
	
		R9 <- "There are no corresponding pairs that can be analyzed on the network."
		return(R9)
		
	}
	
	}
	
}

	
print.sidesplit <- function(x, digits = x$digits, ...) {

  cat("Call:\n")
  print(x$call,row.names=FALSE)
  cat("\n")
  
  cat("Coding:\n", sep = "")
  print(x$coding,row.names=FALSE)
  cat("\n")

  cat("Reference: ", sep = "")
  cat(x[[2]])
  cat("\n")
  cat("\n")
  
  cat("Direct evidence: ", sep = "")
  cat("\n")
  A <- x[[3]]
  ##
  est <- round(A[,1],digits)
  SE <- round(A[,2],digits)
  Lower <- round(A[,3],digits)
  Upper <- round(A[,4],digits)
  pval <- round(A[,5],digits)
  TAB <- cbind(
    "Est." = est,
    "SE" = SE,
	"Lower" = Lower,
	"Upper" = Upper,
    "Pr(>|z|)"  = pval
  )
  rownames(TAB) <- rownames(A)
  print(TAB)
  cat("\n")

  cat("Indirect evidence: ", sep = "")
  cat("\n")
  A <- x[[4]]
  ##
  est <- round(A[,1],digits)
  SE <- round(A[,2],digits)
  Lower <- round(A[,3],digits)
  Upper <- round(A[,4],digits)
  pval <- round(A[,5],digits)
  TAB <- cbind(
    "Est." = est,
    "SE" = SE,
	"Lower" = Lower,
	"Upper" = Upper,
    "Pr(>|z|)"  = pval
  )
  rownames(TAB) <- rownames(A)
  print(TAB)
  cat("\n")

  cat("Difference between direct and indirect evidence with inconsistency test: ", sep = "")
  cat("\n")
  A <- x[[5]]
  ##
  est <- round(A[,1],digits)
  SE <- round(A[,2],digits)
  Lower <- round(A[,3],digits)
  Upper <- round(A[,4],digits)
  pval <- round(A[,5],digits)
  TAB <- cbind(
    "Est." = est,
    "SE" = SE,
	"Lower" = Lower,
	"Upper" = Upper,
    "Pr(>|z|)"  = pval
  )
  rownames(TAB) <- rownames(A)
  print(TAB)
  cat("\n")
  
  if(is.na(x$Warning)==FALSE){
	cat("Note: ")
	cat(x$Warning, sep = "")
	cat("\n")
	cat("\n")
  }
  
  invisible(x)
  
}
