global.ict <- function(x){

	xms <- x$measure

	if(xms=="OR"||xms=="RR"||xms=="RD"){

	study <- x$study
	treat <- x$treat
	d <- x$d
	n <- x$n
	
	study <- as.numeric(factor(study))
	
	y <- x$y			# Contrast-based statistics
	S <- x$S

	MY <- max(y,na.rm=TRUE) - min(y,na.rm=TRUE)

###

N <- dim(y)[1]
p <- max(treat) - 1

des <- n.arm <- numeric(N)
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
	
}

des0 <- sort(unique(des))
L <- length(des0)			# number of design

dof <- sum(nchar(des0)) - L - p			# degree-of-freedom

if(dof>=1){

count <- numeric(p+1)
X1 <- NULL

for(i in 1:p)  X1[[i]] <- t(matrix(numeric(N) + 1))
for(i in 1:p)  rownames(X1[[i]]) <- paste0(i+1,": cons")

for(i in 1:L){

	di <- des0[i]
	wi <- which(des==di)

	ti <- Ti[[wi[1]]]

	count[ti] <- count[ti] + 1
	
	wc2 <- which(count[ti]>=2)
	lc2 <- length(wc2)
	tc2 <- ti[wc2]
	
	if(lc2>=2){
	
		for(j in 2:lc2){

			xi <- numeric(N); xi[wi] <- 1
			xi <- t(matrix(xi)); rownames(xi) <- paste0(tc2[j],": des_",di)
			X1[[tc2[j]-1]] <- rbind(X1[[tc2[j]-1]],xi)
		
		}

	}

}

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

REMLIC <- function(y,S,X1,maxitr=200){
  
  N <- dim(y)[1]
  p <- dim(y)[2]
  
  Q <- 0
  for(i in 1:p) Q <- Q + dim(X1[[i]])[1]
  
  Qp <- numeric(p)
  for(i in 1:p) Qp[i] <- dim(X1[[i]])[1]
  Qc2 <- cumsum(Qp)
  Qc1 <- c(1,(Qc2[1:p-1]+1))
  L1 <- Qc2 - Qc1 + 1

  rnames <- rep(NA,times=Q)
  for(i in 1:p) rnames[Qc1[i]:Qc2[i]] <- rownames(X1[[i]])
    
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
    
    mu <- A1 %*% ginv2(A2)
    g1 <- optimize(LL1, lower = 0, upper = MY)$minimum
    g2 <- 0.5*g1
    
    V.mu <- ginv2(A2)
    
    Qc <- c(mu,g1,g2)
    
    rb <- abs(Qc - Qc0)/abs(Qc0); rb[is.nan(rb)] <- 0
    if(max(rb) < 10^-4) break
        
    Qc0 <- Qc
    
  }
  
  SE <- sqrt(diag(V.mu))
  
  R1 <- as.vector(mu)
  R2 <- as.vector(SE)
  R3 <- as.vector(mu - qnorm(.975)*SE)
  R4 <- as.vector(mu + qnorm(.975)*SE)
  P4 <- 2*(1-pnorm(abs(R1)/R2))
  
  R5 <- cbind(R1,R2,R3,R4,P4); colnames(R5) <- c("Coef.","SE","95%CL","95%CU","P-value"); rownames(R5) <- rnames
  
  R6 <- sqrt(g1)
  R7 <- g2/g1
  
  R8 <- as.numeric(mu[-Qc1]%*%ginv2(V.mu[-Qc1,-Qc1])%*%mu[-Qc1])			# ML-based Wald statistic for the global inconsistency test
  R9 <- length(mu[-Qc1])			# df
  R10 <- 1 - pchisq(R8,df=R9)
  
  R11 <- list("Coefficients"=R5,"Between-studies_SD"=R6,"Between-studies_COR"=R7,"X2-statistic"=R8,"df"=R9,"P-value"=R10)
    
  return(R11)
  
}

C1 <- REMLIC(y,S,X1)

C2 <- list("coding"=x$coding,"reference"=x$reference,"number of studies"=N,designs=des0,"Coefficients of the design-by-treatment interaction model"=C1[[1]],"Between-studies_SD"=C1[[2]],"Between-studies_COR"=C1[[3]],"X2-statistic"=C1[[4]],"df"=C1[[5]],"P-value"=C1[[6]])

return(C2)

}

if(dof==0)	return("The degree-of-freedom of this network is 0.")

}

	if(xms=="MD"||xms=="SMD"){

	study <- x$study
	treat <- x$treat
	m <- x$m
	s <- x$s
	n <- x$n
	
	study <- as.numeric(factor(study))
	
	y <- x$y			# Contrast-based statistics
	S <- x$S

	MY <- max(y,na.rm=TRUE) - min(y,na.rm=TRUE)

###

N <- dim(y)[1]
p <- max(treat) - 1

des <- n.arm <- numeric(N)
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
	
}

des0 <- sort(unique(des))
L <- length(des0)			# number of design

dof <- sum(nchar(des0)) - L - p			# degree-of-freedom

if(dof>=1){

count <- numeric(p+1)
X1 <- NULL

for(i in 1:p)  X1[[i]] <- t(matrix(numeric(N) + 1))
for(i in 1:p)  rownames(X1[[i]]) <- paste0(i+1,": cons")

for(i in 1:L){

	di <- des0[i]
	wi <- which(des==di)

	ti <- Ti[[wi[1]]]

	count[ti] <- count[ti] + 1
	
	wc2 <- which(count[ti]>=2)
	lc2 <- length(wc2)
	tc2 <- ti[wc2]
	
	if(lc2>=2){
	
		for(j in 2:lc2){

			xi <- numeric(N); xi[wi] <- 1
			xi <- t(matrix(xi)); rownames(xi) <- paste0(tc2[j],": des_",di)
			X1[[tc2[j]-1]] <- rbind(X1[[tc2[j]-1]],xi)
		
		}

	}

}

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

REMLIC <- function(y,S,X1,maxitr=200){
  
  N <- dim(y)[1]
  p <- dim(y)[2]
  
  Q <- 0
  for(i in 1:p) Q <- Q + dim(X1[[i]])[1]
  
  Qp <- numeric(p)
  for(i in 1:p) Qp[i] <- dim(X1[[i]])[1]
  Qc2 <- cumsum(Qp)
  Qc1 <- c(1,(Qc2[1:p-1]+1))
  L1 <- Qc2 - Qc1 + 1

  rnames <- rep(NA,times=Q)
  for(i in 1:p) rnames[Qc1[i]:Qc2[i]] <- rownames(X1[[i]])
    
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
    
    mu <- A1 %*% ginv2(A2)
    g1 <- optimize(LL1, lower = 0, upper = MY)$minimum
    g2 <- 0.5*g1
    
    V.mu <- ginv2(A2)
    
    Qc <- c(mu,g1,g2)
    
    rb <- abs(Qc - Qc0)/abs(Qc0); rb[is.nan(rb)] <- 0
    if(max(rb) < 10^-4) break
        
    Qc0 <- Qc
    
  }
  
  SE <- sqrt(diag(V.mu))
  
  R1 <- as.vector(mu)
  R2 <- as.vector(SE)
  R3 <- as.vector(mu - qnorm(.975)*SE)
  R4 <- as.vector(mu + qnorm(.975)*SE)
  P4 <- 2*(1-pnorm(abs(R1)/R2))
  
  R5 <- cbind(R1,R2,R3,R4,P4); colnames(R5) <- c("Coef.","SE","95%CL","95%CU","P-value"); rownames(R5) <- rnames
  
  R6 <- sqrt(g1)
  R7 <- g2/g1
  
  R8 <- as.numeric(mu[-Qc1]%*%ginv2(V.mu[-Qc1,-Qc1])%*%mu[-Qc1])			# ML-based Wald statistic for the global inconsistency test
  R9 <- length(mu[-Qc1])			# df
  R10 <- 1 - pchisq(R8,df=R9)
  
  R11 <- list("Coefficients"=R5,"Between-studies_SD"=R6,"Between-studies_COR"=R7,"X2-statistic"=R8,"df"=R9,"P-value"=R10)
    
  return(R11)
  
}

C1 <- REMLIC(y,S,X1)

C2 <- list("coding"=x$coding,"reference"=x$reference,"number of studies"=N,designs=des0,"Coefficients of the design-by-treatment interaction model"=C1[[1]],"Between-studies_SD"=C1[[2]],"Between-studies_COR"=C1[[3]],"X2-statistic"=C1[[4]],"df"=C1[[5]],"P-value"=C1[[6]])

return(C2)

}

if(dof==0)	return("The degree-of-freedom of this network is 0.")

}

}
