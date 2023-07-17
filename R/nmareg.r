nmareg <- function(x, z, treats){

	xms <- x$measure

	Z <- x$Z
	covariate <- x$covariate

	###

	sz <- substitute(z)

	z1 <- deparse(substitute(z))
	
	z2 <- gsub(" ","",z1)
	
	if(substring(z2,1,2)!="c(")  covname <- z1
	  
	if(substring(z2,1,2)=="c("){

		nz <- nchar(z2)
		z3 <- substring(z2,3,(nz-1))
		covname <- strsplit(z3,",")[[1]]
		
    }
  
	###

	study <- x$study
	treat <- x$treat
	n <- x$n

	###
	
	study <- as.numeric(factor(study))
	
	y <- x$y			# Contrast-based statistics
	S <- x$S

	MY <- max(y,na.rm=TRUE) - min(y,na.rm=TRUE)

	###

	N <- dim(y)[1]
	p <- max(treat) - 1

	q1 <- length(covname)
	c1 <- rep(NA,times=q1)
	for(i in 1:q1)  c1[i] <- which(covariate==covname[i])

	Z1 <- data.frame(Z[,c1])

	z1 <- numeric(N)

	for(i in 1:N){

		wi <- which(study==i)
		z1[i] <- sum(Z1[wi,1]*n[wi])/sum(n[wi])
	
	}
	
	Z2 <- data.frame(z1)
	
	if(q1>=2){
	
		for(j in 2:q1){

			z1 <- numeric(N)

			for(i in 1:N){

				wi <- which(study==i)
				z1[i] <- sum(Z1[wi,j]*n[wi])/sum(n[wi])
	
			}
	
			Z2 <- data.frame(Z2,z1)
		
		}
	
	}
	
	colnames(Z2) <- covname
	

###
	
X1 <- NULL

for(i in 1:p)  X1[[i]] <- t(matrix(numeric(N) + 1))
for(i in 1:p)  rownames(X1[[i]]) <- paste0(i+1,": cons")

###

C <- treats - 1

for(i in C){
		X1[[i]] <- rbind(X1[[i]],t(Z2))
		for(j in 1:q1) rownames(X1[[i]])[j+1] <- paste0(i+1,": ",covname[j])
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

  outcome <- paste0(treats," vs. 1")
  
  R11 <- list("coding"=x$coding,"Covariates"=covname,"Outcomes evaluated the effect modifications"=outcome,"Coefficients"=R5,"Between-studies_SD"=R6,"Between-studies_COR"=R7)
    
  return(R11)
  
}

C1 <- REMLIC(y,S,X1)

return(C1)

}
