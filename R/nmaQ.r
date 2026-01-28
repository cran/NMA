nmaQ <- function(x,digits=3){

	call <- match.call()

	study <- x$study
	treat <- x$treat
	n <- x$n
	
	y <- x$y
	S <- x$S
	dof <- x$df
	
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
		nm[i] <- sum(n[wi])
	
	}
	
	des0 <- sort(unique(des))

	L <- length(des0)
	N.des <- n.des <- arm1 <- Q1 <- df1 <- rep(NA,times=L)
	
	for(i in 1:L){
	
		desi <- des0[i]
		wi <- which(des==desi)

		N.des[i] <- length(wi)
		n.des[i] <- sum(nm[wi])
		arm1[i] <- n.arm[wi][1]
		
		if(length(wi)==1){
			Q1[i] <- 0
			df1[i] <- sum(n.arm[wi]) - length(wi)
		}
		
		if(length(wi)>=2){
			yi <- y[wi,]
			w1 <- which(is.na(yi[1,])==FALSE)
			yi <- yi[,w1]
      
			Si <- S[wi,]
			w2 <- which(is.na(Si[1,])==FALSE)
			Si <- Si[,w2]

			if(length(w1)==1){
				rmai <- rma(yi=yi,vi=Si)
				Q1[i] <- rmai$QE
				df1[i] <- length(yi) - 1
			}
			
			if(length(w1)>=2){
				Q1[i] <- FEM1(yi,Si)$Q
			}

		}
		
	}

	df1 <- N.des * (arm1 - 1) - (arm1 - 1)
	P1 <- 1 - pchisq(Q1,df=df1)
	
	R1 <- data.frame(des0,N.des,n.des,Q1,df1,P1)
	colnames(R1) <- c("design","N","n","Q","df","P-value")
	
	Q2 <- sum(Q1)
	df2 <- sum(df1)
	P2 <- 1 - pchisq(Q2,df=df2)

	R2 <- paste0("Q(df = ",round(df2),") = ",round(Q2,4),", p-value = ",round(P2,4))

	fem0 <- FEM1(y,S)
	Qs <- fem0$Q
	Ps <- 1 - pchisq(Qs,df=dof)
	R3 <- paste0("Q(df = ",round(dof),") = ",round(Qs,4),", p-value = ",round(Ps,4))
	
	Q3 <- Qs - Q2
	df3 <- dof - df2
	P3 <- 1 - pchisq(Q3,df=df3)

	R4 <- paste0("Q(df = ",round(df3),") = ",round(Q3,4),", p-value = ",round(P3,4))

	R6 <- data.frame(c(Q2,Q3,Qs),c(df2,df3,dof),c(P2,P3,Ps))
	colnames(R6) <- c("Q","df","P-value")
	rownames(R6) <- c("Within designs","Between designs","Total")

	R5 <- list("coding"=x$coding,"number of studies"=N,"Within designs (individual designs)"=R1,"Q-statistics"=R6,digits=digits,call=call)

    class(R5) <- "nmaQ"  
	return(R5)
	
}
	
print.nmaQ <- function(x, digits = x$digits, ...) {

  cat("Call:\n")
  print(x$call,row.names=FALSE)
  cat("\n")
  
  cat("Coding:\n", sep = "")
  print(x$coding,row.names=FALSE)
  cat("\n")

  cat("Number of studies: ", sep = "")
  cat(x[[2]])
  cat("\n")
  cat("\n")
  
  cat("Within designs (individual designs): ", sep = "")
  cat("\n")
  A <- x[[3]]
  ##
  design <- A[,1]
  N <- round(A[,2])
  n <- round(A[,3])
  Q <- round(A[,4],digits)
  df <- round(A[,5])
  pval <- round(A[,6],digits)
  TAB <- cbind(
    "N" = N,
	"n" = n,
	"Q" = Q,
	"df" = df,
    "Pr(>Q)"  = pval
  )
  rownames(TAB) <- design
  print(TAB)
  cat("\n")

  cat("Q-statistics: ", sep = "")
  cat("\n")
  A <- x[[4]]
  ##
  Q <- round(A[,1],digits)
  df <- round(A[,2])
  pval <- round(A[,3],digits)
  TAB <- cbind(
    "Q" = Q,
    "df" = df,
    "Pr(>Q)"  = pval
  )
  rownames(TAB) <- rownames(A)
  print(TAB)
  cat("\n")

  invisible(x)
  
}

###

FEM1 <- function(y,S){
  
  N <- dim(y)[1]
  p <- dim(y)[2]
  
  MY <- max(y,na.rm=TRUE) - min(y,na.rm=TRUE)

  mu <- rnorm(p)	# initial values
  
  A1 <- numeric(p)
  A2 <- matrix(numeric(p*p),p)
    
    for(i in 1:N){
      
      yi <- as.vector(y[i,])
      wi <- which(is.na(yi)==FALSE)
      pl <- length(wi)
      
      Si <- vmat(S[i,], p)
      
      yi <- yi[wi]
      Si <- pmat(Si,wi)
      
      Wi <- solve(Si)
      
      A1 <- A1 + ivec(yi %*% Wi, wi, p)
      A2 <- A2 + imat(Wi, wi, p)
      
    }
    
  mu <- A1 %*% solve(A2)
   
  V.mu <- solve(A2)
  
####  

  A3 <- 0

    for(i in 1:N){
      
      yi <- as.vector(y[i,])
      wi <- which(is.na(yi)==FALSE)
      pl <- length(wi)
      
      Si <- vmat(S[i,], p)
      
      yi <- yi[wi]
      Si <- pmat(Si,wi)
      
      Wi <- solve(Si)
	  
	  mui <- mu[wi]
	  
	  A3 <- A3 + t(yi - mui) %*% Wi %*% (yi - mui)
      
    }
    
####

  SE <- sqrt(diag(V.mu))
  
  R1 <- as.vector(mu)
  R2 <- as.vector(SE)
  R3 <- as.vector(mu - qnorm(.975)*SE)
  R4 <- as.vector(mu + qnorm(.975)*SE)
  P4 <- 2*(1-pnorm(abs(R1)/R2))
  
  R5 <- cbind(R1,R2,R3,R4,P4); colnames(R5) <- c("Coef.","SE","95%CL","95%CU","P-value")
  
  #return(R8)
  return(
    list("Coefficients"=R5,"CF"=V.mu,"Q"=A3)
  )
  
}


REML1 <- function(y,S,maxitr=200){
  
  N <- dim(y)[1]
  p <- dim(y)[2]
  
  MY <- max(y,na.rm=TRUE) - min(y,na.rm=TRUE)

  mu <- rnorm(p)	# initial values
  g1 <- 0.2
  g2 <- 0.1
  
  Qc0 <- c(mu,g1,g2)
  
  LL1 <- function(g){
    
    #G <- gmat(g,g2,p)
    G <- gmat(g,(g/2),p)
    
    ll1 <- 0; XWX <- gmat(0,0,p)
    
    for(i in 1:N){
      
      yi <- as.vector(y[i,])
      wi <- which(is.na(yi)==FALSE)
      pl <- length(wi)
      
      Si <- vmat(S[i,], p)
      
      yi <- yi[wi]
      Si <- pmat(Si,wi)
      mui <- mu[wi]
      Gi <- pmat(G,wi)
      
      B1 <- (yi - mui)
      B2 <- solve(Gi + Si)
      
      A1 <- log(det(Gi + Si))
      A2 <- t(B1) %*% B2 %*% B1
      A3 <- pl * log(2*pi)
      
      XWX <- XWX + imat(B2,wi, p)
      
      ll1 <- ll1 + A1 + A2 + A3
      
    }
    
    ll2 <- ll1 + log(det(XWX)) - p*log(2*pi)
    
    return(ll2)
    
  }
  
  LL2 <- function(g){
    
    G <- gmat(g1,g,p)
    
    ll1 <- 0; XWX <- gmat(0,0,p)
    
    for(i in 1:N){
      
      yi <- as.vector(y[i,])
      wi <- which(is.na(yi)==FALSE)
      pl <- length(wi)
      
      Si <- vmat(S[i,], p)
      
      yi <- yi[wi]
      Si <- pmat(Si,wi)
      mui <- mu[wi]
      Gi <- pmat(G,wi)
      
      B1 <- (yi - mui)
      B2 <- solve(Gi + Si)
      
      A1 <- log(det(Gi + Si))
      A2 <- t(B1) %*% B2 %*% B1
      A3 <- pl * log(2*pi)
      
      XWX <- XWX + imat(B2,wi, p)
      
      ll1 <- ll1 + A1 + A2 + A3
      
    }
    
    ll2 <- ll1 + log(det(XWX)) - p*log(2*pi)
    
    return(ll2)
    
  }
  
  for(itr in 1:maxitr){
    
    A1 <- numeric(p)
    A2 <- matrix(numeric(p*p),p)
    
    G <- gmat(g1,g2,p)
    
    for(i in 1:N){
      
      yi <- as.vector(y[i,])
      wi <- which(is.na(yi)==FALSE)
      pl <- length(wi)
      
      Si <- vmat(S[i,], p)
      
      yi <- yi[wi]
      Si <- pmat(Si,wi)
      Gi <- pmat(G,wi)
      
      Wi <- solve(Gi + Si)
      
      A1 <- A1 + ivec(yi %*% Wi, wi, p)
      A2 <- A2 + imat(Wi, wi, p)
      
    }
    
    mu <- A1 %*% solve(A2)
    g1 <- optimize(LL1, lower = 0, upper = MY)$minimum
    g2 <- 0.5*g1
    
    V.mu <- solve(A2)
    
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
  
  R5 <- cbind(R1,R2,R3,R4,P4); colnames(R5) <- c("Coef.","SE","95%CL","95%CU","P-value")
  
  R6 <- sqrt(g1)
  R7 <- g2/g1
  
  #R8 <- list("Coefficients"=R5,"Between-studies_SD"=R6,"Between-studies_COR"=R7)
  
  #return(R8)
  return(
    list("Coefficients"=R5,"Between-studies_SD"=R6,"CR"=V.mu)
  )
  
}

###

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
