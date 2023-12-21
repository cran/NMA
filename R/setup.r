# This is a source code file for several computational functions.

ttrt <- function(x,ref){	

	A1 <- factor(x)
	A2 <- levels(A1)

	w2 <- which(A2==ref)
	
	l <- length(x)
	m <- length(A2)
	
	A3 <- c(A2[w2],setdiff(A2,A2[w2]))

	z <- rep(NA,times=l)

	for(i in 1:m)	z[A1==A3[i]] <- i
	
	A4 <- data.frame(1:m,A3)
	colnames(A4) <- c("code","treatment")

	A5 <- list(coding=A4,treatment=x,code=z)
	return(A5)

}

econtrast1 <- function(X1,X2){		# measure="OR"
  
  N <- dim(X1)[1]
  p <- dim(X1)[2]
  
  y <- matrix(numeric(N*(p-1)),N)
  
  for(i in 2:p){
    
    z1 <- X1[,1]
    z2 <- X2[,1]
    
    x1 <- X1[,i]
    x2 <- X2[,i]
    
    y[,(i-1)] <- log(x1/(x2-x1)) - log(z1/(z2-z1)) 
    
  }
  
  S <- NULL
  
  for(i in 2:(p-1)){
    
    z1 <- X1[,1]
    z2 <- X2[,1]
    
    x1 <- X1[,i]
    x2 <- X2[,i]
    
    s11 <- 1/z1 + 1/(z2-z1) + 1/x1 + 1/(x2-x1)
    
    S <- cbind(S, s11)
    
    for(j in (i+1):p){
      
      w1 <- X1[,j]
      w2 <- X2[,j]
      
      s11 <- 1/z1 + 1/(z2-z1) + (x1 - x1) + (w1 - w1)
      S <- cbind(S, s11)
      
    }
    
  }
  
  x1 <- X1[,p]
  x2 <- X2[,p]
  
  s11 <- 1/z1 + 1/(z2-z1) + 1/x1 + 1/(x2-x1)
  
  S <- cbind(S, s11); colnames(S) <- NULL
  
  mng2 <- list(y=y,S=S)
  
  return(mng2)
  
}


econtrast2 <- function(X1,X2){		# measure="RR"
  
  N <- dim(X1)[1]
  p <- dim(X1)[2]
  
  y <- matrix(numeric(N*(p-1)),N)
  
  for(i in 2:p){
    
    z1 <- X1[,1]
    z2 <- X2[,1]
    
    x1 <- X1[,i]
    x2 <- X2[,i]
    
    y[,(i-1)] <- log(x1/x2) - log(z1/z2)
    
  }
  
  S <- NULL
  
  for(i in 2:(p-1)){
    
    z1 <- X1[,1]
    z2 <- X2[,1]
    
    x1 <- X1[,i]
    x2 <- X2[,i]
    
    s11 <- 1/z1 - 1/z2 + 1/x1 - 1/x2
    
    S <- cbind(S, s11)
    
    for(j in (i+1):p){
      
      w1 <- X1[,j]
      w2 <- X2[,j]
      
      s11 <- 1/z1 - 1/z2 + (x1 - x1) + (w1 - w1)	
      S <- cbind(S, s11)
      
    }
    
  }
  
  x1 <- X1[,p]
  x2 <- X2[,p]
  
  s11 <- 1/z1 - 1/z2 + 1/x1 - 1/x2		
  
  S <- cbind(S, s11); colnames(S) <- NULL
  
  mng2 <- list(y=y,S=S)
  
  return(mng2)
  
}


econtrast3 <- function(X1,X2){		# measure="RD"
  
  N <- dim(X1)[1]
  p <- dim(X1)[2]
  
  y <- matrix(numeric(N*(p-1)),N)
  
  for(i in 2:p){
    
    z1 <- X1[,1]
    z2 <- X2[,1]
    
    x1 <- X1[,i]
    x2 <- X2[,i]
    
    y[,(i-1)] <- x1/x2 - z1/z2
    
  }
  
  S <- NULL
  
  for(i in 2:(p-1)){
    
    z1 <- X1[,1]
    z2 <- X2[,1]
    
    x1 <- X1[,i]
    x2 <- X2[,i]
    
    s11 <- z1*(z2 - z1)/(z2^3) + x1*(x2 - x1)/(x2^3)
    
    S <- cbind(S, s11)
    
    for(j in (i+1):p){
      
      w1 <- X1[,j]
      w2 <- X2[,j]
      
      s11 <- z1*(z2 - z1)/(z2^3) + (x1 - x1) + (w1 - w1)	
      S <- cbind(S, s11)
      
    }
    
  }
  
  x1 <- X1[,p]
  x2 <- X2[,p]
  
  s11 <- z1*(z2 - z1)/(z2^3) + x1*(x2 - x1)/(x2^3)	
  
  S <- cbind(S, s11); colnames(S) <- NULL
  
  mng2 <- list(y=y,S=S)
  
  return(mng2)
  
}


econtrast4 <- function(X1,X2,X3){		# measure="MD"
  
  N <- dim(X1)[1]
  p <- dim(X1)[2]
  
  y <- matrix(numeric(N*(p-1)),N)
  
  for(i in 2:p){
    
    z1 <- X1[,1]
    x1 <- X1[,i]
    
    y[,(i-1)] <- x1 - z1
    
  }
  
  S <- NULL
  
  for(i in 2:(p-1)){
    
    z1 <- X2[,1]
    z2 <- X3[,1]
    
    x1 <- X2[,i]
    x2 <- X3[,i]
    
    s11 <- z1*z1/z2 + x1*x1/x2	
    
    S <- cbind(S, s11)
    
    for(j in (i+1):p){
      
      w1 <- X2[,j]
      w2 <- X3[,j]
      
      s11 <- z1*z1/z2 + (x1 - x1) + (w1 - w1)	
      S <- cbind(S, s11)
      
    }
    
  }
  
  x1 <- X2[,p]
  x2 <- X3[,p]
  
  s11 <- z1*z1/z2 + x1*x1/x2	
  
  S <- cbind(S, s11); colnames(S) <- NULL
  
  mng2 <- list(y=y,S=S)
  
  return(mng2)
  
}


econtrast5 <- function(X1,X2,X3){		# measure="SMD"
  
  N <- dim(X1)[1]
  p <- dim(X1)[2]
  
  y <- matrix(numeric(N*(p-1)),N)
  
  for(i in 2:p){
    
    z1 <- X1[,1]
	z2 <- X2[,1]
	z3 <- X3[,1]
	
    x1 <- X1[,i]
	x2 <- X2[,i]
	x3 <- X3[,i]
	
	#s2 <- (z2 + x2)/2
	s2 <- ((z3-1)*z2*z2 + (x3-1)*x2*x2)/(z3+x3-2)
    
    y[,(i-1)] <- (x1 - z1)/sqrt(s2)
    
  }
  
  S <- NULL
  
  for(i in 2:(p-1)){
    
    z2 <- X3[,1]
    x2 <- X3[,i]
    
    s11 <- 1/z2 + 1/x2		
    
    S <- cbind(S, s11)
    
    for(j in (i+1):p){
      
      w1 <- X2[,j]
      w2 <- X3[,j]
      
      s11 <- 1/z2 + (x2 - x2) + (w2 - w2)	
      S <- cbind(S, s11)
      
    }
    
  }
  
  x2 <- X3[,p]
  
  s11 <- 1/z2 + 1/x2	
  
  S <- cbind(S, s11); colnames(S) <- NULL
  
  mng2 <- list(y=y,S=S)
  
  return(mng2)
  
}


econtrast6 <- function(X1,X2){		# measure="RR"
  
  N <- dim(X1)[1]
  p <- dim(X1)[2]
  
  y <- matrix(numeric(N*(p-1)),N)
  
  for(i in 2:p){
    
    z1 <- X1[,1]
    z2 <- X2[,1]
    
    x1 <- X1[,i]
    x2 <- X2[,i]
    
    y[,(i-1)] <- log(x1/x2) - log(z1/z2)
    
  }
  
  S <- NULL
  
  for(i in 2:(p-1)){
    
    z1 <- X1[,1]
    z2 <- X2[,1]
    
    x1 <- X1[,i]
    x2 <- X2[,i]
    
    s11 <- 1/z1 + 1/x1
    
    S <- cbind(S, s11)
    
    for(j in (i+1):p){
      
      w1 <- X1[,j]
      w2 <- X2[,j]
      
      s11 <- 1/z1 + (x1 - x1) + (w1 - w1)	
      S <- cbind(S, s11)
      
    }
    
  }
  
  x1 <- X1[,p]
  x2 <- X2[,p]
  
  s11 <- 1/z1 + 1/x1
  
  S <- cbind(S, s11); colnames(S) <- NULL
  
  mng2 <- list(y=y,S=S)
  
  return(mng2)
  
}


setup <- function(study,trt,d=NULL,n=NULL,m=NULL,s=NULL,z=NULL,measure,ref,data=NULL){

  data <- data.frame(data)

  sz <- substitute(z)

  if(is.null(sz))  covariate <- NULL

  if(is.null(sz)==FALSE){

	z1 <- deparse(substitute(z))
	
	z2 <- gsub(" ","",z1)
	
	if(substring(z2,1,2)!="c(")  covariate <- z1
	  
	if(substring(z2,1,2)=="c("){

		nz <- nchar(z2)
		z3 <- substring(z2,3,(nz-1))
		covariate <- strsplit(z3,",")[[1]]
		
    }
	
  }
  
  ###

  if(is.null(covariate)) Z <- NULL
  
  if(is.null(covariate)==FALSE){
  
  nc <- length(covariate)
  cn <- colnames(data)
  
  ci <- covariate[1]
  j <- which(cn==ci)
  dj <- data[,j]
  Z <- data.frame(dj)
   
  if(nc>=2){
	for(i in 2:nc){
  
		ci <- covariate[i]
		j <- which(cn==ci)
		dj <- data[,j]
		Z <- data.frame(Z,dj)
  	}
  }
  
  colnames(Z) <- covariate
  
  }
  
  if(measure=="OR"){

	study <- data[, deparse(substitute(study))]
	trt <- data[, deparse(substitute(trt))]
	d <- data[, deparse(substitute(d))]
	n <- data[, deparse(substitute(n))]

	study <- as.numeric(factor(study))

    T1 <- ttrt(trt,ref=ref)
	treat <- T1$code

	N <- max(study)
	p <- max(treat)
	
	L <- length(study)

	X1 <- X2 <- matrix(rep(NA, times=N*p), N)

	for(i in 1:L){
  
		k <- study[i]
		l <- treat[i]
  
		X1[k,l] <- d[i]
		X2[k,l] <- n[i]
  
	}

	p <- p - 1

	for(i in 1:N){
  
	if(length(which(X1[i,]==0))>0){
		X1[i,] <- X1[i,] + 0.5
		X2[i,] <- X2[i,] + 1
	}
  
	if(length(which(X1[i,]==X2[i,]))>0){
		X1[i,] <- X1[i,] + 0.5
		X2[i,] <- X2[i,] + 1
	}
  
	}

	X1[is.na(X1[,1]),1] <- 0.00001		# data autmentation for the reference arms
	X2[is.na(X2[,1]),1] <- 0.0001

	de1 <- econtrast1(X1,X2)				# creating the contrast-based summary statistics
	y <- de1$y
	S <- de1$S
	
	

	n.arm <- numeric(N)

	for(i in 1:N){

		wi <- which(study==i)
		n.arm[i] <- length(wi)
	
	}

	dof <- sum(n.arm) - N - p			# degree-of-freedom for Q-statistic
	
	###

	mng2 <- list(coding=T1$coding,reference=ref,measure=measure,covariate=covariate,N=N,p=p,df=dof,study=study,trt=trt,treat=treat,d=d,n=n,Z=Z,y=y,S=S)
  
	return(mng2)
  
  }
  

  if(measure=="RR"){

	study <- data[, deparse(substitute(study))]
	trt <- data[, deparse(substitute(trt))]
	d <- data[, deparse(substitute(d))]
	n <- data[, deparse(substitute(n))]
	
	study <- as.numeric(factor(study))

    T1 <- ttrt(trt,ref=ref)
	treat <- T1$code

	N <- max(study)
	p <- max(treat)
	
	L <- length(study)

	X1 <- X2 <- matrix(rep(NA, times=N*p), N)

	for(i in 1:L){
  
		k <- study[i]
		l <- treat[i]
  
		X1[k,l] <- d[i]
		X2[k,l] <- n[i]
  
	}

	p <- p - 1

	for(i in 1:N){
  
	if(length(which(X1[i,]==0))>0){
		X1[i,] <- X1[i,] + 0.5
		X2[i,] <- X2[i,] + 1
	}
  
	if(length(which(X1[i,]==X2[i,]))>0){
		X1[i,] <- X1[i,] + 0.5
		X2[i,] <- X2[i,] + 1
	}
  
	}

	X1[is.na(X1[,1]),1] <- 0.00001		# data autmentation for the reference arms
	X2[is.na(X2[,1]),1] <- 0.0001

	de1 <- econtrast2(X1,X2)				# creating the contrast-based summary statistics
	y <- de1$y
	S <- de1$S
	
	

	n.arm <- numeric(N)

	for(i in 1:N){

		wi <- which(study==i)
		n.arm[i] <- length(wi)
	
	}

	dof <- sum(n.arm) - N - p			# degree-of-freedom for Q-statistic
	
	###

	mng2 <- list(coding=T1$coding,reference=ref,measure=measure,covariate=covariate,N=N,p=p,df=dof,study=study,trt=trt,treat=treat,d=d,n=n,Z=Z,y=y,S=S)
  
	return(mng2)
  
  }

    if(measure=="RD"){

	study <- data[, deparse(substitute(study))]
	trt <- data[, deparse(substitute(trt))]
	d <- data[, deparse(substitute(d))]
	n <- data[, deparse(substitute(n))]

	study <- as.numeric(factor(study))

    T1 <- ttrt(trt,ref=ref)
	treat <- T1$code

	N <- max(study)
	p <- max(treat)
	
	L <- length(study)

	X1 <- X2 <- matrix(rep(NA, times=N*p), N)

	for(i in 1:L){
  
		k <- study[i]
		l <- treat[i]
  
		X1[k,l] <- d[i]
		X2[k,l] <- n[i]
  
	}

	p <- p - 1

	for(i in 1:N){
  
	if(length(which(X1[i,]==0))>0){
		X1[i,] <- X1[i,] + 0.5
		X2[i,] <- X2[i,] + 1
	}
  
	if(length(which(X1[i,]==X2[i,]))>0){
		X1[i,] <- X1[i,] + 0.5
		X2[i,] <- X2[i,] + 1
	}
  
	}

	for(i in 1:N){		# data autmentation for the reference arms
	
		P3 <- X1[i,]/X2[i,]
		M3 <- mean(P3,na.rm=TRUE)
		
		if(is.na(X1[i,1])){

			X2[i,1] <- 0.00001
			X1[i,1] <- 0.00001*M3
			
		}
	
	}

	de1 <- econtrast3(X1,X2)				# creating the contrast-based summary statistics
	y <- de1$y
	S <- de1$S
	
	

	n.arm <- numeric(N)

	for(i in 1:N){

		wi <- which(study==i)
		n.arm[i] <- length(wi)
	
	}

	dof <- sum(n.arm) - N - p			# degree-of-freedom for Q-statistic
	
	###

	mng2 <- list(coding=T1$coding,reference=ref,measure=measure,covariate=covariate,N=N,p=p,df=dof,study=study,trt=trt,treat=treat,d=d,n=n,Z=Z,y=y,S=S)
  
	return(mng2)
  
  }
  
  

  if(measure=="MD"){

	study <- data[, deparse(substitute(study))]
	trt <- data[, deparse(substitute(trt))]
	m <- data[, deparse(substitute(m))]
	s <- data[, deparse(substitute(s))]
	n <- data[, deparse(substitute(n))]

	study <- as.numeric(factor(study))

    T1 <- ttrt(trt,ref=ref)
	treat <- T1$code

	N <- max(study)
	p <- max(treat)
	
	L <- length(study)

	X1 <- X2 <- X3 <- matrix(rep(NA, times=N*p), N)

	for(i in 1:L){
  
		k <- study[i]
		l <- treat[i]
  
		X1[k,l] <- m[i]
		X2[k,l] <- s[i]
		X3[k,l] <- n[i]
  
	}

	p <- p - 1

	for(i in 1:N){		# data autmentation for the reference arms
	
		m3 <- mean(X1[i,],na.rm=TRUE)
		s3 <- mean(X2[i,],na.rm=TRUE)
		
		if(is.na(X3[i,1])){

			X1[i,1] <- m3
			X2[i,1] <- s3
			X3[i,1] <- 0.00001
			
		}
	
	}
	
	de1 <- econtrast4(X1,X2,X3)				# creating the contrast-based summary statistics
	y <- de1$y
	S <- de1$S
	
	

	n.arm <- numeric(N)

	for(i in 1:N){

		wi <- which(study==i)
		n.arm[i] <- length(wi)
	
	}

	dof <- sum(n.arm) - N - p			# degree-of-freedom for Q-statistic
	
	###

	mng2 <- list(coding=T1$coding,reference=ref,measure=measure,covariate=covariate,N=N,p=p,df=dof,study=study,trt=trt,treat=treat,m=m,s=s,n=n,Z=Z,y=y,S=S)
  
	return(mng2)
  
  }


  if(measure=="SMD"){

	study <- data[, deparse(substitute(study))]
	trt <- data[, deparse(substitute(trt))]
	m <- data[, deparse(substitute(m))]
	s <- data[, deparse(substitute(s))]
	n <- data[, deparse(substitute(n))]

	study <- as.numeric(factor(study))

    T1 <- ttrt(trt,ref=ref)
	treat <- T1$code

	N <- max(study)
	p <- max(treat)
	
	L <- length(study)

	X1 <- X2 <- X3 <- matrix(rep(NA, times=N*p), N)

	for(i in 1:L){
  
		k <- study[i]
		l <- treat[i]
  
		X1[k,l] <- m[i]
		X2[k,l] <- s[i]
		X3[k,l] <- n[i]
  
	}

	p <- p - 1

	for(i in 1:N){		# data autmentation for the reference arms
	
		m3 <- mean(X1[i,],na.rm=TRUE)
		s3 <- mean(X2[i,],na.rm=TRUE)
		
		if(is.na(X3[i,1])){

			X1[i,1] <- m3
			X2[i,1] <- s3
			X3[i,1] <- 0.00001
			
		}
	
	}
	
	de1 <- econtrast5(X1,X2,X3)				# creating the contrast-based summary statistics
	y <- de1$y
	S <- de1$S
	
	
	###

	n.arm <- numeric(N)

	for(i in 1:N){

		wi <- which(study==i)
		n.arm[i] <- length(wi)
	
	}

	dof <- sum(n.arm) - N - p			# degree-of-freedom for Q-statistic
	
	###

	mng2 <- list(coding=T1$coding,reference=ref,measure=measure,covariate=covariate,N=N,p=p,df=dof,study=study,trt=trt,treat=treat,m=m,s=s,n=n,Z=Z,y=y,S=S)
 
	return(mng2)
  
  }

  if(measure=="HR"){

	study <- data[, deparse(substitute(study))]
	trt <- data[, deparse(substitute(trt))]
	d <- data[, deparse(substitute(d))]
	n <- data[, deparse(substitute(n))]
	
	study <- as.numeric(factor(study))

    T1 <- ttrt(trt,ref=ref)
	treat <- T1$code

	N <- max(study)
	p <- max(treat)
	
	L <- length(study)

	X1 <- X2 <- matrix(rep(NA, times=N*p), N)

	for(i in 1:L){
  
		k <- study[i]
		l <- treat[i]
  
		X1[k,l] <- d[i]
		X2[k,l] <- n[i]
  
	}

	p <- p - 1

	for(i in 1:N){
  
	if(length(which(X1[i,]==0))>0){
		X1[i,] <- X1[i,] + 0.5
		X2[i,] <- X2[i,] + 1
	}
  
	if(length(which(X1[i,]==X2[i,]))>0){
		X1[i,] <- X1[i,] + 0.5
		X2[i,] <- X2[i,] + 1
	}
  
	}

	X1[is.na(X1[,1]),1] <- 0.00001		# data autmentation for the reference arms
	X2[is.na(X2[,1]),1] <- 0.0001

	de1 <- econtrast6(X1,X2)				# creating the contrast-based summary statistics
	y <- de1$y
	S <- de1$S
	
	

	n.arm <- numeric(N)

	for(i in 1:N){

		wi <- which(study==i)
		n.arm[i] <- length(wi)
	
	}

	dof <- sum(n.arm) - N - p			# degree-of-freedom for Q-statistic
	
	###

	mng2 <- list(coding=T1$coding,reference=ref,measure=measure,covariate=covariate,N=N,p=p,df=dof,study=study,trt=trt,treat=treat,d=d,n=n,Z=Z,y=y,S=S)
  
	return(mng2)
  
  }

}

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

fun.I <- function(x){ 
  diag(x)
}
fun.J <- function(x, y = x){ 
  #rep(1, x) %*% t(rep(1, y))
  matrix(1, x, y)
}
fun.e <- function(m, i){
  fun.I(m)[, i]
  #as.numeric((1:m) == i)
}
fun.E <- function(m, i, j){
  fun.e(m, i) %*% t(fun.e(m, j))
}
fun.tilde_P <- function(m){
  (fun.I(m) + fun.J(m)) / 2
}
tr <- function(x){ 
  sum(diag(as.matrix(x)))
  #sum(diag(x))
}
fun.Sum <- function(List){
  rowSums(array(unlist(List), c(dim(List[[1]]), length(List))), dims = 2)
}

ginv2 <- function(A){

	if(is.null(dim(A))==FALSE){
	
		a <- det(A)
	
		if(a>0)	return(solve(A))
		if(a<=0){
			message("Error: The matrix is singular (The inverse matrix does not exist).")
			return(NaN)
		}

	}

	if(is.null(dim(A))==TRUE)	return(solve(A))

}


