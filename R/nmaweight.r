nmaweight <- function(x,digits=3){

	call <- match.call()

	xms <- x$measure	

	if(xms=="OR"||xms=="RR"||xms=="RD"||xms=="HR"||xms=="SPD"){

	study <- x$study
	treat <- x$treat
	n <- x$n

	d <- x$d

	study <- as.numeric(factor(study))
	
	data1 <- data.frame(study,treat,d,n)

    ####

	treat1 <- sort(unique(treat))
	
	N <- length(unique(study))
	p <- max(treat) - 1

	des <- n.arm <- nm <- numeric(N)
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

KR0 <- function(y, S, tau2){
  
  N <- dim(y)[1]
  p <- dim(y)[2]
  
  matrix.indicator <- !(is.na(y))
  
  listlist.y_X_Siinv_tildeP <- lapply(1:N, function(i){
    indicator_i <- matrix.indicator[i, ]
    y_i <- y[i, indicator_i]
    X_i <- fun.I(p)[indicator_i, , drop = FALSE]
    c_i <- length(y_i)
    S_i <- matrix(NA, c_i, c_i)
    S_i[lower.tri(S_i, diag = TRUE)] <- c(na.omit(S[i, ]))
    lower_S_i <- S_i
    S_i <- t(S_i)
    S_i[lower.tri(S_i)] <- lower_S_i[lower.tri(lower_S_i)]
    tilde_P_i <- fun.tilde_P(c_i)
    #Siinv_i <- ginv2(tau2 * tilde_P_i + S_i)
    Siinv_i <- ginv2(tau2 * tilde_P_i + S_i)
    return(list(y_i, X_i, Siinv_i, tilde_P_i))
  })
  
   Sum_xppx <- fun.Sum(lapply(listlist.y_X_Siinv_tildeP, function(z){
    t(z[[2]]) %*% z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[2]]
  }))
  Sum_xppy <- fun.Sum(lapply(listlist.y_X_Siinv_tildeP, function(z){
    t(z[[2]]) %*% z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[1]]
  }))
  Sum_yppy <- fun.Sum(lapply(listlist.y_X_Siinv_tildeP, function(z){
    t(z[[1]]) %*% z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[1]]
  }))
  Sum_xpx <- fun.Sum(lapply(listlist.y_X_Siinv_tildeP, function(z){
    t(z[[2]]) %*% z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[2]]
  }))
  Sum_xpy <- fun.Sum(lapply(listlist.y_X_Siinv_tildeP, function(z){
    t(z[[2]]) %*% z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[1]]
  }))
  Sum_xx <- fun.Sum(lapply(listlist.y_X_Siinv_tildeP, function(z){
    t(z[[2]]) %*% z[[3]] %*% z[[2]]
  }))
  Sum_xy <- fun.Sum(lapply(listlist.y_X_Siinv_tildeP, function(z){
    t(z[[2]]) %*% z[[3]] %*% z[[1]]
  }))
  
  # Phi, P, Q
  #Phi <- ginv2(Sum_xx)
  Phi <- ginv2(Sum_xx)
  P <- - Sum_xpx
  Q <- Sum_xppx
  
  #I_E
  I_E_1 <- sum(unlist(
    lapply(listlist.y_X_Siinv_tildeP, function(z){
      tr(z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[4]])
    })
  ))
  I_E_2 <- tr(2 * (Phi %*% Q) - Phi %*% P %*% Phi %*% P)
  I_E <- (1 / 2) * (I_E_1 - I_E_2)
  
  # I_O
  I_O_1 <- - sum(unlist(
    lapply(listlist.y_X_Siinv_tildeP, function(z){
      sum((z[[3]] %*% z[[4]] %*% z[[3]]) * z[[4]])
    })
  ))
  I_O_2 <- sum(((- Phi %*% P %*% Phi) * P) + (Phi * (2 * Q)))
  I_O_3 <- 2 * c(
    Sum_yppy - t(Sum_xpy) %*% Phi %*% Sum_xpy - 2 * t(Sum_xy) %*% Phi %*% Sum_xppy + 2 * t(Sum_xy) %*% Phi %*% Sum_xpx %*% Phi %*% Sum_xpy + t(Sum_xy) %*% Phi %*% Sum_xppx %*% Phi %*% Sum_xy - t(Sum_xy) %*% Phi %*% Sum_xpx %*% Phi %*% Sum_xpx %*% Phi %*% Sum_xy
  )
  I_O <- (1 / 2) * (I_O_1 + I_O_2 + I_O_3)
  
  # Phi_A_E, Phi_A_O
  Phi_A_E <- Phi + 2 * Phi %*% ((1 / I_E) * (Q - P %*% Phi %*% P)) %*% Phi
  Phi_A_O <- Phi + 2 * Phi %*% ((1 / I_O) * (Q - P %*% Phi %*% P)) %*% Phi
  
  vec.W <- 1 / c(I_E, I_O)
  fun.m_lambda <- function(j){
    L_j <- as.matrix(fun.e(p, j))
    l_j <- 1
    #Theta_j <- L_j %*% ginv2(t(L_j) %*% Phi %*% L_j) %*% t(L_j)
    Theta_j <- L_j %*% ginv2(t(L_j) %*% Phi %*% L_j) %*% t(L_j)
    vec.A1_j <- vec.W * (tr(Theta_j %*% Phi %*% P %*% Phi))^2
    vec.A2_j <- vec.W * tr(Theta_j %*% Phi %*% P %*% Phi %*% Theta_j %*% Phi %*% P %*% Phi)
    vec.g_j <- ((l_j + 1) * vec.A1_j - (l_j + 4) * vec.A2_j) / ((l_j + 2) * vec.A2_j)
    vec.c1_j <- vec.g_j / (3 * l_j + 2 * (1 - vec.g_j))
    vec.c2_j <- (l_j - vec.g_j) / (3 * l_j + 2 * (1 - vec.g_j))
    vec.c3_j <- (l_j + 2 - vec.g_j) / (3 * l_j + 2 * (1 - vec.g_j))
    vec.Estar_j <- 1 / (1 - vec.A2_j / l_j)
    vec.B_j <- (vec.A1_j + 6 * vec.A2_j) / (2 * l_j)
    vec.Vstar_j <- (2 / l_j) * (1 + vec.c1_j * vec.B_j) / ((1 - vec.c2_j * vec.B_j)^2 * (1 - vec.c3_j * vec.B_j))
    vec.rho_j <- vec.Vstar_j / (2 * vec.Estar_j^2)
    vec.m_j <- 4 + (l_j + 2) / (l_j * vec.rho_j - 1)
    vec.lambda_j <- vec.m_j / (vec.Estar_j * (vec.m_j - 2))
    return(c(m_E = vec.m_j[1], lambda_E = vec.lambda_j[1], m_O = vec.m_j[2], lambda_O = vec.lambda_j[2]))
  }
  matrix.m_lambda <- sapply(1:p, fun.m_lambda)
  
  SE_E <- sqrt(pmax(0, diag(Phi_A_E)))
  SE_O <- sqrt(pmax(0, diag(Phi_A_O)))
  df_E <- matrix.m_lambda[1, ]
  lambda_E <- matrix.m_lambda[2, ]
  df_O <- matrix.m_lambda[3, ]
  lambda_O <- matrix.m_lambda[4, ]
  
   if(any(diag(Phi) < 0)){
    warning("At least one of the diagonal elements of the matrix `Phi` is negative. ")
  }
  if(any(diag(Phi_A_E) < 0)){
    warning("At least one of the diagonal elements of the matrix `Phi_A_E` is negative. ")
  }
  if(any(diag(Phi_A_O) < 0)){
    warning("At least one of the diagonal elements of the matrix `Phi_A_O` is negative. ")
  }
  if(any(df_E < 0)){
    warning("At least one element of the vector `df_E` is negative. ")
  }
  if(any(df_O < 0)){
    warning("At least one element of the vector `df_O` is negative. ")
  }
  
  #
  return(
    list(Expected = list(SE = SE_E, df = df_E, lambda = lambda_E), 
         Observed = list(SE = SE_O, df = df_O, lambda = lambda_O),
		 V=Phi_A_E)
  )
  
}

REML <- function(y,S,maxitr=200){
  
  N <- dim(y)[1]
  p <- dim(y)[2]
  
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
      B2 <- ginv2(Gi + Si)
      
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
      B2 <- ginv2(Gi + Si)
      
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
      
      Wi <- ginv2(Gi + Si)
      
      A1 <- A1 + ivec(yi %*% Wi, wi, p)
      A2 <- A2 + imat(Wi, wi, p)
      
    }
    
    mu <- A1 %*% ginv2(A2)
    g1 <- optimize(LL1, lower = 0, upper = 5)$minimum
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
  
  R5 <- cbind(R1,R2,R3,R4); colnames(R5) <- c("Coef.","SE","95%CL","95%CU")
  
  R6 <- sqrt(g1)
  R7 <- g2/g1
  
  R8 <- list("Coefficients"=R5,"Between-studies_SD"=R6,"Between-studies_COR"=R7)
  
  #return(R8)
  return(
    list("Coefficients"=R5,"Between-studies_SD"=R6,mu=mu,V=V.mu)
  )
  
}

####

W1 <- W2 <- W3 <- W4 <- W5 <- cl6 <- NULL

for(k in 1:p){

 T1 <- ttrt(treat, ref=treat1[k])	
 data1$treat1 <- T1$code

 edat <- setup(study=study,trt=treat1,d=d,n=n,measure=xms,ref=1,data=data1)

 y <- edat$y
 S <- edat$S
	
 ###

 result_reml <- REML(y,S)
 tau <- result_reml[[2]]

 ###
 
 G <- gmat(tau*tau,.5*tau*tau,p)

 A2 <- matrix(numeric(p*p),p)
 
 for(i in 1:N){
      
  yi <- as.vector(y[i,])
  wi <- which(is.na(yi)==FALSE)
  pl <- length(wi)
      
  Si <- vmat(S[i,], p)
      
  yi <- yi[wi]
  Si <- pmat(Si,wi)
  Gi <- pmat(G,wi)

  Wi <- ginv2(Gi + Si)
  
  A2 <- A2 + imat(Wi, wi, p)
      
 }

 W <- ginv2(A2)
 
 WM <- matrix(rep(0,times=N*p),N)

 for(i in 1:N){
      
  yi <- as.vector(y[i,])
  wi <- which(is.na(yi)==FALSE)
  pl <- length(wi)
      
  Si <- vmat(S[i,], p)
      
  yi <- yi[wi]
  Si <- pmat(Si,wi)
  Gi <- pmat(G,wi)

  Wi <- ginv2(Gi + Si)
  
  wti <- imat(Wi, wi, p) %*% W   # total weight
  
  WM[i,wi] <- diag(wti)[wi]
        
 }
 
 WD <- matrix(rep(0,times=N*p),N)

 for(i in 1:N){
      
  yi <- as.vector(y[i,])
  wi <- which(is.na(yi)==FALSE)
  pl <- length(wi)
      
  Si <- vmat(S[i,], p)
      
  yi <- yi[wi]
  Si <- pmat(Si,wi)
  Gi <- pmat(G,wi)

  if(pl==1) Wi <- 1/(Gi + Si)
  if(pl>=2) Wi <- ginv2( diag(diag(Gi + Si)) )
  
  wti <- imat(Wi, wi, p) %*% W   # direct weight
  
  WD[i,wi] <- diag(wti)[wi]
        
 }
 
 WI <- WM - WD  # indirect weight

 ###
 
 ifd <- round(apply(WD,2,sum),3)
 ifi <- round(apply(WI,2,sum),3)
 
 ###

 WM <- round(WM[,k:p],3)
 WD <- round(WD[,k:p],3)
 WI <- round(WI[,k:p],3)

 if(k==p){

  WM <- matrix(round(WM,3))
  WD <- matrix(round(WD,3))
  WI <- matrix(round(WI,3))
 
 }
 
 cl <- paste0(k,"-",(k+1):(p+1))
 colnames(WM) <- colnames(WD) <- colnames(WI) <- cl

 ###
 
 W1 <- cbind(W1,WD)
 W2 <- cbind(W2,WI)
 W3 <- cbind(W3,WM)
 
 W4 <- c(W4,ifd[k:p])
 W5 <- c(W5,ifi[k:p])
 
 cl6 <- c(cl6,cl)
 
}

 spi <- data.frame(1:N,nm,des)	
 colnames(spi) <- c("study","n","design")

 colnames(W1) <- cl6
 colnames(W2) <- cl6
 colnames(W3) <- cl6
 
 W2[W2<0] <- 0

 WD <- cbind(spi,W1)
 WI <- cbind(spi,W2)
 WM <- cbind(spi,W3)

 ###
 
 W6 <- rbind(W4,W5)
 colnames(W6) <- cl6
 rownames(W6) <- c("direct","indirect")

 ###
 
 es <- factor(rep(cl6,each=N))
 id <- factor(rep(1:N,times=length(cl6)),levels=N:1)

 ##
 
 #Q3 <- data.frame(id,es,c(W3))
 #colnames(Q3) <- c("study","contrast","weight")

 weight <- c(W3)
 Q3 <- data.frame(id,es,weight)
 
 #ghm3 <- ggplot(Q3, aes(x = contrast, y = study, fill = weight))
 ghm3 <- ggplot(Q3, aes(x = es, y = id, fill = weight))
 ghm3 <- ghm3 + geom_tile()
 ghm3 <- ghm3 + theme_bw()
 ghm3 <- ghm3 + theme(plot.background = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.background = element_blank(),
                   axis.line = element_blank(),
                   axis.ticks = element_blank(),
                   strip.background = element_rect(fill = "white", colour = "white"),
                   axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
 ghm3 <- ghm3 + geom_text(aes(label = round(weight,2), colour="grey"), show.legend = FALSE) + scale_fill_gradient2(low = "white", high = "blue", midpoint = 0.005)
 ghm3 <- ghm3 + xlab("Contrast") + ylab("Study")
 #ghm3   # Making a heatmap by ggplot2

 R1 <- list("Heatmap is created by ggplot2."=ghm3, "coding"=x$coding, "Contribution of direct and indirect information"=W6, "Contribution weights: Direct comparison"=WD, "Contribution weights: Indirect comparison (BoS)"=WI, "Contribution weights: Overall evidence"=WM,digits=digits,call=call)

 # message("Contribution weight matrices for the consistency model:")

 class(R1) <- "nmaweight"
 return(R1)

 }
 
 	if(xms=="MD"||xms=="SMD"){

	study <- x$study
	treat <- x$treat
	n <- x$n

	m <- x$m
	s <- x$s

	study <- as.numeric(factor(study))
	
	data1 <- data.frame(study,treat,m,s,n)

    ####

	treat1 <- sort(unique(treat))
	
	N <- length(unique(study))
	p <- max(treat) - 1

	des <- n.arm <- nm <- numeric(N)
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

KR0 <- function(y, S, tau2){
  
  N <- dim(y)[1]
  p <- dim(y)[2]
  
  matrix.indicator <- !(is.na(y))
  
  listlist.y_X_Siinv_tildeP <- lapply(1:N, function(i){
    indicator_i <- matrix.indicator[i, ]
    y_i <- y[i, indicator_i]
    X_i <- fun.I(p)[indicator_i, , drop = FALSE]
    c_i <- length(y_i)
    S_i <- matrix(NA, c_i, c_i)
    S_i[lower.tri(S_i, diag = TRUE)] <- c(na.omit(S[i, ]))
    lower_S_i <- S_i
    S_i <- t(S_i)
    S_i[lower.tri(S_i)] <- lower_S_i[lower.tri(lower_S_i)]
    tilde_P_i <- fun.tilde_P(c_i)
    #Siinv_i <- ginv2(tau2 * tilde_P_i + S_i)
    Siinv_i <- ginv2(tau2 * tilde_P_i + S_i)
    return(list(y_i, X_i, Siinv_i, tilde_P_i))
  })
  
   Sum_xppx <- fun.Sum(lapply(listlist.y_X_Siinv_tildeP, function(z){
    t(z[[2]]) %*% z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[2]]
  }))
  Sum_xppy <- fun.Sum(lapply(listlist.y_X_Siinv_tildeP, function(z){
    t(z[[2]]) %*% z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[1]]
  }))
  Sum_yppy <- fun.Sum(lapply(listlist.y_X_Siinv_tildeP, function(z){
    t(z[[1]]) %*% z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[1]]
  }))
  Sum_xpx <- fun.Sum(lapply(listlist.y_X_Siinv_tildeP, function(z){
    t(z[[2]]) %*% z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[2]]
  }))
  Sum_xpy <- fun.Sum(lapply(listlist.y_X_Siinv_tildeP, function(z){
    t(z[[2]]) %*% z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[1]]
  }))
  Sum_xx <- fun.Sum(lapply(listlist.y_X_Siinv_tildeP, function(z){
    t(z[[2]]) %*% z[[3]] %*% z[[2]]
  }))
  Sum_xy <- fun.Sum(lapply(listlist.y_X_Siinv_tildeP, function(z){
    t(z[[2]]) %*% z[[3]] %*% z[[1]]
  }))
  
  # Phi, P, Q
  #Phi <- ginv2(Sum_xx)
  Phi <- ginv2(Sum_xx)
  P <- - Sum_xpx
  Q <- Sum_xppx
  
  #I_E
  I_E_1 <- sum(unlist(
    lapply(listlist.y_X_Siinv_tildeP, function(z){
      tr(z[[3]] %*% z[[4]] %*% z[[3]] %*% z[[4]])
    })
  ))
  I_E_2 <- tr(2 * (Phi %*% Q) - Phi %*% P %*% Phi %*% P)
  I_E <- (1 / 2) * (I_E_1 - I_E_2)
  
  # I_O
  I_O_1 <- - sum(unlist(
    lapply(listlist.y_X_Siinv_tildeP, function(z){
      sum((z[[3]] %*% z[[4]] %*% z[[3]]) * z[[4]])
    })
  ))
  I_O_2 <- sum(((- Phi %*% P %*% Phi) * P) + (Phi * (2 * Q)))
  I_O_3 <- 2 * c(
    Sum_yppy - t(Sum_xpy) %*% Phi %*% Sum_xpy - 2 * t(Sum_xy) %*% Phi %*% Sum_xppy + 2 * t(Sum_xy) %*% Phi %*% Sum_xpx %*% Phi %*% Sum_xpy + t(Sum_xy) %*% Phi %*% Sum_xppx %*% Phi %*% Sum_xy - t(Sum_xy) %*% Phi %*% Sum_xpx %*% Phi %*% Sum_xpx %*% Phi %*% Sum_xy
  )
  I_O <- (1 / 2) * (I_O_1 + I_O_2 + I_O_3)
  
  # Phi_A_E, Phi_A_O
  Phi_A_E <- Phi + 2 * Phi %*% ((1 / I_E) * (Q - P %*% Phi %*% P)) %*% Phi
  Phi_A_O <- Phi + 2 * Phi %*% ((1 / I_O) * (Q - P %*% Phi %*% P)) %*% Phi
  
  vec.W <- 1 / c(I_E, I_O)
  fun.m_lambda <- function(j){
    L_j <- as.matrix(fun.e(p, j))
    l_j <- 1
    #Theta_j <- L_j %*% ginv2(t(L_j) %*% Phi %*% L_j) %*% t(L_j)
    Theta_j <- L_j %*% ginv2(t(L_j) %*% Phi %*% L_j) %*% t(L_j)
    vec.A1_j <- vec.W * (tr(Theta_j %*% Phi %*% P %*% Phi))^2
    vec.A2_j <- vec.W * tr(Theta_j %*% Phi %*% P %*% Phi %*% Theta_j %*% Phi %*% P %*% Phi)
    vec.g_j <- ((l_j + 1) * vec.A1_j - (l_j + 4) * vec.A2_j) / ((l_j + 2) * vec.A2_j)
    vec.c1_j <- vec.g_j / (3 * l_j + 2 * (1 - vec.g_j))
    vec.c2_j <- (l_j - vec.g_j) / (3 * l_j + 2 * (1 - vec.g_j))
    vec.c3_j <- (l_j + 2 - vec.g_j) / (3 * l_j + 2 * (1 - vec.g_j))
    vec.Estar_j <- 1 / (1 - vec.A2_j / l_j)
    vec.B_j <- (vec.A1_j + 6 * vec.A2_j) / (2 * l_j)
    vec.Vstar_j <- (2 / l_j) * (1 + vec.c1_j * vec.B_j) / ((1 - vec.c2_j * vec.B_j)^2 * (1 - vec.c3_j * vec.B_j))
    vec.rho_j <- vec.Vstar_j / (2 * vec.Estar_j^2)
    vec.m_j <- 4 + (l_j + 2) / (l_j * vec.rho_j - 1)
    vec.lambda_j <- vec.m_j / (vec.Estar_j * (vec.m_j - 2))
    return(c(m_E = vec.m_j[1], lambda_E = vec.lambda_j[1], m_O = vec.m_j[2], lambda_O = vec.lambda_j[2]))
  }
  matrix.m_lambda <- sapply(1:p, fun.m_lambda)
  
  SE_E <- sqrt(pmax(0, diag(Phi_A_E)))
  SE_O <- sqrt(pmax(0, diag(Phi_A_O)))
  df_E <- matrix.m_lambda[1, ]
  lambda_E <- matrix.m_lambda[2, ]
  df_O <- matrix.m_lambda[3, ]
  lambda_O <- matrix.m_lambda[4, ]
  
   if(any(diag(Phi) < 0)){
    warning("At least one of the diagonal elements of the matrix `Phi` is negative. ")
  }
  if(any(diag(Phi_A_E) < 0)){
    warning("At least one of the diagonal elements of the matrix `Phi_A_E` is negative. ")
  }
  if(any(diag(Phi_A_O) < 0)){
    warning("At least one of the diagonal elements of the matrix `Phi_A_O` is negative. ")
  }
  if(any(df_E < 0)){
    warning("At least one element of the vector `df_E` is negative. ")
  }
  if(any(df_O < 0)){
    warning("At least one element of the vector `df_O` is negative. ")
  }
  
  #
  return(
    list(Expected = list(SE = SE_E, df = df_E, lambda = lambda_E), 
         Observed = list(SE = SE_O, df = df_O, lambda = lambda_O),
		 V=Phi_A_E)
  )
  
}

REML <- function(y,S,maxitr=200){
  
  N <- dim(y)[1]
  p <- dim(y)[2]
  
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
      B2 <- ginv2(Gi + Si)
      
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
      B2 <- ginv2(Gi + Si)
      
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
      
      Wi <- ginv2(Gi + Si)
      
      A1 <- A1 + ivec(yi %*% Wi, wi, p)
      A2 <- A2 + imat(Wi, wi, p)
      
    }
    
    mu <- A1 %*% ginv2(A2)
    g1 <- optimize(LL1, lower = 0, upper = 5)$minimum
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
  
  R5 <- cbind(R1,R2,R3,R4); colnames(R5) <- c("Coef.","SE","95%CL","95%CU")
  
  R6 <- sqrt(g1)
  R7 <- g2/g1
  
  R8 <- list("Coefficients"=R5,"Between-studies_SD"=R6,"Between-studies_COR"=R7)
  
  #return(R8)
  return(
    list("Coefficients"=R5,"Between-studies_SD"=R6,mu=mu,V=V.mu)
  )
  
}

####

W1 <- W2 <- W3 <- W4 <- W5 <- cl6 <- NULL

for(k in 1:p){

 T1 <- ttrt(treat, ref=treat1[k])	
 data1$treat1 <- T1$code

 edat <- setup(study=study,trt=treat1,m=m,s=s,n=n,measure=xms,ref=1,data=data1)

 y <- edat$y
 S <- edat$S
	
 ###

 result_reml <- REML(y,S)
 tau <- result_reml[[2]]

 ###
 
 G <- gmat(tau*tau,.5*tau*tau,p)

 A2 <- matrix(numeric(p*p),p)
 
 for(i in 1:N){
      
  yi <- as.vector(y[i,])
  wi <- which(is.na(yi)==FALSE)
  pl <- length(wi)
      
  Si <- vmat(S[i,], p)
      
  yi <- yi[wi]
  Si <- pmat(Si,wi)
  Gi <- pmat(G,wi)

  Wi <- ginv2(Gi + Si)
  
  A2 <- A2 + imat(Wi, wi, p)
      
 }

 W <- ginv2(A2)
 
 WM <- matrix(rep(0,times=N*p),N)

 for(i in 1:N){
      
  yi <- as.vector(y[i,])
  wi <- which(is.na(yi)==FALSE)
  pl <- length(wi)
      
  Si <- vmat(S[i,], p)
      
  yi <- yi[wi]
  Si <- pmat(Si,wi)
  Gi <- pmat(G,wi)

  Wi <- ginv2(Gi + Si)
  
  wti <- imat(Wi, wi, p) %*% W   # total weight
  
  WM[i,wi] <- diag(wti)[wi]
        
 }
 
 WD <- matrix(rep(0,times=N*p),N)

 for(i in 1:N){
      
  yi <- as.vector(y[i,])
  wi <- which(is.na(yi)==FALSE)
  pl <- length(wi)
      
  Si <- vmat(S[i,], p)
      
  yi <- yi[wi]
  Si <- pmat(Si,wi)
  Gi <- pmat(G,wi)

  if(pl==1) Wi <- 1/(Gi + Si)
  if(pl>=2) Wi <- ginv2( diag(diag(Gi + Si)) )
  
  wti <- imat(Wi, wi, p) %*% W   # direct weight
  
  WD[i,wi] <- diag(wti)[wi]
        
 }
 
 WI <- WM - WD  # indirect weight

 ###
 
 ifd <- round(apply(WD,2,sum),3)
 ifi <- round(apply(WI,2,sum),3)
 
 ###

 WM <- round(WM[,k:p],3)
 WD <- round(WD[,k:p],3)
 WI <- round(WI[,k:p],3)

 if(k==p){

  WM <- matrix(round(WM,3))
  WD <- matrix(round(WD,3))
  WI <- matrix(round(WI,3))
 
 }
 
 cl <- paste0(k,"-",(k+1):(p+1))
 colnames(WM) <- colnames(WD) <- colnames(WI) <- cl

 ###
 
 W1 <- cbind(W1,WD)
 W2 <- cbind(W2,WI)
 W3 <- cbind(W3,WM)
 
 W4 <- c(W4,ifd[k:p])
 W5 <- c(W5,ifi[k:p])
 
 cl6 <- c(cl6,cl)
 
}

 spi <- data.frame(1:N,nm,des)	
 colnames(spi) <- c("study","n","design")
 
 W2[W2<0] <- 0

 colnames(W1) <- cl6
 colnames(W2) <- cl6
 colnames(W3) <- cl6

 WD <- cbind(spi,W1)
 WI <- cbind(spi,W2)
 WM <- cbind(spi,W3)

 ###
 
 W6 <- rbind(W4,W5)
 colnames(W6) <- cl6
 rownames(W6) <- c("direct","indirect")

 ###
 
 es <- factor(rep(cl6,each=N))
 id <- factor(rep(1:N,times=length(cl6)),levels=N:1)

 ##
 
 #Q3 <- data.frame(id,es,c(W3))
 #colnames(Q3) <- c("study","contrast","weight")

 weight <- c(W3)
 Q3 <- data.frame(id,es,weight)
 
 #ghm3 <- ggplot(Q3, aes(x = contrast, y = study, fill = weight))
 ghm3 <- ggplot(Q3, aes(x = es, y = id, fill = weight))
 ghm3 <- ghm3 + geom_tile()
 ghm3 <- ghm3 + theme_bw()
 ghm3 <- ghm3 + theme(plot.background = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.background = element_blank(),
                   axis.line = element_blank(),
                   axis.ticks = element_blank(),
                   strip.background = element_rect(fill = "white", colour = "white"),
                   axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
 ghm3 <- ghm3 + geom_text(aes(label = round(weight,2), colour="grey"), show.legend = FALSE) + scale_fill_gradient2(low = "white", high = "blue", midpoint = 0.005)
 ghm3 <- ghm3 + xlab("Contrast") + ylab("Study")
 #ghm3   # Making a heatmap by ggplot2

 R1 <- list("coding"=x$coding, "Contribution of direct and indirect information"=W6, "Contribution weights: Direct comparison"=WD, "Contribution weights: Indirect comparison (BoS)"=WI, "Contribution weights: Overall evidence"=WM,digits=digits,call=call,"Heatmap is created by ggplot2."=ghm3)

 # message("Contribution weight matrices for the consistency model:")

 class(R1) <- "nmaweight"
 return(R1)

 }
  
}



print.nmaweight <- function(x, digits = x$digits, ...) {

  cat("Call:\n")
  print(x$call,row.names=FALSE)
  cat("\n")
  
  cat("Coding:\n", sep = "")
  print(x$coding,row.names=FALSE)
  cat("\n")

  cat("Contribution of direct and indirect information: ", sep = "")
  cat("\n")
  A <- x[[3]]
  print(A,row.names=FALSE)
  cat("\n")

  cat("Contribution weights: Direct comparison: ", sep = "")
  cat("\n")
  A <- x[[4]]
  print(A,row.names=FALSE)
  cat("\n")

  cat("Contribution weights: Indirect comparison (BoS): ", sep = "")
  cat("\n")
  A <- x[[5]]
  print(A,row.names=FALSE)
  cat("\n")

  cat("Contribution weights: Overall evidence: ", sep = "")
  cat("\n")
  A <- x[[6]]
  print(A,row.names=FALSE)
  cat("\n")

  print(x[[1]])
  
  invisible(x)
  
}
