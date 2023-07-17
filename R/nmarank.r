nmarank <- function(x, B=20000, method="NH",ascending=TRUE){

y <- x$y
S <- x$S
treat <- x$coding[,2]


MY <- max(y,na.rm=TRUE) - min(y,na.rm=TRUE)

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
  
  R5 <- cbind(R1,R2,R3,R4); colnames(R5) <- c("Coef.","SE","95%CL","95%CU")
  
  R6 <- sqrt(g1)
  R7 <- g2/g1
  
  R8 <- list("Coefficients"=R5,"Between-studies_SD"=R6,"Between-studies_COR"=R7)
  
  #return(R8)
  return(
    list("Coefficients"=R5,"Between-studies_SD"=R6,mu=mu,V=V.mu)
  )
  
}

if(method=="NH"){

 result_reml <- REML(y,S)
 tau_reml <- result_reml[[2]]
 R1 <- KR0(y = y, S = S, tau2 = tau_reml^2)
  
 mu <- result_reml$mu
 V <- R1$V
  
 if(ascending==TRUE){

  x1 <- mvrnorm(B, mu, V)
  x1 <- cbind(rep(0,times=B),x1)

  x2 <- apply(x1,1,rank)

  p <- dim(x2)[1]
  
  R1 <- matrix(rep(NA,times=p*p),p)

  for(i in 1:p){
	for(j in 1:p){
		R1[i,j] <- sum(x2[i,]==j)/B
	}
  }

  if(is.null(treat)==FALSE)  rownames(R1) <- treat
  colnames(R1) <- paste0("Probability of rank ",1:p)
  
  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1
    
  if(p<=4) par(mfrow=c(2,2))
  if(p==5||p==6) par(mfrow=c(2,3))
  if(p==7||p==8) par(mfrow=c(2,4))
  if(p==9||p==10) par(mfrow=c(2,5))
  if(p>10&&p<=15) par(mfrow=c(3,5))
  if(p>15) par(mfrow=c(4,5))

  if(p<=20){
    if(is.null(treat))  for(i in 1:p)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=paste0("Treatment ",i),ylim=c(0,1))
    if(is.null(treat)==FALSE)  for(i in 1:p)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=treat[i],ylim=c(0,1))
  }
  
  if(p>20&&p<=40){	
    if(is.null(treat))  for(i in 1:20)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=paste0("Treatment ",i),ylim=c(0,1))
    if(is.null(treat)==FALSE)  for(i in 1:20)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=treat[i],ylim=c(0,1))
	message("The panels will be changed after 10 seconds.")
    Sys.sleep(10)
    if(is.null(treat))  for(i in 21:p)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=paste0("Treatment ",i),ylim=c(0,1))
    if(is.null(treat)==FALSE)  for(i in 21:p)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=treat[i],ylim=c(0,1))
  }

  if(p>40)  message("The ranking probability diagrams cannot be drawn for this setting.")
  
  MEANRANK <- apply(x2,1,mean)

  # R2 <- apply(R1,1,cumsum)

  sp0 <- approxfun(c(1,p),c(0,1),method = "linear")
  sucra <- 1 - sp0(MEANRANK)

  R2 <- cbind(sucra,MEANRANK,R1)
  colnames(R2)[1:2] <- c("SUCRA","MEANRANK")

  ord <- order(sucra,decreasing=TRUE)

  R3 <- list(SUCRA=t(t(R2[ord,1])),MEANRANK=t(t(R2[ord,2])),RANKPROB=R2[ord,3:dim(R2)[2]])
  
  return(R3)

 }
 
 if(ascending==FALSE){

  x1 <- mvrnorm(B, mu, V)
  x1 <- cbind(rep(0,times=B),x1)

  x2 <- apply(x1,1,rank)
  x2 <- max(x2) - x2 + 1

  p <- dim(x2)[1]
  
  R1 <- matrix(rep(NA,times=p*p),p)

  for(i in 1:p){
	for(j in 1:p){
		R1[i,j] <- sum(x2[i,]==j)/B
	}
  }
  
  if(is.null(treat)==FALSE)  rownames(R1) <- treat
  colnames(R1) <- paste0("Probability of rank ",1:p)
  
  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1
    
  if(p<=4) par(mfrow=c(2,2))	
  if(p==5||p==6) par(mfrow=c(2,3))
  if(p==7||p==8) par(mfrow=c(2,4))
  if(p==9||p==10) par(mfrow=c(2,5))
  if(p>10&&p<=15) par(mfrow=c(3,5))
  if(p>15) par(mfrow=c(4,5))

  if(p<=20){
    if(is.null(treat))  for(i in 1:p)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=paste0("Treatment ",i),ylim=c(0,1))
    if(is.null(treat)==FALSE)  for(i in 1:p)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=treat[i],ylim=c(0,1))
  }
  
  if(p>20&&p<=40){		
    if(is.null(treat))  for(i in 1:20)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=paste0("Treatment ",i),ylim=c(0,1))
    if(is.null(treat)==FALSE)  for(i in 1:20)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=treat[i],ylim=c(0,1))
	message("The panels will be changed after 10 seconds.")
    Sys.sleep(10)
    if(is.null(treat))  for(i in 21:p)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=paste0("Treatment ",i),ylim=c(0,1))
    if(is.null(treat)==FALSE)  for(i in 21:p)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=treat[i],ylim=c(0,1))
  }

  if(p>40)  message("The ranking probability diagrams cannot be drawn for this setting.")
  
  MEANRANK <- apply(x2,1,mean)

  # R2 <- apply(R1,1,cumsum)

  sp0 <- approxfun(c(1,p),c(0,1),method = "linear")
  sucra <- 1 - sp0(MEANRANK)

  R2 <- cbind(sucra,MEANRANK,R1)
  colnames(R2)[1:2] <- c("SUCRA","MEANRANK")

  ord <- order(sucra,decreasing=TRUE)

  R3 <- list(SUCRA=t(t(R2[ord,1])),MEANRANK=t(t(R2[ord,2])),RANKPROB=R2[ord,3:dim(R2)[2]])
  
  return(R3)

 }
 
} 

if(method=="REML"){

 R1 <- REML(y,S)
 mu <- R1$mu
 V <- R1$V
  
 if(ascending==TRUE){

  x1 <- mvrnorm(B, mu, V)
  x1 <- cbind(rep(0,times=B),x1)

  x2 <- apply(x1,1,rank)

  p <- dim(x2)[1]
  
  R1 <- matrix(rep(NA,times=p*p),p)

  for(i in 1:p){
	for(j in 1:p){
		R1[i,j] <- sum(x2[i,]==j)/B
	}
  }

  if(is.null(treat)==FALSE)  rownames(R1) <- treat
  colnames(R1) <- paste0("Probability of rank ",1:p)
  
  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1
    
  if(p<=4) par(mfrow=c(2,2))	
  if(p==5||p==6) par(mfrow=c(2,3))
  if(p==7||p==8) par(mfrow=c(2,4))
  if(p==9||p==10) par(mfrow=c(2,5))
  if(p>10&&p<=15) par(mfrow=c(3,5))
  if(p>15) par(mfrow=c(4,5))

  if(p<=20){
    if(is.null(treat))  for(i in 1:p)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=paste0("Treatment ",i),ylim=c(0,1))
    if(is.null(treat)==FALSE)  for(i in 1:p)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=treat[i],ylim=c(0,1))
  }
  
  if(p>20&&p<=40){		
    if(is.null(treat))  for(i in 1:20)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=paste0("Treatment ",i),ylim=c(0,1))
    if(is.null(treat)==FALSE)  for(i in 1:20)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=treat[i],ylim=c(0,1))
	message("The panels will be changed after 10 seconds.")
    Sys.sleep(10)
    if(is.null(treat))  for(i in 21:p)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=paste0("Treatment ",i),ylim=c(0,1))
    if(is.null(treat)==FALSE)  for(i in 21:p)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=treat[i],ylim=c(0,1))
  }

  if(p>40)  message("The ranking probability diagrams cannot be drawn for this setting.")
  
  MEANRANK <- apply(x2,1,mean)

  # R2 <- apply(R1,1,cumsum)

  sp0 <- approxfun(c(1,p),c(0,1),method = "linear")
  sucra <- 1 - sp0(MEANRANK)

  R2 <- cbind(sucra,MEANRANK,R1)
  colnames(R2)[1:2] <- c("SUCRA","MEANRANK")

  ord <- order(sucra,decreasing=TRUE)

  R3 <- list(SUCRA=t(t(R2[ord,1])),MEANRANK=t(t(R2[ord,2])),RANKPROB=R2[ord,3:dim(R2)[2]])
 
  return(R3)

 }
 
 if(ascending==FALSE){

  x1 <- mvrnorm(B, mu, V)
  x1 <- cbind(rep(0,times=B),x1)

  x2 <- apply(x1,1,rank)
  x2 <- max(x2) - x2 + 1

  p <- dim(x2)[1]
  
  R1 <- matrix(rep(NA,times=p*p),p)

  for(i in 1:p){
	for(j in 1:p){
		R1[i,j] <- sum(x2[i,]==j)/B
	}
  }
  
  if(is.null(treat)==FALSE)  rownames(R1) <- treat
  colnames(R1) <- paste0("Probability of rank ",1:p)
  
  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1
    
  if(p<=4) par(mfrow=c(2,2))	
  if(p==5||p==6) par(mfrow=c(2,3))
  if(p==7||p==8) par(mfrow=c(2,4))
  if(p==9||p==10) par(mfrow=c(2,5))
  if(p>10&&p<=15) par(mfrow=c(3,5))
  if(p>15) par(mfrow=c(4,5))

  if(p<=20){
    if(is.null(treat))  for(i in 1:p)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=paste0("Treatment ",i),ylim=c(0,1))
    if(is.null(treat)==FALSE)  for(i in 1:p)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=treat[i],ylim=c(0,1))
  }
  
  if(p>20&&p<=40){		
    if(is.null(treat))  for(i in 1:20)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=paste0("Treatment ",i),ylim=c(0,1))
    if(is.null(treat)==FALSE)  for(i in 1:20)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=treat[i],ylim=c(0,1))
	message("The panels will be changed after 10 seconds.")
    Sys.sleep(10)
    if(is.null(treat))  for(i in 21:p)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=paste0("Treatment ",i),ylim=c(0,1))
    if(is.null(treat)==FALSE)  for(i in 21:p)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=treat[i],ylim=c(0,1))
  }

  if(p>40)  message("The ranking probability diagrams cannot be drawn for this setting.")
  
  MEANRANK <- apply(x2,1,mean)

  # R2 <- apply(R1,1,cumsum)

  sp0 <- approxfun(c(1,p),c(0,1),method = "linear")
  sucra <- 1 - sp0(MEANRANK)

  R2 <- cbind(sucra,MEANRANK,R1)
  colnames(R2)[1:2] <- c("SUCRA","MEANRANK")

  ord <- order(sucra,decreasing=TRUE)

  R3 <- list(SUCRA=t(t(R2[ord,1])),MEANRANK=t(t(R2[ord,2])),RANKPROB=R2[ord,3:dim(R2)[2]])

  return(R3)

 }
 
} 

if(method=="fixed"){

FEM <- function(y,S){
  
  N <- dim(y)[1]
  p <- dim(y)[2]
  
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
      
      Wi <- ginv2(Si)
      
      A1 <- A1 + ivec(yi %*% Wi, wi, p)
      A2 <- A2 + imat(Wi, wi, p)
      
    }
    
  mu <- A1 %*% ginv2(A2)
   
  V.mu <- ginv2(A2)
  
####  

  A3 <- 0

    for(i in 1:N){
      
      yi <- as.vector(y[i,])
      wi <- which(is.na(yi)==FALSE)
      pl <- length(wi)
      
      Si <- vmat(S[i,], p)
      
      yi <- yi[wi]
      Si <- pmat(Si,wi)
      
      Wi <- ginv2(Si)
	  
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
    list("mu"=mu,"V"=V.mu)
  )
  
}

 R1 <- FEM(y,S)
 mu <- R1$mu
 V <- R1$V
  
 if(ascending==TRUE){

  x1 <- mvrnorm(B, mu, V)
  x1 <- cbind(rep(0,times=B),x1)

  x2 <- apply(x1,1,rank)

  p <- dim(x2)[1]
  
  R1 <- matrix(rep(NA,times=p*p),p)

  for(i in 1:p){
	for(j in 1:p){
		R1[i,j] <- sum(x2[i,]==j)/B
	}
  }

  if(is.null(treat)==FALSE)  rownames(R1) <- treat
  colnames(R1) <- paste0("Probability of rank ",1:p)
  
  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1
    
  if(p<=4) par(mfrow=c(2,2))	
  if(p==5||p==6) par(mfrow=c(2,3))
  if(p==7||p==8) par(mfrow=c(2,4))
  if(p==9||p==10) par(mfrow=c(2,5))
  if(p>10&&p<=15) par(mfrow=c(3,5))
  if(p>15) par(mfrow=c(4,5))

  if(p<=20){
    if(is.null(treat))  for(i in 1:p)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=paste0("Treatment ",i),ylim=c(0,1))
    if(is.null(treat)==FALSE)  for(i in 1:p)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=treat[i],ylim=c(0,1))
  }
  
  if(p>20&&p<=40){		
    if(is.null(treat))  for(i in 1:20)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=paste0("Treatment ",i),ylim=c(0,1))
    if(is.null(treat)==FALSE)  for(i in 1:20)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=treat[i],ylim=c(0,1))
	message("The panels will be changed after 10 seconds.")
    Sys.sleep(10)
    if(is.null(treat))  for(i in 21:p)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=paste0("Treatment ",i),ylim=c(0,1))
    if(is.null(treat)==FALSE)  for(i in 21:p)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=treat[i],ylim=c(0,1))
  }

  if(p>40)  message("The ranking probability diagrams cannot be drawn for this setting.")
  
  MEANRANK <- apply(x2,1,mean)

  # R2 <- apply(R1,1,cumsum)

  sp0 <- approxfun(c(1,p),c(0,1),method = "linear")
  sucra <- 1 - sp0(MEANRANK)

  R2 <- cbind(sucra,MEANRANK,R1)
  colnames(R2)[1:2] <- c("SUCRA","MEANRANK")

  ord <- order(sucra,decreasing=TRUE)

  R3 <- list(SUCRA=t(t(R2[ord,1])),MEANRANK=t(t(R2[ord,2])),RANKPROB=R2[ord,3:dim(R2)[2]])
 
  return(R3)

 }
 
 if(ascending==FALSE){

  x1 <- mvrnorm(B, mu, V)
  x1 <- cbind(rep(0,times=B),x1)

  x2 <- apply(x1,1,rank)
  x2 <- max(x2) - x2 + 1

  p <- dim(x2)[1]
  
  R1 <- matrix(rep(NA,times=p*p),p)

  for(i in 1:p){
	for(j in 1:p){
		R1[i,j] <- sum(x2[i,]==j)/B
	}
  }
  
  if(is.null(treat)==FALSE)  rownames(R1) <- treat
  colnames(R1) <- paste0("Probability of rank ",1:p)
  
  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1
    
  if(p<=4) par(mfrow=c(2,2))	
  if(p==5||p==6) par(mfrow=c(2,3))
  if(p==7||p==8) par(mfrow=c(2,4))
  if(p==9||p==10) par(mfrow=c(2,5))
  if(p>10&&p<=15) par(mfrow=c(3,5))
  if(p>15) par(mfrow=c(4,5))

  if(p<=20){
    if(is.null(treat))  for(i in 1:p)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=paste0("Treatment ",i),ylim=c(0,1))
    if(is.null(treat)==FALSE)  for(i in 1:p)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=treat[i],ylim=c(0,1))
  }
  
  if(p>20&&p<=40){		
    if(is.null(treat))  for(i in 1:20)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=paste0("Treatment ",i),ylim=c(0,1))
    if(is.null(treat)==FALSE)  for(i in 1:20)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=treat[i],ylim=c(0,1))
	message("The panels will be changed after 10 seconds.")
    Sys.sleep(10)
    if(is.null(treat))  for(i in 21:p)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=paste0("Treatment ",i),ylim=c(0,1))
    if(is.null(treat)==FALSE)  for(i in 21:p)		plot(R1[i,],type="l",lwd=2,col="blue",xlab="Rank",ylab="Probability",main=treat[i],ylim=c(0,1))
  }

  if(p>40)  message("The ranking probability diagrams cannot be drawn for this setting.")
  
  MEANRANK <- apply(x2,1,mean)

  # R2 <- apply(R1,1,cumsum)

  sp0 <- approxfun(c(1,p),c(0,1),method = "linear")
  sucra <- 1 - sp0(MEANRANK)

  R2 <- cbind(sucra,MEANRANK,R1)
  colnames(R2)[1:2] <- c("SUCRA","MEANRANK")

  ord <- order(sucra,decreasing=TRUE)

  R3 <- list(SUCRA=t(t(R2[ord,1])),MEANRANK=t(t(R2[ord,2])),RANKPROB=R2[ord,3:dim(R2)[2]])

  return(R3)

 }
 
} 

}


