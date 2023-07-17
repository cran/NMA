nma <- function(x, eform=FALSE, method="NH"){

y <- x$y
S <- x$S
dof <- x$df


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
    list(Expected = list(SE = SE_E, df = df_E, lambda = lambda_E, CR = Phi_A_E), 
         Observed = list(SE = SE_O, df = df_O, lambda = lambda_O, CR = Phi_A_O))
  )
  
}

REML <- function(y,S,maxitr=20){
  
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
    list("Coefficients"=R5,"CF"=V.mu,"Q"=A3)
  )
  
}

tPI <- function(y, S){
  
  N <- dim(y)[1]
  p <- dim(y)[2]
  
  A8 <- numeric(p)

  for(i in 1:N){

	vi <- diag(vmat(S[i,],p))
	A8[which(vi<1000)] <- A8[which(vi<1000)] + 1

  }
  
  reml0 <- REML(y,S)
  fem0 <- FEM(y,S)

  e1 <- reml0[[1]][,1]
  v1 <- reml0[[1]][,2]^2
  tau2 <- reml0[[2]]^2
  
  R <- ( det(reml0$CR)/det(fem0$CF) )^(1/(2*p))			# Jackson's multivariate R-statistic
  I2 <- (R^2 - 1)/(R^2)									# Jackson's multivariate I2^statistic

  # I2 <- (det(reml0$CR)^(1/p) - det(fem0$CF)^(1/p))/( det(reml0$CR)^(1/p) )		# alternative formula of I2

  pl <- e1 - qt(0.975,df=N-A8-1)*sqrt(tau2+v1)		# prediction interval
  pu <- e1 + qt(0.975,df=N-A8-1)*sqrt(tau2+v1)

  R3 <- list("Estimates" = reml0[[1]],"Between-studies_SD" = reml0[[2]],"Multivariate I2-statistic" = I2,
					  "95%PI" = cbind(pl,pu))
  
  return(R3)
  
}

KR <- function(y, S){
  
  result_reml <- REML(y,S)
  tau_reml <- result_reml[[2]]
  kr <- KR0(y = y, S = S, tau2 = tau_reml^2)
  
  kr[["Expected"]][["df"]][kr[["Expected"]][["df"]] < 3] <- 3		# truncation of the DF
  
  KR_E <- list("Coefficients" = cbind("Coef." = result_reml[[1]][, 1], 
                                             "SE" = kr[["Expected"]][["SE"]], 
                                             "95%CL" = result_reml[[1]][, 1] - kr[["Expected"]][["SE"]] * qt(0.975, df = kr[["Expected"]][["df"]]), 
                                             "95%CU" = result_reml[[1]][, 1] + kr[["Expected"]][["SE"]] * qt(0.975, df = kr[["Expected"]][["df"]]), 
                                             "df" = kr[["Expected"]][["df"]]), 
                      "Between-studies_SD" = tau_reml)
  
  #
  
  pl <- KR_E[[1]][,1] - qt(0.975,df=KR_E[[1]][,5]-1)*sqrt(KR_E[[1]][,2]^2 + KR_E[[2]]^2)
  pu <- KR_E[[1]][,1] + qt(0.975,df=KR_E[[1]][,5]-1)*sqrt(KR_E[[1]][,2]^2 + KR_E[[2]]^2)

  p <- dim(y)[2]
  fem0 <- FEM(y,S)
  R <- ( det(kr$Expected$CR)/det(fem0$CF) )^(1/(2*p))			# Jackson's multivariate R-statistic
  I2 <- (R^2 - 1)/(R^2)									# Jackson's multivariate I2^statistic

  R3 <- list("Estimates" = cbind("Coef." = result_reml[[1]][, 1], 
                                             "SE" = kr[["Expected"]][["SE"]], 
                                             "95%CL" = result_reml[[1]][, 1] - kr[["Expected"]][["SE"]] * qt(0.975, df = kr[["Expected"]][["df"]]), 
                                             "95%CU" = result_reml[[1]][, 1] + kr[["Expected"]][["SE"]] * qt(0.975, df = kr[["Expected"]][["df"]]), 
											 "P-value" = 2*(1-pt(abs(result_reml[[1]][, 1])/kr[["Expected"]][["SE"]], df=kr[["Expected"]][["df"]]))),
                                             #"df" = kr[["Expected"]][["df"]]), 
                      "Between-studies_SD" = tau_reml,"Multivariate I2-statistic" = I2,
					  "95%PI" = cbind(pl,pu))
  
  return(R3)
  
}

etPI <- function(y, S){
  
  N <- dim(y)[1]
  p <- dim(y)[2]
  
  A8 <- numeric(p)

  for(i in 1:N){

	vi <- diag(vmat(S[i,],p))
	A8[which(vi<1000)] <- A8[which(vi<1000)] + 1

  }
  
  reml0 <- REML(y,S)
  fem0 <- FEM(y,S)

  e1 <- reml0[[1]][,1]
  v1 <- reml0[[1]][,2]^2
  tau2 <- reml0[[2]]^2
  
  R <- ( det(reml0$CR)/det(fem0$CF) )^(1/(2*p))			# Jackson's multivariate R-statistic
  I2 <- (R^2 - 1)/(R^2)									# Jackson's multivariate I2^statistic

  pl <- exp(e1 - qt(0.975,df=N-A8-1)*sqrt(tau2+v1))		# prediction interval
  pu <- exp(e1 + qt(0.975,df=N-A8-1)*sqrt(tau2+v1))

  reml0[[1]][,1] <- exp(reml0[[1]][,1])
  reml0[[1]][,3] <- exp(reml0[[1]][,3])
  reml0[[1]][,4] <- exp(reml0[[1]][,4])

  R3 <- list("Estimates" = reml0[[1]],"Between-studies_SD" = reml0[[2]],"Multivariate I2-statistic" = I2,
					  "95%PI" = cbind(pl,pu))
    
  return(R3)
  
}

eKR <- function(y, S){
  
  result_reml <- REML(y,S)
  tau_reml <- result_reml[[2]]
  kr <- KR0(y = y, S = S, tau2 = tau_reml^2)
  
  kr[["Expected"]][["df"]][kr[["Expected"]][["df"]] < 3] <- 3		# truncation of the DF
  
  KR_E <- list("Coefficients" = cbind("Coef." = result_reml[[1]][, 1], 
                                             "SE" = kr[["Expected"]][["SE"]], 
                                             "95%CL" = result_reml[[1]][, 1] - kr[["Expected"]][["SE"]] * qt(0.975, df = kr[["Expected"]][["df"]]), 
                                             "95%CU" = result_reml[[1]][, 1] + kr[["Expected"]][["SE"]] * qt(0.975, df = kr[["Expected"]][["df"]]), 
                                             "df" = kr[["Expected"]][["df"]]), 
                      "Between-studies_SD" = tau_reml)
  
  #

  p <- dim(y)[2]
  fem0 <- FEM(y,S)
  R <- ( det(kr$Expected$CR)/det(fem0$CF) )^(1/(2*p))			# Jackson's multivariate R-statistic
  I2 <- (R^2 - 1)/(R^2)									# Jackson's multivariate I2^statistic
  
  pl <- exp(KR_E[[1]][,1] - qt(0.975,df=KR_E[[1]][,5]-1)*sqrt(KR_E[[1]][,2]^2 + KR_E[[2]]^2))
  pu <- exp(KR_E[[1]][,1] + qt(0.975,df=KR_E[[1]][,5]-1)*sqrt(KR_E[[1]][,2]^2 + KR_E[[2]]^2))

  R3 <- list("Estimates" = cbind("Coef." = exp(result_reml[[1]][, 1]), 
                                             "SE" = kr[["Expected"]][["SE"]], 
                                             "95%CL" = exp(result_reml[[1]][, 1] - kr[["Expected"]][["SE"]] * qt(0.975, df = kr[["Expected"]][["df"]])), 
                                             "95%CU" = exp(result_reml[[1]][, 1] + kr[["Expected"]][["SE"]] * qt(0.975, df = kr[["Expected"]][["df"]])),
											 "P-value" = 2*(1-pt(abs(result_reml[[1]][, 1])/kr[["Expected"]][["SE"]], df=kr[["Expected"]][["df"]]))),
											 # "df" = kr[["Expected"]][["df"]]), 
                      "Between-studies_SD" = tau_reml,"Multivariate I2-statistic" = I2,
					  "95%PI" = cbind(pl,pu))
    
   return(R3)
  
}

fem0 <- FEM(y,S)
Qs <- fem0$Q
Ps <- 1 - pchisq(Qs,df=dof)
Rs <- paste0("Q(df = ",round(dof),") = ",round(Qs,4),", p-value = ",round(Ps,4))

Hs <- max( as.numeric(Qs/dof), 0)
IH <- max(as.numeric( (Hs-1)/Hs ), 0)

if(eform==FALSE){
	if(method=="NH")	R1 <- KR(y,S)
	if(method=="REML")	R1 <- tPI(y,S)
}

if(eform==TRUE){
	if(method=="NH")	R1 <- eKR(y,S)
	if(method=="REML")	R1 <- etPI(y,S)
}

if(method=="NH"||method=="REML"){

	R2 <- R1[[1]]
	R3 <- R1[[2]]
	R4 <- max(R1[[3]],0)
	R5 <- R1[[4]]

	N <- dim(y)[1]
	p <- dim(y)[2] + 1
	rownames(R2) <- paste0(2:p,": cons")
	rownames(R5) <- paste0(2:p,": cons")

	if(method=="NH") method <- "Noma-Hamura's improved REML-based inference and prediction methods"
	if(method=="REML") method <- "Restricted maximum likelihood (REML) estimation and prediction based on the ordinary t-approximation"

	R5 <- list("coding"=x$coding,"reference"=x$reference,"number of studies"=N,"method"=method,"Coef. (vs. treat 1)"=R2,"tau (Between-studies_SD) estimate"=R3,"tau2 (Between-studies_variance) estimate"=(R3*R3),"Multivariate H2-statistic"=Hs,"Multivariate I2-statistic"=R4,"Test for Heterogeneity"=Rs,"95%PI (vs. treat 1)"=R5)

}

if(method=="fixed"){

	method <- "Fixed-effect model"
	N <- dim(y)[1]
	p <- dim(y)[2] + 1
	rownames(fem0[[1]]) <- paste0(2:p,": cons")
	
	if(eform==TRUE){
		fem0[[1]][,1] <- exp(fem0[[1]][,1])
		fem0[[1]][,3] <- exp(fem0[[1]][,3])
		fem0[[1]][,4] <- exp(fem0[[1]][,4])
	}
	
	R5 <- list("coding"=x$coding,"reference"=x$reference,"number of studies"=N,"method"=method,"Coef. (vs. treat 1)"=fem0[[1]],"Test for Heterogeneity"=Rs)

}

return(R5)

}


