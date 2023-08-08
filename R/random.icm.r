random.icm <- function(x){

	xms <- x$measure
	
	N <- x$N
	p <- x$p

	if(xms=="OR"||xms=="RR"||xms=="RD"){

	study <- x$study
	treat <- x$treat
	d <- x$d
	n <- x$n
	
	study <- as.numeric(factor(study))

	dtx <- data.frame(study,treat,d,n)
	dtx <- dtx[order(treat),]

	studlab <- treat1 <- treat2 <- d1 <- d2 <- n1 <- n2 <- NULL
		
	for(i in 1:N){
	
		wi <- which(dtx$study==i)
		
		if(length(wi)==2){
		
			studlab <- c(studlab,i)
			treat1 <- c(treat1,dtx$treat[wi][1])
			treat2 <- c(treat2,dtx$treat[wi][2])
			d1 <- c(d1,dtx$d[wi][1])
			d2 <- c(d2,dtx$d[wi][2])
			n1 <- c(n1,dtx$n[wi][1])
			n2 <- c(n2,dtx$n[wi][2])
		
		}
		
		if(length(wi)>=3){
		
			for(j in 2:length(wi)){
		
				studlab <- c(studlab,i)
				treat1 <- c(treat1,dtx$treat[wi][1])
				treat2 <- c(treat2,dtx$treat[wi][j])
				d1 <- c(d1,dtx$d[wi][1])
				d2 <- c(d2,dtx$d[wi][j])
				n1 <- c(n1,dtx$n[wi][1])
				n2 <- c(n2,dtx$n[wi][j])

			}
		
		}
			
	}

	jwi <- data.frame(studlab,treat1,treat2,d1,d2,n1,n2)
	colnames(jwi) <- c("study","ref","trt","d1","d2","n1","n2")

	stud <- unique(studlab)
	N <- length(stud)
	
	L <- length(studlab)
		
	contr <- design <- rep(NA,times=L)
	
	for(i in 1:L)	contr[i] <- paste0(treat1[i],"-",treat2[i])
	
	for(i in 1:N){
	
		wi <- which(studlab==stud[i])
		
		desi <- sort(unique(c(treat1[wi],treat2[wi])))
		desj <- desi[1]
		for(j in 2:length(desi)) desj <- paste0(desj,"-",desi[j])
		design[wi] <- desj
		
	}
	
	jwi <- cbind(jwi,contr,design)
	
	###
		
	if(xms=="OR"){

		y <- log(jwi$d1) - log(jwi$n1 - jwi$d1) - log(jwi$d2) + log(jwi$n2 - jwi$d2)
		y[abs(y)==Inf] <- log(jwi$d1[abs(y)==Inf] + .5) - log(jwi$n1[abs(y)==Inf] - jwi$d1[abs(y)==Inf] + 1) - log(jwi$d2[abs(y)==Inf] + .5) + log(jwi$n2[abs(y)==Inf] - jwi$d2[abs(y)==Inf] + 1)
		jwi$y <- y

		VTE <- (jwi$d1)^-1 + (jwi$d2)^-1 + (jwi$n1 - jwi$d1)^-1 + (jwi$n2 - jwi$d2)^-1
		VTE[VTE==Inf] <- (jwi$d1[VTE==Inf]+.5)^-1 + (jwi$d2[VTE==Inf]+.5)^-1 + (jwi$n1[VTE==Inf] - jwi$d1[VTE==Inf] + .5)^-1 + (jwi$n2[VTE==Inf] - jwi$d2[VTE==Inf] +.5)^-1

		Vd <- diag(VTE)

		for(i in 1:N){
	
			wi <- which(jwi$study==i)
		
			if(length(wi)>=2){
		
				w1 <- wi[1]
				C1 <- (jwi$d1[w1])^-1 + (jwi$n1[w1] - jwi$d1[w1])^-1
				if(C1==Inf)	C1 <- (jwi$d1[w1]+.5)^-1 + (jwi$n1[w1] - jwi$d1[w1]+.5)^-1
			
				for(j in 1:length(wi)){
					ai <- 1:length(wi)
					bi <- setdiff(ai,j)
					Vd[wi[j],wi[bi]] <- C1
				}
				
			}
	
		}
	
	}
	
	if(xms=="RR"){

		y <- log(jwi$d1) - log(jwi$n1) - log(jwi$d2) + log(jwi$n2)
		y[abs(y)==Inf] <- log(jwi$d1[abs(y)==Inf] + .5) - log(jwi$n1[abs(y)==Inf] + 1) - log(jwi$d2[abs(y)==Inf] + .5) + log(jwi$n2[abs(y)==Inf] + 1)
		jwi$y <- y

		VTE <- (jwi$d1)^-1 + (jwi$d2)^-1 - (jwi$n1)^-1 - (jwi$n2)^-1
		VTE[VTE==Inf] <- (jwi$d1[VTE==Inf]+.5)^-1 + (jwi$d2[VTE==Inf]+.5)^-1 - (jwi$n1[VTE==Inf] + 1)^-1 - (jwi$n2[VTE==Inf] + 1)^-1

		Vd <- diag(VTE)

		for(i in 1:N){
	
			wi <- which(jwi$study==i)
		
			if(length(wi)>=2){
		
				w1 <- wi[1]
				C1 <- (jwi$d1[w1])^-1 - (jwi$n1[w1])^-1
				if(C1==Inf)	C1 <- (jwi$d1[w1] + .5)^-1 - (jwi$n1[w1] + 1)^-1
			
				for(j in 1:length(wi)){
					ai <- 1:length(wi)
					bi <- setdiff(ai,j)
					Vd[wi[j],wi[bi]] <- C1
				}
				
			}
	
		}
	
	}

	if(xms=="RD"){

		y <- (jwi$d1/jwi$n1) - (jwi$d2/jwi$n2)
		jwi$y <- y

		VTE <- (jwi$d1)*(jwi$n1 - jwi$d1)/(jwi$n1^3) + (jwi$d2)*(jwi$n2 - jwi$d2)/(jwi$n2^3)
		VTE[VTE==Inf] <- (jwi$d1[VTE==Inf] + .5)*(jwi$n1[VTE==Inf] - jwi$d1[VTE==Inf] + .5)/((jwi$n1[VTE==Inf]+1)^3) + (jwi$d2[VTE==Inf] + .5)*(jwi$n2[VTE==Inf] - jwi$d2[VTE==Inf] + .5)/((jwi$n2[VTE==Inf]+1)^3)

		Vd <- diag(VTE)

		for(i in 1:N){
	
			wi <- which(jwi$study==i)
		
			if(length(wi)>=2){
		
				w1 <- wi[1]
				C1 <- (jwi$d1[w1])*(jwi$n1[w1] - jwi$d1[w1])/(jwi$n1[w1]^3)
				if(C1==Inf)	C1 <- (jwi$d1[w1] + .5)*(jwi$n1[w1] - jwi$d1[w1] + .5)/((jwi$n1[w1]+1)^3)
			
				for(j in 1:length(wi)){
					ai <- 1:length(wi)
					bi <- setdiff(ai,j)
					Vd[wi[j],wi[bi]] <- C1
				}
				
			}
	
		}
	
	}
	
	}
	
	if(xms=="MD"||xms=="SMD"){

	study <- x$study
	treat <- x$treat
	m <- x$m
	s <- x$s
	n <- x$n
	
	study <- as.numeric(factor(study))

	dtx <- data.frame(study,treat,m,s,n)
	dtx <- dtx[order(treat),]

	studlab <- treat1 <- treat2 <- m1 <- m2 <- s1 <- s2 <- n1 <- n2 <- NULL
		
	for(i in 1:N){
	
		wi <- which(dtx$study==i)
		
		if(length(wi)==2){
		
			studlab <- c(studlab,i)
			treat1 <- c(treat1,dtx$treat[wi][1])
			treat2 <- c(treat2,dtx$treat[wi][2])
			m1 <- c(m1,dtx$m[wi][1])
			m2 <- c(m2,dtx$m[wi][2])
			s1 <- c(s1,dtx$s[wi][1])
			s2 <- c(s2,dtx$s[wi][2])
			n1 <- c(n1,dtx$n[wi][1])
			n2 <- c(n2,dtx$n[wi][2])
		
		}
		
		if(length(wi)>=3){
		
			for(j in 2:length(wi)){
		
				studlab <- c(studlab,i)
				treat1 <- c(treat1,dtx$treat[wi][1])
				treat2 <- c(treat2,dtx$treat[wi][j])
				m1 <- c(m1,dtx$m[wi][1])
				m2 <- c(m2,dtx$m[wi][j])
				s1 <- c(s1,dtx$s[wi][1])
				s2 <- c(s2,dtx$s[wi][j])
				n1 <- c(n1,dtx$n[wi][1])
				n2 <- c(n2,dtx$n[wi][j])

			}
		
		}
			
	}

	jwi <- data.frame(studlab,treat1,treat2,m1,m2,s1,s2,n1,n2)
	colnames(jwi) <- c("study","ref","trt","m1","m2","s1","s2","n1","n2")

	stud <- unique(studlab)

	N <- length(stud)
	
	L <- length(studlab)
		
	contr <- design <- rep(NA,times=L)
	
	for(i in 1:L)	contr[i] <- paste0(treat1[i],"-",treat2[i])
	
	for(i in 1:N){
	
		wi <- which(studlab==stud[i])
		
		desi <- sort(unique(c(treat1[wi],treat2[wi])))
		desj <- desi[1]
		for(j in 2:length(desi)) desj <- paste0(desj,"-",desi[j])
		design[wi] <- desj
		
	}
	
	jwi <- cbind(jwi,contr,design)
	
	###
		
	if(xms=="MD"){

		y <- jwi$m1 - jwi$m2
		jwi$y <- y

		VTE <- (jwi$s1)^2 / (jwi$n1) + (jwi$s2)^2 / (jwi$n2)

		Vd <- diag(VTE)

		for(i in 1:N){
	
			wi <- which(jwi$study==i)
		
			if(length(wi)>=2){
		
				w1 <- wi[1]
				C1 <- (jwi$s1[w1])^2 / (jwi$n1[w1])
			
				for(j in 1:length(wi)){
					ai <- 1:length(wi)
					bi <- setdiff(ai,j)
					Vd[wi[j],wi[bi]] <- C1
				}
				
			}
	
		}
	
	}
	
	if(xms=="SMD"){

		std <- sqrt( ( (jwi$n1 - 1)*(jwi$s1*jwi$s1) + (jwi$n2 - 1)*(jwi$s2*jwi$s2) ) / (jwi$n1 + jwi$n2 - 2) )
		y <- (jwi$m1 - jwi$m2)/std
		jwi$y <- y

		VTE <- 1/(jwi$n1) + 1/(jwi$n2)

		Vd <- diag(VTE)

		for(i in 1:N){
	
			wi <- which(jwi$study==i)
		
			if(length(wi)>=2){
		
				w1 <- wi[1]
				C1 <- 1 / (jwi$n1[w1])
			
				for(j in 1:length(wi)){
					ai <- 1:length(wi)
					bi <- setdiff(ai,j)
					Vd[wi[j],wi[bi]] <- C1
				}
				
			}
	
		}
	
	}
	
	}
	
	###
			
	contrmat <- function(trt1, trt2, ref) {
		all.lvls <- sort(unique(c(levels(factor(trt1)), levels(factor(trt2)))))
		trt1 <- factor(trt1, levels=all.lvls)
		trt2 <- factor(trt2, levels=all.lvls)
		X <- model.matrix(~ trt2 - 1) - model.matrix(~ trt1 - 1)
		colnames(X) <- all.lvls
		if (missing(ref))
			ref <- all.lvls[1]
		X[, colnames(X) != ref]
	}

	X <- contrmat(jwi$trt, jwi$ref)

	###
	
	# No inconsistency, Fixed-effect
	modF <- rma.mv(y, Vd, mods=X, intercept=FALSE, data=jwi)

	# No inconsistency, Random-effects
	modC <- rma.mv(y, Vd, mods=X, intercept=FALSE, random = ~ contr | study, rho=1/2, data=jwi)
	
	# Random-inconsistency, Random-effects
	modI <- rma.mv(y, Vd, mods=X, intercept=FALSE, random = list(~ contr | study, ~ contr | design), rho=1/2, phi=1/2, data=jwi)

	e1 <- modI$beta
	e2 <- modI$se
	e3 <- modI$ci.lb
	e4 <- modI$ci.ub
	e5 <- modI$pval
	
	R2 <- data.frame(e1,e2,e3,e4,e5)
	rownames(R2) <- paste0(2:(p+1),": cons")
	colnames(R2) <- c("Coef.","SE","95%CL","95%CU","P-value")

	###
	
	LR1 <- anova(modF, modC)
	LR2 <- anova(modI, modC)
	LR3 <- anova(modI, modF)

	df1 <- LR1$parms.f - LR1$parms.r
	A1 <- c(LR1$LRT,df1,LR1$pval)

	df2 <- LR2$parms.f - LR2$parms.r
	A2 <- c(LR2$LRT,df2,LR2$pval)

	df3 <- LR3$parms.f - LR3$parms.r
	A3 <- c(LR3$LRT,df3,LR3$pval)

	R3 <- t(data.frame(A1,A2,A3))
	rownames(R3) <- c("Heterogeneity", "Inconsistency", "Heterogeneity + Inconsistency")
	colnames(R3) <- c("X2-statistic","df","P-value")

	###

	Ra <- ( det(modC$vb)/det(modF$vb) )^(1/(2*p))
	I2_het <- (Ra^2 - 1)/(Ra^2)

	Rb <- ( det(modI$vb)/det(modC$vb) )^(1/(2*p))
	I2_inc <- (Rb^2 - 1)/(Rb^2)

	Rc <- ( det(modI$vb)/det(modF$vb) )^(1/(2*p))
	I2_both <- (Rc^2 - 1)/(Rc^2)
	
	R4 <- t(data.frame(c(Ra,I2_het),c(Rb,I2_inc),c(Rc,I2_both)))
	rownames(R4) <- c("Heterogeneity", "Inconsistency", "Heterogeneity + Inconsistency")
	colnames(R4) <- c("R-statistic","I2-statistic")

	###

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
	L <- length(des0)

	###

	C2 <- list("coding"=x$coding,"reference"=x$reference,"number of studies"=x$N,"number of designs"=L,designs=des0,"Coef. (vs. treat 1)"=R2,"Between-studies_SD"=sqrt(modI$tau2),"Between-designs_SD"=sqrt(modI$gamma2),"Likelihood ratio tests for the variance components"=R3,"Heterogeneity and inconsistency statistics"=R4)

	return(C2)

}

