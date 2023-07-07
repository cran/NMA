rdc <- function(a,digits){

	if(a>=0){
	
	a1 <- round(a,digits)
	a2 <- as.character(floor(a1))
	a3 <- as.character(round(a1-floor(a1),digits))
	a4 <- str_split_fixed(a3, pattern=".", n=2)[2]
	L4 <- digits + 1 - nchar(a4)

	b1 <- a1 - floor(a1)	

	if((L4!=0)&&(L4<(digits+1))){
		for(k in 1:L4)  a4 <- paste0(a4,"0")
	}
		
	if(L4==(digits+1)){
		a4 <- "."
		for(k in 1:(L4-1))  a4 <- paste0(a4,"0")
	}
			
	a5 <- paste0(a2,a4)
	return(a5)

	}
	
	if(a<0){
	
	a <- abs(a)
	
	a1 <- round(a,digits)
	a2 <- as.character(floor(a1))
	a3 <- as.character(round(a1-floor(a1),digits))
	a4 <- str_split_fixed(a3, pattern=".", n=2)[2]
	L4 <- digits + 1 - nchar(a4)

	b1 <- a1 - floor(a1)	

	if((L4!=0)&&(L4<(digits+1))){
		for(k in 1:L4)  a4 <- paste0(a4,"0")
	}
		
	if(L4==(digits+1)){
		a4 <- "."
		for(k in 1:(L4-1))  a4 <- paste0(a4,"0")
	}
			
	a5 <- paste0("-",a2,a4)
	return(a5)

	}


}
