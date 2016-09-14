####################################
# GENERIC USEFUL FUNCTIONS
####################################

# -------------------------- Generic functions ------------------------------------------

# Utility functions that perform a rotation of the elements of a vector
rotate.l <- function(x) {
	n<-length(x)
	if(n==1) x else{
		c(x[2:n],x[1])
	}
}

rotate.r <- function(x) {
	n<-length(x)
	if(n==1) x else{
		c(x[n],x[1:(n-1)])
	}
}

# Return an object of the same length as x, substituting the values from rep, matching elements by name
replace.par <- function(x,rep)
{
	if(length(rep)==0) return(x)
	pos <- sapply(names(rep),function(n) which(names(x)==n))
	replace(x,pos,rep)
}



# ------------------------- Birth pulse function ----------------------------------------

# Periodic Gaussian, using a sin function so that tau = peak time, and normalised so that <B> = 1 
# Parameters: t = time, Bp = list(s(bandwidth ≥0), tau (peak time)), 
# If length(s) == length(tau) == np AND length(t) == 1, return a vector of length np
# WARNING: do not use with length(t) > 1 AND length(s)*length(tau) > 1

Gauss.b <- function(t,p){
	exp( -p$s * (sin(pi*(t-p$tau)))^2 ) / besselI(p$s/2,0,T)
}

# ------------------------ Generate metapopulation matrices -----------------------------

# Generate contact matrices (square matrices of size n with 0 and 1) from a choice of types:
# - "circle": patch i linked to i-1 and i+1, plus a link between 1 and n.
# - "all": every patch linked to each other
meta.matrix <- function(n,type){
	switch(type,
		 circle = if(n==1) matrix(0,1) else if(n==2) matrix(c(0,1,1,0),2) else {
		 	x <- c(c(1),rep(0,n-3),c(1,0))
		 	sapply(1:n, function(i) {x <<- rotate.r(x)})
		 },
		 all = sapply(1:n,function(i) c(rep(1,i-1),0,rep(1,n-i)))
	)
}


