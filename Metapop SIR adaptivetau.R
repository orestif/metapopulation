####################################################
# 	METAPOPULATION MODEL FOR INFECTION DYNAMICS
#		WITH ANNUAL BIRTH PULSES
#	-------------------------------------------
#			DEFINITIONS
####################################################

# Olivier Restif, June 2015

# Stochastic SIR model with seasonal Gaussian birth rate 
# R version using package `adaptivetau`

library(adaptivetau)
library(foreach)
library(parallel)
library(doParallel)
library(deSolve)

# ================================== MODEL DEFINITIONS ===================================

# Variables: {S_i,I_i,R_i} 1≤i≤n 
# These can be handled by names (using string manipulation) or by position
# Order of variables: S_1,..., S_n, I_1, ..., I_n, R_1,... , R_n
# => S_i = var[i], I_i = var[i+n], R_i = var[i+2*n]

# Time unit = year

# ----------------------- Transition matrix -----------------------------------------------
# Create a transition matrix with 3*n rows (variables) and 6*n+3*k columns (events)
# where n=nrow(contact.matrix) is the number of patches and k=sum(contact.matrix) is the number of migration events
# Argument: contact.matrix: square matrix of size n with 0 and 1, M[i,j]=1 <=> migration from i to j
meta.SIR.transitions <- function(contact.matrix){
	n <- nrow(contact.matrix)
	var <- numeric(3*n)
	moves <- which(contact.matrix==1,T) # Matrix of all possible contacts
	cbind(
		# Births
		sapply(1:n,function(i){replace(var,i,+1)}),
		# Deaths (in the order of variables)
		sapply(1:(3*n),function(i){replace(var,i,-1)}),
		# Infection
		sapply(1:n,function(i){replace(var,c(i,i+n),c(-1,+1))}),
		# Recovery
		sapply(1:n,function(i){replace(var,c(i+n,i+2*n),c(-1,+1))}),
		# Migrations (S)
		if(nrow(moves)>0) sapply(1:nrow(moves),function(j){replace(var,moves[j,],c(-1,+1))}) else NULL,
		# Migrations (I)
		if(nrow(moves)>0) sapply(1:nrow(moves),function(j){replace(var,moves[j,]+n,c(-1,+1))}) else NULL,
		# Migrations (R)
		if(nrow(moves)>0) sapply(1:nrow(moves),function(j){replace(var,moves[j,]+2*n,c(-1,+1))}) else NULL
	)
}

# ----------------------- Transition rates -----------------------------------------------
# Calculate rates and return a vector of length 6*n+3*k 
# SIR model with FD transmission within patches, all individuals die and reproduce equally, directional migration
# Arguments:
# - x: vector of variable values, no names required
# - p: list of parameters, with elements:
# 	- n: number of patches
#	- contact.matrix: square matrix of size n with 0 and 1, M[i,j]=1 <=> migration from i to j
#	- B = function(t,Bpar) returning the birth rate per capita
# 	- Bp: list of parameters used by function B()
#	- d: death rate
#	- beta: transmission rate
#	- gamma: recovery rate
#	- mu: migration rate
meta.SIR.Rates <- function(x,p,t){
	moves <- which(p$contact.matrix==1,T) # Matrix of all possible contacts
	return(c(
		# Births: assume all individuals contribute
		# Allow birth rate function to vary among patches
		sapply(1:p$n,function(i){ p$B(t,p$Bp,i)*(x[i]+x[i+p$n]+x[i+2*p$n]) }),
		# Deaths: assume independent of status or patch
		sapply(1:(3*p$n),function(i){ p$d*x[i] }),
		# Infection: assume FD within patch -- test N>0 !!
		sapply(1:p$n,function(i){ 
			Ni = (x[i]+x[i+p$n]+x[i+2*p$n])
			if(Ni>0) p$beta*x[i]*x[i+p$n]/Ni else 0}),
		# Recovery
		sapply(1:p$n,function(i){ p$gamma*x[i+p$n] }),
		# Migrations
		if(nrow(moves)>0){ 
			c(sapply(1:nrow(moves),function(j){ p$mu*x[moves[j,1]] }),
			  sapply(1:nrow(moves),function(j){ p$mu*x[moves[j,1]+p$n] }),
			  sapply(1:nrow(moves),function(j){ p$mu*x[moves[j,1]+2*p$n] }))
		} else NULL
	))
}


# ------------------------- Birth pulse function ----------------------------------------

# Unscaled periodic Gaussian, using a sin function so that tau = peak time 
# Parameters: t = time, Bp = list(s(bandwidth ≥0), tau (peak time)) 
Gauss.B <- function(t,Bp){
	exp(-Bp$s*(sin(pi*(t-Bp$tau)))^2)
}

# Calculate the scaling factor k such that <B> = d. It is a function of d and s only.
B.k <- function(d,s){ if(s>0) d/integrate(Gauss.B,0,1,list(s=s,tau=0),abs.tol=1e-8)$value else d}

# Scaled periodic Gaussian
# Parameters: k (scale), s(bandwidth ≥0), tau (peak time) must all be vectors of length n
scaled.B <- function(t,Bp,i){
	Bp$k[i] * Gauss.B(t,list(s=Bp$s[i],tau=Bp$tau[i]))
}

# ------------------------ Generate metapopulation matrices -----------------------------

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


# ------------------------ Simulations -------------------------------------------------

# Return the list of parameters required from a simulation 
# Use scaled.B() for birth pulse
# mat = metapop contact matrix 
# p = list(span, s[], tau[], R0, IP, mu)
sim.par <- function(mat,p) {
	list(n=nrow(mat), contact.matrix=mat, 
	     B = scaled.B, Bp = list(k=sapply(1:nrow(mat),function(i) B.k(1/p$span,p$s[i])), s=p$s, tau=p$tau), 
	     d=1/p$span, beta=p$R0*(1/p$IP+1/p$span), gamma=1/p$IP, mu=p$mu)
}


# Use ssa.adaptivetau to run multiple simulations
# Parallel version
meta.SIR.ssa.tau <- function(n.simul, mat, par, init, t.end, tau.leap.par=list(epsilon=0.005)){
	foreach(1:n.simul, .inorder = F, .packages='adaptivetau') %dopar% {
		ssa.adaptivetau(init,meta.SIR.transitions(mat),meta.SIR.Rates,sim.par(mat,par),tf=t.end,tl.params = tau.leap.par)
	}
}

# ----------------------- Statistics from stochastic simulations -----------------------------------

# Thin the simulation output at specified time steps
# sim is the output from ssa.adaptivetau(), t is a vector of time steps
ssa.thin <- function(sim,t.vec){
	i <- 2
	t(sapply(t.vec,function(t){
		while(sim[i,1]<t) i <<- i+1
		c(Time=t,sim[i-1,-1])
	}))
}

# Calculate summary statistics about infection dynamics from a single simulation
# n = number of patches
# sim = output from ssa.adaptivetau
meta.SIR.ssa.stats <- function(sim,n){
	I.col <- (n+2):(2*n+1)
	# Matrix of presence-absence
	pres.mat <- sim[,I.col]>0
	# Global extinction
	all.ext.vec <- which(apply(pres.mat,1,sum)==0)
	all.ext.t <- if(length(all.ext.vec)==0) NA else sim[all.ext.vec[1],1]
	# 
	# RETURN
	return(list(pres.mat = pres.mat, ext = all.ext.t))
}

# Run multiple simulations and only save thinned time series
meta.SIR.ssa.thin <- function(n.simul, mat, par, init, t.end, dt = 0.01, tau.leap.par=list(epsilon=0.005)){
	foreach(1:n.simul, .inorder = F, .packages='adaptivetau') %dopar% {
		ssa.thin(ssa.adaptivetau(init,meta.SIR.transitions(mat),meta.SIR.Rates,sim.par(mat,par),tf=t.end,tl.params = tau.leap.par),seq(0,t.end,dt))
	}
}


# ---------------------- Deterministic model -------------------------------------------
# Use transition matrix and rate function to calculate the deterministic solution
# Use deSolve package

# Calculate the derivatives dX/dt = f(t,X,p)
stoch.diff <- function(t,var,par,trans.mat,rate.fun){
	list(trans.mat %*% rate.fun(var,par,t))
}

meta.SIR.ode <- function(meta.mat,par,init,t.end,dt=0.005){
	ode(init,seq(0,t.end,dt),stoch.diff,sim.par(meta.mat,par),trans.mat=meta.SIR.transitions(meta.mat),rate.fun=meta.SIR.Rates)
}

# Use the birth-death deterministic model to calculate the initial population size required to produce a given yearly-average population size
# Run for one year in a disconnected metapopulation with no infection 
meta.SIR.init <- function(par){
	n <- length(par$s)
	meta.mat <- matrix(0,n,n)
	init <- c(rep(1,n),rep(0,2*n))
	S <- meta.SIR.ode(meta.mat,par,init,1,dt=1e-3)[,2:(n+1)]
	if(n==1) 1/mean(S) else 1/apply(S,2,mean)
}

	