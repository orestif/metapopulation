####################################################
# 	METAPOPULATION MODEL FOR INFECTION DYNAMICS
#		   WITH ANNUAL BIRTH PULSES
# 		AND DENSITY-DEPENDENT DEATH RATE
#	------------------------------------------
#	TWO COMPETING STRAINS WITH CROSS-IMMUNITY
#	------------------------------------------
#			DEFINITIONS
####################################################

# Olivier Restif, July 2015

# Stochastic SIR model with seasonal Gaussian birth rate 
# Density-dependent death rate (DDD) introduced to prevent unlimited growth in the metapopulation.
# To allow comparison, a density-independent death rate is also included.

# The single-patch demographic model is:
# dN/dt = [b(t) - d - dd*N(t)] * N(t)/lambda
# Birth rate: b(t)*N(t)/lambda
# Death rate: [d + dd*N(t)] * N(t)/lambda
# where
# - b(t) is the seasonal birthing term with <b(t)> = 1. Annual births per capita = 1/lambda
# - K = 1/dd is the carrying capacity when d=0. Average population size <N> = (1-d)/dd
# - lambda is the average life span when N = <N>

# Infection dynamics in a single patch:
# dI/dt = beta*S*I/N - gamma*I - (d+dd*N)*I/lambda
# should lead to an average R0 = <N>/<S> = beta / (gamma + 1/lambda)

# Two pathogen strains with complete cross-immunity: S Ix Iy R
# Parameters: beta_x, beta_y, R0_x, R0_y, gamma_x, gamma_y

# R version using package `adaptivetau`

library(adaptivetau)
library(foreach)
library(parallel)
library(doParallel)
library(deSolve)

# ================================== MODEL DEFINITIONS ===================================

# Variables: {S_i,Ix_i,Iy_i,R_i} 1≤i≤np 
# These can be handled by names (using string manipulation) or by position
# Order of variables: S_1,..., S_n, Ix_1, ..., Ix_n, Iy_1, ..., Iy_n, R_1,... , R_n
# => S_i = var[i], Ix_i = var[i+np], Iy_i = var[i+2*np], R_i = var[i+3*np]

# Time unit = year

# ----------------------- Transition matrix -----------------------------------------------
# Create a transition matrix with 4*np rows (variables) and 9*np+4*k columns (events)
# where np=nrow(contact.matrix) is the number of patches and k=sum(contact.matrix) is the number of migration events
# Argument: contact.matrix: square matrix of size np with 0 and 1, M[i,j]=1 <=> migration from i to j
meta.DDD.SIR.2.transitions <- function(contact.matrix){
	np <- nrow(contact.matrix)
	var <- numeric(4*np)
	moves <- which(contact.matrix==1,T) # Matrix of all possible contacts
	cbind(
		# Births (np)
		sapply(1:np,function(i){replace(var,i,+1)}),
		# Deaths (4*np)
		sapply(1:(4*np),function(i){replace(var,i,-1)}),
		# Infection with x (np)
		sapply(1:np,function(i){replace(var,c(i,i+np),c(-1,+1))}),
		# Infection with y (np)
		sapply(1:np,function(i){replace(var,c(i,i+2*np),c(-1,+1))}),
		# Recovery from x (np)
		sapply(1:np,function(i){replace(var,c(i+np,i+3*np),c(-1,+1))}),
		# Recovery from y (np)
		sapply(1:np,function(i){replace(var,c(i+2*np,i+3*np),c(-1,+1))}),
		# Migrations (S)
		if(nrow(moves)>0) sapply(1:nrow(moves),function(j){replace(var,moves[j,],c(-1,+1))}) else NULL,
		# Migrations (Ix)
		if(nrow(moves)>0) sapply(1:nrow(moves),function(j){replace(var,moves[j,]+np,c(-1,+1))}) else NULL,
		# Migrations (Iy)
		if(nrow(moves)>0) sapply(1:nrow(moves),function(j){replace(var,moves[j,]+2*np,c(-1,+1))}) else NULL,
		# Migrations (R)
		if(nrow(moves)>0) sapply(1:nrow(moves),function(j){replace(var,moves[j,]+3*np,c(-1,+1))}) else NULL
	)
}

# ----------------------- Transition rates -----------------------------------------------
# Calculate rates and return a vector of length 6*np+3*k 
# SIR model with FD transmission within patches, all individuals die and reproduce equally, directional migration
# Arguments:
# - x: vector of variable values, no names required
# - par: list of parameters, with elements:
# 	- np: number of patches
#	- contact.matrix: square matrix of size np with 0 and 1, M[i,j]=1 <=> migration from i to j
#	- b() = function(t,p) returning a vector of normalised birth rates per capita: <b> = 1
# 	- b.par: list of parameters used by function b(t,p)
#	- d: density-independent death rate
#	- dd[]: density-dependent death rate (1/carrying capacity)
#	- lambda: average life-span
#	- beta[]: transmission rates for each strain
#	- gamma[]: recovery rates
#	- mu: migration rate
meta.DDD.SIR.2.Rates <- function(x,par,t) with(par,{
	moves <- which(contact.matrix==1,T) # Matrix of all possible contacts
	# Calculate population size in each patch 
	N <- x[1:np] + x[(np+1):(2*np)] + x[(2*np+1):(3*np)] + x[(3*np+1):(4*np)] 
	return(c(
		# Births: assume all individuals contribute
		b(t,b.par) * N/lambda,
		# Deaths: rate is independent of status or patch
		sapply(1:(4*np),function(i){
			p <- (i-1)%%np+1
			(d + N[p]*dd[p]) * x[i]/lambda }),
		# Infection: assume FD within patch -- test N>0 !!
		sapply(1:np,function(i){ if(N[i]>0) beta[1]*x[i]*x[i+np]/N[i] else 0}),
		sapply(1:np,function(i){ if(N[i]>0) beta[2]*x[i]*x[i+2*np]/N[i] else 0}),
		# Recovery
		gamma[1]*x[(np+1):(2*np)],
		gamma[2]*x[(2*np+1):(3*np)],
		# Migrations
		if(nrow(moves)>0){ 
			c(sapply(1:nrow(moves),function(j){ mu*x[moves[j,1]] }),
			  sapply(1:nrow(moves),function(j){ mu*x[moves[j,1]+np] }),
			  sapply(1:nrow(moves),function(j){ mu*x[moves[j,1]+2*np] }),
			  sapply(1:nrow(moves),function(j){ mu*x[moves[j,1]+3*np] }))
		} else NULL
	))
})


# ------------------------- Birth pulse function ----------------------------------------

# Periodic Gaussian, using a sin function so that tau = peak time, and normalised so that <B> = 1 
# Parameters: t = time, Bp = list(s(bandwidth ≥0), tau (peak time)), 
# If length(s) == length(tau) == np AND length(t) == 1, return a vector of length np
# WARNING: do not use with length(t) > 1 AND length(s)*length(tau) > 1

Gauss.b <- function(t,p){
	exp( -p$s * (sin(pi*(t-p$tau)))^2 ) / besselI(p$s/2,0,T)
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

# Generate contact matrices (square matrices of size np with 0 and 1) from a choice of types:
# - "circle": patch i linked to i-1 and i+1, plus a link between 1 and np.
# - "all": every patch linked to each other
meta.matrix <- function(np,type){
	switch(type,
		 circle = if(np==1) matrix(0,1) else if(np==2) matrix(c(0,1,1,0),2) else {
		 	x <- c(c(1),rep(0,np-3),c(1,0))
		 	sapply(1:np, function(i) {x <<- rotate.r(x)})
		 },
		 all = sapply(1:np,function(i) c(rep(1,i-1),0,rep(1,np-i)))
	)
}


# ------------------------ Simulations -------------------------------------------------

# Return the list of parameters required from a simulation 
# Use Gauss.b() for birth pulse
# mat = metapop contact matrix 
# p: ecological parameters = list(life.span, d, K[], s[], tau[], R0[], IP[], mu)
# Note: K can be Inf to set dd=0
sim.par <- function(mat,p) {
	list(np=nrow(mat), contact.matrix=mat, 
	     b = Gauss.b, b.par = list(s=p$s, tau=p$tau), 
	     d = p$d, dd = 1/p$K, lambda = p$life.span, 
	     beta = p$R0*(1/p$IP+1/p$life.span), gamma = 1/p$IP, mu = p$mu)
}

# THIS IS THE MAIN FUNCTION FOR THE USER
# Use ssa.adaptivetau to run multiple simulations
# Parallel version
# - n.simul: positive integer
# - mat: square matrix of 0s and 1s for connections between patches
# - par: ecological parameters = list(life.span, d, K[], s[], tau[], R0[], IP[], mu)
# - t.end: positive number
# - thin: save variables at regular time steps. If NA, save all time-points
# Return a list of tables

meta.DDD.SIR.2.ssa.tau <- function(n.simul, mat, par, init, t.end, thin=NA, tau.leap.par=list(epsilon=0.005)){
	foreach(1:n.simul, .inorder = F, .packages='adaptivetau') %dopar% 
	{
		sim <- ssa.adaptivetau(init,meta.DDD.SIR.2.transitions(mat),meta.DDD.SIR.2.Rates,sim.par(mat,par),tf=t.end,tl.params = tau.leap.par)
		if(is.finite(thin)){
			t.vec <- seq(0,t.end,thin)
			i <- 2
			t(sapply(t.vec,function(t){
				while(sim[i,1]<t) i <<- i+1
				c(Time=t,sim[i-1,-1])
			}))
		} else sim
	}
}


# ----------------------- Statistics from stochastic simulations -----------------------------------

# Thin a simulation output at specified time steps
# sim is the output from ssa.adaptivetau(), t is a vector of time steps
ssa.thin <- function(sim,t.vec){
	i <- 2
	t(sapply(t.vec,function(t){
		while(sim[i,1]<t) i <<- i+1
		c(Time=t,sim[i-1,-1])
	}))
}

# Calculate summary statistics about infection dynamics from a single simulation
# np = number of patches
# sim = output from ssa.adaptivetau
meta.DDD.SIR.2.ssa.stats <- function(sim,np){
	I.col <- (np+2):(3*np+1)
	# Matrix of presence-absence
	pres.mat <- sim[,I.col]>0
	# Global extinction
	all.ext.vec <- which(apply(pres.mat,1,sum)==0)
	all.ext.t <- if(length(all.ext.vec)==0) NA else sim[all.ext.vec[1],1]
	# 
	# RETURN
	return(list(pres.mat = pres.mat, ext = all.ext.t))
}


# ---------------------- Deterministic model -------------------------------------------
# Use transition matrix and rate function to calculate the deterministic solution
# Use deSolve package

# Calculate the derivatives dX/dt = f(t,X,p)
stoch.diff <- function(t,var,par,trans.mat,rate.fun){
	list(trans.mat %*% rate.fun(var,par,t))
}

meta.DDD.SIR.2.ode <- function(meta.mat,par,init,t.end,dt=0.005){
	ode(init,seq(0,t.end,dt),stoch.diff,sim.par(meta.mat,par),trans.mat=meta.DDD.SIR.2.transitions(meta.mat),rate.fun=meta.DDD.SIR.2.Rates)
}

# Use the birth-death deterministic model to calculate S(0) on a stable cycle in the absence of infection and migration.
meta.DDD.SIR.2.init <- function(par){
	np <- length(par$s)
	meta.mat <- matrix(0,np,np)
	foreach(i=1:np, .combine=c) %dopar%{
		optimize(function(x){
			init <- rep(0,4*np)
			init[i] <- x
			(meta.DDD.SIR.2.ode(meta.mat,par,init,t.end=1,dt=1)[2,i+1]-x)^2},
			par$K[i]*c(0.5,1.5),tol=1E-4)$minimum
	}
}

