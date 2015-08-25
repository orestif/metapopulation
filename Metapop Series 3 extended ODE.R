####################################################
# 	METAPOPULATION MODEL FOR INFECTION DYNAMICS
#		WITH ANNUAL BIRTH PULSES
#	-------------------------------------------
#	    SERIES 3 - DETERMINISTIC SOLUTIONS
# 		    EXTENDED PARAMETER SET
####################################################

# Note: when this file was executed, some solutions ended early, especially with short IP. 
# Is this a problem with default numerical precision in ode()?

library(foreach)

source('Metapop SIR adaptivetau.R')

par.3 <- list(span=1,s=100,tau=0.5,R0=4,IP=1/12,mu=0.01)

# Explore a range of values for s, tau, mu, IP and span 
meta.try <- c(2,4,6,8)
tau.max.try <- c(0,0.25,0.5)
s.try <- c(0,1,10,100)
mu.try <- c(1E-4,1E-3,1E-2,1E-1,1)
IP.try <- c(7/365, 1/12, 1/2)
span.try <- c(1,2,5,10)

par.try <- expand.grid(meta.try, tau.max.try, s.try, mu.try, IP.try, span.try)
colnames(par.try) <- c('Patches','tau.max','s', 'mu', 'IP', 'span')

# ===================== SOLVE ODES ==================================

ode.series.3.extended <- foreach(i = 1:nrow(par.try)) %dopar% {
	par.i <- par.3
	np <- par.try$Patches[i]

	par.i$span <- par.try$span[i]
	par.i$s <- as.numeric(rep(par.try$s[i],np))
	# varies tau linearly from 0.5 in patch 1 to 0.5-tau.max in the opposite patch
	par.i$tau <- 0.5 - (par.try$tau.max[i] * (1-abs(1-(0:(np-1))*(2/np))))
	par.i$IP <- par.try$IP[i]
	par.i$mu <- par.try$mu[i]
	print(par.try[i,])
	scale.N <- meta.SIR.init(par.i)
	init.i <- c(scale.N, 1E-4, rep(0,np-1), rep(0,np))
	names(init.i) <- c(paste('S',1:np,sep=''),paste('I',1:np,sep=''),paste('R',1:np,sep=''))
	
	# sims <- meta.SIR.ssa.thin(n.simul=5000, meta.matrix(np,'circle'), par.i, init.i, t.end=20, tau.leap.par=list(epsilon=0.02), dt=0.01)
	meta.SIR.ode(meta.matrix(np,'circle'), par.i, init.i, t.end=20, dt=0.02)
}

save(ode.series.3.extended, par.try,file='Series_3_extended_ODE.RData')
