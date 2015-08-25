####################################################
# 	METAPOPULATION MODEL FOR INFECTION DYNAMICS
#		WITH ANNUAL BIRTH PULSES
# 	  AND DENSITY-DEPENDENT DEATH RATE
#	-------------------------------------------
#	    SERIES 1 - DETERMINISTIC SOLUTIONS
# 		    EXTENDED PARAMETER SET
####################################################

# RAN ON OR226 17/07

library(foreach)
registerDoParallel(cores=16)

source('Metapop DDD SIR adaptivetau.R')

def.par <- list(life.span=1,d=0,K=1000,s=100,tau=0.5,R0=4,IP=1/12,mu=0.01)

# Explore a range of values for s, tau, mu, IP and span 
np.try <- c(2,4,6,8)
tau.max.try <- c(0,0.25,0.5)
s.try <- c(0,1,10,100)
mu.try <- c(1E-4,1E-3,1E-2,1E-1)
IP.try <- c(7/365, 1/12, 1/2)
span.try <- c(1,2,5,10)

par.try <- expand.grid(np.try, tau.max.try, s.try, mu.try, IP.try, span.try)
colnames(par.try) <- c('Patches','tau.max','s', 'mu', 'IP', 'life.span')

# ===================== SOLVE ODES ==================================

ode.DDD.series.1 <- foreach(i = 1:nrow(par.try)) %dopar% {
	par.i <- def.par
	np <- par.try$Patches[i]

	par.i$life.span <- par.try$life.span[i]
	par.i$K <- rep(def.par$K,np)
	par.i$s <- as.numeric(rep(par.try$s[i],np))
	# varies tau linearly from 0.5 in patch 1 to 0.5-tau.max in the opposite patch
	par.i$tau <- 0.5 - (par.try$tau.max[i] * (1-abs(1-(0:(np-1))*(2/np))))
	par.i$IP <- par.try$IP[i]
	par.i$mu <- par.try$mu[i]
	print(par.try[i,])
	init.i <- c(meta.DDD.SIR.init(par.i), 1E-4, rep(0,np-1), rep(0,np))
	names(init.i) <- c(paste('S',1:np,sep=''),paste('I',1:np,sep=''),paste('R',1:np,sep=''))
	meta.DDD.SIR.ode(meta.matrix(np,'circle'), par.i, init.i, t.end=20, dt=0.02)
}

save(ode.DDD.series.1, par.try,file='DDD_Series_1_ODE.RData')
