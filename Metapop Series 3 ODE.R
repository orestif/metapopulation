####################################################
# 	METAPOPULATION MODEL FOR INFECTION DYNAMICS
#		WITH ANNUAL BIRTH PULSES
#	-------------------------------------------
#	    SERIES 3 - DETERMINISTIC SOLUTIONS
####################################################



par.3 <- list(span=1,s=100,tau=0.5,R0=4,IP=1/12,mu=0.01)

# Explore a range of values for s, tau, mu and N
meta.try <- c(2,4,6,8)
tau.max.try <- c(0,0.25,0.5)
N.try <- c(1000,5000,10000)

par.try <- expand.grid(meta.try,tau.max.try,N.try)
colnames(par.try) <- c('Patches','tau.max','N')

# ===================== SOLVE ODES ==================================

ode.series.3 <- lapply(1:nrow(par.try),function(i) {
	par.i <- par.3
	np <- par.try[i,'Patches']
	par.i$s <- as.numeric(rep(par.3['s'],np))
	# varies tau linearly from 0.5 in patch 1 to 0.5-tau.max in the opposite patch
	par.i$tau <- 0.5 - (par.try[i,'tau.max'] * (1-abs(1-(0:(np-1))*(2/np))))
	print(par.i)
	scale.N <- meta.SIR.init(par.i)
	init.i <- c(round(scale.N*par.try[i,'N']-1), 1, rep(0,np-1), rep(0,np))
	names(init.i) <- c(paste('S',1:np,sep=''),paste('I',1:np,sep=''),paste('R',1:np,sep=''))
	
	# sims <- meta.SIR.ssa.thin(n.simul=5000, meta.matrix(np,'circle'), par.i, init.i, t.end=20, tau.leap.par=list(epsilon=0.02), dt=0.01)
	meta.SIR.ode(meta.matrix(np,'circle'), par.i, init.i, t.end=20, dt=0.01)
})

save(ode.series.3, par.try,file='Series_3_ODE.RData')
