####################################################
# 	METAPOPULATION MODEL FOR INFECTION DYNAMICS
#		WITH ANNUAL BIRTH PULSES
#	-------------------------------------------
#		    SIMULATION SERIES
####################################################

source('Metapop SIR adaptivetau.R')

# ==================== Numerical Values ===========================

# Run simulations with 2 patches 
# Parameters taken from Peel et al (2014), Fig. 2b, with a CCS around 10,000.
mat.1 <- matrix(c(0,1,1,0),2,2)
init.1 <- c(4990,5000,10,0,0,0)

par.1 <- list(span=1,s=c(10,10),tau=c(0.5,0.5),R0=4,IP=1/12,mu=0.001)

# Explore a range of values for s, tau, mu and N
s.try <- c(0,1,10,100)
tau2.try <- c(0,0.25,0.5)
mu.try <- c(0,1e-4,1e-3,1e-2)
N.try <- c(1000,5000,10000)

par.try <- expand.grid(s.try,tau2.try,mu.try,N.try)
colnames(par.try) <- c('s','tau','mu','N')

# ===================== RUN SIMULATIONS ==================================

registerDoParallel(cores=15) # Register the parallel backend to be used by foreach

# Run series and return extinction times
sim.series.2 <- sapply(1:nrow(par.try),function(i) {
	par.i <- par.1
	par.i$s <- rep(par.try[i,'s'],2)
	par.i$tau <- c(0.5,par.try[i,'tau'])
	par.i$mu <- par.try[i,'mu']
	print(par.try[i,])
	scale.N <- meta.SIR.init(par.i)
	init.i <- c(S1=round(scale.N[1]*par.try[i,'N']-1), S2=round(scale.N[2]*par.try[i,'N']), I1=1, I2=0, R1=0, R2=0)
	
	sims <- meta.SIR.ssa.thin(n.simul=5000, mat.1, par.i, init.i, t.end=20, tau.leap.par=list(epsilon=0.02), dt=0.01)

	save(par.i,init.i,sims,file=paste('Series_2/Sims_',i,'.RData',sep=''))
	sapply(sims, function(sim){
		ext <- which(sim[,4]+sim[,5]==0)
		if(length(ext)==0) NA else sim[ext[1],1]
	})
})

write.csv(sim.series.2,file='Series_2/Extinctions.csv')

save.image('Metapop Series 2.RData')
