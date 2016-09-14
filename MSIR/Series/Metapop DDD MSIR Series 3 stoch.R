####################################################
# 	METAPOPULATION MODEL FOR INFECTION DYNAMICS
# 		   WITH ANNUAL BIRTH PULSES
# 	       DENSITY-DEPENDENT DEATH RATE
#         AND MATERNAL TRANSFER OF ANTIBODIES
#	-------------------------------------------
#		    SIMULATION SERIES 3
####################################################

# In this series I vary:
# - birth pulse tightness s
# - carrying capacity K
# - Probability of MTA rho
# - Duration of maternal immunity MP


# Constants:
# - One patch
# - All simulations initiated with 5 cases to avoid initial fade-outs.
# - life.span=2, tau=0.5, R0=4, d=0, IP=1/12
# Run 1,000 simulations for 10 years.

source('Models/Metapop DDD MSIR adaptivetau.R')

library(foreach)
registerDoParallel(cores=32) # Register the parallel backend to be used by foreach

# Run on serversofia 26/11  


# ==================== Numerical Values ===========================

def.par <- list(life.span=2,d=0,K=1000,s=100,tau=0.5,R0=4,IP=1/12,rho=1,MP=0.1,mu=0)

# Explore a range of values for s, IP, rho, MP and K
par.list <- list(s=c(0,1,2,5,10,20,50,100), rho=c(0,1), MP=c(1/12,1/6,1/3), K=100*2^c(2:8))
par.try <- do.call(expand.grid,par.list)

# Remove duplicate rows with no MTA (rho=0) and variable MP
par.try <- par.try[-which(par.try$rho==0 & par.try$MP>min(par.list$MP)),]


# ===================== RUN SIMULATIONS ==================================

# Run series and return extinction times
DDD.MSIR.stoch.series.3 <- sapply(1:nrow(par.try),function(i) {
	par.i <- replace.par(def.par,par.try[i,])
	np <- 1
	print(par.try[i,])
	init.i <- c(S1=round(meta.DDD.MSIR.init(par.i)-1), I1=5L, R1=0, M1=0)
	sims <- meta.DDD.MSIR.ssa.tau(n.simul=1000, meta.matrix(np,'circle'), par.i, init.i, t.end=10, thin=0.05, tau.leap.par=list(epsilon=0.01))

	save(par.i,init.i,sims,file=paste('DDD_MSIR_Stoch_Series_3/Sims_',i,'.RData',sep=''))
	ext.sim <- sapply(sims, function(sim){
		ext <- which(sim[,3]==0)
		if(length(ext)==0) NA else sim[ext[1],1]
	})
	print(paste("Extinctions:",length(which(is.finite(ext.sim)))))
	return(ext.sim)
})

write.csv(DDD.MSIR.stoch.series.3,file='DDD_MSIR_Stoch_Series_3/Extinctions.csv')

save.image('DDD_MSIR_Stoch_Series_3/DDD_MSIR_Stoch_Series_3.RData')
