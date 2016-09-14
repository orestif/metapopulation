####################################################
# 	METAPOPULATION MODEL FOR INFECTION DYNAMICS
# 		   WITH ANNUAL BIRTH PULSES
# 	       DENSITY-DEPENDENT DEATH RATE
#         AND MATERNAL TRANSFER OF ANTIBODIES
#	-------------------------------------------
#			RUN SOME TESTS
####################################################

source('Models/Metapop DDD MSIR adaptivetau.R')
library(tidyr)
library(dplyr)
library(ggplot2)

# ==================== Numerical values for test ===========================

# Single population
test.mat.1 <- matrix(c(0),1,1)
test.init.1 <- c(S1=990,I1=10,R1=0,M1=0)
test.par.1 <- list(life.span=1, d=0, K=1000, s=100, tau=0.5, R0=2, IP=0.1, mu=0.01, rho=1, MP=0.1)

# Two patches
test.mat.2 <- matrix(c(0,1,1,0),2,2)
test.init.2 <- c(S1=990,S2=1000,I1=10,I2=0,R1=0,R2=0, M1=0, M2=0)
test.par.2 <- list(life.span=1, d=0, K=c(1000,1000), s=c(100,10),tau=c(0.5,0),R0=2,IP=0.1,mu=0.01, rho=1, MP=0.1)


# ===================== RUN SIMULATIONS ==================================
# Run a single simulation 
system.time(test.sim.1 <- ssa.adaptivetau(test.init.1,meta.DDD.MSIR.transitions(test.mat.1),meta.DDD.MSIR.Rates,sim.par.DDD.MSIR(test.mat.1,test.par.1),tf=20,tl.params = list(epsilon=0.001)))

sim.1 <- gather(as.data.frame(test.sim.1),"Variable","N",-1)

ggplot(sim.1,aes(time,N)) + geom_line(aes(color=Variable))

system.time(test.sim.z <- meta.DDD.MSIR.ssa.tau(1,test.mat.1,test.par.1,test.init.1,t.end=10,thin=0.1,tau.leap.par = list(epsilon=0.001)))
sim.z <- gather(as.data.frame(test.sim.z),"Variable","N",-1)
ggplot(sim.z,aes(Time,N)) + geom_line(aes(color=Variable))


# Parallelisation using parallel package
registerDoParallel(cores=16) # Register the parallel backend to be used by foreach

# Run simulations in parallel
system.time(test.sim.2 <- meta.DDD.MSIR.ssa.tau(16,test.mat.1,test.par.1,test.init.1,thin=0.1,t.end=10))

# ======================= PLOTS ==============================

# Plot results
par(mfrow=c(1,1))
plot(0,xlim=c(0,10),xlab='Time (y)', ylim=c(0,200*length(test.sim.2)),ylab='Infected',yaxt='n',type='n')
for(i in 1:length(test.sim.2)){
	abline(h=200*(i-1))
	lines(test.sim.2[[i]][,1],test.sim.2[[i]][,5]+200*(i-1),col='blue')
	lines(test.sim.2[[i]][,1],test.sim.2[[i]][,3]+200*(i-1),col='red')
}

# ================== DETERMINISTIC MODEL =========================

meta.ode.sol.1 <- meta.DDD.MSIR.ode(test.mat.1,test.par.1,test.init.1,t.end=20)
plot(meta.ode.sol.1)

par(mfrow=c(1,1))
plot(meta.ode.sol.1[,4],meta.ode.sol.1[,5],type='l')

meta.DDD.SIR.init(test.par)
