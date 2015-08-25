####################################################
# 	METAPOPULATION MODEL FOR INFECTION DYNAMICS
#		WITH ANNUAL BIRTH PULSES
#	-------------------------------------------
#		  ANALYSIS OF SERIES 3
####################################################

# Olivier Restif, 3 July 2015

source('Metapop SIR adaptivetau.R')
load('Metapop Series 3.RData')
library(dplyr)
library(tidyr)

sim.dir <- '~/Documents/Repository/Metapopulation/Series_3'

registerDoParallel(cores=15) # Register the parallel backend to be used by foreach

# Unless specified otherwise, the code sections below should only be run once, as the results are saved in text files.

# ---------------------- Time to extinction -----------------------------
# DONE - DO NOT RE-RUN

# Extract pathogen extinction times for ach simulation (as this was not correctly assessed at run time)

rm(sim.series.3)

sim.series.3 <- foreach(i=1:nrow(par.try), .combine=rbind) %dopar% {
	np <- par.try[i,'Patches']
	load(file.path(sim.dir,paste('Sims_',i,'.RData',sep='')))
	ex <- sapply(sims, function(sim){
		ext <- which(apply(sim[,(np+2):(2*np+1)],1,sum)==0)
		if(length(ext)==0) NA else sim[ext[1],1]
	})
	rm(sim,par.i,init.i)
	ex
}
write.csv(t(sim.series.3),file='~/Documents/Repository/Metapopulation/Series_3/Extinctions.csv')

# NOTE: The corrected version of sim.series.3 has been added to Metapop Series 3.RData


# ----------------------- Spatiotemporal dynamics --------------------------

# Record average (across simulations) number of patches occupied at every time step
# DONE
series.3.mean.occupancy <- foreach(i=1:nrow(par.try), .combine=rbind) %dopar% {
	np <- par.try[i,'Patches']
	I.col <- (np+2):(2*np+1)
	load(file.path(sim.dir,paste('Sims_',i,'.RData',sep='')))
	ocp <- sapply(sims, function(sim){
		apply(sim[,I.col],1,function(x) length(which(x>0)))
	})
	apply(ocp,1,mean)
}
write.csv(series.3.mean.occupancy,file=file.path(sim.dir,'Mean_Occupancy.csv'),row.names = F)

# Average occupancy conditional on non-extinction at any time
# DONE
series.3.cond.mean.occupancy <- foreach(i=1:nrow(par.try), .combine=rbind) %dopar% {
	np <- par.try[i,'Patches']
	I.col <- (np+2):(2*np+1)
	load(file.path(sim.dir,paste('Sims_',i,'.RData',sep='')))
	ocp <- sapply(sims, function(sim){
		apply(sim[,I.col],1,function(x) length(which(x>0)))
	})
	apply(ocp,1,function(x) mean(x[x>0]))
}
write.csv(series.3.cond.mean.occupancy,file=file.path(sim.dir,'Cond_Mean_Occupancy.csv'),row.names = F)

# Distribution of number of patches occupied - thin time series to 10 points per year
# Done
series.3.dist.occupancy <- foreach(i=1:nrow(par.try), .combine=rbind) %dopar% {
	np <- par.try[i,'Patches']
	I.col <- (np+2):(2*np+1)
	load(file.path(sim.dir,paste('Sims_',i,'.RData',sep='')))
	ocp <- sapply(sims, function(sim){
		apply(sim[seq(1,2001,10),I.col],1,function(x) length(which(x>0)))
	})
	apply(ocp,1,function(x) sapply(0:np, function(i) length(which(x==i))))
}
colnames(series.3.dist.occupancy) <- seq(0,20,0.1)
series.3.dist.heads <- foreach(i=1:nrow(par.try), .combine=rbind) %do% {
	cbind(sapply(par.try[i,], rep, times=par.try$Patches[i]+1),Occupied=0:par.try$Patches[i])
}
write.csv(cbind(series.3.dist.heads,series.3.dist.occupancy),file=file.path(sim.dir,'Dist_Occupancy.csv'),row.names = F)


# Look at individual simulations from tau.max=0.25, 8 patches, N=5000 (index )
# NOTE: the sims object is ~2GB
# DONE
load(file.path(sim.dir,paste('Sims_',which(with(par.try,Patches==8 & tau.max==0.25 & N==5000)),'.RData',sep='')))
# Gather thinned I_n time-series into a single table
sim.I.wide <- data.frame(Patches=8,tau.max=0.25,N=5000,Simulation=rep(1:5000,each=8),p=rep(1:8,5000))
I.col <- (10:17)
thins <- foreach(x=sims, .combine=rbind) %dopar% t(x[seq(1,2001,10),I.col])
colnames(thins) <- seq(0,20,0.1)
sim.I.wide <- cbind(sim.I.wide,thins)
sim.I.table <- gather(sim.I.wide,"Time",'I',6:206,convert = T)
rm(sims,thins,sim.I.wide)
save(sim.I.table,file='Series_3_I_table_8_025_5000.RData')
# Note: the RData file for one set of parameter values is 13 MB


# Generate a single table with 100 thinned time series per parameter combination, 
# for the purpose of easy visualisation of individual simulations, e.g. with Shiny
# DONE (note the parallel code is not optimised)
# However, subsetting from this very large table is slow
sim_3_sample <- foreach(i=1:nrow(par.try),.combine=rbind) %dopar% {
	load(file.path(sim.dir,paste('Sims_',i,'.RData',sep='')))
	# Gather thinned I_n time-series into a single table
	np <- par.try$Patches[i]
	sim.I.wide <- data.frame(Patches=par.try$Patches[i],tau.max=par.try$tau.max[i],N=par.try$N[i],Simulation=rep(1:100,each=np),p=rep(1:np,100))
	I.col <- (np+2):(2*np+1)
	thins <- foreach(x=sims[1:100], .combine=rbind) %dopar% t(x[seq(1,2001,10),I.col])
	colnames(thins) <- seq(0,20,0.1)
	sim.I.wide <- cbind(sim.I.wide,thins)
	gather(sim.I.wide,"Time",'I',6:206,convert = T)
}
save(sim_3_sample,file='Series_3_I_table.RData')
rm(sim_3_sample)


# Save in a list instead, ordered according to par.try
# DONE - MUCH BETTER 
sim_3_list <- foreach(i=1:nrow(par.try)) %dopar% {
	load(file.path(sim.dir,paste('Sims_',i,'.RData',sep='')))
	# Gather thinned I_n time-series into a single table
	np <- par.try$Patches[i]
	sim.I.wide <- data.frame(Simulation=rep(1:100,each=np),p=rep(1:np,100))
	I.col <- (np+2):(2*np+1)
	thins <- foreach(x=sims[1:100], .combine=rbind) %do% t(x[seq(1,2001,10),I.col])
	colnames(thins) <- seq(0,20,0.1)
	sim.I.wide <- cbind(sim.I.wide,thins)
	gather(sim.I.wide,"Time",'I',3:203,convert = T)
}

save(par.try,sim_3_list,file='Series_3_I_list.RData')


# Calculate pairwise correlation coefficients for I across patches
# Save a list: for each parameter combination, create a matrix with Patches rows and Simulations columns
# For each simulation, the column shows correlation coefficients with Patch 1
# calculated for t > 1 year and t < t_ext
# DONE
sim_3_cor <- foreach(i=1:nrow(par.try)) %dopar% {
	load(file.path(sim.dir,paste('Sims_',i,'.RData',sep='')))
	sapply(sims, function(sim){
		np <- par.try$Patches[i]
		I.col <- (np+2):(2*np+1)
		I_mat <- sim[-(1:100),I.col]
		sel <- which(apply(I_mat,1,sum)>0)
		if(length(sel)>10) cor(I_mat[sel,])[,1] else rep(NA,np)
	})
}
save(par.try,sim_3_cor,file='Series_3_cor.RData')

