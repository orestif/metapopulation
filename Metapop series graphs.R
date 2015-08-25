####################################################
# 	METAPOPULATION MODEL FOR INFECTION DYNAMICS
#		WITH ANNUAL BIRTH PULSES
#	-------------------------------------------
#     Visual representation of simulation results
####################################################

library(ggplot2)
library(dplyr)
library(tidyr)


# ======================================= SIR - SERIES 2 =======================================

# TWO PATCHES, R0=4, IP = 1/12 YEARS, LIFE SPAN = 1YEAR

load('Metapop Series 2.RData')

# ---------------- Probability of extinction within 20 years --------------------------

# Data formatting
sim.2.results <- as.data.frame(par.try)
sim.2.results$ext.NA <- apply(sim.series.2,2,function(x) length(which(is.na(x))))
sim.2.results$phase <- 12*(0.5-sim.2.results$tau)

# Heatmap
pdf("~/Documents/Work/Research/Bats/Presentations/Heatmap Series 2.pdf",12,6)
ggplot(sim.2.results,aes(factor(mu),factor(phase))) + geom_tile(aes(fill=1-ext.NA/5000)) + facet_grid(N~s, labeller=label_both) + labs(fill='P(extinction)') + scale_fill_continuous(low=rgb(1,1,0.6),high=rgb(0.9,0,0)) + xlab("Migration rate per year") + ylab("Birth pulse lag (months)") + theme(text=element_text(size=18))
dev.off()

# ------------- Time to extinction -----------------------------

# Data formatting
rownames(sim.series.2) <- 1:nrow(sim.series.2)
colnames(sim.series.2) <- rep('Time',ncol(sim.series.2))
sim.2.all <- cbind(sim.2.results,t(sim.series.2))
sim.2.ext <- gather(sim.2.all,Simulation,Extinction,6:5005,convert=T)

# Boxplots
ggplot(sim.2.ext,aes(interaction(factor(mu),factor(tau)),Extinction)) + geom_violin() + facet_grid(N~s, labeller=label_both)

ggplot(filter(sim.2.ext,mu==0.01 & s>2),aes(factor(tau),Extinction,fill=1-ext.NA/5000)) + geom_violin() + coord_flip(ylim=c(0,5)) + facet_grid(N~s, labeller=label_both)+ labs(fill='P(extinction)') + scale_fill_continuous(low='darkred',high='white')


# ======================================= SIR - SERIES 3 =======================================

rm(list=ls())
load('Metapop Series 3.RData')

# Data formatting
sim.3.results <- as.data.frame(par.try)
sim.3.results$ext.NA <- apply(sim.series.3,2,function(x) length(which(is.na(x))))
rownames(sim.series.3) <- 1:nrow(sim.series.3)
colnames(sim.series.3) <- rep('Time',ncol(sim.series.3))
sim.3.all <- cbind(sim.3.results,t(sim.series.3))
sim.3.ext <- gather(sim.3.all,Simulation,Extinction,5:5004,convert=T)

sim.3.results$Total.Size <- sim.3.results$Patches*sim.3.results$N
sim.3.results$Phase <- 12*sim.3.results$tau.max
sim.3.results$CI.min <- qbinom(0.025,5000,1-sim.3.results$ext.NA/5000)/5000
sim.3.results$CI.max <- qbinom(0.975,5000,1-sim.3.results$ext.NA/5000)/5000

ggplot(filter(sim.3.results,Total.Size==20000),aes(factor(Patches),factor(tau.max))) + geom_tile(aes(fill=1-ext.NA/5000)) + labs(fill='P(extinction)') + scale_fill_continuous(low=rgb(1,1,0.6),high=rgb(0.9,0,0))

pdf("~/Documents/Work/Research/Bats/Presentations/Fragmentation Series 3.pdf",6,6)
ggplot(filter(sim.3.results,Total.Size==20000),aes(factor(Phase),1-ext.NA/5000,fill=factor(Patches))) + geom_bar(stat='identity',position="dodge") + labs(fill='Patches') + ylab("Probability of extinction") + xlab("Birth pulse lag (months)") + theme(text=element_text(size=18)) + geom_errorbar(aes(ymin=CI.min,ymax=CI.max),position=position_dodge(0.9),width=0.2,col="yellow") + scale_fill_manual(values=c("blue","blue4"))
dev.off()

# Individual simulations
load('~/Documents/Work/Research/Bats/Metapopulation/Series_3_I_list.RData')

pdf("~/Documents/Work/Research/Bats/Presentations/Sim Series 3 in phase.pdf",8,4)
ggplot(filter(sim_3_list[[with(par.try,which(Patches==2 & tau.max==0 & N==5000))]], Simulation==12), aes(Time,I)) + geom_line(aes(color=factor(p)),size=1) + scale_color_manual(values=c("orange","purple")) + ylab("Infected individuals") +  theme(text=element_text(size=18)) + coord_cartesian(ylim=c(0,2500),xlim=c(0,16))
dev.off()

pdf("~/Documents/Work/Research/Bats/Presentations/Sim Series 3 out of phase.pdf",9,4)
ggplot(filter(sim_3_list[[with(par.try,which(Patches==2 & tau.max==0.5 & N==5000))]], Simulation==9), aes(Time,I)) + geom_line(aes(color=factor(p)),size=1) + scale_color_manual(values=c("orange","purple")) + ylab("Infected individuals") +  theme(text=element_text(size=18)) + coord_cartesian(ylim=c(0,2500))
dev.off()

ggplot(filter(sim_3_list[[with(par.try,which(Patches==2 & tau.max==0 & N==5000))]], Simulation==12), aes(Time,p)) + geom_tile(aes(fill=factor(ifelse(I>0,1,0))))+ scale_fill_manual(values=c("black","white"))

ggplot(filter(sim_3_list[[with(par.try,which(Patches==2 & tau.max==0.5 & N==5000))]], Simulation==9), aes(Time,p)) + geom_tile(aes(fill=factor(ifelse(I>0,1,0))))+ scale_fill_manual(values=c("black","white"))


# ====================================== DDD SIR 2 SERIES 3 =======================================

load("Metapop DDD SIR 2 Series 3 Extinction.RData")

SIR.2.series.3.P.ext$CI.min <- qbinom(0.025,100,SIR.2.series.3.P.ext$P.ext)/100
SIR.2.series.3.P.ext$CI.max <- qbinom(0.975,100,SIR.2.series.3.P.ext$P.ext)/100
SIR.2.series.3.P.ext$Strain.names <- ifelse(SIR.2.series.3.P.ext$Strain=='x',"Acute","Slow")


pdf("~/Documents/Work/Research/Bats/Presentations/DDD SIR 2 Ext.pdf",6,4)
ggplot(SIR.2.series.3.P.ext, aes(factor(s),P.ext,fill=Strain.names)) + geom_bar(stat="identity",position='dodge') + facet_grid(Patches~., labeller=label_both) + labs(fill="Strain") + ylab("Probability of extinction") + xlab("Birth pulse tightness (s)") + theme(text=element_text(size=18)) + geom_errorbar(aes(ymin=CI.min,ymax=CI.max),position=position_dodge(0.9),width=0.2) + scale_fill_manual(values=c("red","cyan3"))
dev.off()


load("Metapop DDD SIR 2 Series 3 Occupancy.RData")
DDD.SIR.2.series.3.occupancy$Strain.names <- ifelse(DDD.SIR.2.series.3.occupancy$Strain=='x',"Acute","Slow")

pdf("~/Documents/Work/Research/Bats/Presentations/DDD SIR 2 Ocp.pdf",7,5)
ggplot(DDD.SIR.2.series.3.occupancy %>% filter(Patches==12), aes(Time,Avg)) + geom_line(aes(col=Strain.names),size=1.5) + facet_grid(~s, labeller=label_both)+ labs(color="Strain")+ scale_color_manual(values=c("red","cyan3"))+ theme(text=element_text(size=18)) + ylab("Average number of infected patches")
dev.off()


load('DDD_SIR_2_Series_3_I_list.RData')

pdf("~/Documents/Work/Research/Bats/Presentations/DDD SIR 2 Sim s10 p12.pdf",14,5)
ggplot(DDD.SIR.2.series.3.list[[4]] %>% filter(Simulation==11 & p<8),aes(Time,I)) + geom_line(aes(col=Strain),size=1) + facet_grid(.~p) + coord_cartesian(ylim=c(0,2000))+ theme(text=element_text(size=18)) + scale_color_manual(values=c("red","cyan3"))+ylab("Number of infected individuals")
dev.off()

