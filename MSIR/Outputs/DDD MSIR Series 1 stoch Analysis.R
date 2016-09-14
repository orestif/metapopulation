####################################################
# 	METAPOPULATION MODEL FOR INFECTION DYNAMICS
# 		   WITH ANNUAL BIRTH PULSES
# 	       DENSITY-DEPENDENT DEATH RATE
#         AND MATERNAL TRANSFER OF ANTIBODIES
#	-------------------------------------------
#	ANALYSIS OF RESULTS FROM SIMULATION SERIES 1
####################################################

source('Models/Metapop DDD MSIR adaptivetau.R')
load('Outputs/DDD_MSIR_Stoch_Series_1/DDD_MSIR_Stoch_Series_1.RData')

library(dplyr)
library(tidyr)
library(ggplot2)

# --------------- Analyse extinctions ----------------------------

p.ext <- apply(DDD.MSIR.stoch.series.1,2, function(x) length(which(is.finite(x)))/200)
t.ext.mean <- apply(DDD.MSIR.stoch.series.1,2, mean, na.rm=T)
t.ext.median <- apply(DDD.MSIR.stoch.series.1,2, median, na.rm=T)

MSIR.series.ext <- cbind(par.try,p.ext=p.ext,t.ext.mean=t.ext.mean,t.ext.median=t.ext.median)


# Heatmap for P(ext)
ggplot(MSIR.series.ext,aes(factor(s),factor(IP))) + geom_tile(aes(fill=p.ext)) + facet_grid(rho~MP, labeller=label_both) + labs(fill='P(extinction)') + scale_fill_continuous(low=rgb(1,1,0.6),high=rgb(0.9,0,0)) + xlab("s") + ylab("IP") + theme(text=element_text(size=18))

MSIR.series.ext %>% mutate(IP.day=round(IP*365),MP.day=round(MP*365)) %>% ggplot(aes(factor(MP.day),factor(rho))) + geom_tile(aes(fill=p.ext)) + facet_grid(s~IP.day, labeller=label_both) + labs(fill='P(extinction)') + scale_fill_continuous(low=rgb(1,1,0.6),high=rgb(0.9,0,0)) + xlab("Average duration of maternal immunity (days)") + ylab("rho") + theme(text=element_text(size=18))


# Heatmap for <T(ext)>
ggplot(MSIR.series.ext,aes(factor(s),factor(IP))) + geom_tile(aes(fill=t.ext.mean)) + facet_grid(rho~MP, labeller=label_both) + labs(fill='T(extinction)') + scale_fill_continuous(low=rgb(1,1,0.6),high=rgb(0.9,0,0)) + xlab("s") + ylab("IP") + theme(text=element_text(size=18))


