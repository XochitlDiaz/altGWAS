#JUN 2023
#Xochitl Diaz
### Graphs created for the tests



setwd("./Desktop/StatGeno/altGWAS")

#Functions setup

#create summary table
mksum <- function(sim){
  
  cols= c("afreq", "meffect", "sdeffect","skew", "pow.lreg", "pow_rq.t50", 
          "pow_rq.t75", "pow_rq.t95", "pow_rq.t99", "pow.quail", "pow.bartlett",
          "pow.bf", "pow.levene", "pow.bp")
  summ <- data.frame(matrix(nrow = 0,ncol = length(cols)) )
  colnames(summ)= cols
  
  return(summ)
}


### Distribution characterization on different effect sizes ###
#Power on regression models at minor allele frequency = 0.5
plot.summ <- function(summ, title){
  
  test_names= c("LR", "Tau.50", "Tau.75", "Tau.95", "Tau.99", "QUAIL","Bartlett","Levene", "BF", "BP")
  test_col= c("red","#2A595D","#22577A","#56C9D4","#38A3A5","#E8BB00","#C1ED8F",
              "#7A4EBF","#F289CA","#ACBF15")
  plot(summ$sdeffect,summ$pow.lreg, 
       xlab= "Effect on sd", ylab= "Power", type="l", ylim=c(0,1), col=test_col[1], lwd=2,
       main= title)
  
  lines(summ$sdeffect, summ$pow_rq.t50, 
        xlab= "sdeffect", ylab= "tau.50 power", type="l", col= test_col[2],lwd=2)
  
  lines(summ$sdeffect, summ$pow_rq.t75, 
        xlab= "sdeffect", ylab= "tau.75 power",type="l", col=test_col[3],lwd=2)
  
  
  lines(summ$sdeffect, summ$pow_rq.t95, 
        xlab= "sdeffect", ylab= "tau.95 power", col= test_col[4],lwd=2)
  
  lines(summ$sdeffect, summ$pow_rq.t99, 
        xlab= "sdeffect", ylab= "tau.99 power", col= test_col[5],lwd=2)
  
  lines(summ$sdeffect, summ$pow.quail, 
        xlab= "sdeffect", ylab= "Quail power", col= test_col[6],lwd=2)
  
  lines(summ$sdeffect, summ$pow.bartlett, 
        xlab= "sdeffect", ylab= "Bartlett power", col= test_col[7],lwd=2)
  
  lines(summ$sdeffect, summ$pow.levene, 
        xlab= "sdeffect", ylab= "Levene power", col= test_col[8],lwd=2)
  
  lines(summ$sdeffect, summ$pow.bf, 
        xlab= "sdeffect", ylab= "BF power", col= test_col[9],lwd=2)
  
  lines(summ$sdeffect, summ$pow.bp, 
        xlab= "sdeffect", ylab= "BP power", col= test_col[10],lwd=2)
  
  
  lines(x=seq(0, max(summ$sdeffect),0.1), y=rep(0.05,length(seq(0, max(summ$sdeffect),0.1))), lwd=1.5, lty=2)
  
  legend( c(1,0.6),legend = test_names, fill=test_col , box.lty = 0)
}


### !skew meffect=0 ###
sim= read.csv("../batch5/rebatch1/run2/rebatch1/noskewm0/sim_noskewm0.csv")
#summ= mksum(sim = sim)
summ= mksum(sim)
afrq=0.5
m=0
sk=0

#this_model= subset(sim, sim$meffect==m & sim$sdeffect==sd & sim$skew==sk)
this_sub= subset(sim, sim$meffect==m & sim$skew==sk & sim$sdeffect <= 0.8)
for (sd in unique(this_sub$sdeffect)){
  this_model= subset(this_sub, this_sub$sdeffect==sd)
  nsim=nrow(this_model)
  summ[nrow(summ)+1, ]= c(afrq, m, sd , sk,
                          length(which(this_model$lreg_pv < 0.05))/nsim,
                          length(which(this_model$rq_pval.t50 < 0.05))/nsim,
                          length(which(this_model$rq_pval.t75 < 0.05))/nsim,
                          length(which(this_model$rq_pval.t95 < 0.05))/nsim,
                          length(which(this_model$rq_pval.t99 < 0.05))/nsim,
                          length(which(this_model$quail_pval  < 0.05))/nsim,
                          length(which(this_model$bartlett_pval  < 0.05))/nsim,
                          length(which(this_model$bf_pval  < 0.05))/nsim,
                          length(which(this_model$levene_pval  < 0.05))/nsim,
                          length(which(this_model$bp_pval  < 0.05))/nsim
                        )
}


plot.summ(summ, title= "meffect:0 ; skew:no")


### !skew meffect=0.2 ###
sim= read.csv("../batch5/rebatch1/run2/rebatch1/noskewm2/sim_noskewm2.csv")
#summ= mksum(sim = sim)
summ= mksum(sim)
afrq=0.5
m=0.02
sk=0

#this_model= subset(sim, sim$meffect==m & sim$sdeffect==sd & sim$skew==sk)
this_sub= subset(sim, sim$meffect==m & sim$skew==sk & sim$sdeffect <= 0.8)
for (sd in unique(this_sub$sdeffect)){
  this_model= subset(this_sub, this_sub$sdeffect==sd)
  nsim=nrow(this_model)
  summ[nrow(summ)+1, ]= c(afrq, m, sd , sk,
                          length(which(this_model$lreg_pv < 0.05))/nsim,
                          length(which(this_model$rq_pval.t50 < 0.05))/nsim,
                          length(which(this_model$rq_pval.t75 < 0.05))/nsim,
                          length(which(this_model$rq_pval.t95 < 0.05))/nsim,
                          length(which(this_model$rq_pval.t99 < 0.05))/nsim,
                          length(which(this_model$quail_pval  < 0.05))/nsim,
                          length(which(this_model$bartlett_pval  < 0.05))/nsim,
                          length(which(this_model$bf_pval  < 0.05))/nsim,
                          length(which(this_model$levene_pval  < 0.05))/nsim,
                          length(which(this_model$bp_pval  < 0.05))/nsim
  )
}


plot.summ(summ, title= "meffect:0.02 ; skew:no")


### !skew meffect=4 ###
sim= read.csv("../batch5/rebatch1/run2/rebatch1/noskewm4/sim_noskewm4.csv")
#summ= mksum(sim = sim)
summ= mksum(sim)
afrq=0.5
m=0.04
sk=0

#this_model= subset(sim, sim$meffect==m & sim$sdeffect==sd & sim$skew==sk)
this_sub= subset(sim, sim$meffect==m & sim$skew==sk & sim$sdeffect <= 0.8)
for (sd in unique(this_sub$sdeffect)){
  this_model= subset(this_sub, this_sub$sdeffect==sd)
  nsim=nrow(this_model)
  summ[nrow(summ)+1, ]= c(afrq, m, sd , sk,
                          length(which(this_model$lreg_pv < 0.05))/nsim,
                          length(which(this_model$rq_pval.t50 < 0.05))/nsim,
                          length(which(this_model$rq_pval.t75 < 0.05))/nsim,
                          length(which(this_model$rq_pval.t95 < 0.05))/nsim,
                          length(which(this_model$rq_pval.t99 < 0.05))/nsim,
                          length(which(this_model$quail_pval  < 0.05))/nsim,
                          length(which(this_model$bartlett_pval  < 0.05))/nsim,
                          length(which(this_model$bf_pval  < 0.05))/nsim,
                          length(which(this_model$levene_pval  < 0.05))/nsim,
                          length(which(this_model$bp_pval  < 0.05))/nsim
  )
}


plot.summ(summ, title= "meffect:0.04 ; skew:no")


### skew meffect=0 ###
sim= read.csv("../batch5/rebatch1/run2/rebatch1/skewm0/sim_skewm0.csv")
#summ= mksum(sim = sim)
summ= mksum(sim)
afrq=0.5
m=0
sk=3

#this_model= subset(sim, sim$meffect==m & sim$sdeffect==sd & sim$skew==sk)
this_sub= subset(sim, sim$meffect==m & sim$skew==sk & sim$sdeffect <= 0.5)
for (sd in unique(this_sub$sdeffect)){
  this_model= subset(this_sub, this_sub$sdeffect==sd)
  nsim=nrow(this_model)
  summ[nrow(summ)+1, ]= c(afrq, m, sd , sk,
                          length(which(this_model$lreg_pv < 0.05))/nsim,
                          length(which(this_model$rq_pval.t50 < 0.05))/nsim,
                          length(which(this_model$rq_pval.t75 < 0.05))/nsim,
                          length(which(this_model$rq_pval.t95 < 0.05))/nsim,
                          length(which(this_model$rq_pval.t99 < 0.05))/nsim,
                          length(which(this_model$quail_pval  < 0.05))/nsim,
                          length(which(this_model$bartlett_pval  < 0.05))/nsim,
                          length(which(this_model$bf_pval  < 0.05))/nsim,
                          length(which(this_model$levene_pval  < 0.05))/nsim,
                          length(which(this_model$bp_pval  < 0.05))/nsim
  )
}


plot.summ(summ, title= "meffect:0 ; skew:yes")


### skew meffect=0.02 ###
sim= read.csv("../batch5/rebatch1/run2/rebatch1/skewm2/sim_skewm2.csv")
#summ= mksum(sim = sim)
summ= mksum(sim)
afrq=0.5
m=0.02
sk=3

#this_model= subset(sim, sim$meffect==m & sim$sdeffect==sd & sim$skew==sk)
this_sub= subset(sim, sim$meffect==m & sim$skew==sk & sim$sdeffect <= 0.5 )
for (sd in unique(this_sub$sdeffect)){
  this_model= subset(this_sub, this_sub$sdeffect==sd)
  nsim=nrow(this_model)
  summ[nrow(summ)+1, ]= c(afrq, m, sd , sk,
                          length(which(this_model$lreg_pv < 0.05))/nsim,
                          length(which(this_model$rq_pval.t50 < 0.05))/nsim,
                          length(which(this_model$rq_pval.t75 < 0.05))/nsim,
                          length(which(this_model$rq_pval.t95 < 0.05))/nsim,
                          length(which(this_model$rq_pval.t99 < 0.05))/nsim,
                          length(which(this_model$quail_pval  < 0.05))/nsim,
                          length(which(this_model$bartlett_pval  < 0.05))/nsim,
                          length(which(this_model$bf_pval  < 0.05))/nsim,
                          length(which(this_model$levene_pval  < 0.05))/nsim,
                          length(which(this_model$bp_pval  < 0.05))/nsim
  )
}


plot.summ(summ, title= "meffect:0.02 ; skew:yes")


### skew meffect=0.04 ###
sim= read.csv("../batch5/rebatch1/run2/rebatch1/skewm4/sim_skewm4.csv")
#summ= mksum(sim = sim)
summ= mksum(sim)
afrq=0.5
m=0.04
sk=3

#this_model= subset(sim, sim$meffect==m & sim$sdeffect==sd & sim$skew==sk)
this_sub= subset(sim, sim$meffect==m & sim$skew==sk & sim$sdeffect <= 0.5 )
for (sd in unique(this_sub$sdeffect)){
  this_model= subset(this_sub, this_sub$sdeffect==sd)
  nsim=nrow(this_model)
  summ[nrow(summ)+1, ]= c(afrq, m, sd , sk,
                          length(which(this_model$lreg_pv < 0.05))/nsim,
                          length(which(this_model$rq_pval.t50 < 0.05))/nsim,
                          length(which(this_model$rq_pval.t75 < 0.05))/nsim,
                          length(which(this_model$rq_pval.t95 < 0.05))/nsim,
                          length(which(this_model$rq_pval.t99 < 0.05))/nsim,
                          length(which(this_model$quail_pval  < 0.05))/nsim,
                          length(which(this_model$bartlett_pval  < 0.05))/nsim,
                          length(which(this_model$bf_pval  < 0.05))/nsim,
                          length(which(this_model$levene_pval  < 0.05))/nsim,
                          length(which(this_model$bp_pval  < 0.05))/nsim
  )
}

plot.summ(summ, title= "meffect:0.04 ; skew:yes")






### Characterization of power on different minor allele frequencies ###

plot.summ <- function(summ, title){
  
  test_names= c("LR", "tau.50", "tau.75", "tau.95", "tau.99", "QUAIL","Bartlett","Levene", "BF", "BP")
  test_col= c("red","#2A595D","#22577A","#56C9D4","#38A3A5","#E8BB00","#C1ED8F",
              "#7A4EBF","#F289CA","#ACBF15")
  plot(summ$afreq,summ$pow.lreg, 
       xlab= "Minor allele frequency", ylab= "Power", type="l", ylim=c(0,1), col=test_col[1], lwd=2,
       main= title)
  
  lines(summ$afreq, summ$pow_rq.t50, 
        xlab= "afreq", ylab= "tau.50 power", type="l", col= test_col[2],lwd=2)
  
  lines(summ$afreq, summ$pow_rq.t75, 
        xlab= "afreq", ylab= "tau.75 power",type="l", col=test_col[3],lwd=2)
  
  
  lines(summ$afreq, summ$pow_rq.t95, 
        xlab= "afreq", ylab= "tau.95 power", col= test_col[4],lwd=2)
  
  lines(summ$afreq, summ$pow_rq.t99, 
        xlab= "afreq", ylab= "tau.99 power", col= test_col[5],lwd=2)
  
  lines(summ$afreq, summ$pow.quail, 
        xlab= "afreq", ylab= "Quail power", col= test_col[6],lwd=2)
  
  lines(summ$afreq, summ$pow.bartlett, 
        xlab= "afreq", ylab= "Bartlett power", col= test_col[7],lwd=2)
  
  lines(summ$afreq, summ$pow.levene, 
        xlab= "afreq", ylab= "Levene power", col= test_col[8],lwd=2)
  
  lines(summ$afreq, summ$pow.bf, 
        xlab= "afreq", ylab= "BF power", col= test_col[9],lwd=2)
  
  lines(summ$afreq, summ$pow.bp, 
        xlab= "afreq", ylab= "BP power", col= test_col[10],lwd=2)
  
  
  lines(x=seq(0, max(summ$afreq),0.01), y=rep(0.05,length(seq(0, max(summ$afreq),0.01))), lwd=1.5, lty=2)
  
  legend( c(1,0.6),legend = test_names, fill=test_col , box.lty = 0)
}



### skew meffect=0.00 ###
sim= read.csv("../batch5/rebatch3/simulations_afrqchangem0_table.csv")
#summ= mksum(sim = sim)
summ= mksum(sim)
m=0
sk=0
sd=0

#this_model= subset(sim, sim$meffect==m & sim$sdeffect==sd & sim$skew==sk)
#this_sub= subset(sim, sim$meffect==m & sim$skew==sk & sim$sdeffect <= 0.5 )
for (afrq in unique(sim$afreq)){
  this_model= subset(sim, sim$afreq==afrq)
  nsim=nrow(this_model)
  summ[nrow(summ)+1, ]= c(afrq, m, sd , sk,
                          length(which(this_model$lreg_pv < 0.05))/nsim,
                          length(which(this_model$rq_pval.t50 < 0.05))/nsim,
                          length(which(this_model$rq_pval.t75 < 0.05))/nsim,
                          length(which(this_model$rq_pval.t95 < 0.05))/nsim,
                          length(which(this_model$rq_pval.t99 < 0.05))/nsim,
                          length(which(this_model$quail_pval  < 0.05))/nsim,
                          length(which(this_model$bartlett_pval  < 0.05))/nsim,
                          length(which(this_model$bf_pval  < 0.05))/nsim,
                          length(which(this_model$levene_pval  < 0.05))/nsim,
                          length(which(this_model$bp_pval  < 0.05))/nsim
  )
}

plot.summ(summ, title= "meffect:0 ; sdeffect:0 ; skew:no")



### skew meffect=0.2 ###
sim= read.csv("../batch5/rebatch3/simulations_afrqchangem2_table.csv")
#summ= mksum(sim = sim)
summ= mksum(sim)
m=0.2
sk=0
sd=0

#this_model= subset(sim, sim$meffect==m & sim$sdeffect==sd & sim$skew==sk)
#this_sub= subset(sim, sim$meffect==m & sim$skew==sk & sim$sdeffect <= 0.5 )
for (afrq in unique(sim$afreq)){
  this_model= subset(sim, sim$afreq==afrq)
  nsim=nrow(this_model)
  summ[nrow(summ)+1, ]= c(afrq, m, sd , sk,
                          length(which(this_model$lreg_pv < 0.05))/nsim,
                          length(which(this_model$rq_pval.t50 < 0.05))/nsim,
                          length(which(this_model$rq_pval.t75 < 0.05))/nsim,
                          length(which(this_model$rq_pval.t95 < 0.05))/nsim,
                          length(which(this_model$rq_pval.t99 < 0.05))/nsim,
                          length(which(this_model$quail_pval  < 0.05))/nsim,
                          length(which(this_model$bartlett_pval  < 0.05))/nsim,
                          length(which(this_model$bf_pval  < 0.05))/nsim,
                          length(which(this_model$levene_pval  < 0.05))/nsim,
                          length(which(this_model$bp_pval  < 0.05))/nsim
  )
}

plot.summ(summ, title= "meffect:0.2 ; sdeffect:0 ; skew:no")








### P value distribution on a single point QQunif ###
library(gap)

unique(sim$meffect)
unique(sim$sdeffect)
unique(sim$skew)

m= 0
sd= 0
sk= 0

this_model= subset(sim, sim$meffect== m & sim$sdeffect==sd & sim$skew== sk)

par(mfrow=c(2,5), mar=c(6,4,7,1))
qqunif(this_model$lreg_pv, main= "LR", col = test_col[1])
qqunif(this_model$rq_pval.t50, main= "Tau.50", col = test_col[2])
qqunif(this_model$rq_pval.t75, main= "Tau.75", col = test_col[3])
qqunif(this_model$rq_pval.t95, main= "Tau.95", col = test_col[4])
qqunif(this_model$rq_pval.t99, main= "Tau.99", col = test_col[5])
qqunif(this_model$quail_pval, main= "QUAIL", col = test_col[6])
qqunif(this_model$bartlett_pval, main= "Bartlett", col = test_col[7])
qqunif(this_model$bf_pval, main= "BF", col = test_col[8])
qqunif(this_model$levene_pval, main= "Levene", col = test_col[9])
qqunif(this_model$bp_pval, main= "BP", col = test_col[10])
mtext("meffect:0 ; sdeffect:0 ; skew:no", side = 3, line= -2, outer = T, cex= 1.5,font=2)




### Run Time ###
sim= read.csv("../batch5/runtime/simulations_runtime_table.csv")
trunt= read.csv("../batch5/runtime/overall_runtime.csv")
summrunt= data.frame( "Test"= c("LR", "Tau.50", "Tau.75","Tau.95","Tau.99",
                                "QUAIL","Bartlett","Levene","BF","BP"),
                      "runtime"= c(sum(sim$lreg_time), sum(sim$t50_time),
                                 sum(sim$t95_time), sum(sim$t99_time),
                                 sum(sim$t99_time), sum(sim$quail_time),
                                 sum(sim$bartlett_time), sum(sim$levene_time), 
                                 sum(sim$bf_time), sum(sim$bp_time))
                    )


top=max(summrunt$runtime)
bottom=min(summrunt$runtime)
library(RColorBrewer)
library("ggplot2")
library("ggbreak")



plot=ggplot(summrunt, aes(x= factor(Test, test_names), y= runtime/bottom, fill= Test))+ 
  geom_bar(stat = "identity", fill=c("red","#2A595D", "#22577A", "#56C9D4", "#38A3A5", "#E8BB00",
                                     "#C1ED8F", "#7A4EBF", "#F289CA", "#ACBF15")) + 
  geom_text(aes(label= paste(round(runtime, digits = 4), "s", sep = ""), vjust = - 0.6)) +
  theme(text = element_text(size = 12)) +
  ylim(0,392) +
  scale_y_break(c(4.7,390)) +
  ylab("Runtime proportion of LR") +
  xlab("Test")
plot





