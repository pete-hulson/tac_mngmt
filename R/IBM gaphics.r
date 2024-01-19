rm(list=(ls()))
library(ggplot2)
library(lattice)
library(latticeExtra)
library(gridExtra)
require(plotrix); require(PBSmodelling); require(lattice); require(PBSadmb); require(fBasics); 
require(gplots); readADopts("Adopts.txt");require(Hmisc);require(calibrate);require(stringr);require(reshape2);

#setwd('F:/NOAA FILES/Research/Incorporating IBM data/Model/Final Application/Final/Adult and Larval Movement')
#IBM.est<- readList("rs_ibm.rep"); unpackList(IBM.est);
#IBM.est.recruitsE_BM<-recruitsE_BM
#IBM.est.R_aveE_BM<-R_aveE
#IBM.R_aveE1<-c(rep(R_aveE,(nyrs)))
#IBM.est.recruitsW_BM<-recruitsW_BM
#IBM.est.R_aveW_BM<-R_aveW
#IBM.R_aveW1<-c(rep(R_aveW,(nyrs)))
#IBM.est.recruitsE_AM<-recruitsE_AM
#IBM.est.R_aveE_AM<-R_aveE
#IBM.R_aveE1<-c(rep(R_aveE,(nyrs)))
#IBM.est.recruitsW_AM<-recruitsW_AM
#IBM.est.R_aveW_AM<-R_aveW
#IBM.R_aveW1<-c(rep(R_aveW,(nyrs)))
#IBM.abund.movE.W<-recruitsE_BM*TE_W[,1]
#IBM.abund.movW.E<-recruitsW_BM*TW_E[,1]

#setwd('F:/NOAA FILES/Research/Incorporating IBM data/Model/Final Application/Final/Adult Movement')
#spatial.est<- readList("rs_ibm.rep"); unpackList(spatial.est);
#spatial.est.recruitsE_BM<-recruitsE_BM
#spatial.est.R_aveE_BM<-R_aveE
#spatial.R_aveE1<-c(rep(R_aveE,(nyrs)))
#spatial.est.recruitsW_BM<-recruitsW_BM
#spatial.est.R_aveW_BM<-R_aveW
#spatial.R_aveW1<-c(rep(R_aveW,(nyrs)))
#spatial.est.recruitsE_AM<-recruitsE_AM
#spatial.est.R_aveE_AM<-R_aveE
#spatial.R_aveE1<-c(rep(R_aveE,(nyrs)))
#spatial.est.recruitsW_AM<-recruitsW_AM
#spatial.est.R_aveW_AM<-R_aveW
#spatial.R_aveW1<-c(rep(R_aveW,(nyrs)))

#setwd('F:/NOAA FILES/Research/Incorporating IBM data/Model/Final Application/Final/Nonspatial')
#nonspatial<-readList("rs_ibm.rep"); unpackList(nonspatial);
#nonspatial.est.recruitsE_AM<-recruitsE_AM
#nonspatial.est.R_aveE_AM<-R_aveE
#nonspatial.R_aveE1<-c(rep(R_aveE,(nyrs)))
#nonspatial.est.recruitsW_AM<-recruitsW_AM
#nonspatial.est.R_aveW_AM<-R_aveW
#nonspatial.R_aveW1<-c(rep(R_aveW,(nyrs)))
#nonspatial.est.recruitsE_BM<-recruitsE_BM
#nonspatial.est.R_aveE_BM<-R_aveE
#nonspatial.R_aveE1<-c(rep(R_aveE,(nyrs)))
#nonspatial.est.recruitsW_BM<-recruitsW_BM
#nonspatial.est.R_aveW_BM<-R_aveW
#nonspatial.R_aveW1<-c(rep(R_aveW,(nyrs)))



pdf("RS IBM.pdf")

IBM <- readList("rs_ibm.rep"); unpackList(IBM);

ESTIMATES<-readLines("rs_ibm.std",n=-1)
ESTIMATES<-unlist(strsplit(ESTIMATES," "))
ESTIMATES<-matrix(ESTIMATES)
ESTIMATES<-ESTIMATES[ESTIMATES!=""]
ESTIMATES<-ESTIMATES[c(5:length(ESTIMATES))]
ESTIMATES<-matrix(ESTIMATES,byrow=TRUE,ncol=4)
ESTIMATES<-ESTIMATES[,c(2:4)]
ESTIMATES<-data.frame(ESTIMATES)
colnames(ESTIMATES)<-c("Parameter","Value","SD")
ESTIMATES$SD<-as.numeric(as.character(ESTIMATES$SD))
ESTIMATES$Value<-as.numeric(as.character(ESTIMATES$Value))
CV<-rep(NA,times=length(ESTIMATES))
CV<-ESTIMATES$SD/ESTIMATES$Value
ESTIMATES$CV<-CV
Color<-ifelse(abs(ESTIMATES$CV)>0.2,"red","black")
Color<-cbind(Color,Color,Color,Color)
Color<-as.matrix(Color)
ifelse(file.exists("rs_ibm.cor"),converged<-"Yes",converged<-"No")


################################### rss and gradient
Pars<-readLines("rs_ibm.par",n=1)
Pars<-unlist(strsplit(Pars," "))
Pars<-matrix(Pars)
Pars<-Pars[6]

TpenE.larval<-sum(TpenE_larval+MpenE_larval+TpenE_larval_pre_IBM+MpenE_larval_pre_IBM)
TpenW.larval<-sum(TpenW_larval+MpenW_larval+TpenW_larval_pre_IBM+MpenW_larval_pre_IBM)

layout(matrix(c(1,1,1,1,1,1,2,2), 8, 1, byrow = TRUE))
rss.main<-c(rss_landings_HL_E,rss_landings_HL_W,rss_landings_LL_E,rss_landings_LL_W,rss_landings_MRIP_E,rss_landings_MRIP_W,rss_landings_HBT_E,rss_landings_HBT_W,rss_landings_SHR_E,rss_landings_SHR_W,
rss_index_HL_E,rss_index_HL_W,rss_index_MRIP_E,rss_index_MRIP_W,rss_index_HBT_E,rss_index_HBT_W,rss_effort_SHR_E,rss_effort_SHR_W,rss_index_SUMM_E,rss_index_SUMM_W,rss_index_FALL_E,rss_index_FALL_W,
rss_catch_prop_HL_E,rss_catch_prop_HL_W,rss_catch_prop_LL_E,rss_catch_prop_LL_W,rss_catch_prop_MRIP_E,rss_catch_prop_MRIP_W,rss_catch_prop_HBT_E,rss_catch_prop_HBT_W,rss_catch_prop_SHR_E,rss_catch_prop_SHR_W,
rss_index_prop_SUMM_E,rss_index_prop_SUMM_W,rss_index_prop_FALL_E,rss_index_prop_FALL_W,rss_IBM_propE,rss_IBM_propW,
F_devs_pen_HL_E,F_devs_pen_HL_W,F_devs_pen_LL_E,F_devs_pen_LL_W,F_devs_pen_MRIP_E,F_devs_pen_MRIP_W,F_devs_pen_HBT_E,F_devs_pen_HBT_W,F_devs_pen_SHR_E,F_devs_pen_SHR_W,
F_pen_HL_E_sum,F_pen_HL_W_sum,F_pen_LL_E_sum,F_pen_LL_W_sum,F_pen_MRIP_E_sum,F_pen_MRIP_W_sum,F_pen_HBT_E_sum,F_pen_HBT_W_sum,F_pen_SHR_E_sum,F_pen_SHR_W_sum,
rec_devs_penE,rec_devs_penW,R_ave_penE,R_ave_penW,sigma_penE,sigma_penW,rss_recruitE,rss_recruitW,init_abund_penE,init_abund_penW,TpenE,TpenW,TpenE.larval,TpenW.larval)


barplot(rss.main,names.arg=c("HL_E","HL_W","LL_E","LL_W","MRIP_E","MRIP_W","HBT_E","HBT_W","SHR_E","SHR_W",
"HL_E","HL_W","MRIP_E","MRIP_W","HBT_E","HBT_W","SHR_E","SHR_W","SUMM_E","SUMM_W","FALL_E","FALL_W",
"HL_E","HL_W","LL_E","LL_W","MRIP_E","MRIP_W","HBT_E","HBT_W","SHR_E","SHR_W","SUMM_E","SUMM_W","FALL_E","FALL_W","IBM_E","IBM_W",
"F_DEVS_HL_E","F_DEVS_HL_W","F_DEVS_LL_E","F_DEVS_LL_W","F_DEVS_MRIP_E","F_DEVS_MRIP_W","F_DEVS_HBT_E","F_DEVS_HBT_W","F_DEVS_SHR_E","F_DEVS_SHR_W",
"F_PEN_HL_E","F_PEN_HL_W","F_PEN_LL_E","F_PEN_LL_W","F_PEN_MRIP_E","F_PEN_MRIP_W","F_PEN_HBT_E","F_PEN_HBT_W","F_PEN_SHR_E","F_PEN_SHR_W",
"REC_DEVS_E","REC_DEVS_W","R_AVE_E","R_AVE_W","SIGMA_REC_E","SIGMA_REC_W","REC_E","REC_W","PEN_ABUND_E","PEN_ABUND_W","PEN_T_E","PEN_T_W","PEN_T_LARVAL_E","PEN_T_LARVAL_W"),
          col=c(rep("red",times=10),rep("blue",times=6),rep("deeppink",times=2),rep("salmon",times=4),rep("green",times=10),rep("dodgerblue",times=4),rep("magenta",times=2),
rep("orange",times=10),rep("brown",times=10),rep("gray",times=2),rep("pink",times=2),rep("aquamarine",times=2),rep("Yellow",times=2),rep("black",times=2),
rep("Purple",times=4)),cex.names=.65,las=2)
legend("topright",ncol=2, c("Landings","Index (CPUE)","Effort (SHR)","Survey (CPUE)",
               "Catch Proportions","Survey Proportions","IBM Proportions","F_Devs Penalty","Max F Pen","Rec_Devs Penalty","R_AVE Penalty","Recruit Sigma Pen","Recruitment Penalty","Init Abund Dev Pen","Movement Penalty"), 
               col = c("red","blue","deeppink" ,"salmon","green","dodgerblue","magenta","orange","Brown","Gray","Pink","Aquamarine","Yellow","black","purple"),lwd=2,cex=1)
title("Component Contributions to Objective Function (RSS)",line=1)

summary<-matrix(c(max_gradient,f,Pars,converged),ncol=4,byrow=TRUE)  
colnames(summary)<-c("Maximum Gradient","Negative Log-Likelihood","# of Parameters","Model Converged")
summary<-as.table(summary)
textplot(summary,show.rownames=FALSE, cmar=8,valign="top",halign="center")


##################################### correlations ##################


################################# Population Trajectories ##############################################

layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))

matplot(yrs,recruitsE_BM,type='n',main='Estimated Recruitment\n Pre and Post Movement',xlab='Year',
ylab='Recruitment (1000s of fish)',xlim=c(min(yrs)-1,max(yrs)+1),
        ylim=c(0,max(TRUE_recruitsE_BM,TRUE_recruitsW_BM,recruitsE_BM,recruitsW_BM,recruitsE_AM,recruitsW_AM)+.2*max(TRUE_recruitsE_BM,TRUE_recruitsW_BM,recruitsE_BM,recruitsW_BM,recruitsE_AM,recruitsW_AM)))
grid()
matlines(yrs,recruitsE_BM,col='cyan',lty=2,lwd=2);
matlines(yrs,recruitsE_AM,col='blue',lwd=2);
matlines(yrs,recruitsW_BM,col='salmon',,lty=2,lwd=2);
matlines(yrs,recruitsW_AM,col='red',lwd=2);
matlines(yrs,TRUE_recruitsE_BM,col='black',lty=3,lwd=2);
matlines(yrs,TRUE_recruitsW_BM,col='grey',lty=3,lwd=2);
legend("top", ncol=2,c("True_East", "True West","East Before Mov.","West Before Mov.","East After Mov.","West After Mov." ), 
   col = c('black','grey','cyan','salmon','blue','red'),lty=c(3,3,2,2,1,1),lwd=2, cex=.7)

matplot(yrs,ssbE,type='n',main='Estimated Spawning Stock',
       xlab='Year',ylab='SSB (# Eggs)',xlim=c(min(yrs)-1,max(yrs)+1),
        ylim=c(0,max(TRUE_ssbE,TRUE_ssbW,ssbE,ssbW)+.2*max(ssbE,ssbW,TRUE_ssbE,TRUE_ssbW)))
grid()
matlines(yrs,ssbE,col='blue',lwd=2); 
matlines(yrs,ssbW,col='red',lwd=2); 
matlines(yrs,TRUE_ssbE,col='black',lty=3,lwd=2);
matlines(yrs,TRUE_ssbW,col='grey',lty=3,lwd=2);
legend("top", ncol=2,c("TRUE East", "TRUE West","East","West" ), col = c('black','grey','blue','red'),lwd=2, cex=.8)

matplot(yrs,biomassE_BM,type='n',main='Estimated Biomass',
       xlab='Year',ylab='Biomass (mt)',xlim=c(min(yrs)-1,max(yrs)+1),
        ylim=c(0,max(TRUE_biomassE_AM,TRUE_biomassW_AM,biomassE_BM,biomassW_BM,biomassE_AM,biomassW_AM)+.2*max(TRUE_biomassE_AM,TRUE_biomassW_AM,biomassE_BM,biomassW_BM,biomassE_AM,biomassW_AM))) 
grid()
matlines(yrs,biomassE_BM,col='cyan',lty=2,lwd=2);
matlines(yrs,biomassE_AM,col='blue',lwd=2);
matlines(yrs,biomassW_BM,col='salmon',,lty=2,lwd=2);
matlines(yrs,biomassW_AM,col='red',lwd=2);
matlines(yrs,TRUE_biomassE_AM,col='black',,lty=3,lwd=2);
matlines(yrs,TRUE_biomassW_AM,col='gray',lty=3,lwd=2);
legend("top", ncol=2,c("TRUE East","TRUE West","East Before Mov.","West Before Mov.","East After Mov.","West After Mov." ), col = c('black', 'grey','cyan','salmon','blue','red'),lty=c(3,3,2,2,1,1),lwd=2, cex=.8)

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))

biomass_ageE<-abundance_at_ageE_BM*weight_spawnE
biomass_ageW<-abundance_at_ageW_BM*weight_spawnW
biomass_ageE_TRUE<-TRUE_abundance_at_ageE_BM*weight_spawnE
biomass_ageW_TRUE<-TRUE_abundance_at_ageW_BM*weight_spawnW

matplot(ages,biomass_ageE[1,],type='n',main='Initial Biomass at Age East',
       xlab='Age',ylab='Biomass (mt)',ylim=c(0,max(biomass_ageE,biomass_ageE_TRUE)+.2*max(biomass_ageE,biomass_ageE_TRUE))) 
grid()
matlines(ages,biomass_ageE[1,],col='cyan',lty=1,lwd=2);
matlines(ages,biomass_ageE[2,],col='blue',lwd=2);
matlines(ages,biomass_ageE[3,],col='salmon',lty=1,lwd=2);
matlines(ages,biomass_ageE[4,],col='red',lwd=2);
matlines(ages,biomass_ageE[5,],col='green',lwd=2);
matlines(ages,biomass_ageE_TRUE[1,],col='cyan',lty=3,lwd=2);
matlines(ages,biomass_ageE_TRUE[2,],col='blue',lwd=2,lty=3);
matlines(ages,biomass_ageE_TRUE[3,],col='salmon',lty=3,lwd=2);
matlines(ages,biomass_ageE_TRUE[4,],col='red',lwd=2,lty=3);
matlines(ages,biomass_ageE_TRUE[5,],col='green',lwd=2,lty=3);
legend("top", ncol=3,c("1872","1873","1874","1875","1876","TRUE" ), col = c('cyan','blue','salmon','red','green',"black"),lty=c(1,1,1,1,1,3),lwd=2, cex=.8)


matplot(ages,biomass_ageE[128,],type='n',main='Terminal Biomass at Age East',
       xlab='Age',ylab='Biomass (mt)',ylim=c(0,max(biomass_ageE,biomass_ageE_TRUE)+.2*max(biomass_ageE,biomass_ageE_TRUE))) 
grid()
matlines(ages,biomass_ageE[138,],col='cyan',lty=1,lwd=2);
matlines(ages,biomass_ageE[139,],col='blue',lwd=2);
matlines(ages,biomass_ageE[140,],col='salmon',lty=1,lwd=2);
matlines(ages,biomass_ageE[141,],col='red',lwd=2);
matlines(ages,biomass_ageE[142,],col='green',lwd=2);
matlines(ages,biomass_ageE_TRUE[138,],col='cyan',lty=3,lwd=2);
matlines(ages,biomass_ageE_TRUE[139,],col='blue',lwd=2,lty=3);
matlines(ages,biomass_ageE_TRUE[140,],col='salmon',lty=3,lwd=2);
matlines(ages,biomass_ageE_TRUE[141,],col='red',lwd=2,lty=3);
matlines(ages,biomass_ageE_TRUE[142,],col='green',lwd=2,lty=3);
legend("top", ncol=3,c("2009","2010","2011","2012","2013","TRUE" ), col = c('cyan','blue','salmon','red','green','black'),lwd=2, lty=c(1,1,1,1,1,3),cex=.8)

matplot(ages,biomass_ageW[1,],type='n',main='Initial Biomass at Age West',
       xlab='Age',ylab='Biomass (mt)',ylim=c(0,max(biomass_ageW,biomass_ageW_TRUE)+.2*max(biomass_ageW,biomass_ageW_TRUE))) 
grid()
matlines(ages,biomass_ageW[1,],col='cyan',lty=1,lwd=2);
matlines(ages,biomass_ageW[2,],col='blue',lwd=2);
matlines(ages,biomass_ageW[3,],col='salmon',lty=1,lwd=2);
matlines(ages,biomass_ageW[4,],col='red',lwd=2);
matlines(ages,biomass_ageW[5,],col='green',lwd=2);
matlines(ages,biomass_ageW_TRUE[1,],col='cyan',lty=3,lwd=2);
matlines(ages,biomass_ageW_TRUE[2,],col='blue',lwd=2,lty=3);
matlines(ages,biomass_ageW_TRUE[3,],col='salmon',lty=3,lwd=2);
matlines(ages,biomass_ageW_TRUE[4,],col='red',lwd=2,lty=3);
matlines(ages,biomass_ageW_TRUE[5,],col='green',lwd=2,lty=3);
legend("top", ncol=3,c("1872","1873","1874","1875","1876","TRUE" ), col = c('cyan','blue','salmon','red','green',"black"),lty=c(1,1,1,1,1,3),lwd=2, cex=.8)

matplot(ages,biomass_ageW[138,],type='n',main='Terminal Biomass at Age West',
       xlab='Age',ylab='Biomass (mt)',ylim=c(0,max(biomass_ageW,biomass_ageW_TRUE)+.2*max(biomass_ageW,biomass_ageW_TRUE))) 
grid()
matlines(ages,biomass_ageW[138,],col='cyan',lty=1,lwd=2);
matlines(ages,biomass_ageW[139,],col='blue',lwd=2);
matlines(ages,biomass_ageW[140,],col='salmon',lty=1,lwd=2);
matlines(ages,biomass_ageW[141,],col='red',lwd=2);
matlines(ages,biomass_ageW[142,],col='green',lwd=2);
matlines(ages,biomass_ageW_TRUE[138,],col='cyan',lty=3,lwd=2);
matlines(ages,biomass_ageW_TRUE[139,],col='blue',lwd=2,lty=3);
matlines(ages,biomass_ageW_TRUE[140,],col='salmon',lty=3,lwd=2);
matlines(ages,biomass_ageW_TRUE[141,],col='red',lwd=2,lty=3);
matlines(ages,biomass_ageW_TRUE[142,],col='green',lwd=2,lty=3);
legend("top", ncol=3,c("2009","2010","2011","2012","2013","TRUE" ), col = c('cyan','blue','salmon','red','green','black'),lwd=2, lty=c(1,1,1,1,1,3),cex=.8)

##################### MOvement ############################################################################



layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))

matplot(ages,TE_W[1,],type='n',main='Estimated Start Year Annual Movement \nRate by Age',xlab='Age',
ylab='Annual Movement Proportion',xlim=c(min(ages)-1,max(ages)+1),
        ylim=c(0,max(TE_W,TW_E)+.2*max(TE_W,TW_E)))

grid()
matlines(ages,TE_W[1,],col='blue',lwd=2); 
matlines(ages,TW_E[1,],col='red',lwd=2); 
matlines(ages,TRUE_TE_W[1,],col='blue',lwd=2,lty=3); 
matlines(ages,TRUE_TW_E[1,],col='red',lwd=2,lty=3); 
legend("top", ncol=2,c("East to West","West to East","TRUE" ), col = c('blue','red','black'),lty=c(1,1,3),lwd=2, cex=.8)


matplot(ages,TE_W[nyrs,],type='n',main='Estimated End Year Annual Movement \nRate by Age',xlab='Age',
ylab='Annual Movement Proportion',xlim=c(min(ages)-1,max(ages)+1),
        ylim=c(0,max(TE_W,TW_E)+.2*max(TE_W,TW_E)))

grid()
matlines(ages,TE_W[nyrs,],col='blue',lwd=2); 
matlines(ages,TW_E[nyrs,],col='red',lwd=2); 
matlines(ages,TRUE_TE_W[nyrs,],col='blue',lwd=2,lty=3); 
matlines(ages,TRUE_TW_E[nyrs,],col='red',lwd=2,lty=3); 
legend("top", ncol=2,c("East to West","West to East","TRUE" ), col = c('blue','red','black'),lty=c(1,1,3),lwd=2, cex=.8)

bio_move_recE<-abund_moveE[,1]*weight_spawnE[,1]
bio_move_recW<-abund_moveW[,1]*weight_spawnW[,1]
bio_move_recE_TRUE<-TRUE_abund_moveE[,1]*weight_spawnE[,1]
bio_move_recW_TRUE<-TRUE_abund_moveW[,1]*weight_spawnW[,1]

matplot(yrs,bio_moveE,type='n',main='Biomass Emigrating',
       xlab='Year',ylab='Biomass (mt)',xlim=c(min(yrs),max(yrs)+1),
ylim=c(0,max(bio_moveE,TRUE_bio_moveE,TRUE_bio_moveW,bio_moveW)
+.4*max(bio_moveE,bio_moveW)))
grid()
matlines(yrs,bio_moveE,col="blue",lwd=2,lty=1); 
matlines(yrs,bio_moveW,col="red",lwd=2,lty=1); 
matlines(yrs,bio_move_recE,col="cyan",lwd=2,lty=2); 
matlines(yrs,bio_move_recW,col="salmon",lwd=2,lty=2); 
matlines(yrs,TRUE_bio_moveE,col="blue",lwd=2,lty=3); 
matlines(yrs,TRUE_bio_moveW,col="red",lwd=2,lty=3); 
matlines(yrs,bio_move_recE_TRUE,col="cyan",lwd=2,lty=3); 
matlines(yrs,bio_move_recW_TRUE,col="salmon",lwd=2,lty=3); 
legend("top", ncol=2,c("East to West (Larval)","West to East (Larval)","East to West (Total)","West to East (Total)" ), col = c('cyan','salmon','blue','red'),lwd=2,lty=c(2,2,1,1), cex=.8)



############################ F and Selectivity ##########################################
fmult.HL.E<-c(rep(0,times=(nyrs-nyrs_landings_HL_E)),fmult_HL_E)
fmult.HL.W<-c(rep(0,times=(nyrs-nyrs_landings_HL_W)),fmult_HL_W)
fmult.LL.E<-c(rep(0,times=(nyrs-nyrs_landings_LL_E)),fmult_LL_E)
fmult.LL.W<-c(rep(0,times=(nyrs-nyrs_landings_LL_W)),fmult_LL_W)
fmult.MRIP.E<-c(rep(0,times=(nyrs-nyrs_landings_MRIP_E)),fmult_MRIP_E)
fmult.MRIP.W<-c(rep(0,times=(nyrs-nyrs_landings_MRIP_W)),fmult_MRIP_W)
fmult.HBT.E<-c(rep(0,times=(nyrs-nyrs_landings_HBT_E)),fmult_HBT_E)
fmult.HBT.W<-c(rep(0,times=(nyrs-nyrs_landings_HBT_W)),fmult_HBT_W)
fmult.SHR.E<-c(rep(0,times=(nyrs-nyrs_landings_SHR_E)),fmult_SHR_E)
fmult.SHR.W<-c(rep(0,times=(nyrs-nyrs_landings_SHR_W)),fmult_SHR_W)

TRUE_fmult.HL.E<-c(rep(0,times=(nyrs-nyrs_landings_HL_E)),TRUE_fmult_HL_E)
TRUE_fmult.HL.W<-c(rep(0,times=(nyrs-nyrs_landings_HL_W)),TRUE_fmult_HL_W)
TRUE_fmult.LL.E<-c(rep(0,times=(nyrs-nyrs_landings_LL_E)),TRUE_fmult_LL_E)
TRUE_fmult.LL.W<-c(rep(0,times=(nyrs-nyrs_landings_LL_W)),TRUE_fmult_LL_W)
TRUE_fmult.MRIP.E<-c(rep(0,times=(nyrs-nyrs_landings_MRIP_E)),TRUE_fmult_MRIP_E)
TRUE_fmult.MRIP.W<-c(rep(0,times=(nyrs-nyrs_landings_MRIP_W)),TRUE_fmult_MRIP_W)
TRUE_fmult.HBT.E<-c(rep(0,times=(nyrs-nyrs_landings_HBT_E)),TRUE_fmult_HBT_E)
TRUE_fmult.HBT.W<-c(rep(0,times=(nyrs-nyrs_landings_HBT_W)),TRUE_fmult_HBT_W)
TRUE_fmult.SHR.E<-c(rep(0,times=(nyrs-nyrs_landings_SHR_E)),TRUE_fmult_SHR_E)
TRUE_fmult.SHR.W<-c(rep(0,times=(nyrs-nyrs_landings_SHR_W)),TRUE_fmult_SHR_W)

pred_landings.HL.E<-c(rep(0,times=(nyrs-nyrs_landings_HL_E)),pred_landings_HL_E)
pred_landings.HL.W<-c(rep(0,times=(nyrs-nyrs_landings_HL_W)),pred_landings_HL_W)
pred_landings.LL.E<-c(rep(0,times=(nyrs-nyrs_landings_LL_E)),pred_landings_LL_E)
pred_landings.LL.W<-c(rep(0,times=(nyrs-nyrs_landings_LL_W)),pred_landings_LL_W)
pred_landings.MRIP.E<-c(rep(0,times=(nyrs-nyrs_landings_MRIP_E)),pred_landings_MRIP_E)
pred_landings.MRIP.W<-c(rep(0,times=(nyrs-nyrs_landings_MRIP_W)),pred_landings_MRIP_W)
pred_landings.HBT.E<-c(rep(0,times=(nyrs-nyrs_landings_HBT_E)),pred_landings_HBT_E)
pred_landings.HBT.W<-c(rep(0,times=(nyrs-nyrs_landings_HBT_W)),pred_landings_HBT_W)
pred_landings.SHR.E<-c(rep(0,times=(nyrs-nyrs_landings_SHR_E)),pred_landings_SHR_E)
pred_landings.SHR.W<-c(rep(0,times=(nyrs-nyrs_landings_SHR_W)),pred_landings_SHR_W)

OBS_landings.HL.E<-c(rep(0,times=(nyrs-nyrs_landings_HL_E)),OBS_landings_HL_E)
OBS_landings.HL.W<-c(rep(0,times=(nyrs-nyrs_landings_HL_W)),OBS_landings_HL_W)
OBS_landings.LL.E<-c(rep(0,times=(nyrs-nyrs_landings_LL_E)),OBS_landings_LL_E)
OBS_landings.LL.W<-c(rep(0,times=(nyrs-nyrs_landings_LL_W)),OBS_landings_LL_W)
OBS_landings.MRIP.E<-c(rep(0,times=(nyrs-nyrs_landings_MRIP_E)),OBS_landings_MRIP_E)
OBS_landings.MRIP.W<-c(rep(0,times=(nyrs-nyrs_landings_MRIP_W)),OBS_landings_MRIP_W)
OBS_landings.HBT.E<-c(rep(0,times=(nyrs-nyrs_landings_HBT_E)),OBS_landings_HBT_E)
OBS_landings.HBT.W<-c(rep(0,times=(nyrs-nyrs_landings_HBT_W)),OBS_landings_HBT_W)
OBS_landings.SHR.E<-c(rep(0,times=(nyrs-nyrs_landings_SHR_E)),OBS_landings_SHR_E)
OBS_landings.SHR.W<-c(rep(0,times=(nyrs-nyrs_landings_SHR_W)),OBS_landings_SHR_W)

layout(matrix(c(1,2), 2, 1, byrow = TRUE))
matplot(yrs,fmult.HL.E,type='n',main='Fully Selected Fishing Mortality by Fleet',
     xlab='Year',ylab='F',xlim=c(min(yrs)-1,max(yrs)+1),
ylim=c(0,(max(fmult_HL_E,fmult_HL_W,fmult_LL_E,fmult_LL_W,fmult_MRIP_E,fmult_MRIP_W,fmult_HBT_E,fmult_HBT_W,fmult_SHR_E,fmult_SHR_W)+.1*max(fmult_HL_E,fmult_HL_W,fmult_LL_E,fmult_LL_W,fmult_MRIP_E,fmult_MRIP_W,fmult_HBT_E,fmult_HBT_W,fmult_SHR_E,fmult_SHR_W)))) 
grid()
matlines(yrs,fmult.HL.E,col='cyan',lwd=2);
matlines(yrs,fmult.HL.W,col='green',lwd=2);
matlines(yrs,fmult.LL.E,col='red',lwd=2); 
matlines(yrs,fmult.LL.W,col='yellow',lwd=2);
matlines(yrs,fmult.MRIP.E,col='orange',lwd=2);
matlines(yrs,fmult.MRIP.W,col='purple',lwd=2);
matlines(yrs,fmult.HBT.E,col='orchid',lwd=2);
matlines(yrs,fmult.HBT.W,col='blue',lwd=2);
matlines(yrs,fmult.SHR.E,col='black',lwd=2);
matlines(yrs,fmult.SHR.W,col='gray',lwd=2);
legend("top",ncol=5, c("HL_E","HL_W", "LL_E","LL_W","MRIP_E","MRIP_W","HBT_E","HBT_W","SHR_E","SHR_W"
               ), col = c('cyan','green','red','yellow','orange','purple','orchid','blue','black','gray'),lwd=2,cex=.45)

matplot(ages,selectivity_HL_E[1,],type='n',main='Fishery Selectivity',
xlab='Age',ylab='Fishery Selectivity',xlim=c(ages[1],nages),ylim=c(0,1.2))
grid() 
matlines(ages,selectivity_HL_E[1,],col='cyan',lwd=2);
matlines(ages,selectivity_HL_W[1,],col='green',lwd=2);
matlines(ages,selectivity_LL_E[1,],col='red',lwd=2); 
matlines(ages,selectivity_LL_W[1,],col='yellow',lwd=2);
matlines(ages,selectivity_MRIP_E[1,],col='orange',lwd=2);
matlines(ages,selectivity_MRIP_W[1,],col='purple',lwd=2);
matlines(ages,selectivity_HBT_E[1,],col='orchid',lwd=2);
matlines(ages,selectivity_HBT_W[1,],col='blue',lwd=2);
matlines(ages,selectivity_SHR_E[1,],col='black',lwd=2);
matlines(ages,selectivity_SHR_W[1,],col='gray',lwd=2);
legend("top",ncol=5, c("HL_E","HL_W", "LL_E","LL_W","MRIP_E","MRIP_W","HBT_E","HBT_W","SHR_E","SHR_W"
               ), col = c('cyan','green','red','yellow','orange','purple','orchid','blue','black','gray'),lwd=2,cex=.45)



############################ Surveys ##########################################

index.SUMM.E<-c(rep(0,times=(nyrs-nyrs_index_SUMM_E)),Index_SUMM_E)
index.SUMM.W<-c(rep(0,times=(nyrs-nyrs_index_SUMM_W)),Index_SUMM_W)
index.FALL.E<-c(rep(0,times=(nyrs-nyrs_index_FALL_E)),Index_FALL_E)
index.FALL.W<-c(rep(0,times=(nyrs-nyrs_index_FALL_W)),Index_FALL_W)

OBS.index.SUMM.E<-c(rep(0,times=(nyrs-nyrs_index_SUMM_E)),OBS_index_SUMM_E)
OBS.index.SUMM.W<-c(rep(0,times=(nyrs-nyrs_index_SUMM_W)),OBS_index_SUMM_W)
OBS.index.FALL.E<-c(rep(0,times=(nyrs-nyrs_index_FALL_E)),OBS_index_FALL_E)
OBS.index.FALL.W<-c(rep(0,times=(nyrs-nyrs_index_FALL_W)),OBS_index_FALL_W)

layout(matrix(c(1,2), 2, 1, byrow = TRUE))
matplot(yrs,index.SUMM.E,type='n',main='Predicted and Observed Fishery Independent Surveys',
     xlab='Year',ylab='CPUE',xlim=c(min(yrs)-1,max(yrs)+1),
ylim=c(0,(max(index.SUMM.E,index.SUMM.W,index.FALL.E,index.FALL.W,OBS.index.SUMM.E,OBS.index.SUMM.W,OBS.index.FALL.E,OBS.index.FALL.W)+.1*max(index.SUMM.E,index.SUMM.W,index.FALL.E,index.FALL.W,OBS.index.SUMM.E,OBS.index.SUMM.W,OBS.index.FALL.E,OBS.index.FALL.W)))) 
grid()
matlines(yrs,index.SUMM.E,col='cyan',lwd=2);
matlines(yrs,index.SUMM.W,col='green',lwd=2);
matlines(yrs,index.FALL.E,col='red',lwd=2); 
matlines(yrs,index.FALL.W,col='yellow',lwd=2);
matpoints(yrs,OBS.index.SUMM.E,col='cyan',pch=15);
matpoints(yrs,OBS.index.SUMM.W,col='green',pch=16);
matpoints(yrs,OBS.index.FALL.E,col='red',pch=17); 
matpoints(yrs,OBS.index.FALL.W,col='yellow',pch=18);
legend("top",ncol=5, c("SUMM_E","SUMM_W", "FALL_E","FALL_W","OBS_SUMM_E","OBS_SUMM_W", "OBS_FALL_E","OBS_FALL_W"),
 col = c('cyan','green','red','yellow','cyan','green','red','yellow'),lty=c(1,1,1,1,-1,-1,-1,-1),pch=c(-1,-1,-1,-1,15,16,17,18),lwd=2,cex=.45)


matplot(ages,selectivity_SUMM_E[1,],type='n',main='Survey Selectivity',
xlab='Age',ylab='Survey Selectivity',xlim=c(ages[1],nages),ylim=c(0,1.2))
grid() 
matlines(ages,selectivity_SUMM_E[1,],col='cyan',lwd=2);
matlines(ages,selectivity_SUMM_W[1,],col='green',lwd=2);
matlines(ages,selectivity_FALL_E[1,],col='red',lwd=2); 
matlines(ages,selectivity_FALL_W[1,],col='yellow',lwd=2);

legend("top",ncol=5, c("SUMM_E","SUMM_W", "FALL_E","FALL_W" ), col = c('cyan','green','red','yellow'),lwd=2,cex=.45)


############################ F and Selectivity by fishery ##########################################

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))

matplot(yrs,pred_landings.HL.E,type='n',main='Landings Handlines',
     xlab='Year',ylab='Landings (mt)',xlim=c(min(yrs)-1,max(yrs)+1),
ylim=c(0,(max(OBS_landings.HL.E,OBS_landings.HL.W,pred_landings_HL_E,pred_landings_HL_W)+.1*max(OBS_landings.HL.E,OBS_landings.HL.W,pred_landings_HL_E,pred_landings_HL_W)))) 
grid()
matlines(yrs,pred_landings.HL.E,col='cyan',lwd=2);
matlines(yrs,pred_landings.HL.W,col='green',lwd=2);
matlines(yrs,OBS_landings.HL.E,col='cyan',lwd=2,lty=3);
matlines(yrs,OBS_landings.HL.W,col='green',lwd=2,lty=3);
legend("topright",ncol=1, c("HL_E","HL_W"
               ), col = c('cyan','green'),lwd=2,cex=.45)

matplot(yrs,fmult.HL.E,type='n',main='Fully Selected Fishing Mortality \nHandlines',
     xlab='Year',ylab='F',xlim=c(min(yrs)-1,max(yrs)+1),
ylim=c(0,(max(TRUE_fmult.HL.E,TRUE_fmult.HL.W,fmult_HL_E,fmult_HL_W)+.1*max(TRUE_fmult.HL.E,TRUE_fmult.HL.W,fmult_HL_E,fmult_HL_W)))) 
grid()
matlines(yrs,fmult.HL.E,col='cyan',lwd=2);
matlines(yrs,fmult.HL.W,col='green',lwd=2);
matlines(yrs,TRUE_fmult.HL.E,col='cyan',lwd=2,lty=3);
matlines(yrs,TRUE_fmult.HL.W,col='green',lwd=2,lty=3);
legend("topright",ncol=1, c("HL_E","HL_W"
               ), col = c('cyan','green'),lwd=2,cex=.45)

matplot(ages,selectivity_HL_E[1,],type='n',main='Beginning Year Handline Selectivity',
xlab='Age',ylab='Fishery Selectivity',xlim=c(ages[1],nages),ylim=c(0,1.1))
grid() 
matlines(ages,selectivity_HL_E[1,],col='cyan',lwd=2);
matlines(ages,selectivity_HL_W[1,],col='green',lwd=2);
matlines(ages,TRUE_selectivity_HL_E[1,],col='cyan',lwd=2,lty=3);
matlines(ages,TRUE_selectivity_HL_W[1,],col='green',lwd=2,lty=3);
legend("topright",ncol=1, c("HL_E","HL_W"
               ), col = c('cyan','green'),lwd=2,cex=.45)

matplot(ages,selectivity_HL_E[nages,],type='n',main='End Year Handline Selectivity',
xlab='Age',ylab='Fishery Selectivity',xlim=c(ages[1],nages),ylim=c(0,1.1))
grid() 
matlines(ages,selectivity_HL_E[nages,],col='cyan',lwd=2);
matlines(ages,selectivity_HL_W[nages,],col='green',lwd=2);
matlines(ages,TRUE_selectivity_HL_E[nages,],col='cyan',lwd=2,lty=3);
matlines(ages,TRUE_selectivity_HL_W[nages,],col='green',lwd=2,lty=3);
legend("topright",ncol=1, c("HL_E","HL_W"
               ), col = c('cyan','green'),lwd=2,cex=.45)



layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))

matplot(yrs,pred_landings.LL.E,type='n',main='Landings Longlines',
     xlab='Year',ylab='Landings (mt)',xlim=c(min(yrs)-1,max(yrs)+1),
ylim=c(0,(max(OBS_landings.LL.E,OBS_landings.LL.W,pred_landings_LL_E,pred_landings_LL_W)+.1*max(OBS_landings.LL.E,OBS_landings.LL.W,pred_landings_LL_E,pred_landings_LL_W)))) 
grid()
matlines(yrs,pred_landings.LL.E,col='red',lwd=2);
matlines(yrs,pred_landings.LL.W,col='yellow',lwd=2);
matlines(yrs,OBS_landings.LL.E,col='red',lwd=2,lty=3);
matlines(yrs,OBS_landings.LL.W,col='yellow',lwd=2,lty=3);
legend("topright",ncol=1, c("LL_E","LL_W"
               ), col = c('red','yellow'),lwd=2,cex=.45)

matplot(yrs,fmult.LL.E,type='n',main='Fully Selected Fishing Mortality \nLonglines',
     xlab='Year',ylab='F',xlim=c(min(yrs)-1,max(yrs)+1),
ylim=c(0,(max(TRUE_fmult.LL.E,TRUE_fmult.LL.W,fmult_LL_E,fmult_LL_W)+.1*max(TRUE_fmult.LL.E,TRUE_fmult.LL.W,fmult_LL_E,fmult_LL_W)))) 
grid()
matlines(yrs,fmult.LL.E,col='red',lwd=2);
matlines(yrs,fmult.LL.W,col='yellow',lwd=2);
matlines(yrs,TRUE_fmult.LL.E,col='red',lwd=2,lty=3);
matlines(yrs,TRUE_fmult.LL.W,col='yellow',lwd=2,lty=3);
legend("topright",ncol=1, c("LL_E","LL_W"
               ), col = c('red','yellow'),lwd=2,cex=.45)

matplot(ages,selectivity_LL_E[1,],type='n',main='Beginning Year Longline Selectivity',
xlab='Age',ylab='Fishery Selectivity',xlim=c(ages[1],nages),ylim=c(0,1.1))
grid() 
matlines(ages,selectivity_LL_E[1,],col='red',lwd=2);
matlines(ages,selectivity_LL_W[1,],col='yellow',lwd=2);
matlines(ages,TRUE_selectivity_LL_E[1,],col='red',lwd=2,lty=3);
matlines(ages,TRUE_selectivity_LL_W[1,],col='yellow',lwd=2,lty=3);
legend("topright",ncol=1, c("LL_E","LL_W"
               ), col = c('red','yellow'),lwd=2,cex=.45)

matplot(ages,selectivity_LL_E[nages,],type='n',main='End Year Longline Selectivity',
xlab='Age',ylab='Fishery Selectivity',xlim=c(ages[1],nages),ylim=c(0,1.1))
grid() 
matlines(ages,selectivity_LL_E[nages,],col='red',lwd=2);
matlines(ages,selectivity_LL_W[nages,],col='yellow',lwd=2);
matlines(ages,TRUE_selectivity_LL_E[nages,],col='red',lwd=2,lty=3);
matlines(ages,TRUE_selectivity_LL_W[nages,],col='yellow',lwd=2,lty=3);
legend("topright",ncol=1, c("LL_E","LL_W"
               ), col = c('red','yellow'),lwd=2,cex=.45)



layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))

matplot(yrs,pred_landings.MRIP.E,type='n',main='Landings MRIP',
     xlab='Year',ylab='Landings (1000s of Fish)',xlim=c(min(yrs)-1,max(yrs)+1),
ylim=c(0,(max(OBS_landings.MRIP.W,OBS_landings.MRIP.E,pred_landings_MRIP_E,pred_landings_MRIP_W)+.1*max(OBS_landings.MRIP.W,OBS_landings.MRIP.E,pred_landings_MRIP_E,pred_landings_MRIP_W)))) 
grid()
matlines(yrs,pred_landings.MRIP.E,col='orange',lwd=2);
matlines(yrs,pred_landings.MRIP.W,col='purple',lwd=2);
matlines(yrs,OBS_landings.MRIP.E,col='orange',lwd=2,lty=3);
matlines(yrs,OBS_landings.MRIP.W,col='purple',lwd=2,lty=3);
legend("topright",ncol=1, c("MRIP_E","MRIP_W"
               ), col = c('orange','purple'),lwd=2,cex=.45)

matplot(yrs,fmult.MRIP.E,type='n',main='Fully Selected Fishing Mortality \nMRIP',
     xlab='Year',ylab='F',xlim=c(min(yrs)-1,max(yrs)+1),
ylim=c(0,(max(fmult_MRIP_E,fmult_MRIP_W,TRUE_fmult.MRIP.E,TRUE_fmult.MRIP.W)+.1*max(fmult_MRIP_E,TRUE_fmult.MRIP.E,TRUE_fmult.MRIP.W,fmult_MRIP_W)))) 
grid()
matlines(yrs,fmult.MRIP.E,col='orange',lwd=2);
matlines(yrs,fmult.MRIP.W,col='purple',lwd=2);
matlines(yrs,TRUE_fmult.MRIP.E,col='orange',lwd=2,lty=3);
matlines(yrs,TRUE_fmult.MRIP.W,col='purple',lwd=2,lty=3);
legend("topright",ncol=1, c("MRIP_E","MRIP_W"
               ), col = c('orange','purple'),lwd=2,cex=.45)

matplot(ages,selectivity_MRIP_E[1,],type='n',main='Beginning Year MRIP Selectivity',
xlab='Age',ylab='Fishery Selectivity',xlim=c(ages[1],nages),ylim=c(0,1.1))
grid() 
matlines(ages,selectivity_MRIP_E[1,],col='orange',lwd=2);
matlines(ages,selectivity_MRIP_W[1,],col='purple',lwd=2);
matlines(ages,TRUE_selectivity_MRIP_E[1,],col='orange',lwd=2,lty=3);
matlines(ages,TRUE_selectivity_MRIP_W[1,],col='purple',lwd=2,lty=3);
legend("topright",ncol=1, c("MRIP_E","MRIP_W"
               ), col = c('orange','purple'),lwd=2,cex=.45)

matplot(ages,selectivity_MRIP_E[nages,],type='n',main='End Year MRIP Selectivity',
xlab='Age',ylab='Fishery Selectivity',xlim=c(ages[1],nages),ylim=c(0,1.1))
grid() 
matlines(ages,selectivity_MRIP_E[nages,],col='orange',lwd=2);
matlines(ages,selectivity_MRIP_W[nages,],col='purple',lwd=2);
matlines(ages,TRUE_selectivity_MRIP_E[nages,],col='orange',lwd=2,lty=3);
matlines(ages,TRUE_selectivity_MRIP_W[nages,],col='purple',lwd=2,lty=3);
legend("topright",ncol=1, c("MRIP_E","MRIP_W"
               ), col = c('orange','purple'),lwd=2,cex=.45)


layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))

matplot(yrs,pred_landings.HBT.E,type='n',main='Landings Headboat',
     xlab='Year',ylab='Landings (1000s of Fish)',xlim=c(min(yrs)-1,max(yrs)+1),
ylim=c(0,(max(OBS_landings.HBT.E,OBS_landings.HBT.W,pred_landings_HBT_E,pred_landings_HBT_W)+.1*max(OBS_landings.HBT.E,OBS_landings.HBT.W,pred_landings_HBT_E,pred_landings_HBT_W)))) 
grid()
matlines(yrs,pred_landings.HBT.E,col='orchid',lwd=2);
matlines(yrs,pred_landings.HBT.W,col='blue',lwd=2);
matlines(yrs,OBS_landings.HBT.E,col='orchid',lwd=2,lty=3);
matlines(yrs,OBS_landings.HBT.W,col='blue',lwd=2,lty=3);
legend("topright",ncol=1, c("HBT_E","HBT_W"
               ), col = c('orchid','blue'),lwd=2,cex=.45)

matplot(yrs,fmult.HBT.E,type='n',main='Fully Selected Fishing Mortality \nHeadboat',
     xlab='Year',ylab='F',xlim=c(min(yrs)-1,max(yrs)+1),
ylim=c(0,(max(TRUE_fmult.HBT.W,TRUE_fmult.HBT.E,fmult_HBT_E,fmult_HBT_W)+.1*max(TRUE_fmult.HBT.W,TRUE_fmult.HBT.E,fmult_HBT_E,fmult_HBT_W)))) 
grid()
matlines(yrs,fmult.HBT.E,col='orchid',lwd=2);
matlines(yrs,fmult.HBT.W,col='blue',lwd=2);
matlines(yrs,TRUE_fmult.HBT.E,col='orchid',lwd=2,lty=3);
matlines(yrs,TRUE_fmult.HBT.W,col='blue',lwd=2,lty=3);
legend("topright",ncol=1, c("HBT_E","HBT_W"
               ), col = c('orchid','blue'),lwd=2,cex=.45)

matplot(ages,selectivity_HBT_E[1,],type='n',main='Beginning Year Headboat Selectivity',
xlab='Age',ylab='Fishery Selectivity',xlim=c(ages[1],nages),ylim=c(0,1.1))
grid() 
matlines(ages,selectivity_HBT_E[1,],col='orchid',lwd=2);
matlines(ages,selectivity_HBT_W[1,],col='blue',lwd=2);
matlines(ages,TRUE_selectivity_HBT_E[1,],col='orchid',lwd=2,lty=3);
matlines(ages,TRUE_selectivity_HBT_W[1,],col='blue',lwd=2,lty=3);
legend("topright",ncol=1, c("HBT_E","HBT_W"
               ), col = c('orchid','blue'),lwd=2,cex=.45)

matplot(ages,selectivity_HBT_E[nages,],type='n',main='End Year Headboat Selectivity',
xlab='Age',ylab='Fishery Selectivity',xlim=c(ages[1],nages),ylim=c(0,1.1))
grid() 
matlines(ages,selectivity_HBT_E[nages,],col='orchid',lwd=2);
matlines(ages,selectivity_HBT_W[nages,],col='blue',lwd=2);
matlines(ages,TRUE_selectivity_HBT_E[nages,],col='orchid',lwd=2,lty=3);
matlines(ages,TRUE_selectivity_HBT_W[nages,],col='blue',lwd=2,lty=3);
legend("topright",ncol=1, c("HBT_E","HBT_W"
               ), col = c('orchid','blue'),lwd=2,cex=.45)


layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))

matplot(yrs,pred_landings.SHR.E,type='n',main='Bycatch Shrimp',
     xlab='Year',ylab='Bycatch (1000s of Fish)',xlim=c(min(yrs)-1,max(yrs)+1),
ylim=c(0,(max(OBS_landings.SHR.E,OBS_landings.SHR.W,pred_landings_SHR_E,pred_landings_SHR_W)+.1*max(OBS_landings.SHR.E,OBS_landings.SHR.W,pred_landings_SHR_E,pred_landings_SHR_W)))) 
grid()
matlines(yrs,pred_landings.SHR.E,col='black',lwd=2);
matlines(yrs,pred_landings.SHR.W,col='gray',lwd=2);
matlines(yrs,OBS_landings.SHR.E,col='black',lwd=2,lty=3);
matlines(yrs,OBS_landings.SHR.W,col='gray',lwd=2,lty=3);
legend("topright",ncol=1, c("SHR_E","SHR_W","Mean_E","Mean_W","OBS_Mean_E","OBS_Mean_W"
               ), col = c('black','gray','black','gray','blue','red'),lwd=2,cex=.45,lty=c(1,1,2,2,3,3))

matplot(yrs,fmult.SHR.E,type='n',main='Fully Selected Fishing Mortality \nShrimp',
     xlab='Year',ylab='F',xlim=c(min(yrs)-1,max(yrs)+1),
ylim=c(0,(max(TRUE_fmult.SHR.E,TRUE_fmult.SHR.W,fmult_SHR_E,fmult_SHR_W)+.1*max(TRUE_fmult.SHR.E,TRUE_fmult.SHR.W,fmult_SHR_E,fmult_SHR_W)))) 
grid()
matlines(yrs,fmult.SHR.E,col='black',lwd=2);
matlines(yrs,fmult.SHR.W,col='gray',lwd=2);
matlines(yrs,TRUE_fmult.SHR.E,col='black',lwd=2,lty=3);
matlines(yrs,TRUE_fmult.SHR.W,col='gray',lwd=2,lty=3);
legend("topright",ncol=1, c("SHR_E","SHR_W"
               ), col = c('black','gray'),lwd=2,cex=.45)

matplot(ages,selectivity_SHR_E[1,],type='n',main='Beginning Year Shrimp Selectivity',
xlab='Age',ylab='Fishery Selectivity',xlim=c(ages[1],nages),ylim=c(0,1.1))
grid() 
matlines(ages,selectivity_SHR_E[1,],col='black',lwd=2);
matlines(ages,selectivity_SHR_W[1,],col='gray',lwd=2);
matlines(ages,TRUE_selectivity_SHR_E[1,],col='black',lwd=2,lty=3);
matlines(ages,TRUE_selectivity_SHR_W[1,],col='gray',lwd=2,lty=3);
legend("topright",ncol=1, c("SHR_E","SHR_W"
               ), col = c('black','gray'),lwd=2,cex=.45)

matplot(ages,selectivity_SHR_E[nages,],type='n',main='End Year Shrimp Selectivity',
xlab='Age',ylab='Fishery Selectivity',xlim=c(ages[1],nages),ylim=c(0,1.1))
grid() 
matlines(ages,selectivity_SHR_E[nages,],col='black',lwd=2);
matlines(ages,selectivity_SHR_W[nages,],col='gray',lwd=2);
matlines(ages,TRUE_selectivity_SHR_E[nages,],col='black',lwd=2,lty=3);
matlines(ages,TRUE_selectivity_SHR_W[nages,],col='gray',lwd=2,lty=3);
legend("topright",ncol=1, c("SHR_E","SHR_W"
               ), col = c('black','gray'),lwd=2,cex=.45)


#################### Recruitment ##############################################################
#ssb.BH.E<-seq(1,(max(ssbE)+.5*max(ssbE)),by=1000000000)
#ssb.BH.W<-seq(1,(max(ssbW)+.5*max(ssbW)),by=100000000000)

#BH_REC.E=NULL
#BH_REC.W=NULL
#for(i in 1:length(ssb.BH.E))
#{
#BH_REC.E[i]<-(4*h_E*R_aveE*ssb.BH.E[i])/(SSB_zeroE*(1-h_E)+ssb.BH.E[i]*(5*h_E-1))
#}
#for(i in 1:length(ssb.BH.W))
#{
#BH_REC.W[i]<-(4*h_W*R_aveW*ssb.BH.W[i])/(SSB_zeroW*(1-h_W)+ssb.BH.W[i]*(5*h_W-1))
#}

recruitsE1<-recruitsE_BM[2:nyrs]
ssbE1<-ssbE[1:(nyrs-1)]

layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE))
par(oma=c(0,0,0,1.4))
plot((yrs[1:(nyrs-1)]),ssbE1,type="h",lwd=22,ylim=c(0,max(ssbE1)+.2*max(ssbE1)),col="red4",
xlab='Year',ylab='SSB (# Eggs)',main='SSB and Recruitment in the East\nStock Specific (Before Movement)',
sub='Recruitment in Year+1 Overlayed with Causal SSB in Current Year')
grid()
par(new=T)
plot(recruitsE1,type="o",pch=16,cex=1.5,lwd=4,ylim=c(0,max(recruitsE1)+.2*max(recruitsE1)),
axes=FALSE,xlab="",ylab="",col="blue4")
axis(4,ylim=c(0,max(ssbE1)+.2*max(ssbE1)))
mtext('Recruitment (1000s of Fish)',side=4,line=2)
legend("top",ncol=2, c(" SSB", " Recruits"), lty=c(1,1),lwd=2,
pch=c(-18,18),col = c("red4","blue4"),cex=.7, merge=FALSE)

R_aveE1<-c(rep(R_aveE,(nyrs)))
TRUE_R_aveE1<-c(rep(TRUE_R_aveE,(nyrs)))
matplot((yrs[1:(nyrs)]),recruitsE_BM,type='n',main='Estimated Recruitment East',xlab='Year',
ylab='Recruitment (1000s of fish)',xlim=c(min(yrs)-1,max(yrs)+1),
        ylim=c(0,max(recruitsE_BM,TRUE_recruitsE_BM)+.2*max(recruitsE_BM,TRUE_recruitsE_BM)))
grid()
matlines((yrs[1:(nyrs)]),recruitsE_BM,col="blue",type='b',pch=18,lwd=2);
matlines((yrs[1:(nyrs)]),R_aveE1,col="blue",lwd=2,lty=3);
matlines((yrs[1:(nyrs)]),TRUE_recruitsE_BM,col="red",type='b',pch=18,lwd=2);
matlines((yrs[1:(nyrs)]),TRUE_R_aveE1,col="red",lwd=2,lty=3);
legend("top",ncol=2, c("Recruits", "R_Ave","True Recruits", "True R_Ave"), lty=c(1,3,1,3),lwd=2,
pch=c(18,-18,18,-18),col = c("blue","blue","red","red"),cex=.6, merge=FALSE)

#matplot(ssb.BH.E,BH_REC.E,type='n',main='Recruitment vs. SSB East'
#,xlab='SSB (# Eggs)',
#ylab='Recruitment (1000s of fish)',ylim=c(0,max(recruitsE_BM,BH_REC.E)+.25*max(recruitsE_BM,BH_REC.E))) 
#grid()
#matpoints(ssbE1,recruitsE1,col="red",pch=18);
#matlines(ssb.BH.E,BH_REC.E,col="blue",lwd=2);
#legend("top",ncol=2, c("Estimate", "BH"), lty=c(-1,1),lwd=2,
#pch=c(18,-18),col = c("red","blue"),cex=.6, merge=FALSE)





recruitsW1<-recruitsW_BM[2:nyrs]
ssbW1<-ssbW[1:(nyrs-1)]

layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE))
par(oma=c(0,0,0,1.4))
plot((yrs[1:(nyrs-1)]),ssbW1,type="h",lwd=22,ylim=c(0,max(ssbW1)+.2*max(ssbW1)),col="red4",
xlab='Year',ylab='SSB (# Eggs)',main='SSB and Recruitment in the West\nStock Specific (Before Movement)',
sub='Recruitment in Year+1 Overlayed with Causal SSB in Current Year')
grid()
par(new=T)
plot(recruitsW1,type="o",pch=16,cex=1.5,lwd=4,ylim=c(0,max(recruitsW1)+.2*max(recruitsW1)),
axes=FALSE,xlab="",ylab="",col="blue4")
axis(4,ylim=c(0,max(ssbW1)+.2*max(ssbW1)))
mtext('Recruitment (1000s of Fish)',side=4,line=2)
legend("top",ncol=2, c(" SSB", " Recruits"), lty=c(1,1),lwd=2,
pch=c(-18,18),col = c("red4","blue4"),cex=.7, merge=FALSE)

R_aveW1<-c(rep(R_aveW,(nyrs)))
TRUE_R_aveW1<-c(rep(TRUE_R_aveW,(nyrs)))
matplot((yrs[1:(nyrs)]),recruitsW_BM,type='n',main='Estimated Recruitment West',xlab='Year',
ylab='Recruitment (1000s of fish)',xlim=c(min(yrs)-1,max(yrs)+1),
        ylim=c(0,max(recruitsW_BM,TRUE_recruitsW_BM)+.2*max(recruitsW_BM,TRUE_recruitsW_BM)))
grid()
matlines((yrs[1:(nyrs)]),recruitsW_BM,col="blue",type='b',pch=18,lwd=2);
matlines((yrs[1:(nyrs)]),R_aveW1,col="blue",lwd=2,lty=3);
matlines((yrs[1:(nyrs)]),TRUE_recruitsW_BM,col="red",type='b',pch=18,lwd=2);
matlines((yrs[1:(nyrs)]),TRUE_R_aveW1,col="red",lwd=2,lty=3);
legend("top",ncol=2, c("Recruits", "R_Ave","True Recruits", "True R_Ave"), lty=c(1,3,1,3),lwd=2,
pch=c(18,-18,18,-18),col = c("blue","blue","red","red"),cex=.6, merge=FALSE)



layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE))


matplot((yrs[1:(nyrs)]),recruitsE_BM,type='n',main='Estimated Recruitment East',xlab='Year',
ylab='Recruitment (1000s of fish)',xlim=c(min(yrs)-1,max(yrs)+1),
        ylim=c(0,max(recruitsE_BM,TRUE_recruitsE_BM)+.2*max(recruitsE_BM,TRUE_recruitsE_BM)))
grid()
matlines((yrs[1:(nyrs)]),recruitsE_BM,col="blue",type='b',pch=18,lwd=2);
matlines((yrs[1:(nyrs)]),R_aveE1,col="blue",lwd=2,lty=3);
matlines((yrs[1:(nyrs)]),TRUE_recruitsE_BM,col="red",type='b',pch=18,lwd=2);
matlines((yrs[1:(nyrs)]),TRUE_R_aveE1,col="red",lwd=2,lty=3);
legend("top",ncol=2, c("Recruits", "R_Ave","True Recruits", "True R_Ave"), lty=c(1,3,1,3),lwd=2,
pch=c(18,-18,18,-18),col = c("blue","blue","red","red"),cex=.6, merge=FALSE)

matplot((yrs[1:(nyrs)]),recruitsE_AM,type='n',main='Estimated Recruitment After Larval Movement East',xlab='Year',
ylab='Recruitment (1000s of fish)',xlim=c(min(yrs)-1,max(yrs)+1),
        ylim=c(0,max(recruitsE_AM,TRUE_recruitsE_AM)+.2*max(recruitsE_AM,TRUE_recruitsE_AM)))
grid()
matlines((yrs[1:(nyrs)]),recruitsE_AM,col="blue",type='b',pch=18,lwd=2);
matlines((yrs[1:(nyrs)]),TRUE_recruitsE_AM,col="red",type='b',pch=18,lwd=2);
legend("top",ncol=2, c("Recruits", "True Recruits" ), lty=c(1,3),lwd=2,
pch=c(18,18),col = c("blue","red"),cex=.6, merge=FALSE)


layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE))

matplot((yrs[109:(nyrs)]),recruitsE_BM[109:(nyrs)],type='n',main='Estimated Recruitment East',xlab='Year',
ylab='Recruitment (1000s of fish)',xlim=c(yrs[109]-1,max(yrs)+1),
        ylim=c(0,max(recruitsE_BM[109:(nyrs)],TRUE_recruitsE_BM[109:(nyrs)])+.2*max(recruitsE_BM[109:(nyrs)],TRUE_recruitsE_BM[109:(nyrs)])))
grid()
matlines((yrs[109:(nyrs)]),recruitsE_BM[109:(nyrs)],col="blue",type='b',pch=18,lwd=2);
matlines((yrs[109:(nyrs)]),R_aveE1[109:(nyrs)],col="blue",lwd=2,lty=3);
matlines((yrs[109:(nyrs)]),TRUE_recruitsE_BM[109:(nyrs)],col="red",type='b',pch=18,lwd=2);
matlines((yrs[109:(nyrs)]),TRUE_R_aveE1[109:(nyrs)],col="red",lwd=2,lty=3);
legend("top",ncol=2, c("Recruits", "R_Ave","True Recruits", "True R_Ave"), lty=c(1,3,1,3),lwd=2,
pch=c(18,-18,18,-18),col = c("blue","blue","red","red"),cex=.6, merge=FALSE)

matplot((yrs[109:(nyrs)]),recruitsE_AM[109:(nyrs)],type='n',main='Estimated Recruitment After Larval Movement East',xlab='Year',
ylab='Recruitment (1000s of fish)',xlim=c(yrs[109]-1,max(yrs)+1),
        ylim=c(0,max(recruitsE_AM[109:(nyrs)],TRUE_recruitsE_AM[109:(nyrs)])+.2*max(recruitsE_AM[109:(nyrs)],TRUE_recruitsE_AM[109:(nyrs)])))
grid()
matlines((yrs[109:(nyrs)]),recruitsE_AM[109:(nyrs)],col="blue",type='b',pch=18,lwd=2);
matlines((yrs[109:(nyrs)]),TRUE_recruitsE_AM[109:(nyrs)],col="red",type='b',pch=18,lwd=2);
legend("top",ncol=2, c("Recruits", "True Recruits" ), lty=c(1,3),lwd=2,
pch=c(18,18),col = c("blue","red"),cex=.6, merge=FALSE)



layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE))


matplot((yrs[1:(nyrs)]),recruitsW_BM,type='n',main='Estimated Recruitment West',xlab='Year',
ylab='Recruitment (1000s of fish)',xlim=c(min(yrs)-1,max(yrs)+1),
        ylim=c(0,max(recruitsW_BM,TRUE_recruitsW_BM)+.2*max(recruitsW_BM,TRUE_recruitsW_BM)))
grid()
matlines((yrs[1:(nyrs)]),recruitsW_BM,col="blue",type='b',pch=18,lwd=2);
matlines((yrs[1:(nyrs)]),R_aveW1,col="blue",lwd=2,lty=3);
matlines((yrs[1:(nyrs)]),TRUE_recruitsW_BM,col="red",type='b',pch=18,lwd=2);
matlines((yrs[1:(nyrs)]),TRUE_R_aveW1,col="red",lwd=2,lty=3);
legend("top",ncol=2, c("Recruits", "R_Ave","True Recruits", "True R_Ave"), lty=c(1,3,1,3),lwd=2,
pch=c(18,-18,18,-18),col = c("blue","blue","red","red"),cex=.6, merge=FALSE)

matplot((yrs[1:(nyrs)]),recruitsW_AM,type='n',main='Estimated Recruitment After Larval Movement West',xlab='Year',
ylab='Recruitment (1000s of fish)',xlim=c(min(yrs)-1,max(yrs)+1),
        ylim=c(0,max(recruitsW_AM,TRUE_recruitsW_AM)+.2*max(recruitsW_AM,TRUE_recruitsW_AM)))
grid()
matlines((yrs[1:(nyrs)]),recruitsW_AM,col="blue",type='b',pch=18,lwd=2);
matlines((yrs[1:(nyrs)]),TRUE_recruitsW_AM,col="red",type='b',pch=18,lwd=2);
legend("top",ncol=2, c("Recruits", "True Recruits" ), lty=c(1,3),lwd=2,
pch=c(18,18),col = c("blue","red"),cex=.6, merge=FALSE)


layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE))

matplot((yrs[109:(nyrs)]),recruitsW_BM[109:(nyrs)],type='n',main='Estimated Recruitment West',xlab='Year',
ylab='Recruitment (1000s of fish)',xlim=c(yrs[109]-1,max(yrs)+1),
        ylim=c(0,max(recruitsW_BM[109:(nyrs)],TRUE_recruitsW_BM[109:(nyrs)])+.2*max(recruitsW_BM[109:(nyrs)],TRUE_recruitsW_BM[109:(nyrs)])))
grid()
matlines((yrs[109:(nyrs)]),recruitsW_BM[109:(nyrs)],col="blue",type='b',pch=18,lwd=2);
matlines((yrs[109:(nyrs)]),R_aveW1[109:(nyrs)],col="blue",lwd=2,lty=3);
matlines((yrs[109:(nyrs)]),TRUE_recruitsW_BM[109:(nyrs)],col="red",type='b',pch=18,lwd=2);
matlines((yrs[109:(nyrs)]),TRUE_R_aveW1[109:(nyrs)],col="red",lwd=2,lty=3);
legend("top",ncol=2, c("Recruits", "R_Ave","True Recruits", "True R_Ave"), lty=c(1,3,1,3),lwd=2,
pch=c(18,-18,18,-18),col = c("blue","blue","red","red"),cex=.6, merge=FALSE)

matplot((yrs[109:(nyrs)]),recruitsW_AM[109:(nyrs)],type='n',main='Estimated Recruitment After Larval Movement West',xlab='Year',
ylab='Recruitment (1000s of fish)',xlim=c(yrs[109]-1,max(yrs)+1),
        ylim=c(0,max(recruitsW_AM[109:(nyrs)],TRUE_recruitsW_AM[109:(nyrs)])+.2*max(recruitsW_AM[109:(nyrs)],TRUE_recruitsW_AM[109:(nyrs)])))
grid()
matlines((yrs[109:(nyrs)]),recruitsW_AM[109:(nyrs)],col="blue",type='b',pch=18,lwd=2);
matlines((yrs[109:(nyrs)]),TRUE_recruitsW_AM[109:(nyrs)],col="red",type='b',pch=18,lwd=2);
legend("top",ncol=2, c("Recruits", "True Recruits" ), lty=c(1,3),lwd=2,
pch=c(18,18),col = c("blue","red"),cex=.6, merge=FALSE)




layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))

matplot(yrs,TE_W[,1],type='n',main='Estimated Larval Movement',xlab='Year',
ylab='Proportion of Population Moving',
        ylim=c(0,1))
grid()
matlines(yrs,TE_W[,1],col="blue",lwd=2);
matlines(yrs,TW_E[,1],col="red",lwd=2);
legend("top",ncol=2, c("East to West", "West to East"), 
col = c("blue","red"),cex=.6, merge=FALSE,lty=1)

M_larval_E1<-rep(M_larval_E,times=nyrs)
M_larval_W1<-rep(M_larval_W,times=nyrs)
matplot(yrs,TE_W[,1],type='n',main='Estimated Larval Mortality',xlab='Year',
ylab='Yearly Mortality',
        ylim=c(0,1))
grid()
matlines(yrs,M_larval_E,col="blue",lwd=2);
matlines(yrs,M_larval_W,col="red",lwd=2);
legend("top",ncol=2, c("East", "West"), 
col = c("blue","red"),cex=.6, lty=1,merge=FALSE)



matplot(yrs,recruitsE_BM,type='n',main='Estimated Recruitment Pre and Post Movement',xlab='Year',
ylab='Recruitment (1000s of fish)',xlim=c(min(yrs)-1,max(yrs)+1),
        ylim=c(0,max(recruitsE_BM,recruitsW_BM,recruitsE_AM,recruitsW_AM)+.2*max(recruitsE_BM,recruitsW_BM,recruitsE_AM,recruitsW_AM)))
grid()
matlines(yrs,recruitsW_BM,col='salmon',,lty=2,lwd=2);
matlines(yrs,recruitsW_AM,col='red',lwd=2);
matlines(yrs,recruitsE_BM,col='cyan',lty=2,lwd=2);
matlines(yrs,recruitsE_AM,col='blue',lwd=2);
legend("top", ncol=2,c("East Before Mov.","West Before Mov.","East After Mov.","West After Mov." ), col = c('cyan','salmon','blue','red'),lty=c(2,2,1,1),lwd=2, cex=.8)




layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE))



matplot((yrs),TE_W[,1],type='n',main='Estimated Larval Movement ',xlab='Age',
ylab='Annual Movement Proportion',xlim=c((yrs[1])-1,(yrs[142])+1),
        ylim=c(0,max(TRUE_TE_W[,1],TRUE_TW_E[,1],TE_W,TW_E)+.2*max(TRUE_TW_E[,1],TRUE_TE_W[,1],TE_W,TW_E)))

grid()
matlines((yrs),TE_W[,1],col='blue',lwd=2); 
matlines((yrs),TW_E[,1],col='red',lwd=2); 
matlines((yrs),TRUE_TE_W[,1],col='blue',lwd=2,lty=3); 
matlines((yrs),TRUE_TW_E[,1],col='red',lwd=2,lty=3); 
legend("top", ncol=2,c("East to West","West to East" ), col = c('blue','red'),lwd=2, cex=.8)

abund.movE.W<-recruitsE_BM*TE_W[,1]
abund.movW.E<-recruitsW_BM*TW_E[,1]
TRUE_abund.movE.W<-TRUE_recruitsE_BM*TRUE_TE_W[,1]
TRUE_abund.movW.E<-TRUE_recruitsW_BM*TRUE_TW_E[,1]

matplot((yrs),abund.movE.W,type='n',main='Number of Recruits Moving',xlab='Year',
ylab='1000s of fish',xlim=c((yrs[1])-1,(yrs[142])+1),
        ylim=c(0,max(TRUE_abund.movE.W,TRUE_abund.movW.E,abund.movE.W,abund.movW.E)+.2*max(TRUE_abund.movE.W,TRUE_abund.movW.E,abund.movE.W,abund.movW.E)))
grid()
matlines((yrs),abund.movE.W,col="blue",lwd=2);
matlines((yrs),abund.movW.E,col="red",lwd=2);
matlines((yrs),TRUE_abund.movE.W,col="blue",lwd=2,lty=3);
matlines((yrs),TRUE_abund.movW.E,col="red",lwd=2,lty=3);
legend("top",ncol=2, c(" West_East", "East_West"), lty=c(1,1),lwd=2,
col = c("red","blue"),cex=.6, merge=FALSE)



################################### Abundance at age ############################################


layout(matrix(c(1), 1, byrow = TRUE))
years <- yrs
nyears <- length(years)
nindices <- nages
abundanceE <- abundance_at_ageE_BM 
resid.pos.color <- "blue"
resid.neg.color <- "red"
resid.col=matrix(NA, nrow=nyears, ncol=nindices)   # set color for residual bubbles
resid.col <- ifelse(abundanceE > 0.0,resid.pos.color, resid.neg.color) 
# use a scaler to increase (>1) or decrease (0<val<1) the size of the bubbles to make them look good
scale.bubble <- 10/(max(abundanceE)) 
abundanceE1 <- abundanceE*scale.bubble
# first set up a blank space for the plot
plot(seq(0,nindices-1), rep(years[1],nindices), ylim=c(years[nyears],years[1]), xlab = "Age",
     ylab = "Year",  type = "n", axes=F,main='Abundance at Age East Before Movement',xlim=c(-1,21),
sub=bquote("(Max Abundance"==.(max(abundance_at_ageE_BM))~")"))


axis(1, at= seq(0,nindices-1), lab=seq(0,nindices-1))
axis(2, at = seq(years[nyears],years[1], by=-5), lab = seq(years[nyears],years[1], by=-5),las=1)
box()
abline(h=seq( years[1], years[nyears]),col = "lightgray", lty = 1)
abline(a=1855,b=1,col = "lightgray", lty = 1)
abline(a=1860,b=1,col = "lightgray", lty = 1)
abline(a=1865,b=1,col = "lightgray", lty = 1)
abline(a=1870,b=1,col = "lightgray", lty = 1)
abline(a=1875,b=1,col = "lightgray", lty = 1)
abline(a=1880,b=1,col = "lightgray", lty = 1)
abline(a=1885,b=1,col = "lightgray", lty = 1)
abline(a=1890,b=1,col = "lightgray", lty = 1)
abline(a=1895,b=1,col = "lightgray", lty = 1)
abline(a=1900,b=1,col = "lightgray", lty = 1)
abline(a=1905,b=1,col = "lightgray", lty = 1)
abline(a=1910,b=1,col = "lightgray", lty = 1)
abline(a=1915,b=1,col = "lightgray", lty = 1)
abline(a=1920,b=1,col = "lightgray", lty = 1)
abline(a=1925,b=1,col = "lightgray", lty = 1)
abline(a=1930,b=1,col = "lightgray", lty = 1)
abline(a=1935,b=1,col = "lightgray", lty = 1)
abline(a=1940,b=1,col = "lightgray", lty = 1)
abline(a=1945,b=1,col = "lightgray", lty = 1)
abline(a=1950,b=1,col = "lightgray", lty = 1)
abline(a=1955,b=1,col = "lightgray", lty = 1)
abline(a=1960,b=1,col = "lightgray", lty = 1)
abline(a=1965,b=1,col = "lightgray", lty = 1)
abline(a=1970,b=1,col = "lightgray", lty = 1)
abline(a=1975,b=1,col = "lightgray", lty = 1)
abline(a=1980,b=1,col = "lightgray", lty = 1)
abline(a=1985,b=1,col = "lightgray", lty = 1)
abline(a=1990,b=1,col = "lightgray", lty = 1)
abline(a=1995,b=1,col = "lightgray", lty = 1)
abline(a=2000,b=1,col = "lightgray", lty = 1)
abline(a=2005,b=1,col = "lightgray", lty = 1)
abline(a=2010,b=1,col = "lightgray", lty = 1)
# fill in the bubbles using points command, cex determines size, bg determines color
for (i in 1:nyears){
  points(seq(0, nindices-1), rep((years[i]), nindices),  
  cex =abs(abundanceE1[i,]) ,  pch = 21, bg = resid.col[i,], col='black')
}


layout(matrix(c(1), 1, byrow = TRUE))
years <- yrs
nyears <- length(years)
nindices <- nages
abundanceE <- abundance_at_ageE_AM 
resid.pos.color <- "blue"
resid.neg.color <- "red"
resid.col=matrix(NA, nrow=nyears, ncol=nindices)   # set color for residual bubbles
resid.col <- ifelse(abundanceE > 0.0,resid.pos.color, resid.neg.color) 
# use a scaler to increase (>1) or decrease (0<val<1) the size of the bubbles to make them look good
scale.bubble <- 10/(max(abundanceE)) 
abundanceE1 <- abundanceE*scale.bubble
# first set up a blank space for the plot
plot(seq(0,nindices-1), rep(years[1],nindices), ylim=c(years[nyears],years[1]), xlab = "Age",
     ylab = "Year",  type = "n", axes=F,main='Abundance at Age East After Movement',xlim=c(-1,21),
sub=bquote("(Max Abundance"==.(max(abundance_at_ageE_AM))~")"))


axis(1, at= seq(0,nindices-1), lab=seq(0,nindices-1))
axis(2, at = seq(years[nyears],years[1], by=-5), lab = seq(years[nyears],years[1]-2, by=-5),las=1)
box()
abline(h=seq( years[1], years[nyears]),col = "lightgray", lty = 1)
abline(a=1855,b=1,col = "lightgray", lty = 1)
abline(a=1860,b=1,col = "lightgray", lty = 1)
abline(a=1865,b=1,col = "lightgray", lty = 1)
abline(a=1870,b=1,col = "lightgray", lty = 1)
abline(a=1875,b=1,col = "lightgray", lty = 1)
abline(a=1880,b=1,col = "lightgray", lty = 1)
abline(a=1885,b=1,col = "lightgray", lty = 1)
abline(a=1890,b=1,col = "lightgray", lty = 1)
abline(a=1895,b=1,col = "lightgray", lty = 1)
abline(a=1900,b=1,col = "lightgray", lty = 1)
abline(a=1905,b=1,col = "lightgray", lty = 1)
abline(a=1910,b=1,col = "lightgray", lty = 1)
abline(a=1915,b=1,col = "lightgray", lty = 1)
abline(a=1920,b=1,col = "lightgray", lty = 1)
abline(a=1925,b=1,col = "lightgray", lty = 1)
abline(a=1930,b=1,col = "lightgray", lty = 1)
abline(a=1935,b=1,col = "lightgray", lty = 1)
abline(a=1940,b=1,col = "lightgray", lty = 1)
abline(a=1945,b=1,col = "lightgray", lty = 1)
abline(a=1950,b=1,col = "lightgray", lty = 1)
abline(a=1955,b=1,col = "lightgray", lty = 1)
abline(a=1960,b=1,col = "lightgray", lty = 1)
abline(a=1965,b=1,col = "lightgray", lty = 1)
abline(a=1970,b=1,col = "lightgray", lty = 1)
abline(a=1975,b=1,col = "lightgray", lty = 1)
abline(a=1980,b=1,col = "lightgray", lty = 1)
abline(a=1985,b=1,col = "lightgray", lty = 1)
abline(a=1990,b=1,col = "lightgray", lty = 1)
abline(a=1995,b=1,col = "lightgray", lty = 1)
abline(a=2000,b=1,col = "lightgray", lty = 1)
abline(a=2005,b=1,col = "lightgray", lty = 1)
abline(a=2010,b=1,col = "lightgray", lty = 1)
# fill in the bubbles using points command, cex determines size, bg determines color
for (i in 1:nyears){
  points(seq(0, nindices-1), rep((years[i]), nindices),  
  cex =abs(abundanceE1[i,]) ,  pch = 21, bg = resid.col[i,], col='black')
}


layout(matrix(c(1), 1, byrow = TRUE))
years <- yrs
nyears <- length(years)
nindices <- nages
abundanceW <- abundance_at_ageW_BM 
resid.pos.color <- "red"
resid.neg.color <- "blue"
resid.col=matrix(NA, nrow=nyears, ncol=nindices)   # set color for residual bubbles
resid.col <- ifelse(abundanceW > 0.0,resid.pos.color, resid.neg.color) 
# use a scaler to increase (>1) or decrease (0<val<1) the size of the bubbles to make them look good
scale.bubble <- 10/(max(abundanceW)) 
abundanceW1 <- abundanceW*scale.bubble
# first set up a blank space for the plot
plot(seq(0,nindices-1), rep(years[1],nindices), ylim=c(years[nyears],years[1]-2), xlab = "Age",
     ylab = "Year",  type = "n", axes=F,main='Abundance at Age West Before Movement',xlim=c(-1,21),
sub=bquote("(Max Abundance"==.(max(abundance_at_ageW_BM))~")"))


axis(1, at= seq(0,nindices-1), lab=seq(0,nindices-1))
axis(2, at = seq(years[nyears],years[1], by=-5), lab = seq(years[nyears],years[1], by=-5),las=1)
box()
abline(h=seq( years[1], years[nyears]),col = "lightgray", lty = 1)
abline(a=1855,b=1,col = "lightgray", lty = 1)
abline(a=1860,b=1,col = "lightgray", lty = 1)
abline(a=1865,b=1,col = "lightgray", lty = 1)
abline(a=1870,b=1,col = "lightgray", lty = 1)
abline(a=1875,b=1,col = "lightgray", lty = 1)
abline(a=1880,b=1,col = "lightgray", lty = 1)
abline(a=1885,b=1,col = "lightgray", lty = 1)
abline(a=1890,b=1,col = "lightgray", lty = 1)
abline(a=1895,b=1,col = "lightgray", lty = 1)
abline(a=1900,b=1,col = "lightgray", lty = 1)
abline(a=1905,b=1,col = "lightgray", lty = 1)
abline(a=1910,b=1,col = "lightgray", lty = 1)
abline(a=1915,b=1,col = "lightgray", lty = 1)
abline(a=1920,b=1,col = "lightgray", lty = 1)
abline(a=1925,b=1,col = "lightgray", lty = 1)
abline(a=1930,b=1,col = "lightgray", lty = 1)
abline(a=1935,b=1,col = "lightgray", lty = 1)
abline(a=1940,b=1,col = "lightgray", lty = 1)
abline(a=1945,b=1,col = "lightgray", lty = 1)
abline(a=1950,b=1,col = "lightgray", lty = 1)
abline(a=1955,b=1,col = "lightgray", lty = 1)
abline(a=1960,b=1,col = "lightgray", lty = 1)
abline(a=1965,b=1,col = "lightgray", lty = 1)
abline(a=1970,b=1,col = "lightgray", lty = 1)
abline(a=1975,b=1,col = "lightgray", lty = 1)
abline(a=1980,b=1,col = "lightgray", lty = 1)
abline(a=1985,b=1,col = "lightgray", lty = 1)
abline(a=1990,b=1,col = "lightgray", lty = 1)
abline(a=1995,b=1,col = "lightgray", lty = 1)
abline(a=2000,b=1,col = "lightgray", lty = 1)
abline(a=2005,b=1,col = "lightgray", lty = 1)
abline(a=2010,b=1,col = "lightgray", lty = 1)
# fill in the bubbles using points command, cex determines size, bg determines color
for (i in 1:nyears){
  points(seq(0, nindices-1), rep((years[i]), nindices),  
  cex =abs(abundanceW1[i,]) ,  pch = 21, bg = resid.col[i,], col='black')
}

layout(matrix(c(1), 1, byrow = TRUE))
years <- yrs
nyears <- length(years)
nindices <- nages
abundanceW <- abundance_at_ageW_AM 
resid.pos.color <- "red"
resid.neg.color <- "blue"
resid.col=matrix(NA, nrow=nyears, ncol=nindices)   # set color for residual bubbles
resid.col <- ifelse(abundanceW > 0.0,resid.pos.color, resid.neg.color) 
# use a scaler to increase (>1) or decrease (0<val<1) the size of the bubbles to make them look good
scale.bubble <- 10/(max(abundanceW)) 
abundanceW1 <- abundanceW*scale.bubble
# first set up a blank space for the plot
plot(seq(0,nindices-1), rep(years[1],nindices), ylim=c(years[nyears],years[1]-2), xlab = "Age",
     ylab = "Year",  type = "n", axes=F,main='Abundance at Age West After Movement',xlim=c(-1,21),
sub=bquote("(Max Abundance"==.(max(abundance_at_ageW_AM))~")"))


axis(1, at= seq(0,nindices-1), lab=seq(0,nindices-1))
axis(2, at = seq(years[nyears],years[1], by=-5), lab = seq(years[nyears],years[1], by=-5),las=1)
box()
abline(h=seq( years[1], years[nyears]),col = "lightgray", lty = 1)
abline(a=1855,b=1,col = "lightgray", lty = 1)
abline(a=1860,b=1,col = "lightgray", lty = 1)
abline(a=1865,b=1,col = "lightgray", lty = 1)
abline(a=1870,b=1,col = "lightgray", lty = 1)
abline(a=1875,b=1,col = "lightgray", lty = 1)
abline(a=1880,b=1,col = "lightgray", lty = 1)
abline(a=1885,b=1,col = "lightgray", lty = 1)
abline(a=1890,b=1,col = "lightgray", lty = 1)
abline(a=1895,b=1,col = "lightgray", lty = 1)
abline(a=1900,b=1,col = "lightgray", lty = 1)
abline(a=1905,b=1,col = "lightgray", lty = 1)
abline(a=1910,b=1,col = "lightgray", lty = 1)
abline(a=1915,b=1,col = "lightgray", lty = 1)
abline(a=1920,b=1,col = "lightgray", lty = 1)
abline(a=1925,b=1,col = "lightgray", lty = 1)
abline(a=1930,b=1,col = "lightgray", lty = 1)
abline(a=1935,b=1,col = "lightgray", lty = 1)
abline(a=1940,b=1,col = "lightgray", lty = 1)
abline(a=1945,b=1,col = "lightgray", lty = 1)
abline(a=1950,b=1,col = "lightgray", lty = 1)
abline(a=1955,b=1,col = "lightgray", lty = 1)
abline(a=1960,b=1,col = "lightgray", lty = 1)
abline(a=1965,b=1,col = "lightgray", lty = 1)
abline(a=1970,b=1,col = "lightgray", lty = 1)
abline(a=1975,b=1,col = "lightgray", lty = 1)
abline(a=1980,b=1,col = "lightgray", lty = 1)
abline(a=1985,b=1,col = "lightgray", lty = 1)
abline(a=1990,b=1,col = "lightgray", lty = 1)
abline(a=1995,b=1,col = "lightgray", lty = 1)
abline(a=2000,b=1,col = "lightgray", lty = 1)
abline(a=2005,b=1,col = "lightgray", lty = 1)
abline(a=2010,b=1,col = "lightgray", lty = 1)
# fill in the bubbles using points command, cex determines size, bg determines color
for (i in 1:nyears){
  points(seq(0, nindices-1), rep((years[i]), nindices),  
  cex =abs(abundanceW1[i,]) ,  pch = 21, bg = resid.col[i,], col='black')
}
############################### Weight At Age  ###############################################

layout(matrix(c(1,2), 2, 1, byrow = TRUE))

matplot(yrs,weight_spawnE,type='l',main='Spawning Weight-at-Age East',lty=1,col=c(1:21),
xlab='Year',ylab='Weight(kg)',ylim=c(0,max(weight_spawnE)+.00001*max(weight_spawnE)),xlim=c(yrs[1],yrs[nyrs]+1)) 
grid()
textxy(rep((yrs[nyrs]+1),times=nages),weight_spawnE[nyrs,],c(0:(nages-1)),cex=.3)

matplot(yrs,weight_spawnW,type='l',main='Spawning Weight-at-Age West',lty=1,col=c(1:21),
xlab='Year',ylab='Weight(kg)',ylim=c(0,max(weight_spawnW)+.00001*max(weight_spawnW)),xlim=c(yrs[1],yrs[nyrs]+1)) 
grid()
textxy(rep((yrs[nyrs]+1),times=nages),weight_spawnW[nyrs,],c(0:(nages-1)),cex=.3)

################################ fit to IBM proportions ################################################

if(IBM_switch==1)
{
if(IBM_full_switch==0)
{
layout(matrix(c(1,1), 1, 1, byrow = TRUE))

matplot(yrs_IBM,IBM_recap_prop_totalE[,2],type='n',main='Fit to IBM Proportions East',xlab='Year',
ylab='Proportion of Fish in Each State',xlim=c(yrs_IBM[1]-1,yrs_IBM[nyrs_IBM]+1),
        ylim=c(0,1.1))
grid()
matpoints(yrs_IBM,OBS_IBM_recap_propE[,2],col="red",pch=16);
matlines(yrs_IBM,IBM_recap_prop_totalE[,2],col="blue",lwd=2);
matpoints(yrs_IBM,OBS_IBM_recap_propE[,1],col="red",pch=15);
matlines(yrs_IBM,IBM_recap_prop_totalE[,1],col="blue",lwd=2,lty=3);
matpoints(yrs_IBM,OBS_IBM_recap_propE[,3],col="red",pch=17);
matlines(yrs_IBM,IBM_recap_prop_totalE[,3],col="blue",lwd=2,lty=2);
legend("top",ncol=3, c("Resident (E_E)", "OBS Resident (E_E)", "Move (E_W)","OBS Move (E_W)","Dead","OBS Dead"), lty=c(3,-1,1,-1,2,-1),lwd=2,
col = c("blue","red","blue","red","blue","red"),pch=c(-1,15,-1,16,-1,17),cex=.6, merge=FALSE)


matplot(yrs_IBM,IBM_recap_prop_totalW[,2],type='n',main='Fit to IBM Proportions West',xlab='Year',
ylab='Proportion of Fish in Each State',xlim=c(yrs_IBM[1]-1,yrs_IBM[nyrs_IBM]+1),
        ylim=c(0,1.1))
grid()
matpoints(yrs_IBM,OBS_IBM_recap_propW[,2],col="red",pch=16);
matlines(yrs_IBM,IBM_recap_prop_totalW[,2],col="blue",lwd=2);
matpoints(yrs_IBM,OBS_IBM_recap_propW[,1],col="red",pch=15);
matlines(yrs_IBM,IBM_recap_prop_totalW[,1],col="blue",lwd=2,lty=3);
matpoints(yrs_IBM,OBS_IBM_recap_propW[,3],col="red",pch=17);
matlines(yrs_IBM,IBM_recap_prop_totalW[,3],col="blue",lwd=2,lty=2);
legend("top",ncol=3, c("Resident (W_W)", "OBS Resident (W_W)", "Move (W_E)","OBS Move (W_E)","Dead","OBS Dead"), lty=c(3,-1,1,-1,2,-1),lwd=2,
col = c("blue","red","blue","red","blue","red"),pch=c(-1,15,-1,16,-1,17),cex=.6, merge=FALSE)
}
if(IBM_full_switch==1)
{
layout(matrix(c(1,1), 1, 1, byrow = TRUE))

matplot(yrs,IBM_recap_prop_totalE_full[,2],type='n',main='Fit to IBM Proportions East',xlab='Year',
ylab='Proportion of Fish in Each State',xlim=c(yrs[1]-1,yrs[nyrs]+1),
        ylim=c(0,1.1))
grid()
matpoints(yrs,OBS_IBM_recap_propE_full[,2],col="red",pch=16);
matlines(yrs,IBM_recap_prop_totalE_full[,2],col="blue",lwd=2);
matpoints(yrs,OBS_IBM_recap_propE_full[,1],col="red",pch=15);
matlines(yrs,IBM_recap_prop_totalE_full[,1],col="blue",lwd=2,lty=3);
matpoints(yrs,OBS_IBM_recap_propE_full[,3],col="red",pch=17);
matlines(yrs,IBM_recap_prop_totalE_full[,3],col="blue",lwd=2,lty=2);
legend("top",ncol=3, c("Resident (E_E)", "OBS Resident (E_E)", "Move (E_W)","OBS Move (E_W)","Dead","OBS Dead"), lty=c(3,-1,1,-1,2,-1),lwd=2,
col = c("blue","red","blue","red","blue","red"),pch=c(-1,15,-1,16,-1,17),cex=.6, merge=FALSE)


matplot(yrs,IBM_recap_prop_totalW_full[,2],type='n',main='Fit to IBM Proportions West',xlab='Year',
ylab='Proportion of Fish in Each State',xlim=c(yrs[1]-1,yrs[nyrs]+1),
        ylim=c(0,1.1))
grid()
matpoints(yrs,OBS_IBM_recap_propW_full[,2],col="red",pch=16);
matlines(yrs,IBM_recap_prop_totalW_full[,2],col="blue",lwd=2);
matpoints(yrs,OBS_IBM_recap_propW_full[,1],col="red",pch=15);
matlines(yrs,IBM_recap_prop_totalW_full[,1],col="blue",lwd=2,lty=3);
matpoints(yrs,OBS_IBM_recap_propW_full[,3],col="red",pch=17);
matlines(yrs,IBM_recap_prop_totalW_full[,3],col="blue",lwd=2,lty=2);
legend("top",ncol=3, c("Resident (W_W)", "OBS Resident (W_W)", "Move (W_E)","OBS Move (W_E)","Dead","OBS Dead"), lty=c(3,-1,1,-1,2,-1),lwd=2,
col = c("blue","red","blue","red","blue","red"),pch=c(-1,15,-1,16,-1,17),cex=.6, merge=FALSE)
}
}

############################ Survey Catchability ##############################################

layout(matrix(c(1), 1, byrow = TRUE))
qsurvey<-as.data.frame(cbind(q_HL_E,q_HL_W,q_MRIP_E,q_MRIP_W,q_HBT_E,q_HBT_W,q_SHR_E,q_SHR_W,q_SUMM_E,q_SUMM_W,q_FALL_E,q_FALL_W))
qsurvey2<-stack(qsurvey)
names(qsurvey2)<-c("q","names")

matplot(1:nrow(qsurvey2),qsurvey2$q,col="blue",main='Catchability Coefficients',xaxt="n", xlab="",
ylab='q',ylim=c(0,max(qsurvey)+.01*max(qsurvey)),xlim=c(0,nrow(qsurvey2)+1),pch=16)
grid()
textxy(1:nrow(qsurvey2),qsurvey2$q,qsurvey2$q)
axis(1,at=1:nrow(qsurvey2),labels=qsurvey2$names,cex.axis=.7,tick=TRUE,las=2)


####################### survey fit and resids

layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow = TRUE))

matplot(yrs_index_HL_E,Index_HL_E,type='n',main='CPUE Handline East',
     xlab='Year',ylab='CPUE',xlim=c(min(yrs_index_HL_E),max(yrs_index_HL_E)),
ylim=c(0,(max(Index_HL_E,OBS_index_HL_E)+.1*max(Index_HL_E,OBS_index_HL_E)))) 
grid()
matlines(yrs_index_HL_E,Index_HL_E,col='blue',lwd=2);
matpoints(yrs_index_HL_E,OBS_index_HL_E,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)

matplot(yrs_index_HL_W,Index_HL_W,type='n',main='CPUE Handline West',
     xlab='Year',ylab='CPUE',xlim=c(min(yrs_index_HL_W),max(yrs_index_HL_W)),
ylim=c(0,(max(Index_HL_W,OBS_index_HL_W)+.1*max(Index_HL_W,OBS_index_HL_W)))) 
grid()
matlines(yrs_index_HL_W,Index_HL_W,col='blue',lwd=2);
matpoints(yrs_index_HL_W,OBS_index_HL_W,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)

matplot(yrs_index_MRIP_E,Index_MRIP_E,type='n',main='CPUE MRIP East',
     xlab='Year',ylab='CPUE',xlim=c(min(yrs_index_MRIP_E),max(yrs_index_MRIP_E)),
ylim=c(0,(max(Index_MRIP_E,OBS_index_MRIP_E)+.1*max(Index_MRIP_E,OBS_index_MRIP_E)))) 
grid()
matlines(yrs_index_MRIP_E,Index_MRIP_E,col='blue',lwd=2);
matpoints(yrs_index_MRIP_E,OBS_index_MRIP_E,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)

matplot(yrs_index_MRIP_W,Index_MRIP_W,type='n',main='CPUE MRIP West',
     xlab='Year',ylab='CPUE',xlim=c(min(yrs_index_MRIP_W),max(yrs_index_MRIP_W)),
ylim=c(0,(max(Index_MRIP_W,OBS_index_MRIP_W)+.1*max(Index_MRIP_W,OBS_index_MRIP_W)))) 
grid()
matlines(yrs_index_MRIP_W,Index_MRIP_W,col='blue',lwd=2);
matpoints(yrs_index_MRIP_W,OBS_index_MRIP_W,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)


matplot(yrs_index_HBT_E,Index_HBT_E,type='n',main='CPUE Headboat East',
     xlab='Year',ylab='CPUE',xlim=c(min(yrs_index_HBT_E),max(yrs_index_HBT_E)),
ylim=c(0,(max(Index_HBT_E,OBS_index_HBT_E)+.1*max(Index_HBT_E,OBS_index_HBT_E)))) 
grid()
matlines(yrs_index_HBT_E,Index_HBT_E,col='blue',lwd=2);
matpoints(yrs_index_HBT_E,OBS_index_HBT_E,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)

matplot(yrs_index_HBT_W,Index_HBT_W,type='n',main='CPUE Headboat West',
     xlab='Year',ylab='CPUE',xlim=c(min(yrs_index_HBT_W),max(yrs_index_HBT_W)),
ylim=c(0,(max(Index_HBT_W,OBS_index_HBT_W)+.1*max(Index_HBT_W,OBS_index_HBT_W)))) 
grid()
matlines(yrs_index_HBT_W,Index_HBT_W,col='blue',lwd=2);
matpoints(yrs_index_HBT_W,OBS_index_HBT_W,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)

matplot(yrs_effort_SHR_E,effort_SHR_E,type='n',main='Effort Shrimp East',
     xlab='Year',ylab='Effort',xlim=c(min(yrs_effort_SHR_E),max(yrs_effort_SHR_E)),
ylim=c(0,(max(effort_SHR_E,OBS_effort_SHR_E)+.1*max(effort_SHR_E,OBS_effort_SHR_E)))) 
grid()
matlines(yrs_effort_SHR_E,effort_SHR_E,col='blue',lwd=2);
matpoints(yrs_effort_SHR_E,OBS_effort_SHR_E,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)

matplot(yrs_effort_SHR_W,effort_SHR_W,type='n',main='Effort Shrimp West',
     xlab='Year',ylab='Effort',xlim=c(min(yrs_effort_SHR_W),max(yrs_effort_SHR_W)),
ylim=c(0,(max(effort_SHR_W,OBS_effort_SHR_W)+.1*max(effort_SHR_W,OBS_effort_SHR_W)))) 
grid()
matlines(yrs_effort_SHR_W,effort_SHR_W,col='blue',lwd=2);
matpoints(yrs_effort_SHR_W,OBS_effort_SHR_W,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)



layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))

matplot(yrs_index_SUMM_E,Index_SUMM_E,type='n',main='CPUE Summer East',
     xlab='Year',ylab='CPUE',xlim=c(min(yrs_index_SUMM_E),max(yrs_index_SUMM_E)),
ylim=c(0,(max(Index_SUMM_E,OBS_index_SUMM_E)+.1*max(Index_SUMM_E,OBS_index_SUMM_E)))) 
grid()
matlines(yrs_index_SUMM_E,Index_SUMM_E,col='blue',lwd=2);
matpoints(yrs_index_SUMM_E,OBS_index_SUMM_E,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)

matplot(yrs_index_SUMM_W,Index_SUMM_W,type='n',main='CPUE Summer West',
     xlab='Year',ylab='CPUE',xlim=c(min(yrs_index_SUMM_W),max(yrs_index_SUMM_W)),
ylim=c(0,(max(Index_SUMM_W,OBS_index_SUMM_W)+.1*max(Index_SUMM_W,OBS_index_SUMM_W)))) 
grid()
matlines(yrs_index_SUMM_W,Index_SUMM_W,col='blue',lwd=2);
matpoints(yrs_index_SUMM_W,OBS_index_SUMM_W,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)

matplot(yrs_index_FALL_E,Index_FALL_E,type='n',main='CPUE Fall East',
     xlab='Year',ylab='CPUE',xlim=c(min(yrs_index_FALL_E),max(yrs_index_FALL_E)),
ylim=c(0,(max(Index_FALL_E,OBS_index_FALL_E)+.1*max(Index_FALL_E,OBS_index_FALL_E)))) 
grid()
matlines(yrs_index_FALL_E,Index_FALL_E,col='blue',lwd=2);
matpoints(yrs_index_FALL_E,OBS_index_FALL_E,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)

matplot(yrs_index_FALL_W,Index_FALL_W,type='n',main='CPUE Fall West',
     xlab='Year',ylab='CPUE',xlim=c(min(yrs_index_FALL_W),max(yrs_index_FALL_W)),
ylim=c(0,(max(Index_FALL_W,OBS_index_FALL_W)+.1*max(Index_FALL_W,OBS_index_FALL_W)))) 
grid()
matlines(yrs_index_FALL_W,Index_FALL_W,col='blue',lwd=2);
matpoints(yrs_index_FALL_W,OBS_index_FALL_W,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)

layout(matrix(c(1,2), 2, 1, byrow = TRUE),
   widths=c(1,1), heights=c(1,1))

matplot(yrs_index_HL_E,Index_HL_E,type='n',main='CPUE Handline East',
     xlab='Year',ylab='CPUE',xlim=c(min(yrs_index_HL_E),max(yrs_index_HL_E)),
ylim=c(0,(max(Index_HL_E,OBS_index_HL_E)+.1*max(Index_HL_E,OBS_index_HL_E)))) 
grid()
matlines(yrs_index_HL_E,Index_HL_E,col='blue',lwd=2);
matpoints(yrs_index_HL_E,OBS_index_HL_E,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)


ln_RMSE_index_HL_E<-sqrt(mean((log(Index_HL_E)-log(OBS_index_HL_E))^2)) #lognormal dist. so resids are ln(obs)-ln(pred) and so want SD of difference in log space, also sd function does not account for n in r so use RMSE
stan_res_index_HL_E<-res_index_HL_E/ln_RMSE_index_HL_E

matplot(yrs_index_HL_E,stan_res_index_HL_E,type='n',xlab='Year',sub=bquote("(RMSE"==.(ln_RMSE_index_HL_E)~")"),
        ylab='[ln(Obs)-ln(Pred)]/SD',xlim=c(min(yrs_index_HL_E)-1,max(yrs_index_HL_E)+1),
        ylim=c(min(stan_res_index_HL_E)-.05*abs(min(stan_res_index_HL_E)),max(stan_res_index_HL_E)+.05*abs(max(stan_res_index_HL_E)))) 
grid()
matlines(x=((min(yrs_index_HL_E)-1):(max(yrs_index_HL_E)+1)),y=rep(0,times=nyrs_index_HL_E+2),col="blue",lwd=2);
matpoints(yrs_index_HL_E,stan_res_index_HL_E,col="red",pch=19);
lines(x=c(yrs_index_HL_E,yrs_index_HL_E),y=c(stan_res_index_HL_E,rep(0,nyrs_index_HL_E)),type='h');
title("CPUE Standardized Residuals Handline East",line=1)



matplot(yrs_index_HL_W,Index_HL_W,type='n',main='CPUE Handline West',
     xlab='Year',ylab='CPUE',xlim=c(min(yrs_index_HL_W),max(yrs_index_HL_W)),
ylim=c(0,(max(Index_HL_W,OBS_index_HL_W)+.1*max(Index_HL_W,OBS_index_HL_W)))) 
grid()
matlines(yrs_index_HL_W,Index_HL_W,col='blue',lwd=2);
matpoints(yrs_index_HL_W,OBS_index_HL_W,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)


ln_RMSE_index_HL_W<-sqrt(mean((log(Index_HL_W)-log(OBS_index_HL_W))^2)) #lognormal dist. so resids are ln(obs)-ln(pred) and so want SD of difference in log space, also sd function does not account for n in r so use RMSE
stan_res_index_HL_W<-res_index_HL_W/ln_RMSE_index_HL_W

matplot(yrs_index_HL_W,stan_res_index_HL_W,type='n',xlab='Year',sub=bquote("(RMSE"==.(ln_RMSE_index_HL_W)~")"),
        ylab='[ln(Obs)-ln(Pred)]/SD',xlim=c(min(yrs_index_HL_W)-1,max(yrs_index_HL_W)+1),
        ylim=c(min(stan_res_index_HL_W)-.05*abs(min(stan_res_index_HL_W)),max(stan_res_index_HL_W)+.05*abs(max(stan_res_index_HL_W)))) 
grid()
matlines(x=((min(yrs_index_HL_W)-1):(max(yrs_index_HL_W)+1)),y=rep(0,times=nyrs_index_HL_W+2),col="blue",lwd=2);
matpoints(yrs_index_HL_W,stan_res_index_HL_W,col="red",pch=19);
lines(x=c(yrs_index_HL_W,yrs_index_HL_W),y=c(stan_res_index_HL_W,rep(0,nyrs_index_HL_W)),type='h');
title("CPUE Standardized Residuals Handline West",line=1)


matplot(yrs_index_MRIP_E,Index_MRIP_E,type='n',main='CPUE MRIP East',
     xlab='Year',ylab='CPUE',xlim=c(min(yrs_index_MRIP_E),max(yrs_index_MRIP_E)),
ylim=c(0,(max(Index_MRIP_E,OBS_index_MRIP_E)+.1*max(Index_MRIP_E,OBS_index_MRIP_E)))) 
grid()
matlines(yrs_index_MRIP_E,Index_MRIP_E,col='blue',lwd=2);
matpoints(yrs_index_MRIP_E,OBS_index_MRIP_E,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)


ln_RMSE_index_MRIP_E<-sqrt(mean((log(Index_MRIP_E)-log(OBS_index_MRIP_E))^2)) #lognormal dist. so resids are ln(obs)-ln(pred) and so want SD of difference in log space, also sd function does not account for n in r so use RMSE
stan_res_index_MRIP_E<-res_index_MRIP_E/ln_RMSE_index_MRIP_E

matplot(yrs_index_MRIP_E,stan_res_index_MRIP_E,type='n',xlab='Year',sub=bquote("(RMSE"==.(ln_RMSE_index_MRIP_E)~")"),
        ylab='[ln(Obs)-ln(Pred)]/SD',xlim=c(min(yrs_index_MRIP_E)-1,max(yrs_index_MRIP_E)+1),
        ylim=c(min(stan_res_index_MRIP_E)-.05*abs(min(stan_res_index_MRIP_E)),max(stan_res_index_MRIP_E)+.05*abs(max(stan_res_index_MRIP_E)))) 
grid()
matlines(x=((min(yrs_index_MRIP_E)-1):(max(yrs_index_MRIP_E)+1)),y=rep(0,times=nyrs_index_MRIP_E+2),col="blue",lwd=2);
matpoints(yrs_index_MRIP_E,stan_res_index_MRIP_E,col="red",pch=19);
lines(x=c(yrs_index_MRIP_E,yrs_index_MRIP_E),y=c(stan_res_index_MRIP_E,rep(0,nyrs_index_MRIP_E)),type='h');
title("CPUE Standardized Residuals MRIP East",line=1)



matplot(yrs_index_MRIP_W,Index_MRIP_W,type='n',main='CPUE MRIP West',
     xlab='Year',ylab='CPUE',xlim=c(min(yrs_index_MRIP_W),max(yrs_index_MRIP_W)),
ylim=c(0,(max(Index_MRIP_W,OBS_index_MRIP_W)+.1*max(Index_MRIP_W,OBS_index_MRIP_W)))) 
grid()
matlines(yrs_index_MRIP_W,Index_MRIP_W,col='blue',lwd=2);
matpoints(yrs_index_MRIP_W,OBS_index_MRIP_W,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)


ln_RMSE_index_MRIP_W<-sqrt(mean((log(Index_MRIP_W)-log(OBS_index_MRIP_W))^2)) #lognormal dist. so resids are ln(obs)-ln(pred) and so want SD of difference in log space, also sd function does not account for n in r so use RMSE
stan_res_index_MRIP_W<-res_index_MRIP_W/ln_RMSE_index_MRIP_W

matplot(yrs_index_MRIP_W,stan_res_index_MRIP_W,type='n',xlab='Year',sub=bquote("(RMSE"==.(ln_RMSE_index_MRIP_W)~")"),
        ylab='[ln(Obs)-ln(Pred)]/SD',xlim=c(min(yrs_index_MRIP_W)-1,max(yrs_index_MRIP_W)+1),
        ylim=c(min(stan_res_index_MRIP_W)-.05*abs(min(stan_res_index_MRIP_W)),max(stan_res_index_MRIP_W)+.05*abs(max(stan_res_index_MRIP_W)))) 
grid()
matlines(x=((min(yrs_index_MRIP_W)-1):(max(yrs_index_MRIP_W)+1)),y=rep(0,times=nyrs_index_MRIP_W+2),col="blue",lwd=2);
matpoints(yrs_index_MRIP_W,stan_res_index_MRIP_W,col="red",pch=19);
lines(x=c(yrs_index_MRIP_W,yrs_index_MRIP_W),y=c(stan_res_index_MRIP_W,rep(0,nyrs_index_MRIP_W)),type='h');
title("CPUE Standardized Residuals MRIP West",line=1)



matplot(yrs_index_HBT_E,Index_HBT_E,type='n',main='CPUE Headboat East',
     xlab='Year',ylab='CPUE',xlim=c(min(yrs_index_HBT_E),max(yrs_index_HBT_E)),
ylim=c(0,(max(Index_HBT_E,OBS_index_HBT_E)+.1*max(Index_HBT_E,OBS_index_HBT_E)))) 
grid()
matlines(yrs_index_HBT_E,Index_HBT_E,col='blue',lwd=2);
matpoints(yrs_index_HBT_E,OBS_index_HBT_E,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)


ln_RMSE_index_HBT_E<-sqrt(mean((log(Index_HBT_E)-log(OBS_index_HBT_E))^2)) #lognormal dist. so resids are ln(obs)-ln(pred) and so want SD of difference in log space, also sd function does not account for n in r so use RMSE
stan_res_index_HBT_E<-res_index_HBT_E/ln_RMSE_index_HBT_E

matplot(yrs_index_HBT_E,stan_res_index_HBT_E,type='n',xlab='Year',sub=bquote("(RMSE"==.(ln_RMSE_index_HBT_E)~")"),
        ylab='[ln(Obs)-ln(Pred)]/SD',xlim=c(min(yrs_index_HBT_E)-1,max(yrs_index_HBT_E)+1),
        ylim=c(min(stan_res_index_HBT_E)-.05*abs(min(stan_res_index_HBT_E)),max(stan_res_index_HBT_E)+.05*abs(max(stan_res_index_HBT_E)))) 
grid()
matlines(x=((min(yrs_index_HBT_E)-1):(max(yrs_index_HBT_E)+1)),y=rep(0,times=nyrs_index_HBT_E+2),col="blue",lwd=2);
matpoints(yrs_index_HBT_E,stan_res_index_HBT_E,col="red",pch=19);
lines(x=c(yrs_index_HBT_E,yrs_index_HBT_E),y=c(stan_res_index_HBT_E,rep(0,nyrs_index_HBT_E)),type='h');
title("CPUE Standardized Residuals Headboat East",line=1)



matplot(yrs_index_HBT_W,Index_HBT_W,type='n',main='CPUE Headboat West',
     xlab='Year',ylab='CPUE',xlim=c(min(yrs_index_HBT_W),max(yrs_index_HBT_W)),
ylim=c(0,(max(Index_HBT_W,OBS_index_HBT_W)+.1*max(Index_HBT_W,OBS_index_HBT_W)))) 
grid()
matlines(yrs_index_HBT_W,Index_HBT_W,col='blue',lwd=2);
matpoints(yrs_index_HBT_W,OBS_index_HBT_W,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)


ln_RMSE_index_HBT_W<-sqrt(mean((log(Index_HBT_W)-log(OBS_index_HBT_W))^2)) #lognormal dist. so resids are ln(obs)-ln(pred) and so want SD of difference in log space, also sd function does not account for n in r so use RMSE
stan_res_index_HBT_W<-res_index_HBT_W/ln_RMSE_index_HBT_W

matplot(yrs_index_HBT_W,stan_res_index_HBT_W,type='n',xlab='Year',sub=bquote("(RMSE"==.(ln_RMSE_index_HBT_W)~")"),
        ylab='[ln(Obs)-ln(Pred)]/SD',xlim=c(min(yrs_index_HBT_W)-1,max(yrs_index_HBT_W)+1),
        ylim=c(min(stan_res_index_HBT_W)-.05*abs(min(stan_res_index_HBT_W)),max(stan_res_index_HBT_W)+.05*abs(max(stan_res_index_HBT_W)))) 
grid()
matlines(x=((min(yrs_index_HBT_W)-1):(max(yrs_index_HBT_W)+1)),y=rep(0,times=nyrs_index_HBT_W+2),col="blue",lwd=2);
matpoints(yrs_index_HBT_W,stan_res_index_HBT_W,col="red",pch=19);
lines(x=c(yrs_index_HBT_W,yrs_index_HBT_W),y=c(stan_res_index_HBT_W,rep(0,nyrs_index_HBT_W)),type='h');
title("CPUE Standardized Residuals Headboat West",line=1)


matplot(yrs_effort_SHR_E,effort_SHR_E,type='n',main='Effort Shrimp East',
     xlab='Year',ylab='Effort',xlim=c(min(yrs_effort_SHR_E),max(yrs_effort_SHR_E)),
ylim=c(0,(max(effort_SHR_E,OBS_effort_SHR_E)+.1*max(effort_SHR_E,OBS_effort_SHR_E)))) 
grid()
matlines(yrs_effort_SHR_E,effort_SHR_E,col='blue',lwd=2);
matpoints(yrs_effort_SHR_E,OBS_effort_SHR_E,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)


ln_RMSE_effort_SHR_E<-sqrt(mean((log(effort_SHR_E)-log(OBS_effort_SHR_E))^2)) #lognormal dist. so resids are ln(obs)-ln(pred) and so want SD of difference in log space, also sd function does not account for n in r so use RMSE
stan_res_effort_SHR_E<-res_effort_SHR_E/ln_RMSE_effort_SHR_E

matplot(yrs_effort_SHR_E,stan_res_effort_SHR_E,type='n',xlab='Year',sub=bquote("(RMSE"==.(ln_RMSE_effort_SHR_E)~")"),
        ylab='[ln(Obs)-ln(Pred)]/SD',xlim=c(min(yrs_effort_SHR_E)-1,max(yrs_effort_SHR_E)+1),
        ylim=c(min(stan_res_effort_SHR_E)-.05*abs(min(stan_res_effort_SHR_E)),max(stan_res_effort_SHR_E)+.05*abs(max(stan_res_effort_SHR_E)))) 
grid()
matlines(x=((min(yrs_effort_SHR_E)-1):(max(yrs_effort_SHR_E)+1)),y=rep(0,times=nyrs_effort_SHR_E+2),col="blue",lwd=2);
matpoints(yrs_effort_SHR_E,stan_res_effort_SHR_E,col="red",pch=19);
lines(x=c(yrs_effort_SHR_E,yrs_effort_SHR_E),y=c(stan_res_effort_SHR_E,rep(0,nyrs_effort_SHR_E)),type='h');
title("Effort Standardized Residuals Shrimp East",line=1)

matplot(yrs_effort_SHR_W,effort_SHR_W,type='n',main='Effort Shrimp West',
     xlab='Year',ylab='Effort',xlim=c(min(yrs_effort_SHR_W),max(yrs_effort_SHR_W)),
ylim=c(0,(max(effort_SHR_W,OBS_effort_SHR_W)+.1*max(effort_SHR_W,OBS_effort_SHR_W)))) 
grid()
matlines(yrs_effort_SHR_W,effort_SHR_W,col='blue',lwd=2);
matpoints(yrs_effort_SHR_W,OBS_effort_SHR_W,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)


ln_RMSE_effort_SHR_W<-sqrt(mean((log(effort_SHR_W)-log(OBS_effort_SHR_W))^2)) #lognormal dist. so resids are ln(obs)-ln(pred) and so want SD of difference in log space, also sd function does not account for n in r so use RMSE
stan_res_effort_SHR_W<-res_effort_SHR_W/ln_RMSE_effort_SHR_W

matplot(yrs_effort_SHR_W,stan_res_effort_SHR_W,type='n',xlab='Year',sub=bquote("(RMSE"==.(ln_RMSE_effort_SHR_W)~")"),
        ylab='[ln(Obs)-ln(Pred)]/SD',xlim=c(min(yrs_effort_SHR_W)-1,max(yrs_effort_SHR_W)+1),
        ylim=c(min(stan_res_effort_SHR_W)-.05*abs(min(stan_res_effort_SHR_W)),max(stan_res_effort_SHR_W)+.05*abs(max(stan_res_effort_SHR_W)))) 
grid()
matlines(x=((min(yrs_effort_SHR_W)-1):(max(yrs_effort_SHR_W)+1)),y=rep(0,times=nyrs_effort_SHR_W+2),col="blue",lwd=2);
matpoints(yrs_effort_SHR_W,stan_res_effort_SHR_W,col="red",pch=19);
lines(x=c(yrs_effort_SHR_W,yrs_effort_SHR_W),y=c(stan_res_effort_SHR_W,rep(0,nyrs_effort_SHR_W)),type='h');
title("Effort Standardized Residuals Shrimp West",line=1)



layout(matrix(c(1,2), 2, 1, byrow = TRUE),
   widths=c(1,1), heights=c(1,1))

matplot(yrs_index_SUMM_E,Index_SUMM_E,type='n',main='CPUE Summer East',
     xlab='Year',ylab='CPUE',xlim=c(min(yrs_index_SUMM_E),max(yrs_index_SUMM_E)),
ylim=c(0,(max(Index_SUMM_E,OBS_index_SUMM_E)+.1*max(Index_SUMM_E,OBS_index_SUMM_E)))) 
grid()
matlines(yrs_index_SUMM_E,Index_SUMM_E,col='blue',lwd=2);
matpoints(yrs_index_SUMM_E,OBS_index_SUMM_E,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)


ln_RMSE_index_SUMM_E<-sqrt(mean((log(Index_SUMM_E)-log(OBS_index_SUMM_E))^2)) #lognormal dist. so resids are ln(obs)-ln(pred) and so want SD of difference in log space, also sd function does not account for n in r so use RMSE
stan_res_index_SUMM_E<-res_index_SUMM_E/ln_RMSE_index_SUMM_E

matplot(yrs_index_SUMM_E,stan_res_index_SUMM_E,type='n',xlab='Year',sub=bquote("(RMSE"==.(ln_RMSE_index_SUMM_E)~")"),
        ylab='[ln(Obs)-ln(Pred)]/SD',xlim=c(min(yrs_index_SUMM_E)-1,max(yrs_index_SUMM_E)+1),
        ylim=c(min(stan_res_index_SUMM_E)-.05*abs(min(stan_res_index_SUMM_E)),max(stan_res_index_SUMM_E)+.05*abs(max(stan_res_index_SUMM_E)))) 
grid()
matlines(x=((min(yrs_index_SUMM_E)-1):(max(yrs_index_SUMM_E)+1)),y=rep(0,times=nyrs_index_SUMM_E+2),col="blue",lwd=2);
matpoints(yrs_index_SUMM_E,stan_res_index_SUMM_E,col="red",pch=19);
lines(x=c(yrs_index_SUMM_E,yrs_index_SUMM_E),y=c(stan_res_index_SUMM_E,rep(0,nyrs_index_SUMM_E)),type='h');
title("CPUE Standardized Residuals Summer East",line=1)

layout(matrix(c(1,2), 2, 1, byrow = TRUE),
   widths=c(1,1), heights=c(1,1))

matplot(yrs_index_SUMM_W,Index_SUMM_W,type='n',main='CPUE Summer West',
     xlab='Year',ylab='CPUE',xlim=c(min(yrs_index_SUMM_W),max(yrs_index_SUMM_W)),
ylim=c(0,(max(Index_SUMM_W,OBS_index_SUMM_W)+.1*max(Index_SUMM_W,OBS_index_SUMM_W)))) 
grid()
matlines(yrs_index_SUMM_W,Index_SUMM_W,col='blue',lwd=2);
matpoints(yrs_index_SUMM_W,OBS_index_SUMM_W,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)


ln_RMSE_index_SUMM_W<-sqrt(mean((log(Index_SUMM_W)-log(OBS_index_SUMM_W))^2)) #lognormal dist. so resids are ln(obs)-ln(pred) and so want SD of difference in log space, also sd function does not account for n in r so use RMSE
stan_res_index_SUMM_W<-res_index_SUMM_W/ln_RMSE_index_SUMM_W

matplot(yrs_index_SUMM_W,stan_res_index_SUMM_W,type='n',xlab='Year',sub=bquote("(RMSE"==.(ln_RMSE_index_SUMM_W)~")"),
        ylab='[ln(Obs)-ln(Pred)]/SD',xlim=c(min(yrs_index_SUMM_W)-1,max(yrs_index_SUMM_W)+1),
        ylim=c(min(stan_res_index_SUMM_W)-.05*abs(min(stan_res_index_SUMM_W)),max(stan_res_index_SUMM_W)+.05*abs(max(stan_res_index_SUMM_W)))) 
grid()
matlines(x=((min(yrs_index_SUMM_W)-1):(max(yrs_index_SUMM_W)+1)),y=rep(0,times=nyrs_index_SUMM_W+2),col="blue",lwd=2);
matpoints(yrs_index_SUMM_W,stan_res_index_SUMM_W,col="red",pch=19);
lines(x=c(yrs_index_SUMM_W,yrs_index_SUMM_W),y=c(stan_res_index_SUMM_W,rep(0,nyrs_index_SUMM_W)),type='h');
title("CPUE Standardized Residuals Summer West",line=1)




layout(matrix(c(1,2), 2, 1, byrow = TRUE),
   widths=c(1,1), heights=c(1,1))

matplot(yrs_index_FALL_E,Index_FALL_E,type='n',main='CPUE Fall East',
     xlab='Year',ylab='CPUE',xlim=c(min(yrs_index_FALL_E),max(yrs_index_FALL_E)),
ylim=c(0,(max(Index_FALL_E,OBS_index_FALL_E)+.1*max(Index_FALL_E,OBS_index_FALL_E)))) 
grid()
matlines(yrs_index_FALL_E,Index_FALL_E,col='blue',lwd=2);
matpoints(yrs_index_FALL_E,OBS_index_FALL_E,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)


ln_RMSE_index_FALL_E<-sqrt(mean((log(Index_FALL_E)-log(OBS_index_FALL_E))^2)) #lognormal dist. so resids are ln(obs)-ln(pred) and so want SD of difference in log space, also sd function does not account for n in r so use RMSE
stan_res_index_FALL_E<-res_index_FALL_E/ln_RMSE_index_FALL_E

matplot(yrs_index_FALL_E,stan_res_index_FALL_E,type='n',xlab='Year',sub=bquote("(RMSE"==.(ln_RMSE_index_FALL_E)~")"),
        ylab='[ln(Obs)-ln(Pred)]/SD',xlim=c(min(yrs_index_FALL_E)-1,max(yrs_index_FALL_E)+1),
        ylim=c(min(stan_res_index_FALL_E)-.05*abs(min(stan_res_index_FALL_E)),max(stan_res_index_FALL_E)+.05*abs(max(stan_res_index_FALL_E)))) 
grid()
matlines(x=((min(yrs_index_FALL_E)-1):(max(yrs_index_FALL_E)+1)),y=rep(0,times=nyrs_index_FALL_E+2),col="blue",lwd=2);
matpoints(yrs_index_FALL_E,stan_res_index_FALL_E,col="red",pch=19);
lines(x=c(yrs_index_FALL_E,yrs_index_FALL_E),y=c(stan_res_index_FALL_E,rep(0,nyrs_index_FALL_E)),type='h');
title("CPUE Standardized Residuals Fall East",line=1)

layout(matrix(c(1,2), 2, 1, byrow = TRUE),
   widths=c(1,1), heights=c(1,1))

matplot(yrs_index_FALL_W,Index_FALL_W,type='n',main='CPUE Fall West',
     xlab='Year',ylab='CPUE',xlim=c(min(yrs_index_FALL_W),max(yrs_index_FALL_W)),
ylim=c(0,(max(Index_FALL_W,OBS_index_FALL_W)+.1*max(Index_FALL_W,OBS_index_FALL_W)))) 
grid()
matlines(yrs_index_FALL_W,Index_FALL_W,col='blue',lwd=2);
matpoints(yrs_index_FALL_W,OBS_index_FALL_W,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)


ln_RMSE_index_FALL_W<-sqrt(mean((log(Index_FALL_W)-log(OBS_index_FALL_W))^2)) #lognormal dist. so resids are ln(obs)-ln(pred) and so want SD of difference in log space, also sd function does not account for n in r so use RMSE
stan_res_index_FALL_W<-res_index_FALL_W/ln_RMSE_index_FALL_W

matplot(yrs_index_FALL_W,stan_res_index_FALL_W,type='n',xlab='Year',sub=bquote("(RMSE"==.(ln_RMSE_index_FALL_W)~")"),
        ylab='[ln(Obs)-ln(Pred)]/SD',xlim=c(min(yrs_index_FALL_W)-1,max(yrs_index_FALL_W)+1),
        ylim=c(min(stan_res_index_FALL_W)-.05*abs(min(stan_res_index_FALL_W)),max(stan_res_index_FALL_W)+.05*abs(max(stan_res_index_FALL_W)))) 
grid()
matlines(x=((min(yrs_index_FALL_W)-1):(max(yrs_index_FALL_W)+1)),y=rep(0,times=nyrs_index_FALL_W+2),col="blue",lwd=2);
matpoints(yrs_index_FALL_W,stan_res_index_FALL_W,col="red",pch=19);
lines(x=c(yrs_index_FALL_W,yrs_index_FALL_W),y=c(stan_res_index_FALL_W,rep(0,nyrs_index_FALL_W)),type='h');
title("CPUE Standardized Residuals Fall West",line=1)


################### landings fit and resids


layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow = TRUE))

matplot(yrs_landings_HL_E,pred_landings_HL_E,type='n',main='Landings Handline East',
     xlab='Year',ylab='Landings (mt)',xlim=c(min(yrs_landings_HL_E),max(yrs_landings_HL_E)),
ylim=c(0,(max(pred_landings_HL_E,OBS_landings_HL_E)+.1*max(pred_landings_HL_E,OBS_landings_HL_E)))) 
grid()
matlines(yrs_landings_HL_E,pred_landings_HL_E,col='blue',lwd=2);
matpoints(yrs_landings_HL_E,OBS_landings_HL_E,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)

matplot(yrs_landings_HL_W,pred_landings_HL_W,type='n',main='Landings Handline West',
     xlab='Year',ylab='Landings (mt)',xlim=c(min(yrs_landings_HL_W),max(yrs_landings_HL_W)),
ylim=c(0,(max(pred_landings_HL_W,OBS_landings_HL_W)+.1*max(pred_landings_HL_W,OBS_landings_HL_W)))) 
grid()
matlines(yrs_landings_HL_W,pred_landings_HL_W,col='blue',lwd=2);
matpoints(yrs_landings_HL_W,OBS_landings_HL_W,col='red',pch=16);
legend("bottom",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)

matplot(yrs_landings_LL_E,pred_landings_LL_E,type='n',main='Landings Longline East',
     xlab='Year',ylab='Landings (mt)',xlim=c(min(yrs_landings_LL_E),max(yrs_landings_LL_E)),
ylim=c(0,(max(pred_landings_LL_E,OBS_landings_LL_E)+.1*max(pred_landings_LL_E,OBS_landings_LL_E)))) 
grid()
matlines(yrs_landings_LL_E,pred_landings_LL_E,col='blue',lwd=2);
matpoints(yrs_landings_LL_E,OBS_landings_LL_E,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)

matplot(yrs_landings_LL_W,pred_landings_LL_W,type='n',main='Landings Longline West',
     xlab='Year',ylab='Landings (mt)',xlim=c(min(yrs_landings_LL_W),max(yrs_landings_LL_W)),
ylim=c(0,(max(pred_landings_LL_W,OBS_landings_LL_W)+.1*max(pred_landings_LL_W,OBS_landings_LL_W)))) 
grid()
matlines(yrs_landings_LL_W,pred_landings_LL_W,col='blue',lwd=2);
matpoints(yrs_landings_LL_W,OBS_landings_LL_W,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)

matplot(yrs_landings_MRIP_E,pred_landings_MRIP_E,type='n',main='Landings MRIP East',
     xlab='Year',ylab='Landings (1000s of Fish)',xlim=c(min(yrs_landings_MRIP_E),max(yrs_landings_MRIP_E)),
ylim=c(0,(max(pred_landings_MRIP_E,OBS_landings_MRIP_E)+.1*max(pred_landings_MRIP_E,OBS_landings_MRIP_E)))) 
grid()
matlines(yrs_landings_MRIP_E,pred_landings_MRIP_E,col='blue',lwd=2);
matpoints(yrs_landings_MRIP_E,OBS_landings_MRIP_E,col='red',pch=16);
legend("bottom",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)

matplot(yrs_landings_MRIP_W,pred_landings_MRIP_W,type='n',main='Landings MRIP West',
     xlab='Year',ylab='Landings (1000s of Fish)',xlim=c(min(yrs_landings_MRIP_W),max(yrs_landings_MRIP_W)),
ylim=c(0,(max(pred_landings_MRIP_W,OBS_landings_MRIP_W)+.1*max(pred_landings_MRIP_W,OBS_landings_MRIP_W)))) 
grid()
matlines(yrs_landings_MRIP_W,pred_landings_MRIP_W,col='blue',lwd=2);
matpoints(yrs_landings_MRIP_W,OBS_landings_MRIP_W,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)

matplot(yrs_landings_HBT_E,pred_landings_HBT_E,type='n',main='Landings Headboat East',
     xlab='Year',ylab='Landings (1000s of Fish)',xlim=c(min(yrs_landings_HBT_E),max(yrs_landings_HBT_E)),
ylim=c(0,(max(pred_landings_HBT_E,OBS_landings_HBT_E)+.1*max(pred_landings_HBT_E,OBS_landings_HBT_E)))) 
grid()
matlines(yrs_landings_HBT_E,pred_landings_HBT_E,col='blue',lwd=2);
matpoints(yrs_landings_HBT_E,OBS_landings_HBT_E,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)

matplot(yrs_landings_HBT_W,pred_landings_HBT_W,type='n',main='Landings Headboat West',
     xlab='Year',ylab='Landings (1000s of Fish)',xlim=c(min(yrs_landings_HBT_W),max(yrs_landings_HBT_W)),
ylim=c(0,(max(pred_landings_HBT_W,OBS_landings_HBT_W)+.1*max(pred_landings_HBT_W,OBS_landings_HBT_W)))) 
grid()
matlines(yrs_landings_HBT_W,pred_landings_HBT_W,col='blue',lwd=2);
matpoints(yrs_landings_HBT_W,OBS_landings_HBT_W,col='red',pch=16);
legend("topright",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)


layout(matrix(c(1,2), 2, 1, byrow = TRUE),
   widths=c(1,1), heights=c(1,1))

matplot(yrs_landings_HL_E,pred_landings_HL_E,type='n',main='Landings Handline East',
     xlab='Year',ylab='Landings (mt)',xlim=c(min(yrs_landings_HL_E),max(yrs_landings_HL_E)),
ylim=c(0,(max(pred_landings_HL_E,OBS_landings_HL_E)+.1*max(pred_landings_HL_E,OBS_landings_HL_E)))) 
grid()
matlines(yrs_landings_HL_E,pred_landings_HL_E,col='blue',lwd=2);
matpoints(yrs_landings_HL_E,OBS_landings_HL_E,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)


ln_RMSE_landings_HL_E<-sqrt(mean((log(pred_landings_HL_E)-log(OBS_landings_HL_E))^2)) #lognormal dist. so resids are ln(obs)-ln(pred) and so want SD of difference in log space, also sd function does not account for n in r so use RMSE
stan_res_landings_HL_E<-res_landings_HL_E/ln_RMSE_landings_HL_E

matplot(yrs_landings_HL_E,stan_res_landings_HL_E,type='n',xlab='Year',sub=bquote("(RMSE"==.(ln_RMSE_landings_HL_E)~")"),
        ylab='[ln(Obs)-ln(Pred)]/SD',xlim=c(min(yrs_landings_HL_E)-1,max(yrs_landings_HL_E)+1),
        ylim=c(min(stan_res_landings_HL_E)-.05*abs(min(stan_res_landings_HL_E)),max(stan_res_landings_HL_E)+.05*abs(max(stan_res_landings_HL_E)))) 
grid()
matlines(x=((min(yrs_landings_HL_E)-1):(max(yrs_landings_HL_E)+1)),y=rep(0,times=nyrs_landings_HL_E+2),col="blue",lwd=2);
matpoints(yrs_landings_HL_E,stan_res_landings_HL_E,col="red",pch=19);
lines(x=c(yrs_landings_HL_E,yrs_landings_HL_E),y=c(stan_res_landings_HL_E,rep(0,nyrs_landings_HL_E)),type='h');
title("Landings Standardized Residuals Handline East",line=1)



matplot(yrs_landings_HL_W,pred_landings_HL_W,type='n',main='Landings Handline West',
     xlab='Year',ylab='Landings (mt)',xlim=c(min(yrs_landings_HL_W),max(yrs_landings_HL_W)),
ylim=c(0,(max(pred_landings_HL_W,OBS_landings_HL_W)+.1*max(pred_landings_HL_W,OBS_landings_HL_W)))) 
grid()
matlines(yrs_landings_HL_W,pred_landings_HL_W,col='blue',lwd=2);
matpoints(yrs_landings_HL_W,OBS_landings_HL_W,col='red',pch=16);
legend("topleft",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)


ln_RMSE_landings_HL_W<-sqrt(mean((log(pred_landings_HL_W)-log(OBS_landings_HL_W))^2)) #lognormal dist. so resids are ln(obs)-ln(pred) and so want SD of difference in log space, also sd function does not account for n in r so use RMSE
stan_res_landings_HL_W<-res_landings_HL_W/ln_RMSE_landings_HL_W

matplot(yrs_landings_HL_W,stan_res_landings_HL_W,type='n',xlab='Year',sub=bquote("(RMSE"==.(ln_RMSE_landings_HL_W)~")"),
        ylab='[ln(Obs)-ln(Pred)]/SD',xlim=c(min(yrs_landings_HL_W)-1,max(yrs_landings_HL_W)+1),
        ylim=c(min(stan_res_landings_HL_W)-.05*abs(min(stan_res_landings_HL_W)),max(stan_res_landings_HL_W)+.05*abs(max(stan_res_landings_HL_W)))) 
grid()
matlines(x=((min(yrs_landings_HL_W)-1):(max(yrs_landings_HL_W)+1)),y=rep(0,times=nyrs_landings_HL_W+2),col="blue",lwd=2);
matpoints(yrs_landings_HL_W,stan_res_landings_HL_W,col="red",pch=19);
lines(x=c(yrs_landings_HL_W,yrs_landings_HL_W),y=c(stan_res_landings_HL_W,rep(0,nyrs_landings_HL_W)),type='h');
title("Landings Standardized Residuals Handline West",line=1)


matplot(yrs_landings_LL_E,pred_landings_LL_E,type='n',main='Landings Longline East',
     xlab='Year',ylab='Landings (mt)',xlim=c(min(yrs_landings_LL_E),max(yrs_landings_LL_E)),
ylim=c(0,(max(pred_landings_LL_E,OBS_landings_LL_E)+.1*max(pred_landings_LL_E,OBS_landings_LL_E)))) 
grid()
matlines(yrs_landings_LL_E,pred_landings_LL_E,col='blue',lwd=2);
matpoints(yrs_landings_LL_E,OBS_landings_LL_E,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)


ln_RMSE_landings_LL_E<-sqrt(mean((log(pred_landings_LL_E)-log(OBS_landings_LL_E))^2)) #lognormal dist. so resids are ln(obs)-ln(pred) and so want SD of difference in log space, also sd function does not account for n in r so use RMSE
stan_res_landings_LL_E<-res_landings_LL_E/ln_RMSE_landings_LL_E

matplot(yrs_landings_LL_E,stan_res_landings_LL_E,type='n',xlab='Year',sub=bquote("(RMSE"==.(ln_RMSE_landings_LL_E)~")"),
        ylab='[ln(Obs)-ln(Pred)]/SD',xlim=c(min(yrs_landings_LL_E)-1,max(yrs_landings_LL_E)+1),
        ylim=c(min(stan_res_landings_LL_E)-.05*abs(min(stan_res_landings_LL_E)),max(stan_res_landings_LL_E)+.05*abs(max(stan_res_landings_LL_E)))) 
grid()
matlines(x=((min(yrs_landings_LL_E)-1):(max(yrs_landings_LL_E)+1)),y=rep(0,times=nyrs_landings_LL_E+2),col="blue",lwd=2);
matpoints(yrs_landings_LL_E,stan_res_landings_LL_E,col="red",pch=19);
lines(x=c(yrs_landings_LL_E,yrs_landings_LL_E),y=c(stan_res_landings_LL_E,rep(0,nyrs_landings_LL_E)),type='h');
title("Landings Standardized Residuals Longline East",line=1)



matplot(yrs_landings_LL_W,pred_landings_LL_W,type='n',main='Landings Longline West',
     xlab='Year',ylab='Landings (mt)',xlim=c(min(yrs_landings_LL_W),max(yrs_landings_LL_W)),
ylim=c(0,(max(pred_landings_LL_W,OBS_landings_LL_W)+.1*max(pred_landings_LL_W,OBS_landings_LL_W)))) 
grid()
matlines(yrs_landings_LL_W,pred_landings_LL_W,col='blue',lwd=2);
matpoints(yrs_landings_LL_W,OBS_landings_LL_W,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)


ln_RMSE_landings_LL_W<-sqrt(mean((log(pred_landings_LL_W)-log(OBS_landings_LL_W))^2)) #lognormal dist. so resids are ln(obs)-ln(pred) and so want SD of difference in log space, also sd function does not account for n in r so use RMSE
stan_res_landings_LL_W<-res_landings_LL_W/ln_RMSE_landings_LL_W

matplot(yrs_landings_LL_W,stan_res_landings_LL_W,type='n',xlab='Year',sub=bquote("(RMSE"==.(ln_RMSE_landings_LL_W)~")"),
        ylab='[ln(Obs)-ln(Pred)]/SD',xlim=c(min(yrs_landings_LL_W)-1,max(yrs_landings_LL_W)+1),
        ylim=c(min(stan_res_landings_LL_W)-.05*abs(min(stan_res_landings_LL_W)),max(stan_res_landings_LL_W)+.05*abs(max(stan_res_landings_LL_W)))) 
grid()
matlines(x=((min(yrs_landings_LL_W)-1):(max(yrs_landings_LL_W)+1)),y=rep(0,times=nyrs_landings_LL_W+2),col="blue",lwd=2);
matpoints(yrs_landings_LL_W,stan_res_landings_LL_W,col="red",pch=19);
lines(x=c(yrs_landings_LL_W,yrs_landings_LL_W),y=c(stan_res_landings_LL_W,rep(0,nyrs_landings_LL_W)),type='h');
title("Landings Standardized Residuals Longline West",line=1)



matplot(yrs_landings_MRIP_E,pred_landings_MRIP_E,type='n',main='Landings MRIP East',
     xlab='Year',ylab='Landings (1000s of Fish)',xlim=c(min(yrs_landings_MRIP_E),max(yrs_landings_MRIP_E)),
ylim=c(0,(max(pred_landings_MRIP_E,OBS_landings_MRIP_E)+.1*max(pred_landings_MRIP_E,OBS_landings_MRIP_E)))) 
grid()
matlines(yrs_landings_MRIP_E,pred_landings_MRIP_E,col='blue',lwd=2);
matpoints(yrs_landings_MRIP_E,OBS_landings_MRIP_E,col='red',pch=16);
legend("topleft",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)


ln_RMSE_landings_MRIP_E<-sqrt(mean((log(pred_landings_MRIP_E)-log(OBS_landings_MRIP_E))^2)) #lognormal dist. so resids are ln(obs)-ln(pred) and so want SD of difference in log space, also sd function does not account for n in r so use RMSE
stan_res_landings_MRIP_E<-res_landings_MRIP_E/ln_RMSE_landings_MRIP_E

matplot(yrs_landings_MRIP_E,stan_res_landings_MRIP_E,type='n',xlab='Year',sub=bquote("(RMSE"==.(ln_RMSE_landings_MRIP_E)~")"),
        ylab='[ln(Obs)-ln(Pred)]/SD',xlim=c(min(yrs_landings_MRIP_E)-1,max(yrs_landings_MRIP_E)+1),
        ylim=c(min(stan_res_landings_MRIP_E)-.05*abs(min(stan_res_landings_MRIP_E)),max(stan_res_landings_MRIP_E)+.05*abs(max(stan_res_landings_MRIP_E)))) 
grid()
matlines(x=((min(yrs_landings_MRIP_E)-1):(max(yrs_landings_MRIP_E)+1)),y=rep(0,times=nyrs_landings_MRIP_E+2),col="blue",lwd=2);
matpoints(yrs_landings_MRIP_E,stan_res_landings_MRIP_E,col="red",pch=19);
lines(x=c(yrs_landings_MRIP_E,yrs_landings_MRIP_E),y=c(stan_res_landings_MRIP_E,rep(0,nyrs_landings_MRIP_E)),type='h');
title("Landings Standardized Residuals MRIP East",line=1)



matplot(yrs_landings_MRIP_W,pred_landings_MRIP_W,type='n',main='Landings MRIP West',
     xlab='Year',ylab='Landings (1000s of Fish)',xlim=c(min(yrs_landings_MRIP_W),max(yrs_landings_MRIP_W)),
ylim=c(0,(max(pred_landings_MRIP_W,OBS_landings_MRIP_W)+.1*max(pred_landings_MRIP_W,OBS_landings_MRIP_W)))) 
grid()
matlines(yrs_landings_MRIP_W,pred_landings_MRIP_W,col='blue',lwd=2);
matpoints(yrs_landings_MRIP_W,OBS_landings_MRIP_W,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)


ln_RMSE_landings_MRIP_W<-sqrt(mean((log(pred_landings_MRIP_W)-log(OBS_landings_MRIP_W))^2)) #lognormal dist. so resids are ln(obs)-ln(pred) and so want SD of difference in log space, also sd function does not account for n in r so use RMSE
stan_res_landings_MRIP_W<-res_landings_MRIP_W/ln_RMSE_landings_MRIP_W

matplot(yrs_landings_MRIP_W,stan_res_landings_MRIP_W,type='n',xlab='Year',sub=bquote("(RMSE"==.(ln_RMSE_landings_MRIP_W)~")"),
        ylab='[ln(Obs)-ln(Pred)]/SD',xlim=c(min(yrs_landings_MRIP_W)-1,max(yrs_landings_MRIP_W)+1),
        ylim=c(min(stan_res_landings_MRIP_W)-.05*abs(min(stan_res_landings_MRIP_W)),max(stan_res_landings_MRIP_W)+.05*abs(max(stan_res_landings_MRIP_W)))) 
grid()
matlines(x=((min(yrs_landings_MRIP_W)-1):(max(yrs_landings_MRIP_W)+1)),y=rep(0,times=nyrs_landings_MRIP_W+2),col="blue",lwd=2);
matpoints(yrs_landings_MRIP_W,stan_res_landings_MRIP_W,col="red",pch=19);
lines(x=c(yrs_landings_MRIP_W,yrs_landings_MRIP_W),y=c(stan_res_landings_MRIP_W,rep(0,nyrs_landings_MRIP_W)),type='h');
title("Landings Standardized Residuals MRIP West",line=1)



matplot(yrs_landings_HBT_E,pred_landings_HBT_E,type='n',main='Landings Headboat East',
     xlab='Year',ylab='Landings (1000s of Fish)',xlim=c(min(yrs_landings_HBT_E),max(yrs_landings_HBT_E)),
ylim=c(0,(max(pred_landings_HBT_E,OBS_landings_HBT_E)+.1*max(pred_landings_HBT_E,OBS_landings_HBT_E)))) 
grid()
matlines(yrs_landings_HBT_E,pred_landings_HBT_E,col='blue',lwd=2);
matpoints(yrs_landings_HBT_E,OBS_landings_HBT_E,col='red',pch=16);
legend("top",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)


ln_RMSE_landings_HBT_E<-sqrt(mean((log(pred_landings_HBT_E)-log(OBS_landings_HBT_E))^2)) #lognormal dist. so resids are ln(obs)-ln(pred) and so want SD of difference in log space, also sd function does not account for n in r so use RMSE
stan_res_landings_HBT_E<-res_landings_HBT_E/ln_RMSE_landings_HBT_E

matplot(yrs_landings_HBT_E,stan_res_landings_HBT_E,type='n',xlab='Year',sub=bquote("(RMSE"==.(ln_RMSE_landings_HBT_E)~")"),
        ylab='[ln(Obs)-ln(Pred)]/SD',xlim=c(min(yrs_landings_HBT_E)-1,max(yrs_landings_HBT_E)+1),
        ylim=c(min(stan_res_landings_HBT_E)-.05*abs(min(stan_res_landings_HBT_E)),max(stan_res_landings_HBT_E)+.05*abs(max(stan_res_landings_HBT_E)))) 
grid()
matlines(x=((min(yrs_landings_HBT_E)-1):(max(yrs_landings_HBT_E)+1)),y=rep(0,times=nyrs_landings_HBT_E+2),col="blue",lwd=2);
matpoints(yrs_landings_HBT_E,stan_res_landings_HBT_E,col="red",pch=19);
lines(x=c(yrs_landings_HBT_E,yrs_landings_HBT_E),y=c(stan_res_landings_HBT_E,rep(0,nyrs_landings_HBT_E)),type='h');
title("Landings Standardized Residuals Headboat East",line=1)



matplot(yrs_landings_HBT_W,pred_landings_HBT_W,type='n',main='Landings Headboat West',
     xlab='Year',ylab='Landings (1000s of Fish)',xlim=c(min(yrs_landings_HBT_W),max(yrs_landings_HBT_W)),
ylim=c(0,(max(pred_landings_HBT_W,OBS_landings_HBT_W)+.1*max(pred_landings_HBT_W,OBS_landings_HBT_W)))) 
grid()
matlines(yrs_landings_HBT_W,pred_landings_HBT_W,col='blue',lwd=2);
matpoints(yrs_landings_HBT_W,OBS_landings_HBT_W,col='red',pch=16);
legend("topright",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)


ln_RMSE_landings_HBT_W<-sqrt(mean((log(pred_landings_HBT_W)-log(OBS_landings_HBT_W))^2)) #lognormal dist. so resids are ln(obs)-ln(pred) and so want SD of difference in log space, also sd function does not account for n in r so use RMSE
stan_res_landings_HBT_W<-res_landings_HBT_W/ln_RMSE_landings_HBT_W

matplot(yrs_landings_HBT_W,stan_res_landings_HBT_W,type='n',xlab='Year',sub=bquote("(RMSE"==.(ln_RMSE_landings_HBT_W)~")"),
        ylab='[ln(Obs)-ln(Pred)]/SD',xlim=c(min(yrs_landings_HBT_W)-1,max(yrs_landings_HBT_W)+1),
        ylim=c(min(stan_res_landings_HBT_W)-.05*abs(min(stan_res_landings_HBT_W)),max(stan_res_landings_HBT_W)+.05*abs(max(stan_res_landings_HBT_W)))) 
grid()
matlines(x=((min(yrs_landings_HBT_W)-1):(max(yrs_landings_HBT_W)+1)),y=rep(0,times=nyrs_landings_HBT_W+2),col="blue",lwd=2);
matpoints(yrs_landings_HBT_W,stan_res_landings_HBT_W,col="red",pch=19);
lines(x=c(yrs_landings_HBT_W,yrs_landings_HBT_W),y=c(stan_res_landings_HBT_W,rep(0,nyrs_landings_HBT_W)),type='h');
title("Landings Standardized Residuals Headboat West",line=1)


matplot(yrs_landings_SHR_E,pred_landings_SHR_E,type='n',main='Bycatch Shrimp East',
     xlab='Year',ylab='Bycatch (1000s of Fish)',xlim=c(min(yrs_landings_SHR_E),max(yrs_landings_SHR_E)),
ylim=c(0,(max(pred_landings_SHR_E,OBS_landings_SHR_E)+.1*max(pred_landings_SHR_E,OBS_landings_SHR_E)))) 
grid()
matlines(yrs_landings_SHR_E,pred_landings_SHR_E,col='blue',lwd=2);
matpoints(yrs_landings_SHR_E,OBS_landings_SHR_E,col='red',pch=16);
legend("topright",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)

matplot(yrs_landings_SHR_W,pred_landings_SHR_W,type='n',main='Bycatch Shrimp West',
     xlab='Year',ylab='Bycatch (1000s of Fish)',xlim=c(min(yrs_landings_SHR_W),max(yrs_landings_SHR_W)),
ylim=c(0,(max(pred_landings_SHR_W,OBS_landings_SHR_W)+.1*max(pred_landings_SHR_W,OBS_landings_SHR_W)))) 
grid()
matlines(yrs_landings_SHR_W,pred_landings_SHR_W,col='blue',lwd=2);
matpoints(yrs_landings_SHR_W,OBS_landings_SHR_W,col='red',pch=16);
legend("topright",ncol=2, c("Pred","Obs"
               ), col = c('blue','red'),lty=c(1,-1),pch=c(-1,16),lwd=2,cex=.65)





###################### catch prop fit and resids



pred.catch.propHL_E<- as.data.frame(pred_catch_prop_HL_E)
pred.catch.propHL_E2<-stack(pred.catch.propHL_E)
pred.catch.propHL_E2$age<-rep(0:(nages-1),each=nyrs_age_comps_HL_E)
pred.catch.propHL_E2$year<-rep(yrs_age_comps_HL_E,times=nages)
pred.catch.propHL_E2$type<-rep("Predicted",times=nages*nyrs_age_comps_HL_E)

OBS.catch.propHL_E<- as.data.frame(OBS_catch_prop_HL_E)
OBS.catch.propHL_E2<-stack(OBS.catch.propHL_E)
OBS.catch.propHL_E2$age<-rep(0:(nages-1),each=nyrs_age_comps_HL_E)
OBS.catch.propHL_E2$year<-rep(yrs_age_comps_HL_E,times=nages)
OBS.catch.propHL_E2$type<-rep("Observed",times=nages*nyrs_age_comps_HL_E)


catch.propHL_E<-as.data.frame(cbind(OBS.catch.propHL_E2,pred.catch.propHL_E2$values,
              pred.catch.propHL_E2$type))

print(xyplot(pred.catch.propHL_E2$values+values~age|paste(year),data=catch.propHL_E, 
      scales=list(alternating=c(1,2)),
      main='Handline East Catch Age Proportions',xlab='Ages', ylab='Proportion at Age',
      outer=FALSE,distribute.type = TRUE, col=c("blue","red"),type=c('l','o'),pch=c(0,16),lty=c(1,0),lwd=c(2,0),
      auto.key = list ("topright",border=TRUE,padding.text=2.5,columns =2,text=c("Pred","Obs")), 
      par.settings = simpleTheme (pch=c(15,16), lwd=c(8,0),lty=c(1,0),col=c("blue","red")))
)




RMSE_catch_prop_HL_E<-sqrt(mean((pred_catch_prop_HL_E-OBS_catch_prop_HL_E)^2))
stan_res_catch_prop_HL_E<-(OBS_catch_prop_HL_E-pred_catch_prop_HL_E)/RMSE_catch_prop_HL_E

layout(matrix(c(1), 1, byrow = TRUE))
years <- (yrs_age_comps_HL_E[1]-2):(yrs_age_comps_HL_E[nyrs_age_comps_HL_E])
nyears <- length(years)
nindices <- nages
res.catch.prop_HL_E <- stan_res_catch_prop_HL_E #matrix(rnorm(nyears*nindices), ncol=nindices, nrow=nyears) 
resid.pos.color <- "blue"
resid.neg.color <- "red"
resid.col=matrix(NA, nrow=nyears, ncol=nindices)   # set color for residual bubbles
resid.col <- ifelse(res.catch.prop_HL_E > 0.0,resid.pos.color, resid.neg.color) 
# use a scaler to increase (>1) or decrease (0<val<1) the size of the bubbles to make them look good
res.catch.prop_HL_E1 <- res.catch.prop_HL_E
# first set up a blank space for the plot
plot(seq(0,(nindices-1)), rep(years[1],nindices), ylim=c(years[nyears],years[1]), xlab = "Age",
     ylab = "Year",  type = "n", axes=F,main='Handline East Catch-at-Age Proportion Standardized\nResiduals [(OBS-PRED)/SD]',
sub=bquote("(Absolute Value Max Residual"==.(max(abs(stan_res_catch_prop_HL_E)))~", RMSE"==.(RMSE_catch_prop_HL_E)~")"))


axis(1, at= seq(0,(nindices-1)), lab=seq(0,(nindices-1)))
axis(2, at = seq(years[nyears],years[1], by=-3), lab = seq(years[nyears],years[1], by=-3),las=1)
box()
abline(h=seq( years[3], years[nyears]),col = "lightgray", lty = 1)
abline(a=1965,b=1,col = "lightgray", lty = 1)
abline(a=1970,b=1,col = "lightgray", lty = 1)
abline(a=1975,b=1,col = "lightgray", lty = 1)
abline(a=1980,b=1,col = "lightgray", lty = 1)
abline(a=1985,b=1,col = "lightgray", lty = 1)
abline(a=1990,b=1,col = "lightgray", lty = 1)
abline(a=1995,b=1,col = "lightgray", lty = 1)
abline(a=2000,b=1,col = "lightgray", lty = 1)
abline(a=2005,b=1,col = "lightgray", lty = 1)
abline(a=2010,b=1,col = "lightgray", lty = 1)

# fill in the bubbles using points command, cex determines size, bg determines color
for (i in 3:nyears){
  points(seq(0, (nindices-1)), rep((years[i]), nindices),  
  cex =abs(res.catch.prop_HL_E1[i-2,]) ,  pch = 21, bg = resid.col[i-2,], col='black')
}
legend("top",ncol=2, c(" Negative", " Positive"),pch=c(16,16),
col = c("red","blue")
,cex=.8, merge=FALSE)


RMSE_catch_prop_HL_E<-sqrt(mean((pred_catch_prop_HL_E-OBS_catch_prop_HL_E)^2))
stan_res_catch_prop_HL_E<-(OBS_catch_prop_HL_E-pred_catch_prop_HL_E)/RMSE_catch_prop_HL_E

res.catch.prop_HL_E<- as.data.frame(stan_res_catch_prop_HL_E)
res.catch.prop_HL_E2<-stack(res.catch.prop_HL_E)
res.catch.prop_HL_E2$age<-rep(0:(nages-1),each=nyrs_age_comps_HL_E)
res.catch.prop_HL_E2$year<-rep(yrs_age_comps_HL_E,times=nages)

print(xyplot(values~age|paste(year),
      data=res.catch.prop_HL_E2,scales=list(alternating=c(1,2)),
      main='Handline East Catch-at-Age Proportion \nStandardized Residuals',xlab='Ages', ylab='Residuals',
      outer=FALSE,distribute.type = TRUE, type=c("b",'g'),lwd=1,pch=16,col='blue2')
)




pred.catch.propHL_W<- as.data.frame(pred_catch_prop_HL_W)
pred.catch.propHL_W2<-stack(pred.catch.propHL_W)
pred.catch.propHL_W2$age<-rep(0:(nages-1),each=nyrs_age_comps_HL_W)
pred.catch.propHL_W2$year<-rep(yrs_age_comps_HL_W,times=nages)
pred.catch.propHL_W2$type<-rep("Predicted",times=nages*nyrs_age_comps_HL_W)

OBS.catch.propHL_W<- as.data.frame(OBS_catch_prop_HL_W)
OBS.catch.propHL_W2<-stack(OBS.catch.propHL_W)
OBS.catch.propHL_W2$age<-rep(0:(nages-1),each=nyrs_age_comps_HL_W)
OBS.catch.propHL_W2$year<-rep(yrs_age_comps_HL_W,times=nages)
OBS.catch.propHL_W2$type<-rep("Observed",times=nages*nyrs_age_comps_HL_W)


catch.propHL_W<-as.data.frame(cbind(OBS.catch.propHL_W2,pred.catch.propHL_W2$values,
              pred.catch.propHL_W2$type))

print(xyplot(pred.catch.propHL_W2$values+values~age|paste(year),data=catch.propHL_W, 
      scales=list(alternating=c(1,2)),
      main='Handline West Catch Age Proportions',xlab='Ages', ylab='Proportion at Age',
      outer=FALSE,distribute.type = TRUE, col=c("blue","red"),type=c('l','o'),pch=c(0,16),lty=c(1,0),lwd=c(2,0),
      auto.key = list ("topright",border=TRUE,padding.text=2.5,columns =2,text=c("Pred","Obs")), 
      par.settings = simpleTheme (pch=c(15,16), lwd=c(8,0),lty=c(1,0),col=c("blue","red")))
)


RMSE_catch_prop_HL_W<-sqrt(mean((pred_catch_prop_HL_W-OBS_catch_prop_HL_W)^2))
stan_res_catch_prop_HL_W<-(OBS_catch_prop_HL_W-pred_catch_prop_HL_W)/RMSE_catch_prop_HL_W

layout(matrix(c(1), 1, byrow = TRUE))
years <- (yrs_age_comps_HL_W[1]-2):(yrs_age_comps_HL_W[nyrs_age_comps_HL_W])
nyears <- length(years)
nindices <- nages
res.catch.prop_HL_W <- stan_res_catch_prop_HL_W #matrix(rnorm(nyears*nindices), ncol=nindices, nrow=nyears) 
resid.pos.color <- "blue"
resid.neg.color <- "red"
resid.col=matrix(NA, nrow=nyears, ncol=nindices)   # set color for residual bubbles
resid.col <- ifelse(res.catch.prop_HL_W > 0.0,resid.pos.color, resid.neg.color) 
# use a scaler to increase (>1) or decrease (0<val<1) the size of the bubbles to make them look good
res.catch.prop_HL_W1 <- res.catch.prop_HL_W
# first set up a blank space for the plot
plot(seq(0,(nindices-1)), rep(years[1],nindices), ylim=c(years[nyears],years[1]), xlab = "Age",
     ylab = "Year",  type = "n", axes=F,main='Handline West Catch-at-Age Proportion Standardized\nResiduals [(OBS-PRED)/SD]',
sub=bquote("(Absolute Value Max Residual"==.(max(abs(stan_res_catch_prop_HL_W)))~", RMSE"==.(RMSE_catch_prop_HL_W)~")"))


axis(1, at= seq(0,(nindices-1)), lab=seq(0,(nindices-1)))
axis(2, at = seq(years[nyears],years[1], by=-3), lab = seq(years[nyears],years[1], by=-3),las=1)
box()
abline(h=seq( years[3], years[nyears]),col = "lightgray", lty = 1)
abline(a=1965,b=1,col = "lightgray", lty = 1)
abline(a=1970,b=1,col = "lightgray", lty = 1)
abline(a=1975,b=1,col = "lightgray", lty = 1)
abline(a=1980,b=1,col = "lightgray", lty = 1)
abline(a=1985,b=1,col = "lightgray", lty = 1)
abline(a=1990,b=1,col = "lightgray", lty = 1)
abline(a=1995,b=1,col = "lightgray", lty = 1)
abline(a=2000,b=1,col = "lightgray", lty = 1)
abline(a=2005,b=1,col = "lightgray", lty = 1)
abline(a=2010,b=1,col = "lightgray", lty = 1)

# fill in the bubbles using points command, cex determines size, bg determines color
for (i in 3:nyears){
  points(seq(0, (nindices-1)), rep((years[i]), nindices),  
  cex =abs(res.catch.prop_HL_W1[i-2,]) ,  pch = 21, bg = resid.col[i-2,], col='black')
}
legend("top",ncol=2, c(" Negative", " Positive"),pch=c(16,16),
col = c("red","blue")
,cex=.8, merge=FALSE)


RMSE_catch_prop_HL_W<-sqrt(mean((pred_catch_prop_HL_W-OBS_catch_prop_HL_W)^2))
stan_res_catch_prop_HL_W<-(OBS_catch_prop_HL_W-pred_catch_prop_HL_W)/RMSE_catch_prop_HL_W

res.catch.prop_HL_W<- as.data.frame(stan_res_catch_prop_HL_W)
res.catch.prop_HL_W2<-stack(res.catch.prop_HL_W)
res.catch.prop_HL_W2$age<-rep(0:(nages-1),each=nyrs_age_comps_HL_W)
res.catch.prop_HL_W2$year<-rep(yrs_age_comps_HL_W,times=nages)

print(xyplot(values~age|paste(year),
      data=res.catch.prop_HL_W2,scales=list(alternating=c(1,2)),
      main='Handline West Catch-at-Age Proportion \nStandardized Residuals',xlab='Ages', ylab='Residuals',
      outer=FALSE,distribute.type = TRUE, type=c("b",'g'),lwd=1,pch=16,col='blue2')
)







pred.catch.propLL_E<- as.data.frame(pred_catch_prop_LL_E)
pred.catch.propLL_E2<-stack(pred.catch.propLL_E)
pred.catch.propLL_E2$age<-rep(0:(nages-1),each=nyrs_age_comps_LL_E)
pred.catch.propLL_E2$year<-rep(yrs_age_comps_LL_E,times=nages)
pred.catch.propLL_E2$type<-rep("Predicted",times=nages*nyrs_age_comps_LL_E)

OBS.catch.propLL_E<- as.data.frame(OBS_catch_prop_LL_E)
OBS.catch.propLL_E2<-stack(OBS.catch.propLL_E)
OBS.catch.propLL_E2$age<-rep(0:(nages-1),each=nyrs_age_comps_LL_E)
OBS.catch.propLL_E2$year<-rep(yrs_age_comps_LL_E,times=nages)
OBS.catch.propLL_E2$type<-rep("Observed",times=nages*nyrs_age_comps_LL_E)


catch.propLL_E<-as.data.frame(cbind(OBS.catch.propLL_E2,pred.catch.propLL_E2$values,
              pred.catch.propLL_E2$type))

print(xyplot(pred.catch.propLL_E2$values+values~age|paste(year),data=catch.propLL_E, 
      scales=list(alternating=c(1,2)),
      main='Longline East Catch Age Proportions',xlab='Ages', ylab='Proportion at Age',
      outer=FALSE,distribute.type = TRUE, col=c("blue","red"),type=c('l','o'),pch=c(0,16),lty=c(1,0),lwd=c(2,0),
      auto.key = list ("topright",border=TRUE,padding.text=2.5,columns =2,text=c("Pred","Obs")), 
      par.settings = simpleTheme (pch=c(15,16), lwd=c(8,0),lty=c(1,0),col=c("blue","red")))
)




RMSE_catch_prop_LL_E<-sqrt(mean((pred_catch_prop_LL_E-OBS_catch_prop_LL_E)^2))
stan_res_catch_prop_LL_E<-(OBS_catch_prop_LL_E-pred_catch_prop_LL_E)/RMSE_catch_prop_LL_E

layout(matrix(c(1), 1, byrow = TRUE))
years <- (yrs_age_comps_LL_E[1]-2):(yrs_age_comps_LL_E[nyrs_age_comps_LL_E])
nyears <- length(years)
nindices <- nages
res.catch.prop_LL_E <- stan_res_catch_prop_LL_E #matrix(rnorm(nyears*nindices), ncol=nindices, nrow=nyears) 
resid.pos.color <- "blue"
resid.neg.color <- "red"
resid.col=matrix(NA, nrow=nyears, ncol=nindices)   # set color for residual bubbles
resid.col <- ifelse(res.catch.prop_LL_E > 0.0,resid.pos.color, resid.neg.color) 
# use a scaler to increase (>1) or decrease (0<val<1) the size of the bubbles to make them look good
res.catch.prop_LL_E1 <- res.catch.prop_LL_E
# first set up a blank space for the plot
plot(seq(0,(nindices-1)), rep(years[1],nindices), ylim=c(years[nyears],years[1]), xlab = "Age",
     ylab = "Year",  type = "n", axes=F,main='Longline East Catch-at-Age Proportion Standardized\nResiduals [(OBS-PRED)/SD]',
sub=bquote("(Absolute Value Max Residual"==.(max(abs(stan_res_catch_prop_LL_E)))~", RMSE"==.(RMSE_catch_prop_LL_E)~")"))


axis(1, at= seq(0,(nindices-1)), lab=seq(0,(nindices-1)))
axis(2, at = seq(years[nyears],years[1], by=-3), lab = seq(years[nyears],years[1], by=-3),las=1)
box()
abline(h=seq( years[3], years[nyears]),col = "lightgray", lty = 1)
abline(a=1965,b=1,col = "lightgray", lty = 1)
abline(a=1970,b=1,col = "lightgray", lty = 1)
abline(a=1975,b=1,col = "lightgray", lty = 1)
abline(a=1980,b=1,col = "lightgray", lty = 1)
abline(a=1985,b=1,col = "lightgray", lty = 1)
abline(a=1990,b=1,col = "lightgray", lty = 1)
abline(a=1995,b=1,col = "lightgray", lty = 1)
abline(a=2000,b=1,col = "lightgray", lty = 1)
abline(a=2005,b=1,col = "lightgray", lty = 1)
abline(a=2010,b=1,col = "lightgray", lty = 1)

# fill in the bubbles using points command, cex determines size, bg determines color
for (i in 3:nyears){
  points(seq(0, (nindices-1)), rep((years[i]), nindices),  
  cex =abs(res.catch.prop_LL_E1[i-2,]) ,  pch = 21, bg = resid.col[i-2,], col='black')
}
legend("top",ncol=2, c(" Negative", " Positive"),pch=c(16,16),
col = c("red","blue")
,cex=.8, merge=FALSE)


RMSE_catch_prop_LL_E<-sqrt(mean((pred_catch_prop_LL_E-OBS_catch_prop_LL_E)^2))
stan_res_catch_prop_LL_E<-(OBS_catch_prop_LL_E-pred_catch_prop_LL_E)/RMSE_catch_prop_LL_E

res.catch.prop_LL_E<- as.data.frame(stan_res_catch_prop_LL_E)
res.catch.prop_LL_E2<-stack(res.catch.prop_LL_E)
res.catch.prop_LL_E2$age<-rep(0:(nages-1),each=nyrs_age_comps_LL_E)
res.catch.prop_LL_E2$year<-rep(yrs_age_comps_LL_E,times=nages)

print(xyplot(values~age|paste(year),
      data=res.catch.prop_LL_E2,scales=list(alternating=c(1,2)),
      main='Longline East Catch-at-Age Proportion \nStandardized Residuals',xlab='Ages', ylab='Residuals',
      outer=FALSE,distribute.type = TRUE, type=c("b",'g'),lwd=1,pch=16,col='blue2')
)




pred.catch.propLL_W<- as.data.frame(pred_catch_prop_LL_W)
pred.catch.propLL_W2<-stack(pred.catch.propLL_W)
pred.catch.propLL_W2$age<-rep(0:(nages-1),each=nyrs_age_comps_LL_W)
pred.catch.propLL_W2$year<-rep(yrs_age_comps_LL_W,times=nages)
pred.catch.propLL_W2$type<-rep("Predicted",times=nages*nyrs_age_comps_LL_W)

OBS.catch.propLL_W<- as.data.frame(OBS_catch_prop_LL_W)
OBS.catch.propLL_W2<-stack(OBS.catch.propLL_W)
OBS.catch.propLL_W2$age<-rep(0:(nages-1),each=nyrs_age_comps_LL_W)
OBS.catch.propLL_W2$year<-rep(yrs_age_comps_LL_W,times=nages)
OBS.catch.propLL_W2$type<-rep("Observed",times=nages*nyrs_age_comps_LL_W)


catch.propLL_W<-as.data.frame(cbind(OBS.catch.propLL_W2,pred.catch.propLL_W2$values,
              pred.catch.propLL_W2$type))

print(xyplot(pred.catch.propLL_W2$values+values~age|paste(year),data=catch.propLL_W, 
      scales=list(alternating=c(1,2)),
      main='Longline West Catch Age Proportions',xlab='Ages', ylab='Proportion at Age',
      outer=FALSE,distribute.type = TRUE, col=c("blue","red"),type=c('l','o'),pch=c(0,16),lty=c(1,0),lwd=c(2,0),
      auto.key = list ("topright",border=TRUE,padding.text=2.5,columns =2,text=c("Pred","Obs")), 
      par.settings = simpleTheme (pch=c(15,16), lwd=c(8,0),lty=c(1,0),col=c("blue","red")))
)


RMSE_catch_prop_LL_W<-sqrt(mean((pred_catch_prop_LL_W-OBS_catch_prop_LL_W)^2))
stan_res_catch_prop_LL_W<-(OBS_catch_prop_LL_W-pred_catch_prop_LL_W)/RMSE_catch_prop_LL_W

layout(matrix(c(1), 1, byrow = TRUE))
years <- (yrs_age_comps_LL_W[1]-2):(yrs_age_comps_LL_W[nyrs_age_comps_LL_W])
nyears <- length(years)
nindices <- nages
res.catch.prop_LL_W <- stan_res_catch_prop_LL_W #matrix(rnorm(nyears*nindices), ncol=nindices, nrow=nyears) 
resid.pos.color <- "blue"
resid.neg.color <- "red"
resid.col=matrix(NA, nrow=nyears, ncol=nindices)   # set color for residual bubbles
resid.col <- ifelse(res.catch.prop_LL_W > 0.0,resid.pos.color, resid.neg.color) 
# use a scaler to increase (>1) or decrease (0<val<1) the size of the bubbles to make them look good
res.catch.prop_LL_W1 <- res.catch.prop_LL_W
# first set up a blank space for the plot
plot(seq(0,(nindices-1)), rep(years[1],nindices), ylim=c(years[nyears],years[1]), xlab = "Age",
     ylab = "Year",  type = "n", axes=F,main='Longline West Catch-at-Age Proportion Standardized\nResiduals [(OBS-PRED)/SD]',
sub=bquote("(Absolute Value Max Residual"==.(max(abs(stan_res_catch_prop_LL_W)))~", RMSE"==.(RMSE_catch_prop_LL_W)~")"))


axis(1, at= seq(0,(nindices-1)), lab=seq(0,(nindices-1)))
axis(2, at = seq(years[nyears],years[1], by=-3), lab = seq(years[nyears],years[1], by=-3),las=1)
box()
abline(h=seq( years[3], years[nyears]),col = "lightgray", lty = 1)
abline(a=1965,b=1,col = "lightgray", lty = 1)
abline(a=1970,b=1,col = "lightgray", lty = 1)
abline(a=1975,b=1,col = "lightgray", lty = 1)
abline(a=1980,b=1,col = "lightgray", lty = 1)
abline(a=1985,b=1,col = "lightgray", lty = 1)
abline(a=1990,b=1,col = "lightgray", lty = 1)
abline(a=1995,b=1,col = "lightgray", lty = 1)
abline(a=2000,b=1,col = "lightgray", lty = 1)
abline(a=2005,b=1,col = "lightgray", lty = 1)
abline(a=2010,b=1,col = "lightgray", lty = 1)

# fill in the bubbles using points command, cex determines size, bg determines color
for (i in 3:nyears){
  points(seq(0, (nindices-1)), rep((years[i]), nindices),  
  cex =abs(res.catch.prop_LL_W1[i-2,]) ,  pch = 21, bg = resid.col[i-2,], col='black')
}
legend("top",ncol=2, c(" Negative", " Positive"),pch=c(16,16),
col = c("red","blue")
,cex=.8, merge=FALSE)


RMSE_catch_prop_LL_W<-sqrt(mean((pred_catch_prop_LL_W-OBS_catch_prop_LL_W)^2))
stan_res_catch_prop_LL_W<-(OBS_catch_prop_LL_W-pred_catch_prop_LL_W)/RMSE_catch_prop_LL_W

res.catch.prop_LL_W<- as.data.frame(stan_res_catch_prop_LL_W)
res.catch.prop_LL_W2<-stack(res.catch.prop_LL_W)
res.catch.prop_LL_W2$age<-rep(0:(nages-1),each=nyrs_age_comps_LL_W)
res.catch.prop_LL_W2$year<-rep(yrs_age_comps_LL_W,times=nages)

print(xyplot(values~age|paste(year),
      data=res.catch.prop_LL_W2,scales=list(alternating=c(1,2)),
      main='Longline West Catch-at-Age Proportion \nStandardized Residuals',xlab='Ages', ylab='Residuals',
      outer=FALSE,distribute.type = TRUE, type=c("b",'g'),lwd=1,pch=16,col='blue2')
)





pred.catch.propMRIP_E<- as.data.frame(pred_catch_prop_MRIP_E)
pred.catch.propMRIP_E2<-stack(pred.catch.propMRIP_E)
pred.catch.propMRIP_E2$age<-rep(0:(nages-1),each=nyrs_age_comps_MRIP_E)
pred.catch.propMRIP_E2$year<-rep(yrs_age_comps_MRIP_E,times=nages)
pred.catch.propMRIP_E2$type<-rep("Predicted",times=nages*nyrs_age_comps_MRIP_E)

OBS.catch.propMRIP_E<- as.data.frame(OBS_catch_prop_MRIP_E)
OBS.catch.propMRIP_E2<-stack(OBS.catch.propMRIP_E)
OBS.catch.propMRIP_E2$age<-rep(0:(nages-1),each=nyrs_age_comps_MRIP_E)
OBS.catch.propMRIP_E2$year<-rep(yrs_age_comps_MRIP_E,times=nages)
OBS.catch.propMRIP_E2$type<-rep("Observed",times=nages*nyrs_age_comps_MRIP_E)


catch.propMRIP_E<-as.data.frame(cbind(OBS.catch.propMRIP_E2,pred.catch.propMRIP_E2$values,
              pred.catch.propMRIP_E2$type))

print(xyplot(pred.catch.propMRIP_E2$values+values~age|paste(year),data=catch.propMRIP_E, 
      scales=list(alternating=c(1,2)),
      main='MRIP East Catch Age Proportions',xlab='Ages', ylab='Proportion at Age',
      outer=FALSE,distribute.type = TRUE, col=c("blue","red"),type=c('l','o'),pch=c(0,16),lty=c(1,0),lwd=c(2,0),
      auto.key = list ("topright",border=TRUE,padding.text=2.5,columns =2,text=c("Pred","Obs")), 
      par.settings = simpleTheme (pch=c(15,16), lwd=c(8,0),lty=c(1,0),col=c("blue","red")))
)




RMSE_catch_prop_MRIP_E<-sqrt(mean((pred_catch_prop_MRIP_E-OBS_catch_prop_MRIP_E)^2))
stan_res_catch_prop_MRIP_E<-(OBS_catch_prop_MRIP_E-pred_catch_prop_MRIP_E)/RMSE_catch_prop_MRIP_E

layout(matrix(c(1), 1, byrow = TRUE))
years <- (yrs_age_comps_MRIP_E[1]-2):(yrs_age_comps_MRIP_E[nyrs_age_comps_MRIP_E])
nyears <- length(years)
nindices <- nages
res.catch.prop_MRIP_E <- stan_res_catch_prop_MRIP_E #matrix(rnorm(nyears*nindices), ncol=nindices, nrow=nyears) 
resid.pos.color <- "blue"
resid.neg.color <- "red"
resid.col=matrix(NA, nrow=nyears, ncol=nindices)   # set color for residual bubbles
resid.col <- ifelse(res.catch.prop_MRIP_E > 0.0,resid.pos.color, resid.neg.color) 
# use a scaler to increase (>1) or decrease (0<val<1) the size of the bubbles to make them look good
res.catch.prop_MRIP_E1 <- res.catch.prop_MRIP_E
# first set up a blank space for the plot
plot(seq(0,(nindices-1)), rep(years[1],nindices), ylim=c(years[nyears],years[1]), xlab = "Age",
     ylab = "Year",  type = "n", axes=F,main='MRIP East Catch-at-Age Proportion Standardized\nResiduals [(OBS-PRED)/SD]',
sub=bquote("(Absolute Value Max Residual"==.(max(abs(stan_res_catch_prop_MRIP_E)))~", RMSE"==.(RMSE_catch_prop_MRIP_E)~")"))


axis(1, at= seq(0,(nindices-1)), lab=seq(0,(nindices-1)))
axis(2, at = seq(years[nyears],years[1], by=-3), lab = seq(years[nyears],years[1], by=-3),las=1)
box()
abline(h=seq( years[3], years[nyears]),col = "lightgray", lty = 1)
abline(a=1965,b=1,col = "lightgray", lty = 1)
abline(a=1970,b=1,col = "lightgray", lty = 1)
abline(a=1975,b=1,col = "lightgray", lty = 1)
abline(a=1980,b=1,col = "lightgray", lty = 1)
abline(a=1985,b=1,col = "lightgray", lty = 1)
abline(a=1990,b=1,col = "lightgray", lty = 1)
abline(a=1995,b=1,col = "lightgray", lty = 1)
abline(a=2000,b=1,col = "lightgray", lty = 1)
abline(a=2005,b=1,col = "lightgray", lty = 1)
abline(a=2010,b=1,col = "lightgray", lty = 1)

# fill in the bubbles using points command, cex determines size, bg determines color
for (i in 3:nyears){
  points(seq(0, (nindices-1)), rep((years[i]), nindices),  
  cex =abs(res.catch.prop_MRIP_E1[i-2,]) ,  pch = 21, bg = resid.col[i-2,], col='black')
}
legend("top",ncol=2, c(" Negative", " Positive"),pch=c(16,16),
col = c("red","blue")
,cex=.8, merge=FALSE)


RMSE_catch_prop_MRIP_E<-sqrt(mean((pred_catch_prop_MRIP_E-OBS_catch_prop_MRIP_E)^2))
stan_res_catch_prop_MRIP_E<-(OBS_catch_prop_MRIP_E-pred_catch_prop_MRIP_E)/RMSE_catch_prop_MRIP_E

res.catch.prop_MRIP_E<- as.data.frame(stan_res_catch_prop_MRIP_E)
res.catch.prop_MRIP_E2<-stack(res.catch.prop_MRIP_E)
res.catch.prop_MRIP_E2$age<-rep(0:(nages-1),each=nyrs_age_comps_MRIP_E)
res.catch.prop_MRIP_E2$year<-rep(yrs_age_comps_MRIP_E,times=nages)

print(xyplot(values~age|paste(year),
      data=res.catch.prop_MRIP_E2,scales=list(alternating=c(1,2)),
      main='MRIP East Catch-at-Age Proportion \nStandardized Residuals',xlab='Ages', ylab='Residuals',
      outer=FALSE,distribute.type = TRUE, type=c("b",'g'),lwd=1,pch=16,col='blue2')
)




pred.catch.propMRIP_W<- as.data.frame(pred_catch_prop_MRIP_W)
pred.catch.propMRIP_W2<-stack(pred.catch.propMRIP_W)
pred.catch.propMRIP_W2$age<-rep(0:(nages-1),each=nyrs_age_comps_MRIP_W)
pred.catch.propMRIP_W2$year<-rep(yrs_age_comps_MRIP_W,times=nages)
pred.catch.propMRIP_W2$type<-rep("Predicted",times=nages*nyrs_age_comps_MRIP_W)

OBS.catch.propMRIP_W<- as.data.frame(OBS_catch_prop_MRIP_W)
OBS.catch.propMRIP_W2<-stack(OBS.catch.propMRIP_W)
OBS.catch.propMRIP_W2$age<-rep(0:(nages-1),each=nyrs_age_comps_MRIP_W)
OBS.catch.propMRIP_W2$year<-rep(yrs_age_comps_MRIP_W,times=nages)
OBS.catch.propMRIP_W2$type<-rep("Observed",times=nages*nyrs_age_comps_MRIP_W)


catch.propMRIP_W<-as.data.frame(cbind(OBS.catch.propMRIP_W2,pred.catch.propMRIP_W2$values,
              pred.catch.propMRIP_W2$type))

print(xyplot(pred.catch.propMRIP_W2$values+values~age|paste(year),data=catch.propMRIP_W, 
      scales=list(alternating=c(1,2)),
      main='MRIP West Catch Age Proportions',xlab='Ages', ylab='Proportion at Age',
      outer=FALSE,distribute.type = TRUE, col=c("blue","red"),type=c('l','o'),pch=c(0,16),lty=c(1,0),lwd=c(2,0),
      auto.key = list ("topright",border=TRUE,padding.text=2.5,columns =2,text=c("Pred","Obs")), 
      par.settings = simpleTheme (pch=c(15,16), lwd=c(8,0),lty=c(1,0),col=c("blue","red")))
)


RMSE_catch_prop_MRIP_W<-sqrt(mean((pred_catch_prop_MRIP_W-OBS_catch_prop_MRIP_W)^2))
stan_res_catch_prop_MRIP_W<-(OBS_catch_prop_MRIP_W-pred_catch_prop_MRIP_W)/RMSE_catch_prop_MRIP_W

layout(matrix(c(1), 1, byrow = TRUE))
years <- (yrs_age_comps_MRIP_W[1]-2):(yrs_age_comps_MRIP_W[nyrs_age_comps_MRIP_W])
nyears <- length(years)
nindices <- nages
res.catch.prop_MRIP_W <- stan_res_catch_prop_MRIP_W #matrix(rnorm(nyears*nindices), ncol=nindices, nrow=nyears) 
resid.pos.color <- "blue"
resid.neg.color <- "red"
resid.col=matrix(NA, nrow=nyears, ncol=nindices)   # set color for residual bubbles
resid.col <- ifelse(res.catch.prop_MRIP_W > 0.0,resid.pos.color, resid.neg.color) 
# use a scaler to increase (>1) or decrease (0<val<1) the size of the bubbles to make them look good
res.catch.prop_MRIP_W1 <- res.catch.prop_MRIP_W
# first set up a blank space for the plot
plot(seq(0,(nindices-1)), rep(years[1],nindices), ylim=c(years[nyears],years[1]), xlab = "Age",
     ylab = "Year",  type = "n", axes=F,main='MRIP West Catch-at-Age Proportion Standardized\nResiduals [(OBS-PRED)/SD]',
sub=bquote("(Absolute Value Max Residual"==.(max(abs(stan_res_catch_prop_MRIP_W)))~", RMSE"==.(RMSE_catch_prop_MRIP_W)~")"))


axis(1, at= seq(0,(nindices-1)), lab=seq(0,(nindices-1)))
axis(2, at = seq(years[nyears],years[1], by=-3), lab = seq(years[nyears],years[1], by=-3),las=1)
box()
abline(h=seq( years[3], years[nyears]),col = "lightgray", lty = 1)
abline(a=1965,b=1,col = "lightgray", lty = 1)
abline(a=1970,b=1,col = "lightgray", lty = 1)
abline(a=1975,b=1,col = "lightgray", lty = 1)
abline(a=1980,b=1,col = "lightgray", lty = 1)
abline(a=1985,b=1,col = "lightgray", lty = 1)
abline(a=1990,b=1,col = "lightgray", lty = 1)
abline(a=1995,b=1,col = "lightgray", lty = 1)
abline(a=2000,b=1,col = "lightgray", lty = 1)
abline(a=2005,b=1,col = "lightgray", lty = 1)
abline(a=2010,b=1,col = "lightgray", lty = 1)

# fill in the bubbles using points command, cex determines size, bg determines color
for (i in 3:nyears){
  points(seq(0, (nindices-1)), rep((years[i]), nindices),  
  cex =abs(res.catch.prop_MRIP_W1[i-2,]) ,  pch = 21, bg = resid.col[i-2,], col='black')
}
legend("top",ncol=2, c(" Negative", " Positive"),pch=c(16,16),
col = c("red","blue")
,cex=.8, merge=FALSE)


RMSE_catch_prop_MRIP_W<-sqrt(mean((pred_catch_prop_MRIP_W-OBS_catch_prop_MRIP_W)^2))
stan_res_catch_prop_MRIP_W<-(OBS_catch_prop_MRIP_W-pred_catch_prop_MRIP_W)/RMSE_catch_prop_MRIP_W

res.catch.prop_MRIP_W<- as.data.frame(stan_res_catch_prop_MRIP_W)
res.catch.prop_MRIP_W2<-stack(res.catch.prop_MRIP_W)
res.catch.prop_MRIP_W2$age<-rep(0:(nages-1),each=nyrs_age_comps_MRIP_W)
res.catch.prop_MRIP_W2$year<-rep(yrs_age_comps_MRIP_W,times=nages)

print(xyplot(values~age|paste(year),
      data=res.catch.prop_MRIP_W2,scales=list(alternating=c(1,2)),
      main='MRIP West Catch-at-Age Proportion \nStandardized Residuals',xlab='Ages', ylab='Residuals',
      outer=FALSE,distribute.type = TRUE, type=c("b",'g'),lwd=1,pch=16,col='blue2')
)



pred.catch.propHBT_E<- as.data.frame(pred_catch_prop_HBT_E)
pred.catch.propHBT_E2<-stack(pred.catch.propHBT_E)
pred.catch.propHBT_E2$age<-rep(0:(nages-1),each=nyrs_age_comps_HBT_E)
pred.catch.propHBT_E2$year<-rep(yrs_age_comps_HBT_E,times=nages)
pred.catch.propHBT_E2$type<-rep("Predicted",times=nages*nyrs_age_comps_HBT_E)

OBS.catch.propHBT_E<- as.data.frame(OBS_catch_prop_HBT_E)
OBS.catch.propHBT_E2<-stack(OBS.catch.propHBT_E)
OBS.catch.propHBT_E2$age<-rep(0:(nages-1),each=nyrs_age_comps_HBT_E)
OBS.catch.propHBT_E2$year<-rep(yrs_age_comps_HBT_E,times=nages)
OBS.catch.propHBT_E2$type<-rep("Observed",times=nages*nyrs_age_comps_HBT_E)


catch.propHBT_E<-as.data.frame(cbind(OBS.catch.propHBT_E2,pred.catch.propHBT_E2$values,
              pred.catch.propHBT_E2$type))

print(xyplot(pred.catch.propHBT_E2$values+values~age|paste(year),data=catch.propHBT_E, 
      scales=list(alternating=c(1,2)),
      main='Headboat East Catch Age Proportions',xlab='Ages', ylab='Proportion at Age',
      outer=FALSE,distribute.type = TRUE, col=c("blue","red"),type=c('l','o'),pch=c(0,16),lty=c(1,0),lwd=c(2,0),
      auto.key = list ("topright",border=TRUE,padding.text=2.5,columns =2,text=c("Pred","Obs")), 
      par.settings = simpleTheme (pch=c(15,16), lwd=c(8,0),lty=c(1,0),col=c("blue","red")))
)




RMSE_catch_prop_HBT_E<-sqrt(mean((pred_catch_prop_HBT_E-OBS_catch_prop_HBT_E)^2))
stan_res_catch_prop_HBT_E<-(OBS_catch_prop_HBT_E-pred_catch_prop_HBT_E)/RMSE_catch_prop_HBT_E

layout(matrix(c(1), 1, byrow = TRUE))
years <- (yrs_age_comps_HBT_E[1]-2):(yrs_age_comps_HBT_E[nyrs_age_comps_HBT_E])
nyears <- length(years)
nindices <- nages
res.catch.prop_HBT_E <- stan_res_catch_prop_HBT_E #matrix(rnorm(nyears*nindices), ncol=nindices, nrow=nyears) 
resid.pos.color <- "blue"
resid.neg.color <- "red"
resid.col=matrix(NA, nrow=nyears, ncol=nindices)   # set color for residual bubbles
resid.col <- ifelse(res.catch.prop_HBT_E > 0.0,resid.pos.color, resid.neg.color) 
# use a scaler to increase (>1) or decrease (0<val<1) the size of the bubbles to make them look good
res.catch.prop_HBT_E1 <- res.catch.prop_HBT_E
# first set up a blank space for the plot
plot(seq(0,(nindices-1)), rep(years[1],nindices), ylim=c(years[nyears],years[1]), xlab = "Age",
     ylab = "Year",  type = "n", axes=F,main='Headboat East Catch-at-Age Proportion Standardized\nResiduals [(OBS-PRED)/SD]',
sub=bquote("(Absolute Value Max Residual"==.(max(abs(stan_res_catch_prop_HBT_E)))~", RMSE"==.(RMSE_catch_prop_HBT_E)~")"))


axis(1, at= seq(0,(nindices-1)), lab=seq(0,(nindices-1)))
axis(2, at = seq(years[nyears],years[1], by=-3), lab = seq(years[nyears],years[1], by=-3),las=1)
box()
abline(h=seq( years[3], years[nyears]),col = "lightgray", lty = 1)
abline(a=1965,b=1,col = "lightgray", lty = 1)
abline(a=1970,b=1,col = "lightgray", lty = 1)
abline(a=1975,b=1,col = "lightgray", lty = 1)
abline(a=1980,b=1,col = "lightgray", lty = 1)
abline(a=1985,b=1,col = "lightgray", lty = 1)
abline(a=1990,b=1,col = "lightgray", lty = 1)
abline(a=1995,b=1,col = "lightgray", lty = 1)
abline(a=2000,b=1,col = "lightgray", lty = 1)
abline(a=2005,b=1,col = "lightgray", lty = 1)
abline(a=2010,b=1,col = "lightgray", lty = 1)

# fill in the bubbles using points command, cex determines size, bg determines color
for (i in 3:nyears){
  points(seq(0, (nindices-1)), rep((years[i]), nindices),  
  cex =abs(res.catch.prop_HBT_E1[i-2,]) ,  pch = 21, bg = resid.col[i-2,], col='black')
}
legend("top",ncol=2, c(" Negative", " Positive"),pch=c(16,16),
col = c("red","blue")
,cex=.8, merge=FALSE)


RMSE_catch_prop_HBT_E<-sqrt(mean((pred_catch_prop_HBT_E-OBS_catch_prop_HBT_E)^2))
stan_res_catch_prop_HBT_E<-(OBS_catch_prop_HBT_E-pred_catch_prop_HBT_E)/RMSE_catch_prop_HBT_E

res.catch.prop_HBT_E<- as.data.frame(stan_res_catch_prop_HBT_E)
res.catch.prop_HBT_E2<-stack(res.catch.prop_HBT_E)
res.catch.prop_HBT_E2$age<-rep(0:(nages-1),each=nyrs_age_comps_HBT_E)
res.catch.prop_HBT_E2$year<-rep(yrs_age_comps_HBT_E,times=nages)

print(xyplot(values~age|paste(year),
      data=res.catch.prop_HBT_E2,scales=list(alternating=c(1,2)),
      main='Headboat East Catch-at-Age Proportion \nStandardized Residuals',xlab='Ages', ylab='Residuals',
      outer=FALSE,distribute.type = TRUE, type=c("b",'g'),lwd=1,pch=16,col='blue2')
)




pred.catch.propHBT_W<- as.data.frame(pred_catch_prop_HBT_W)
pred.catch.propHBT_W2<-stack(pred.catch.propHBT_W)
pred.catch.propHBT_W2$age<-rep(0:(nages-1),each=nyrs_age_comps_HBT_W)
pred.catch.propHBT_W2$year<-rep(yrs_age_comps_HBT_W,times=nages)
pred.catch.propHBT_W2$type<-rep("Predicted",times=nages*nyrs_age_comps_HBT_W)

OBS.catch.propHBT_W<- as.data.frame(OBS_catch_prop_HBT_W)
OBS.catch.propHBT_W2<-stack(OBS.catch.propHBT_W)
OBS.catch.propHBT_W2$age<-rep(0:(nages-1),each=nyrs_age_comps_HBT_W)
OBS.catch.propHBT_W2$year<-rep(yrs_age_comps_HBT_W,times=nages)
OBS.catch.propHBT_W2$type<-rep("Observed",times=nages*nyrs_age_comps_HBT_W)


catch.propHBT_W<-as.data.frame(cbind(OBS.catch.propHBT_W2,pred.catch.propHBT_W2$values,
              pred.catch.propHBT_W2$type))

print(xyplot(pred.catch.propHBT_W2$values+values~age|paste(year),data=catch.propHBT_W, 
      scales=list(alternating=c(1,2)),
      main='Headboat West Catch Age Proportions',xlab='Ages', ylab='Proportion at Age',
      outer=FALSE,distribute.type = TRUE, col=c("blue","red"),type=c('l','o'),pch=c(0,16),lty=c(1,0),lwd=c(2,0),
      auto.key = list ("topright",border=TRUE,padding.text=2.5,columns =2,text=c("Pred","Obs")), 
      par.settings = simpleTheme (pch=c(15,16), lwd=c(8,0),lty=c(1,0),col=c("blue","red")))
)


RMSE_catch_prop_HBT_W<-sqrt(mean((pred_catch_prop_HBT_W-OBS_catch_prop_HBT_W)^2))
stan_res_catch_prop_HBT_W<-(OBS_catch_prop_HBT_W-pred_catch_prop_HBT_W)/RMSE_catch_prop_HBT_W

layout(matrix(c(1), 1, byrow = TRUE))
years <- (yrs_age_comps_HBT_W[1]-2):(yrs_age_comps_HBT_W[nyrs_age_comps_HBT_W])
nyears <- length(years)
nindices <- nages
res.catch.prop_HBT_W <- stan_res_catch_prop_HBT_W #matrix(rnorm(nyears*nindices), ncol=nindices, nrow=nyears) 
resid.pos.color <- "blue"
resid.neg.color <- "red"
resid.col=matrix(NA, nrow=nyears, ncol=nindices)   # set color for residual bubbles
resid.col <- ifelse(res.catch.prop_HBT_W > 0.0,resid.pos.color, resid.neg.color) 
# use a scaler to increase (>1) or decrease (0<val<1) the size of the bubbles to make them look good
res.catch.prop_HBT_W1 <- res.catch.prop_HBT_W
# first set up a blank space for the plot
plot(seq(0,(nindices-1)), rep(years[1],nindices), ylim=c(years[nyears],years[1]), xlab = "Age",
     ylab = "Year",  type = "n", axes=F,main='Headboat West Catch-at-Age Proportion Standardized\nResiduals [(OBS-PRED)/SD]',
sub=bquote("(Absolute Value Max Residual"==.(max(abs(stan_res_catch_prop_HBT_W)))~", RMSE"==.(RMSE_catch_prop_HBT_W)~")"))


axis(1, at= seq(0,(nindices-1)), lab=seq(0,(nindices-1)))
axis(2, at = seq(years[nyears],years[1], by=-3), lab = seq(years[nyears],years[1], by=-3),las=1)
box()
abline(h=seq( years[3], years[nyears]),col = "lightgray", lty = 1)
abline(a=1965,b=1,col = "lightgray", lty = 1)
abline(a=1970,b=1,col = "lightgray", lty = 1)
abline(a=1975,b=1,col = "lightgray", lty = 1)
abline(a=1980,b=1,col = "lightgray", lty = 1)
abline(a=1985,b=1,col = "lightgray", lty = 1)
abline(a=1990,b=1,col = "lightgray", lty = 1)
abline(a=1995,b=1,col = "lightgray", lty = 1)
abline(a=2000,b=1,col = "lightgray", lty = 1)
abline(a=2005,b=1,col = "lightgray", lty = 1)
abline(a=2010,b=1,col = "lightgray", lty = 1)

# fill in the bubbles using points command, cex determines size, bg determines color
for (i in 3:nyears){
  points(seq(0, (nindices-1)), rep((years[i]), nindices),  
  cex =abs(res.catch.prop_HBT_W1[i-2,]) ,  pch = 21, bg = resid.col[i-2,], col='black')
}
legend("top",ncol=2, c(" Negative", " Positive"),pch=c(16,16),
col = c("red","blue")
,cex=.8, merge=FALSE)


RMSE_catch_prop_HBT_W<-sqrt(mean((pred_catch_prop_HBT_W-OBS_catch_prop_HBT_W)^2))
stan_res_catch_prop_HBT_W<-(OBS_catch_prop_HBT_W-pred_catch_prop_HBT_W)/RMSE_catch_prop_HBT_W

res.catch.prop_HBT_W<- as.data.frame(stan_res_catch_prop_HBT_W)
res.catch.prop_HBT_W2<-stack(res.catch.prop_HBT_W)
res.catch.prop_HBT_W2$age<-rep(0:(nages-1),each=nyrs_age_comps_HBT_W)
res.catch.prop_HBT_W2$year<-rep(yrs_age_comps_HBT_W,times=nages)

print(xyplot(values~age|paste(year),
      data=res.catch.prop_HBT_W2,scales=list(alternating=c(1,2)),
      main='Headboat West Catch-at-Age Proportion \nStandardized Residuals',xlab='Ages', ylab='Residuals',
      outer=FALSE,distribute.type = TRUE, type=c("b",'g'),lwd=1,pch=16,col='blue2')
)



pred.catch.propSHR_E<- as.data.frame(pred_catch_prop_SHR_E)
pred.catch.propSHR_E2<-stack(pred.catch.propSHR_E)
pred.catch.propSHR_E2$age<-rep(0:(nages-1),each=nyrs_age_comps_SHR_E)
pred.catch.propSHR_E2$year<-rep(yrs_age_comps_SHR_W,times=nages)
pred.catch.propSHR_E2$type<-rep("Predicted",times=nages*nyrs_age_comps_SHR_E)

OBS.catch.propSHR_E<- as.data.frame(OBS_catch_prop_SHR_E)
OBS.catch.propSHR_E2<-stack(OBS.catch.propSHR_E)
OBS.catch.propSHR_E2$age<-rep(0:(nages-1),each=nyrs_age_comps_SHR_E)
OBS.catch.propSHR_E2$year<-rep(yrs_age_comps_SHR_E,times=nages)
OBS.catch.propSHR_E2$type<-rep("Observed",times=nages*nyrs_age_comps_SHR_E)


catch.propSHR_E<-as.data.frame(cbind(OBS.catch.propSHR_E2,pred.catch.propSHR_E2$values,
              pred.catch.propSHR_E2$type))

print(xyplot(pred.catch.propSHR_E2$values+values~age|paste(year),data=catch.propSHR_E, 
      scales=list(alternating=c(1,2)),
      main='Shrimp East Catch Age Proportions',xlab='Ages', ylab='Proportion at Age',
      outer=FALSE,distribute.type = TRUE, col=c("blue","red"),type=c('l','o'),pch=c(0,16),lty=c(1,0),lwd=c(2,0),
      auto.key = list ("topright",border=TRUE,padding.text=2.5,columns =2,text=c("Pred","Obs")), 
      par.settings = simpleTheme (pch=c(15,16), lwd=c(8,0),lty=c(1,0),col=c("blue","red")))
)


RMSE_catch_prop_SHR_E<-sqrt(mean((pred_catch_prop_SHR_E-OBS_catch_prop_SHR_E)^2))
stan_res_catch_prop_SHR_E<-(OBS_catch_prop_SHR_E-pred_catch_prop_SHR_E)/RMSE_catch_prop_SHR_E

layout(matrix(c(1), 1, byrow = TRUE))
years <- (yrs_age_comps_SHR_E[1]-2):(yrs_age_comps_SHR_E[nyrs_age_comps_SHR_E])
nyears <- length(years)
nindices <- nages
res.catch.prop_SHR_E <- stan_res_catch_prop_SHR_E #matrix(rnorm(nyears*nindices), ncol=nindices, nrow=nyears) 
resid.pos.color <- "blue"
resid.neg.color <- "red"
resid.col=matrix(NA, nrow=nyears, ncol=nindices)   # set color for residual bubbles
resid.col <- ifelse(res.catch.prop_SHR_E > 0.0,resid.pos.color, resid.neg.color) 
# use a scaler to increase (>1) or decrease (0<val<1) the size of the bubbles to make them look good
res.catch.prop_SHR_E1 <- res.catch.prop_SHR_E
# first set up a blank space for the plot
plot(seq(0,(nindices-1)), rep(years[1],nindices), ylim=c(years[nyears],years[1]), xlab = "Age",
     ylab = "Year",  type = "n", axes=F,main='Shrimp East Catch-at-Age Proportion Standardized\nResiduals [(OBS-PRED)/SD]',
sub=bquote("(Absolute Value Max Residual"==.(max(abs(stan_res_catch_prop_SHR_E)))~", RMSE"==.(RMSE_catch_prop_SHR_E)~")"))


axis(1, at= seq(0,(nindices-1)), lab=seq(0,(nindices-1)))
axis(2, at = seq(years[nyears],years[1], by=-3), lab = seq(years[nyears],years[1], by=-3),las=1)
box()
abline(h=seq( years[3], years[nyears]),col = "lightgray", lty = 1)
abline(a=1965,b=1,col = "lightgray", lty = 1)
abline(a=1970,b=1,col = "lightgray", lty = 1)
abline(a=1975,b=1,col = "lightgray", lty = 1)
abline(a=1980,b=1,col = "lightgray", lty = 1)
abline(a=1985,b=1,col = "lightgray", lty = 1)
abline(a=1990,b=1,col = "lightgray", lty = 1)
abline(a=1995,b=1,col = "lightgray", lty = 1)
abline(a=2000,b=1,col = "lightgray", lty = 1)
abline(a=2005,b=1,col = "lightgray", lty = 1)
abline(a=2010,b=1,col = "lightgray", lty = 1)

# fill in the bubbles using points command, cex determines size, bg determines color
for (i in 3:nyears){
  points(seq(0, (nindices-1)), rep((years[i]), nindices),  
  cex =abs(res.catch.prop_SHR_E1[i-2,]) ,  pch = 21, bg = resid.col[i-2,], col='black')
}
legend("top",ncol=2, c(" Negative", " Positive"),pch=c(16,16),
col = c("red","blue")
,cex=.8, merge=FALSE)


RMSE_catch_prop_SHR_E<-sqrt(mean((pred_catch_prop_SHR_E-OBS_catch_prop_SHR_E)^2))
stan_res_catch_prop_SHR_E<-(OBS_catch_prop_SHR_E-pred_catch_prop_SHR_E)/RMSE_catch_prop_SHR_E

res.catch.prop_SHR_E<- as.data.frame(stan_res_catch_prop_SHR_E)
res.catch.prop_SHR_E2<-stack(res.catch.prop_SHR_E)
res.catch.prop_SHR_E2$age<-rep(0:(nages-1),each=nyrs_age_comps_SHR_E)
res.catch.prop_SHR_E2$year<-rep(yrs_age_comps_SHR_E,times=nages)

print(xyplot(values~age|paste(year),
      data=res.catch.prop_SHR_E2,scales=list(alternating=c(1,2)),
      main='Shrimp East Catch-at-Age Proportion \nStandardized Residuals',xlab='Ages', ylab='Residuals',
      outer=FALSE,distribute.type = TRUE, type=c("b",'g'),lwd=1,pch=16,col='blue2')
)






pred.catch.propSHR_W<- as.data.frame(pred_catch_prop_SHR_W)
pred.catch.propSHR_W2<-stack(pred.catch.propSHR_W)
pred.catch.propSHR_W2$age<-rep(0:(nages-1),each=nyrs_age_comps_SHR_W)
pred.catch.propSHR_W2$year<-rep(yrs_age_comps_SHR_W,times=nages)
pred.catch.propSHR_W2$type<-rep("Predicted",times=nages*nyrs_age_comps_SHR_W)

OBS.catch.propSHR_W<- as.data.frame(OBS_catch_prop_SHR_W)
OBS.catch.propSHR_W2<-stack(OBS.catch.propSHR_W)
OBS.catch.propSHR_W2$age<-rep(0:(nages-1),each=nyrs_age_comps_SHR_W)
OBS.catch.propSHR_W2$year<-rep(yrs_age_comps_SHR_W,times=nages)
OBS.catch.propSHR_W2$type<-rep("Observed",times=nages*nyrs_age_comps_SHR_W)


catch.propSHR_W<-as.data.frame(cbind(OBS.catch.propSHR_W2,pred.catch.propSHR_W2$values,
              pred.catch.propSHR_W2$type))

print(xyplot(pred.catch.propSHR_W2$values+values~age|paste(year),data=catch.propSHR_W, 
      scales=list(alternating=c(1,2)),
      main='Shrimp West Catch Age Proportions',xlab='Ages', ylab='Proportion at Age',
      outer=FALSE,distribute.type = TRUE, col=c("blue","red"),type=c('l','o'),pch=c(0,16),lty=c(1,0),lwd=c(2,0),
      auto.key = list ("topright",border=TRUE,padding.text=2.5,columns =2,text=c("Pred","Obs")), 
      par.settings = simpleTheme (pch=c(15,16), lwd=c(8,0),lty=c(1,0),col=c("blue","red")))
)


RMSE_catch_prop_SHR_W<-sqrt(mean((pred_catch_prop_SHR_W-OBS_catch_prop_SHR_W)^2))
stan_res_catch_prop_SHR_W<-(OBS_catch_prop_SHR_W-pred_catch_prop_SHR_W)/RMSE_catch_prop_SHR_W

layout(matrix(c(1), 1, byrow = TRUE))
years <- (yrs_age_comps_SHR_W[1]-2):(yrs_age_comps_SHR_W[nyrs_age_comps_SHR_W])
nyears <- length(years)
nindices <- nages
res.catch.prop_SHR_W <- stan_res_catch_prop_SHR_W #matrix(rnorm(nyears*nindices), ncol=nindices, nrow=nyears) 
resid.pos.color <- "blue"
resid.neg.color <- "red"
resid.col=matrix(NA, nrow=nyears, ncol=nindices)   # set color for residual bubbles
resid.col <- ifelse(res.catch.prop_SHR_W > 0.0,resid.pos.color, resid.neg.color) 
# use a scaler to increase (>1) or decrease (0<val<1) the size of the bubbles to make them look good
res.catch.prop_SHR_W1 <- res.catch.prop_SHR_W
# first set up a blank space for the plot
plot(seq(0,(nindices-1)), rep(years[1],nindices), ylim=c(years[nyears],years[1]), xlab = "Age",
     ylab = "Year",  type = "n", axes=F,main='Shrimp West Catch-at-Age Proportion Standardized\nResiduals [(OBS-PRED)/SD]',
sub=bquote("(Absolute Value Max Residual"==.(max(abs(stan_res_catch_prop_SHR_W)))~", RMSE"==.(RMSE_catch_prop_SHR_W)~")"))


axis(1, at= seq(0,(nindices-1)), lab=seq(0,(nindices-1)))
axis(2, at = seq(years[nyears],years[1], by=-3), lab = seq(years[nyears],years[1], by=-3),las=1)
box()
abline(h=seq( years[3], years[nyears]),col = "lightgray", lty = 1)
abline(a=1965,b=1,col = "lightgray", lty = 1)
abline(a=1970,b=1,col = "lightgray", lty = 1)
abline(a=1975,b=1,col = "lightgray", lty = 1)
abline(a=1980,b=1,col = "lightgray", lty = 1)
abline(a=1985,b=1,col = "lightgray", lty = 1)
abline(a=1990,b=1,col = "lightgray", lty = 1)
abline(a=1995,b=1,col = "lightgray", lty = 1)
abline(a=2000,b=1,col = "lightgray", lty = 1)
abline(a=2005,b=1,col = "lightgray", lty = 1)
abline(a=2010,b=1,col = "lightgray", lty = 1)

# fill in the bubbles using points command, cex determines size, bg determines color
for (i in 3:nyears){
  points(seq(0, (nindices-1)), rep((years[i]), nindices),  
  cex =abs(res.catch.prop_SHR_W1[i-2,]) ,  pch = 21, bg = resid.col[i-2,], col='black')
}
legend("top",ncol=2, c(" Negative", " Positive"),pch=c(16,16),
col = c("red","blue")
,cex=.8, merge=FALSE)


RMSE_catch_prop_SHR_W<-sqrt(mean((pred_catch_prop_SHR_W-OBS_catch_prop_SHR_W)^2))
stan_res_catch_prop_SHR_W<-(OBS_catch_prop_SHR_W-pred_catch_prop_SHR_W)/RMSE_catch_prop_SHR_W

res.catch.prop_SHR_W<- as.data.frame(stan_res_catch_prop_SHR_W)
res.catch.prop_SHR_W2<-stack(res.catch.prop_SHR_W)
res.catch.prop_SHR_W2$age<-rep(0:(nages-1),each=nyrs_age_comps_SHR_W)
res.catch.prop_SHR_W2$year<-rep(yrs_age_comps_SHR_W,times=nages)

print(xyplot(values~age|paste(year),
      data=res.catch.prop_SHR_W2,scales=list(alternating=c(1,2)),
      main='Shrimp West Catch-at-Age Proportion \nStandardized Residuals',xlab='Ages', ylab='Residuals',
      outer=FALSE,distribute.type = TRUE, type=c("b",'g'),lwd=1,pch=16,col='blue2')
)










pred.catch.propSUMM_E<- as.data.frame(Index_prop_SUMM_E)
pred.catch.propSUMM_E2<-stack(pred.catch.propSUMM_E)
pred.catch.propSUMM_E2$age<-rep(0:(nages-1),each=nyrs_age_comps_SUMM_E)
pred.catch.propSUMM_E2$year<-rep(yrs_age_comps_SUMM_E,times=nages)
pred.catch.propSUMM_E2$type<-rep("Predicted",times=nages*nyrs_age_comps_SUMM_E)

OBS.catch.propSUMM_E<- as.data.frame(OBS_index_prop_SUMM_E)
OBS.catch.propSUMM_E2<-stack(OBS.catch.propSUMM_E)
OBS.catch.propSUMM_E2$age<-rep(0:(nages-1),each=nyrs_age_comps_SUMM_E)
OBS.catch.propSUMM_E2$year<-rep(yrs_age_comps_SUMM_E,times=nages)
OBS.catch.propSUMM_E2$type<-rep("Observed",times=nages*nyrs_age_comps_SUMM_E)


catch.propSUMM_E<-as.data.frame(cbind(OBS.catch.propSUMM_E2,pred.catch.propSUMM_E2$values,
              pred.catch.propSUMM_E2$type))

print(xyplot(pred.catch.propSUMM_E2$values+values~age|paste(year),data=catch.propSUMM_E, 
      scales=list(alternating=c(1,2)),
      main='Summer East Catch Age Proportions',xlab='Ages', ylab='Proportion at Age',
      outer=FALSE,distribute.type = TRUE, col=c("blue","red"),type=c('l','o'),pch=c(0,16),lty=c(1,0),lwd=c(2,0),
      auto.key = list ("topright",border=TRUE,padding.text=2.5,columns =2,text=c("Pred","Obs")), 
      par.settings = simpleTheme (pch=c(15,16), lwd=c(8,0),lty=c(1,0),col=c("blue","red")))
)




RMSE_catch_prop_SUMM_E<-sqrt(mean((Index_prop_SUMM_E-OBS_index_prop_SUMM_E)^2))
stan_res_catch_prop_SUMM_E<-(OBS_index_prop_SUMM_E-Index_prop_SUMM_E)/RMSE_catch_prop_SUMM_E

layout(matrix(c(1), 1, byrow = TRUE))
years <- (yrs_age_comps_SUMM_E[1]-2):(yrs_age_comps_SUMM_E[nyrs_age_comps_SUMM_E])
nyears <- length(years)
nindices <- nages
res.catch.prop_SUMM_E <- stan_res_catch_prop_SUMM_E #matrix(rnorm(nyears*nindices), ncol=nindices, nrow=nyears) 
resid.pos.color <- "blue"
resid.neg.color <- "red"
resid.col=matrix(NA, nrow=nyears, ncol=nindices)   # set color for residual bubbles
resid.col <- ifelse(res.catch.prop_SUMM_E > 0.0,resid.pos.color, resid.neg.color) 
# use a scaler to increase (>1) or decrease (0<val<1) the size of the bubbles to make them look good
res.catch.prop_SUMM_E1 <- res.catch.prop_SUMM_E
# first set up a blank space for the plot
plot(seq(0,(nindices-1)), rep(years[1],nindices), ylim=c(years[nyears],years[1]), xlab = "Age",
     ylab = "Year",  type = "n", axes=F,main='Summer East Catch-at-Age Proportion Standardized\nResiduals [(OBS-PRED)/SD]',
sub=bquote("(Absolute Value Max Residual"==.(max(abs(stan_res_catch_prop_SUMM_E)))~", RMSE"==.(RMSE_catch_prop_SUMM_E)~")"))


axis(1, at= seq(0,(nindices-1)), lab=seq(0,(nindices-1)))
axis(2, at = seq(years[nyears],years[1], by=-3), lab = seq(years[nyears],years[1], by=-3),las=1)
box()
abline(h=seq( years[3], years[nyears]),col = "lightgray", lty = 1)
abline(a=1965,b=1,col = "lightgray", lty = 1)
abline(a=1970,b=1,col = "lightgray", lty = 1)
abline(a=1975,b=1,col = "lightgray", lty = 1)
abline(a=1980,b=1,col = "lightgray", lty = 1)
abline(a=1985,b=1,col = "lightgray", lty = 1)
abline(a=1990,b=1,col = "lightgray", lty = 1)
abline(a=1995,b=1,col = "lightgray", lty = 1)
abline(a=2000,b=1,col = "lightgray", lty = 1)
abline(a=2005,b=1,col = "lightgray", lty = 1)
abline(a=2010,b=1,col = "lightgray", lty = 1)

# fill in the bubbles using points command, cex determines size, bg determines color
for (i in 3:nyears){
  points(seq(0, (nindices-1)), rep((years[i]), nindices),  
  cex =abs(res.catch.prop_SUMM_E1[i-2,]) ,  pch = 21, bg = resid.col[i-2,], col='black')
}
legend("top",ncol=2, c(" Negative", " Positive"),pch=c(16,16),
col = c("red","blue")
,cex=.8, merge=FALSE)


RMSE_catch_prop_SUMM_E<-sqrt(mean((Index_prop_SUMM_E-OBS_index_prop_SUMM_E)^2))
stan_res_catch_prop_SUMM_E<-(OBS_index_prop_SUMM_E-Index_prop_SUMM_E)/RMSE_catch_prop_SUMM_E

res.catch.prop_SUMM_E<- as.data.frame(stan_res_catch_prop_SUMM_E)
res.catch.prop_SUMM_E2<-stack(res.catch.prop_SUMM_E)
res.catch.prop_SUMM_E2$age<-rep(0:(nages-1),each=nyrs_age_comps_SUMM_E)
res.catch.prop_SUMM_E2$year<-rep(yrs_age_comps_SUMM_E,times=nages)

print(xyplot(values~age|paste(year),
      data=res.catch.prop_SUMM_E2,scales=list(alternating=c(1,2)),
      main='Summer East Catch-at-Age Proportion \nStandardized Residuals',xlab='Ages', ylab='Residuals',
      outer=FALSE,distribute.type = TRUE, type=c("b",'g'),lwd=1,pch=16,col='blue2')
)





pred.catch.propSUMM_W<- as.data.frame(Index_prop_SUMM_W)
pred.catch.propSUMM_W2<-stack(pred.catch.propSUMM_W)
pred.catch.propSUMM_W2$age<-rep(0:(nages-1),each=nyrs_age_comps_SUMM_W)
pred.catch.propSUMM_W2$year<-rep(yrs_age_comps_SUMM_W,times=nages)
pred.catch.propSUMM_W2$type<-rep("Predicted",times=nages*nyrs_age_comps_SUMM_W)

OBS.catch.propSUMM_W<- as.data.frame(OBS_index_prop_SUMM_W)
OBS.catch.propSUMM_W2<-stack(OBS.catch.propSUMM_W)
OBS.catch.propSUMM_W2$age<-rep(0:(nages-1),each=nyrs_age_comps_SUMM_W)
OBS.catch.propSUMM_W2$year<-rep(yrs_age_comps_SUMM_W,times=nages)
OBS.catch.propSUMM_W2$type<-rep("Observed",times=nages*nyrs_age_comps_SUMM_W)


catch.propSUMM_W<-as.data.frame(cbind(OBS.catch.propSUMM_W2,pred.catch.propSUMM_W2$values,
              pred.catch.propSUMM_W2$type))

print(xyplot(pred.catch.propSUMM_W2$values+values~age|paste(year),data=catch.propSUMM_W, 
      scales=list(alternating=c(1,2)),
      main='Summer West Catch Age Proportions',xlab='Ages', ylab='Proportion at Age',
      outer=FALSE,distribute.type = TRUE, col=c("blue","red"),type=c('l','o'),pch=c(0,16),lty=c(1,0),lwd=c(2,0),
      auto.key = list ("topright",border=TRUE,padding.text=2.5,columns =2,text=c("Pred","Obs")), 
      par.settings = simpleTheme (pch=c(15,16), lwd=c(8,0),lty=c(1,0),col=c("blue","red")))
)




RMSE_catch_prop_SUMM_W<-sqrt(mean((Index_prop_SUMM_W-OBS_index_prop_SUMM_W)^2))
stan_res_catch_prop_SUMM_W<-(OBS_index_prop_SUMM_W-Index_prop_SUMM_W)/RMSE_catch_prop_SUMM_W

layout(matrix(c(1), 1, byrow = TRUE))
years <- (yrs_age_comps_SUMM_W[1]-2):(yrs_age_comps_SUMM_W[nyrs_age_comps_SUMM_W])
nyears <- length(years)
nindices <- nages
res.catch.prop_SUMM_W <- stan_res_catch_prop_SUMM_W #matrix(rnorm(nyears*nindices), ncol=nindices, nrow=nyears) 
resid.pos.color <- "blue"
resid.neg.color <- "red"
resid.col=matrix(NA, nrow=nyears, ncol=nindices)   # set color for residual bubbles
resid.col <- ifelse(res.catch.prop_SUMM_W > 0.0,resid.pos.color, resid.neg.color) 
# use a scaler to increase (>1) or decrease (0<val<1) the size of the bubbles to make them look good
res.catch.prop_SUMM_W1 <- res.catch.prop_SUMM_W
# first set up a blank space for the plot
plot(seq(0,(nindices-1)), rep(years[1],nindices), ylim=c(years[nyears],years[1]), xlab = "Age",
     ylab = "Year",  type = "n", axes=F,main='Summer West Catch-at-Age Proportion Standardized\nResiduals [(OBS-PRED)/SD]',
sub=bquote("(Absolute Value Max Residual"==.(max(abs(stan_res_catch_prop_SUMM_W)))~", RMSE"==.(RMSE_catch_prop_SUMM_W)~")"))


axis(1, at= seq(0,(nindices-1)), lab=seq(0,(nindices-1)))
axis(2, at = seq(years[nyears],years[1], by=-3), lab = seq(years[nyears],years[1], by=-3),las=1)
box()
abline(h=seq( years[3], years[nyears]),col = "lightgray", lty = 1)
abline(a=1965,b=1,col = "lightgray", lty = 1)
abline(a=1970,b=1,col = "lightgray", lty = 1)
abline(a=1975,b=1,col = "lightgray", lty = 1)
abline(a=1980,b=1,col = "lightgray", lty = 1)
abline(a=1985,b=1,col = "lightgray", lty = 1)
abline(a=1990,b=1,col = "lightgray", lty = 1)
abline(a=1995,b=1,col = "lightgray", lty = 1)
abline(a=2000,b=1,col = "lightgray", lty = 1)
abline(a=2005,b=1,col = "lightgray", lty = 1)
abline(a=2010,b=1,col = "lightgray", lty = 1)

# fill in the bubbles using points command, cex determines size, bg determines color
for (i in 3:nyears){
  points(seq(0, (nindices-1)), rep((years[i]), nindices),  
  cex =abs(res.catch.prop_SUMM_W1[i-2,]) ,  pch = 21, bg = resid.col[i-2,], col='black')
}
legend("top",ncol=2, c(" Negative", " Positive"),pch=c(16,16),
col = c("red","blue")
,cex=.8, merge=FALSE)


RMSE_catch_prop_SUMM_W<-sqrt(mean((Index_prop_SUMM_W-OBS_index_prop_SUMM_W)^2))
stan_res_catch_prop_SUMM_W<-(OBS_index_prop_SUMM_W-Index_prop_SUMM_W)/RMSE_catch_prop_SUMM_W

res.catch.prop_SUMM_W<- as.data.frame(stan_res_catch_prop_SUMM_W)
res.catch.prop_SUMM_W2<-stack(res.catch.prop_SUMM_W)
res.catch.prop_SUMM_W2$age<-rep(0:(nages-1),each=nyrs_age_comps_SUMM_W)
res.catch.prop_SUMM_W2$year<-rep(yrs_age_comps_SUMM_W,times=nages)

print(xyplot(values~age|paste(year),
      data=res.catch.prop_SUMM_W2,scales=list(alternating=c(1,2)),
      main='Summer West Catch-at-Age Proportion \nStandardized Residuals',xlab='Ages', ylab='Residuals',
      outer=FALSE,distribute.type = TRUE, type=c("b",'g'),lwd=1,pch=16,col='blue2')
)



pred.catch.propFALL_E<- as.data.frame(Index_prop_FALL_E)
pred.catch.propFALL_E2<-stack(pred.catch.propFALL_E)
pred.catch.propFALL_E2$age<-rep(0:(nages-1),each=nyrs_age_comps_FALL_E)
pred.catch.propFALL_E2$year<-rep(yrs_age_comps_FALL_E,times=nages)
pred.catch.propFALL_E2$type<-rep("Predicted",times=nages*nyrs_age_comps_FALL_E)

OBS.catch.propFALL_E<- as.data.frame(OBS_index_prop_FALL_E)
OBS.catch.propFALL_E2<-stack(OBS.catch.propFALL_E)
OBS.catch.propFALL_E2$age<-rep(0:(nages-1),each=nyrs_age_comps_FALL_E)
OBS.catch.propFALL_E2$year<-rep(yrs_age_comps_FALL_E,times=nages)
OBS.catch.propFALL_E2$type<-rep("Observed",times=nages*nyrs_age_comps_FALL_E)


catch.propFALL_E<-as.data.frame(cbind(OBS.catch.propFALL_E2,pred.catch.propFALL_E2$values,
              pred.catch.propFALL_E2$type))

print(xyplot(pred.catch.propFALL_E2$values+values~age|paste(year),data=catch.propFALL_E, 
      scales=list(alternating=c(1,2)),
      main='Fall East Catch Age Proportions',xlab='Ages', ylab='Proportion at Age',
      outer=FALSE,distribute.type = TRUE, col=c("blue","red"),type=c('l','o'),pch=c(0,16),lty=c(1,0),lwd=c(2,0),
      auto.key = list ("topright",border=TRUE,padding.text=2.5,columns =2,text=c("Pred","Obs")), 
      par.settings = simpleTheme (pch=c(15,16), lwd=c(8,0),lty=c(1,0),col=c("blue","red")))
)



RMSE_catch_prop_FALL_E<-sqrt(mean((Index_prop_FALL_E-OBS_index_prop_FALL_E)^2))
stan_res_catch_prop_FALL_E<-(OBS_index_prop_FALL_E-Index_prop_FALL_E)/RMSE_catch_prop_FALL_E

layout(matrix(c(1), 1, byrow = TRUE))
years <- (yrs_age_comps_FALL_E[1]-2):(yrs_age_comps_FALL_E[nyrs_age_comps_FALL_E])
nyears <- length(years)
nindices <- nages
res.catch.prop_FALL_E <- stan_res_catch_prop_FALL_E #matrix(rnorm(nyears*nindices), ncol=nindices, nrow=nyears) 
resid.pos.color <- "blue"
resid.neg.color <- "red"
resid.col=matrix(NA, nrow=nyears, ncol=nindices)   # set color for residual bubbles
resid.col <- ifelse(res.catch.prop_FALL_E > 0.0,resid.pos.color, resid.neg.color) 
# use a scaler to increase (>1) or decrease (0<val<1) the size of the bubbles to make them look good
res.catch.prop_FALL_E1 <- res.catch.prop_FALL_E
# first set up a blank space for the plot
plot(seq(0,(nindices-1)), rep(years[1],nindices), ylim=c(years[nyears],years[1]), xlab = "Age",
     ylab = "Year",  type = "n", axes=F,main='Fall East Catch-at-Age Proportion Standardized\nResiduals [(OBS-PRED)/SD]',
sub=bquote("(Absolute Value Max Residual"==.(max(abs(stan_res_catch_prop_FALL_E)))~", RMSE"==.(RMSE_catch_prop_FALL_E)~")"))


axis(1, at= seq(0,(nindices-1)), lab=seq(0,(nindices-1)))
axis(2, at = seq(years[nyears],years[1], by=-3), lab = seq(years[nyears],years[1], by=-3),las=1)
box()
abline(h=seq( years[3], years[nyears]),col = "lightgray", lty = 1)
abline(a=1965,b=1,col = "lightgray", lty = 1)
abline(a=1970,b=1,col = "lightgray", lty = 1)
abline(a=1975,b=1,col = "lightgray", lty = 1)
abline(a=1980,b=1,col = "lightgray", lty = 1)
abline(a=1985,b=1,col = "lightgray", lty = 1)
abline(a=1990,b=1,col = "lightgray", lty = 1)
abline(a=1995,b=1,col = "lightgray", lty = 1)
abline(a=2000,b=1,col = "lightgray", lty = 1)
abline(a=2005,b=1,col = "lightgray", lty = 1)
abline(a=2010,b=1,col = "lightgray", lty = 1)

# fill in the bubbles using points command, cex determines size, bg determines color
for (i in 3:nyears){
  points(seq(0, (nindices-1)), rep((years[i]), nindices),  
  cex =abs(res.catch.prop_FALL_E1[i-2,]) ,  pch = 21, bg = resid.col[i-2,], col='black')
}
legend("top",ncol=2, c(" Negative", " Positive"),pch=c(16,16),
col = c("red","blue")
,cex=.8, merge=FALSE)


RMSE_catch_prop_FALL_E<-sqrt(mean((Index_prop_FALL_E-OBS_index_prop_FALL_E)^2))
stan_res_catch_prop_FALL_E<-(OBS_index_prop_FALL_E-Index_prop_FALL_E)/RMSE_catch_prop_FALL_E

res.catch.prop_FALL_E<- as.data.frame(stan_res_catch_prop_FALL_E)
res.catch.prop_FALL_E2<-stack(res.catch.prop_FALL_E)
res.catch.prop_FALL_E2$age<-rep(0:(nages-1),each=nyrs_age_comps_FALL_E)
res.catch.prop_FALL_E2$year<-rep(yrs_age_comps_FALL_E,times=nages)

print(xyplot(values~age|paste(year),
      data=res.catch.prop_FALL_E2,scales=list(alternating=c(1,2)),
      main='Fall East Catch-at-Age Proportion \nStandardized Residuals',xlab='Ages', ylab='Residuals',
      outer=FALSE,distribute.type = TRUE, type=c("b",'g'),lwd=1,pch=16,col='blue2')
)





pred.catch.propFALL_W<- as.data.frame(Index_prop_FALL_W)
pred.catch.propFALL_W2<-stack(pred.catch.propFALL_W)
pred.catch.propFALL_W2$age<-rep(0:(nages-1),each=nyrs_age_comps_FALL_W)
pred.catch.propFALL_W2$year<-rep(yrs_age_comps_FALL_W,times=nages)
pred.catch.propFALL_W2$type<-rep("Predicted",times=nages*nyrs_age_comps_FALL_W)

OBS.catch.propFALL_W<- as.data.frame(OBS_index_prop_FALL_W)
OBS.catch.propFALL_W2<-stack(OBS.catch.propFALL_W)
OBS.catch.propFALL_W2$age<-rep(0:(nages-1),each=nyrs_age_comps_FALL_W)
OBS.catch.propFALL_W2$year<-rep(yrs_age_comps_FALL_W,times=nages)
OBS.catch.propFALL_W2$type<-rep("Observed",times=nages*nyrs_age_comps_FALL_W)


catch.propFALL_W<-as.data.frame(cbind(OBS.catch.propFALL_W2,pred.catch.propFALL_W2$values,
              pred.catch.propFALL_W2$type))

print(xyplot(pred.catch.propFALL_W2$values+values~age|paste(year),data=catch.propFALL_W, 
      scales=list(alternating=c(1,2)),
      main='Fall West Catch Age Proportions',xlab='Ages', ylab='Proportion at Age',
      outer=FALSE,distribute.type = TRUE, col=c("blue","red"),type=c('l','o'),pch=c(0,16),lty=c(1,0),lwd=c(2,0),
      auto.key = list ("topright",border=TRUE,padding.text=2.5,columns =2,text=c("Pred","Obs")), 
      par.settings = simpleTheme (pch=c(15,16), lwd=c(8,0),lty=c(1,0),col=c("blue","red")))
)




RMSE_catch_prop_FALL_W<-sqrt(mean((Index_prop_FALL_W-OBS_index_prop_FALL_W)^2))
stan_res_catch_prop_FALL_W<-(OBS_index_prop_FALL_W-Index_prop_FALL_W)/RMSE_catch_prop_FALL_W

layout(matrix(c(1), 1, byrow = TRUE))
years <- (yrs_age_comps_FALL_W[1]-2):(yrs_age_comps_FALL_W[nyrs_age_comps_FALL_W])
nyears <- length(years)
nindices <- nages
res.catch.prop_FALL_W <- stan_res_catch_prop_FALL_W #matrix(rnorm(nyears*nindices), ncol=nindices, nrow=nyears) 
resid.pos.color <- "blue"
resid.neg.color <- "red"
resid.col=matrix(NA, nrow=nyears, ncol=nindices)   # set color for residual bubbles
resid.col <- ifelse(res.catch.prop_FALL_W > 0.0,resid.pos.color, resid.neg.color) 
# use a scaler to increase (>1) or decrease (0<val<1) the size of the bubbles to make them look good
res.catch.prop_FALL_W1 <- res.catch.prop_FALL_W
# first set up a blank space for the plot
plot(seq(0,(nindices-1)), rep(years[1],nindices), ylim=c(years[nyears],years[1]), xlab = "Age",
     ylab = "Year",  type = "n", axes=F,main='Fall West Catch-at-Age Proportion Standardized\nResiduals [(OBS-PRED)/SD]',
sub=bquote("(Absolute Value Max Residual"==.(max(abs(stan_res_catch_prop_FALL_W)))~", RMSE"==.(RMSE_catch_prop_FALL_W)~")"))


axis(1, at= seq(0,(nindices-1)), lab=seq(0,(nindices-1)))
axis(2, at = seq(years[nyears],years[1], by=-3), lab = seq(years[nyears],years[1], by=-3),las=1)
box()
abline(h=seq( years[3], years[nyears]),col = "lightgray", lty = 1)
abline(a=1965,b=1,col = "lightgray", lty = 1)
abline(a=1970,b=1,col = "lightgray", lty = 1)
abline(a=1975,b=1,col = "lightgray", lty = 1)
abline(a=1980,b=1,col = "lightgray", lty = 1)
abline(a=1985,b=1,col = "lightgray", lty = 1)
abline(a=1990,b=1,col = "lightgray", lty = 1)
abline(a=1995,b=1,col = "lightgray", lty = 1)
abline(a=2000,b=1,col = "lightgray", lty = 1)
abline(a=2005,b=1,col = "lightgray", lty = 1)
abline(a=2010,b=1,col = "lightgray", lty = 1)

# fill in the bubbles using points command, cex determines size, bg determines color
for (i in 3:nyears){
  points(seq(0, (nindices-1)), rep((years[i]), nindices),  
  cex =abs(res.catch.prop_FALL_W1[i-2,]) ,  pch = 21, bg = resid.col[i-2,], col='black')
}
legend("top",ncol=2, c(" Negative", " Positive"),pch=c(16,16),
col = c("red","blue")
,cex=.8, merge=FALSE)


RMSE_catch_prop_FALL_W<-sqrt(mean((Index_prop_FALL_W-OBS_index_prop_FALL_W)^2))
stan_res_catch_prop_FALL_W<-(OBS_index_prop_FALL_W-Index_prop_FALL_W)/RMSE_catch_prop_FALL_W

res.catch.prop_FALL_W<- as.data.frame(stan_res_catch_prop_FALL_W)
res.catch.prop_FALL_W2<-stack(res.catch.prop_FALL_W)
res.catch.prop_FALL_W2$age<-rep(0:(nages-1),each=nyrs_age_comps_FALL_W)
res.catch.prop_FALL_W2$year<-rep(yrs_age_comps_FALL_W,times=nages)

print(xyplot(values~age|paste(year),
      data=res.catch.prop_FALL_W2,scales=list(alternating=c(1,2)),
      main='Fall West Catch-at-Age Proportion \nStandardized Residuals',xlab='Ages', ylab='Residuals',
      outer=FALSE,distribute.type = TRUE, type=c("b",'g'),lwd=1,pch=16,col='blue2')
)




for(i in 1:ceiling(nrow(ESTIMATES)/50))
{
if(50*i<nrow(ESTIMATES))
{
textplot(ESTIMATES[c((50*i-49):(50*i)),],mar=(c(.5,.5,.5,.5)),show.rownames=FALSE,col.data=Color[c((50*i-49):(50*i)),])
}
if(50*i>nrow(ESTIMATES))
{
textplot(ESTIMATES[c((50*i-49):nrow(ESTIMATES)),],mar=(c(.5,.5,.5,.5)),show.rownames=FALSE,col.data=Color[c((50*i-49):nrow(Color)),])
}
}


dev.off()

  shell("copy rs_ibm.cor rs_ibm_correlation.cor") #rename correlation file
  shell("del rs_ibm.cor") #delete original correlation file to prevent the file_exists command above from giving false positive about convergence


