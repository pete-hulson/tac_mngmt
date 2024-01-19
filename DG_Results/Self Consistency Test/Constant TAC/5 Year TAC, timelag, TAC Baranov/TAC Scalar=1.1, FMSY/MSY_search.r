# script to iterate over F values by stock/region/fleet in search for FMSY

wd<-wd

optimize.by.fleet<-'TRUE'             # optimize MSY by fleets within areas (=TRUE) or simply by using an F multiplier within an area (=FALSE)....latter reduces dimensionality of F permuations when number of fleets is too large to do full analysis
F.name<-"Fmsy_input"                                           # F values from .dat file to replace when iterating over F

update=readLines("SIM_TAC.dat",n=-1)

nyrs<-as.numeric(update[(grep("nyrs",update)+1)])
nfleets<-as.numeric(update[(grep("nfleets",update)+1)])  

F.seq<-seq(F.start,F.end,interval)                                                      # Create the desired F values to iterate over
if(optimize.by.fleet=='TRUE')
{
  permutation<-permutations(length(F.seq),sum(nfleets),F.seq,repeats.allowed=TRUE)  # determine all permutations of F in each stock/region/fleet
}

ntrials<-length(permutation[,1])
seq.length<-length(F.seq)
dir.create(paste0(wd,"/MSY Results",sep=""))
dir.create(paste0(wd,"/MSY Results/Figures",sep=""))
dir.create(paste0(wd,"/MSY Results/Report Files",sep=""))

Roundup <- function(from,to) ceiling(from/to)*to

MSY<-function(WD,ntrial,F.grab,perm)
 {
  setwd(WD)

  dir.create(paste0(WD,"/MSY Results/Run",ntrial,sep=""))
  
  invisible(file.copy(from=paste0(WD,"/SIM_TAC.exe",sep=""),to=paste0(WD,"/MSY Results/Run",ntrial,"/SIM_TAC.exe",sep="")))
  invisible(file.copy(from=paste0(WD,"/SIM_TAC.dat",sep=""),to=paste0(WD,"/MSY Results/Run",ntrial,"/SIM_TAC.dat",sep="")))
  invisible(file.copy(from=paste0(WD,"/SIM_TAC.tpl",sep=""),to=paste0(WD,"/MSY Results/Run",ntrial,"/SIM_TAC.tpl",sep="")))

  setwd(paste0(WD,"/MSY Results/Run",ntrial,sep=""))

  update=readLines("SIM_TAC.dat",n=-1)
  
  nyrs<-as.numeric(update[(grep("nyrs",update)+1)])

    F.region<-(perm[ntrial,])

  update[(grep(F.grab,update)+1)]=F.region
  update[(grep("yr_fishing_start",update)+1)]=2
  update[(grep("sigma_landings",update)+1)]=paste0(as.character(rep(sigma_landings[1],times=nyrs-1)),collapse=" ") 

  writeLines(update,"SIM_TAC.dat")
  
  invisible(shell("SIM_TAC -nohess",wait=T)) #show.output.on.console=FALSE))  ### Run ADMB with update F
  
  invisible(file.remove(paste0(WD,"/MSY Results/Run",ntrial,"/SIM_TAC.exe",sep="")))
  
  invisible(file.copy(from=paste0(WD,"/MSY Results/Run",ntrial,"/SIM_TAC.rep",sep=""),to=paste0(WD,"/MSY Results/Report Files/Report",ntrial,".rep",sep="")))
  
  out =readList("SIM_TAC.rep")
  par_names=c('SIM_ssb','SIM_recruits','SIM_biomass','SIM_F','TRUE_landings')
  result=out[par_names]
  
  bio_total<-out$SIM_biomass
  yield_total<-out$TRUE_landings
  SSB_total<-out$SIM_ssb
  F_total<-out$SIM_F
  Recruits_total<-out$SIM_recruits

    c(i,F.region,bio_total[1],bio_total[nyrs],yield_total[1],yield_total[nyrs-1],SSB_total[1],SSB_total[nyrs],F_total[1],F_total[nyrs-1],
      Recruits_total[1],Recruits_total[nyrs])
}

    
cl <- makeSOCKcluster(detectCores())
clusterExport(cl, c("wd","ntrials","F.name","permutation"))
registerDoSNOW(cl)
pb <- winProgressBar(paste("MSY Search Progress Bar"), label=paste("Simulation Run 0 of ",ntrials,sep=""),max=100)
progress<-function(n) setWinProgressBar(pb,(n/ntrials*100),label=paste("Simulation Run", n,"of", ntrials,"Completed"))
opts<-list(progress=progress)

t<- foreach(i=1:ntrials,.combine=rbind,.options.snow=opts,.packages=c('PBSmodelling','matrixStats','TeachingDemos','snowfall','parallel')
) %dopar% {
  MSY(wd,i,F.name,permutation) 
}

colnames(t)=c("trial","F","bio_st","bio_end","yield_st","yield_end","SSB_st","SSB_end","F_st","F_end","Rec_st","Rec_end")

write.csv(t, file=paste0(wd,"/MSY Results/Figures/Output Quantities.csv",sep=""))

close(pb) 
stopCluster(cl)
closeAllConnections()
output3<-gc()


graph<-read.csv(paste0(wd,"/MSY Results/Figures/Output Quantities.csv",sep=""))

MSY1<-graph$yield_end
SPR<-graph$SSB_end/graph$SSB_st
bio<-graph$bio_end

stock<-graph[which.max(graph$yield_end),]

MSY.SPR<-data.frame(cbind(SPR,MSY1))
MSY.SPR<-data.table(MSY.SPR,key="SPR")

MSY.bio<-cbind(bio,MSY1)
MSY.bio<-data.table(MSY.bio,key="bio")


MSY.values<-graph[which.max(graph$yield_end),]
write.csv(MSY.values, file=paste0(wd,"/MSY Results/Figures/MSY Outputs.csv",sep=""))

png(file=paste0(wd,"/MSY Results/Figures/MSY Curve (TOTAL).png",sep=""))
par(mfrow=c(2,1)) #,  mgp=c(2.5,.6,.5), mar=c(5, 4, 3, .5) + 0.1)

par(las=1, mar=c(5,5,3,1))
plot(MSY.SPR$SPR,MSY.SPR$MSY, typ='p', xlab = 'Equilibrium SSB Ratio', ylab = 'Yield', lwd=2,lty=1, main='MSY vs. SSB Ratio',cex.main=1.15,
        ylim=c(0,Roundup(max(MSY.SPR$MSY),1000)),
        cex.axis=1., cex.lab=1.,xlim=c(0.0, Roundup(max(MSY.SPR$SPR),.1)) ,axes=FALSE)
axis( 2, pos=0,seq(0,Roundup(max(MSY.SPR$MSY),1000),Roundup(max(MSY.SPR$MSY),1000)/5),cex.axis=.8)
axis( 1,pos=0,seq(0, Roundup(max(MSY.SPR$SPR),.1), Roundup(max(MSY.SPR$SPR),.1)/5),cex.axis=.8)

segments(0, max(MSY.SPR$MSY), MSY.SPR$SPR[which.max(MSY.SPR$MSY)], max(MSY.SPR$MSY), lty=6, col='grey70', lwd=2)
segments(MSY.SPR$SPR[which.max(MSY.SPR$MSY)],0, MSY.SPR$SPR[which.max(MSY.SPR$MSY)], max(MSY.SPR$MSY), lty=6, col='grey70', lwd=2)
matpoints(MSY.SPR$SPR[which.max(MSY.SPR$MSY)],max(MSY.SPR$MSY),pch=19,col='grey70',cex=2)


par(las=1, mar=c(5,5,3,1))
plot(MSY.bio$bio,MSY.bio$MSY, typ='p', xlab = 'Equilibrium Biomass', ylab = 'Yield', lwd=2,lty=1, main='MSY vs. Biomass',cex.main=1.15,
     ylim=c(0,Roundup(max(MSY.bio$MSY),1000)),
     cex.axis=1., cex.lab=1.,xlim=c(0.0, Roundup(max(MSY.bio$bio),.1)) ,axes=FALSE)
axis( 2, pos=0,seq(0,Roundup(max(MSY.bio$MSY),1000),Roundup(max(MSY.bio$MSY),1000)/5),cex.axis=.8)
axis( 1,pos=0,seq(0, Roundup(max(MSY.bio$bio),1000), Roundup(max(MSY.bio$bio),1000)/5),cex.axis=.8)

segments(0, max(MSY.bio$MSY), MSY.bio$bio[which.max(MSY.bio$MSY)], max(MSY.bio$MSY), lty=6, col='grey70', lwd=2)
segments(MSY.bio$bio[which.max(MSY.bio$MSY)],0, MSY.bio$bio[which.max(MSY.bio$MSY)], max(MSY.bio$MSY), lty=6, col='grey70', lwd=2)
matpoints(MSY.bio$bio[which.max(MSY.bio$MSY)],max(MSY.bio$MSY),pch=19,col='grey70',cex=2);

dev.off()

