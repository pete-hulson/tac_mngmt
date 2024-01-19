
rm(list=(ls()))
wd<-getwd()
setwd(wd)

############### INPUTS ##############################################################################################################################
################ ALL FLAGS MUST BE IN ALL ****CAPS********########################################################
#### FMSY Sim Inputs ####################################
Perform.FMSY.SIM<-'TRUE'                                     # 'TRUE' or 'FALSE'
FMSY.SIM.Recruit.Type<-'DETERMINISTIC'                       # 'DETERMINISTIC' or 'STOCHASTIC'
F.start=0.25                                                 # Starting F value that will search from for Fmsy
F.end=.35                                               # Terminal F value that will search to for Fmsy 
interval=0.001                                               # Step increase in Fmsy search 
#######################################################################

nloops <-50                                                 # How many simulations to perform

############# Model Flags/Switches ###################################################################
################ ALL FLAGS MUST BE IN ALL ****CAPS********########################################################

SIM.Recruit.Type<-'STOCHASTIC'                               # 'DETERMINISTIC' or 'STOCHASTIC' recruitment for simulations

SIM.Survey.Sel.Type<-'LOGISTIC'                              # 'LOGISTIC' or 'BY_AGE'                      
Survey.Sel.Type<-'LOGISTIC'
SIM.Fishery.Sel.Type<-'LOGISTIC'
Fishery.Sel.Type<-'LOGISTIC'

Init.Abund.Est.Type<-'ESTIMATE'                     # 'EQUIL_OFFSET' or 'ESTIMATE', if equil then assume pop in equil with no fishing based on an initial Recruitment offset from R0 by R0_offset, if estimate then estimate abund at age

##################### Types of SIM Model Runs ########################################################################################################

SIM.PRE.TAC.F.LEVEL<-'FMSY'                                   # 'LOW', 'HIGH', 'FMSY', 'ALL', should pre-TAC F be simulated below, at or above FMSY (or all 3)...FMSY.scalar is used to scale the average F before the TAC is implemented (the true F is then determined by the random F_devs)
FMSY.SCALAR.LOW<-0.5                                         # Value for low FMSY Scalar
FMSY.SCALAR.HIGH<-1.5                                         # Value for HI FMSY Scalar

TAC.Constant<-'BOTH'                                         # 'TRUE' 'FALSE' or 'BOTH', if true TAC is held constant at level in first year, if false TAC is reset each year based on estimated abundance and FMSY, if both computes results for both

####################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################

DATA.FILE=readLines("SIM_TAC.dat",n=-1)

if(identical(Perform.FMSY.SIM,'TRUE')==TRUE)
{
  DATA.FILE[(grep("SIM_Fmsy_switch",DATA.FILE)+1)]=1
  if(identical(FMSY.SIM.Recruit.Type,'DETERMINISTIC')==TRUE)
  {
   DATA.FILE[(grep("SIM_deterministic_recruit",DATA.FILE)+1)]=1
  }
  if(identical(FMSY.SIM.Recruit.Type,'STOCHASTIC')==TRUE)
  {
    DATA.FILE[(grep("SIM_deterministic_recruit",DATA.FILE)+1)]=0
  }
}
if(identical(Perform.FMSY.SIM,'FALSE')==TRUE)
{
  DATA.FILE[(grep("SIM_Fmsy_switch",DATA.FILE)+1)]=0
}
#####################################################################
if(identical(SIM.Survey.Sel.Type,'LOGISTIC')==TRUE)
{
  DATA.FILE[(grep("SIM_survey_sel_switch",DATA.FILE)+1)]=2
}
if(identical(SIM.Survey.Sel.Type,'BY_AGE')==TRUE)
{
  DATA.FILE[(grep("SIM_survey_sel_switch",DATA.FILE)+1)]=1
}
if(identical(SIM.Fishery.Sel.Type,'LOGISTIC')==TRUE)
{
  DATA.FILE[(grep("SIM_sel_switch",DATA.FILE)+1)]=2
}
if(identical(SIM.Fishery.Sel.Type,'BY_AGE')==TRUE)
{
  DATA.FILE[(grep("SIM_sel_switch",DATA.FILE)+1)]=1
}
################################################################
if(identical(Survey.Sel.Type,'LOGISTIC')==TRUE)
{
  DATA.FILE[(grep("ASS_survey_sel_switch",DATA.FILE)+1)]=2
  DATA.FILE[(grep("phase_beta1_index",DATA.FILE)+1)]=4
  DATA.FILE[(grep("phase_beta2_index",DATA.FILE)+1)]=4
  DATA.FILE[(grep("phase_ln_sel_age_index",DATA.FILE)+1)]=-4
}
if(identical(Survey.Sel.Type,'BY_AGE')==TRUE)
{
  DATA.FILE[(grep("ASS_survey_sel_switch",DATA.FILE)+1)]=1
  DATA.FILE[(grep("phase_beta1_index",DATA.FILE)+1)]=-4
  DATA.FILE[(grep("phase_beta2_index",DATA.FILE)+1)]=-4
  DATA.FILE[(grep("phase_ln_sel_age_index",DATA.FILE)+1)]=4
}
if(identical(Fishery.Sel.Type,'LOGISTIC')==TRUE)
{
  DATA.FILE[(grep("ASS_sel_switch",DATA.FILE)+1)]=2
  DATA.FILE[(grep("phase_beta1_fish",DATA.FILE)+1)]=4
  DATA.FILE[(grep("phase_beta2_fish",DATA.FILE)+1)]=4
  DATA.FILE[(grep("phase_ln_sel_age_fish",DATA.FILE)+1)]=-4
}
if(identical(Fishery.Sel.Type,'BY_AGE')==TRUE)
{
  DATA.FILE[(grep("ASS_sel_switch",DATA.FILE)+1)]=1
  DATA.FILE[(grep("phase_beta1_fish",DATA.FILE)+1)]=-4
  DATA.FILE[(grep("phase_beta2_fish",DATA.FILE)+1)]=-4
  DATA.FILE[(grep("phase_ln_sel_age_fish",DATA.FILE)+1)]=4
}
######################################################################
######################################################################
if(identical(Init.Abund.Est.Type,'EQUIL_OFFSET')==TRUE)
{
  DATA.FILE[(grep("init_abund_switch",DATA.FILE)+1)]=0
  DATA.FILE[(grep("phase_ln_R0_offset",DATA.FILE)+1)]=1
  DATA.FILE[(grep("phase_ln_init_abund",DATA.FILE)+1)]=-1
}
if(identical(Init.Abund.Est.Type,'ESTIMATE')==TRUE)
{
  DATA.FILE[(grep("init_abund_switch",DATA.FILE)+1)]=1
  DATA.FILE[(grep("phase_ln_R0_offset",DATA.FILE)+1)]=-1
  DATA.FILE[(grep("phase_ln_init_abund",DATA.FILE)+1)]=1
}
#######################################################################
  writeLines(DATA.FILE,"SIM_TAC.dat")
############################################################################################################################  
  
#--------------------------------------------------------------------------------------------------------------------------------------------------
suppressWarnings(suppressMessages(require(PBSmodelling)))
suppressWarnings(suppressMessages(require(matrixStats)))
suppressWarnings(suppressMessages(require(TeachingDemos)))
suppressWarnings(suppressMessages(require(snowfall)))
suppressWarnings(suppressMessages(library(parallel)))
suppressWarnings(suppressMessages(library('snow')))
suppressWarnings(suppressMessages(library(foreach)))
suppressWarnings(suppressMessages(library(doSNOW)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(gtools)))
suppressWarnings(suppressMessages(library(spatstat)))
suppressWarnings(suppressMessages(library(alphahull)))
suppressWarnings(suppressMessages(library(beanplot)))
#--------------------------------------------------------------------------------------------------------------------------------------------------


##########################################################################################################################################
############ Run Fmsy Simulation ###########################################################################################################

update1=readLines("SIM_TAC.dat",n=-1)
nyrs<-as.numeric(update1[(grep("nyrs",update1)+1)])
yr_fishing_start<-as.numeric(update1[(grep("yr_fishing_start",update1)+1)])
SIM_Fmsy_Switch<-as.numeric(update1[(grep("SIM_Fmsy_switch",update1)+1)])
sigma_landings<-as.numeric(unlist(strsplit(update1[(grep("sigma_landings",update1)+1)]," ")))

if(SIM_Fmsy_Switch==1)
{                                            
  source(file="MSY_search.r")
  Fmsy.temp<-read.csv(paste0(wd,"/MSY Results/Figures/MSY Outputs.csv",sep=""))
  Fmsy<-Fmsy.temp$F
  update2=readLines("SIM_TAC.dat",n=-1)
  update2[(grep("Fmsy_input",update2)+1)]=Fmsy
  update2[(grep("SIM_Fmsy_switch",update2)+1)]=0
  if(identical(SIM.Recruit.Type,'DETERMINISTIC')==TRUE)
  {
    update2[(grep("SIM_deterministic_recruit",update2)+1)]=1
  }
  if(identical(SIM.Recruit.Type,'STOCHASTIC')==TRUE)
  {
    update2[(grep("SIM_deterministic_recruit",update2)+1)]=0
  }
  update2[(grep("yr_fishing_start",update2)+1)]=yr_fishing_start
  update2[(grep("sigma_landings",update2)+1)]=paste0(as.character(rep(sigma_landings[1],times=(nyrs-yr_fishing_start+1))),collapse=" ")
  writeLines(update2,"SIM_TAC.dat")
}

######################################################################################################################

run.SIM<-function(ntrial,WD1,WD,TAC.CNST,PRE.TAC.F.SCALAR)
{
  if(ntrial==1)
  {

  SIM.DATA =readLines(paste0(WD1,"/SIM_TAC.dat",sep=""),n=-1)
  
  if(TAC.CNST==1)
  {
    SIM.DATA[(grep("TAC_CNST_switch",SIM.DATA)+1)]=1
  }
  if(TAC.CNST==0)
  {
    SIM.DATA[(grep("TAC_CNST_switch",SIM.DATA)+1)]=0
  }    
  SIM.DATA[(grep("Fmsy_scalar",SIM.DATA)+1)]=PRE.TAC.F.SCALAR
  writeLines(SIM.DATA,paste0(WD1,"/SIM_TAC.dat",sep=""))
  }
  dir.create(paste0(WD,"/Simulation Results/Run",ntrial,sep=""))
  dir.create(paste0(WD,"/Assessment Results/Run",ntrial,sep=""))
  invisible(file.copy(from=paste0(WD1,"/SIM_TAC.exe",sep=""),to=paste0(WD,"/Simulation Results/Run",ntrial,"/SIM_TAC.exe",sep="")))
  invisible(file.copy(from=paste0(WD1,"/SIM_TAC.tpl",sep=""),to=paste0(WD,"/Simulation Results/Run",ntrial,"/SIM_TAC.tpl",sep="")))
  invisible(file.copy(from=paste0(WD1,"/SIM_TAC.dat",sep=""),to=paste0(WD,"/Simulation Results/Run",ntrial,"/SIM_TAC.dat",sep="")))
  
  setwd(paste0(WD,"/Simulation Results/Run",ntrial,sep=""))
  
  SIM.DAT=readLines("SIM_TAC.dat",n=-1)
  SIM.DAT[(grep("myseed",SIM.DAT)+1)]=ntrial
  writeLines(SIM.DAT,"SIM_TAC.dat")
  
  system(paste0("SIM_TAC -nox -nohess SIM_TAC.dat",sep=""),wait=TRUE,show.output.on.console=FALSE)  #-nox keeps from showing the ADMB gradient outputs on screen
  
  invisible(shell(paste("copy SIM_TAC.par SIM_TAC",ntrial,".par", sep="")))
  invisible(shell(paste("copy SIM_TAC.rep SIM_TAC",ntrial,".rep", sep="")))

  invisible(file.copy(from=paste0(WD,"/Simulation Results/Run",ntrial,"/SIM_TAC",ntrial,".rep",sep=""),to=paste0(WD,"/Assessment Results/Run",ntrial,"/TAC_est",ntrial,".dat",sep="")))

  invisible(shell(paste("del SIM_TAC.par", sep="")))
  invisible(shell(paste("del SIM_TAC.rep", sep="")))
  invisible(shell("del SIM_TAC.exe")) 
  invisible(shell(paste("del SIM_TAC.bar",sep="")))
  invisible(shell(paste("del SIM_TAC.log",sep="")))
  invisible(shell("del fmin.log"))
  invisible(shell("del variance"))
  
}

run.est<-function(ntrial,WD1,WD)
{

  invisible(file.copy(from=paste0(WD1,"/TAC_est.exe",sep=""),to=paste0(WD,"/Assessment Results/Run",ntrial,"/TAC_est.exe",sep="")))
  invisible(file.copy(from=paste0(WD1,"/TAC_est.tpl",sep=""),to=paste0(WD,"/Assessment Results/Run",ntrial,"/TAC_est.tpl",sep="")))

  setwd(paste0(WD,"/Assessment Results/Run",ntrial,sep=""))

  system(paste0("TAC_est -nox -ind TAC_est",ntrial,".dat",sep=""),wait=TRUE,show.output.on.console=FALSE)  #-nox keeps from showing the ADMB gradient outputs on screen
  
  invisible(shell(paste("copy TAC_est.par TAC_est",ntrial,".par", sep="")))
  invisible(shell(paste("copy TAC_est.rep TAC_est",ntrial,".rep", sep="")))

  if(file.exists("TAC_est.cor")==TRUE)
  {
    print(paste0("Model Run ",ntrial," Converged"))
    invisible(shell(paste("copy TAC_est.std TAC_est",ntrial,".std", sep="")))
    invisible(shell(paste("copy TAC_est.cor TAC_est",ntrial,".cor", sep="")))
    invisible(file.copy(from=paste0(WD,"/Assessment Results/Run",ntrial,"/TAC_est",ntrial,".rep",sep=""),to=paste0(WD,"/Assessment Results/Report Files/TAC_est",ntrial,".rep",sep="")))
    invisible(file.copy(from=paste0(WD,"/Assessment Results/Run",ntrial,"/TAC_est",ntrial,".cor",sep=""),to=paste0(WD,"/Assessment Results/Report Files/TAC_est",ntrial,".cor",sep="")))
    invisible(file.copy(from=paste0(WD,"/Assessment Results/Run",ntrial,"/TAC_est",ntrial,".std",sep=""),to=paste0(WD,"/Assessment Results/Report Files/TAC_est",ntrial,".std",sep="")))
    invisible(file.copy(from=paste0(WD,"/Assessment Results/Run",ntrial,"/TAC_est",ntrial,".par",sep=""),to=paste0(WD,"/Assessment Results/Report Files/TAC_est",ntrial,".par",sep="")))
    invisible(shell(paste("del TAC_est.std", sep="")))
    invisible(shell(paste("del TAC_est.cor", sep="")))  
    invisible(shell("del admodel.cov"))
  }
  if(file.exists("TAC_est.cor")==FALSE)
  {
    print(paste0("Model Run ",ntrial," Did Not Converge"))
    invisible(file.copy(from=paste0(WD,"/Assessment Results/Run",ntrial,"/TAC_est",ntrial,".rep",sep=""),to=paste0(WD,"/Assessment Results/Report Files/TAC_est",ntrial,".rep",sep="")))
  }
  
  #if(file.exists("IBM.cor"))
 # {
  #  source(file="IBM graphics.r")
  #  shell(paste("copy RS_IBM.pdf IBM",z,".pdf", sep="")) 
 # }
  
  invisible(shell(paste("del TAC_est.par", sep="")))
  invisible(shell(paste("del TAC_est.rep", sep="")))
  invisible(shell("del TAC_est.exe")) 
  invisible(shell("del admodel.dep"))
  invisible(shell("del admodel.hes")) 
  invisible(shell(paste("del TAC_est.b01",sep="")))
  invisible(shell(paste("del TAC_est.P01",sep="")))
  invisible(shell(paste("del TAC_est.R01",sep="")))
  invisible(shell(paste("del TAC_est.b02",sep="")))
  invisible(shell(paste("del TAC_est.P02",sep="")))
  invisible(shell(paste("del TAC_est.R02",sep="")))
  invisible(shell(paste("del TAC_est.b03",sep="")))
  invisible(shell(paste("del TAC_est.P03",sep="")))
  invisible(shell(paste("del TAC_est.R03",sep="")))
  invisible(shell(paste("del TAC_est.bar",sep="")))
  invisible(shell(paste("del TAC_est.eva",sep="")))
  invisible(shell(paste("del TAC_est.log",sep="")))
  invisible(shell("del fmin.log"))
  invisible(shell("del variance"))
  
}

plot.beans.combined<-function(parameter.name,data.frame.name,Converge,percent.converge,...)
{
  par(mfrow=c(1,1),  mgp=c(2.5,.6,.5), mar=c(5, 4, 3, .5) + 0.1)
  
  beanplot(as.numeric(unlist(data.frame.name)),bw="nrd0",boxwex=1.2,what=c(0,1,0,0),col=c('grey','black','black','white'), #ylim=c(min(data.frame.name,na.rm=TRUE),max(data.frame.name,na.rm=TRUE)),
           beanlines="median",main=paste0(parameter.name," Percent Bias"),  mgp=c(2.75,.75,.5),
           sub=paste0(Converge, " runs converged out of ", nloops," (", percent.converge,"% Convergence)"),las=1,log="",
           cex.axis=.5,na.rm=TRUE,...)
  abline(h=0,lty=1,lwd=2)
  abline(h=median(as.numeric(unlist(data.frame.name)),na.rm=TRUE),lty=2,col="red",lwd=2)
  legend("top", c(paste0("Median Bias ", signif(median(as.numeric(unlist(data.frame.name)),na.rm=TRUE),3),"%"),
                  paste0("Max Bias ", signif(max(abs(as.numeric(unlist(data.frame.name))),na.rm=TRUE),3),"%")),lty=c(4,-1),cex=.5,ncol=2,bg='white')
 
}   

plot.beans=function(parameter.name,data.frame.name,Converge,percent.converge,...)
{
  par(mfrow=c(1,1),  mgp=c(2.5,.6,.5), mar=c(5, 4, 3, .5) + 0.1)
  
  beanplot(data.frame.name,bw="nrd0",boxwex=1.2,what=c(0,1,0,0),col=c('grey','black','black','white'),#ylim=c(min(data.frame.name,na.rm=TRUE),max(data.frame.name,na.rm=TRUE)),
           beanlines="median",main=paste0(parameter.name," Percent Bias"),  mgp=c(2.75,.75,.5),
           sub=paste0(Converge, " runs converged out of ", nloops," (", percent.converge,"% Convergence)"),las=1,log="",
           cex.axis=.5,na.rm=TRUE,...)
  abline(h=0,lty=1,lwd=2)
  abline(h=median(as.numeric(unlist(data.frame.name)),na.rm=TRUE),lty=2,col="red",lwd=2)
  legend("top", c(paste0("Median Bias ", signif(median(as.numeric(unlist(data.frame.name)),na.rm=TRUE),3),"%"),
                  paste0("Max Bias ", signif(max(abs(as.numeric(unlist(data.frame.name))),na.rm=TRUE),3),"%")),lty=c(4,-1),cex=.5,ncol=2,bg='white')
  
}   

plot.beans.value=function(parameter.name,data.frame.name,Converge,percent.converge,...)
{
  par(mfrow=c(1,1),  mgp=c(2.5,.6,.5), mar=c(5, 4, 3, .5) + 0.1)
  
  beanplot(data.frame.name,bw="nrd0",boxwex=1.2,what=c(0,1,0,0),col=c('grey','black','black','white'),#ylim=c(min(data.frame.name,na.rm=TRUE),max(data.frame.name,na.rm=TRUE)),
           beanlines="median",main=paste0(parameter.name," Observed and Predicted Value"),  mgp=c(2.75,.75,.5),
           sub=paste0(Converge, " runs converged out of ", nloops," (", percent.converge,"% Convergence)"),las=1,log="",
           cex.axis=.5,na.rm=TRUE,...)
  abline(h=0,lty=1,lwd=2)
  legend("top",c('True Value'),lty=c(1),pch=1,bg='white')
  
}  

result=function(WD,iteration)
{
  setwd(paste0(WD,"/Assessment Results/Report Files"))
  
  if(file.exists(paste0(WD,"/Assessment Results/Report Files/TAC_est",iteration,".cor",sep=""))==TRUE)
  {
    output =readList(paste0(WD,"/Assessment Results/Report Files/TAC_est",iteration,".rep",sep=""))
    par_names=c('Obj_Func','max_gradient','SIM_h','SIM_R_ave','SIM_F','SIM_ssb','SIM_recruits','h','R_ave','fmult','ssb','recruits')
    results=output[par_names]
    c(rep=iteration,
      objective=results$Obj_Func,
      grad=results$max_gradient,
      h_est=results$h,
      R0_est=results$R_ave,
      F_est=as.vector(results$fmult),
      ssb_est=as.vector(results$ssb),
      recruits_est=as.vector(results$recruits)
    )
  }
  else
  {
    output =readList(paste0(WD,"/Assessment Results/Report Files/TAC_est",iteration,".rep",sep=""))
    par_names=c('Obj_Func','max_gradient','SIM_h','SIM_R_ave','SIM_F','SIM_ssb','SIM_recruits','h','R_ave','fmult','ssb','recruits')
    results=output[par_names]
    c(rep=iteration,
      objective=NA,
      grad=NA,
      h_est=NA,
      R0_est=NA,
      F_est=as.vector(rep(NA,times=length(results$fmult))),
      ssb_est=as.vector(rep(NA,times=length(results$ssb))),
      recruits_est=as.vector(rep(NA,times=length(results$recruits)))
    )
  }
}

true_value=function(WD,iteration)
{
  setwd(paste0(WD,"/Assessment Results/Report Files"))
  output =readList(paste0(WD,"/Assessment Results/Report Files/TAC_est",iteration,".rep",sep=""))
  par_names=c('Obj_Func','max_gradient','SIM_h','SIM_R_ave','SIM_F','SIM_ssb','SIM_recruits','h','R_ave','fmult','ssb','recruits')
  results=output[par_names]
  c(rep=iteration,
    h_true=results$SIM_h,
    R0_true=results$SIM_R_ave,
    F_true=as.vector(results$SIM_F),
    ssb_true=as.vector(results$SIM_ssb),
    recruits_true=as.vector(results$SIM_recruits)
  )
}

calc_bias=function(WD,iteration)
{
  est=cbind(sapply(1:iteration,function(i)result(WD,i)))
  true_temp=cbind(sapply(1:iteration,function(i)true_value(WD,i)))
  
  Converged<-(iteration-length(which(is.na(est[2,]))))
  percent.converged<-(Converged/iteration)*100
  
  est1<-est[-c(1:3),]
  write.csv(est1,file=paste0(WD,"/Bias/Estimated Parameters.csv",sep=""))
  
  setwd(paste0(WD,"/Assessment Results/Report Files"))
  report.file=readList(paste0(WD,"/Assessment Results/Report Files/TAC_est1.rep",sep=""))
  nyrs<-as.numeric(report.file['nyrs_SIM'])
  nyrs_ASS<-as.numeric(report.file['nyrs_assessment'])
  
  true=true_temp[c(2:(nyrs_ASS+3),(nyrs_ASS+3+1+(nyrs-nyrs_ASS)):(nyrs_ASS+3+1+(nyrs)-1),(nyrs_ASS+3+1+(nyrs)+(nyrs-nyrs_ASS)):(nyrs_ASS+3+1+(nyrs)-1+(nyrs))),]
  write.csv(true_temp,file=paste0(WD,"/Bias/Simulated Parameters(ALL).csv",sep=""))
  write.csv(true,file=paste0(WD,"/Bias/Simulated Parameters(Assess).csv",sep=""))
  
  setwd(WD)
  bias<-matrix(NA,nrow=length(true[,1]),ncol=iteration)
  percent_bias<-matrix(NA,nrow=length(true[,1]),ncol=iteration)
  
  for(i in 1:length(true[,1]))  
  {
    for(j in 1:iteration)
    {
      if(true[i,j]==0 || is.na(est1[i,j]))
      {
        bias[i,j]=NA
        percent_bias[i,j]=NA
      }
      else{
        bias[i,j]=est1[i,j]-true[i,j]
        percent_bias[i,j]=bias[i,j]/true[i]*100
      }
    }
  }
  
  row.names(percent_bias)<-c('Steep','R0',rep('F',times=nyrs_ASS),rep('SSB',times=nyrs_ASS),rep('RECR',times=nyrs_ASS))
  row.names(bias)<-c('Steep','R0',rep('F',times=nyrs_ASS),rep('SSB',times=nyrs_ASS),rep('RECR',times=nyrs_ASS))
  write.csv(percent_bias,file=paste0(WD,"/Bias/Percent_Bias.csv",sep=""))
  write.csv(bias,file=paste0(WD,"/Bias/Bias.csv",sep=""))
  
  pdf(file=paste0(WD,"/Bias/Percent_Bias.pdf",sep=""))
  plot.beans.combined('Steepness',percent_bias[1,],Converged,percent.converged,ylab='% Bias')
  plot.beans.combined('R0',percent_bias[2,],Converged,percent.converged,ylab='% Bias')
  plot.beans.combined('F (Combined Across Years)',percent_bias[c(3:(nyrs_ASS+2)),],Converged,percent.converged,ylab='% Bias')
  plot.beans('F (By Year)',as.data.frame(t(percent_bias[c(3:(nyrs_ASS+2)),])),Converged,percent.converged,names=c(1:75),xlab='Assessment Year',ylab='% Bias')

  plot.beans.value('F (By Year)',as.data.frame(t(est1[c(3:(nyrs_ASS+2)),])),Converged,percent.converged,names=c(1:75),xlab='Assessment Year',ylab='F')
  matplot(1:75,true[c(3:(nyrs_ASS+2)),2],type='l',lwd=2,add=T)
  matpoints(1:75,true[c(3:(nyrs_ASS+2)),2],pch=21,col='black',cex=0.75,bg='white')
  
  plot.beans.combined('SSB (Combined Across Years)',percent_bias[c((nyrs_ASS+2+1):(nyrs_ASS+2+1+(nyrs_ASS)-1)),],Converged,percent.converged,ylab='% Bias')
  plot.beans('SSB (By Year)',as.data.frame(t(percent_bias[c((nyrs_ASS+2+1):(nyrs_ASS+2+1+(nyrs_ASS)-1)),])),Converged,percent.converged,names=c(1:75),xlab='Assessment Year',ylab='% Bias')

  plot.beans.value('SSB (By Year)',as.data.frame(t(est1[c((nyrs_ASS+2+1):(nyrs_ASS+2+1+(nyrs_ASS)-1)),])),Converged,percent.converged,names=c(1:75),xlab='Assessment Year',ylab='SSB')
  matplot(1:75,true[c((nyrs_ASS+2+1):(nyrs_ASS+2+1+(nyrs_ASS)-1)),2],type='l',lwd=2,add=T)
  matpoints(1:75,true[c((nyrs_ASS+2+1):(nyrs_ASS+2+1+(nyrs_ASS)-1)),2],pch=21,col='black',cex=0.75,bg='white')
  
  plot.beans.combined('Recruits (Combined Across Years)',percent_bias[c((nyrs_ASS+2+1+(nyrs_ASS)):(nyrs_ASS+2+1+(nyrs_ASS)+nyrs_ASS-1)),],Converged,percent.converged,ylab='% Bias')
  plot.beans('Recruits (By Year)',as.data.frame(t(percent_bias[c((nyrs_ASS+2+1+(nyrs_ASS)):(nyrs_ASS+2+1+(nyrs_ASS)+nyrs_ASS-1)),])),Converged,percent.converged,names=c(1:75),xlab='Assessment Year',ylab='% Bias')

  plot.beans.value('Recruits (By Year)',as.data.frame(t(est1[c((nyrs_ASS+2+1+(nyrs_ASS)):(nyrs_ASS+2+1+(nyrs_ASS)+nyrs_ASS-1)),])),Converged,percent.converged,names=c(1:75),xlab='Assessment Year',ylab='Recruits')
  matplot(1:75,true[c((nyrs_ASS+2+1+(nyrs_ASS)):(nyrs_ASS+2+1+(nyrs_ASS)+nyrs_ASS-1)),2],type='l',lwd=2,add=T)
  matpoints(1:75,true[c((nyrs_ASS+2+1+(nyrs_ASS)):(nyrs_ASS+2+1+(nyrs_ASS)+nyrs_ASS-1)),2],pch=21,col='black',cex=0.75,bg='white')
  
  dev.off()
}


Model_Outputs=function(WD,iteration)
{
 result(WD,iteration)
 true_value(WD,iteration)
 calc_bias(WD,iteration)
}




###############################################################################################################################################
###############################################################################################################################################
################ Perform Simulations and Estimation ###########################################################################################
###############################################################################################################################################
wd1<-getwd()
setwd(wd1)

if(identical(SIM.PRE.TAC.F.LEVEL,'ALL')==TRUE & identical(TAC.Constant,'BOTH')==TRUE)
{
  ######### CNST TAC, low F ###############################################################
  
  label<-'CNST TAC, LOW F'
  dir.create(paste0(wd1,"/",label,sep=""))
  
  wd<-paste0(wd1,"/",label,sep="")
  F.SCALAR<-FMSY.SCALAR.LOW
  TAC_switch<-1
  dir.create(paste0(wd,"/Simulation Results",sep=""))
  dir.create(paste0(wd,"/Bias",sep=""))
  dir.create(paste0(wd,"/Assessment Results",sep=""))
  dir.create(paste0(wd,"/Assessment Results/Report Files",sep=""))
  
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd1","wd","run.SIM","run.est","F.SCALAR","TAC_switch"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste0("Simulation Model Progress Bar (",label,")",sep=""), label=paste("Simulation Run 0 of ",nloops,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/nloops*100),label=paste("Simulation Run", n,"of", nloops,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:nloops,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
  ) %dopar% {
    run.SIM(i,wd1,wd,TAC_switch,F.SCALAR) 
    run.est(i,wd1,wd) 
  }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  Model_Outputs(wd,nloops)

  
  ######### CNST TAC, High F ###############################################################
  
  
  label<-'CNST TAC, High F'
  dir.create(paste0(wd1,"/",label,sep=""))
  
  wd<-paste0(wd1,"/",label,sep="")
  F.SCALAR<-FMSY.SCALAR.HIGH
  TAC_switch<-1
  dir.create(paste0(wd,"/Simulation Results",sep=""))
  dir.create(paste0(wd,"/Bias",sep=""))
  dir.create(paste0(wd,"/Assessment Results",sep=""))
  dir.create(paste0(wd,"/Assessment Results/Report Files",sep=""))
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd1","wd","run.SIM","run.est","F.SCALAR","TAC_switch"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste0("Simulation Model Progress Bar (",label,")",sep=""), label=paste("Simulation Run 0 of ",nloops,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/nloops*100),label=paste("Simulation Run", n,"of", nloops,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:nloops,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
  ) %dopar% {
    run.SIM(i,wd1,wd,TAC_switch,F.SCALAR) 
    run.est(i,wd1,wd) 
  }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  Model_Outputs(wd,nloops)
  
  ######### CNST TAC, Fmsy ###############################################################
  
  
  label<-'CNST TAC, Fmsy'
  dir.create(paste0(wd1,"/",label,sep=""))
  
  wd<-paste0(wd1,"/",label,sep="")
  F.SCALAR<-1.0
  TAC_switch<-1
  dir.create(paste0(wd,"/Simulation Results",sep=""))
  dir.create(paste0(wd,"/Bias",sep=""))
  dir.create(paste0(wd,"/Assessment Results",sep=""))
  dir.create(paste0(wd,"/Assessment Results/Report Files",sep=""))
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd1","wd","run.SIM","run.est","F.SCALAR","TAC_switch"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste0("Simulation Model Progress Bar (",label,")",sep=""), label=paste("Simulation Run 0 of ",nloops,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/nloops*100),label=paste("Simulation Run", n,"of", nloops,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:nloops,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
  ) %dopar% {
    run.SIM(i,wd1,wd,TAC_switch,F.SCALAR) 
    run.est(i,wd1,wd) 
  }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  Model_Outputs(wd,nloops)
  
  ######### Dynamic TAC, low F ###############################################################
  
  
  label<-'Dynamic TAC, LOW F'
  dir.create(paste0(wd1,"/",label,sep=""))
  
  wd<-paste0(wd1,"/",label,sep="")
  F.SCALAR<-FMSY.SCALAR.LOW
  TAC_switch<-0
  dir.create(paste0(wd,"/Simulation Results",sep=""))
  dir.create(paste0(wd,"/Bias",sep=""))
  dir.create(paste0(wd,"/Assessment Results",sep=""))
  dir.create(paste0(wd,"/Assessment Results/Report Files",sep=""))
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd1","wd","run.SIM","run.est","F.SCALAR","TAC_switch"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste0("Simulation Model Progress Bar (",label,")",sep=""), label=paste("Simulation Run 0 of ",nloops,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/nloops*100),label=paste("Simulation Run", n,"of", nloops,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:nloops,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
  ) %dopar% {
    run.SIM(i,wd1,wd,TAC_switch,F.SCALAR) 
    run.est(i,wd1,wd) 
  }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  Model_Outputs(wd,nloops)
  
  
  ######### Dynamic TAC, High F ###############################################################
  
  
  label<-'Dynamic TAC, High F'
  dir.create(paste0(wd1,"/",label,sep=""))
  
  wd<-paste0(wd1,"/",label,sep="")
  F.SCALAR<-FMSY.SCALAR.HIGH
  TAC_switch<-0
  dir.create(paste0(wd,"/Simulation Results",sep=""))
  dir.create(paste0(wd,"/Bias",sep=""))
  dir.create(paste0(wd,"/Assessment Results",sep=""))
  dir.create(paste0(wd,"/Assessment Results/Report Files",sep=""))
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd1","wd","run.SIM","run.est","F.SCALAR","TAC_switch"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste0("Simulation Model Progress Bar (",label,")",sep=""), label=paste("Simulation Run 0 of ",nloops,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/nloops*100),label=paste("Simulation Run", n,"of", nloops,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:nloops,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
  ) %dopar% {
    run.SIM(i,wd1,wd,TAC_switch,F.SCALAR) 
    run.est(i,wd1,wd) 
  }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  Model_Outputs(wd,nloops)
  
  ######### Dynamic TAC, Fmsy ###############################################################

  
  label<-'Dynamic TAC, Fmsy'
  dir.create(paste0(wd1,"/",label,sep=""))
  
  wd<-paste0(wd1,"/",label,sep="")
  F.SCALAR<-1.0
  TAC_switch<-0
  dir.create(paste0(wd,"/Simulation Results",sep=""))
  dir.create(paste0(wd,"/Bias",sep=""))
  dir.create(paste0(wd,"/Assessment Results",sep=""))
  dir.create(paste0(wd,"/Assessment Results/Report Files",sep=""))
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd1","wd","run.SIM","run.est","F.SCALAR","TAC_switch"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste0("Simulation Model Progress Bar (",label,")",sep=""), label=paste("Simulation Run 0 of ",nloops,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/nloops*100),label=paste("Simulation Run", n,"of", nloops,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:nloops,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
  ) %dopar% {
    run.SIM(i,wd1,wd,TAC_switch,F.SCALAR) 
    run.est(i,wd1,wd) 
  }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  Model_Outputs(wd,nloops)
}

if(identical(SIM.PRE.TAC.F.LEVEL,'ALL')==TRUE & identical(TAC.Constant,'TRUE')==TRUE)
{
  ######### CNST TAC, low F ###############################################################

  
  label<-'CNST TAC, LOW F'
  dir.create(paste0(wd1,"/",label,sep=""))
  
  wd<-paste0(wd1,"/",label,sep="")
  F.SCALAR<-FMSY.SCALAR.LOW
  TAC_switch<-1
  dir.create(paste0(wd,"/Simulation Results",sep=""))
  dir.create(paste0(wd,"/Bias",sep=""))
  dir.create(paste0(wd,"/Assessment Results",sep=""))
  dir.create(paste0(wd,"/Assessment Results/Report Files",sep=""))
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd1","wd","run.SIM","run.est","F.SCALAR","TAC_switch"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste0("Simulation Model Progress Bar (",label,")",sep=""), label=paste("Simulation Run 0 of ",nloops,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/nloops*100),label=paste("Simulation Run", n,"of", nloops,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:nloops,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
  ) %dopar% {
    run.SIM(i,wd1,wd,TAC_switch,F.SCALAR) 
    run.est(i,wd1,wd) 
  }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  Model_Outputs(wd,nloops)
  
  
  ######### CNST TAC, High F ###############################################################
  
  label<-'CNST TAC, High F'
  dir.create(paste0(wd1,"/",label,sep=""))
  
  wd<-paste0(wd1,"/",label,sep="")
  F.SCALAR<-FMSY.SCALAR.HIGH
  TAC_switch<-1
  dir.create(paste0(wd,"/Simulation Results",sep=""))
  dir.create(paste0(wd,"/Bias",sep=""))
  dir.create(paste0(wd,"/Assessment Results",sep=""))
  dir.create(paste0(wd,"/Assessment Results/Report Files",sep=""))
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd1","wd","run.SIM","run.est","F.SCALAR","TAC_switch"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste0("Simulation Model Progress Bar (",label,")",sep=""), label=paste("Simulation Run 0 of ",nloops,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/nloops*100),label=paste("Simulation Run", n,"of", nloops,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:nloops,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
  ) %dopar% {
    run.SIM(i,wd1,wd,TAC_switch,F.SCALAR) 
    run.est(i,wd1,wd) 
  }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  Model_Outputs(wd,nloops)
  
  ######### CNST TAC, Fmsy ###############################################################

  
  label<-'CNST TAC, Fmsy'
  dir.create(paste0(wd1,"/",label,sep=""))
  
  wd<-paste0(wd1,"/",label,sep="")
  F.SCALAR<-1.0
  TAC_switch<-1
  dir.create(paste0(wd,"/Simulation Results",sep=""))
  dir.create(paste0(wd,"/Bias",sep=""))
  dir.create(paste0(wd,"/Assessment Results",sep=""))
  dir.create(paste0(wd,"/Assessment Results/Report Files",sep=""))
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd1","wd","run.SIM","run.est","F.SCALAR","TAC_switch"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste0("Simulation Model Progress Bar (",label,")",sep=""), label=paste("Simulation Run 0 of ",nloops,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/nloops*100),label=paste("Simulation Run", n,"of", nloops,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:nloops,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
  ) %dopar% {
    run.SIM(i,wd1,wd,TAC_switch,F.SCALAR) 
    run.est(i,wd1,wd) 
  }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  Model_Outputs(wd,nloops)

}


if(identical(SIM.PRE.TAC.F.LEVEL,'ALL')==TRUE & identical(TAC.Constant,'FALSE')==TRUE)
{

  ######### Dynamic TAC, low F ###############################################################
  
  label<-'Dynamic TAC, LOW F'
  dir.create(paste0(wd1,"/",label,sep=""))
  
  wd<-paste0(wd1,"/",label,sep="")
  F.SCALAR<-FMSY.SCALAR.LOW
  TAC_switch<-0
  dir.create(paste0(wd,"/Simulation Results",sep=""))
  dir.create(paste0(wd,"/Bias",sep=""))
  dir.create(paste0(wd,"/Assessment Results",sep=""))
  dir.create(paste0(wd,"/Assessment Results/Report Files",sep=""))
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd1","wd","run.SIM","run.est","F.SCALAR","TAC_switch"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste0("Simulation Model Progress Bar (",label,")",sep=""), label=paste("Simulation Run 0 of ",nloops,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/nloops*100),label=paste("Simulation Run", n,"of", nloops,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:nloops,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
  ) %dopar% {
    run.SIM(i,wd1,wd,TAC_switch,F.SCALAR) 
    run.est(i,wd1,wd) 
  }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  Model_Outputs(wd,nloops)
  
  
  ######### Dynamic TAC, High F ###############################################################
  
  label<-'Dynamic TAC, High F'
  dir.create(paste0(wd1,"/",label,sep=""))
  
  wd<-paste0(wd1,"/",label,sep="")
  F.SCALAR<-FMSY.SCALAR.HIGH
  TAC_switch<-0
  dir.create(paste0(wd,"/Simulation Results",sep=""))
  dir.create(paste0(wd,"/Bias",sep=""))
  dir.create(paste0(wd,"/Assessment Results",sep=""))
  dir.create(paste0(wd,"/Assessment Results/Report Files",sep=""))
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd1","wd","run.SIM","run.est","F.SCALAR","TAC_switch"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste0("Simulation Model Progress Bar (",label,")",sep=""), label=paste("Simulation Run 0 of ",nloops,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/nloops*100),label=paste("Simulation Run", n,"of", nloops,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:nloops,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
  ) %dopar% {
    run.SIM(i,wd1,wd,TAC_switch,F.SCALAR) 
    run.est(i,wd1,wd) 
  }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  Model_Outputs(wd,nloops)
  
  ######### Dynamic TAC, Fmsy ###############################################################
  
  label<-'Dynamic TAC, Fmsy'
  dir.create(paste0(wd1,"/",label,sep=""))
  
  wd<-paste0(wd1,"/",label,sep="")
  F.SCALAR<-1.0
  TAC_switch<-0
  dir.create(paste0(wd,"/Simulation Results",sep=""))
  dir.create(paste0(wd,"/Bias",sep=""))
  dir.create(paste0(wd,"/Assessment Results",sep=""))
  dir.create(paste0(wd,"/Assessment Results/Report Files",sep=""))
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd1","wd","run.SIM","run.est","F.SCALAR","TAC_switch"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste0("Simulation Model Progress Bar (",label,")",sep=""), label=paste("Simulation Run 0 of ",nloops,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/nloops*100),label=paste("Simulation Run", n,"of", nloops,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:nloops,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
  ) %dopar% {
    run.SIM(i,wd1,wd,TAC_switch,F.SCALAR) 
    run.est(i,wd1,wd) 
  }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  Model_Outputs(wd,nloops)
}


if(identical(SIM.PRE.TAC.F.LEVEL,'LOW')==TRUE & identical(TAC.Constant,'BOTH')==TRUE)
{
  ######### CNST TAC, low F ###############################################################
  
  label<-'CNST TAC, LOW F'
  dir.create(paste0(wd1,"/",label,sep=""))
  
  wd<-paste0(wd1,"/",label,sep="")
  F.SCALAR<-FMSY.SCALAR.LOW
  TAC_switch<-1
  dir.create(paste0(wd,"/Simulation Results",sep=""))
  dir.create(paste0(wd,"/Bias",sep=""))
  dir.create(paste0(wd,"/Assessment Results",sep=""))
  dir.create(paste0(wd,"/Assessment Results/Report Files",sep=""))
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd1","wd","run.SIM","run.est","F.SCALAR","TAC_switch"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste0("Simulation Model Progress Bar (",label,")",sep=""), label=paste("Simulation Run 0 of ",nloops,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/nloops*100),label=paste("Simulation Run", n,"of", nloops,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:nloops,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
  ) %dopar% {
    run.SIM(i,wd1,wd,TAC_switch,F.SCALAR) 
    run.est(i,wd1,wd) 
  }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  Model_Outputs(wd,nloops)
  
  ######### Dynamic TAC, low F ###############################################################
  
  label<-'Dynamic TAC, LOW F'
  dir.create(paste0(wd1,"/",label,sep=""))
  
  wd<-paste0(wd1,"/",label,sep="")
  F.SCALAR<-FMSY.SCALAR.LOW
  TAC_switch<-0
  dir.create(paste0(wd,"/Simulation Results",sep=""))
  dir.create(paste0(wd,"/Bias",sep=""))
  dir.create(paste0(wd,"/Assessment Results",sep=""))
  dir.create(paste0(wd,"/Assessment Results/Report Files",sep=""))
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd1","wd","run.SIM","run.est","F.SCALAR","TAC_switch"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste0("Simulation Model Progress Bar (",label,")",sep=""), label=paste("Simulation Run 0 of ",nloops,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/nloops*100),label=paste("Simulation Run", n,"of", nloops,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:nloops,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
  ) %dopar% {
    run.SIM(i,wd1,wd,TAC_switch,F.SCALAR) 
    run.est(i,wd1,wd) 
  }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  Model_Outputs(wd,nloops)
 
}

if(identical(SIM.PRE.TAC.F.LEVEL,'LOW')==TRUE & identical(TAC.Constant,'TRUE')==TRUE)
{
  ######### CNST TAC, low F ###############################################################
  
  label<-'CNST TAC, LOW F'
  dir.create(paste0(wd1,"/",label,sep=""))
  
  wd<-paste0(wd1,"/",label,sep="")
  F.SCALAR<-FMSY.SCALAR.LOW
  TAC_switch<-1
  dir.create(paste0(wd,"/Simulation Results",sep=""))
  dir.create(paste0(wd,"/Bias",sep=""))
  dir.create(paste0(wd,"/Assessment Results",sep=""))
  dir.create(paste0(wd,"/Assessment Results/Report Files",sep=""))
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd1","wd","run.SIM","run.est","F.SCALAR","TAC_switch"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste0("Simulation Model Progress Bar (",label,")",sep=""), label=paste("Simulation Run 0 of ",nloops,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/nloops*100),label=paste("Simulation Run", n,"of", nloops,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:nloops,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
  ) %dopar% {
    run.SIM(i,wd1,wd,TAC_switch,F.SCALAR) 
    run.est(i,wd1,wd) 
  }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  Model_Outputs(wd,nloops)

}

if(identical(SIM.PRE.TAC.F.LEVEL,'LOW')==TRUE & identical(TAC.Constant,'FALSE')==TRUE)
{
 
  ######### Dynamic TAC, low F ###############################################################
  
  label<-'Dynamic TAC, LOW F'
  dir.create(paste0(wd1,"/",label,sep=""))
  
  wd<-paste0(wd1,"/",label,sep="")
  F.SCALAR<-FMSY.SCALAR.LOW
  TAC_switch<-0
  dir.create(paste0(wd,"/Simulation Results",sep=""))
  dir.create(paste0(wd,"/Bias",sep=""))
  dir.create(paste0(wd,"/Assessment Results",sep=""))
  dir.create(paste0(wd,"/Assessment Results/Report Files",sep=""))
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd1","wd","run.SIM","run.est","F.SCALAR","TAC_switch"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste0("Simulation Model Progress Bar (",label,")",sep=""), label=paste("Simulation Run 0 of ",nloops,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/nloops*100),label=paste("Simulation Run", n,"of", nloops,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:nloops,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
  ) %dopar% {
    run.SIM(i,wd1,wd,TAC_switch,F.SCALAR) 
    run.est(i,wd1,wd) 
  }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  Model_Outputs(wd,nloops)
 
}

if(identical(SIM.PRE.TAC.F.LEVEL,'HIGH')==TRUE & identical(TAC.Constant,'BOTH')==TRUE)
{

  ######### CNST TAC, High F ###############################################################
  
  label<-'CNST TAC, High F'
  dir.create(paste0(wd1,"/",label,sep=""))
  
  wd<-paste0(wd1,"/",label,sep="")
  F.SCALAR<-FMSY.SCALAR.HIGH
  TAC_switch<-1
  dir.create(paste0(wd,"/Simulation Results",sep=""))
  dir.create(paste0(wd,"/Bias",sep=""))
  dir.create(paste0(wd,"/Assessment Results",sep=""))
  dir.create(paste0(wd,"/Assessment Results/Report Files",sep=""))
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd1","wd","run.SIM","run.est","F.SCALAR","TAC_switch"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste0("Simulation Model Progress Bar (",label,")",sep=""), label=paste("Simulation Run 0 of ",nloops,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/nloops*100),label=paste("Simulation Run", n,"of", nloops,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:nloops,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
  ) %dopar% {
    run.SIM(i,wd1,wd,TAC_switch,F.SCALAR) 
    run.est(i,wd1,wd) 
  }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  Model_Outputs(wd,nloops)
  
  ######### Dynamic TAC, High F ###############################################################
  
  label<-'Dynamic TAC, High F'
  dir.create(paste0(wd1,"/",label,sep=""))
  
  wd<-paste0(wd1,"/",label,sep="")
  F.SCALAR<-FMSY.SCALAR.HIGH
  TAC_switch<-0
  dir.create(paste0(wd,"/Simulation Results",sep=""))
  dir.create(paste0(wd,"/Bias",sep=""))
  dir.create(paste0(wd,"/Assessment Results",sep=""))
  dir.create(paste0(wd,"/Assessment Results/Report Files",sep=""))
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd1","wd","run.SIM","run.est","F.SCALAR","TAC_switch"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste0("Simulation Model Progress Bar (",label,")",sep=""), label=paste("Simulation Run 0 of ",nloops,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/nloops*100),label=paste("Simulation Run", n,"of", nloops,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:nloops,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
  ) %dopar% {
    run.SIM(i,wd1,wd,TAC_switch,F.SCALAR) 
    run.est(i,wd1,wd) 
  }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  Model_Outputs(wd,nloops)
  
}

if(identical(SIM.PRE.TAC.F.LEVEL,'HIGH')==TRUE & identical(TAC.Constant,'TRUE')==TRUE)
{

  ######### CNST TAC, High F ###############################################################
  
  label<-'CNST TAC, High F'
  dir.create(paste0(wd1,"/",label,sep=""))
  
  wd<-paste0(wd1,"/",label,sep="")
  F.SCALAR<-FMSY.SCALAR.HIGH
  TAC_switch<-1
  dir.create(paste0(wd,"/Simulation Results",sep=""))
  dir.create(paste0(wd,"/Bias",sep=""))
  dir.create(paste0(wd,"/Assessment Results",sep=""))
  dir.create(paste0(wd,"/Assessment Results/Report Files",sep=""))
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd1","wd","run.SIM","run.est","F.SCALAR","TAC_switch"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste0("Simulation Model Progress Bar (",label,")",sep=""), label=paste("Simulation Run 0 of ",nloops,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/nloops*100),label=paste("Simulation Run", n,"of", nloops,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:nloops,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
  ) %dopar% {
    run.SIM(i,wd1,wd,TAC_switch,F.SCALAR) 
    run.est(i,wd1,wd) 
  }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  Model_Outputs(wd,nloops)
  
}


if(identical(SIM.PRE.TAC.F.LEVEL,'HIGH')==TRUE & identical(TAC.Constant,'FALSE')==TRUE)
{
 
  ######### Dynamic TAC, High F ###############################################################
  
  label<-'Dynamic TAC, High F'
  dir.create(paste0(wd1,"/",label,sep=""))
  
  wd<-paste0(wd1,"/",label,sep="")
  F.SCALAR<-FMSY.SCALAR.HIGH
  TAC_switch<-0
  dir.create(paste0(wd,"/Simulation Results",sep=""))
  dir.create(paste0(wd,"/Bias",sep=""))
  dir.create(paste0(wd,"/Assessment Results",sep=""))
  dir.create(paste0(wd,"/Assessment Results/Report Files",sep=""))
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd1","wd","run.SIM","run.est","F.SCALAR","TAC_switch"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste0("Simulation Model Progress Bar (",label,")",sep=""), label=paste("Simulation Run 0 of ",nloops,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/nloops*100),label=paste("Simulation Run", n,"of", nloops,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:nloops,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
  ) %dopar% {
    run.SIM(i,wd1,wd,TAC_switch,F.SCALAR) 
    run.est(i,wd1,wd) 
  }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  Model_Outputs(wd,nloops)

}

if(identical(SIM.PRE.TAC.F.LEVEL,'FMSY')==TRUE & identical(TAC.Constant,'BOTH')==TRUE)
{

  ######### CNST TAC, Fmsy ###############################################################
  
  label<-'CNST TAC, Fmsy'
  dir.create(paste0(wd1,"/",label,sep=""))
  
  wd<-paste0(wd1,"/",label,sep="")
  F.SCALAR<-1.0
  TAC_switch<-1
  dir.create(paste0(wd,"/Simulation Results",sep=""))
  dir.create(paste0(wd,"/Bias",sep=""))
  dir.create(paste0(wd,"/Assessment Results",sep=""))
  dir.create(paste0(wd,"/Assessment Results/Report Files",sep=""))
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd1","wd","run.SIM","run.est","F.SCALAR","TAC_switch"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste0("Simulation Model Progress Bar (",label,")",sep=""), label=paste("Simulation Run 0 of ",nloops,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/nloops*100),label=paste("Simulation Run", n,"of", nloops,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:nloops,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
  ) %dopar% {
    run.SIM(i,wd1,wd,TAC_switch,F.SCALAR) 
    run.est(i,wd1,wd) 
  }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  Model_Outputs(wd,nloops)
  
  ######### Dynamic TAC, Fmsy ###############################################################
  
  label<-'Dynamic TAC, Fmsy'
  dir.create(paste0(wd1,"/",label,sep=""))
  
  wd<-paste0(wd1,"/",label,sep="")
  F.SCALAR<-1.0
  TAC_switch<-0
  dir.create(paste0(wd,"/Simulation Results",sep=""))
  dir.create(paste0(wd,"/Bias",sep=""))
  dir.create(paste0(wd,"/Assessment Results",sep=""))
  dir.create(paste0(wd,"/Assessment Results/Report Files",sep=""))
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd1","wd","run.SIM","run.est","F.SCALAR","TAC_switch"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste0("Simulation Model Progress Bar (",label,")",sep=""), label=paste("Simulation Run 0 of ",nloops,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/nloops*100),label=paste("Simulation Run", n,"of", nloops,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:nloops,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
  ) %dopar% {
    run.SIM(i,wd1,wd,TAC_switch,F.SCALAR) 
    run.est(i,wd1,wd) 
  }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  Model_Outputs(wd,nloops)
}

if(identical(SIM.PRE.TAC.F.LEVEL,'FMSY')==TRUE & identical(TAC.Constant,'TRUE')==TRUE)
{
 
  ######### CNST TAC, Fmsy ###############################################################
  
  label<-'CNST TAC, Fmsy'
  dir.create(paste0(wd1,"/",label,sep=""))
  
  wd<-paste0(wd1,"/",label,sep="")
  F.SCALAR<-1.0
  TAC_switch<-1
  dir.create(paste0(wd,"/Simulation Results",sep=""))
  dir.create(paste0(wd,"/Bias",sep=""))
  dir.create(paste0(wd,"/Assessment Results",sep=""))
  dir.create(paste0(wd,"/Assessment Results/Report Files",sep=""))
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd1","wd","run.SIM","run.est","F.SCALAR","TAC_switch"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste0("Simulation Model Progress Bar (",label,")",sep=""), label=paste("Simulation Run 0 of ",nloops,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/nloops*100),label=paste("Simulation Run", n,"of", nloops,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:nloops,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
  ) %dopar% {
    run.SIM(i,wd1,wd,TAC_switch,F.SCALAR) 
    run.est(i,wd1,wd) 
  }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  Model_Outputs(wd,nloops)
  
}

if(identical(SIM.PRE.TAC.F.LEVEL,'FMSY')==TRUE & identical(TAC.Constant,'FALSE')==TRUE)
{
 
  ######### Dynamic TAC, Fmsy ###############################################################
  
  label<-'Dynamic TAC, Fmsy'
  dir.create(paste0(wd1,"/",label,sep=""))
  
  wd<-paste0(wd1,"/",label,sep="")
  F.SCALAR<-1.0
  TAC_switch<-0
  dir.create(paste0(wd,"/Simulation Results",sep=""))
  dir.create(paste0(wd,"/Bias",sep=""))
  dir.create(paste0(wd,"/Assessment Results",sep=""))
  dir.create(paste0(wd,"/Assessment Results/Report Files",sep=""))
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd1","wd","run.SIM","run.est","F.SCALAR","TAC_switch"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste0("Simulation Model Progress Bar (",label,")",sep=""), label=paste("Simulation Run 0 of ",nloops,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/nloops*100),label=paste("Simulation Run", n,"of", nloops,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:nloops,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
  ) %dopar% {
    run.SIM(i,wd1,wd,TAC_switch,F.SCALAR) 
    run.est(i,wd1,wd) 
  }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  Model_Outputs(wd,nloops)
}







