
rm(list=(ls()))
wd<-getwd()
setwd(wd)

############### INPUTS ##############################################################################################################################
################ ALL FLAGS MUST BE IN ALL ****CAPS********########################################################
#### FMSY Sim Inputs ####################################
Perform.FMSY.SIM<-'FALSE'                                     # 'TRUE' or 'FALSE'
FMSY.SIM.Recruit.Type<-'DETERMINISTIC'                       # 'DETERMINISTIC' or 'STOCHASTIC'
F.start=0.25                                                 # Starting F value that will search from for Fmsy
F.end=.35                                               # Terminal F value that will search to for Fmsy 
interval=0.001                                               # Step increase in Fmsy search 
#######################################################################

nloops <-1                                                 # How many simulations to perform

############# Model Flags/Switches ###################################################################
################ ALL FLAGS MUST BE IN ALL ****CAPS********########################################################

SIM.Recruit.Type<-'STOCHASTIC'                               # 'DETERMINISTIC' or 'STOCHASTIC' recruitment for simulations

SIM.Survey.Sel.Type<-'LOGISTIC'                              # 'LOGISTIC' or 'BY_AGE'                      
Survey.Sel.Type<-'LOGISTIC'
SIM.Fishery.Sel.Type<-'LOGISTIC'
Fishery.Sel.Type<-'LOGISTIC'

Init.Abund.Est.Type<-'ESTIMATE'                     # 'EQUIL_OFFSET' or 'ESTIMATE', if equil then assume pop in equil with no fishing based on an initial Recruitment offset from R0 by R0_offset, if estimate then estimate abund at age

##################### Types of SIM Model Runs ########################################################################################################

PRE.TAC.F.SCALAR<-1

TAC.CNST<-0                                         # 1='TRUE' 0='FALSE' if true TAC is held constant at level in first year, if false TAC is reset each year based on estimated abundance and FMSY, if both computes results for both

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


    SIM.DATA =readLines(paste0(wd,"/SIM_TAC.dat",sep=""),n=-1)
    
    if(TAC.CNST==1)
    {
      SIM.DATA[(grep("TAC_CNST_switch",SIM.DATA)+1)]=1
    }
    if(TAC.CNST==0)
    {
      SIM.DATA[(grep("TAC_CNST_switch",SIM.DATA)+1)]=0
    }    
    SIM.DATA[(grep("Fmsy_scalar",SIM.DATA)+1)]=PRE.TAC.F.SCALAR
    writeLines(SIM.DATA,paste0(wd,"/SIM_TAC.dat",sep=""))
  
  setwd(wd)
  
  SIM.DAT=readLines("SIM_TAC.dat",n=-1)
  SIM.DAT[(grep("myseed",SIM.DAT)+1)]=1
  writeLines(SIM.DAT,"SIM_TAC.dat")
  
  system(paste0("SIM_TAC -nox -nohess SIM_TAC.dat",sep=""),wait=TRUE,show.output.on.console=FALSE)  #-nox keeps from showing the ADMB gradient outputs on screen
  
  invisible(file.copy(from=paste0(wd,"/SIM_TAC.rep",sep=""),to=paste0(wd,"/TAC_est.dat",sep="")))
  
  invisible(shell(paste("del SIM_TAC.par", sep="")))
  invisible(shell(paste("del SIM_TAC.rep", sep="")))
  invisible(shell(paste("del SIM_TAC.bar",sep="")))
  invisible(shell(paste("del SIM_TAC.log",sep="")))
  invisible(shell("del fmin.log"))
  invisible(shell("del variance"))
  
  setwd(wd)
  
  system(paste0("TAC_est -nox -ind TAC_est.dat",sep=""),wait=TRUE,show.output.on.console=FALSE)  #-nox keeps from showing the ADMB gradient outputs on screen

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
  

  output =readList(paste0(wd,"/TAC_est.rep",sep=""))
  par_names=c('Obj_Func','max_gradient','SIM_h','SIM_R_ave','SIM_F','SIM_ssb','SIM_recruits','h','R_ave','fmult','ssb','recruits')
  results=output[par_names]
  EST<-c(
    h_est=results$h,
    R0_est=results$R_ave,
    F_est=as.vector(results$fmult),
    ssb_est=as.vector(results$ssb),
    recruits_est=as.vector(results$recruits)
  )



true_values =readList(paste0(wd,"/TAC_est.rep",sep=""))
par_names=c('Obj_Func','max_gradient','SIM_h','SIM_R_ave','SIM_F','SIM_ssb','SIM_recruits','h','R_ave','fmult','ssb','recruits')
results1=true_values[par_names]
true<-c(
  h_true=results1$SIM_h,
  R0_true=results1$SIM_R_ave,
  F_true=as.vector(results1$SIM_F),
  ssb_true=as.vector(results1$SIM_ssb),
  recruits_true=as.vector(results1$SIM_recruits))



report.file=readList(paste0(wd,"/TAC_est.rep",sep=""))

nyrs<-as.numeric(report.file['nyrs_SIM'])
nyrs_ASS<-as.numeric(report.file['nyrs_assessment'])

true=true[c(1:(nyrs_ASS+2),(nyrs_ASS+2+1+(nyrs-nyrs_ASS)):(nyrs_ASS+2+1+(nyrs)-1),(nyrs_ASS+2+1+(nyrs)+(nyrs-nyrs_ASS)):(nyrs_ASS+2+1+(nyrs)-1+(nyrs)))]


pdf("Consistency Run Results.pdf")

matplot(1:75,true[c(3:(nyrs_ASS+2))],type='none',lwd=2, ylab="F", xlab='Year')
matpoints(1:75,true[c(3:(nyrs_ASS+2))],pch=16,col='black',cex=0.75)
matlines(1:75,EST[c(3:(nyrs_ASS+2))],lty=1,col='black')
abline(h=0,lty=1,lwd=2)
legend("top",c('True Value','Pred Value'),lty=c(-1,1),pch=c(16,-1),bg='white',ncol=2,cex=0.6)


matplot(1:75,true[c((nyrs_ASS+2+1):(nyrs_ASS+2+1+(nyrs_ASS)-1))],type='n',lwd=2, ylab="SSB", xlab='Year')
matpoints(1:75,true[c((nyrs_ASS+2+1):(nyrs_ASS+2+1+(nyrs_ASS)-1))],pch=16,col='black',cex=0.75)
matlines(1:75,EST[c((nyrs_ASS+2+1):(nyrs_ASS+2+1+(nyrs_ASS)-1))],lty=1,col='black')
abline(h=0,lty=1,lwd=2)
legend("top",c('True Value','Pred Value'),lty=c(-1,1),pch=c(16,-1),bg='white',ncol=2,cex=0.6)


matplot(1:75,true[c((nyrs_ASS+2+1+(nyrs_ASS)):(nyrs_ASS+2+1+(nyrs_ASS)+nyrs_ASS-1))],type='none',lwd=2, ylab="Recruits", xlab='Year')
matpoints(1:75,true[c((nyrs_ASS+2+1+(nyrs_ASS)):(nyrs_ASS+2+1+(nyrs_ASS)+nyrs_ASS-1))],pch=16,col='black',cex=0.75)
matlines(1:75,EST[c((nyrs_ASS+2+1+(nyrs_ASS)):(nyrs_ASS+2+1+(nyrs_ASS)+nyrs_ASS-1))],lty=1,col='black')
abline(h=0,lty=1,lwd=2)
legend("top",c('True Value','Pred Value'),lty=c(-1,1),pch=c(16,-1),bg='white',ncol=2,cex=0.6)



dev.off()