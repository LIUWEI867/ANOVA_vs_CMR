setwd('C:/Users/Liuwe/Desktop/Sta_vs_Anova/STACMR-R')
source("STAvsANOVAfunctions.R")
source("staCMRsetup.R")
resultsDf=data.frame(niter=0,#number of iterations
                     nANOVA=0,#anova reject or not
                     nANOVA_AB=0,#anova reject or not for AB inter only
                     nboot=0,#bootanova reject or not for AB inter only
                     p_AinterB=0,#pvalue of AB interaction in ANOVA
                     p_AinterC=0,#pvalue of AC interaction in ANOVA
                     p_AinterBinterC=0,#pvalue of ABC interaction in ANOVA
                     p_boot_AB=0,#pvalue of AB interaction in boot-ANOVA
                     nSTA_NP=0,#STA without partial order reject or not
                     nSTA_LP=0,#STA with looser partial order reject or not
                     nSTA=0,#STA with partial order reject or not
                     p_STA_NP=0,#pvalue of STA without partial order
                     p_STA_LP=0,#pvalue of STA with looser partial order
                     p_STA=0,#pvalue of STA with partial order
                     fit_STA_NP=0,fit_STA_LP=0,fit_STA=0,
                     MR_fit_NP=0,CMR_fit_NP=0,MR_fit_LP=0,CMR_fit_LP=0,MR_fit=0,CMR_fit=0,
                     nSTAMR=0,p_STAMR=0,fit_STAMR=0
)
se=4
n=30
modelType=c("Linear","NonLinear1LV","NonLinear2LV")[1]#select a model
filename=paste0("results/" ,modelType,"_SE=",se,"_n=",n,".csv")
for (iter in 1:1000) {
  Thisdata=genLinData(se=se,n=n,a=c(-1,1),b=c(-0.5,0.5),c=c(-2,-1,0,1,2),modelType=modelType)
  anovaResults=ResultsofAnova(df=Thisdata[[1]])
  STARestults=ResutlsofSTA(staDF=Thisdata[[2]])
  resultsDf[iter,]=c(iter,anovaResults,STARestults)
  nalpha_ANOVA=sum(resultsDf[3])
  nalpha_ANOVAboot=sum(resultsDf[4])
  nalpha_STA_NP=sum(resultsDf[9])
  nalpha_STA_LP=sum(resultsDf[10])
  nalpha_STA=sum(resultsDf[11])
  cat('iteration:',iter,'\n','number of rejection for ANOVA:',nalpha_ANOVA,'\n',
      'number of rejection for bootANOVA:',nalpha_ANOVAboot,'\n',
      'number of rejection for STA without partial order:',nalpha_STA_NP,'\n',
      'number of rejection for STA with looser partial order:',nalpha_STA_NP,'\n',
      'number of rejection for STA with partial order:',nalpha_STA,'\n'
      )
  if (iter %% 20==0) {
    write.csv(resultsDf,file = filename)
  }
}


