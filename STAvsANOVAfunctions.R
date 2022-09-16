genLinData=function(se,n,mu=5,a=c(-1,1),b=c(-0.5,0.5),c=c(-2,-1,0,1,2),#function for generating data
                    ab=array(c(0,0,0,0),dim=c(2,2)),modelType="Linear"){
  na=n*length(b)*length(c)#samplesize for each level of A
  nb=n*length(c)#samplesize for each level of B at each level of A
  m=length(a)*length(b)*length(c)#number of conditions
  idx=1
  error_R=rnorm((m*n),mean = 0,sd=se)#random error
  Y=NULL#to store simulated values
  for (i in 1:length(a)) {#simulated values
    for (j in 1:length(b)) {
      for (k in 1:length(c)) {
        for (l in 1:n) {
          if (modelType=="Linear") {Y[idx]=mu+a[i]+b[j]+c[k]+ab[i,j]+error_R[idx]}
          if (modelType=="NonLinear1LV") {
            LV=b[j]+c[k]#the single latent variable
            if(i==1)Y[idx]=LV^3+error_R[idx]#observed values in a1 level
            if(i==2)Y[idx]=2^LV+error_R[idx]#observed values in a2 level
          }
          if (modelType=="NonLinear2LV") {
            LV1=b[j]+c[k]#the first latent variable
            LV2=b[j]-c[k]#the second latent variable
            if(i==1)Y[idx]=a[i]+LV1^3+error_R[idx]#observed values in a1 level
            if(i==2)Y[idx]=a[i]+2^LV2+error_R[idx]#observed values in a2 level
          }
          idx=idx+1
        }
      }
    }
  }
  anovaDF=data.frame(A=c(rep(1,na),rep(2,na)),#dataset for anova
                     B=c(rep(1,nb),rep(2,nb),rep(1,nb),rep(2,nb)),
                     C=rep(c(rep(1,n),rep(2,n),rep(3,n),rep(4,n),rep(5,n)),4),
                     anovaY=Y)
  staDF=data.frame(SNo=c(1:(n*m)),#dataset for STA. SNO represent number of subjects, not used for analysis
                   Cno=rep(rep(c(1:10),each=n),2),#10conbinations of B and C at each level of A. 
                                                   #1-5 represent C1-C5 at B1；6-10 represent C1-C5 at B2
                   DV=rep(c(1,2),each=na),
                   staY=Y)
  return(list(anovaDF,staDF))
}

ResultsofAnova=function(df,nboot=10000){#ANOVA interaction bootstrap 
  anovaResults=aov(anovaY~factor(A)*factor(B)*factor(C),df)
  Anova=summary(anovaResults)
  if((Anova[[1]]$`Pr(>F)`[4]<=0.05)||
     (Anova[[1]]$`Pr(>F)`[5]<=0.05)||
     (Anova[[1]]$`Pr(>F)`[7]<=0.05)){nANOVA=1} else {nANOVA=0}#at least one interaction involving A is significant
  if (Anova[[1]]$`Pr(>F)`[4]<=0.05) {nANOVA_AB=1} else {nANOVA_AB=0}
  CF=Anova[[1]]$`F value`[4]#A by B interaction F value as bootstrap critical F
  anoP=Anova[[1]]$`Pr(>F)`[4]#p from anova
  agMean=aggregate(df$anovaY, 
                   by=list(df$A,df$B),FUN = "mean")#means of each A,B combinations
  Delta=agMean$x[4]-agMean$x[3]-agMean$x[2]+agMean$x[1]#different of interaction
  null11=df$anovaY[which(df$A==1 & df$B==1)]-Delta/4#Null hypothesis AB11
  null12=df$anovaY[which(df$A==1 & df$B==2)]+Delta/4#Null hypothesis AB12
  null21=df$anovaY[which(df$A==2 & df$B==1)]+Delta/4#Null hypothesis AB21
  null22=df$anovaY[which(df$A==2 & df$B==2)]-Delta/4#Null hypothesis AB22
  dfresample=df#resample dataframe to be completed
  bootFvalue=c()#to store the resampled F value of AB interaction
  for (i in 1:nboot) {
    resample11=sample(null11,length(null11),replace = TRUE)#resample from AB11
    resample12=sample(null12,length(null12),replace = TRUE)#resample from AB12
    resample21=sample(null21,length(null21),replace = TRUE)#resample from AB21
    resample22=sample(null22,length(null22),replace = TRUE)#resample from AB22
    dfresample$anovaY=c(resample11,resample12,resample21,resample22)#resample dataframe completed
    resampleResult=aov(anovaY~factor(A)*factor(B)*factor(C),dfresample)
    sumTR=summary(resampleResult)
    bootFvalue[i]=sumTR[[1]]$`F value`[4]#The F value of resampled data
  }
  bootp=sum(bootFvalue>CF)/length(bootFvalue)#bootstrap p value
  if (bootp<=0.05) {nboot=1}else{nboot=0}
  return(c(nANOVA,nANOVA_AB,nboot,
           Anova[[1]]$`Pr(>F)`[4],Anova[[1]]$`Pr(>F)`[5],Anova[[1]]$`Pr(>F)`[7],
           bootp))
}

ResutlsofSTA=function(staDF,
                      partialOrderS=list(c(1:5),c(6:10),c(1,6),c(2,7),c(3,8),c(4,9),c(5,10)),
                      partialOrderL=list(c(1:5),c(6:10)),
                      nsample=10000){
  MR_NP_results=staMR(staDF,partial=list())#fit of MR
  CMR_NP_results=staCMR(staDF,partial=list())#fit of CMR
  staResult_NP=staCMRFIT(staDF,partial=list(),nsample)#no partial order
  if(staResult_NP$p<=0.05){nSTA_NP=1}else{nSTA_NP=0}
  MR_LP_results=staMR(staDF,partial=partialOrderL)#fit of MR_LP
  CMR_LP_results=staCMR(staDF,partial=partialOrderL)#fit of CMR_LP
  staResult_LP=staCMRFIT(staDF,partial=partialOrderL,nsample)#loose partial order constraint
  if(staResult_LP$p<=0.05){nSTA_LP=1}else{nSTA_LP=0}
  MR_results=staMR(staDF,partial=partialOrderS)#fit of MR
  CMR_results=staCMR(staDF,partial=partialOrderS)#fit of CMR
  staResultMR=staMRFIT(staDF,partial=partialOrderS,nsample)
  if(staResultMR$p<=0.05){nSTAMR=1}else{nSTAMR=0}
  staResult=staCMRFIT(staDF,partial = partialOrderS,nsample)#strict partial order constraint
  if(staResult$p<=0.05){nSTA=1}else{nSTA=0}
  return(c(nSTA_NP,nSTA_LP,nSTA,staResult_NP$p,staResult_LP$p,staResult$p,
           staResult_NP$datafit,staResult_LP$datafit,staResult$datafit,
           MR_NP_results$fval,CMR_NP_results$fval,
           MR_LP_results$fval,CMR_LP_results$fval,
           MR_results$fval,CMR_results$fval,
           nSTAMR,staResultMR$p,staResultMR$datafit))
}
