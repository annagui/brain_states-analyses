

#### R script for ERP analyses
### Gui et al., Diminished engagement of attentive brain states predicts later Autism Spectrum Disorder
## 06 March 2020




library(erp.easy)
library(ggplot2)
library(pastecs)
library(dplyr)
library(tidyr)
library(plyr)
library(moments)
library(Rmisc)
library(car)
library(nlme)
library(DescTools)
library(ggpubr)




#load phenotypic info

phenInfo<-read.csv("brain_states_dataset.csv", colClasses = 'character', header=T, stringsAsFactors = F)
#select phenotypic information: Outcome groups (), age in days, sex, Mullen Early Learning Composite at 8 months, VABS Socialization Score at 3 years, VABS Motor Skills Score at 3 years
phenInfo<-subset(phenInfo, select=c(ID, Outcome, phase, CA_8mo_days, sex, M_ELC_SS_1, V_Soc_SS_4, V_Mot_SS_4)) 

SummaryTable_phase <- ddply(phenInfo, "phase", summarise,  N= length(unique(ID)))
SummaryTable_phase
SummaryTable_sex_out<-ddply(phenInfo, c("Outcome", "sex"), summarise,  N= length(unique(ID)))
SummaryTable_sex_out
SummaryTable_out<-ddply(phenInfo, "Outcome", summarise,  N= length(unique(ID)))
SummaryTable_out




###CHECK DIFFERENCE IN NUMBER OF TRIALS

nTrials<-read.csv("NrTrials.csv", header=T, stringsAsFactors=F)

ddply(nTrials, "Outcome", summarise,  N= length(unique(ID)), mean=mean(NvalidTrialsFD), sd=sd(NvalidTrialsFD))
summary(aov(NvalidTrialsFD~Outcome, data=nTrials))
etaSquared(aov(NvalidTrialsFD~Outcome, data=nTrials))

ddply(nTrials, "Outcome", summarise,  N= length(unique(ID)), mean=mean(NvalidTrialsNoise), sd=sd(NvalidTrialsNoise))
summary(aov(NvalidTrialsNoise~Outcome, data=nTrials))
etaSquared(aov(NvalidTrialsNoise~Outcome, data=nTrials))

ddply(nTrials, "Outcome", summarise,  N= length(unique(ID)), mean=mean(NvalidTrialsFA), sd=sd(NvalidTrialsFA))
summary(aov(NvalidTrialsFA~Outcome, data=nTrials))
etaSquared(aov(NvalidTrialsFA~Outcome, data=nTrials))





############################################## prepare datasets ###############################################

#load the ERP data
FaceDirect<-load.data(path = "~/ERPs_folder/", condition = "FaceDirect", 
                      num.subs= 131, epoch.st=-200, epoch.end=794, header=F)
FaceAverted<-load.data(path = "~/ERPs_folder/", condition = "FaceAverted", 
                       num.subs= 131, epoch.st=-200, epoch.end=794, header=F)
Noise<-load.data(path = "~/ERPs_folder/", condition = "Noise", 
                 num.subs= 131, epoch.st=-200, epoch.end=794, header=F)
combData<-rbind.data.frame(FaceDirect,Noise) #combine the data 





##### extract ERPs from selected channels as in Webb et al., 2011 




### dataset for 3-way ANOVAs (Relation to categorical outcome of ASD)


## Mean Amplitude 

left<-combData[,c("Subject","Stimulus", "Time","V34","V35","V27","V28","V23","V24","V29","V20")]
mid<-combData[,c("Subject","Stimulus", "Time","V18","V16","V10","V19","V11","V4","V12","V5")]
right<-combData[,c("Subject","Stimulus", "Time","V3","V123","V124","V117","V116","V110","V118","V111")]

left$ELmean<-rowSums(combData[,c("V34","V35","V27","V28","V23","V24","V29","V20")])/8 #mean of amplitudes left region
left_Nc<-left[left$Time>299,] #only keep timepoints from 300 ms from the stimulus onset
left_Nc_meanAmp<-ddply(left_Nc[,c(1,2,12)],.(Subject,Stimulus),summarize,MeanAmp=mean(ELmean)) #compute mean amplitude per subject by stimulus
left_Nc_meanAmp$region<-rep("left", nrow(left_Nc_meanAmp))

mid$ELmean<-rowSums(combData[,c("V18","V16","V10","V19","V11","V4","V12","V5")])/8 #mean of amplitudes central region
mid_Nc<-mid[mid$Time>299,] #only keep timepoints from 300 ms from the stimulus onset
mid_Nc_meanAmp<-ddply(mid_Nc[,c(1,2,12)],.(Subject,Stimulus),summarize,MeanAmp=mean(ELmean)) #compute mean amplitude per subject by stimulus
mid_Nc_meanAmp$region<-rep("mid", nrow(mid_Nc_meanAmp))

right$ELmean<-rowSums(combData[,c("V3","V123","V124","V117","V116","V110","V118","V111")])/8 #mean of amplitudes right region
right_Nc<-right[right$Time>299,] #only keep timepoints from 300 ms from the stimulus onset
right_Nc_meanAmp<-ddply(right_Nc[,c(1,2,12)],.(Subject,Stimulus),summarize,MeanAmp=mean(ELmean)) #compute mean amplitude per subject by stimulus
right_Nc_meanAmp$region<-rep("right", nrow(right_Nc_meanAmp))

meanAmp_left_mid<-rbind.data.frame(left_Nc_meanAmp,mid_Nc_meanAmp)
meanAmp_all_L<-rbind.data.frame(meanAmp_left_mid,right_Nc_meanAmp)


## Peak Latency

NcPeaks_left<-p.measures(combData, electrodes=c("V34","V35","V27","V28","V23","V24","V29","V20"),window=c(300,794),num.pts = 0, pol="neg")
NcPeaks_mid<-p.measures(combData, electrodes=c("V18","V16","V10","V19","V11","V4","V12","V5"),window=c(300,794),num.pts = 0, pol="neg")
NcPeaks_right<-p.measures(combData, electrodes=c("V3","V123","V124","V117","V116","V110","V118","V111"),window=c(300,794),num.pts = 0, pol="neg")


#compute latency difference between face and noise for left and right regions
LatDiff_left_FD<-NcPeaks_left[NcPeaks_left$`Trial Type`=='FaceDirect',] 
LatDiff_left_FD$`Trial Type`<-droplevels(LatDiff_left_FD$`Trial Type`)
LatDiff_left_Noise<-NcPeaks_left[NcPeaks_left$`Trial Type`=='Noise',]
LatDiff_left_Noise$`Trial Type`<-droplevels(LatDiff_left_Noise$`Trial Type`)
LatDiff_left_FvsN<-merge(LatDiff_left_FD,LatDiff_left_Noise, by.x=c('Subject'), by.y=c('Subject'))
LatDiff_left_FvsN$`Trial Type.x`<-NULL
LatDiff_left_FvsN$`Trial Type.y`<-NULL
LatDiff_left_FvsN$`Peak Latency.x`->LatDiff_left_FvsN$PeakLat_FD
LatDiff_left_FvsN$`Peak Latency.x`<-NULL
LatDiff_left_FvsN$`Peak Amplitude.x`->LatDiff_left_FvsN$PeakAmp_FD
LatDiff_left_FvsN$`Peak Amplitude.x`<-NULL
LatDiff_left_FvsN$`Peak Latency.y`->LatDiff_left_FvsN$PeakLat_Noise
LatDiff_left_FvsN$`Peak Latency.y`<-NULL
LatDiff_left_FvsN$`Peak Amplitude.y`->LatDiff_left_FvsN$PeakAmp_Noise
LatDiff_left_FvsN$`Peak Amplitude.y`<-NULL
LatDiff_left_FvsN$LatDiff<-LatDiff_left_FvsN$PeakLat_FD-LatDiff_left_FvsN$PeakLat_Noise

LatDiff_mid_FD<-NcPeaks_mid[NcPeaks_mid$`Trial Type`=='FaceDirect',] 
LatDiff_mid_FD$`Trial Type`<-droplevels(LatDiff_mid_FD$`Trial Type`)
LatDiff_mid_Noise<-NcPeaks_mid[NcPeaks_mid$`Trial Type`=='Noise',]
LatDiff_mid_Noise$`Trial Type`<-droplevels(LatDiff_mid_Noise$`Trial Type`)
LatDiff_mid_FvsN<-merge(LatDiff_mid_FD,LatDiff_mid_Noise, by.x=c('Subject'), by.y=c('Subject'))
LatDiff_mid_FvsN$`Trial Type.x`<-NULL
LatDiff_mid_FvsN$`Trial Type.y`<-NULL
LatDiff_mid_FvsN$`Peak Latency.x`->LatDiff_mid_FvsN$PeakLat_FD
LatDiff_mid_FvsN$`Peak Latency.x`<-NULL
LatDiff_mid_FvsN$`Peak Amplitude.x`->LatDiff_mid_FvsN$PeakAmp_FD
LatDiff_mid_FvsN$`Peak Amplitude.x`<-NULL
LatDiff_mid_FvsN$`Peak Latency.y`->LatDiff_mid_FvsN$PeakLat_Noise
LatDiff_mid_FvsN$`Peak Latency.y`<-NULL
LatDiff_mid_FvsN$`Peak Amplitude.y`->LatDiff_mid_FvsN$PeakAmp_Noise
LatDiff_mid_FvsN$`Peak Amplitude.y`<-NULL
LatDiff_mid_FvsN$LatDiff<-LatDiff_mid_FvsN$PeakLat_FD-LatDiff_mid_FvsN$PeakLat_Noise

LatDiff_right_FD<-NcPeaks_right[NcPeaks_right$`Trial Type`=='FaceDirect',]
LatDiff_right_FD$`Trial Type`<-droplevels(LatDiff_right_FD$`Trial Type`)
LatDiff_right_Noise<-NcPeaks_right[NcPeaks_right$`Trial Type`=='Noise',]
LatDiff_right_Noise$`Trial Type`<-droplevels(LatDiff_right_Noise$`Trial Type`)
LatDiff_right_FvsN<-merge(LatDiff_right_FD,LatDiff_right_Noise, by.x=c('Subject'), by.y=c('Subject'))
LatDiff_right_FvsN$`Trial Type.x`<-NULL
LatDiff_right_FvsN$`Trial Type.y`<-NULL
LatDiff_right_FvsN$`Peak Latency.x`->LatDiff_right_FvsN$PeakLat_FD
LatDiff_right_FvsN$`Peak Latency.x`<-NULL
LatDiff_right_FvsN$`Peak Amplitude.x`->LatDiff_right_FvsN$PeakAmp_FD
LatDiff_right_FvsN$`Peak Amplitude.x`<-NULL
LatDiff_right_FvsN$`Peak Latency.y`->LatDiff_right_FvsN$PeakLat_Noise
LatDiff_right_FvsN$`Peak Latency.y`<-NULL
LatDiff_right_FvsN$`Peak Amplitude.y`->LatDiff_right_FvsN$PeakAmp_Noise
LatDiff_right_FvsN$`Peak Amplitude.y`<-NULL
LatDiff_right_FvsN$LatDiff<-LatDiff_right_FvsN$PeakLat_FD-LatDiff_right_FvsN$PeakLat_Noise

LatDiff_left_FvsN$region<-rep("left", nrow(LatDiff_left_FvsN))
LatDiff_mid_FvsN$region<-rep("mid", nrow(LatDiff_mid_FvsN))
LatDiff_right_FvsN$region<-rep("right", nrow(LatDiff_right_FvsN))

Lat_left_mid<-rbind.data.frame(LatDiff_left_FvsN,LatDiff_mid_FvsN)
Lat_all<-rbind.data.frame(Lat_left_mid,LatDiff_right_FvsN)

Lat_all_L <- gather(Lat_all[,-c(3,5)], condition, PeakLat, c(PeakLat_Noise,PeakLat_FD), factor_key=TRUE)
Lat_all_L$condition<-as.character(Lat_all_L$condition) #rename conditions
Lat_all_L$condition[Lat_all_L$condition=='PeakLat_Noise']<-'Noise'
Lat_all_L$condition[Lat_all_L$condition=='PeakLat_FD']<-'FaceDirect'
head(Lat_all_L) 


all_Nc_reg_L<-merge(meanAmp_all_L, Lat_all_L, by.x=c("Subject", "Stimulus", "region"), by.y=c("Subject", "condition", "region"), all=F)
all_Nc_reg_L$Subject<-as.numeric(as.character(all_Nc_reg_L$Subject))
all_Nc_reg_L_phen<-merge(all_Nc_reg_L, phenInfo, by.x="Subject", by.y="ID", all=F)

all_Nc_reg_L_phen$Subject<-as.factor(all_Nc_reg_L_phen$Subject)
all_Nc_reg_L_phen$CA_8mo_days<-as.numeric(all_Nc_reg_L_phen$CA_8mo_days)
all_Nc_reg_L_phen$Outcome<-as.factor(all_Nc_reg_L_phen$Outcome)
all_Nc_reg_L_phen$region<-as.factor(all_Nc_reg_L_phen$region)
all_Nc_reg_L_phen$Stimulus<-as.factor(all_Nc_reg_L_phen$Stimulus)
all_Nc_reg_L_phen$sex<-as.factor(all_Nc_reg_L_phen$sex)





#### dataset for regression (Relation to VABS Socialization Scores)

###select all frontal electrodes
all<-combData[,c("Subject","Stimulus", "Time","V18","V16","V10","V19","V11","V4","V12","V5","V34","V35","V27","V28","V23","V24","V29","V20","V3","V123","V124","V117","V116","V110","V118","V111")]
all$ELmean<-rowSums(combData[,c("V18","V16","V10","V19","V11","V4","V12","V5","V34","V35","V27","V28","V23","V24","V29","V20","V3","V123","V124","V117","V116","V110","V118","V111")])/24 #mean of amplitudes 
all_Nc<-all[all$Time>299,] #from 300 ms from the stimulus onset
all_Nc_meanAmp<-ddply(all_Nc[,c(1,2,28)],.(Subject,Stimulus),summarize,MeanAmp=mean(ELmean)) #compute mean amplitude per subject by stimulus


###compute latencies (and peak amplitude) by stimulus

#all anterior electrodes together
NcPeaks<-p.measures(combData, electrodes=c("V18","V16","V10","V19","V11","V4","V12","V5","V34","V35","V27","V28","V23","V24","V29","V20","V3","V123","V124","V117","V116","V110","V118","V111"),window=c(300,794),num.pts = 0, pol="neg") #peak latency

#compute latency difference between face and noise
FD<-NcPeaks[NcPeaks$`Trial Type`=='FaceDirect',]
FD$`Trial Type`<-droplevels(FD$`Trial Type`)
FA<-NcPeaks[NcPeaks$`Trial Type`=='FaceAverted',]
FA$`Trial Type`<-droplevels(FA$`Trial Type`)
Noise<-NcPeaks[NcPeaks$`Trial Type`=='Noise',]
Noise$`Trial Type`<-droplevels(Noise$`Trial Type`)
FvsN<-merge(FD,Noise, by.x=c('Subject'), by.y=c('Subject'))

#rename columns
FvsN$`Trial Type.x`<-NULL
FvsN$`Trial Type.y`<-NULL
FvsN$`Peak Latency.x`->FvsN$PeakLat_FD
FvsN$`Peak Latency.x`<-NULL
FvsN$`Peak Amplitude.x`->FvsN$PeakAmp_FD
FvsN$`Peak Amplitude.x`<-NULL
FvsN$`Peak Latency.y`->FvsN$PeakLat_Noise
FvsN$`Peak Latency.y`<-NULL
FvsN$`Peak Amplitude.y`->FvsN$PeakAmp_Noise
FvsN$`Peak Amplitude.y`<-NULL
FvsN$LatDiff<-FvsN$PeakLat_FD-FvsN$PeakLat_Noise

#make long-format file for plot Latency Face Direct Gaze vs Noise (region and stimulus as within-subject variables)
FvsN_L <- gather(FvsN[,-c(3,5)], condition, PeakLat, c(PeakLat_Noise,PeakLat_FD), factor_key=TRUE)
FvsN_L$condition<-as.character(FvsN_L$condition) #rename conditions
FvsN_L$condition[FvsN_L$condition=='PeakLat_Noise']<-'Noise'
FvsN_L$condition[FvsN_L$condition=='PeakLat_FD']<-'FaceDirect'
head(FvsN_L) 


#merge datasets mean amplitude and peak latency
all_Nc_meanAmp$Subject<-as.numeric(as.character(all_Nc_meanAmp$Subject))
FvsN_L$Subject<-as.numeric(as.character(FvsN_L$Subject))
FvsN_L<-FvsN_L[complete.cases(FvsN_L$Subject),] #remove grand averages
Nc_AmpLat_L<-merge(all_Nc_meanAmp,FvsN_L, by.x=c("Subject", "Stimulus"), by.y=c("Subject","condition"), all=T)
#merge phenotypic info
Nc_AmpLat_L_phen<-merge(Nc_AmpLat_L, phenInfo, by.x='Subject', by.y='ID', all.x =T) #with condition as wsv






###################################### Relation to categorical outcome of ASD ###############################################


# check assumptions for ANOVA Mean Amplitude
Amp_W<-spread(Nc_AmpLat_L[,c(1,2,3)], key=Stimulus, value= MeanAmp)
Amp_W_phen<-merge(Amp_W, phenInfo, by.x='Subject', by.y='ID', all=F)

by(Amp_W_phen$FaceDirect,Amp_W_phen$Outcome,  shapiro.test)
by(Amp_W_phen$Noise,Amp_W_phen$Outcome,  shapiro.test)

leveneTest(Amp_W_phen$FaceDirect,Amp_W_phen$Outcome)
leveneTest(Amp_W_phen$Noise,Amp_W_phen$Outcome)

# check assumptions for ANOVA Latency
Lat_W<-spread(Nc_AmpLat_L[,c(1,2,5)], key=Stimulus, value= PeakLat)
Lat_W_phen<-merge(Lat_W, phenInfo, by.x='Subject', by.y='ID', all=F)

by(Lat_W_phen$FaceDirect,Lat_W_phen$Outcome,  shapiro.test)
kurtosis(Lat_W_phen[Lat_W_phen$Outcome=='1',]$FaceDirect)
skewness(Lat_W_phen[Lat_W_phen$Outcome=='1',]$FaceDirect)

by(Lat_W_phen$Noise,Lat_W_phen$Outcome,  shapiro.test)
kurtosis(Lat_W_phen[Lat_W_phen$Outcome=='0',]$FaceDirect)
skewness(Lat_W_phen[Lat_W_phen$Outcome=='0',]$FaceDirect)
kurtosis(Lat_W_phen[Lat_W_phen$Outcome=='1',]$FaceDirect)
skewness(Lat_W_phen[Lat_W_phen$Outcome=='1',]$FaceDirect)

leveneTest(Lat_W_phen$FaceDirect,Lat_W_phen$Outcome)
leveneTest(Lat_W_phen$Noise,Lat_W_phen$Outcome)



##### ANOVA MEAN AMPLITUDE #####

all_Nc_reg_L_phen$CA_8mo_days<-as.numeric(all_Nc_reg_L_phen$CA_8mo_days)
all_Nc_reg_L_phen$M_ELC_SS_1<-as.numeric(all_Nc_reg_L_phen$M_ELC_SS_1)

aov_MeanAmp<-aov(MeanAmp~Stimulus*region*Outcome + sex + CA_8mo_days + M_ELC_SS_1 + Error(Subject), data=all_Nc_reg_L_phen[complete.cases(all_Nc_reg_L_phen$M_ELC_SS_1),]) #control for effect of sex
summary(aov_MeanAmp)
EtaSq(aov_MeanAmp, type=1, anova=TRUE)

#without covariates
aov_MeanAmp_wc<-aov(MeanAmp~Outcome*region*Stimulus + Error(Subject), data=all_Nc_reg_L_phen) 
summary(aov_MeanAmp_wc)
EtaSq(aov_MeanAmp_wc, type=1, anova=TRUE)

ddply(all_Nc_reg_L_phen, "Stimulus", summarise, N=length(unique(Subject)), mean=mean(MeanAmp), sd=sd(MeanAmp))

#plot significant effects of covariates
plot(all_Nc_reg_L_phen$MeanAmp~all_Nc_reg_L_phen$CA_8mo_days)
abline(lm(all_Nc_reg_L_phen$MeanAmp~all_Nc_reg_L_phen$CA_8mo_days))
plot(all_Nc_reg_L_phen$MeanAmp~all_Nc_reg_L_phen$ M_ELC_SS_1)
abline(lm(all_Nc_reg_L_phen$MeanAmp~all_Nc_reg_L_phen$ M_ELC_SS_1))

##explore significant effect of stimulus by group
aov_MeanAmp_TL<-aov(MeanAmp~Stimulus+CA_8mo_days + M_ELC_SS_1 + sex +Error(Subject), data=all_Nc_reg_L_phen[all_Nc_reg_L_phen$Outcome=='0'&complete.cases(all_Nc_reg_L_phen$M_ELC_SS_1),])
summary(aov_MeanAmp_TL)
EtaSq(aov_MeanAmp_TL, type=1, anova=TRUE)

aov_MeanAmp_ELnoASD<-aov(MeanAmp~Stimulus+CA_8mo_days+ M_ELC_SS_1 + sex+ Error(Subject), data=all_Nc_reg_L_phen[all_Nc_reg_L_phen$Outcome=='1'&complete.cases(all_Nc_reg_L_phen$M_ELC_SS_1),])
summary(aov_MeanAmp_ELnoASD)
EtaSq(aov_MeanAmp_ELnoASD, type=1, anova=TRUE)

aov_MeanAmp_ELASD<-aov(MeanAmp~Stimulus+CA_8mo_days+ M_ELC_SS_1 + sex+ Error(Subject), data=all_Nc_reg_L_phen[all_Nc_reg_L_phen$Outcome=='3'&complete.cases(all_Nc_reg_L_phen$M_ELC_SS_1),])
summary(aov_MeanAmp_ELASD)
EtaSq(aov_MeanAmp_ELASD, type=1, anova=TRUE)

ddply(all_Nc_reg_L_phen, c("Stimulus","Outcome"), summarise, N=length(unique(Subject)), mean=mean(MeanAmp), sd=sd(MeanAmp))

#explore significant effect of outcome by stimulus
aov_MeanAmp_FD<-aov(MeanAmp~Outcome +CA_8mo_days + M_ELC_SS_1 + sex+ Error(Subject), data=all_Nc_reg_L_phen[complete.cases(all_Nc_reg_L_phen$M_ELC_SS_1)& all_Nc_reg_L_phen$Stimulus=='FaceDirect',])
summary(aov_MeanAmp_FD)
EtaSq(aov_MeanAmp_FD, type=1, anova=TRUE)

aov_MeanAmp_Noise<-aov(MeanAmp~Outcome +CA_8mo_days + M_ELC_SS_1 + sex + Error(Subject), data=all_Nc_reg_L_phen[complete.cases(all_Nc_reg_L_phen$M_ELC_SS_1)& all_Nc_reg_L_phen$Stimulus=='Noise',])
summary(aov_MeanAmp_Noise)
EtaSq(aov_MeanAmp_Noise, type=1, anova=TRUE)



##### ANOVA PEAK LATENCY #####

aov_PeakLat<-aov(PeakLat~Stimulus*region*Outcome+sex+CA_8mo_days+ M_ELC_SS_1 + Error(Subject), data=all_Nc_reg_L_phen[complete.cases(all_Nc_reg_L_phen$M_ELC_SS_1),]) #control for effect of sex
summary(aov_PeakLat)
EtaSq(aov_PeakLat, type=1, anova=TRUE)

#without covariates
aov_PeakLat_wc<-aov(PeakLat~Stimulus*region*Outcome + Error(Subject), data=all_Nc_reg_L_phen) 
summary(aov_PeakLat_wc)
EtaSq(aov_PeakLat_wc, type=1, anova=TRUE)

TukeyHSD(aov(PeakLat~Outcome , data=all_Nc_reg_L_phen[complete.cases(all_Nc_reg_L_phen$M_ELC_SS_1),]))
ddply(all_Nc_reg_L_phen, "Outcome", summarise, N=length(unique(Subject)), mean=mean(PeakLat), sd=sd(PeakLat))


##explore significant effect of stimulus by group
aov_PeakLat_TL<-aov(PeakLat~Stimulus+CA_8mo_days + M_ELC_SS_1 + sex +Error(Subject), data=all_Nc_reg_L_phen[all_Nc_reg_L_phen$Outcome=='0'&complete.cases(all_Nc_reg_L_phen$M_ELC_SS_1),])
summary(aov_PeakLat_TL)
EtaSq(aov_PeakLat_TL, type=1, anova=TRUE)

aov_PeakLat_ELnoASD<-aov(PeakLat~Stimulus+CA_8mo_days+ M_ELC_SS_1 + sex+ Error(Subject), data=all_Nc_reg_L_phen[all_Nc_reg_L_phen$Outcome=='1'&complete.cases(all_Nc_reg_L_phen$M_ELC_SS_1),])
summary(aov_PeakLat_ELnoASD)
EtaSq(aov_PeakLat_ELnoASD, type=1, anova=TRUE)

aov_PeakLat_ELASD<-aov(PeakLat~Stimulus+CA_8mo_days+ M_ELC_SS_1 + sex+ Error(Subject), data=all_Nc_reg_L_phen[all_Nc_reg_L_phen$Outcome=='3'&complete.cases(all_Nc_reg_L_phen$M_ELC_SS_1),])
summary(aov_PeakLat_ELASD)
EtaSq(aov_PeakLat_ELASD, type=1, anova=TRUE)

#explore significant effect of outcome by stimulus
aov_PeakLat_FD<-aov(PeakLat~Outcome +CA_8mo_days + M_ELC_SS_1 + sex+ Error(Subject), data=all_Nc_reg_L_phen[complete.cases(all_Nc_reg_L_phen$M_ELC_SS_1)& all_Nc_reg_L_phen$Stimulus=='FaceDirect',])
summary(aov_PeakLat_FD)
EtaSq(aov_PeakLat_FD, type=1, anova=TRUE)
TukeyHSD(aov(PeakLat~Outcome, data=all_Nc_reg_L_phen[complete.cases(all_Nc_reg_L_phen$M_ELC_SS_1)& all_Nc_reg_L_phen$Stimulus=='FaceDirect',]))

aov_PeakLat_Noise<-aov(PeakLat~Outcome +CA_8mo_days + M_ELC_SS_1 + sex + Error(Subject), data=all_Nc_reg_L_phen[complete.cases(all_Nc_reg_L_phen$M_ELC_SS_1)& all_Nc_reg_L_phen$Stimulus=='Noise',])
summary(aov_PeakLat_Noise)
EtaSq(aov_PeakLat_Noise, type=1, anova=TRUE)
TukeyHSD(aov(PeakLat~Outcome, data=all_Nc_reg_L_phen[complete.cases(all_Nc_reg_L_phen$M_ELC_SS_1)& all_Nc_reg_L_phen$Stimulus=='Noise',]))

ddply(all_Nc_reg_L_phen, c("Stimulus","Outcome"), summarise, N=length(unique(Subject)), mean=mean(PeakLat), sd=sd(PeakLat))


################ PLOTS ############################
#plot mean amplitude in face and noise by group
sumSE_NcAmpDiff <- summarySEwithin(Nc_AmpLat_L_phen, measurevar="MeanAmp", betweenvars="Outcome", withinvars='Stimulus', idvar='Subject')

#jpeg(file="FvsN_amp_Outcome.jpg",width = 4, height = 4, units = 'in', res = 300)
FvsN_amp_Outcome<-ggplot(sumSE_NcAmpDiff, aes(x=Outcome, y=MeanAmp, shape=Stimulus, fill=Outcome)) + 
  xlab("Outcome groups") +
  ggtitle("Average frontal regions") +
  scale_x_discrete(breaks=c("0", "1", "3"),labels=c("TL", "EL-noASD", "EL-ASD")) +
  ylab(expression(Nc~mean~amplitude~(mu~V))) +
  geom_point(data=sumSE_NcAmpDiff, aes(x=Outcome, y=MeanAmp, shape=Stimulus,fill=Outcome), size=8, stat="identity",position=position_dodge(.8), color="black")+
  geom_errorbar(aes(ymin=MeanAmp-se, ymax=MeanAmp+se, color=Outcome), width=.2, position=position_dodge(.8) ) +
  scale_colour_manual(breaks=c("0", "1", "3"),values=c('springgreen4', 'royalblue', 'firebrick'), labels=c("TL", "EL-noASD", "EL-ASD")) + 
  scale_fill_manual(breaks=c("0", "1", "3"),values=c('springgreen4', 'royalblue',  'firebrick'), labels=c("TL", "EL-noASD", "EL-ASD")) + 
  scale_shape_manual(values=c(23,22),labels=c("Face with Direct Gaze", "Noise")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="bottom", plot.title = element_text(hjust = 0.5)) + guides(color=F,fill=F)
print(FvsN_amp_Outcome)
dev.off()



#plot latency in face and noise by group

sumSE_NcLatDiff <- summarySEwithin(Nc_AmpLat_L_phen, measurevar="PeakLat", betweenvars="Outcome", withinvars='Stimulus', idvar='Subject')

#jpeg(file="FvsN_lat_Outcome.jpg",width = 4, height = 4, units = 'in', res = 300)
FvsN_lat_Outcome<-ggplot(sumSE_NcLatDiff, aes(x=Outcome, y=PeakLat, shape=Stimulus, fill=Outcome)) + 
  ggtitle("Average frontal regions") +
  xlab("Outcome groups") +
  scale_x_discrete(breaks=c("0", "1","3"),labels=c("TL", "EL-noASD", "EL-ASD")) +
  ylab("Nc peak latency (ms)") +
  geom_point(data=sumSE_NcLatDiff, aes(x=Outcome, y=PeakLat, shape=Stimulus,fill=Outcome), size=8, stat="identity",position=position_dodge(.8))+
  geom_errorbar(aes(ymin=PeakLat-se, ymax=PeakLat+se, color=Outcome), width=0.2, position=position_dodge(.8) ) +
  scale_colour_manual(breaks=c("0", "1", "3"),values=c('springgreen4', 'royalblue', 'firebrick'), labels=c("TL", "EL-noASD", "EL-ASD")) + 
  scale_fill_manual(breaks=c("0", "1", "3"),values=c('springgreen4', 'royalblue', 'firebrick'), labels=c("TL", "EL-noASD", "EL-ASD")) + 
  scale_shape_manual(values=c(23,22),labels=c("Face with Direct Gaze", "Noise")) +  theme_bw() +
  theme(axis.line = element_line(colour = "black"), text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="bottom", plot.title = element_text(hjust = 0.5)) + guides(color=F,fill=F)
print(FvsN_lat_Outcome)
dev.off()


# combine plots

jpeg('/Users/annagui/Library/Mobile Documents/com~apple~CloudDocs/Birkbeck/ERPs/Files_for_paper/TP_review/ANOVAeffects_Nc_Outcome.jpeg', units='in', width=12, height = 6, res = 300)
ANOVAeffects_Nc_Outcome<-ggarrange(FvsN_amp_Outcome, FvsN_lat_Outcome, labels = c("A", "B"),
                                    ncol = 2, nrow = 1,
                                    common.legend = TRUE, legend = "bottom")
print(ANOVAeffects_Nc_Outcome)
dev.off()






##################################### Relation to VABS Socialization scores #############################################

Amp_W_phen$AmpDiff<-Amp_W_phen$FaceDirect-Amp_W_phen$Noise #amplitude difference score between FD and Noise
Lat_W_phen$LatDiff<-Lat_W_phen$FaceDirect-Lat_W_phen$Noise #latency difference score between FD and Noise


#### REGRESSION MEAN AMPLITUDE DIFFERENCE SCORE ####
Amp_W_phen$V_Soc_SS_4<-as.numeric(Amp_W_phen$V_Soc_SS_4)
Amp_W_phen$CA_8mo_days<-as.numeric(Amp_W_phen$CA_8mo_days)

summary(lm(V_Soc_SS_4 ~ AmpDiff + sex + CA_8mo_days + as.numeric(M_ELC_SS_1), data=Amp_W_phen)) #OLS control for effect of sex
#summary(lm(V_Soc_SS_4 ~ AmpDiff , data=Amp_W_phen)) #no covariates


##check assumptions
lmMA<-lm(V_Soc_SS_4 ~ AmpDiff + sex + CA_8mo_days + as.numeric(M_ELC_SS_1), data=Amp_W_phen)
mean(lmMA$residuals) #mean residuals = 0

MA_Vabs<-Amp_W_phen[complete.cases(Amp_W_phen$V_Soc_SS_4)&complete.cases(as.numeric(Amp_W_phen$M_ELC_SS_1)),]$AmpDiff
cor.test(MA_Vabs, as.vector(lmMA$residuals)) #no correlation between X and residuals

lmtest::dwtest(lmMA) #no autocorrelation of the residuals
vif(lmMA) #no multicollinearity
gvlma::gvlma(lmMA) #normality and homoschedasticity

#plot diagnostics to check homoschedasticity, normality of residuals, heteroschedasticity and outliers
jpeg('lmMA_assumptions.jpeg', units='in', width=6, height = 6, res = 300)
par(mfrow=c(2,2)) 
plot(lmMA) 
print(plot(lmMA))
dev.off()


#check if effect is significant for non-social scale
Amp_W_phen$V_Mot_SS_4<-as.numeric(Amp_W_phen$V_Mot_SS_4)
summary(lm(V_Mot_SS_4 ~ AmpDiff + sex + CA_8mo_days + as.numeric(M_ELC_SS_1), data=Amp_W_phen)) 
#summary(lm(V_Mot_SS_4 ~ AmpDiff, data=Amp_W_phen)) #no covariares


#### REGRESSION PEAK LATENCY DIFFERENCE SCORE ####

Lat_W_phen$V_Soc_SS_4<-as.numeric(Lat_W_phen$V_Soc_SS_4)
Lat_W_phen$CA_8mo_days<-as.numeric(Lat_W_phen$CA_8mo_days)
summary(lm(V_Soc_SS_4 ~ LatDiff + sex + CA_8mo_days + as.numeric(M_ELC_SS_1), data=Lat_W_phen)) 
#summary(lm(V_Soc_SS_4 ~ LatDiff, data=Lat_W_phen)) #no covariates

#check assumptions
lmPL<-lm(V_Soc_SS_4 ~ LatDiff + sex + CA_8mo_days + as.numeric(M_ELC_SS_1), data=Lat_W_phen)
mean(lmPL$residuals) #mean residuals = 0

PL_Vabs<-Lat_W_phen[complete.cases(Lat_W_phen$V_Soc_SS_4)&complete.cases(as.numeric(Lat_W_phen$M_ELC_SS_1)),]$LatDiff
cor.test(PL_Vabs, as.vector(lmPL$residuals)) #no correlation between X and residuals

lmtest::dwtest(lmPL) #no autocorrelation of the residuals
vif(lmPL) #no multicollinearity
gvlma::gvlma(lmPL) #normality and homoschedasticity

#plot diagnostics to check homoschedasticity, normality of residuals, heteroschedasticity and outliers
jpeg('lmPL_assumptions.jpeg', units='in', width=6, height = 6, res = 300)
par(mfrow=c(2,2)) 
plot(lmPL) 
print(plot(lmPL))
dev.off()

#check if effect is significant for non-social scale
Lat_W_phen$V_Mot_SS_4<-as.numeric(Lat_W_phen$V_Mot_SS_4)
summary(lm(V_Mot_SS_4 ~ LatDiff + sex + CA_8mo_days + as.numeric(M_ELC_SS_1), data=Lat_W_phen)) 
#summary(lm(V_Mot_SS_4 ~ LatDiff, data=Lat_W_phen)) # no covariates




#plot relationship VABS - mean amplitude difference score

#jpeg(file="regrVABS_FvsN_amp_Outcome.jpeg",width = 7, height = 4, units = 'in', res = 300)
regrVABS_FvsN_amp_Outcome<-ggplot() +
  geom_smooth(method='lm', fullrange=T, data=Amp_W_phen, aes(y=V_Soc_SS_4, x= AmpDiff),colour="black") +
  geom_point(data= Amp_W_phen, aes(y=V_Soc_SS_4, x=AmpDiff, color=Outcome, shape=Outcome), size=2)  +
  xlab("Nc mean amplitude difference score Face with Direct Gaze - Noise") +
  ylab("VABS - Socialization Score at 3 years") +
  scale_color_manual(values=c("0"='springgreen4', "1"='royalblue', "3"= 'firebrick'), name="Outcome groups", labels=c("TL", "EL-noASD", "EL-ASD"), breaks=levels(factor(Amp_W_phen$Outcome))) +
  scale_shape_manual(name="Outcome groups", values=c(15,17,18), labels=c("TL", "EL-noASD", "EL-ASD")) +
  theme_bw() + 
  theme(axis.title.x = element_text(vjust=-1, size=10, hjust=0.5), axis.line = element_line(colour = "black"), text = element_text(size=12),  panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  guides(fill=guide_legend(override.aes=list(shape=21))) 
print(regrVABS_FvsN_amp_Outcome)
dev.off()


#plot relationship VABS - peak latency difference score

#jpeg(file="regrVABS_FvsN_lat_Outcome.jpeg",width = 7, height = 4, units = 'in', res = 300)
regrVABS_FvsN_lat_Outcome<-ggplot() +
  geom_smooth(method='lm', fullrange=T, data=Lat_W_phen, aes(y=V_Soc_SS_4, x= LatDiff),colour="black") +
  geom_point(data= Lat_W_phen, aes(y=V_Soc_SS_4, x=LatDiff, color=Outcome, shape=Outcome), size=2)  +
  xlab("Nc latency difference score Face with Direct Gaze - Noise") +
  ylab("VABS - Socialization Score at 3 years") +
  scale_color_manual(values=c("0"='springgreen4', "1"='royalblue',"3"= 'firebrick'), name="Outcome groups", labels=c("TL", "EL-noASD", "EL-ASD"), breaks=levels(factor(Amp_W_phen$Outcome))) +
  scale_shape_manual(name="Outcome groups", values=c(15,17,18), labels=c("TL", "EL-noASD", "EL-ASD")) +
  theme_bw() + 
  theme(axis.title.x = element_text(vjust=-1, size=10, hjust=0.5), axis.line = element_line(colour = "black"), text = element_text(size=12),  panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  guides(fill=guide_legend(override.aes=list(shape=21))) 
print(regrVABS_FvsN_lat_Outcome)
dev.off()




# combine plots 


jpeg('two_plots_NcXvabs_Outcome.jpeg', units='in', width=12, height = 5, res = 300)
two_plots_NcXvabs_Outcome<-ggarrange(regrVABS_FvsN_amp_Outcome, regrVABS_FvsN_lat_Outcome, labels = c("A", "B"),
                                      ncol = 2, 
                                      common.legend = TRUE, legend = "bottom")
print(two_plots_NcXvabs_Outcome)
dev.off()



