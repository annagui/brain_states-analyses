
#### R script to select the optimal number of microstates
### Gui et al., Diminished engagememnt of attentive brain states predicts later Autism Spectrum Disorder
## 06 March 2020

EV<-read.csv("ExplainedVariance_microstates_FD.csv", header=F, stringsAsFactors=F) 


EV<-as.data.frame(EV)
colnames(EV)<-EV[1,]
EV<-EV[2:nrow(EV),]
rownames(EV)<-EV[,1]
EV<-EV[,2:ncol(EV)]


EV<-t(EV)
EV<-as.data.frame(EV)
EV_L<-gather(EV, NrMs, ExpVar, n1_ms:n10_ms, factor_key = TRUE)

glmNrMs<-glm(ExpVar~NrMs, data=EV_L, family='gaussian')
summary(glmNrMs)
aovExpVar<-aov(ExpVar~NrMs, data=EV_L)
summary(aov(ExpVar~NrMs, data=EV_L))


pairwise.t.test(EV_L$ExpVar,EV_L$NrMs ,'bonf' , alternative='greater')


