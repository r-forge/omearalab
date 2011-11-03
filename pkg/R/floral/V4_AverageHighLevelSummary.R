#need to remove the second deltaAIC and AICweight
setwd("../Summaries")
load("Highlevel.dataframe.withrates.Rsave")
hl<-highlevel.dataframe[,-3] #removes second deltaAIC
hl<-highlevel.dataframe[,-3] #removes second AICweight
highlevel.dataframe<-hl
save(file="Highlevel.dataframe.withrates.fixedAIC.Rsave",compress=TRUE)
parameter.dataframe<-highlevel.dataframe[,11:dim(highlevel.dataframe)[2] ]

multiplyVector<-function(x,weights) {
	return(x*weights)
}

parameter.dataframe<-apply(parameter.dataframe,2,multiplyVector,weights=highlevel.dataframe$AICweight)


save(file="Highlevel.dataframe.withrates.fixedAIC.Rsave",compress=TRUE)
averageVector<-colSums(parameter.dataframe)
save(file="Highlevel.dataframe.withrates.fixedAIC.Rsave",compress=TRUE)
names(averageVector)<-names(highlevel.dataframe)[11:dim(highlevel.dataframe)[2] ]
save(file="Highlevel.dataframe.withrates.fixedAIC.Rsave",compress=TRUE)
print(averageVector)

