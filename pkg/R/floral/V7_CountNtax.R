library(RColorBrewer)
library(stringr)
setwd("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/SimsMay2013FINAL")
source("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral/V7_NewSimulator.R")
source("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral/V7_StochasticSSASims_Functions.R")
source("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral/V7_UtilityFns.R")
source("/Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral/V7_NewSimulatorJustGetRates.R")

GetFirstAppearance<-function(counts) {
  age<-137
  nonzeros<-which(counts!=0)
  if(length(nonzeros)>0) {
    age<-nonzeros[1]-1 
  }
  return(age)
}

in.progress<-read.csv("~/Dropbox/ResultsPeekPruned.csv")
in.progress$ntax.rescale<-as.numeric(str_extract(str_extract(in.progress$file, "scale_\\d+e*\\+*\\d*"), "\\d+e*\\+*\\d*"))
in.progress$rate.est<-log(in.progress$ntax/2)/in.progress$time
in.progress$ntax.predicted<-2*exp(136*in.progress$rate.est)

filter.v<-c("full")
results.df<-c()
for (filter.index in sequence(length(filter.v))) {
  filter=filter.v[filter.index]
  dirs<-system(paste("ls -1 | grep old | grep ", filter, " | grep -v RData | grep -v '\\.R' | grep -iv pdf", sep=""),intern=TRUE)
  #mypalette<-brewer.pal(8,"Dark2")
  mypalette<-heat.colors(12,alpha=0.3)
  original.data<-read.csv("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/RunsJan2012/SourceData/Stebbins_prunenoper25i2012BCO.csv",stringsAsFactors=FALSE)
  observed.ntax<-dim(original.data)[1]
  all.histories.list<-list()
  combo.names<-c("0x00xx","0x01xx","0x10xx","0x11xx","1x00xx","1x01xx","1x10xx","1x11xx")
  
  
  observed.proportions<-rep(NA, length(combo.names))
  for (i in sequence(length(combo.names))) {
    chars<-strsplit(combo.names[i],"")[[1]]
    small.dataset<-original.data[which(original.data[,1+1]==chars[1]),] #look at the relevant chars
    small.dataset<-small.dataset[which(small.dataset[,3+1]==chars[3]),]
    small.dataset<-small.dataset[which(small.dataset[,4+1]==chars[4]),]
    observed.proportions[i]<-dim(small.dataset)[1]/dim(original.data)[1]
  }
  names(observed.proportions)<-combo.names
  
  system("rm *.pdf")
  
  for (dir.index in sequence(length(dirs))) {
    ntax.rescale<-as.numeric(str_extract(str_extract(dirs[dir.index], "scale_\\d+e*\\+*\\d*"), "\\d+e*\\+*\\d*"))
    proportion.df<-c()
    #print(dirs[dir.index])
    setwd(dirs[dir.index])
    ntax.vector<-c()
    ages.df<-c()
    load("Rates.Rsave")
    tf<-136
    all.matrix<-OMearaSSASave(x0, q.vector=q.means, lambda.vector=lambda.means, mu.vector=mu.means, tf=136, verbose=FALSE, print.freq=100, full.history=FALSE, rescale.species=250000, yule.scale=0, history.steps.to.save=seq(from=1,to=floor(tf),length.out=floor(tf)), t.rescale=tf, x0.rescale=x0.rescale, ntax.old.scale=ntax.rescale)
    lambda.actual<-all.matrix[,9]
    names(lambda.actual)<-paste("lambda", rownames(all.matrix), sep="_")
    mu.actual<-all.matrix[,10]
    names(mu.actual)<-paste("mu", rownames(all.matrix), sep="_")
    div.actual<-all.matrix[,9]-all.matrix[,10]
    names(div.actual)<-paste("div", rownames(all.matrix), sep="_")
    turnover.actual<-all.matrix[,9]+all.matrix[,10]
    names(turnover.actual)<-paste("turnover", rownames(all.matrix), sep="_")
    q.actual<-c()
    q.actual.names<-c()
    for (row.index in sequence(8)) {
      for (col.index in sequence(8)) {
        q.actual<-append(q.actual, all.matrix[row.index, col.index])
        q.actual.names<-append(q.actual.names, paste("q",rownames(all.matrix)[row.index], "_", colnames(all.matrix)[col.index], sep=""))
      }
    }
    names(q.actual)<-q.actual.names
    #    file.list<-system("ls -1 | grep total151.RSave",intern=TRUE)
    #print(round(all.matrix, 2))
    file.list<-system("ls -1 | grep RSave",intern=TRUE)
    if(grepl("old", dirs[dir.index])) {
      file.list<-system("ls -1 | grep RSave | grep -v app1", intern=TRUE)
    }
    if (length(file.list)>0) {
      file.list<-paste("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/SimsMay2013FINAL/", dirs[dir.index], "/", file.list, sep="")
    }
    #  if(!grepl("stacey", dirs[dir.index])) {
    #    setwd("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/SimsMay2013Additional")
    #    setwd(dirs[dir.index])
    #    file.list2<-system("ls -1 | grep total137.RSave",intern=TRUE)
    #    if (length(file.list2)>0) {
    #      file.list2<-paste("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/SimsMay2013Additional/", dirs[dir.index], "/", file.list2, sep="")
    #    }
    #    file.list<-c(file.list, file.list2)
    #  }
    #print(file.list)
    setwd("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/SimsMay2013FINAL")
    setwd(dirs[dir.index])
    distances<-c()
    for (i in sequence(length(file.list))) {
      load(file.list[i])
      #print(tail(history))
      #  print(file.list[i])
      proportional.history<-history
      for (row.index in sequence(dim(history)[1])) {
        proportional.history[row.index, 2:9]<-history[row.index, 2:9]/sum( history[row.index, 2:9])
      }
      distances<-append(distances, dist(as.matrix(rbind(observed.proportions, proportional.history[137, 2:9])))[1]/sqrt(2)) #sqrt(2) is the max distance
      #mat.for.dist<<-as.matrix(rbind(observed.proportions, proportional.history[137, 2:9]))
      #dput(as.matrix(rbind(observed.proportions, proportional.history[137, 2:9])))
      #print(dist(as.matrix(rbind(observed.proportions, proportional.history[137, 2:9])))[1])
      #print(dist(as.matrix(rbind(proportional.history[137, 2:9], proportional.history[137, 2:9])))[1])
      #print(rbind(observed.proportions, proportional.history[137, 2:9]))
      if (i==1) {
        combo.names<-colnames(history)[2:9]
        for (j in 2:9) {
          all.histories.list[[j-1]]<-proportional.history[,c(1,j)] 
        }
      }
      else {
        for (j in 2:9) {
          all.histories.list[[j-1]]<-cbind( all.histories.list[[j-1]], proportional.history[,j])
        }
      }
      ages.appearance<-rep(NA,8)
      try(ages.appearance<-apply(history[,2:9], 2, GetFirstAppearance))
      names(ages.appearance)<-paste("age.first_",combo.names,sep="")
      #print(ages.appearance)
      if(!is.na(ages.appearance[1])) {
      ages.df<-rbind(ages.df, ages.appearance)
      } else {
	print("problem with ages.appearance")
	print(ages.appearance)
      }
      #print(proportional.history)
      # print(paste("ntax present = ",sum(history[137,2:9])))
      ntax.vector<-append(ntax.vector, sum(history[137,2:9]))
      print(paste(ntax.vector[length(ntax.vector)], file.list[i]))
      proportion.df<-rbind(proportion.df, proportional.history[137, 2:9])
      rm(history)
      rm(proportional.history)
    }
    if (length(file.list)>0) {
      proportional.average<-colMeans(proportion.df)
      names(proportional.average)<-paste("proport", combo.names, sep="_")
     # ntax.vector<-unique(ntax.vector)
      ntax.include.predictions.vector<-c(ntax.vector, in.progress$ntax.predicted[which(in.progress$ntax.rescale==ntax.rescale)])
      ntax.include.partials.vector<-c(ntax.vector, in.progress$ntax[which(in.progress$ntax.rescale==ntax.rescale)])
      ages.leave.root<-NA
	if(dim(ages.df)[1]>0) {
      		try(ages.leave.root<-apply(ages.df[,2:8], 1, min))
	}
      results.df<-rbind(results.df, data.frame(model=filter, dir=dirs[dir.index], ntax.rescale=ntax.rescale, ntax.actual=median(ntax.vector), ntax.actual.025=quantile(ntax.vector, 0.025), ntax.actual.975=quantile(ntax.vector, 0.975), ntax.vector.string=paste(ntax.vector, collapse="_"), ntax.predicted=median(ntax.include.predictions.vector), ntax.predicted.025=quantile(ntax.include.predictions.vector, 0.025), ntax.predicted.975=quantile(ntax.include.predictions.vector, 0.975), ntax.partials=median(ntax.include.partials.vector), ntax.partials.025=quantile(ntax.include.partials.vector, 0.025), ntax.partials.975=quantile(ntax.include.partials.vector, 0.975), nruns=length(ntax.vector), nruns.with.predictions=length(ntax.include.predictions.vector), state1.ages=paste(ages.df[,1], collapse="_"), state2.ages=paste(ages.df[,2], collapse="_"), state3.ages=paste(ages.df[,3], collapse="_"), state4.ages=paste(ages.df[,4], collapse="_"), state5.ages=paste(ages.df[,5], collapse="_"), state6.ages=paste(ages.df[,6], collapse="_"), state7.ages=paste(ages.df[,7], collapse="_"), state8.ages=paste(ages.df[,8], collapse="_"),  leave.root.ages=paste(ages.leave.root, collapse="_"), t(c(div.actual, proportional.average, lambda.actual, mu.actual, turnover.actual, q.actual))))
      #print(c(dirs[dir.index], ntax.rescale, median(ntax.vector)))
    }
    setwd("..")
    
  }
  #system("open All.PDF")
}
#print(results.df)
write.csv(results.df, file=paste("~/Dropbox/ntax.Allfilter.results.csv", sep=""))
save(results.df, file="~/Dropbox/ntax.Allfilter.results.Rdata")
pdf(file="~/Dropbox/Regression.pdf")
par(mfcol=c(1,2))
ntax.ideal<-250000
colors<-rainbow(length(results.df$nruns))[order(runif(length(results.df$nruns)))]
plot(x=range(results.df$ntax.rescale), y=range(c(results.df$ntax.actual.025, results.df$ntax.predicted.975)), type="n", bty="n", log="xy",  main="finished")
text(x=results.df$ntax.rescale, y=results.df$ntax.actual, labels=results.df$nruns, col=colors)
for (i in sequence(dim(results.df)[1])) {
  lines(x=rep(results.df$ntax.rescale[i], 2), y=c(results.df$ntax.actual.025[i], results.df$ntax.actual.975[i]), col=colors[i]) 
}
abline(h=ntax.ideal, lty="dotted")
abline(v=158500)
plot(x=range(results.df$ntax.rescale), y=range(c(results.df$ntax.actual.025, results.df$ntax.predicted.975)), type="n", bty="n", log="xy",  main="finished + in progress")
text(x=results.df$ntax.rescale, y=results.df$ntax.predicted, labels=results.df$nruns.with.predictions, col=colors)
for (i in sequence(dim(results.df)[1])) {
  lines(x=rep(results.df$ntax.rescale[i], 2), y=c(results.df$ntax.predicted.025[i], results.df$ntax.predicted.975[i]), col=colors[i]) 
}
abline(h=ntax.ideal, lty="dotted")
#abline(v=158500)
y<-log(results.df$ntax.predicted)
x<-log(results.df$ntax.rescale)
regression1<-lm(y~x)
lines(results.df$ntax.rescale, exp(predict(regression1, data.frame(x=log(results.df$ntax.rescale)))), pch=16)
#plot(x,y, ylim=c(0,25))
regression2<-lm(x~y)
prediction2<-predict(regression2, data.frame(y=log(250000)))
lines(y=results.df$ntax.predicted, x=predict(regression2, data.frame(y=results.df$ntax.predicted)), lwd=3)
#exp(prediction2)
abline(v=exp(prediction2), lwd=4)

dev.off()

best.results.df<-results.df[which(results.df$ntax.rescale==11043),]
wait.times<-data.frame(cbind( state2=as.numeric(strsplit(as.character(best.results.df$state2.ages), "_")[[1]]), state3=as.numeric(strsplit(as.character(best.results.df$state3.ages), "_")[[1]]), state4=as.numeric(strsplit(as.character(best.results.df$state4.ages), "_")[[1]]), state5=as.numeric(strsplit(as.character(best.results.df$state5.ages), "_")[[1]]), state6=as.numeric(strsplit(as.character(best.results.df$state6.ages), "_")[[1]]), state7=as.numeric(strsplit(as.character(best.results.df$state7.ages), "_")[[1]]), state8=as.numeric(strsplit(as.character(best.results.df$state8.ages), "_")[[1]])))
wait.times<-cbind(wait.times, first.leavings=apply(wait.times, 1, min))
wait.times<-cbind(wait.times, ntax=as.numeric(strsplit(as.character(best.results.df$ntax.vector.string), "_")[[1]]))
wait.times<-cbind(wait.times, ntax.dev=wait.times$ntax/250000)
wait.times$ntax.dev[which(wait.times$ntax.dev>1)]<-1/wait.times$ntax.dev[which(wait.times$ntax.dev>1)]

library(RColorBrewer)
mypalette<-brewer.pal(8,"Dark2")
par(mfcol=c(2,1))
combo.namesNot1<-combo.names[-1]

plot(x=c(-20,136),y=c(0,8),type="n", xlab="MY from root", yaxt="n", xaxt="n",ylab="", bty="n", main="unweighted")
axis(side=1, at=c(0, 30, 60, 90, 120, 136))
scaling<-10
for (i in 1:7) {
	density.val<-density(wait.times[,i], from=0, to=136, adjust=.7)
	lines(density.val$x, i+scaling*density.val$y, col=mypalette[i])
	lines(density.val$x, i-scaling*density.val$y, col=mypalette[i])
	focus<-quantile(wait.times[,i], 0.5)
	lines(rep(focus,2), i+c(-1, 1) * approx(density.val$x, y = scaling*density.val$y, xout=focus)$y, lwd=3)
	focus<-quantile(wait.times[,i], 0.025)
	lines(rep(focus,2), i+c(-1, 1) * approx(density.val$x, y = scaling*density.val$y, xout=focus)$y, lwd=2)
	focus<-quantile(wait.times[,i], 0.975)
	lines(rep(focus,2), i+c(-1, 1) * approx(density.val$x, y = scaling*density.val$y, xout=focus)$y, lwd=2)

	points(wait.times[,i],rep(i,length(wait.times[,i])), col=mypalette[i], pch=".")
	text(x=-10, y=i, combo.namesNot1[i], col=mypalette[i])
}


plot(x=c(-20,136),y=c(0,8),type="n", xlab="MY from root", yaxt="n", xaxt="n",ylab="", bty="n", main="weighted")
axis(side=1, at=c(0, 30, 60, 90, 120, 136))
scaling<-10
for (i in 1:7) {
	samples<-apply(rmultinom(1000, size=length(wait.times[,i]), prob=wait.times$ntax.dev),1, sum)
	sampled.times<-c()
	for (run in sequence(length(samples))) {
		sampled.times<-c(sampled.times, rep(wait.times[run,i], samples[run]))	
	}
	density.val<-density(sampled.times, from=0, to=136)
	lines(density.val$x, i+scaling*density.val$y, col=mypalette[i])
	lines(density.val$x, i-scaling*density.val$y, col=mypalette[i])
	focus<-quantile(sampled.times, 0.5)
	lines(rep(focus,2), i+c(-1, 1) * approx(density.val$x, y = scaling*density.val$y, xout=focus)$y, lwd=3)
	focus<-quantile(sampled.times, 0.025)
	lines(rep(focus,2), i+c(-1, 1) * approx(density.val$x, y = scaling*density.val$y, xout=focus)$y, lwd=2)
	focus<-quantile(sampled.times, 0.975)
	lines(rep(focus,2), i+c(-1, 1) * approx(density.val$x, y = scaling*density.val$y, xout=focus)$y, lwd=2)

	symbols(wait.times[,i],rep(i,length(wait.times[,i])), circles=samples/(max(samples)), bg=mypalette[i], fg=NA, add=TRUE, inches=FALSE)
	text(x=-10, y=i, combo.namesNot1[i], col=mypalette[i])
}


# par(mfrow=c(2,2))
# plot(x=range(c(results.df$ntax.rescale, results.df$ntax.actual)), y=c(0,1), type="n", bty="n", ylab="proportion 0x00xx", xlab="ntax actual", log="x")
# text(x=results.df[which(results.df$model=="divonly"),]$ntax.actual, y=results.df[which(results.df$model=="divonly"),]$proport_0x00xx, labels="D", col="blue")
# text(x=results.df[which(results.df$model=="transonly"),]$ntax.actual, y=results.df[which(results.df$model=="transonly"),]$proport_0x00xx, labels="T", col="red")
# text(x=results.df[which(results.df$model=="full"),]$ntax.actual, y=results.df[which(results.df$model=="full"),]$proport_0x00xx, labels="f", col="orchid")
# abline(v=250000, lty="dotted")
# 
# plot(x=-range(c(results.df$ntax.rescale, results.df$ntax.actual)), y=c(0,1), type="n", bty="n", ylab="proportion 0x00xx", xlab="neg ntax rescale param")
# text(x=-results.df[which(results.df$model=="divonly"),]$ntax.rescale, y=results.df[which(results.df$model=="divonly"),]$proport_0x00xx, labels="D", col="blue")
# text(x=-results.df[which(results.df$model=="transonly"),]$ntax.rescale, y=results.df[which(results.df$model=="transonly"),]$proport_0x00xx, labels="T", col="red")
# text(x=-results.df[which(results.df$model=="full"),]$ntax.rescale, y=results.df[which(results.df$model=="full"),]$proport_0x00xx, labels="f", col="orchid")
# 
# 
# plot(x=range(c(results.df$ntax.rescale, results.df$ntax.actual)), y=c(0,1), type="n", bty="n", ylab="proportion 0x11xx", xlab="ntax actual", log="x")
# text(x=results.df[which(results.df$model=="divonly"),]$ntax.actual, y=results.df[which(results.df$model=="divonly"),]$proport_0x11xx, labels="D", col="blue")
# text(x=results.df[which(results.df$model=="transonly"),]$ntax.actual, y=results.df[which(results.df$model=="transonly"),]$proport_0x11xx, labels="T", col="red")
# text(x=results.df[which(results.df$model=="full"),]$ntax.actual, y=results.df[which(results.df$model=="full"),]$proport_0x11xx, labels="f", col="orchid")
# abline(v=250000, lty="dotted")
# 
# plot(x=-range(c(results.df$ntax.rescale, results.df$ntax.actual)), y=c(0,1), type="n", bty="n", ylab="proportion 0x11xx", xlab="neg ntax rescale param")
# text(x=-results.df[which(results.df$model=="divonly"),]$ntax.rescale, y=results.df[which(results.df$model=="divonly"),]$proport_0x11xx, labels="D", col="blue")
# text(x=-results.df[which(results.df$model=="transonly"),]$ntax.rescale, y=results.df[which(results.df$model=="transonly"),]$proport_0x11xx, labels="T", col="red")
# text(x=-results.df[which(results.df$model=="full"),]$ntax.rescale, y=results.df[which(results.df$model=="full"),]$proport_0x11xx, labels="f", col="orchid")
# 
# par(mfrow=c(2,1))
# plot(x=range(c(results.df$ntax.rescale)), y=range(c(results.df$div_0x00xx, results.df$div_0x11xx)), xlab="ntax rescale param", ylab="div. rate", log="xy", bty="n", type="n")
# text(x=results.df[which(results.df$model=="divonly"),]$ntax.rescale, y=results.df[which(results.df$model=="divonly"),]$div_0x00xx, labels="D", col="blue")
# text(x=results.df[which(results.df$model=="divonly"),]$ntax.rescale, y=results.df[which(results.df$model=="divonly"),]$div_0x11xx, labels="d", col="blue")
# text(x=results.df[which(results.df$model=="transonly"),]$ntax.rescale, y=results.df[which(results.df$model=="transonly"),]$div_0x00xx, labels="T", col="red")
# text(x=results.df[which(results.df$model=="transonly"),]$ntax.rescale, y=results.df[which(results.df$model=="transonly"),]$div_0x11xx, labels="t", col="red")
# text(x=results.df[which(results.df$model=="full"),]$ntax.rescale, y=results.df[which(results.df$model=="full"),]$div_0x00xx, labels="F", col="orchid")
# text(x=results.df[which(results.df$model=="full"),]$ntax.rescale, y=results.df[which(results.df$model=="full"),]$div_0x11xx, labels="f", col="orchid")
# 
# 
# q.indices<-which(grepl("q", colnames(results.df)))
# plot(x=range(results.df[,q.indices]), y=range(0:5), type="n", bty="n", yaxt="n", ylab="", xlab="q")
# unique.plot<-function(df, col, y) {
#   vals<-unique(df)
#   names(vals)<-NULL
#   vals<-unique(as.numeric(vals))
#   val.counts<-rep(NA, length(vals))
#   for (val.index in sequence(length(vals))) {
#     val.counts[val.index]<-sum(df==vals[val.index])
#   }
#   zero.index<-which(vals==0)
#   vals<-vals[-zero.index]
#   val.counts<-val.counts[-zero.index]
#   val.counts<-val.counts/sum(val.counts)
#   for (val.index in sequence(length(vals))) {
#     symbols(x=vals[val.index], y=y, circles=min(val.counts[val.index]*0.0002, 0.00005), fg=col, bg=col, add=TRUE, inches=FALSE)
#   }
# }
# unique.plot(results.df[which(results.df$model=="divonly"),q.indices], "blue", 1)
# unique.plot(results.df[which(results.df$model=="transonly"),q.indices], "red", 2)
# unique.plot(results.df[which(results.df$model=="full"),q.indices], "orchid", 3)
