data<-read.csv("~/Dropbox/ntax.Allfilter.results.csv")
data<-data[which(data$nruns.with.predictions>2),]
#plot(data$ntax.rescale, data$ntax.predicted, log="xy")
x.points<-seq(from=6, to=13, length.out=50)

y<-log(data$ntax.predicted)
x<-log(data$ntax.rescale)
plot(x,y, ylim=c(0,25), xlim=range(x.points))
points(x, y=log(data$ntax.actual), pch="-")
points(x, y=log(data$ntax.partials), pch="x")
text(x, y=1+y, data$nruns, col="blue")
lines(x[order(x)], log(data$ntax.predicted.025)[order(x)], lty="dotted")
lines(x[order(x)], log(data$ntax.predicted.975)[order(x)], lty="dotted")
lines(x[order(x)], log(data$ntax.partials.025)[order(x)], lty="dashed")
lines(x[order(x)], log(data$ntax.partials.975)[order(x)], lty="dashed")

#abline(lm(y~x))
regression1<-lm(y~x)
lines(x.points, (predict(regression1, data.frame(x=x.points))), lwd=4)
#data<-data[which(data$nruns>10),]
abline(h=log(250000))
data.full<-data[which(data$nruns.with.predictions - data$nruns <6),]
y<-log(data.full$ntax.predicted)
x<-log(data.full$ntax.rescale)
points(x, y=log(data.full$ntax.actual), col="red")
regression.full<-lm(y ~ x)
lines(x.points, (predict(regression.full, data.frame(x=x.points))), lwd=1, col="red")

y<-log(data$ntax.partials)
x<-log(data$ntax.rescale)
regression.partials<-lm(y ~ x)
lines(x.points, (predict(regression.partials, data.frame(x=x.points))), lwd=1, col="purple")

ToOptimize<-function(x, regression, goal) {
  return(abs(goal - predict(regression, data.frame(x=x)))) 
}

optimize.result<-optimize(f=ToOptimize, interval=c(5,9), regression=regression.partials, goal=log(250000))
abline(v=optimize.result$minimum, col="purple")
print(paste("optimal scaling partials is ",round(exp(optimize.result$minimum))))

optimize.result<-optimize(f=ToOptimize, interval=c(5,11), regression= regression1, goal=log(250000))
abline(v=optimize.result$minimum, col="black")
print(paste("optimal scaling all is ",round(exp(optimize.result$minimum))))

optimize.local<-data$ntax.rescale[which.min(abs(data$ntax.predicted-250000))]

print(paste("optimal scaling based on local is ", optimize.local))
abline(v=log(optimize.local), col="green")