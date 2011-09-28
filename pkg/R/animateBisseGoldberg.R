library(animation)
rm(list=ls())
#0=compatible
#1=incompatible
bisse.ani<-function(q0=0.0,q1=0.555,b0=5.199,b1=2.933,d0=5.433,d1=2.319,starting0=0,starting1=1,colorVector=c(rgb(0,1,.2,0.9),rgb(0.3,0,1,0.7),rgb(0,1,.2,0.2),rgb(0.3,0,1,0.2)),radius=0.1,drawtomutatecount=2,allowExtinction=FALSE,subsampleAt=100) {
	colorMatrix<-matrix(nrow=2,ncol=2) #first row is colors for state 0 and 1 living, second is colors for state 0 and 1, dead
	colorMatrix[1,]=colorVector[1:2]
	colorMatrix[2,]=colorVector[3:4]
	qrate<-c(q0,q1)/10 #rescaling for better animation
	qcounts<-c(0,0)
	brate<-c(b0,b1)/10
	bcounts<-c(0,0)
	drate<-c(d0,d1)/10
	dcounts<-c(0,0)
	labelText<-c("c","i")
	drawtomutatecount=round(drawtomutatecount)
	makePoint<-function(x=0,y=0,state=0,pointAge=0,stateAge=0,extant=1) {
		pointVector<-c(x,y,x,y,state,pointAge,stateAge,extant)
		return(pointVector)
	}
	movePoint<-function(point, pointsMatrix=pointsMatrix,flee.scaling=1) {
		x.velocity=0
		y.velocity=0
		
		#flee other points
		 x.me<-pointsMatrix[point,1]
		 y.me<-pointsMatrix[point,2]
		 for (otherpoint in 1:dim(pointsMatrix)[1]) {
			 if (point!=otherpoint) {
				 x.other<-pointsMatrix[otherpoint,1]
				 y.other<-pointsMatrix[otherpoint,2]
				 distance<-min(1000,sqrt(((x.other-x.me)^2) + ((y.other-y.me)^2)),na.rm=TRUE)
				 if (distance<=radius*2) {
				 	x.velocity<-x.velocity+(x.me-x.other)/distance
				 	y.velocity<-y.velocity+(y.me-y.other)/distance
				 }
			 }
		 }
	 
		 if (pointsMatrix[point,8]==1) { 		#go towards optimum if alive
			 x.me<-pointsMatrix[point,1]
			 y.me<-pointsMatrix[point,2]
			 x.goal<- -2 #goal if state=0
			 if (pointsMatrix[point,5]==1) {
				 x.goal<-2 #goal if state=1
			 }
			 y.goal=0
			 distance<-sqrt(((x.me-x.goal)^2) + ((y.me-y.goal)^2))
			 x.velocity<-x.velocity-(x.me-x.goal)*distance
			 y.velocity<-y.velocity-(y.me-y.goal)*distance
		 }
		 else {    #drift towards bottom if dead
			 y.me<-max(-5,pointsMatrix[point,2],na.rm=TRUE) 
			 y.goal<- -4 #bottom
			 if (y.me>y.goal) {
				 y.velocity<-y.velocity - 1
			 }
		 }
		 pointsMatrix[point,3]<-x.velocity*flee.scaling
		 pointsMatrix[point,4]<-y.velocity*flee.scaling
		 return(pointsMatrix)
 	}
	
	moveAllPoints<-function(pointsMatrix=pointsMatrix,velocity.max=0.2) {
		for (point in 1:dim(pointsMatrix)[1]) {
			pointsMatrix<-movePoint(point,pointsMatrix)
		}
		#scalingFactor=velocity.max/max(abs(c(pointsMatrix[,3],pointsMatrix[,4])))
		for (point in 1:dim(pointsMatrix)[1]) {
			pointsMatrix[point,3]<-sign(pointsMatrix[point,3])*min(abs(velocity.max),abs(pointsMatrix[point,3]))
			pointsMatrix[point,4]<-sign(pointsMatrix[point,4])*min(abs(velocity.max),abs(pointsMatrix[point,4]))
			pointsMatrix[point,1]<-pointsMatrix[point,3]+pointsMatrix[point,1]
			pointsMatrix[point,2]<-pointsMatrix[point,4]+pointsMatrix[point,2]
		}
		for (point in 1:dim(pointsMatrix)[1]) {
		      pointsMatrix[point,2]<-max(-4,pointsMatrix[point,2])
		      pointsMatrix[point,1]<-max(-4,pointsMatrix[point,1])
		      pointsMatrix[point,1]<-min(4,pointsMatrix[point,1])
		}
		return(pointsMatrix)
	}
	
	plotPoints<-function(pointsMatrix=pointsMatrix) {
		#print(colorMatrix)
		#print(pointsMatrix)
		plot(x=c(-4,4),y=c(-4,4),bty="n",xlab="",ylab="",type="n",xlim=c(-4,4),ylim=c(-4,4),xaxt="n",yaxt="n")
		for (point in 1:dim(pointsMatrix)[1]) {
			symbols(x=pointsMatrix[point,1],y=pointsMatrix[point,2],fg="black",bg=colorMatrix[abs(2-pointsMatrix[point,8]), 1+pointsMatrix[point,5] ], circles=radius,add=TRUE,inches=FALSE)
			text(x=pointsMatrix[point,1],y=pointsMatrix[point,2],labels=labelText[1+pointsMatrix[point,5] ],col="white")
		}
		text(x=-2,y=4,labels=c(paste("0=SC\nq01=",q0,"\nb0=",b0,"\nd0=",d0,sep="")),adj=c(0.5,1))
		#text(x=-1,y=4,labels=c(paste("0->1=",qcounts[1],"\n#birth0=",bcounts[1],"\n#death0=",dcounts[1],sep="")),adj=c(0.5,1))
		text(x=2,y=4,labels=c(paste("1=SI\nq10=",q1,"\nb1=",b1,"\nd1=",d1,sep="")),adj=c(0.5,1))
		#text(x=3,y=4,labels=c(paste("1->0=",qcounts[2],"\n#birth1=",bcounts[2],"\n#death1=",dcounts[2],sep="")),adj=c(0.5,1))

		if(max(pointsMatrix[,8])==0) {
			text(x=0,y=0,labels=c("Everything\nExtinct!"),adj=c(0.5,0.5),cex=2)
		}
	}
	
	mutatePoints<-function(pointsMatrix=pointsMatrix) {
		for (point in 1:dim(pointsMatrix)[1]) {
			if (pointsMatrix[point,8]==1) { #if point alive
				point.state<-pointsMatrix[point,5]
				point.reversestate<-0
				if (point.state==0) {
					point.reversestate<-1	
				}
				if (runif(1,0,1)<qrate[(point.state+1)]) {
					pointsMatrix[point,5]<-point.reversestate
					pointsMatrix[point,7]=0
					qcounts[(point.state+1)]<-qcounts[(point.state+1)]+1			
				}
			}
		}
		return(pointsMatrix)
	}
	
	birthdeathPoints<-function(pointsMatrix=pointsMatrix) {
		for (point in 1:dim(pointsMatrix)[1]) {
			if (pointsMatrix[point,8]==1) { #if point alive
				point.state<-pointsMatrix[point,5]
				randomdraw<-runif(1,0,1)
				#print(paste("randomdraw",randomdraw,"brate[(point.state+1)]",brate[(point.state+1)]))
				if (randomdraw<brate[(point.state+1)]) { #birth!
					pointsMatrix[point,6]=0 #new age
					pointsMatrix<-rbind(pointsMatrix,makePoint(x=(pointsMatrix[point,1]+runif(1,-0.1,.1)), y=(pointsMatrix[point,2]+runif(1,-0.1,.1)), state=point.state))
					bcounts[(point.state+1)]<-bcounts[(point.state+1)]+1
				}
				else {
					if (randomdraw<(brate[(point.state+1)]+drate[(point.state+1)])) {
						if (sum(pointsMatrix[,8])>1 || allowExtinction) {
							pointsMatrix[point,8]=0 #death
							dcounts[(point.state+1)]<-dcounts[(point.state+1)]+1
						}
					}
				}
			}
		}
		return(pointsMatrix)
	}
	
	i <- 1
	timesteps<-c()
	#make a matrix containing points
	pointsMatrix<-matrix(nrow=starting0+starting1,ncol=8)
	colnames(pointsMatrix)<-c("x_now","y_now","x_velocity","y_velocity","state","pointAge","stateAge","extant")
	for (point in 1:starting0) {
		pointsMatrix[point,]<-makePoint(x=runif(1,-2,0),y=runif(1,-1,1),state=0)
	}
	for (point in (starting0+1):(starting0+starting1)) {
		pointsMatrix[point,]<-makePoint(x=runif(1,0,2),y=runif(1,-1,1),state=1)
	}
	for (initialmove in 1:10) {
		pointsMatrix<-moveAllPoints(pointsMatrix)
	}
	
   while (i <= ani.options("nmax")) {
   #while(i<=10000) {
   		if (floor(i/drawtomutatecount)==ceiling(i/drawtomutatecount)) {
 			pointsMatrix<-mutatePoints(pointsMatrix)
 			pointsMatrix<-birthdeathPoints(pointsMatrix)
 		}
 		#print(dcounts)
       pointsMatrix<-moveAllPoints(pointsMatrix)
       if (dim(pointsMatrix)[1]>subsampleAt) {
       		pointsMatrix<-pointsMatrix[-1*seq(1,dim(pointsMatrix)[1],2),]
       }
      	print(i)
       plotPoints(pointsMatrix)
        Sys.sleep(ani.options("interval"))
        #Sys.sleep(1)
        i <- i + 1
    }
}
#par(ann=FALSE)
#ani.start(nmax=1000,interval=1,title="Diversification and Character Evolution")
#bisse.ani()
#ani.stop()
oopt=ani.options(nmax=500,interval=0.1)
#quartz()
#bisse.ani()
saveMovie(bisse.ani(),convert='convert',outdir="/Users/bomeara/Desktop/",movietype = "gif",width=600, height=600,interval=0.5)
ani.options(oopt)
#ani.stop()