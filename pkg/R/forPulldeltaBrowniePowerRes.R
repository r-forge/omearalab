#setwd("/data/cichlids/")
counter<-0

ntax.vector=c(2^4, 2^5, 2^6, 2^7, 2^8)
shape.vector=c("balanced", "right")
rate.vector=c(0.1, 0.5, 1, 2, 10)
change.position.vector=c("root", "quarter", "cherry")
reps.per.combination=10

numberRuns=length(ntax.vector)*length(shape.vector)*length(rate.vector)*length(change.position.vector)*(reps.per.combination)

all.a<-vector("list", numberRuns)



for(rep in 1:reps.per.combination){
	for(ntax.index in 1:length(ntax.vector)){
		ntax=ntax.vector[ntax.index]
	
		for(shape.index in 1:length(shape.vector)){
			shape=shape.vector[shape.index]
		
			for(rate.index in 1:length(rate.vector)){
				rate=rate.vector[rate.index]
				
				for(change.position.index in 1:length(change.position.vector)){
					change=change.position.vector[change.position.index]
					
					counter<-counter+1
					
					#results.ntax.128.shape.balanced.rate.0.5.change.root.rep.1.data.frame
					
					fileNameRoot<-paste("ntax.",ntax,".shape.",shape,".rate.",rate,".change.",change,".rep.",rep,sep="",collapse="")
					names(all.a[[counter]])<-c(fileNameRoot) #paste(fileNameRoot)
					
					if(system(paste("ls results* | grep -c results.",fileNameRoot,".data.frame",sep="",collapse=""),intern=TRUE)==1) {#list all files 
						load(paste("results.",fileNameRoot,".data.frame",sep="",collapse=""))
						all.a$fileNameRoot<-paste("results.", fileNameRoot, sep="")

					}
				}
			}
		}		
	}
}




