

ntax.vector=c(2^4, 2^5, 2^6, 2^7, 2^8)
shape.vector=c("balanced", "right")
rate.vector=c(0.1, 0.5, 1, 2, 10)
change.position.vector=c("root", "quarter", "cherry")
reps.per.combination=1





for(rep in 1:reps.per.combination){
	for(ntax.index in 1:length(ntax.vector)){
		ntax=ntax.vector[ntax.index]
	
		for(shape.index in 1:length(shape.vector)){
			shape=shape.vector[shape.index]
		
			for(rate.index in 1:length(rate.vector)){
				rate=rate.vector[rate.index]
				
				for(change.position.index in 1:length(change.position.vector)){
					change=change.position.vector[change.position.index]
					
					fileNameRoot<-paste("ntax.",ntax,".shape.",shape,".rate.",rate,".change.",change,".rep.",rep,sep="",collapse="")
					batchFileName<-paste(fileNameRoot,".R",sep="",collapse="")
					cat("library(geiger)\n","library(RBrownie)\n",sep="",file=batchFileName,append=FALSE)
					cat("ntax=",ntax,"\nshape='",shape,"'\nrate=",rate,"\nchange='",change,"'\nrep=",rep,"\n",sep="",file=batchFileName,append=TRUE)
					cat("source('deltaBrownie.R')","\n",sep="",file=batchFileName,append=TRUE)
					cat("streeBrlen(ntax, type=shape)->phy\n",sep="",file=batchFileName,append=TRUE)
					cat("blMultiplier(phy, rate, change, shape )->phy1\n",sep="",file=batchFileName,append=TRUE )
					cat("model.matrix<-matrix(1)","\n", sep="",file=batchFileName,append=TRUE )
					cat("data<-sim.char(phy1, model.matrix, nsims = 1, model = 'brownian', root.state = 1 )", "\n",sep="",file=batchFileName,append=TRUE)
					cat("data<-c(data[,,1])", "\n",sep="",file=batchFileName,append=TRUE)
					cat("iterateNonCensored(phy, data)->results.",fileNameRoot, "\n",sep="",file=batchFileName,append=TRUE)
					cat("save(results.",fileNameRoot, ", file=paste('results.",fileNameRoot,".data.frame'",',sep=""), ascii=TRUE, compress=FALSE)\n', sep="",file= batchFileName, append=TRUE) 
										
					#for plotting the pdf:
					
					cat("pdf('phy1", fileNameRoot,".pdf')", "\n", sep="",file=batchFileName,append=TRUE)
					cat("plot(phy1)", "\n",sep="",file=batchFileName,append=TRUE)
					cat("title(main='", fileNameRoot, "', col.main='red', font.main=2)", "\n",sep="",file=batchFileName,append=TRUE)
					cat("dev.off()", "\n",sep="",file=batchFileName,append=TRUE)
					
					#TO DO: KEEP ADDING COMMANDS FOR THE BATCH FILE 
						#where is the phy that streeBrlen is outputting with each file (i.e. ntax256shaperightrate10changecherryrep8.R), this object needs to be fed in the write() below marked "HERE!"
						
						
						
						
					#MAKE SURE EACH DATA.FRAME IS SAVED AS AN R OBJECT IN ANOTHER FILE DONW
					
					#TO DO: SYSTEM(EZSUB R CMD BATCH BATCHFILENAME)  DONE
					#intern=TRUE captures the command as an R character vector; all other TRUE/FALSE arguements are left with their defaults, except for the command argument wich pastes the directory
					lsString=paste("ls -l results.", fileNameRoot,".data.frame",'|grep -c data.frame', sep="", collapse="")
					outputCount=suppressWarnings(as.numeric(system(lsString, intern=TRUE)))
					if(outputCount==0){
						system(command=paste("condor_submitR ", batchFileName, " deltaBrownie.R", sep=""), intern=TRUE, ignore.stderr=FALSE, wait=TRUE, input=NULL, show.output.on.console=FALSE)
					}
					Sys.sleep(5)
					
				}
			}
		}
	}	
}





