library(compiler)
enableJIT(3)
tf=136
t.stop=136+15
while(1<2) { #keep looping
  
  setwd("/Users/bomeara/Documents/MyDocuments/Active/FloralAssembly/SimsMay2013")
  system("cp /Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral/V7*R .")
  source("V7_NewSimulator.R")
  files<-system("ls */*RSave", intern=TRUE)
  if(length(files)>0) {
    files.total<-rep(0, length(files))
    for (i in sequence(length(files))) {
      total.time<-0
      try(total.time<-as.numeric(regmatches(files[i],regexec("total(\\d+)", files[i]))[[1]][2]))
      if (!is.na(total.time)) {
        files.total[i]<-total.time
      }
    }
    files<-files[which(files.total<t.stop)] #don't want to keep running forever.
    files.total<-files.total[which(files.total<t.stop)]
    files<-files[order(files.total)] #do the ones that have run the least amount of time first
    files.total<-files[order(files.total)]
    print(cbind(files, files.total[order(files.total)]))
    Sys.sleep(30) #give time to finish transferring any files
    for (i in sequence(length(files))) {
      system(paste("svn add ", files[i]))
      dir<-strsplit(files[i], "/")[[1]][1]
      file<-strsplit(files[i], "/")[[1]][2]
      setwd(dir)
      system("cp /Users/bomeara/Documents/MyDocuments/Active/OMearaLabR/pkg/R/floral/V7*R .")
      file.root<-strsplit(gsub("_app\\d+total\\d+","", file), "\\.")[[1]][1]
      
      system(paste("svn mv ", file, " ", file, ".staged", sep=""))
      system("svn commit -m'staging'")
      t.additional<-1
      if(files.total[i]<100) {
        t.additional<-5
      }
      
      if(files.total[i]<80) {
        t.additional<-10 
      }
      
      cat(paste('source("V7_UtilityFns.R")
load("',file,'.staged")
source("V7_NewSimulator.R")
source("V7_StochasticSSASims_Functions.R")
appendParallelSSA(x0=x0, q.means=q.means, lambda.means=lambda.means, mu.means=mu.means, prev.history=history, t.additional=', t.additional, ', maxWallTime=Inf, verbose=F, file.string="', file.root, '", print.freq=100, rescale.species=rescale.species, yule.scale=0, t.rescale=t.rescale, x0.rescale=x0.rescale)
', sep=""), file=paste(file.root, ".append.R", sep=""), append=FALSE)
      
      
      
      cat(paste('#! /bin/sh


echo "I am process id $$ on" `hostname`
/usr/bin/R CMD BATCH --no-save ', paste(file.root, ".append.R", sep=""), ' ', paste(file.root, ".append.R", sep=""), '.$$.`hostname`.Rout
ls
cat *Rout
pwd
echo "Now finished"
', sep=""), file=paste(file.root, ".myProg", sep=""))
      
      
      system(paste("chmod u+x ", paste(file.root, ".myProg", sep=""), sep=""))
      
      
      cat(paste('executable=', paste(file.root, ".myProg", sep=""), '
universe=vanilla
arguments=', paste(file.root, ".append.R", sep=""), '$(Cluster).$(Process) 1
requirements = Memory >= 0
output=results.output.$(Process)
error=results.error.$(Process)
transfer_input_files=', paste(file.root, ".append.R", sep=""),',V7_NewSimulator.R,V7_StochasticSSASims_Functions.R,V7_UtilityFns.R,', paste(file, ".staged", sep=""),'
log=results.log
notification=never
should_transfer_files=YES
when_to_transfer_output = ON_EXIT_OR_EVICT
queue 1

', sep=""), file=paste(file.root, ".submit", sep=""))
      
      
      system(paste("/condor/condor-installed/bin//condor_submit ", paste(file.root, ".submit", sep=""), sep=""))
      
      
      setwd("..")
    }
    system("svn commit -m'more runs'")
  }
  Sys.sleep(10)
}
