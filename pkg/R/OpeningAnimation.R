setwd("/Users/bomeara/Documents/MyDocuments/Ongoing/Talks/2014vi24_Evolution_Stebbins")
#library(extrafont)
library(animation)
#font_import()
#loadfonts()
lead <- c("Brian O'Meara", "Stacey Smith")
middle <- c("Scott Armbruster", "Lawrence Harder", "Chris Hardy", "Lena Hileman", "Larry Hufford", "Amy Litt", "Susana Magallon", "Stephen Smith", "Peter Stevens")
senior <- c("Charles Fenster", "Pamela Diggle")
#saveVideo({
#ani.options(nmax = 30)
#for (i in sequence(30)) {
#plot(x=c(.5,2.5), y=c(1.5,8.5), bty="n", type="n", xlab="", xaxt="n", ylab="", yaxt="n", mar=c(0,0,0,0))
#rect(0.5,7.5,2.5,8.5, col="lightyellow", border=NA)
#rect(0.5,1.5,2.5,2.5, col="lightyellow", border=NA)
#text(1.5, 7.5, "Equal lead authors", family="serif", cex=0.9, pos=3)
#text(c(1,2), c(8, 8), lead[1+c((i+1)%%2, i%%2)], family="serif", cex=2.4)
#text(rep(c(1,2),5), rep(c(7:3),2), c(sample(middle, size=9, replace=FALSE), ""),  family="serif", cex=2.4)
#text(c(1,2), c(2, 2), senior[1+c((i+1)%%2, i%%2)], family="serif", cex=2.4)
#text(1.5, 1.5, "Equal senior authors", family="serif", cex=0.9, pos=3)

#}
#}, interval = 3, video.name = "authorsWhiteBG.mp4", ani.width = 700, ani.height = 700)

saveVideo({
ani.options(nmax = 30)
for (i in sequence(10)) {
plot(x=c(.5,2.5), y=c(1.5,8.5), bty="n", type="n", xlab="", xaxt="n", ylab="", yaxt="n", mar=c(0,0,0,0))
rect(-1, -1, 10, 10, col="black", border=NA)
rect(0.5,7.5,2.5,8.5, col="gray26", border=NA)
rect(0.5,1.5,2.5,2.5, col="gray26", border=NA)
text(1.5, 7.5, "Equal lead authors", family="serif", cex=1.0, pos=3, col="white")
text(c(1,2), c(8, 8), (lead[1+c((i+1)%%2, i%%2)]), family="serif", cex=2.4, col="white")
text(rep(c(1,2),5), rep(c(7:3),2), (c(sample(middle, size=9, replace=FALSE), "")),  family="serif", cex=2.4, col="white")
text(c(1,2), c(2, 2), (senior[1+c((i+1)%%2, i%%2)]), family="serif", cex=2.4, col="white")
text(1.5, 1.5, "Equal senior authors", family="serif", cex=1.0, pos=3, col="white")

}
}, interval = 3, video.name = "authorsBlackBG.mp4", ani.width = 700, ani.height = 700)

system(paste("cp ", tempdir(), "/authors* .", sep=""))