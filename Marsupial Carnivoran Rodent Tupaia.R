###04/17/2013

##DFA on osteological indices(<0.001)

setwd("/Users/Big Bear/Documents/Projects/Doctoral Dissertation/Chapter 1/Statistical Analysis/LinearDiscriminantAnalysis/Sig0001")

## library used for the plots
library(MASS)

##Reading the data
mydata.new <- read.table("RatioMasterW_oNaVaule_1Species_Sorted_Select_Indices_P0001_18April2013.txt", head=T, sep=",")

## Discriminant Functional Analysis on Significant Indices(<0.001)
sigdata<-mydata.new[,4:34]
rowname<-mydata.new[,1]

Strudata<-sigdata[,2:31]

siglda.mod <- lda(Locomotion ~ ., data = sigdata,
									prior=c(0.125,0.125,0.125,0.125,
													0.125,0.125,0.125,0.125))

## Testing for significant differences between groups
sigindex <- as.matrix(sigdata[,c(2:31)])
sigoutput <- manova(sigindex~sigdata$Locomotion)
sigsum.wilks<-summary(sigoutput, test="Wilks")
sigsum.avo<-summary.aov(sigoutput)

siglda.res <- predict(siglda.mod, sigdata)

## First two coefficient of the LDA of the osteological indices
coemtr<-cbind(siglda.mod$scaling[,1],siglda.mod$scaling[,2])

## Calculate the centroid of each locomotion

## Arboreal
ar.ld1=sum(siglda.res$x[,1]*(sigdata$Locomotion=="Arboreal"))/
	sum(sigdata$Locomotion=="Arboreal")
ar.ld2=sum(siglda.res$x[,2]*(sigdata$Locomotion=="Arboreal"))/
	sum(sigdata$Locomotion=="Arboreal")
ar.ld3=sum(siglda.res$x[,3]*(sigdata$Locomotion=="Arboreal"))/
	sum(sigdata$Locomotion=="Arboreal")

## Fossorial
f.ld1=sum(siglda.res$x[,1]*(sigdata$Locomotion=="Fossorial"))/
	sum(sigdata$Locomotion=="Fossorial")
f.ld2=sum(siglda.res$x[,2]*(sigdata$Locomotion=="Fossorial"))/
	sum(sigdata$Locomotion=="Fossorial")
f.ld3=sum(siglda.res$x[,3]*(sigdata$Locomotion=="Fossorial"))/
	sum(sigdata$Locomotion=="Fossorial")

## Gliding
g.ld1=sum(siglda.res$x[,1]*(sigdata$Locomotion=="Gliding"))/
	sum(sigdata$Locomotion=="Gliding")
g.ld2=sum(siglda.res$x[,2]*(sigdata$Locomotion=="Gliding"))/
	sum(sigdata$Locomotion=="Gliding")
g.ld3=sum(siglda.res$x[,3]*(sigdata$Locomotion=="Gliding"))/
	sum(sigdata$Locomotion=="Gliding")

## Saltatorial
s.ld1=sum(siglda.res$x[,1]*(sigdata$Locomotion=="Saltatorial"))/
	sum(sigdata$Locomotion=="Saltatorial")
s.ld2=sum(siglda.res$x[,2]*(sigdata$Locomotion=="Saltatorial"))/
	sum(sigdata$Locomotion=="Saltatorial")
s.ld3=sum(siglda.res$x[,3]*(sigdata$Locomotion=="Saltatorial"))/
	sum(sigdata$Locomotion=="Saltatorial")

## Scansorial
sc.ld1=sum(siglda.res$x[,1]*(sigdata$Locomotion=="Scansorial"))/
	sum(sigdata$Locomotion=="Scansorial")
sc.ld2=sum(siglda.res$x[,2]*(sigdata$Locomotion=="Scansorial"))/
	sum(sigdata$Locomotion=="Scansorial")
sc.ld3=sum(siglda.res$x[,3]*(sigdata$Locomotion=="Scansorial"))/
	sum(sigdata$Locomotion=="Scansorial")

## Semi-aquatic
sa.ld1=sum(siglda.res$x[,1]*(sigdata$Locomotion=="Semi-aquatic"))/
	sum(sigdata$Locomotion=="Semi-aquatic")
sa.ld2=sum(siglda.res$x[,2]*(sigdata$Locomotion=="Semi-aquatic"))/
	sum(sigdata$Locomotion=="Semi-aquatic")
sa.ld3=sum(siglda.res$x[,3]*(sigdata$Locomotion=="Semi-aquatic"))/
	sum(sigdata$Locomotion=="Semi-aquatic")

## Semi-fossorial
sf.ld1=sum(siglda.res$x[,1]*(sigdata$Locomotion=="Semi-fossorial"))/
	sum(sigdata$Locomotion=="Semi-fossorial")
sf.ld2=sum(siglda.res$x[,2]*(sigdata$Locomotion=="Semi-fossorial"))/
	sum(sigdata$Locomotion=="Semi-fossorial")
sf.ld3=sum(siglda.res$x[,3]*(sigdata$Locomotion=="Semi-fossorial"))/
	sum(sigdata$Locomotion=="Semi-fossorial")


## Terrestrial
t.ld1=sum(siglda.res$x[,1]*(sigdata$Locomotion=="Terrestrial"))/
	sum(sigdata$Locomotion=="Terrestrial")
t.ld2=sum(siglda.res$x[,2]*(sigdata$Locomotion=="Terrestrial"))/
	sum(sigdata$Locomotion=="Terrestrial")
t.ld3=sum(siglda.res$x[,3]*(sigdata$Locomotion=="Terrestrial"))/
	sum(sigdata$Locomotion=="Terrestrial")

centroid1<-matrix(c(ar.ld1,ar.ld2,f.ld1,f.ld2,g.ld1,g.ld2,
										s.ld1,s.ld2,sc.ld1,sc.ld2,sa.ld1,sa.ld2,
										sf.ld1,sf.ld2,t.ld1,t.ld2),nrow=8,ncol=2,byrow=TRUE)

centroid2<-matrix(c(ar.ld1,ar.ld3,f.ld1,f.ld3,g.ld1,g.ld3,
										s.ld1,s.ld3,sc.ld1,sc.ld3,sa.ld1,sa.ld3,
										sf.ld1,sf.ld3,t.ld1,t.ld3),nrow=8,ncol=2,byrow=TRUE)


cols=c(rep("firebrick1",30),          ## Arboreal
			 rep("grey36",12),              ## Fossorial
			 rep("dodgerblue3",3),          ## Gliding
			 rep("magenta",5),              ## Saltatorial
			 rep("green3",9),               ## Scansorial
			 rep("deepskyblue",9),          ## Semi-aquatic
			 rep("yellow3",9),              ## Semi-fossorial
			 rep("darkorange",30))          ## Terrestrial


sybs=c(rep(1,30),          ## Arboreal
			 rep(15,12),              ## Fossorial
			 rep(2,3),          ## Gliding
			 rep(4,5),              ## Saltatorial
			 rep(5,9),               ## Scansorial
			 rep(18,9),          ## Semiaquatic
			 rep(19,9),              ## Semifossorial
			 rep(0,30))          ## Terrestrial

cols.centroid=c(
	rep("firebrick1",1),          ## Arboreal
	rep("grey36",1),              ## Fossorial
	rep("dodgerblue3",1),          ## Gliding
	rep("magenta",1),              ## Saltatorial
	rep("green3",1),               ## Scansorial
	rep("deepskyblue",1),          ## Semi-aquatic
	rep("yellow3",1),              ## Semi-fossorial
	rep("darkorange",1))          ## Terrestrial


## set up plot area
par(mfrow=c(2,2),oma=c(2,2,0,0),
		mar=c(2,2,0,0),mgp=c(3,0.5,0),
		tck=-0.005)

## Plot all data for Marsupialia
plot(siglda.res$x[,1],
		 siglda.res$x[,2],
		 ann=F,
		 pch=sybs,
		 col="black",
		 las=1,
		 frame=F,
		 cex=1.5,
		 xaxt="none")
points(centroid1, 
			 pch=8,
			 col="black",
			 cex=2)
abline(v=0,lty=3)
abline(h=0,lty=3)
## text the x and y axes
mtext(outer=F, side=2, line=2, text="Canonical Function 2 (25.83%)",cex=0.7)
mtext(outer=F, side=3, lin=-1, text="Marsupialia", cex=0.8,at=6)

## Plot all data for Carnivora
plot(siglda.res$x[,1],
		 siglda.res$x[,2],
		 ann=F,
		 axes=F,
		 pch=sybs,
		 col="black",
		 las=1,
		 frame=F,
		 cex=1.5)
points(centroid1, 
			 pch=8,
			 col="black",
			 cex=2)
abline(v=0,lty=3)
abline(h=0,lty=3)
## text the x and y axes
mtext(outer=F, side=3, lin=-1, text="Carnivora", cex=0.8,at=6)

## Plot all data for Rodentia
plot(siglda.res$x[,1],
		 siglda.res$x[,2],
		 ann=F,
		 pch=sybs,
		 col="black",
		 las=1,
		 frame=F,cex=1.5)
points(centroid1, 
			 pch=8,
			 col="black",
			 cex=2)
abline(v=0,lty=3)
abline(h=0,lty=3)
## text the x and y axes
mtext(outer=F, side=1, line=2, text="Canonical Function 1 (49.13%)",cex=0.7)
mtext(outer=F, side=2, line=2, text="Canonical Function 2 (25.83%)",cex=0.7)
mtext(outer=F, side=3, lin=-1, text="Rodentia", cex=0.8,at=6)

## Plot all data for Tupaia
plot(siglda.res$x[,1],
		 siglda.res$x[,2],
		 ann=F,
		 pch=sybs,
		 col="black",
		 las=1,
		 frame=F,
		 cex=1.5,
		 yaxt="none")
points(centroid1, 
			 pch=8,
			 col="black",
			 cex=2)
abline(v=0,lty=3)
abline(h=0,lty=3)
## text the x and y axes
mtext(outer=F, side=1, line=2, text="Canonical Function 1 (49.13%)",cex=0.7)
mtext(outer=F, side=3, lin=-1, text="Tupaia",cex=0.8,at=6,font=3)
