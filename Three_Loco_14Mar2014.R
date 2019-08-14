###04/17/2013

##DFA on osteological indices(<0.001)

setwd("/Users/Big Bear/Documents/Projects/Doctoral Dissertation/Chapter 1/Statistical Analysis/LinearDiscriminantAnalysis/Sig0001")

##Libraries and packages
##Sourcing a package to color dendrograms
#source("http://bioinformatics.mdanderson.org/OOMPA/oompaLite.R")
#need to install all packages indenpendently


# library(vegan)
library(scatterplot3d)
library(rgl)
library(permute)
library(MASS)
library(calibrate)
library(vcd)
library(grid)
library(colorspace)
source("biostats.R")##helps with data screening, compiled by K. McGarigal at UMass

##Reading the data
mydata <- read.table("RatioMasterW_oNaVaule_1Species_Sorted_Select_Indices_three_P0001_18April2013.txt", head=TRUE, sep=",")

## Discriminant Functional Analysis on Significant Indices(<0.01)
sigdata<-mydata[,4:34]
rowname<-mydata[,1]

Strudata<-sigdata[,2:31]

siglda.mod <- lda(Locomotion ~ ., data = sigdata, prior=c(1/3,1/3,1/3))

## Testing for significant differences between groups
sigindex <- as.matrix(sigdata[,c(2:31)])
sigoutput <- manova(sigindex~sigdata$Locomotion)
sigsum.wilks<-summary(sigoutput, test="Wilks")
sigsum.avo<-summary.aov(sigoutput)

siglda.res <- predict(siglda.mod, sigdata)

## Calculate the centroid of each locomotion

## Arboreal
ar.ld1=sum(siglda.res$x[,1]*(sigdata$Locomotion=="Arboreal"))/
	sum(sigdata$Locomotion=="Arboreal")
ar.ld2=sum(siglda.res$x[,2]*(sigdata$Locomotion=="Arboreal"))/
	sum(sigdata$Locomotion=="Arboreal")


## Scansorial
sc.ld1=sum(siglda.res$x[,1]*(sigdata$Locomotion=="Scansorial"))/
	sum(sigdata$Locomotion=="Scansorial")
sc.ld2=sum(siglda.res$x[,2]*(sigdata$Locomotion=="Scansorial"))/
	sum(sigdata$Locomotion=="Scansorial")


## Terrestrial
t.ld1=sum(siglda.res$x[,1]*(sigdata$Locomotion=="Terrestrial"))/
	sum(sigdata$Locomotion=="Terrestrial")
t.ld2=sum(siglda.res$x[,2]*(sigdata$Locomotion=="Terrestrial"))/
	sum(sigdata$Locomotion=="Terrestrial")

centroid1<-matrix(c(ar.ld1,ar.ld2,sc.ld1,sc.ld2,
										t.ld1,t.ld2),nrow=8,ncol=2,byrow=TRUE)


cols=c(rep("firebrick1",30),          ## Arboreal
			 rep("green3",9),               ## Scansorial
			 rep("darkorange",30))          ## Terrestrial

cols.centroid=c(
	rep("firebrick1",1),          ## Arboreal
	rep("green3",1),               ## Scansorial
	rep("darkorange",1))          ## Terrestrial

sybs=c(rep(1,30),          ## Arboreal
			 rep(5,9),               ## Scansorial
			 rep(0,30))          ## Terrestrial

##DF1vsDF2
par(oma=c(2,2,0,0),mar=c(2,2,2,2),mgp=c(3,0.5,0),tck=-0.005)
layout(matrix(c(1,1,1,0,
								2,2,2,3,
								2,2,2,3,
								2,2,2,3),
							nrow=4,ncol=4),
			 width=c(0.5,1,1,1),
			 heights=c(1,1,1,0.5))


## Structual Matrix of DFA 
StruCol<-lda.structure(siglda.res$x,Strudata) 
ld1<-as.matrix(StruCol$LD1)
ld2<-as.matrix(StruCol$LD2)

## Plot all data
plot(siglda.res$x[,1],
		 siglda.res$x[,2],
		 ann=F,
		 pch=sybs,
		 col=cols,
		 las=1,
		 frame=F,cex=1.5)
points(centroid1, 
			 pch=8,
			 col=cols.centroid,
			 cex=2)
mtext(outer=F, side=1, line=2, text="Discriminant Function 1 (55.93%)",cex=0.7)
mtext(outer=F, side=2, line=2, text="Discriminant Function 2 (44.07%)",cex=0.7)

textxy(siglda.res$x[,1], siglda.res$x[,2], labs = rowname, cex = 1)

## Plot LD1
# barplot(t(ld1),ylim=c(-1,1),col="grey80",border="black",
# 				xaxt="none",space=0.9)
# abline(h=0)
# abline(h=-0.5,lty=2,col="grey30")
# abline(h=0.5,lty=2,col="grey30")

###### New plots of data
###### 
StruCol<-lda.structure(siglda.res$x,Strudata) 
head<-c("Sl:Hl","Hsw:Hl","Hpw:Hl","Hdw:Hl","Hsw:Hpw","HHl:Hl",
				"HHw:Hpw","Hdcw:Hpw","Hdcw:Hsw","Hdcw:Hdw","Ul:Hl","Uol:Ul",
				"Uol:Hl","Rl:Hl","Rl:Ul","Uol:Rl","Mcl:(Hl.Rl)","Ppw:Ppl",
				"Ipw:Ipl","(Ppl+Ipl):Mcl","(Mcl+Ppl+Ipl+Dpl):(Hl+Rl)",
				"Il:Pel","FGh:Fl","Fsw:Fl","Tl:Fl","(Hl+Rl):(Tl+Fl)",
				"Tmw:Tl","Cal:Cl","Ctl:Cl","Cal:Ctl")
S<-as.matrix(cbind(ld1,ld2,ld3))

par(mfrow=c(1,2),oma=c(3,3,1,1),mar=c(0,1,0,1),mgp=c(3,0.5,0),tck=-.005)
## Plot CV1 vs CV2
plot(siglda.res$x[,1],
		 siglda.res$x[,2],
		 ann=F,
		 pch=sybs,
		 col=cols,
		 las=1,
		 frame=F,
		 cex=1.5,
		 cex.axis=.8)
points(centroid1, 
			 pch=8,
			 col=cols.centroid,
			 cex=2)
abline(h=0,lty=2,col="grey40")
abline(v=0,lty=2,col="grey40")
text(siglda.res$x[,1],siglda.res$x[,2],mydata$Taxon,cex=.5)
mtext(outer=F, side=1, line=2, text="Canonical function 1 (55.93%)",cex=0.7)
mtext(outer=F, side=2, line=2, text="Canonical function 2 (44.07%)",cex=0.7)

## Coefficient CV1 vs CV2
plot(S, pch=NA_integer_,las=1,frame=F,cex.axis=.8, ann=F,cex.axis=0.6)
for (i in 1:30){
	arrows(0,0,as.numeric(S[i,1]),as.numeric(S[i,2]), 
				 lty=1, 
				 lwd=1, 
				 col="chocolate1",
				 length=0.1)
}
abline(h=0,lty=2,col="grey40")
abline(v=0,lty=2,col="grey40")
text(as.numeric(S[,1]), as.numeric(S[,2]),head,cex=.5)
mtext(outer=F, side=1, line=2, text="Canonical function 1",cex=0.7)
mtext(outer=F, side=2, line=2, text="Canonical function 2",cex=0.7)

