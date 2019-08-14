###04/17/2013

##DFA on osteological indices(<0.001)

setwd("/Users/mengchen/Documents/Research Projects/1 Doctoral Dissertation/Chapter 2/Statistical Analysis/LinearDiscriminantAnalysis/Sig0001")


##Libraries and packages
##Sourcing a package to color dendrograms
#source("http://bioinformatics.mdanderson.org/OOMPA/oompaLite.R")
#need to install all packages indenpendently


library(vegan)
library(scatterplot3d)
library(rgl)
library(permute)
library(MASS)
library(calibrate)
library(grid)
library(vcd)
library(colorspace)
library(ggplot2)
library(car)
library("ggthemes")
library(shape)
library(plot3D)
library(plot3Drgl)
source("biostats.R")##helps with data screening, compiled by K. McGarigal at UMass

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

siglda.res.all <- cbind(mydata.new[,1:2], siglda.res$x[,1:2])
colnames(siglda.res.all) <- c("Species No", "Taxon", "DF1", "DF2")
write.csv(siglda.res.all, file = "DF1_2_Results_15Jan2017.csv")


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


## Structual Matrix of DFA 
StruCol<-lda.structure(siglda.res$x,Strudata) 
ld1<-as.matrix(StruCol$LD1)
ld2<-as.matrix(StruCol$LD2)
ld3<-as.matrix(StruCol$LD3)

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

## text the x and y axes
abline(v=0,lty=3)
abline(h=0,lty=3)
mtext(outer=F, side=1, line=2, text="Canonical function 1 (49.13%)",cex=0.7)
mtext(outer=F, side=2, line=2, text="Canonical function 2 (25.83%)",cex=0.7)
textxy(siglda.res$x[,1], siglda.res$x[,2], labs = rowname, cex = 1)

# ## Plot LD1
# barplot(t(ld1),ylim=c(-1,1),col="grey80",border="black",
# 				xaxt="none",space=0.9)
# abline(h=0)
# abline(h=-0.5,lty=2,col="grey30")
# abline(h=0.5,lty=2,col="grey30")

##DF1 vs DF3
# par(oma=c(2,2,0,0),mar=c(2,2,2,2),mgp=c(3,0.5,0),tck=-0.005)
# layout(matrix(c(1,1,1,0,
# 								2,2,2,3,
# 								2,2,2,3,
# 								2,2,2,3),
# 							nrow=4,ncol=4),
# 			 width=c(0.5,1,1,1),
# 			 heights=c(1,1,1,0.5))

## Plot LD3
# barplot(t(ld3),xlim=c(-1,1),col="grey80",border="black",
# 				horiz=T,yaxt="none", space=0.9)
# abline(v=0)
# abline(v=-0.5,lty=2,col="grey30")
# abline(v=0.5,lty=2,col="grey30")
# mtext(outer=F, side=1, line=7,text="Structual \ncorrelation",cex=0.7)

## Plot all data
plot(siglda.res$x[,1],
		 siglda.res$x[,3],
		 ann=F,
		 pch=sybs,
		 col=cols,
		 las=1,
		 frame=F,cex=1.5)
points(centroid2, 
			 pch=8,
			 col=cols.centroid,
			 cex=2)
mtext(outer=F, side=1, line=2, text="Discriminant Function 1 (49.13%)",cex=0.7)
mtext(outer=F, side=2, line=2, text="Discriminant Function 3 (10.54%)",cex=0.7)

# ## Plot LD1
# barplot(t(ld1),ylim=c(-1,1),col="grey80",border="black",
# 				xaxt="none",space=0.9)
# abline(h=0)
# abline(h=-0.5,lty=2,col="grey30")
# abline(h=0.5,lty=2,col="grey30")


#### Plot for the coefficients
head<-c("Sl:Hl","Hsw:Hl","Hpw:Hl","Hdw:Hl","Hsw:Hpw","HHl:Hl",
				"HHw:Hpw","Hdcw:Hpw","Hdcw:Hsw","Hdcw:Hdw","Ul:Hl","Uol:Ul",
				"Uol:Hl","Rl:Hl","Rl:Ul","Uol:Rl","Mcl:(Hl.Rl)","Ppw:Ppl",
				"Ipw:Ipl","(Ppl+Ipl):Mcl","(Mcl+Ppl+Ipl+Dpl):(Hl+Rl)",
				"Il:Pel","FGh:Fl","Fsw:Fl","Tl:Fl","(Hl+Rl):(Tl+Fl)",
				"Tmw:Tl","Cal:Cl","Ctl:Cl","Cal:Ctl")

StruCol<-lda.structure(siglda.res$x,Strudata) 
S<-as.matrix(cbind(ld1,ld2,ld3))

##Structure matrix plot CF1 and CF2
plot(S[,1:2], pch=NA_integer_,las=1,frame=F,cex.axis=.8, ann=F,cex.axis=0.6)
for (i in 1:30){
	arrows(0,0,as.numeric(S[i,1]),as.numeric(S[i,2]), 
				 lty=1, 
				 lwd=1, 
				 col="chocolate1",
				 length=0.1)
}
abline(h=0,lty=2,col="grey40")
abline(v=0,lty=2,col="grey40")
textxy(as.numeric(S[,1]), as.numeric(S[,2]),head,cex=.5)
mtext(outer=F, side=1, line=2, text="Coefficient",las=1,cex=.7)

#######
#######
par(mfrow=c(1,2),oma=c(3,3,1,1),mar=c(0,1,0,1),
		mgp=c(3,0.5,0),tck=-0.005,cex=0.8,cex.axis=0.7)
## Plot CF1 vs CF3
plot(siglda.res$x[,1],
		 siglda.res$x[,3],
		 ann=F,
		 pch=sybs,
		 col=cols,
		 las=1,
		 frame=F)
points(centroid2, 
			 pch=8,
			 col=cols.centroid,
			 cex=2)
mtext(outer=F, side=1, line=2, text="Canonical Function 1 (49.13%)",cex=0.7)
mtext(outer=F, side=2, line=2, text="Canonical Function 3 (10.54%)",cex=0.7)
##Structure matrix plot CF1 and CF3
plot(S[,c(1,3)], pch=NA_integer_,las=1,frame=F,cex.axis=.8, ann=F,cex.axis=0.6)
for (i in 1:30){
	arrows(0,0,as.numeric(S[i,1]),as.numeric(S[i,3]), 
				 lty=1, 
				 lwd=1, 
				 col="chocolate1",
				 length=0.1)
}
abline(h=0,lty=2,col="grey40")
abline(v=0,lty=2,col="grey40")
textxy(as.numeric(S[,1]), as.numeric(S[,3]),head,cex=.5)
mtext(outer=F, side=1, line=2, text="Canonical function 1",las=1,cex=.7)
mtext(outer=F, side=2, line=2, text="Canonical function 3",cex=0.7)


# ## Old codes
# plot(siglda.res$x[,1],
# 		 siglda.res$x[,3],
# 		 main="Eight Locomotor Mode (p<0.001)",
# 		 xlab="DF1(45.13%)",
# 		 ylab="DF3(10.54%)",
# 		 pch=19,col=cols)
# 
# points(centroid2, pch=2,col=cols.centroid)
# 
# textxy(siglda.res$x[,1],
# 			 siglda.res$x[,3],
# 			 rowname,
# 			 cx = 0.5, 
# 			 dcol = "red")

s3d<-scatterplot3d(siglda.res$x[,1],
									 siglda.res$x[,2],
									 siglda.res$x[,3],
									 xlab="DF1(32.52%)",
									 ylab="DF2(35.33%)",
									 zlab="DF3(12.66%)",
									 pch=sybs,color=cols)
s3d.coords<-s3d$xyz.convert(siglda.res$x[,1],
														siglda.res$x[,2],
														siglda.res$x[,3])
text(s3d.coords$x, s3d.coords$y,
		 labels=rowname,cex=0.5,pos=4, col = "red")

par()
scatter3D(siglda.res$x[,1],
       siglda.res$x[,2],
       siglda.res$x[,3],
       xlab="LD1",
       ylab="LD2",
       zlab="LD3",
			 pch=sybs,
       col=cols)

StruCol<-lda.structure(siglda.res$x,Strudata) 

table(sigdata$Locomotion, siglda.res$class)

