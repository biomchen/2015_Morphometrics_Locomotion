# created on 04/17/2013
# reorganized on 24Oct2019

# Analyses on morphometrics data

setwd("you directory")

# Libraries and packages
library(vegan)
library(scatterplot3d)
library(rgl)
library(permute)
library(MASS)
library(calibrate)
library(car)
library(shape)
source("biostats.R")##helps with data screening, compiled by K. McGarigal at UMass

# reading the data
mydata.new <- read.table("morphmetrics data", head=T, sep=",")
sigdata<-mydata.new[,4:34]

######################################################
# testing for significant differences between groups #
######################################################

sigindex <- as.matrix(sigdata[,c(2:31)])
sigoutput <- manova(sigindex~sigdata$Locomotion)
sigsum.wilks <- summary(sigoutput, test="Wilks")
sigsum.avo <- summary.aov(sigoutput)

######################################
# DFA on Significant Indices(<0.001) #
######################################

rowname <- mydata.new[,1]
Strudata <- sigdata[,2:31]
siglda.mod <- lda(Locomotion ~ .,
                  data=sigdata,
                  prior=c(0.125,0.125,0.125,0.125, 0.125,0.125,0.125,0.125))
siglda.res <- predict(siglda.mod, sigdata)
siglda.res.all <- cbind(mydata.new[,1:2], siglda.res$x[,1:2])
colnames(siglda.res.all) <- c("Species No", "Taxon", "DF1", "DF2")
# First two coefficientS of the DFA
coemtr<-cbind(siglda.mod$scaling[,1],siglda.mod$scaling[,2])

# Calculate the centroid of each locomotion
# Arboreal
ar.ld1=sum(siglda.res$x[,1]*(sigdata$Locomotion=="Arboreal"))/
	sum(sigdata$Locomotion=="Arboreal")
ar.ld2=sum(siglda.res$x[,2]*(sigdata$Locomotion=="Arboreal"))/
	sum(sigdata$Locomotion=="Arboreal")
ar.ld3=sum(siglda.res$x[,3]*(sigdata$Locomotion=="Arboreal"))/
	sum(sigdata$Locomotion=="Arboreal")
# Fossorial
f.ld1=sum(siglda.res$x[,1]*(sigdata$Locomotion=="Fossorial"))/
	sum(sigdata$Locomotion=="Fossorial")
f.ld2=sum(siglda.res$x[,2]*(sigdata$Locomotion=="Fossorial"))/
	sum(sigdata$Locomotion=="Fossorial")
f.ld3=sum(siglda.res$x[,3]*(sigdata$Locomotion=="Fossorial"))/
	sum(sigdata$Locomotion=="Fossorial")
# Gliding
g.ld1=sum(siglda.res$x[,1]*(sigdata$Locomotion=="Gliding"))/
	sum(sigdata$Locomotion=="Gliding")
g.ld2=sum(siglda.res$x[,2]*(sigdata$Locomotion=="Gliding"))/
	sum(sigdata$Locomotion=="Gliding")
g.ld3=sum(siglda.res$x[,3]*(sigdata$Locomotion=="Gliding"))/
	sum(sigdata$Locomotion=="Gliding")
# Saltatorial
s.ld1=sum(siglda.res$x[,1]*(sigdata$Locomotion=="Saltatorial"))/
	sum(sigdata$Locomotion=="Saltatorial")
s.ld2=sum(siglda.res$x[,2]*(sigdata$Locomotion=="Saltatorial"))/
	sum(sigdata$Locomotion=="Saltatorial")
s.ld3=sum(siglda.res$x[,3]*(sigdata$Locomotion=="Saltatorial"))/
	sum(sigdata$Locomotion=="Saltatorial")
# Scansorial
sc.ld1=sum(siglda.res$x[,1]*(sigdata$Locomotion=="Scansorial"))/
	sum(sigdata$Locomotion=="Scansorial")
sc.ld2=sum(siglda.res$x[,2]*(sigdata$Locomotion=="Scansorial"))/
	sum(sigdata$Locomotion=="Scansorial")
sc.ld3=sum(siglda.res$x[,3]*(sigdata$Locomotion=="Scansorial"))/
	sum(sigdata$Locomotion=="Scansorial")
# Semi-aquatic
sa.ld1=sum(siglda.res$x[,1]*(sigdata$Locomotion=="Semi-aquatic"))/
	sum(sigdata$Locomotion=="Semi-aquatic")
sa.ld2=sum(siglda.res$x[,2]*(sigdata$Locomotion=="Semi-aquatic"))/
	sum(sigdata$Locomotion=="Semi-aquatic")
sa.ld3=sum(siglda.res$x[,3]*(sigdata$Locomotion=="Semi-aquatic"))/
	sum(sigdata$Locomotion=="Semi-aquatic")
# Semi-fossorial
sf.ld1=sum(siglda.res$x[,1]*(sigdata$Locomotion=="Semi-fossorial"))/
	sum(sigdata$Locomotion=="Semi-fossorial")
sf.ld2=sum(siglda.res$x[,2]*(sigdata$Locomotion=="Semi-fossorial"))/
	sum(sigdata$Locomotion=="Semi-fossorial")
sf.ld3=sum(siglda.res$x[,3]*(sigdata$Locomotion=="Semi-fossorial"))/
	sum(sigdata$Locomotion=="Semi-fossorial")
# Terrestrial
t.ld1=sum(siglda.res$x[,1]*(sigdata$Locomotion=="Terrestrial"))/
	sum(sigdata$Locomotion=="Terrestrial")
t.ld2=sum(siglda.res$x[,2]*(sigdata$Locomotion=="Terrestrial"))/
	sum(sigdata$Locomotion=="Terrestrial")
t.ld3=sum(siglda.res$x[,3]*(sigdata$Locomotion=="Terrestrial"))/
	sum(sigdata$Locomotion=="Terrestrial")
# centroids
centroid1<-matrix(c(ar.ld1,ar.ld2,f.ld1,f.ld2,g.ld1,g.ld2,
                    s.ld1,s.ld2,sc.ld1,sc.ld2,sa.ld1,sa.ld2,
                    sf.ld1,sf.ld2,t.ld1,t.ld2),
                  nrow=8,ncol=2,byrow=TRUE)

centroid2<-matrix(c(ar.ld1,ar.ld3,f.ld1,f.ld3,g.ld1,g.ld3,
                    s.ld1,s.ld3,sc.ld1,sc.ld3,sa.ld1,sa.ld3,
                    sf.ld1,sf.ld3,t.ld1,t.ld3),
                  nrow=8,ncol=2,byrow=TRUE)

# Plots the results
# colors
cols=c(rep("firebrick1",30), # Arboreal
       rep("grey36",12),     # Fossorial
	   rep("dodgerblue3",3), # Gliding
       rep("magenta",5),     # Saltatorial
       rep("green3",9),      # Scansorial
       rep("deepskyblue",9), # Semi-aquatic
       rep("yellow3",9),     # Semi-fossorial
       rep("darkorange",30)) # Terrestrial
# symbols
sybs=c(rep(1,30),  # Arboreal
       rep(15,12), # Fossorial
       rep(2,3),   # Gliding
       rep(4,5),   # Saltatorial
       rep(5,9),   # Scansorial
       rep(18,9),  # Semiaquatic
       rep(19,9),  # Semifossorial
       rep(0,30))  # Terrestrial
# centroid colors
cols.centroid=c(rep("firebrick1",1),  # Arboreal
                rep("grey36",1),      # Fossorial
                rep("dodgerblue3",1), # Gliding
                rep("magenta",1),     # Saltatorial
                rep("green3",1),      # Scansorial
                rep("deepskyblue",1), # Semi-aquatic
                rep("yellow3",1),     # Semi-fossorial
                rep("darkorange",1))  # Terrestrial

# DF1 vs DF2
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
abline(v=0,lty=3)
abline(h=0,lty=3)
mtext(outer=F, side=1, line=2, text="Canonical function 1 (49.13%)",cex=0.7)
mtext(outer=F, side=2, line=2, text="Canonical function 2 (25.83%)",cex=0.7)
textxy(siglda.res$x[,1], siglda.res$x[,2], labs=rowname, cex=1)

# DF1 vs DF3
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
textxy(siglda.res$x[,1], siglda.res$x[,3], labs=rowname, cex=1)

# DF2 vs DF3
plot(siglda.res$x[,2],
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
textxy(siglda.res$x[,2], siglda.res$x[,3], labs=rowname, cex=1)

# 3D plot
s3d <- scatterplot3d(siglda.res$x[,1],
                     siglda.res$x[,2],
                     siglda.res$x[,3],
                     xlab="DF1(32.52%)",
                     ylab="DF2(35.33%)",
                     zlab="DF3(12.66%)",
                     pch=sybs,color=cols)
s3d.coords <- s3d$xyz.convert(siglda.res$x[,1],
                              siglda.res$x[,2],
                              siglda.res$x[,3])
text(s3d.coords$x, s3d.coords$y, labels=rowname, cex=0.5, pos=4, col="red")

par()
scatter3D(siglda.res$x[,1],
          siglda.res$x[,2],
          siglda.res$x[,3],
          xlab="LD1",
          ylab="LD2",
          zlab="LD3",
          pch=sybs,
          col=cols)

######################################
# structural matrix coefficient plot #
######################################

# header
header <- c("Sl:Hl","Hsw:Hl","Hpw:Hl","Hdw:Hl","Hsw:Hpw","HHl:Hl",
            "HHw:Hpw","Hdcw:Hpw","Hdcw:Hsw","Hdcw:Hdw","Ul:Hl","Uol:Ul",
            "Uol:Hl","Rl:Hl","Rl:Ul","Uol:Rl","Mcl:(Hl.Rl)","Ppw:Ppl",
            "Ipw:Ipl","(Ppl+Ipl):Mcl","(Mcl+Ppl+Ipl+Dpl):(Hl+Rl)",
            "Il:Pel","FGh:Fl","Fsw:Fl","Tl:Fl","(Hl+Rl):(Tl+Fl)",
            "Tmw:Tl","Cal:Cl","Ctl:Cl","Cal:Ctl")

# structual Matrix of DFA
StruCol <- lda.structure(siglda.res$x, Strudata)
ld1 <- as.matrix(StruCol$LD1)
ld2 <- as.matrix(StruCol$LD2)
ld3 <- as.matrix(StruCol$LD3)
S <- as.matrix(cbind(ld1,ld2,ld3))

# set up parameter of layout
par(mfrow=c(1,2), oma=c(3,3,1,1),
    mar=c(0,1,0,1), mgp=c(3,0.5,0),
    tck=-0.005, cex=0.8, cex.axis=0.7)

# DF1 vs DF2
plot(S[,1:2], pch=NA_integer_,las=1, frame=F, cex.axis=.8, ann=F, cex.axis=0.6)
for (i in 1:30){
    arrows(0, 0,
           as.numeric(S[i,1]), as.numeric(S[i,2]),
           lty=1, lwd=1,
           col="chocolate1", length=0.1)
}
abline(h=0, lty=2, col="grey40")
abline(v=0, lty=2, col="grey40")
textxy(as.numeric(S[,1]), as.numeric(S[,2]), header, cex=.5)
mtext(outer=F, side=1, line=2, text="Coefficient", las=1, cex=.7)
# DF1 vs DF3
plot(S[,c(1,3)], pch=NA_integer_, las=1, frame=F, cex.axis=.8, ann=F, cex.axis=0.6)
for (i in 1:30){
	arrows(0, 0,
           as.numeric(S[i,1]), as.numeric(S[i,3]),
           lty=1, lwd=1,
           col="chocolate1", length=0.1)
}
abline(h=0, lty=2, col="grey40")
abline(v=0, lty=2, col="grey40")
textxy(as.numeric(S[,1]), as.numeric(S[,3]), header, cex=.5)
mtext(outer=F, side=1, line=2, text="Canonical function 1", las=1, cex=.7)
mtext(outer=F, side=2, line=2, text="Canonical function 3", cex=0.7)

