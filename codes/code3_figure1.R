#3dplot for the 3 h-experiments (Fig1c)
library(rgl)
library(tidyverse)

d<-data.frame(read_tsv("data3_gene_table.txt"))
d[,9:15]<-log10(d[,9:15])
d<-d[is.na(d$insertion_num)==FALSE,]
na.omit(d)

ns3h<-d[is.na(d$res_3h)==T,][,c(1,9,10,11)]
dep3h<-subset(d,res_3h=="Depleted")[,c(1,9,10,11)]
enr3h<-subset(d,res_3h=="Enriched")[,c(1,9,10,11)]

plot3d(ns[,2],ns[,4],ns[,3],size=6,box=F,xlim=c(-1.1,2.5),ylim=c(-1.1,2.5),zlim=c(-1.1,2.5),
       xlab="",ylab="",zlab="")
idx <- grid3d('x')      
idy <- grid3d('y')
idz <- grid3d('z')                 
gridx <- rgl.attrib(idx[1], "vertices")[1,3]
gridy <- rgl.attrib(idy[1], "vertices")[1,3] 
gridz <- rgl.attrib(idz[1], "vertices")[1,3]  
save <- par3d(ignoreExtent = TRUE)          
with(ns, points3d(ns[,2], gridy, ns[,3], col = "gray70",size=5))
with(ns, points3d(gridx, ns[,4], ns[,3], col = "gray70",size=5))
with(ns, points3d(ns[,2],ns[,4], gridz, col = "gray70",size=5))
with(ns, points3d(pos[,2],pos[,4],pos[,3],col="blue",size=5))
with(ns, points3d(neg[,2],neg[,4],neg[,3],col="red",size=5))
with(ns, points3d(pos[,2],gridy,pos[,3],col="#c8c8ff",size=5))
with(ns, points3d(neg[,2],gridy,neg[,3],col="#ffc8c8",size=5))
with(ns, points3d(gridx,pos[,4],pos[,3],col="#c8c8ff",size=5))
with(ns, points3d(gridx,neg[,4],neg[,3],col="#ffc8c8",size=5))
with(ns, points3d(pos[,2],pos[,4],gridz,col="#c8c8ff",size=5))
with(ns, points3d(neg[,2],neg[,4],gridz,col="#ffc8c8",size=5))
par3d(save)        

rgl.viewpoint(30,15,1)
rgl.postscript("3dplot3h.pdf")


#3dplot for the 7 d-experiments (Fig1c)

library(rgl)
library(tidyverse)

d<-data.frame(read_tsv("data3_gene_table.txt"))
d[,9:15]<-log10(d[,9:15])
d<-d[is.na(d$insertion_num)==FALSE,]
na.omit(d)

ns3h<-d[is.na(d$res_3h)==T,][,c(1,12,13,14)]
dep3h<-subset(d,res_3h=="Depleted")[,c(1,12,13,14)]
enr3h<-subset(d,res_3h=="Enriched")[,c(1,12,13,14)]

plot3d(ns[,2],ns[,4],ns[,3],size=6,box=F,xlim=c(-1.3,5.9),ylim=c(-1.3,5.9),zlim=c(-1.3,5.9),
       xlab="",ylab="",zlab="")
idx <- grid3d('x')      
idy <- grid3d('y')
idz <- grid3d('z')                 
gridx <- rgl.attrib(idx[1], "vertices")[1,3]
gridy <- rgl.attrib(idy[1], "vertices")[1,3] 
gridz <- rgl.attrib(idz[1], "vertices")[1,3]  
save <- par3d(ignoreExtent = TRUE)          
with(ns, points3d(ns[,2], gridy, ns[,3], col = "gray70",size=5))
with(ns, points3d(gridx, ns[,4], ns[,3], col = "gray70",size=5))
with(ns, points3d(ns[,2],ns[,4], gridz, col = "gray70",size=5))
with(ns, points3d(pos[,2],pos[,4],pos[,3],col="blue",size=5))
with(ns, points3d(neg[,2],neg[,4],neg[,3],col="red",size=5))
with(ns, points3d(pos[,2],gridy,pos[,3],col="#c8c8ff",size=5))
with(ns, points3d(neg[,2],gridy,neg[,3],col="#ffc8c8",size=5))
with(ns, points3d(gridx,pos[,4],pos[,3],col="#c8c8ff",size=5))
with(ns, points3d(gridx,neg[,4],neg[,3],col="#ffc8c8",size=5))
with(ns, points3d(pos[,2],pos[,4],gridz,col="#c8c8ff",size=5))
with(ns, points3d(neg[,2],neg[,4],gridz,col="#ffc8c8",size=5))
par3d(save)        

rgl.viewpoint(30,15,1)
rgl.postscript("fig1c_7d.pdf")


#Comparison of 3 h vs 7 d experiment (Fig1d)
library(tidyverse)

gene_table<-read_tsv("data3_gene_table.txt")
d<-subset(gene_table,insertion_num>=2)[,c(17,21,20,24)]
d$res<-paste(d$res_3h,d$res_7d,sep="")

win.graph(10,10)
par(lwd=2,pty="s")
plot(d[,1:2],xlim=c(-4.9,9.7),ylim=c(-4.9,9.7),type="n",ylab="",xlab="",cex.axis=1.5)
par(lwd=1)
points(subset(d,res=="NANA")[,1:2],pch=1,col="grey45",cex=1.5)
par(lwd=2.2)
points(subset(d,res=="DepletedDepleted")[,1:2],pch=7,col="#0000FF",cex=1.7,lwd=2)
points(subset(d,res=="NADepleted")[,1:2],pch=0,col="#0000FF",cex=1.6)
points(subset(d,res=="DepletedNA")[,1:2],pch=4,col="#0000FF",cex=1.6)
points(subset(d,res=="EnrichedEnriched")[,1:2],pch=7,col="#FF0000",cex=1.7,lwd=2)
points(subset(d,res=="NAEnriched")[,1:2],pch=0,col="#FF0000",cex=1.6)
points(subset(d,res=="EnrichedNA")[,1:2],pch=4,col="#FF0000",cex=1.6)
points(subset(d,res=="DepletedEnriched")[,1:2],pch=7,col="#990099",cex=1.7,lwd=2)
points(subset(d,res=="EnrichedDepleted")[,1:2],pch=7,col="#CC9911",cex=1.7,lwd=2)

abline(h=0,lty=5,lwd=1.5)
abline(v=0,lty=5,lwd=1.5)
legend("topright",
       legend=c("3 h-Depleted","3 h-Enriched","7 d-Depleted","7 d-Enriched","3 h-Depleted","7 d-Enriched","3 h-Enriched","7 d-Depleted"),
       col=c("#0000FF","#FF0000","#0000FF","#FF0000","#990099","#FFFFFF","#CC9911","#FFFFFF"),
       pch=c(4,4,0,0,7,7,7,7),
       cex=1.5,pt.cex=2.5)

