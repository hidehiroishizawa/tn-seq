#Location of insertion sites (Fig. S1a)
library(plotly)
library(tidyverse)

adjust<-2.3
d<-data.frame(read_tsv("data2_site_table.txt"))
d$Library<-log(d$Library,10)

x <- c()
y <- c()
xend <- c()
yend <- c()

for(ctr in 1:nrow(d)){
  a <- d$Library[ctr] + adjust
  theta <- d$nc[ctr]*2*pi/4852135
  x[ctr] <- adjust * cos(theta)
  y[ctr] <- adjust * sin(theta)
  xend[ctr] <- a * cos(theta)
  yend[ctr] <- a * sin(theta)
}

plot.df <- data.frame(x, y, xend, yend)
p1 <- plot_ly(plot.df, 
              x = ~x, y = ~y,
              xend = ~xend, yend = ~yend) %>% 
  add_segments(line = list(width = 0.3))

p1<-layout(p1, xaxis = list(domain = c(0, 0.5), title = "", showgrid = F, zeroline = F, showticklabels = F),
           yaxis = list(title = "", showgrid = F, zeroline = F, showticklabels = F))
p1

#Histgram showing the number of insertion sites per gene (Fig. S1b)
library(tidyverse)
d<-data.frame(read_tsv("data3_gene_table.txt"))
win.graph(15,10)
hist(d$insertion_num,breaks=seq(0,100,1),col="#5d87c1",xlab="Number of insertion sites per gene",ylab="Number of genes",cex.lab=1.5,main="")


#Comparison of detection frequencies with Library sample (Fig. S2)
library(tidyverse)
for(i in 9:14){
d<-data.frame(read_tsv("data3_gene_table.txt"))
d<-cbind(d$Library,d[,i])
d<-log10(d)
par(mfrow=c(1,1),bg="white",mar=c(5,5,1,1))
win.graph(10,10)
if(i<=11){j<-4}else{j<-6}
plot(0,0,type="n",xlim=c(-3,j),ylim=c(-3,j),cex.lab=2.4,cex.axis=1.5,xlab="",ylab="")
points(d[,1],d[,2],col="grey20",pch=20,cex=1)
}

#Volcano plot of the 3 h-experiment (Fig. S3)
library(tidyverse)
d<-data.frame(read_tsv("data3_gene_table.txt"))
d<-data.frame("fc"=d$log2FC_3h,"q"=-log(d$qval_3h,10))

win.graph(10,10)
plot(d$fc,d$q,xlim=c(-6,6),xlab="Log2 Fold Change (Plant / Control)",ylab="－Log10 q-value",type="n")
abline(h=-log(0.05,10),lty=2) 
abline(v=1,lty=2)
abline(v=-1,lty=2)

ben<-subset(subset(d,q>log(0.05,10)*-1),fc<1*-1)
det<-subset(subset(d,q>log(0.05,10)*-1),fc>1)
points(d$fc,d$q,col="black",lwd=1,cex=1.3)
points(ben$fc,ben$q,col="blue",lwd=1.5,cex=1.3)
points(det$fc,det$q,col="red",lwd=1.5,cex=1.3)

library(tidyverse)
d<-data.frame(read_tsv("data3_gene_table.txt"))
d<-data.frame("fc"=d$log2FC_7d,"q"=-log(d$qval_7d,10))
plot(d$fc,d$q,xlim=c(-10.5,10.5),xlab="Log2 Fold Change (Plant / Control)",ylab="－Log10 q-value",type="n")
abline(h=-log(0.05,10),lty=2) 
abline(v=1,lty=2)
abline(v=-1,lty=2)

ben<-subset(subset(d,q>log(0.05,10)*-1),fc<1*-1)
det<-subset(subset(d,q>log(0.05,10)*-1),fc>1)
points(d$fc,d$q,col="black",lwd=1,cex=1.3)
points(ben$fc,ben$q,col="blue",lwd=1.5,cex=1.3)
points(det$fc,det$q,col="red",lwd=1.5,cex=1.3)

#Comparison of 3 h and 7 d experiments per COG category (Fig. S4)
library(tidyverse)
library(ggplot2)
library(ggrepel)

gene_table<-read_tsv("data3_gene_table.txt")
gene_table$COG_category<-substr(gene_table$COG_category,1,1)
d<-subset(gene_table,insertion_num>=2)[,c(17,21,20,24,8,7)]
d$res<-paste(d$res_3h,d$res_7d,sep="")

colors<-c("grey65","#0000FF","#0000FF","#FF0000","#990099","#CC9911","#FF0000","#0000FF","#FF0000")
names(colors)<-unique(d$res)
pchs<-c(1,7,0,0,7,7,4,4,7)
names(pchs)<-unique(d$res)

d$COG_category[d$COG_category=="E"]<-"Amino acid transport\nand metabolism (E)"
d$COG_category[d$COG_category=="K"]<-"Transcription (K)"
d$COG_category[d$COG_category=="T"]<-"Signal transduction\nmechanisms (T)"
d$COG_category[d$COG_category=="C"]<-"Energy production\nand conversion (C)"
d$COG_category[d$COG_category=="M"]<-"Cell wall/membrane/envelope\nbiogenesis (M)"
d$COG_category[d$COG_category=="N"]<-"Cell motility (N)"
d$COG_category[d$COG_category=="I"]<-"Lipid transport and\nmetabolism (I)"
d$COG_category[d$COG_category=="G"]<-"Carbohydrate transport\nand metabolism (G)"
d$COG_category[d$COG_category=="P"]<-"Inorganic ion transport and\nmetabolism (P)"
d$COG_category[d$COG_category=="O"]<-"Post-translational modification,\nprotein turnover, chaperones (O)"
d$COG_category[d$COG_category=="L"]<-"Replication, recombination\nand repair (L)"
d$COG_category[d$COG_category=="B"]<-"Minor categories (BDFHJQUV)"
d$COG_category[d$COG_category=="D"]<-"Minor categories (BDFHJQUV)"
d$COG_category[d$COG_category=="F"]<-"Minor categories (BDFHJQUV)"
d$COG_category[d$COG_category=="H"]<-"Minor categories (BDFHJQUV)"
d$COG_category[d$COG_category=="J"]<-"Minor categories (BDFHJQUV)"
d$COG_category[d$COG_category=="Q"]<-"Minor categories (BDFHJQUV)"
d$COG_category[d$COG_category=="U"]<-"Minor categories (BDFHJQUV)"
d$COG_category[d$COG_category=="V"]<-"Minor categories (BDFHJQUV)"
d$COG_category[d$COG_category=="R"]<-"General function prediction only (R)"
d$COG_category[d$COG_category=="S"]<-"Function unknown (S)"
d$COG_category[d$COG_category==NA]<-"No COG annotation"

tiff("figS4.tiff", units="in", width=20, height=16, res=300)
  g<-ggplot(NULL)+
    geom_point(data=d,aes(x=log2FC_3h,y=log2FC_7d,color=res,shape=res),size=2.2,stroke=1.1)+
    scale_color_manual(values=colors)+
    scale_shape_manual(values=pchs)+
    scale_x_continuous(limits=c(-5,10))+
    scale_y_continuous(limits=c(-5,10))+
    labs(x=expression(paste("3 h ",{Log[2]}," FC (Plant/Control)")),y=expression(paste("7 d ",{Log[2]}," FC (Plant/Control)")))+
    coord_fixed(ratio=1)+
    theme_bw(base_size=20)+
    geom_hline(yintercept=0,size=0.5,color="grey50",linetype=2)+
    geom_vline(xintercept=0,size=0.5,color="grey50",linetype=2)+
    geom_text_repel(data=subset(d,res!="NANA"),aes(x=log2FC_3h,y=log2FC_7d,label=gene_name),size=3)+
    theme(legend.position="none")+
    facet_wrap(~COG_category,ncol=5)+
    theme(
      strip.text=element_text(size=13,face="bold",margin=margin(1.5,0,1.5,0,"mm"))
    )
  g
    dev.off()
    