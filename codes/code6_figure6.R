library(ggplot2)
library(agricolae)
library(tidyverse)

d<-read.csv("data6_fig6source.csv")
cor.test(subset(d,time=="3 h")$tn_log2fc,subset(d,time=="3 h")$perc_mut,method="spearman",exact=F)
cor.test(subset(d,time=="7 d")$tn_log2fc,subset(d,time=="7 d")$perc_mut,method="spearman",exact=F)

shapes<-c(21,22,23,24,25,23,22)
names(shapes)<-unique(d$strain)

d1<-subset(d,rep==1)
d2<-subset(d,rep==2)
d3<-subset(d,rep==3)

g<-ggplot(NULL)+
  geom_point(data=d1,aes(x=tn_log2fc,y=perc_mut,fill=time,shape=strain),size=6,stroke=1.3)+
  geom_point(data=subset(d1,strain=="ΔwecC"),aes(x=tn_log2fc,y=perc_mut),shape=3,size=5,stroke=1.3)+
  geom_point(data=subset(d1,strain=="Δrtn"),aes(x=tn_log2fc,y=perc_mut),shape=4,size=5,stroke=1.3)+
  geom_point(data=d2,aes(x=tn_log2fc,y=perc_mut,fill=time,shape=strain),size=6,stroke=1.3)+
  geom_point(data=subset(d2,strain=="ΔwecC"),aes(x=tn_log2fc,y=perc_mut),shape=3,size=5,stroke=1.3)+
  geom_point(data=subset(d2,strain=="Δrtn"),aes(x=tn_log2fc,y=perc_mut),shape=4,size=5,stroke=1.3)+
  geom_point(data=d3,aes(x=tn_log2fc,y=perc_mut,fill=time,shape=strain),size=6,stroke=1.3)+
  geom_point(data=subset(d3,strain=="ΔwecC"),aes(x=tn_log2fc,y=perc_mut),shape=3,size=5,stroke=1.3)+
  geom_point(data=subset(d3,strain=="Δrtn"),aes(x=tn_log2fc,y=perc_mut),shape=4,size=5,stroke=1.3)+
  scale_y_continuous(limits=c(0,100))+
  scale_x_continuous(limits=c(-6,6))+
  scale_shape_manual(values=shapes,guide=F)+
  scale_fill_manual(values=c("#f9a39d","#55d4d8"),guide=F)+
  labs(x=expression(paste({Log[2]},"FC (Plant/Control) in Tn-seq")),y="Colonization ratio of mutants (%)")+
  theme_minimal()+
  theme(
    aspect.ratio=1,
    axis.line=element_line(size=1.2),
    axis.text=element_text(size=20,color="black"),
    axis.ticks=element_line(size=1.2),
    axis.ticks.length=unit(3,"mm"),
    axis.title=element_text(size=26)
  )+
  geom_hline(yintercept=50,size=1,color="grey50",linetype=2)+
  geom_vline(xintercept=0,size=1,color="grey50",linetype=2)+
  coord_cartesian(ylim=c(2,98),xlim=c(-5.7,5.7))

dev.off()
win.graph(10,8)
g
