library(ggplot2)
library(agricolae)
library(tidyverse)

d<-read.csv("data5_fig5csource.csv")
d$system<-paste(d$strain,d$condition,sep="_")
d3<-subset(d,time=="3h");d7<-subset(d,time=="7d")

hsd_3h<-HSD.test(aov(perc_mut~system,data=d3),"system",group=T)$groups
hsd_7d<-HSD.test(aov(perc_mut~system,data=d7),"system",group=T)$groups
hsd_3h$system<-rownames(hsd_3h);hsd_7d$system<-rownames(hsd_7d)
d3<-left_join(d3,hsd_3h,by="system");d7<-left_join(d7,hsd_7d,by="system")


g1<-ggplot(NULL)+
  geom_jitter(data=d3,aes(x=reorder(system,as.numeric(rownames(d3))),y=perc_mut.x),size=6,shape=21,width=0.02,stroke=1.3,fill="#f9a39d")+
  scale_y_continuous(limits=c(0,100))+
  scale_x_discrete(labels=subset(d3,rep==1)$strain)+
  labs(y="",x="")+
  geom_text(data=d3,aes(x=system,y=perc_mut.y+15,label=groups),size=7)+
  geom_text(aes(x=1,y=95,label="3 h"),size=14,fontface="bold")+
  theme_minimal(base_size=18)+
  theme(
    axis.line=element_line(size=1.2),
    axis.text.x=element_text(size=25,color="black",angle=90,hjust=1,vjust=0.5),
    axis.text.y=element_text(size=25,color="black"),
    axis.ticks=element_line(size=1.2),
    axis.ticks.length=unit(2.7,"mm")
  )+
  coord_cartesian(ylim=c(4,96),xlim=c(1,10))
win.graph(10,7)
g1
g2<-ggplot(NULL)+
  geom_jitter(data=d7,aes(x=reorder(system,as.numeric(rownames(d7))),y=perc_mut.x),size=6,shape=21,width=0.02,stroke=1.3,fill="#55d4d8")+
  scale_y_continuous(limits=c(0,100))+
  scale_x_discrete(labels=subset(d7,rep==1)$strain)+
  labs(y="",x="")+
  geom_text(data=d7,aes(x=system,y=perc_mut.y+15,label=groups),size=7)+
  geom_text(aes(x=1,y=95,label="7 d"),size=14,fontface="bold")+
  theme_minimal(base_size=18)+
  theme(
    axis.line=element_line(size=1.2),
    axis.text.x=element_text(size=25,color="black",angle=90,hjust=1,vjust=0.5),
    axis.text.y=element_text(size=25,color="black"),
    axis.ticks=element_line(size=1.2),
    axis.ticks.length=unit(2.7,"mm")
  )+
  coord_cartesian(ylim=c(4,96),xlim=c(1,10))
win.graph(10,7)
g2
