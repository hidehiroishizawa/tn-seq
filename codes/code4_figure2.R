library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggrepel)
library(reshape2)
source("GeomSplitViolin.R")

d<-read.csv("data4_fig2source.csv")

g1<-ggplot(d,aes(x=as.factor(loc),y=FC,fill=time,colour=time,width=2.5))+
  coord_fixed(ratio=0.35)+
  theme_bw(base_size=15)+
  theme(legend.position="none",
        axis.text.x=element_text(angle=30,hjust=1,colour="black",size=11),
        axis.title.y=element_text(size=16),
        axis.text=element_text(size=11,colour="black"))+
  geom_split_violin(trim=F,scale="count",linetype="blank",alpha=I(2/3),adjust=1.3)+
  scale_y_continuous(breaks=seq(-6,10,2),limits=c(-6.2,9.5))+
  scale_x_discrete(labels=c("Amino acid transport\nand metabolism (E)","Transcription (K)","Signal transduction\nmechanisms (T)",
                            "Energy production\nand conversion (C)","Cell wall/membrane/envelope\nbiogenesis (M)",
                            "Carbohydrate transport\nand metabolism (G)","Inorganic ion transport and\nmetabolism (P)",
                            "Cell motility (N)","Lipid transport and\nmetabolism (I)","Post-translational modification,\nprotein turnover, chaperones (O)",
                            "Replication, recombination\nand repair (L)","Minor categories (BDFHJQUV)",
                            "General function prediction only (R)","Function unknown (S)","No COG annotation"))+
  geom_hline(yintercept=c(-1,1),size=0.5,color="grey50",linetype=2)+
  geom_point(data=d %>% filter(sig !="",time=="3h"),aes(x=loc-0.12,y=FC),size=1.5,shape=21,colour="grey20")+
  geom_point(data=d %>% filter(sig !="",time=="7d"),aes(x=loc+0.12,y=FC),size=1.5,shape=21,colour="grey20")+
  geom_text_repel(aes(x=loc-0.1,y=FC,label=name,fontface=ifelse(res=="","italic","bold.italic"),colour=res),
                  data=d %>% filter(sig !="",time=="3h"),
                  direction="y",nudge_x=-0.06,hjust=1,size=2.7,box.padding=0,point.padding=0,
                  segment.size=0,min.segment.length = 10,force=0.8)+
  geom_text_repel(aes(x=loc+0.1,y=FC,label=name,fontface=ifelse(res=="","italic","bold.italic"),colour=res),
                  data=d %>% filter(sig !="",time=="7d"),
                  direction="y",nudge_x=+0.07,hjust=0,size=2.7,box.padding=0,point.padding=0,
                  segment.size=0,min.segment.length = 10,force=0.8)+
  scale_colour_manual(values=c("grey10","#990099","#CC9911","white","white","grey10"))+
  labs(x="",y=expression(paste({Log[2]}," Fold Change (Plant / Control)"))) +
  geom_point(aes(x=14.8,y=-5.3,colour="#f9a39d"),shape=15,size=7.2,colour="#f9a39d")+
  geom_point(aes(x=14.8,y=-6.2,colour="#55d4d8"),shape=15,size=7.2,colour="#55d4d8")+
  geom_text(aes(x=15.3,y=-5.3,label="3 h"),colour="black",size=6.8)+
  geom_text(aes(x=15.3,y=-6.2,label="7 d"),colour="black",size=6.8)+
  annotate("text",x=seq(1,15,1)-0.18,y=8,label=data.frame(subset(d,time=="3h") %>% group_by(loc) %>% tally())$n,colour="black",size=2.9,vjust=-7.15)+
  annotate("text",x=seq(1,15,1)+0.18,y=8,label=data.frame(subset(d,time=="7d") %>% group_by(loc) %>% tally())$n,colour="black",size=2.9,vjust=-7.15)+
  annotate("text",x=14.9,y=-3.2,label="Depleted genes ←",colour="black",size=4.2,angle=90,vjust=6.5)+
  annotate("text",x=14.9,y=3.2,label="→ Enriched genes",colour="Black",size=4.2,angle=90,vjust=6.5)+
  coord_cartesian(xlim=c(0.9,15.4))

win.graph(70,40)
g1
ggsave(file="fig2.tiff",plot=g1,dpi=300)
dev.off()
