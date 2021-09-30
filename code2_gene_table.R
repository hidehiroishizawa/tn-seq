setwd("C:/Users/uadgw/デスクトップ/WorkingFiles/Tnseq原稿/CommunBiol投稿/Dryad")
library(reshape2)
library(tidyverse)
library(qvalue)
site_table<-data.frame(read_tsv("data2_site_table.txt"))
cds<-data.frame(read_tsv("data0_CDS_info.txt"))

site_table$sum3h<-site_table$X3hControl+site_table$X3hPlant
site_table$sum7d<-site_table$X7dControl+site_table$X7dPlant

#Calculate log2 Plant/Control using insertion site data >=30 read counts
genes<-na.omit(unique(site_table$Locus_tag))
results<-data.frame("Locus_tag"=genes,
                    "log2FC_3h"=rep(NA,length(genes)),"pval_3h"=rep(NA,length(genes)),"qval_3h"=rep(NA,length(genes)),"res_3h"=rep(NA,length(genes)),
                    "log2FC_7d"=rep(NA,length(genes)),"pval_7d"=rep(NA,length(genes)),"qval_7d"=rep(NA,length(genes)),"res_7d"=rep(NA,length(genes)))
freq3h<-subset(site_table,sum3h>=30)[,c(9,2,4)]
freq7d<-subset(site_table,sum7d>=30)[,c(9,5,7)]
freq3h[,2]<-(freq3h[,2]+1)/median(freq3h[,2]+2);freq3h[,3]<-(freq3h[,3]+1)/median(freq3h[,3]+2)
freq7d[,2]<-(freq7d[,2]+1)/median(freq7d[,2]+2);freq7d[,3]<-(freq7d[,3]+1)/median(freq7d[,3]+2)

for (i in genes){
  mutants3h<-subset(freq3h,Locus_tag==i)
  mutants7d<-subset(freq7d,Locus_tag==i)
  results$log2FC_3h[match(i,genes)]<-log(sum(mutants3h[,3])/sum(mutants3h[,2]),2)
  results$log2FC_7d[match(i,genes)]<-log(sum(mutants7d[,3])/sum(mutants7d[,2]),2)
  mutants3h$fc<-log(mutants3h[,2]/mutants3h[,3],2)
  mutants7d$fc<-log(mutants7d[,2]/mutants7d[,3],2)
  if(nrow(mutants3h)>=2 && var(mutants3h$fc)>1E-10){
    results$pval_3h[match(i,genes)]<-t.test(mutants3h$fc,mu=0)$p.value
  }
  if(nrow(mutants7d)>=2 && var(mutants7d$fc)>1E-10){
    results$pval_7d[match(i,genes)]<-t.test(mutants7d$fc,mu=0)$p.value
  }
}
results$qval_3h<-qvalue(results$pval_3h)$qvalues
results$qval_7d<-qvalue(results$pval_7d)$qvalues

#test for normal distribution
ks.test(x=na.omit(results$log2FC_3h),y="pnorm",mean=mean(na.omit(results$log2FC_3h)),sd=sd(na.omit(results$log2FC_3h)))
ks.test(x=na.omit(results$log2FC_7d),y="pnorm",mean=mean(na.omit(results$log2FC_7d)),sd=sd(na.omit(results$log2FC_7d)))


#Determine Enriched and Depleted genes (FC>2, q<0.05)
for(i in 1:nrow(results)){
  if(is.na(results$qval_3h[i])==FALSE){
    if (results$qval_3h[i]<0.05){
      if(results$log2FC_3h[i] > 1){results$res_3h[i]<-"Enriched"}
      if(results$log2FC_3h[i] < -1){results$res_3h[i]<-"Depleted"}
    }}
  if(is.na(results$qval_7d[i])==FALSE){
    if (results$qval_7d[i]<0.05){
      if(results$log2FC_7d[i] > 1){results$res_7d[i]<-"Enriched"}
      if(results$log2FC_7d[i] < -1){results$res_7d[i]<-"Depleted"}
    }}
}

#Median normalization,
for(j in 2:8){
  site_table[,j]<-(site_table[,j]+1)/median(site_table[,j]+1)
}

#Re-summarize site table into gene table
d<-melt(site_table,id.vars=c("nc","Locus_tag"),
        measure.vars=c("X3hControl","X3hMedium","X3hPlant","X7dControl","X7dMedium","X7dPlant","Library"))
insertion_num<-dcast(d,Locus_tag~variable,length)[,2]
gene_table<-cbind(dcast(d,Locus_tag~variable,sum)[,],insertion_num)
gene_table<-cbind(gene_table[-nrow(gene_table),],results[,-1])
gene_table<-left_join(cds,gene_table)

write.table(gene_table,"data3_gene_table.txt",sep="\t",row.names=F)

