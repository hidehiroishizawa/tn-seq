library(tidyverse)

cds<-read_tsv("data0_CDS_info.txt")
data<-read_csv("data1_insertion_sites.csv")

#Remove insertion sites with <10 total read counts
d<-data[rowSums(data[,-1])>=10,]

#Remove insertion sites with unnaturally different Left to Right ratio
d<-cbind(d,rep("",nrow(d)))
for(i in 1:nrow(d)){
  for(j in 1:7){
  if(sum(d[i,2*j],d[i,2*j+1])>0){
    p<-binom.test(min(d[i,2*j],d[i,2*j+1]),sum(d[i,2*j],d[i,2*j+1]),0.1)
    if(p$p.value<0.05 && p$estimate<0.1){
    d[i,16]<-"out"
    break}
  }}
}
d<-d[d[,16]!="out",-16]

#Assemble left- and right-side read counts
d<-data.frame("nc"=d$nc,"3hControl"=d[,2]+d[,3],"3hMedium"=d[,4]+d[,5],"3hPlant"=d[,6]+d[,7],
              "7dControl"=d[,8]+d[,9],"7dMedium"=d[,10]+d[,11],"7dPlant"=d[,12]+d[,13],"Library"=d[,14]+d[,15])

#Combine insertion site data with gene information
site_table<-numeric()
d$nc<-as.numeric(d$nc)
cds<-rbind(cds,rep(0,8))
k<-1
for(i in 1:nrow(d)){
  hit<-0
  for(j in k:(nrow(cds)-1)){
     if(d$nc[i]>cds$start[j] && d$nc[i]<=cds$stop[j]){
       site_table<-rbind(site_table,c(d[i,],cds[j,]))
       hit<-1;k<-j
       if(d$nc[i]>cds$start[j+1] && d$nc[i]<=cds$stop[j+1]){
         site_table<-rbind(site_table,c(d[i,],cds[j+1,]))
       }
     break} 
    }
  if(hit==0){site_table<-rbind(site_table,c(d[i,],rep("",8)))}
}

#Save site table
write.table(site_table,"data2_site_table.txt",sep="\t",row.names=F)
