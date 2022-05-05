#Dependencies and colors
library(vegan)
library(dplyr)
library(reshape2)
library(rary(ggsci)
library(RColorBrewer)
library(scales)
library(data.table)
mypal =pal_d3("category20")(20)
########################################################################
#tea genetic dissimilarity from PLINK output
        
tea.dist<- read.table("D:/2.tea/tea_resequencing/SNP-202012/plink.mdist", quote="\"", comment.char="")
#tea cultivar IDs
        
tea_names<-read.delim("D:/2.tea/tea_resequencing/SNP-202012/plink.mdist.id", header=FALSE)
column<-tea_names$V1
row<-tea_names$V2
dimnames(tea.dist)=lrow,column)
###############################################

#Cluster genetic dist into 2 and 3 groupcust(tea.dist,method = "complete")
clclust<-hclust(tea.dist,method = "complete")
clclust.id = cutree(clclust, k=2)
clclust.id3 =cutree(clclust, k=3)
clclust.id<-data.frame(Cultivar=rownames(tea_pca.dist),group=clclust.id,group3=clclust.id3)

#Admixture Cross Validation errors#
cv<-fread("2.admixture/adm.csv",header = T)


p1<-ggplot(cv, aes(x =K,y = CV_Error))+
  geom_line(linetype = "dashed",color="#BEBEBE")+
  geom_point(size=3,shape=17,color="black")+
  theme_bw()+  xlab("K")+
  ylab("CV Error")+
  theme(legend.position  = "right",panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.background = element_blank(),
        axis.line=element_line(color = "black"),
        axis.text = element_text(colour = "black")
  )+
  scale_x_continuous(limits = c(1,10),breaks = seq(1,10,1))

#Admixture barplots
ta2<-read.table("2.admixture/tea_reseq2021.2.Q")
ta3<-read.table("2.admixture/tea_reseq2021.3.Q")

head(ta2)
ta2.1<-mutate(ta2,Cultivar=clclust.id$Cultivar,cluster2=clclust.id$group)
ta3.1<-mutate(ta3,Cultivar=clclust.id$Cultivar,cluster3=clclust.id$group3)

ta2.2<-melt(ta2.1)
ta2.2<-melt(ta2.1)

names(ta3.2)<-c("Cultivar","Ancestry","value")
names(ta2.2)<-c("Cultivar","Ancestry","value")




ta3.2<-ta3.2[-301:-400,]
ta2.2<-ta2.2[-301:-400,]

#Group3
  bar.plot<-ggplot(ta3.2,aes(x=factor(Cultivar,levels = rank$Cultivar),y=value*100,fill=factor(Ancestry)))+
    geom_bar(stat = "identity",width = 1) +    
    theme_minimal()+
    xlab("Cultivar")+
    ylab("Ancestry")+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          panel.background = element_blank(),panel.border = element_blank(),legend.position = "none")+
    theme(axis.line = element_line(color = "black"),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank())+
    scale_fill_manual(values = c(mypal[1],mypal[5],mypal[2]),name="")+
    scale_y_continuous(expand = c(0,0))+
    coord_fixed(ratio=1/4)+
    geom_vline(xintercept =61.5,color="white")+
    geom_vline(xintercept =72.5,color="white")
  
  
  
  #Group2
  bar.plot1<-ggplot(ta2.2,aes(x=factor(Cultivar,levels = rank$Cultivar),y=value*100,fill=factor(Ancestry)))+
    geom_bar(stat = "identity",width = 1) +    
    theme_minimal()+
    xlab("Cultivar")+
    ylab("Ancestry")+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          panel.background = element_blank(),panel.border = element_blank(),legend.position = "none")+
    theme(axis.line = element_line(color = "black"),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank())+
    scale_fill_manual(values = c(mypal[1],mypal[5],mypal[2]),name="")+
    scale_y_continuous(expand = c(0,0))+
    coord_fixed(ratio=1/4)+
    geom_vline(xintercept =89.5,color="white")
  
  
  #################################################################
  #PCA for tea genetics
  library(dplyr)
  tea.pca <- read.table("D:/2.tea/tea_resCA.eigenvec", quote="\"", comment.char="")
  tea.pca<-tea.pca[order(tea.pca$V1),]
  
  pca.eig <- read.table("D:/2.tea/tea_resCA.eigenval", quote="\"", comment.char="")
  
  
  #  % Variation explained
  
  pca.eig<-mutate(pca.eig,V2=round(V1*100/sum(V1),2))


  tea.pca<-(tea.pca,Group=tea.group$group2)
  tea.pca1<-mutate(tea.pca1,group.hclust2=clclust.id$group,group.hclust3=clclust.id$group3)

ylab=paste("PC2(",pca.eig[2,2],"%)",sep="")
zlab=paste("PC3(",pca.eig[3,2],"%)",sep="")


pca<-ggplot(tea.pca1, aes(x =V3,
                        y = V4,color=factor(group.hclust3)))+
  geom_point(size=2,alpha=0.5)+
  scale_color_manual(values = c(mypal[1],mypal[5],mypal[2]))+
  xlab(xlab)+
  ylab(ylab)+
  theme(legend.position  = "top",panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color="black"))+ 
  scale_y_continuous(expand = c(0.05,0.05))+ coord_fixed(ratio=2/2.8)