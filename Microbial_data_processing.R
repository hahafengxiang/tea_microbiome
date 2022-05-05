#Dependencies and colors

library(phyloseq)
library(vegan)
library(dplyr)
library(tidyr)
library(ggsci)
library(agricolae)
library(scales)
library(ggplot2)
library(ape)

pal2<-c("#BEBEBE","#DAA520","#0072B2","#3C6950","#B83F2F","#FFF148","#2E3092")
mypal =pal_d3("category20")(20)
########################################################################################################################
#16S sequencing data processing

source("sintax_info.R")
taxa.tax <- read.table("otu.sintax", sep = "\t")

deteleotu <-taxa.tax[taxa.tax$Class=="c: Chloroplast",]
taxa.tax$ID <- row.names(taxa.tax)
taxa.tax[taxa.tax$ID %in% row.names(deteleotu),]<-NA
taxa.tax<-as.matrix(taxa.tax)
taxa1<-na.omit(taxa.tax)
taxa1 <-taxa1[,-7]

otutab_bac <- read.delim("otutab.txt",row.names = 1)
otutab_bac$ID <- row.names(otutab_bac)
otutab_bac[otutab_bac$ID %in% row.names(deteleotu),]<-NA
otutab_bac1<-na.omit(otutab_bac)


#2. make phyloseq object
library(phyloseq)
bac.phylo <- phyloseq(otu_table(bac.otu3,
                                taxa_are_rows = TRUE),
                      tax_table(taxa1),
                      sample_data(env400),
                      phy_tree(bac.tree))

library(DESeq2)

bac.deseq <- phyloseq_to_deseq2(physeq = bac.phylo,
                                design = ~position)

bac.deseq.wald <- DESeq(bac.deseq,
                        fitType = "parametric",
                        test = "Wald")

bac.norm <- phyloseq(otu_table(counts(bac.deseq.wald, 
                                      normalized = TRUE),
                               taxa_are_rows = TRUE),
                     tax_table(taxa1),
                     sample_data(env400),
                     phy_tree(bac.tree))

##########################################################################################################################

#Fig. 1a PCoA for all compartments
#6. pcoa
bac.pcoa <- ordinate(bac.norm,method = c("PCoA"), distance = "wunifrac")
ord.df <- data.frame(bac.pcoa$vectors[,1:3],sample_names(bac.norm),sample_data(bac.norm))
result.bray <-bac.pcoa$values[,"Relative_eig"]
PCo1 = as.numeric(sprintf("%.4f",result.bray[1]))*100
PCo2 = as.numeric(sprintf("%.4f",result.bray[2]))*100
PCo3 = as.numeric(sprintf("%.4f",result.bray[3]))*100
xlab=paste("PCOA1(",PCo1,"%)",sep="") 
ylab=paste("PCOA2(",PCo2,"%)",sep="")
zlab=paste("PCOA3(",PCo3,"%)",sep="")


pca<-ggplot(ord.df, aes(x = Axis.1,
                        y = Axis.2,
                        color =position))+
  geom_point(size = 2,alpha=0.4)+
  theme_bw()+  xlab(xlab)+
  ylab(ylab)+  theme(legend.position = c("top"))
pca
pca+
  scale_y_continuous(expand = c(0,0.005))+scale_color_manual(values = pal2[c(1,4,2,3)])+
  theme(panel.background = element_blank(),panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  theme(legend.position = "top")

#Zoom-in view for Rhizosphere, bulk soil and root endosphere
ggplot(ord.df[101:400,], aes(x = Axis.1,
                             y = Axis.2,
                             color =position))+
  geom_point(size = 2,alpha=0.6)+
  theme_bw()+  xlab(xlab)+
  ylab(ylab)+  theme(legend.position = c("top"))+
  scale_y_continuous(expand = c(0,0.005))+scale_color_manual(values = pal2[c(1,2,3)])+
  theme(panel.background = element_blank(),panel.grid = element_blank(),panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  theme(legend.position = "top")
##########################################################################################################################
#Fig. 1b ANOSIM between compartments
total.dist<-distance(bac.norm, method="wunifrac", type="samples")
total.anosim<-with(env400, anosim(total.dist, factor(position)),permutations=999)
summary(total.anosim)
anosim.df<-data.frame(rank=total.anosim$dis.rank,group=total.anosim$class.vec)
anosim.df<-filter(anosim.df,group!="Between")
result=paste("R=",round(total.anosim$statistic,2),"P=", total.anosim$signif)

p1<-ggplot(anosim.df, aes(x=factor(group,levels = c("Bulk","Rhizo","Root","Phyllosphere")), y=rank, color=group)) +
  geom_boxplot(notch = TRUE)+
  xlab("Compartments")+ylab("Weighted-UniFrac Rank")+
  ggtitle("ANOSIM R = 0.46, P = 0.001 ")
p1+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.background = element_blank(),panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(colour = "black",size = 8,angle = 45),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black",size = 8,angle = 45))+ 
  scale_color_manual(values = pal2[c(1,4,2,3)])+
  theme(legend.position = "none")+
  scale_y_continuous(expand = c(0.1,0))

#Pairwise anosim
ss<-mutate(env400,LRRB=c(rep("L",100),rep("RRB",300)))
anosim.lrrb<-with(ss, anosim(total.dist, factor(LRRB)),permutations=999)
summary(anosim.lrrb)
##########################################################################################################################
#Fig. 1c Alpha diversity
#shannon
div.df <- diversity(t(bac.otu))
div.df<-data.frame(position=env.df,total=shannon.total)
names(div.df)<-c("Compartments","Total")
div.df1<-aggregate(.~cultivar,mean,data = div.df)

#one-way anova
shannon.aov<-aov(Total~Compartments,data = div.df1)
result_1 <- HSD.test(shannon.aov, "Compartments", group = T)
print(result_1)

p<-ggplot(div.df1, aes(x=factor(Compartments,levels = c("BulkSoil","Rhizosphere","Root","Phyllosphere")),
                       y=Total)) +
  geom_boxplot(aes(colour=Compartments),notch = TRUE,outlier.colour = NA)+
  xlab("Compartments")+ylab("Shannon Index")+
  geom_jitter(size=2,shape=16, position=position_jitter(0.1),aes(colour=Compartments),alpha=0.5)
p+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.background = element_blank(),panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(colour = "black",size = 8,angle = 45),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.y = element_text(colour = "black",size = 8))+ 
  scale_color_manual(values = pal2[c(1,4,2,3)])+
  theme(legend.position = "none")+
  scale_y_continuous(expand = c(0.1,0)) 
##########################################################################################################################
#Fig. 1 d-g compartment specific microbial composition at phylum level
leaf.phylum<-filter(phylum.df,position=="Phyllosphere")
root.phylum<-filter(phylum.df,position=="Root")
rhizo.phylum<-filter(phylum.df,position=="Rhizo")
bulk.phylum<-filter(phylum.df,position=="Bulk")


#cultivars were ranked by Proteobacteria abundance, take phyllosphere as an example
aaa<-filter(leaf.phylum,phylum=="Proteobacteria")
rank<-arrange(aaa,value)

phylum.plot<-ggplot(leaf.phylum,aes(x=factor(cultivar,levels = rank$cultivar),
                                    y=value*100,fill=phylum))+
  geom_bar(stat = "identity",width = 0.8)+
  xlab("100 Cultivars")+
  ylab("Relative abundance (%)")
phylum.plot+scale_fill_manual(values = c(pal2[c(5,2,4,7,6)],mypal[16],pal2[c(3,1)]),name="")+
  theme(axis.text.x = element_text(hjust = 0.5,vjust = 1))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.background = element_blank(),panel.border = element_blank(),
        axis.line = element_line(color = "black"),axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(colour = "black",size = 8))+
  scale_y_continuous(expand = c(0,0))+
  ylab("Relative abundance (%)")

##########################################################################################################################
#Supplementary Fig. 1: phylum composition compared between compartments
phylum.comp<-ggplot(phylum.df,aes(x=factor(position,
                                            levels = c("Bulk","Rhizo","Root","Phyllosphere")),
                                   y=value*100))+
  geom_boxplot(aes(colour=position),notch = F)+
  xlab("")+
  ylab("Relative abundance (%)")+facet_wrap(~factor(phylum),scales = "free")

phylum.comp
phylum.comp+scale_color_manual(values = c(pal2[c(1,4,2,3)],mypal[16],pal2[c(3,1)]),name="")+
  theme(axis.text.x = element_text(hjust = 0.5,vjust = 1))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.background = element_blank(),panel.border = element_blank(),
        axis.line = element_line(color = "black"),axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(colour = "black",size = 8))+
  scale_y_continuous(expand = c(0,0))+
  ylab("Relative abundance (%)")

##########################################################################################################################
#Fig. 2c microbial beta diversity between tea subpopulations

#ANONISM
set.seed(315)
summary(with(tea.clust, anosim(leaf.dist, factor(group3)),permutations=999))
set.seed(315)
summary(with(tea.clust, anosim(rhizo.dist, factor(group3)),permutations=999))
set.seed(315)
summary(with(tea.clust, anosim(root.dist, factor(group3)),permutations=999))

#Compartment specific PCoA: root rhizo phyllosphere, take phyllosphere as an example

pcoa<-leaf.pcoa
phylo<-leaf.phylo1
pcoa.df <- data.frame(tea.clust,pcoa$vectors[,1:4],
                      sample_names(phylo),sample_data(phylo))


result.bray <-pcoa$values[,"Relative_eig"]
PCo1 = as.numeric(sprintf("%.4f",result.bray[1]))*100
PCo2 = as.numeric(sprintf("%.4f",result.bray[2]))*100
PCo3 = as.numeric(sprintf("%.4f",result.bray[3]))*100
xlab=paste("PCOA1(",PCo1,"%)",sep="") 
ylab=paste("PCOA2(",PCo2,"%)",sep="")
zlab=paste("PCOA3(",PCo3,"%)",sep="")


pca<-ggplot(pcoa.df, aes(x = Axis.1,
                         y = Axis.2,
                         color =factor(group3)))+
  geom_point(size = 2)+
  theme_bw()+xlab(xlab)+ylab(ylab)+ 
  stat_bag(prop = 1, alpha = 0, show.legend = F,aes(color = factor(group3)),linetype=2)
pca+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.background = element_blank(),
        axis.line = element_line(color = "black"),axis.text = element_text(size=8,color = "black"),
        axis.title = element_text(size = 8),
        axis.ticks.length = unit(0.05, "cm"),
        legend.position = "none")+
  scale_color_manual(values = c(mypal[1],mypal[2],mypal[5]))


