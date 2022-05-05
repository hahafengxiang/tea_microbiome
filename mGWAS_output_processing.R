setwd("D:/2.tea/tea_BRRL_Total/5.gemma_202104/1.gemma_out/")
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)

#annotation files for tea reference genome
tea.annot<-readRDS("D:/2.tea/tea_resequencing/SNP-202012/tea_reseq2021_SNPs_Annot.rds")
tea.genes<-fread("tea_genes.csv",header = T)
tea.func<-fread("D:/2.tea/tea_resequencing/SNP-202012/0.genome/CSS_ChrLev_20200506_Function.fas")


##################################################################################################
# mGWAS visualization

#mGWAS for OTUs
# phyllosphere Otu2029, root Otu956, and rhizosphere Otu4847, take Otu4877 as example
library(rMVP)
o4877<-fread("Otu4877.txt.assoc.txt",header=T)
o4877<-o4877[,c(2,1,3,12)]
#Manhattan plot
MVP.Report(o4877, plot.type="m", 
           col=c("#DAA520","#B83F2F"), LOG10=TRUE, ylim=NULL,
           threshold=1.95e-12, threshold.lty=2, threshold.lwd=1,
           threshold.col=c("black", "grey"), amplify=F,chr.den.col=NULL,
           signal.cex=1,signal.pch=19,
           file.type="TIFF",memo="",dpi=600)
#qq plot

MVP.Report(o4877,plot.type="q",
           conf.int=F,box=TRUE,file.type="TIFF",memo="",dpi=600)

#######################################################
#mGWAS for microbial PCOs, take phyllosphere as example
a = list.files(pattern = ".txt.assoc.txt")                                     
leaf.multi<-fread(a[1],header=T) 
leaf.multi<-leaf.multi[,c(2,1,3,12)]
ss<-basename(a[1])
ss<-gsub('.txt.assoc.txt','',ss)
names(leaf.multi)[4]<-ss

for (i in 2:length(a)){
  new.data = fread(a[i],header=T)
  new.data<-new.data[,12]
  ss<-basename(a[i])
  ss<-gsub('.txt.assoc.txt','',ss)
  names(new.data)<-ss
  leaf.multi = cbind(leaf.multi,new.data)
}

MVP.Report(leaf.multi, plot.type="m", 
           multracks=TRUE,
           threshold=3.81e-9,threshold.lty=2, 
           threshold.lwd=c(1,1), threshold.col="grey", amplify=TRUE,bin.size=1e6,
           chr.den.col=NULL,col=c("#DAA520","#0072B2","#009E73","#CC79A7","#BEBEBE","#00A4DE","#2E3092"),
           file.type="tiff",memo="rhizo",dpi=600)

#multi trait qq plot
MVP.Report(leaf.multi,plot.type="q",col=c("#DAA520","#0072B2","#009E73","#CC79A7","#BEBEBE","#00A4DE","#2E3092"),
           signal.pch=19,signal.cex=1.5,signal.col="red",
           conf.int=T,conf.int.col="grey",box=FALSE,
           multracks= TRUE,
           file.type="tiff",memo="",dpi=300)

#Single trait visualization for leaf PCO1
leaf.pc1.assoc<-fread("leaf.pc1.txt.assoc.txt",header=T)
leaf.pc1.assoc<-leaf.pc1.assoc[,c(2,1,3,12)]
MVP.Report(leaf.pc1.assoc, plot.type="m", LOG10=TRUE, ylim=NULL,
           threshold=3.81e-9, threshold.lty=2, threshold.lwd=1,
           threshold.col="grey", amplify=F,chr.den.col=NULL,
           signal.cex=1,signal.pch=19,
           file.type="TIFF",memo="",dpi=600)

#qq plot
MVP.Report(leaf.pc1.assoc,plot.type="q",
           conf.int=F,box=TRUE,file.type="TIFF",memo="",dpi=600)


###################################################################################################
#mGWAS signal SNPs

#gather all signigicant outputs from mGWAS

a = list.files(pattern = ".txt.assoc.txt")                                     

assoc.data<-fread(a[1],header=T) 
assoc.data<-filter(assoc.data,p_wald<=1E-7) 
ss<-basename(a[1])
ss<-gsub('.txt.assoc.txt','',ss)
assoc.data<-mutate(assoc.data,OTU=rep(ss,length(assoc.data[,1])))

for (i in 2:length(a)){
  new.data = fread(a[i],header=T)
  new.data<-filter(new.data,p_wald<=1E-7) 
  ss<-basename(a[i])
  ss<-gsub('.txt.assoc.txt','',ss)
  new.data<-mutate(new.data,OTU=rep(ss,length(new.data[,1])))
  
  
  assoc.data = rbind(assoc.data,new.data)
}

otu.gemma<-rbind(rhizo.assoc.data,root.assoc.data,leaf.assoc.data)

otu.gemma1<-filter(otu.gemma,p_wald<=1.95E-12)

keep.rows<-which(tea.annot$RS%in%otu.gemma1$RS)
otu.annot<-tea.annot[keep.rows,]
otu.snps<-merge(otu.gemma1,otu.annot,by="RS")

#
library(tidyr)


#annotate the snps to tea genes
keep.rows<-which(tea.genes$gene%in%otu.snps$candidate)
otu.genes<-tea.genes[keep.rows,]

#add function annotations
keep.rows<-which(tea.func$GeneID%in%otu.genes$cds)
otu.go<-tea.func[keep.rows,]

names(otu.genes)<-c("GeneID","candidate")
otu.go<-merge(otu.genes,otu.go,by="GeneID")
otu.go<-merge(otu.go,otu.snps,by="candidate")


otu.go<-filter(otu.go,gPosType!="intronic")
otu.go<-filter(otu.go,dist<=21100)




###################################################################################################
#Taxon-Function network

aa<-aggregate(otu.go$GO,list(otu.go$OTU,otu.go$GO,otu.go$ASP),length)
names(aa)<-c("OTU","GO","ASP","Count")
aa<-filter(aa,GO=="Biological Process")
library(igraph)
otu.el1<-as.matrix(aa[,c(7,3)])
otu.el1<-as.matrix(otu.el1[!duplicated(otu.el1),])
otu.graph<-graph_from_edgelist(otu.el1,directed = F)
write.graph(otu.graph,file = "otu_otu_GO.graphml",format = "graphml")


##########################################################################
#GO Enrichment analysis

leaf.enrich<-fread("Leaf_TBtools.enrich.temp.GO.Enrichment.final.xls",header = T)
rhizo.enrich<-fread("Rhizo_TBtools.enrich.temp.GO.Enrichment.final.xls",header = T)
root.enrich<-fread("Root_TBtools.enrich.temp.GO.Enrichment.final.xls",header = T)

leaf.enrich1<-filter(leaf.enrich,Class=="Biological process")
rhizo.enrich1<-filter(rhizo.enrich,Class=="Biological process")
root.enrich1<-filter(root.enrich,Class=="Biological process")

leaf.enrich1<-arrange(leaf.enrich1,P_value)
rhizo.enrich1<-arrange(rhizo.enrich1,P_value)
root.enrich1<-arrange(root.enrich1,P_value)


#plot
rank1<-arrange(leaf.enrich1,P_value)

bar.plot<-ggplot(leaf.GO.enrich,aes(y=factor(GO_Name,levels = rev(rank1$GO_Name)),
                                    x=-log10(P_value),fill=Class))+
  geom_bar(stat  = "identity",width = 0.7)+
  geom_text(aes(label=factor(GO_Name,levels = rev(rank1$GO_Name))), vjust=0.3,hjust=-0.05, size=3.5)+
  theme_bw()+
  ylab("GO Term")+
  xlab("-log10(P-value)")+
  scale_y_discrete(guide = guide_axis(position = "left"),expand = c(0.03,0.03))+
  scale_x_continuous(expand = c(0,0.01),guide = guide_axis(position = "bottom"))+
  scale_fill_manual(values = pal1[c(4,2,1)])+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border = element_blank())+
  theme(panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black",angle = 0, hjust = 0.3,vjust = 0.3),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.length = unit(0.05, "cm"))
######################################################################################################################################################################################################
#Genotypes of SNPs in association with phyllosphere Otu2029, root Otu956, and rhizosphere Otu4847
#vcfR
library(vcfR)
a<-read.vcfR( "otu_gemma_peaks_snp.recode.vcf")

a1<-data.frame(a@fix)
a1<-mutate(a1,RS=paste0(CHROM,"_",POS))
a2<-data.frame(a@gt)
a2<-a2[,-1]

a2<-as.data.frame(lapply(a2, as.numeric))
rownames(a2)<-a1$RS

a3<-extract.gt(
  a,
  element = "GT",
  mask = FALSE,
  as.numeric = F,
  return.alleles = T,
  IDtoRowNames = TRUE,
  extract = TRUE,
  convertNA = TRUE
)
peak.gt<-t(a3)