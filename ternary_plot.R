
#Dependencies

library(ggtern)
library(tidyverse)
library(edgeR)
####################################################################################

#Data reformation function
data_clean <- function(otu, design, type=c("relative", "absolute"), threshold=0.0001, times=100){
  
  # ????????????
  # library(amplicon)
  # otu=otutab
  # metadata$SampleID=rownames(metadata)
  # design=metadata[,c("SampleID","Group")]
  # type="absolute"
  # threshold=0.0005
  # times=100
  
  if (type == "absolute"){
    otu_relative <- apply(otu, 2, function(x){x/sum(x)})
  }else {otu_relative <- otu}
  
  idx <- rowSums(otu_relative > threshold) >= 1
  otu_threshold <- as.data.frame(otu_relative[idx, ])
  otu_threshold$OTUs <- row.names(otu_threshold)
  
  otu_longer <- pivot_longer(data=otu_threshold, 
                             cols=-OTUs,
                             names_to="SampleID", 
                             values_to="value")
  
  merge_data <- merge(otu_longer, design, by ="SampleID")
  # otu <- subset(merge_data, select=-SampleID)
  otu <- subset(merge_data, select=c("Group","OTUs","value"))
  otu_mean <- otu %>% group_by(OTUs, Group) %>% 
    summarise(value=mean(value))
  otu_tern <- otu_mean %>%
    group_by(Group, OTUs) %>%
    mutate(index=row_number()) %>%
    pivot_wider(names_from=Group,values_from=value) %>%
    select(-index)

  otu_tern$size <- (apply(otu_tern[2:4], 1, mean))*times   
  return(otu_tern)
}





####################################################################################

# enrich index calculation function
enrich_data <- function(otu, design, p.value=0.05, adjust.method="bonferroni"){
  
  # library(amplicon)
  # otu=otutab
  # metadata$SampleID=rownames(metadata)
  # design=metadata[,c("SampleID","Group")]
  # p.value=0.05
  # adjust.method="fdr" bonferroni
  
  dge_list <- DGEList(counts=otu, group=design$Group)
  # Remove the lower abundance/(cpm, rpkm)
  keep <- rowSums(dge_list$counts) >= 0
  dge_keep <- dge_list[keep, ,keep.lib.sizes=F]
  # scale the raw library sizes dgelist
  dge <- calcNormFactors(dge_keep)
  # fit the GLM
  design.mat <- model.matrix(~ 0 + dge$samples$group)
  d2 <- estimateGLMCommonDisp(dge, design.mat)
  d2 <- estimateGLMTagwiseDisp(d2, design.mat)
  fit <- glmFit(d2, design.mat)
  #######
  # if (missing(adjust.method))
  #   adjust.method="fdr"
  # if (missing(p.value))
  #   p.value=0.05
  group_index <- as.character(design$Group[!duplicated(design$Group)])
  # enrich groups
  lrt_1_2 <- glmLRT(fit, contrast=c(1, -1, 0))
  lrt_1_3 <- glmLRT(fit, contrast=c(1, 0, -1))
  
  de_1_2 <- decideTestsDGE(lrt_1_2, adjust.method=adjust.method, 
                           p.value=p.value)
  de_1_3 <- decideTestsDGE(lrt_1_3, adjust.method=adjust.method, 
                           p.value=p.value)
  
  rich_1 <- rownames(otu)[de_1_2 == 1 & de_1_3 == 1]
  enrich_1 <- data.frame(OTUs=rich_1, 
                         enrich=rep(group_index[1], length(rich_1)))
  ###############################
  lrt_2_3 <- glmLRT(fit, contrast=c(0, 1, -1))
  lrt_2_1 <- glmLRT(fit, contrast=c(-1, 1, 0))
  
  de_2_3 <- decideTestsDGE(lrt_2_3, adjust.method=adjust.method, 
                           p.value=p.value)
  de_2_1 <- decideTestsDGE(lrt_2_1, adjust.method=adjust.method, 
                           p.value=p.value)
  
  rich_2 <- rownames(otu)[de_2_3 == 1 & de_2_1 == 1]
  enrich_2 <- data.frame(OTUs=rich_2, 
                         enrich=rep(group_index[2], length(rich_2)))
  ###################
  lrt_3_1 <- glmLRT(fit, contrast=c(-1, 0, 1))
  lrt_3_2 <- glmLRT(fit, contrast=c(0, -1, 1))
  
  de_3_1 <- decideTestsDGE(lrt_3_1, adjust.method=adjust.method, 
                           p.value=p.value)
  de_3_2 <- decideTestsDGE(lrt_3_2, adjust.method=adjust.method, 
                           p.value=p.value)
  
  rich_3 <- rownames(otu)[de_3_1 == 1 & de_3_2 == 1]
  enrich_3 <- data.frame(OTUs=rich_3, 
                         enrich=rep(group_index[3], length(rich_3)))
  enrich_index <- rbind(enrich_1, enrich_2, enrich_3)
  return(enrich_index)
}


########################################################################################################


#genus ternary

#Rhizo-root-phyllosphere
g.rrl<-bac.genus.norm[,1:300]

g.rrl<-g.rrl[,order(colnames(g.rrl))]
g.rrl<-g.rrl[order(rownames(g.rrl)),]

g_tern<- data_clean(g.brr, design, type="absolute", threshold=0.0001, times=100)

enrich_index <- enrich_data(g.brr, design, p.value=0.05)
plot_data <- merge(g_tern, enrich_index, by="OTUs", all.x=T)
library(stringr)
aaa<-plot_data[,6]
enrich<-str_replace_na(aaa, replacement = "NA")
plot_data<-data.frame(plot_data[,-6],enrich)
plot_data<-plot_data[order(plot_data$Genus),]

#Bulk-rhizo-root
g.brr<-bac.genus.norm[,101:400]

g.brr<-g.brr[,order(colnames(g.brr))]
g.brr<-g.brr[order(rownames(g.brr)),]

g_tern1<- data_clean(g.rrl, design1, type="absolute", threshold=0.0001, times=100)

enrich_index1 <- enrich_data(g.rrl, design1, p.value=0.05)
plot_data1 <- merge(g_tern1, enrich_index1, by="Genus", all.x=T)
library(stringr)
aaa<-plot_data1[,6]
enrich1<-str_replace_na(aaa, replacement = "NA")
plot_data1<-data.frame(plot_data1[,-6],enrich1)
plot_data1<-plot_data1[order(plot_data1$Genus),]
####################################################################################


#Ternary Plot, plot_data,plot_data1
p.1<- ggtern(data=plot_data,
              aes(x=Phyllosphere, y=Root,z=Rhizo))+
  geom_mask() +
  geom_point(aes(size=size, 
                 color=Phylum),alpha=0.8) + 
  scale_color_manual(values = colors) +
  theme_bw() +
  theme(legend.title = element_blank(),axis.text=element_blank(), 
        axis.ticks=element_blank())+
  theme_void(base_size = 12)
p.1
ggsave("RRL_genus_tern_plot_e0.01.pdf", p.1, width=8, height=8,useDingbats=F)



####################################################################################
#Top enriched genus

#heatmap , plot_data,plot_data1
library(pheatmap)
g<-pheatmap(main = "Enrichment",plot_data,scale = "none",
            cluster_row = F, cluster_col = F, border_color = "white",
            display_numbers = F,
            fontsize_number = 7, number_color = "white",
            cellwidth = 30, cellheight =10,color=mycolor,
            annotation_row= annot_data)
dev.off()
ggsave(g,filename = "RRL_10Genus_heatmap.pdf",width = 7,height = 5,useDingbats=F)

##############
