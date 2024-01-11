library(tidyverse)
library(ggrepel)
library(ggstatsplot)
library(ggpubr)
library(org.Hs.eg.db)
library(UniProt.ws)
library(factoextra)
library(ComplexHeatmap)
library(circlize)

dat <- read.csv(file = "../processed_data/NAguideR_seqknn.csv")
dat_mat <- dat[,-1] %>% as.matrix()
rownames(dat_mat) <- dat$X
sample_info <- read.csv("../processed_data/sample_data.csv")


# -------------------------------- PCA ----------------------------------------

pca_dat <- prcomp(x=t(dat_mat),rank.=10,scale. = TRUE)

fviz_eig(pca_dat, addlabels = T,
                linecolor="#E64B35FF",
                barfill="#4DBBD5FF",
                barcolor="#4DBBD5FF",
                ggtheme=ggplot2::theme(
                  axis.title = element_text(face="bold"),
                  axis.text = element_text(face="bold",color = "black"),
                  axis.line = element_line(size=0.8,color="black"),
                  axis.ticks= element_line(size=0.8,colour = "black"),
                  panel.grid =element_blank(),
                  panel.background = element_blank(),
                  title = element_blank(),
                ),
                ylim=c(0,30)
)



pca.var <- pca_dat$sdev^2  
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
data_pca <- data.frame(
  Group=factor(sample_info$Label,levels = c(0,1),
               labels = c("Non-responders","Responders")),
  PC1 = pca_dat$x[,1],
  PC2 = pca_dat$x[,2],
  sample = 1:28 )

p_pca <- ggplot(data=data_pca,aes(x=PC1,y=PC2,color=Group,shape=Group))+
  geom_point(size=3)+
  geom_vline(xintercept = 0, linetype=2)+
  geom_hline(yintercept = 0, linetype=2)+
  #geom_text(aes(label=sample))+
  scale_color_manual(values = c("Non-responders"= "#BDB76B",
                                "Responders"= "#483D8B"))+
  xlab(paste0("PC1(",pca.var.per[1],"%"," ","variance)"))+
  ylab(paste0("PC2(",pca.var.per[2],"%"," ","variance)"))+
  theme_bw() +
  theme(legend.position = "top")

png(filename = "../results/PCA_plot.png",
    width =10 ,height =10 ,units = "cm",res = 300)
p_pca
dev.off()

# --------------------------  Difference analysis -----------------------------

prog_dat <- dat_mat[,sample_info$Label == 0]
resp_dat <- dat_mat[,sample_info$Label == 1]


diff_fun <- function(d1,d2){
  ratio <- c()
  pvalue <- c()
  for(i in 1:nrow(d1)){
    x1 <- as.numeric(d1[i,])
    x2 <- as.numeric(d2[i,])
    x1 <- x1[x1!=0]
    x2 <- x2[x2!=0]
    if(length(x1)<2 | length(x2)<2){
      r <- NA
      p <- NA
    }else{
      r <- mean(x1)/mean(x2)
      p <- t.test(log2(x1),log2(x2),var.equal = T)$p.value
    }
    ratio <- c(ratio,r)
    pvalue <- c(pvalue,p)
  }
  return(data.frame(ratio,pvalue))
}
diff_res <- diff_fun(prog_dat,resp_dat)
diff_res$fdr <- p.adjust(diff_res$pvalue,method ="BH")

up <- UniProt.ws()
genes <- select(up, keys=rownames(dat_mat), keytype="UniProtKB", column="gene_primary", multiVals="first")

res <- data.frame(Protein = rownames(dat_mat), Gene = genes$Gene.Names..primary. ,diff_res)
write.csv(res,file="../results/diff_result.csv",row.names = F)

# ----------------------------- Volcano plot -----------------------------------

res <- read.csv(file="../results/diff_result.csv")
plot_dat <- res[complete.cases(res),]
plot_dat$fill <- "NS"
plot_dat$fill[plot_dat$pvalue < 0.05 & plot_dat$ratio < (1/1.2)] <- "Down"
plot_dat$fill[plot_dat$pvalue < 0.05 & plot_dat$ratio > 1.2] <- "Up"

plot_txt <- plot_dat[plot_dat$fill != "NS",]

res_plot <- ggplot(data=plot_dat,aes(x=log2(ratio),y=-log10(pvalue),fill = fill))+
  geom_point(alpha=1, size=3,shape = 21,show.legend = F) +
  scale_fill_manual(values = c(Up = c("#EE6A50"), NS = c("#BABABA") , Down = c("#009ACD")))+
  labs(x="log2(Fold change)",
       y="-log10 (P value)")+
  theme_bw()


png(filename = "../results/volcano plot.png",width =10 ,height = 10,units = "cm",res = 300)
res_plot
dev.off()

# ------------------------------ Heatmap plot --------------------------------



res_filter <- res[res$pvalue< 0.01 & (res$ratio > 1.5 | res$ratio < 1/1.5),]

heatplot_dat <- dat_mat[match(res_filter$Protein,rownames(dat_mat)),] %>% 
  log2() %>% t() %>% scale %>% t()
rownames(heatplot_dat) <- res_filter$Gene

a_top <- HeatmapAnnotation(`Clinical outcomes` = c(rep("Non-responders",11),rep("Responders",17)),
                           col = list(`Clinical outcomes` = c("Non-responders"= c("#7A8B8B"), "Responders"=c("#6E8B3D"))),
                           show_legend = F,show_annotation_name = F
)


col_fun = colorRamp2(c(-2, 0, 2),c("#009ACD", "#FFFFFF", "#EE6363")) 

heat_plot <- Heatmap(heatplot_dat,
                     #name = "Relative abundance",
                     na_col = "gray",
                     col = col_fun,
                     heatmap_legend_param = list(col_fun = col_fun, title = "Relative abundance", 
                                                 at = c(-2,0,2),
                                                 labels = c("low", "median", "high"), 
                                                 direction = "horizontal",
                                                 heatmap_legend_side = "bottom"
                     ),
                     top_annotation = a_top,
                     row_names_gp = gpar(fontsize = 8),
                     cluster_columns = T,
                     cluster_rows = T,
                     show_column_names = F,
                     column_split = c(rep("Non-responders",11),rep("Responders",17))
)


png(filename = "../results/heatmap plot.png",width = 16,height = 12,units = "cm",res = 300)
draw(heat_plot, heatmap_legend_side = "bottom")
dev.off()



