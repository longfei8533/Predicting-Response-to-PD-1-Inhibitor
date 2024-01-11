library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(UniProt.ws)
library(enrichplot)
library(ComplexHeatmap)
library(circlize)

res <- read.csv(file="../results/diff_result.csv")
res_filter <- res[res$pvalue< 0.05,]

up_pro <- res_filter[res_filter$ratio > 1.2,]
down_pro <- res_filter[res_filter$ratio < 1/1.2,]

# ------------------------------------- KEGG --------------------------------

kegg_all <- enrichKEGG(c(up_pro$Protein,down_pro$Protein),
                       organism = "hsa",
                       keyType = "uniprot")

xx<- kegg_all@result$geneID[1] %>% strsplit(split="/") %>% .[[1]]
xx_up <- xx[xx %in% up_pro$Protein]
xx_down <- xx[xx %in% down_pro$Protein]

kegg_plot <- function(filename,data){
  data <- data[1:10,]
  data$Ratio <- sapply(data$GeneRatio,function(x){
    n <- strsplit(x,split = "[/]") %>% unlist()
    r <- as.numeric(n[1])/as.numeric(n[2])
    r
  })
  data$log2_P <- -log10(data$pvalue)
  data$Description <- 
    factor(data$Description,levels = data$Description[order(data$Ratio)])
  
  p <- ggplot(data,aes(x = Ratio , y = Description)) +
    geom_point(aes(size = Count, color = log2_P))+
    scale_color_gradient(low = "#339999",high = "#fc2e00")+
    labs(x = "Ratio",
         y = "",
         size = "Gene Count",
         color = expression(paste(-log[10]," P value"))
    )+
    theme_bw()+
    theme(axis.text = element_text(color = "black")
    )
  pdf(filename,width = 6,height =3.2)
  print(p)
  dev.off()
}

kegg_plot("../results/kegg_all.pdf",kegg_all@result)

# ------------------------------------- GO --------------------------------
ego <- enrichGO(c(up_pro$Gene,down_pro$Gene),
                   OrgDb = "org.Hs.eg.db" ,
                   keyType = "SYMBOL",
                   ont = "ALL"
                   )
ego_bp <- ego@result %>% dplyr::filter(ONTOLOGY=="BP") %>% 
  arrange(pvalue) %>% .[1:10,] %>% arrange(desc(pvalue)) %>% 
  dplyr::mutate(ONTOLOGY = "Biological process")

ego_bp$Description <- stringr::str_to_title(ego_bp$Description)
ego_bp$Description <- factor(ego_bp$Description,levels = ego_bp$Description)

go_bp_plot <- ggplot(data=ego_bp,aes(y = Description, x = -log10(pvalue))) +
  geom_col(aes(fill=Count))+
  scale_fill_gradient(low = c("#1874CD"),high = "#fc2e00")+
  labs(x = expression(paste(-log[10]," P value")), y = "",fill = "Gene Count")+
  theme_bw()

pdf("../results/Go_bp.pdf",width = 7,height =3)
print(go_bp_plot)
dev.off()

# ------------------------------ KEGG GSEA  GO --------------------------------

gsea_dat <- read.csv("../results/diff_result.csv")
gsea_dat$rank_value <- -(log2(gsea_dat$ratio) * log10(gsea_dat$pvalue))
gsea_dat <- dplyr::arrange(gsea_dat,desc(rank_value))
geneList <- gsea_dat$rank_value
names(geneList) <- gsea_dat$Gene

kegg_gmt <- read.gmt("../data/GSEA/c2.cp.kegg.v2022.1.Hs.symbols.gmt")
kegg_gsea_res <- GSEA(geneList,TERM2GENE = kegg_gmt)

pdf("../results/kegg gsea.pdf",width = 8,height =6)
ridgeplot(kegg_gsea_res,10)
dev.off()

#  --------------------------------- CCCP heatmap -----------------------------
dat <- read.csv(file = "../processed_data/NAguideR_seqknn.csv")
dat_mat <- dat[,-1] %>% as.matrix()
rownames(dat_mat) <- dat$X
sample_info <- read.csv("../processed_data/sample_data.csv")

ccc_protein <- kegg_all@result$geneID[1] %>% strsplit(split="/") %>% .[[1]]

xx <- select(org.Hs.eg.db,keys = ccc_protein, keytype = "UNIPROT", column = "ENTREZID" ) 
cat(xx$ENTREZID)


ccc_gene <- select(org.Hs.eg.db,keys = ccc_protein, keytype = "UNIPROT", column = "SYMBOL" ) 


heatplot_dat <- dat_mat[match(ccc_gene$UNIPROT,rownames(dat_mat)),] %>% 
  log2() %>% t() %>% scale %>% t()
rownames(heatplot_dat) <- ccc_gene$SYMBOL

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

png(filename = "../results/CCC_heatplot.png",width = 12,height = 10,units = "cm",res = 300)
draw(heat_plot, heatmap_legend_side = "bottom")
dev.off()




