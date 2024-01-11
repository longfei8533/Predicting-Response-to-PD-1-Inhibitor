library(tidyverse)
library(org.Hs.eg.db)
library(UniProt.ws)
library(GSVA)
library(pROC)
library(data.table)


cell_dat <- readLines("../data/CellReports.txt")
geneSet <- list()
for(i in cell_dat){
  itms <- strsplit(i,"\t") %>% .[[1]]
  geneSet[[itms[1]]] <- itms[3:length(itms)]
  print(itms[1])
}

kegg_dat <- readLines("../data/GSEA/c2.cp.kegg.v2022.1.Hs.symbols.gmt")
kegg_geneSet <- list()
for(i in kegg_dat){
  itms <- strsplit(i,"\t") %>% .[[1]]
  kegg_geneSet[[itms[1]]] <- itms[3:length(itms)]
  print(itms[1])
}

geneSet[["KEGG_COMPLEMENT_AND_COAGULATION_CASCADES"]] <- 
  kegg_geneSet[["KEGG_COMPLEMENT_AND_COAGULATION_CASCADES"]]


## --------------------------------- GSVA -------------------------------------
anno_dat <- fread(file = "../data/Annotation_combine.txt") %>% as.data.frame()
sampleDat <- read.csv("../processed_data/sample_data.csv")
dat <- fread("../processed_data/NAguideR_seqknn.csv") %>% as.data.frame()
genes <- anno_dat$`Gene name`[match(dat[,1],anno_dat$`Protein accession`)]
dat <- dat[,-1]
dat <- log2(dat)
dat <- data.frame(genes = genes,dat)

ssdat <- dat[!duplicated(dat$genes),]
ssmat <- ssdat[,-1] %>% as.matrix()
rownames(ssmat) <- ssdat$genes
ss_res <- gsva(ssmat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=T)

# val

val_sample <- read.csv("../processed_data/validation_sample_data.csv")
val_dat <- read.csv("../processed_data/validation_data_seqknn.csv") 
val_dat <- separate_longer_delim(val_dat,X,";")
val_dat_mat <- val_dat[,-1] %>% as.matrix() %>% log2()
row.names(val_dat_mat) <- val_dat$X

up <- UniProt.ws()
genes <- select(up, keys=row.names(val_dat_mat), 
                keytype="UniProtKB", 
                column="gene_primary", 
                multiVals="first")

row.names(val_dat_mat) <- genes$Gene.Names..primary.
val_ss_res <- gsva(val_dat_mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=T)


# ------------------ COMPLEMENT_AND_COAGULATION_CASCADES (CCC) ----------------

plot_dat <- data.frame("Score" = ss_res["KEGG_COMPLEMENT_AND_COAGULATION_CASCADES",], 
                       "Label" = sampleDat$Label,
                       check.names = F
)
plot_dat$Data <- "Our data"
val_plot_dat <- data.frame("Score" = val_ss_res["KEGG_COMPLEMENT_AND_COAGULATION_CASCADES",], 
                           "Label" = val_sample$Label,
                           check.names = F
)
val_plot_dat$Data <- "Validation data"

roc_res <- roc(plot_dat$Label, plot_dat$Score)
val_roc_res <- roc(val_plot_dat$Label, val_plot_dat$Score)

# roc plot

png(filename = "../results/ccc_roc.png",
    width =12 ,height =12 ,units = "cm",res = 300)
plot(roc_res, main = "Complement and coagulation cascades",
     col = "#3CB371", 
     lwd = 3, 
     print.auc=TRUE,
     print.auc.y=0.4
)
plot(val_roc_res,add=TRUE,
     col= "#9A32CD",
     lwd = 3,
     print.auc = TRUE,
     print.auc.y=0.5
)

legend("bottomright", legend=c("Our data","Validation data"),
       col=c("#3CB371","#9A32CD"),lwd=3)
dev.off()

# box plot
box_plot_dat <- rbind(plot_dat,val_plot_dat)
box_plot_dat$Label <- ifelse(box_plot_dat$Label,"Responders","Non-responders")
box_plot_dat$Label <- factor(box_plot_dat$Label)

t.test(box_plot_dat$Score[box_plot_dat$Label == "Responders" & box_plot_dat$Data == "Our data"],
       box_plot_dat$Score[box_plot_dat$Label == "Non-responders" & box_plot_dat$Data == "Our data"]
       )
t.test(box_plot_dat$Score[box_plot_dat$Label == "Responders" & box_plot_dat$Data == "Validation data"],
       box_plot_dat$Score[box_plot_dat$Label == "Non-responders" & box_plot_dat$Data == "Validation data"]
)

png(filename = "../results/ccc_box.png",
    width =9 ,height =10 ,units = "cm",res = 300)
ggplot(box_plot_dat,aes(x = Data, y = Score))+
  geom_boxplot(aes(fill = Label),position = position_dodge(0.5), width = 0.5, outlier.alpha = 0) +
  scale_fill_manual(values = c("#CDAA7D", "#6495ED"))+
  labs(x = "", fill = "")+
  theme_bw()+
  theme( axis.text.x = element_text(angle = 45, hjust = 1 ),
         legend.position = "top"
  )

dev.off()

## ----------------------- Tumor microenvironment (TME) ------------------------

plot_dat <- data.frame("Score" = ss_res["Activated CD8 T cell",], 
                       "Label" = sampleDat$Label,
                       check.names = F
                       )
plot_dat$Data <- "Our data"
val_plot_dat <- data.frame("Score" = val_ss_res["Activated CD8 T cell",], 
                       "Label" = val_sample$Label,
                       check.names = F
)
val_plot_dat$Data <- "Validation data"

roc_res <- roc(plot_dat$Label, plot_dat$Score)
val_roc_res <- roc(val_plot_dat$Label, val_plot_dat$Score)

# roc plot
png(filename = "../results/tme_roc.png",
    width =12 ,height =12 ,units = "cm",res = 300)
plot(roc_res, main = "Activated CD8 T cell",
     col = "#3CB371", 
     lwd = 3, 
     print.auc=TRUE,
     print.auc.y=0.4
)
plot(val_roc_res,add=TRUE,
     col= "#9A32CD",
     lwd = 3,
     print.auc = TRUE,
     print.auc.y=0.5
)

legend("bottomright", legend=c("Our data","Validation data"),
       col=c("#3CB371","#9A32CD"),lwd=3)
dev.off()

# box plot
box_plot_dat <- rbind(plot_dat,val_plot_dat)
box_plot_dat$Label <- ifelse(box_plot_dat$Label,"Responders","Non-responders")
box_plot_dat$Label <- factor(box_plot_dat$Label)

t.test(box_plot_dat$Score[box_plot_dat$Label == "Responders" & box_plot_dat$Data == "Our data"],
       box_plot_dat$Score[box_plot_dat$Label == "Non-responders" & box_plot_dat$Data == "Our data"]
)
t.test(box_plot_dat$Score[box_plot_dat$Label == "Responders" & box_plot_dat$Data == "Validation data"],
       box_plot_dat$Score[box_plot_dat$Label == "Non-responders" & box_plot_dat$Data == "Validation data"]
)

png(filename = "../results/tme_box.png",
    width =9 ,height =10 ,units = "cm",res = 300)
ggplot(box_plot_dat,aes(x = Data, y = Score))+
  geom_boxplot(aes(fill = Label),position = position_dodge(0.5), width = 0.5, outlier.alpha = 0) +
  scale_fill_manual(values = c("#CDAA7D", "#6495ED"))+
  #coord_cartesian(ylim = c(-3,3) )+
  labs(x = "", fill = "")+
  #facet_grid(Data~.) +
  theme_bw()+
  theme( axis.text.x = element_text(angle = 45, hjust = 1 ),
         legend.position = "top"
  )

dev.off()


# ----------------------------- Single gene ------------------------------------

my_dat <- read.csv(file = "../processed_data/my_dat_ml_diff.csv") %>% as.matrix()
val_dat <- read.csv(file = "../processed_data/val_dat_ml_diff.csv") %>% as.matrix()

my_sample <- read.csv(file = "../processed_data/sample_data.csv")
val_sample <- read.csv(file = "../processed_data/validation_sample_data.csv")

res <- read.csv(file="../results/diff_result.csv")
res_filter <- res[res$pvalue< 0.05 & (res$ratio > 2 | res$ratio < 0.5),]
res_filter <- res_filter[res_filter$Protein %in% colnames(val_dat),]


for (i in 1:nrow(res_filter)){
  roc_res <- roc(my_sample$Label,my_dat[,res_filter$Protein[i]])
  val_roc_res <- roc(val_sample$Label,val_dat[,res_filter$Protein[i]])
  png(filename = paste0("../results/gene_roc/",res_filter$Gene[i],".png"),
      width =12 ,height =12 ,units = "cm",res = 300)
  plot(roc_res, main = res_filter$Gene[i],
       col = "#3CB371", 
       lwd = 3, 
       print.auc=TRUE,
       print.auc.y=0.4
  )
  plot(val_roc_res,add=TRUE,
       col= "#9A32CD",
       lwd = 3,
       print.auc = TRUE,
       print.auc.y=0.5
  )
  legend("bottomright", legend=c("Our data","Validation data"),
         col=c("#3CB371","#9A32CD"),lwd=3)
  dev.off()
}













