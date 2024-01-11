library(data.table)
library(tidyverse)
library(GSVA)
library(UniProt.ws)

cell_dat <- readLines("../data/CellReports.txt")

geneSet <- list()
for(i in cell_dat){
  itms <- strsplit(i,"\t") %>% .[[1]]
  geneSet[[itms[1]]] <- itms[3:length(itms)]
  print(itms[1])
}

# ------------------------------------- ssGSEA --------------------------------

anno_dat <- fread(file = "../data/Annotation_combine.txt") %>% as.data.frame()
sample <- read.csv("../processed_data/sample_data.csv")
dat <- fread("../processed_data/NAguideR_seqknn.csv") %>% as.data.frame()
genes <- anno_dat$`Gene name`[match(dat[,1],anno_dat$`Protein accession`)]
dat <- dat[,-1]
dat <- log2(dat)
dat <- data.frame(genes = genes,dat)

ssdat <- dat[!duplicated(dat$genes),]
ssmat <- ssdat[,-1] %>% as.matrix()
rownames(ssmat) <- ssdat$genes

ssdat <- gsva(ssmat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=T)


# --------------------------------- val ssGSEA --------------------------------

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

val_ssdat <- gsva(val_dat_mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=T)

# -------------------------------- t-test and boxplot -------------------------

label_fun <- function(ss,samp){
  pvalue_less <- c()
  pvalue_greater <- c()
  for(i in 1:nrow(ss)){
    tRes_less <- t.test(ss[i,][samp$Label == 0],
                        ss[i,][samp$Label == 1],alternative = "less")
    tRes_greater <- t.test(ss[i,][samp$Label == 0],
                           ss[i,][samp$Label == 1],alternative = "greater")
    pvalue_less <- c(pvalue_less,tRes_less$p.value)
    pvalue_greater <- c(pvalue_greater,tRes_greater$p.value)
  }
  
  label_ <- c()
  
  for(i in 1:length(pvalue_greater)){
    l = ""
    if(pvalue_less[i] < 0.05 | pvalue_greater[i] < 0.05){l = "*"}
    if(pvalue_less[i] < 0.01 | pvalue_greater[i] < 0.01){l = "**"}
    if(pvalue_less[i] < 0.001 | pvalue_greater[i] < 0.001){l = "***"}
    label_ <- c(label_,l)
  }
  return(label_)
}

label_fun(val_ssdat,val_sample)
label_fun(ssdat,sample)

ssplot_dat <- ssdat %>% t() %>% scale() %>% as.data.frame()
ssplot_dat$Group <- ifelse(sample$Label,"Responders","Non-responders")
ssplot_dat <- tidyr::pivot_longer(ssplot_dat, -Group,
                                  names_to = "Cell type", values_to = "Score")
ssplot_dat$Data <- "Our data"

val_ssplot_dat <- val_ssdat %>% t() %>% scale() %>% as.data.frame()
val_ssplot_dat$Group <- ifelse(val_sample$Label,"Responders","Non-responders")
val_ssplot_dat <- tidyr::pivot_longer(val_ssplot_dat, -Group,
                                      names_to = "Cell type", values_to = "Score")
val_ssplot_dat$Data <- "Validation data"

plot_dat <- rbind(ssplot_dat, val_ssplot_dat)

# plot
p_box <- ggplot(plot_dat,aes(x = `Cell type`, y = Score))+
  geom_boxplot(aes(fill = Group),position = position_dodge(0.5), width = 0.5, outlier.alpha = 0) +
  scale_fill_manual(values = c("#FF8C00", "#6E8B3D"))+
  coord_cartesian(ylim = c(-3,3) )+
  labs(x = "")+
  facet_grid(Data~.) +
  theme_bw()+
  theme( axis.text.x = element_text(angle = 45, hjust = 1 ),
         legend.position = "top"
  )

png(filename = "../results/tme_boxplot.png",
    width =20 ,height =16 ,units = "cm",res = 300)
p_box
dev.off()

