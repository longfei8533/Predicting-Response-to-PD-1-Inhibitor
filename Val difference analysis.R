library(data.table)
library(tidyverse)
library(org.Hs.eg.db)
library(UniProt.ws)


# ------------------------------ Maxquant -----------------------------------

sam_dat <- fread("../data/IPX0004819000/proteinGroups.txt") %>% 
  dplyr::select(`Protein IDs`,starts_with("LFQ intensity"))

sam_dat <- sam_dat[-c(unique(c(grep("CON_",sam_dat$`Protein IDs`),
                                  grep("REV_",sam_dat$`Protein IDs`)
                                  ))),]
sam_dat = sam_dat[!rowSums(sam_dat[,-1]) ==0,]
fwrite(sam_dat,file = "../processed_data/validation_data.csv")


# ------------------------------------ NAguideR -------------------------------
# Use NAguideR to fill in the missing values
# Web url: http://www.omicsolution.org/wukong/NAguideR/
# output: ../processed_data/validation_data_seqknn.csv


# ------------------------------ Difference analysis ---------------------------

dat <- read.csv("../processed_data/validation_data_seqknn.csv") 
dat <- separate_longer_delim(dat,X,";")

dat_mat <- dat[,-1]
rownames(dat_mat) <- dat$X

prog_dat <- dplyr::select(dat_mat,ends_with("_N"))
resp_dat <- dplyr::select(dat_mat,ends_with("_"))

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
write.csv(res,file="../results/diff_val_result.csv",row.names = F)


# ------------------------- Change direction correlation -----------------------

res <- read.csv(file="../results/diff_val_result.csv")
my_res <-  fread(file="../results/diff_result.csv")

my_res_ <- my_res[my_res$pvalue < 0.05,]
res_ <- res[res$pvalue < 0.05,]

com_dat <- merge(my_res_,res_,by = "Protein")

com_dat$point_type <- 
  ifelse(log2(com_dat$ratio.x) * log2(com_dat$ratio.y) <0 , "Inconsistent","Consistent") %>% 
  as.factor()


cor.test(log2(com_dat$ratio.x), log2(com_dat$ratio.y))

p <- ggplot(data = com_dat,aes(x = log2(ratio.x), y = log2(ratio.y)))+
  geom_point(aes(color = point_type), shape = 19, size = 3 ) + 
  geom_smooth(method = "lm") + 
  annotate("text", x = -0.5 , y = 2,label = "Cor = 0.89",colour="black") + 
  scale_color_manual(values = c("Inconsistent" = "gray","Consistent" = "#1abc9c" )) +
  labs(x = "log2(Fold change) (Our data)", y = "log2(Fold change) (Validation data)" , color = "Direction of change") + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw()+
  theme(legend.position = "top")


png(filename = "../results/direction corr.png",width = 10,height = 10,units = "cm",res = 300)
print(p)
dev.off()
