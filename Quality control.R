
library(tidyverse)
library(readr)
library(ggrepel)
library(ggstatsplot)
library(ggpubr)
library(ggplot2)
library(corrplot)


dat <- read.csv(file = "../processed_data/protein_intensity.csv")
dat_mat <- dat[,-1] %>% as.matrix()
rownames(dat_mat) <- dat$Uniprot.id

sample_info <- read.csv("../processed_data/sample_data.csv")


#-------------------------------可定量蛋白-------------------------------------
his_dat_ <- apply(dat_mat,2,FUN = function(x){sum(x>0)}) 
his_dat <- data.frame("sample" = factor(sample_info$Sample,
                                        levels = sample_info$Sample, 
                                        labels = substr(sample_info$Sample,8,9)),
                      "number" = his_dat_
                      )
png(filename = "../results/protein number.png",
    width =12 ,height =8 ,units = "cm",res = 300)
ggplot(data = his_dat,aes(x = sample, y = number)) +
  geom_col(fill = c("#36648B"))+
  theme_bw()+
  labs(x = "Sample",
       y = "The number of proteins"
       )

dev.off()

#---------------------------- Protein location ---------------------------------


## Method
# This file contains the large-scale whole human proteome predictions by Hum-mPLoc 3.0 (total 20197 proteins). 
# 
# Format:
# AC
# ID	
# SEQ
# Predicted Subcellular location by Hum-mPLoc 3.0	
# Outputted prediction scores for each location from Hum-mPLoc 3.0: Centrosome,Cytoplasm,Cytoskeleton,Endoplasmic reticulum,Endosome,Extracellular,Golgi apparatus,Lysosome,Mitochondrion,Nucleus,Peroxisome,Plasma membrane
# //
# 
# Online search is also supported at: 
# http://www.csbio.sjtu.edu.cn/bioinf/Hum-mPLoc3/WPP.html
# 
# Contact: hbshen@sjtu.edu.cn, zhouhang@sjtu.edu.cn



p_dat <- readLines("../data/protein subcellular localization prediction.txt")

u <- which(p_dat=="//")

uid <- u[-length(u)]
uniprot_id <- p_dat[uid+1]

ul <- u[-1]
loc <- p_dat[ul-2]

loc_name <- lapply(loc,function(x){
  strsplit(x,split=".  ") %>% unlist()
}) %>% unlist() %>% unique()


ll <- lapply(loc,function(x){
  l <- strsplit(x,split=".  ") %>% unlist()
  v <- c(rep(0,12))
  v[which(loc_name %in% l)] <- 1
  return(v)
})

loc_mat <- data.frame(ll) %>% t()
colnames(loc_mat) <- loc_name
rownames(loc_mat) <- uniprot_id

loc_dat <- colSums(loc_mat)
loc_dat <- data.frame(loc = factor(names(loc_dat),
                                   levels = names(loc_dat)[order(loc_dat,decreasing = T)]),
                      number = loc_dat
)
loc_dat$Group <- "Reference"



my_loc_mat <- loc_mat[uniprot_id %in% rownames(dat_mat),]
my_loc_dat <- colSums(my_loc_mat)
my_loc_dat <- data.frame(loc = factor(names(my_loc_dat),
                                       levels = names(my_loc_dat)[order(my_loc_dat,decreasing = T)]),
                          number = my_loc_dat
)
my_loc_dat$Group <- "Our"
my_loc_dat$number <- my_loc_dat$number*2


all_dat <- rbind(loc_dat,my_loc_dat)
all_dat$Group <- factor(all_dat$Group,levels = c("Reference","Our"))

p <- ggplot(data = all_dat,aes(x = loc, y = number)) +
  geom_col(aes(fill = Group),position = "dodge")+
  labs(x = "Subcellular location",
       y = "Reference number"
  )+
  scale_fill_manual(values = c("Reference" = "#36648B", "Our"="#458B00"))+
  theme_bw()+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1)
  )
  

p <- p + scale_y_continuous(
  sec.axis = sec_axis(~./2, 
                      name = "Our number"))

png(filename = "../results/protein loc.png",
    width =12 ,height =10 ,units = "cm",res = 300)
print(p)
dev.off()

#------------------------- Pairwise correlation between samples ----------------

cor_dat <- dat_mat

M <- cor(cor_dat)

corr_value <- (sum(M)-28)/(28*28-28)
print(corr_value)

png(filename = "../results/pairwise corr.png",
    width =12 ,height =10 ,units = "cm",res = 300)
corrplot(M,type = "full", method = "square",
         tl.pos = "l",
         col = rev(COL2('RdBu', 200)),
         tl.col = "black" , tl.cex = 0.5
)
dev.off()

# -------------------------------  Distribution of proteins intensity  -------------------------------------
pp_dat <- tidyr::gather(as.data.frame(log2(dat_mat+1)),key="cases",value="value")
pp_dat$cases <- factor(sample_info$Sample,levels = sample_info$Sample, labels = substr(sample_info$Sample,8,9))
pp_density <- ggplot(pp_dat,aes(x=value))+
  geom_density(aes(color=cases),show.legend = F)+
  labs(x="",y="Density")+
  theme_bw()+
  theme(
    legend.position = "none"
  )
pp_violin <- ggplot (pp_dat,aes(x=cases,y=value))+
  geom_violin(aes(fill=cases),show.legend = F)+
  stat_summary(fun = "median", fill = "white", size = 2, geom = "point", shape = 23) +
  coord_cartesian(ylim = c(5,30))+
  labs(y=substitute(paste(log[2],"(Protein intensity)")),x="Sample")+
  theme_bw()
png(filename = "../results/violin_plot.png",
    width =12 ,height =8 ,units = "cm",res = 300)
pp_violin
dev.off()


