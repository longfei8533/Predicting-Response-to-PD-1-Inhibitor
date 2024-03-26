library(e1071)
library(tidyverse)
library(ggrepel)
library(ggstatsplot)
library(ggpubr)
library(data.table)
library(caret)
library(pROC)
library(Boruta)

# load data

my_dat <- read.csv(file = "../processed_data/my_dat_ml_diff.csv") %>% as.matrix()
val_dat <- read.csv(file = "../processed_data/val_dat_ml_diff.csv") %>% as.matrix()

my_sample <- read.csv(file = "../processed_data/sample_data.csv")
val_sample <- read.csv(file = "../processed_data/validation_sample_data.csv")

# metric

metric_fun <- function(pred,true){
  confusion_matrix <- confusionMatrix(pred, true)
  acc <- confusion_matrix$overall["Accuracy"]
  pre <- confusion_matrix$byClass["Precision"]
  rec <- confusion_matrix$byClass["Recall"]
  f1 <- confusion_matrix$byClass["F1"]
  auroc <- roc(true, attr(pred,"probabilities")[,1])
}

## ----------------------------- select features ------------------------------

set.seed(2023)
boruta <- Boruta(x = my_dat,y = factor(my_sample$Label), doTrace = 2, maxRuns = 500)

imp_dat <- boruta$ImpHistory %>% as.data.frame() %>% 
  tidyr::pivot_longer(everything(),names_to = "key",values_to = "value",cols_vary = "fastest")
imp_dat$type <- rep(c(as.character(boruta$finalDecision),"shadowMax","shadowMean","shadowMin"),
                    nrow(boruta$ImpHistory))

imp_dat <- imp_dat[is.finite(imp_dat$value),]

sort_by <- group_by(imp_dat,key) %>%
  summarise(median = median(value)) %>% 
  arrange(median)

imp_dat$key <- factor(imp_dat$key,levels = sort_by$key)


png(file = "../results/attribute importance.png",width = 24 ,height =16 ,units = "cm",res = 300)
ggplot(data = imp_dat,aes(x = key,y = value))+
  geom_boxplot(aes(fill = type)) +
  labs(x = "Attributes", y = "Importance", fill = "Type") +
  theme_bw() +
  theme( axis.text.x = element_text(angle = 90 ),
         legend.position = "top"
  )
dev.off()

sel_my_dat <- my_dat[,boruta$finalDecision == "Confirmed"]
sel_val_dat <- val_dat[,boruta$finalDecision == "Confirmed"]

## ---------------------------- Building model ---------------------------------

for(i in c("linear","polynomia","sigmoid","radial")){
  mymodel <- svm(x = sel_my_dat,y = as.character(my_sample$Label),type = "C",kernel = i)
  pred <- predict(mymodel,sel_my_dat)
  met <- metric_fun(pred = pred,factor(my_sample$Label))
  print(i)
  print(met)
}

# linear model
mymodel <- svm(x = sel_my_dat,y = as.character(my_sample$Label),
               type = "C",kernel = "linear",probability = TRUE)

# support vectors
support_vectors <- mymodel$SV

# coefficients
coefficients <- t(mymodel$coefs) %*% mymodel$SV

# intercept
intercept <- -mymodel$rho

write.table(support_vectors,file = "../results/support_vectors.csv",sep = "\t",row.names = F)


# --------------------------------- performance -------------------------------

pred1 <- predict(mymodel,sel_my_dat,probability = TRUE)
pred2 <- predict(mymodel,sel_val_dat,probability = TRUE)
met1 <- metric_fun(pred = pred1,factor(my_sample$Label))
met2 <- metric_fun(pred = pred2,factor(val_sample$Label))
print(met1)
print(met2)

# ROC-AUC

plot_dat <- data.frame("Score" = attr(pred1,"probabilities")[,1], 
                       "Label" = my_sample$Label,
                       check.names = F)
plot_dat$Data <- "Our data"
val_plot_dat <- data.frame("Score" = attr(pred2,"probabilities")[,1], 
                           "Label" = val_sample$Label,
                           check.names = F)
val_plot_dat$Data <- "Validation data"

roc_res <- roc(plot_dat$Label, plot_dat$Score)
val_roc_res <- roc(val_plot_dat$Label, val_plot_dat$Score)

png(filename = "../results/svm_roc.png",
    width =12 ,height =12 ,units = "cm",res = 300)
plot(roc_res, main = "Machine learning",
     col = "#3CB371", 
     lwd = 3, 
     print.auc=TRUE,
     print.auc.y=0.4
     )
plot(val_roc_res,add=TRUE,
     col="#9A32CD",
     lwd = 3,
     print.auc = TRUE,
     print.auc.y=0.5
     )

legend("bottomright", legend=c("Training","Validation"),
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

png(filename = "../results/svm_box.png",
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









