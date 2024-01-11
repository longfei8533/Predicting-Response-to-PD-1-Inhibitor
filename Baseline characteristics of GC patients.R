library(tidyverse)
library(tableone) 
library(knitr)


info_dat <- read.csv("../data/clinical_information.csv",
                     check.names = F)

a <- CreateTableOne(vars=c("Sex", "Age",
                           "Degree of differentiation", "Lauren's criteria", 
                           "T", "N", "M"), 
                    data =info_dat,
                    strata="Clinical outcomes", 
                    factorVars=c("Sex", 
                                 "Degree of differentiation", 
                                 "Lauren's criteria", 
                                 "T", "N", "M")) 

a_csv<- print(a, nonnormal = c("Age"),
              exact =c("Sex","Degree of differentiation", "Lauren's criteria", "T", "N", "M"),
              showAllLevels = TRUE,
              quote = FALSE, 
              noSpaces = TRUE, 
              printToggle = FALSE)

kable(a_csv,align = 'c')

write.csv(x = a_csv,file = "../results/baseline.csv")

