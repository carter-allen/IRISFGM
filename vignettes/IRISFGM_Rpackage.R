## ----setwd, eval =TRUE, echo = TRUE-------------------------------------------
# dir.create("your working directory",recursive = TRUE)
# setwd("your working directory")
library(IRISFGM)

## ----txt, eval= TRUE, echo = TRUE---------------------------------------------
InputMatrix <- read.table(url("https://bmbl.bmi.osumc.edu/downloadFiles/Yan_expression.txt"),
                          header = TRUE, 
                          row.names = 1,
                          check.names = FALSE)

## ----create_object, eval= TRUE, echo = TRUE,message=TRUE----------------------
set.seed(123)
seed_idx <- sample(1:nrow(InputMatrix),3000)
InputMatrix_sub <- InputMatrix[seed_idx,]
object <- CreateIRISFGMObject(InputMatrix_sub)

## ----add_metadata, eval= TRUE, echo = TRUE------------------------------------
my_meta <- read.table(url("https://bmbl.bmi.osumc.edu/downloadFiles/Yan_cell_label.txt"),header = TRUE,row.names = 1)
object <- AddMeta(object, meta.info = my_meta)

## ----plot_metadata,eval= TRUE, echo = TRUE------------------------------------
PlotMeta(object)

## ----subset_data,eval= TRUE, echo =  TRUE-------------------------------------
object <- SubsetData(object , nFeature.upper=2000,nFeature.lower=250)

