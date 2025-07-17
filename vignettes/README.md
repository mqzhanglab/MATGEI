# Introduction
Here we present a demo version of our real data analysis. This demo dataset contains 100 individuals from four kinds of skin cancers(basal cell carcinoma, squamous cell carcinoma,malignant melanoma,non-melanoma skin cancer) with each cancer subtype has 25 individuals.The four kinds of cancer subtypes are further collapsed into a binary factor carcinoma(TRUE for basal cell carcinoma and squamous cell carcinoma) for the comparison of the logistic regression.

This dataset also contains the sex, age of diagnose, and two environmental factor: uv_protect(binary) and sun_hour_summer(continuous), and the first 10 PCs from genetics.

This dataset contains 10 variants.

# Initial Setup

R packages required in this pipeline:
```{r load_libraries,eval=runcode}
library(BEDMatrix)
library(nnet)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr) 
library(scales)
library(dplyr)
library(SKAT)
library(ACAT)
library(MATGEI)

```

# Computation



## Basic Settings

Senario settings in our analysis are as follows: we set the column UV_protect as our environmental factors.


```{r basic_settings,eval=runcode}
#we assign column names below.
env_factor="uv_protect" #column name of df_pheno. The column for environmental factors.
response="cancer"       #column name of df_pheno. Response level for categorical regresion.
response_g="carcinoma"  #column name of df_pheno. Response level for logistic regression.

X_item=c("sex","age",
         "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
E_item=c(env_factor)
G_item=c("geno")
GE_item=c(paste0(env_factor,":geno"))
Y_item=c(response)
Yg_item=c(response_g)


```

## Data Preparation
```{R data_load,eval=runcode}
###########input###############

df_pheno=readRDS("../data/df_pheno.rds")
df_geno=readRDS("../data/df_geno.rds")
df_pc10=readRDS("../data/df_pc10.rds")

```

After read in the data, we run MATGEI for testing.

```{r run_MATGEI,eval=runcode,message = FALSE,warning = FALSE,include = FALSE,results = FALSE}

#for G test
res_G=MATGEI::run_invariant_pval(df_pheno=df_pheno,df_geno=df_geno,df_pc10=df_pc10,
                         X_term=c(X_item,E_item),G_term=c(G_item),Y_term=Y_item)
                   

#for GE test
res_GE=MATGEI::run_invariant_pval(df_pheno=df_pheno,df_geno=df_geno,df_pc10=df_pc10,
                          X_term=c(X_item,E_item,G_item),G_term=c(GE_item),Y_term=Y_item)

#for G+GE joint test
res_GGE=MATGEI::run_invariant_pval(df_pheno=df_pheno,df_geno=df_geno,df_pc10=df_pc10,
                           X_term=c(X_item,E_item),G_term=c(G_item,GE_item),Y_term=Y_item)



```


Then we run logistic regression for comparison.
```{r run_logistic,eval=FALSE}

#for G test
glm_G=MATGEI::run_glm_pval(df_pheno=df_pheno,df_geno=df_geno,df_pc10=df_pc10,
                   X_term=c(X_item,E_item),G_term=c(G_item),Y_term=Yg_item)

#for GE test
glm_GE=MATGEI::run_glm_pval(df_pheno=df_pheno,df_geno=df_geno,df_pc10=df_pc10,
                    X_term=c(X_item,E_item,G_item),G_term=c(GE_item),Y_term=Yg_item)

#for G+GE joint test
glm_GGE=MATGEI::run_glm_pval(df_pheno=df_pheno,df_geno=df_geno,df_pc10=df_pc10,
                     X_term=c(X_item,E_item),G_term=c(G_item,GE_item),Y_term=Yg_item)

```


## Extract Results

We can get our results from the result list. For example, the G test of the 1st variant.

```{r sig_genes,eval=runcode}

res_G[[1]]

```


# Session Information
```{r ,eval=runcode}
sessionInfo()
```


# Reference
