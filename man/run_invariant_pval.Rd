% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{run_invariant_pval}
\alias{run_invariant_pval}
\title{run multinomial regression and calc the pvalues of invariant}
\usage{
run_invariant_pval(df_pheno, df_geno, df_pc10 = NA, X_term, G_term, Y_term)
}
\arguments{
\item{df_pheno}{the matrix that each row is a sample and each column is a feature.}

\item{df_geno}{the matrix of geno information, the additive model with each column is a sample, each row is a SNP.}

\item{df_pc10}{the matrix of first 10 pcs from geno, in the 10 columns.}

\item{G_term}{column names in df_pheno that used for all Geno parts GxE parts(interest items to be test)}

\item{Y_term}{column names in df_pheno that used for interested response.}

\item{X_tem}{column names in df_pheno that used for all covariates(non-interest items)}
}
\description{
this function  calculate the core pvalues of invariant. This function is wraped by run_invariant_pval
This function construct the input for function run_invariant_pval
}
\examples{
X_item=c("sex","age",
"PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
E_item=c("sun_hour_summer")
G_item=c("geno")
GE_item=c("sun_hour_summer:geno")
Y_item=c("response")
run_invariant_pval(df_pheno=cur_pheno,df_geno=geno,df_pc10=pc10_geno,
               X_term=c(X_item,E_item),G_term=c(G_item,GE_item),Y_term=Y_item)

}
