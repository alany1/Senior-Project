library(TwoSampleMR)
library(ggplot2)
#insomnia Exposure Data
insomnia_file <- "C:/Users/alany/Documents/insomnia.txt.xls"

#Read outcome data
outcome_file <- "C:/Users/alany/Documents/TRICL_insomnia.txt.xls"
tmp <- read.table(outcome_file,header=1,sep="\t")

#Add beta value column
tmp$beta <- log(tmp$OR_fixed)
write.csv(tmp,"C:/Users/alany/Documents/TRICL_insomnia2.txt.csv",row.names = FALSE)

for (i in 1:4) {
  #Read in insomnia exposure data
  insomnia_exp_dat <- read_exposure_data(filename = insomnia_file, 
                                         sep='\t',
                                         snp_col = "SNP",
                                         beta_col = paste("Beta",i,sep=""),
                                         se_col = paste("SE",i,sep=""),
                                         effect_allele_col = "Eff",
                                         other_allele_col = "Alt",
                                         eaf_col = "EAF",
                                         pval_col = paste("p.val",i,sep=""),
                                         gene_col = "NearestGene")
  #Name Exposure
  insomnia_exp_dat$exposure <- "Insomnia"
  
  #Clump data: ensure that the instruments (SNPs) are independent
  insomnia_exp_dat <- clump_data(insomnia_exp_dat)
  
  #Extract exposure SNPs from outcome data
  outcome_dat <- read_outcome_data(snps = insomnia_exp_dat$SNP,
                                   filename = "TRICL_insomnia2.txt.csv",
                                   sep = ",",
                                   snp_col = "rs_number",
                                   beta_col = "beta",
                                   se_col = "StdError_fixed",
                                   effect_allele_col = "effect_allele",
                                   other_allele_col = "reference_allele",
                                   eaf_col = "EAF",
                                   pval_col = "Pvalue_fixed")
  #Rename outcome to lung cancer
  outcome_dat$outcome <- "Lung Cancer"
  #Harmonize exposure and outcome data
  dat <- harmonise_data(exposure_dat = insomnia_exp_dat,
                        outcome_dat = outcome_dat)
  #Perform Mendelian Randomization
  res <- mr(dat,method_list=c("mr_ivw"))
  p1<-mr_scatter_plot(res,dat)
  ggsave(p1[[1]],file=paste("SNPInsomnia",i,".png",sep=""))
  res_single<- mr_singlesnp(dat,all_method=c("mr_ivw"))
  p2<-mr_forest_plot(res_single)
  ggsave(p2[[1]],file=paste("ForestInsomnia",i,".png",sep=""))
}
  


