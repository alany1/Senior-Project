library(TwoSampleMR)
library(glue)
library(ggplot2)
#Chronotype exposure Data
chronotype_file <- "C:/Users/alany/Documents/Chronotype_mod.txt.xls"

#Read in chronotype exposure data
chronotype_exp_dat <- read_exposure_data(filename = chronotype_file, 
                                    sep='\t',
                                    snp_col = "LeadVariant",
                                    beta_col = "MetalnOR",
                                    se_col = "SE.corrected",
                                    effect_allele_col = "Allele1", #TODO: check if effect allele is 1 or 2?
                                    other_allele_col = "Allele2",
                                    eaf_col = "Freq1",
                                    pval_col = "MetaPCorrected")
                                 

#Name Exposure
chronotype_exp_dat$exposure <- "Chronotype"

#Print data head
print(head(chronotype_exp_dat))

#Clump data: ensure that the instruments (SNPs) are independent
chronotype_exp_dat <- clump_data(chronotype_exp_dat)

#Read outcome data
outcome_file <- "C:/Users/alany/Documents/TRICL_chronotype.txt.xls"
tmp <- read.table(outcome_file,header=1,sep="\t")

#Add beta value column
tmp$beta <- log(tmp$OR_fixed)
write.csv(tmp,"C:/Users/alany/Documents/TRICL_chronotype2.txt.csv", row.names = FALSE)

#Extract exposure SNPs from outcome data
outcome_dat <- read_outcome_data(snps = chronotype_exp_dat$SNP,
                                 filename = "TRICL_chronotype2.txt.csv",
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
dat <- harmonise_data(exposure_dat = chronotype_exp_dat,
                      outcome_dat = outcome_dat)

#Perform Mendelian Randomization
res <- mr(dat)

#Plot SNP effect on outcome vs exposure
#Plotting ONLY IVW
res <- mr(dat, method_list=c("mr_ivw"))
p1 <- mr_scatter_plot(res,dat)

#Save to WD
ggsave(p1[[1]],file="SNPchronotype.png",width=7,height=7)

#Forest plot
res_single <- mr_singlesnp(dat,all_method=c("mr_ivw"))
p2 <- mr_forest_plot(res_single)
ggsave(p2[[1]], file="ForestChronotype.png",width=12,height=12)

