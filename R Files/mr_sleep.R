library(TwoSampleMR)
library(ggplot2)
#Sleep Exposure Data
sleep_file <- "C:/Users/alany/Documents/sleep.txt.xls"

#Read in sleep exposure data
sleep_exp_dat <- read_exposure_data(filename = sleep_file, 
                                    sep='\t',
                                    snp_col = "SNP",
                                    beta_col = "Beta",
                                    se_col = "SE",
                                    effect_allele_col = "Effect",
                                    other_allele_col = "other",
                                    eaf_col = "EAF",
                                    pval_col = "P",
                                    gene_col = "Gene")

#Name Exposure
sleep_exp_dat$exposure <- "Sleep"

#Print data head
print(head(sleep_exp_dat))

#Clump data: ensure that the instruments (SNPs) are independent
sleep_exp_dat <- clump_data(sleep_exp_dat)

#Read outcome data
outcome_file <- "C:/Users/alany/Documents/TRICL_sleepf.txt.xls"
tmp <- read.table(outcome_file,header=1,sep="\t")

#Add beta value column
tmp$beta <- log(tmp$OR_fixed)
write.csv(tmp,"C:/Users/alany/Documents/TRICL_sleepf2.txt.csv", row.names = FALSE)

#Extract exposure SNPs from outcome data
outcome_dat <- read_outcome_data(snps = sleep_exp_dat$SNP,
                                 filename = "TRICL_sleepf2.txt.csv",
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
dat <- harmonise_data(exposure_dat = sleep_exp_dat,
                      outcome_dat = outcome_dat)

#Perform Mendelian Randomization
res <- mr(dat)
res

#Plot SNP effect on outcome vs exposure
#Plotting ONLY IVW
res <- mr(dat, method_list=c("mr_ivw"))
p1 <- mr_scatter_plot(res,dat)

#Save to WD
ggsave(p1[[1]],file="SNPsleep.png",width=7,height=7)

#Forest plot
res_single <- mr_singlesnp(dat,all_method=c("mr_ivw"))
p2 <- mr_forest_plot(res_single)
ggsave(p2[[1]], file="ForestSleep.png")

