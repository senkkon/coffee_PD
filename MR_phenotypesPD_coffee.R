# Load required libraries with necessary checks
library(TwoSampleMR)
library(ggplot2)
library(dplyr)
library(MRPRESSO)
library(data.table)

file0 <- c("~/Downloads/MR_diplom/gwas_ready_UK.txt")
df1 <- fread(file0)
colnames(df1)
df1$pval.exposure


file2 <- c("~/Downloads/MR_diplom/new_full_coffee_exposure.csv")

exp_data <- read_exposure_data(file0, sep='\t', snp_col = "SNP", beta_col = "beta.exposure",  eaf_col = "eaf.exposure", se_col="se.exposure", effect_allele_col = "effect_allele.exposure", other_allele_col= "other_allele.exposure", pval_col = "pval.exposure", samplesize_col = "samplesize.exposure", clump=TRUE)

nrow(exp_data)
out_data2 <- read_outcome_data(snps = exp_data$SNP, filename = "~/Downloads/Phenotype/Data_phenotype/base_SLEEP.txt_merged_with_ref.txt", sep="\t", 
                               snp_col = "RSID", beta_col = "BETA", se_col = "SE",
                               eaf_col = "MAF", effect_allele_col = "REF", 
                               other_allele_col = "ALT", pval_col = "P",
                               samplesize_col="N")

dat <- harmonise_data(exposure_dat = exp_data, outcome_dat = out_data2, action = 2)
res <- mr(dat)



# Adjust units if needed
dat$units.outcome <- "log odds"
dat$units.exposure <- "log odds"
dat1<-subset(dat, dat$eaf.exposure!="NA")
dat1$rsq.exposure<- get_r_from_pn(dat1$pval.exposure, dat1$samplesize.exposure)
dat1$rsq.outcome<- get_r_from_pn(dat1$pval.outcome, dat1$samplesize.outcome)
colnames(dat1)

#p_val_steiger <- cor.test(dat2$rsq.exposure, dat2$rsq.outcome, method = "pearson")$p.value
dat1$effective_n.outcome <- dat1$samplesize.outcome
dat1$effective_n.exposure <- dat1$samplesize.exposure

st <- psych::r.test(n = dat1$effective_n.exposure, n2 = dat1$effective_n.outcome, 
                    r12 = sqrt(dat1$rsq.exposure), r34 = sqrt(dat1$rsq.outcome))
dat1$steiger_dir <- dat1$rsq.exposure > dat1$rsq.outcome
dat1$steiger_pval <- pnorm(-abs(st$z)) * 2


# Perform Steiger filtering
steiger <- steiger_filtering(dat1)
sig <- subset(steiger, steiger$steiger_dir == TRUE)

mr_heterogeneity(sig)
mr_pleiotropy_test(sig)
directionality_test(sig)
# Perform MR PRESSO analysis
presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                    SdOutcome = "se.outcome", SdExposure = "se.exposure",
                    OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat,
                    NbDistribution = 1000, SignifThreshold = 0.05) 

# Save MR PRESSO results to a file
capture.output(print(presso), file = "~/Downloads/Phenotype/presso_results2.txt")

# Save statistics to separate files
capture.output(print(R2), file = "~/Downloads/Phenotype/mean_R2.txt")
capture.output(print(F), file = "~/Downloads/Phenotype/F_statistic.txt")

mr(sig)
R2<-mean(sig$rsq.exposure)
capture.output(print(R2), file = paste( "r2.txt", sep = "/"))
n<-mean(sig$samplesize.exposure)
k <- nrow(subset(sig, sig$ambiguous == FALSE))
F<-(R2*(n-1-k))/((1-R2)*k)
capture.output(print(F), file = paste( "f.txt", sep = "/"))
mr_report(sig) 

capture.output(print(F), file = "~/Downloads/Phenotype/mr_res_with_OR.txt")
mr_res <- mr(sig)
mr_res_with_OR <- generate_odds_ratios(mr_res)
write.table(mr_res_with_OR, file="~/Downloads/Phenotype/mr_res_with_OR.txt", row.names=FALSE, quote=FALSE, sep = "\t")



# Generate MR single SNP forest plot
res_single <- mr_singlesnp(sig)
p5 <- mr_forest_plot(res_single)
p5[[1]]
ggsave(p5[[1]], file = "~/Downloads/plot2.jpg", width = 7, height = 12)


