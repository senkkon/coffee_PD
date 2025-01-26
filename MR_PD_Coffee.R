library(TwoSampleMR)
library(ggplot2)
library(dplyr)
library(MRPRESSO)
library(data.table)

file2 <- c("~/Downloads/Caffeine_from_Coffee_GWAS_UKBiobank_MAFandINFOfiltered_Said_et_al 3")

df1 <- fread(file2)
colnames(df1)

file2_1 <- df1 %>% filter(df1[["pval.exposure"]] < 5e-08 & df1[["pval.exposure"]] != 0)
nrow(df1)
nrow(file2_1)


file2_1$samplesize.exposure <- 362316

write.csv(file2_1, "~/Downloads/MR_diplom/new_file_exposure.csv", row.names = FALSE)


exp_data <- read_exposure_data("~/Downloads/MR_diplom/gwas_ready_UK.txt", sep="\t", snp_col = "SNP", beta_col = "beta.exposure",  eaf_col = "eaf.exposure", se_col="se.exposure", effect_allele_col = "effect_allele.exposure", other_allele_col= "other_allele.exposure", pval_col = "pval.exposure", samplesize_col="samplesize.exposure", clump=TRUE)

nrow(out_data2)
out_data2 <- read_outcome_data(snps = exp_data$SNP, filename = "~/Downloads/MR_diplom/new.csv", sep=",", 
                               snp_col = "SNP", beta_col = "b", se_col = "se",
                               eaf_col = "Freq1", effect_allele_col = "Allele1", 
                               other_allele_col = "Allele2", pval_col = "P-value",
                               ncase_col = "ncase", ncontrol_col = "ncontrol")

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

# Create an empty steiger_pval column in dat1
#dat1$steiger_pval <- NA  # Or initialize it with appropriate values
#str(dat1)
#colnames(dat)


# Perform Steiger filtering
steiger <- steiger_filtering(dat1)
sig <- subset(steiger, steiger$steiger_dir == TRUE)


# Perform MR PRESSO analysis
presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                    SdOutcome = "se.outcome", SdExposure = "se.exposure",
                    OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = sig,
                    NbDistribution = 1000, SignifThreshold = 0.05) 

# Save MR PRESSO results to a file
capture.output(print(presso), file = "~/Downloads/Phenotype/PD/presso_results2.txt")

# Save statistics to separate files
capture.output(print(R2), file = "~/Downloads/Phenotype/PD/mean_R2.txt")
capture.output(print(F), file = "~/Downloads/Phenotype/PD/F_statistic.txt")

mr(sig)
mr_res_with_OR <- generate_odds_ratios(mr_res)
write.table(mr_res_with_OR, file=file.path("~/Downloads/Phenotype/PD", "mr_res_with_OR.txt"), row.names=FALSE, quote=FALSE, sep = "\t")
R2<-mean(sig$rsq.exposure)
capture.output(print(R2), file = paste( "r2.txt", sep = "/"))
n<-mean(sig$samplesize.exposure)
k <- nrow(subset(sig, sig$ambiguous == FALSE))
F<-(R2*(n-1-k))/((1-R2)*k)
capture.output(print(F), file = paste( "f.txt", sep = "/"))
mr_report(sig) 



# Generate MR single SNP forest plot
res_single <- mr_singlesnp(sig)
p5 <- mr_forest_plot(res_single)
p5[[1]]
ggsave(p5[[1]], file = "~/Downloads/Phenotype/PD/plot2.jpg", width = 7, height = 12)

