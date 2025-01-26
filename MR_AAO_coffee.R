install.packages("plyr")
install.packages("devtools")
devtools::install_github("rondolab/MR-PRESSO")

install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR") 

library(ggplot2)
library(dplyr)
library(MRPRESSO)
library(data.table)
library(googledrive)
library(TwoSampleMR)
 
file3 <- c("~/Downloads/MR_diplom/gwas_ready_UK.txt")

df1 <- fread(file3)
colnames(df1)

df1$samplesize.exposure <- 392083


file2_1 <- df1 %>% filter(df1[["pval.exposure"]] < 5e-8 & df1[["pval.exposure"]] != 0)


nrow(file2_1)
class(file2_1)

file2_1$samplesize.exposure <- 373522 
write.csv(df1,"~/Downloads/MR_diplom/gwas_ready_UK.txt",  row.names = FALSE)
write.table(df1, "~/Downloads/MR_diplom/gwas_ready_UK2.txt", sep = "\t", quote = FALSE, row.names = FALSE)


file0 <- c("~/Downloads/MR_diplom/gwas_ready_UK2.txt")
df1 <- fread(file0)
colnames(df1)
exp_data <- read_exposure_data(file0, sep="\t", snp_col = "SNP", beta_col = "beta.exposure",  eaf_col = "eaf.exposure", se_col="se.exposure", effect_allele_col = "effect_allele.exposure", other_allele_col= "other_allele.exposure", pval_col = "pval.exposure", samplesize_col = "samplesize.exposure", clump=TRUE)

exp_data$pval.exposure
nrow(exp_data)
# adding chr:pos manualy
chromosome =    c( "chr15", "chr7",  "chr7",  "chr7",  "chr18", "chr16", "chr2", 
"chr19", "chr2",  "chr12", "chr20", "chr15",
"chr1",  "chr6",  "chr5",   "chr1", "chr8",  "chr5",  "chr20",
"chr17", "chr17", "chr6",  "chr4",   "chr15", "chr16", "chr18", 
"chr3",  "chr22")
  
base_pair_location =   c( 75027880,17284577,75615006,73037956,  57808978, 53828066, 27730940, 
                          41353107, 637498,11316437, 62892739,75174251,
                          50576710,31840021,7381260,177889025,109156532,87949158, 45840459,
                          17845800, 46155786,51174232,2933031, 91437388, 70927078, 55080437, 
                          50429876,24844948)
  
  
chrpos <-  c(  "chr15:75027880", "chr7:17284577",  "chr7:75615006",  "chr7:73037956",  "chr18:57808978", "chr16:53828066", "chr2:27730940", 
               "chr19:41353107", "chr2:637498",  "chr12:11316437", "chr20:62892739", "chr15:75174251",
               "chr1:50576710",  "chr6:31840021",  "chr5:7381260",   "chr1:177889025", "chr8:109156532",  "chr5:87949158",  "chr20:45840459",
               "chr17:17845800", "chr17:46155786", "chr6:51174232",  "chr4:2933031",   "chr15:91437388", "chr16:70927078", "chr18:55080437", 
               "chr3:50429876",  "chr22:24844948")  
  

exp_data$chrpos <- chrpos
exp_data$chromosome <- chromosome
exp_data$base_pair_location <- base_pair_location
exp_data$SNP
nrow(exp_data)
Coffee_chrpos_asID <- format_data(exp_data, type = "exposure", snps = NULL,
                              header = TRUE, snp_col = "chrpos", beta_col = "beta.exposure",
                              se_col = "se.exposure",  eaf_col = "eaf.exposure", effect_allele_col = "effect_allele.exposure",
                              other_allele_col = "other_allele.exposure", pval_col = "pval.exposure", log_pval = FALSE,
                              samplesize_col = "samplesize.exposure", chr_col = "chromosome", pos_col = "base_pair_location")


write.csv(Coffee_chrpos_asID, "~/Downloads/MR_diplom/coffee_exposure_with_chr_position.csv", row.names = FALSE)

file2 <- "~/Downloads/IPDGC_AAO_GWAS_sumstats_april_2018_rsid.txt"
df1 <- fread(file2)
colnames(df1)
df1$rsid
out_data2 <- read_outcome_data(
  snps = exp_data$SNP,
  sep = "\t",
  filename = "~/Downloads/IPDGC_AAO_GWAS_sumstats_april_2018_rsid.txt",
  snp_col = "rsid",
  beta_col = "Effect",
  se_col = "FreqSE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "Freq1",
  pval_col = "P-value")
out_data2$SNP
out_data2$pval.outcome

nrow(out_data2)
out_data2$samplesize.outcome <- 34498

dat <- harmonise_data(exposure_dat = exp_data , outcome_dat = out_data2, action = 2)
res <- mr(dat)
dat$mr_keep.exposure
dat$mr_keep.outcome
dat$mr_keep

dat <- subset(dat, dat$mr_keep == TRUE)

# Adjust units if needed
dat$units.outcome <- "log odds"
dat$units.exposure <- "log odds"
str(dat)
dat2<-subset(dat, dat$eaf.exposure!="NA")
dat2$rsq.exposure<- get_r_from_pn(dat2$pval.exposure, dat2$samplesize.exposure)
dat2$rsq.outcome<- get_r_from_pn(dat2$pval.outcome, dat2$samplesize.outcome )


# Calculate Steiger's Z-test
dat2$effective_n.outcome <- dat2$samplesize.outcome
dat2$effective_n.exposure <- dat2$samplesize.exposure

st <- psych::r.test(n = dat2$effective_n.exposure, n2 = dat2$effective_n.outcome, 
                    r12 = sqrt(dat2$rsq.exposure), r34 = sqrt(dat2$rsq.outcome))
dat2$steiger_dir <- dat2$rsq.exposure > dat2$rsq.outcome
dat2$steiger_pval <- pnorm(-abs(st$z)) * 2
nrow(dat2)

# Perform Steiger filtering
steiger <- steiger_filtering(dat2)
nrow(steiger)
sig <- subset(steiger, steiger$steiger_dir == TRUE)
nrow(sig)
sig$SNP

# Perform MR PRESSO analysis
presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                    SdOutcome = "se.outcome", SdExposure = "se.exposure",
                    OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = sig,
                    NbDistribution = 1000, SignifThreshold = 0.05) 

filtered_dat <- sig[-c( 2  ,5,  6,  7,  8 , 9, 14, 15, 21, 25, 27), ]

sig <- filtered_dat

filtered_dat <- sig[-c( 8,14,16), ]

sig <- filtered_dat

# Perform MR PRESSO analysis again
presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                    SdOutcome = "se.outcome", SdExposure = "se.exposure",
                    OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = sig,
                    NbDistribution = 1000, SignifThreshold = 0.05) 

# Save MR PRESSO results to a file
capture.output(print(presso), file = "~/Downloads/Phenotype/presso_results2.txt")

# Save statistics to separate files
capture.output(print(R2), file = "~/Downloads/Phenotype/mean_R2.txt")
capture.output(print(F), file = "~/Downloads/Phenotype/F_statistic.txt")

sig$mr_keep

mr_res <- mr(sig)
mr_heterogeneity(sig)
mr_pleiotropy_test(sig)
directionality_test(sig)
mr_leaveoneout(sig)


mr_res_with_OR <- generate_odds_ratios(mr_res)
write.table(mr_res_with_OR, file=file.path("~/Downloads/Phenotype/", "mr_res_with_OR.txt"), row.names=FALSE, quote=FALSE, sep = "\t")
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
ggsave(p5[[1]], file = "~/Downloads/Phenotype/plot2.jpg", width = 7, height = 12)

