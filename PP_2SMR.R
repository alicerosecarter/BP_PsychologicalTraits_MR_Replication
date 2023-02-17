################################################################################
##          Performing 2-sample MR using the TwoSampleMR R package            ##
################################################################################

setwd("")

## Install and load packages used in script
################################################################################

#install.packages("rlang")
library(rlang)

# Install packages
#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")
#remotes::install_github("MRCIEU/ieugwasr")
#remotes::install_github("MRCIEU/MRInstruments")
#remotes::install_github("n-mounier/MRlap")
#install.packages("data.table")
#install.packages("devtools")
#install.packages("dplyr")
#install.packages("xlsx")
#install.packages("readr)
#install.packages("stringr")
#install.packages("tidyverse")

# Load Packages
library(TwoSampleMR)
library(ieugwasr)
library(MendelianRandomization)
library(MRInstruments)
library(MRlap)
library(dplyr)
library(xlsx)
library(readr)
library(stringr)

ieugwasr::check_access_token()

################################################################################
                      # Pule pressure exposure #
################################################################################
# Read in pulse pressure GWAS (GWAS not available in IEU open GWAS - summ stats downloaded from EMBL)

pp_GWAS <- read.delim("Data/PP/Evangelou_30224653_PP.txt.gz", sep=" ")

# Split marker ID variable into chormosome, position and type (i.e., SNP)
# this will create a 3 columns matrix, for before, middle and after :
split_chr_pos <- str_split_fixed(pp_GWAS$b,fixed(":"),3)

# Merge chromosome, position and SNP with summary stats
pp_GWAS_pos <- (cbind(pp_GWAS, split_chr_pos))
 
# Remove columns with completely missing data
pp_GWAS_pos = subset(pp_GWAS_pos, select = -c(TotalSampleSize,N_effective) )

# Rename variables to be useful (and match up with the columns, since when reading in these get jumbled)
pp_GWAS <- rename(pp_GWAS_pos, c(
  "MarkerName"="b",
  "Allele1"="a",
  "Allele2"="pMarkerName",
  "Freq1"="Allele1",
  "Effect"="Allele2",
  "StdErr"="Freq1",
  "P"="Effect",
  "TotalSampleSize"="StdErr",
  "N_effective"="P",
  "Chr"="1",
  "position"="2",
  "type"="3"))

# Open summary stats for SBP GWAS - same consortia as PP - from IEU open GWAS to merge on Marker IDs and identify RSIDs in PP GWAS

# Need to double check the reference panels used are the same (should be since from same original GWAS)

gwasinfo("ieu-b-38")
sbp_GWAS <- tophits(id="ieu-b-38", clump=0, pval =1)

# Merge two GWAS together based on Marker IDs, keeping only RSIDs from the SBP GWAS
merged_GWAS <- merge(x = pp_GWAS, y =sbp_GWAS [ , c("position", "rsid")], by = "position", all.x=TRUE)

# CHECK THIS DOESN'T ALSO NEED TO MERGE ON CHROMOSOME AS WELL AS POSITION

merged_pp_gwas <- rename(merged_GWAS, c(
  "SNP"="rsid",
  "effect_allele.exposure"="Allele1",
  "other_allele.exposure"="Allele2",
  "beta.exposure"="Effect",
  "se.exposure"="StdErr",
  "eaf.exposure"="Freq1",
  "pval.exposure"="P",
  "N"="N_effective"))

merged_pp_gwas$exposure = "PP"

# Subset to only keep genome-wide significant SNPs
merged_pp_sig <- subset(merged_pp_gwas, pval.exposure<5e-08)

# Subset data to remove any missing SNPs (used in reverse direction analyses)
merged_pp_matchonly <- subset(merged_pp_gwas, SNP!="NA")

# Clump SNPs (make sure to ONLY do this after restricting to significant SNPs)
# Check that merging is adequate and all SNPs for analysis contain an RSID
pp_liberal <- clump_data(merged_pp_sig, clump_kb = 10000, clump_r2 = 0.05,)
pp_stringent <- clump_data(merged_pp_sig, clump_kb = 10000, clump_r2 = 0.001,)

pp_liberal$id.exposure = "PP - liberal clump"
pp_stringent$id.exposure = "PP - stringent/standard clump"

# Calculate R2 and F stat for exposure data
# Liberal PP F stat
pp_liberal$r2 <- (2 * (pp_liberal$beta.exposure^2) * pp_liberal$eaf.exposure * (1 - pp_liberal$eaf.exposure)) /
  (2 * (pp_liberal$beta.exposure^2) * pp_liberal$eaf.exposure * (1 - pp_liberal$eaf.exposure) +
     2 * pp_liberal$N * pp_liberal$eaf.exposure * (1 - pp_liberal$eaf.exposure) * pp_liberal$se.exposure^2)

pp_liberal$F <- pp_liberal$r2 * (pp_liberal$N - 2) / (1 - pp_liberal$r2)
pp_liberal_meanF <- mean(pp_liberal$F)

# Stringent pp F stat
pp_stringent$r2 <- (2 * (pp_stringent$beta.exposure^2) * pp_stringent$eaf.exposure * (1 - pp_stringent$eaf.exposure)) /
  (2 * (pp_stringent$beta.exposure^2) * pp_stringent$eaf.exposure * (1 - pp_stringent$eaf.exposure) +
     2 * pp_stringent$N * pp_stringent$eaf.exposure * (1 - pp_stringent$eaf.exposure) * pp_stringent$se.exposure^2)

pp_stringent$F <- pp_stringent$r2 * (pp_stringent$N - 2) / (1 - pp_stringent$r2)
pp_stringent_meanF <- mean(pp_stringent$F)

# Find SNPs in outcome data
ao <- available_outcomes()

anxiety_pp_liberal <- extract_outcome_data(snps = pp_liberal$SNP, outcomes = "ukb-b-11311")
anxiety_pp_stringent <- extract_outcome_data(snps = pp_stringent$SNP, outcomes = "ukb-b-11311")

swb_pp_liberal <- extract_outcome_data(snps = pp_liberal$SNP, outcomes = "ieu-a-1009")
swb_pp_stringent <- extract_outcome_data(snps = pp_stringent$SNP, outcomes = "ieu-a-1009")

neuroticism_pp_liberal <- extract_outcome_data(snps = pp_liberal$SNP, outcomes = "ieu-a-1007")
neuroticism_pp_stringent <- extract_outcome_data(snps = pp_stringent$SNP, outcomes = "ieu-a-1007")

depressive_symp_pp_liberal <- extract_outcome_data(snps = pp_liberal$SNP, outcomes = "ieu-a-1000")
depressive_symp_pp_stringent <- extract_outcome_data(snps = pp_stringent$SNP, outcomes = "ieu-a-1000")

# Harmonise datasets
# Note a large number of "incompatible alleles" are removed - possibly just because the data hasn't gone through the same checking as those on openGWAS
harm_anxiety_pp_liberal <- harmonise_data(exposure_dat = pp_liberal, outcome_dat = anxiety_pp_liberal)
harm_anxiety_pp_stringent <- harmonise_data(exposure_dat = pp_stringent, outcome_dat = anxiety_pp_stringent)

harm_swb_pp_liberal <- harmonise_data(exposure_dat = pp_liberal, outcome_dat = swb_pp_liberal)
harm_swb_pp_stringent <- harmonise_data(exposure_dat = pp_stringent, outcome_dat = swb_pp_stringent)

harm_neuroticism_pp_liberal <- harmonise_data(exposure_dat = pp_liberal, outcome_dat = neuroticism_pp_liberal)
harm_neuroticism_pp_stringent <- harmonise_data(exposure_dat = pp_stringent, outcome_dat = neuroticism_pp_stringent)

harm_depressive_symp_pp_liberal <- harmonise_data(exposure_dat = pp_liberal, outcome_dat = depressive_symp_pp_liberal)
harm_depressive_symp_pp_stringent <- harmonise_data(exposure_dat = pp_stringent, outcome_dat = depressive_symp_pp_stringent)


# F stats for each analysis
anxiety_pp_liberal_meanF <- mean(harm_anxiety_pp_liberal$F)
anxiety_pp_stringent_meanF <- mean(harm_anxiety_pp_stringent$F)

depressive_symp_pp_liberal_meanF <- mean(harm_depressive_symp_pp_liberal$F)
depressive_symp_pp_stringent_meanF <- mean(harm_depressive_symp_pp_stringent$F)

neuroticism_pp_liberal_meanF <- mean(harm_neuroticism_pp_liberal$F)
neuroticism_pp_stringent_meanF <- mean(harm_neuroticism_pp_stringent$F)

swb_pp_liberal_meanF <- mean(harm_swb_pp_liberal$F)
swb_pp_stringent_meanF <- mean(harm_swb_pp_stringent$F)

# Run MR
anxiety_pp_liberal_results <- mr(dat = harm_anxiety_pp_liberal)
anxiety_pp_stringent_results <- mr(dat = harm_anxiety_pp_stringent)

swb_pp_liberal_results <- mr(dat = harm_swb_pp_liberal)
swb_pp_stringent_results <- mr(dat = harm_swb_pp_stringent)

neuroticism_pp_liberal_results <- mr(dat = harm_neuroticism_pp_liberal)
neuroticism_pp_stringent_results <- mr(dat = harm_neuroticism_pp_stringent)

depressive_symp_pp_liberal_results <- mr(dat = harm_depressive_symp_pp_liberal)
depressive_symp_pp_stringent_results <- mr(dat = harm_depressive_symp_pp_stringent)

# Send results to excel

results_summary <- rbind(anxiety_pp_liberal_results , 
                         anxiety_pp_stringent_results,
                         swb_pp_liberal_results,
                         swb_pp_stringent_results,
                         neuroticism_pp_liberal_results,
                         neuroticism_pp_stringent_results, 
                         depressive_symp_pp_liberal_results,
                         depressive_symp_pp_stringent_results)

results_summary <- results_summary %>%
  mutate(id.outcome = recode(id.outcome, 'ukb-b-11311' = 'Anxiety', 'ieu-a-1009' = 'SWB', 'ieu-a-1007' = 'Neuroticism', 'ieu-a-1000' ='Depressive symptoms' ))

write.xlsx(results_summary, file = "Results/PP_exposure.xlsx")

################################################################################
                        # Psychiatric traits exposure #
################################################################################

# Format pulse pressure outcome data
merged_pp_outcome <- rename(merged_pp_gwas, c(
  "effect_allele.outcome"="effect_allele.exposure",
  "other_allele.outcome"="other_allele.exposure",
  "beta.outcome"="beta.exposure",
  "se.outcome"="se.exposure",
  "eaf.outcome"="eaf.exposure",
  "pval.outcome"="pval.exposure",))
merged_pp_outcome$outcome = "Pulse Pressure"
merged_pp_outcome$id.outcome = "Pulse Pressure"

available_outcomes(access_token = ieugwasr::check_access_token())

# Anxiety exposure
# No SNPs available using default parameters. By relaxing the P value threshold I get the error "Warning message:Unknown or uninitialised column: `trait`."
# Therefore extract SNPs manually from the Open GWAS catalogue using a relaxed P value threshold, but standard clumping thresholds
#anxiety_exposure <- extract_instruments(outcomes = "ukb-b-11311",  p1 = 5e-06)

anxiety_tophits <- tophits(id="ukb-b-11311", pval = 5e-06)

anxiety_exposure <- rename(anxiety_tophits, c(
  "SNP"="rsid",
  "effect_allele.exposure"="ea",
  "other_allele.exposure"="nea",
  "beta.exposure"="beta",
  "se.exposure"="se",
  "eaf.exposure"="eaf",
  "pval.exposure"="p",
  "N"="n"))

anxiety_exposure$exposure = "anxiety"
anxiety_exposure$id.exposure = "anxiety"

harmonised_data_anxiety_pp <- harmonise_data(exposure_dat = anxiety_exposure, outcome_dat = merged_pp_matchonly, action=2)
anxiety_pp_results <- mr(dat = harmonised_data_anxiety_pp)

# Subjective wellbeing exposure - 
# Single SNP for SWB, which isn't available in the PP GWAS 
swb_exposure <- extract_instruments(outcomes = "ieu-a-1009")

# Impute RSID for SWB SNP
merged_pp_outcome <- rename(merged_pp_outcome, c(
  "pos.exposure"="position",))
merged_pp_outcome <- merge(x = merged_pp_outcome, y =swb_exposure [ , c("pos.exposure", "SNP")], by = "pos.exposure", all.x=TRUE)
merged_pp_outcome <- rename(merged_pp_outcome, c(
  "position"="pos.exposure",  
  "SNP"="SNP.y"))

merged_pp_outcome = subset(merged_pp_outcome, select = -c(SNP.x,exposure) )

harmonised_data_swb_pp <- harmonise_data(exposure_dat = swb_exposure, outcome_dat = merged_pp_outcome, action=2)
swb_pp_results <- mr(dat = harmonised_data_swb_pp)

# Neuroticism exposure - 
neuroticism_exposure <- extract_instruments(outcomes = "ieu-a-1007")

# Impute RSID for neuroticism SNP
merged_pp_outcome <- rename(merged_pp_outcome, c(
  "pos.exposure"="position",))
merged_pp_outcome <- merge(x = merged_pp_outcome, y =neuroticism_exposure [ , c("pos.exposure", "SNP")], by = "pos.exposure", all.x=TRUE)
merged_pp_outcome <- rename(merged_pp_outcome, c(
  "position"="pos.exposure",  
  "SNP"="SNP.y"))

merged_pp_outcome = subset(merged_pp_outcome, select = -c(SNP.x) )

harmonised_data_neuroticism_pp <- harmonise_data(exposure_dat = neuroticism_exposure, outcome_dat = merged_pp_outcome, action=2)
neuroticism_pp_results <- mr(dat = harmonised_data_neuroticism_pp)

# Depressive symtpoms exposure - 
depressive_symp_exposure <- extract_instruments(outcomes = "ieu-a-1000")

# Impute RSID for depression SNPs

merged_pp_outcome <- rename(merged_pp_outcome, c(
  "pos.exposure"="position",))
merged_pp_outcome <- merge(x = merged_pp_outcome, y =depressive_symp_exposure [ , c("pos.exposure", "SNP")], by = "pos.exposure", all.x=TRUE)
merged_pp_outcome <- rename(merged_pp_outcome, c(
  "position"="pos.exposure",  
  "SNP"="SNP.y"))

merged_pp_outcome = subset(merged_pp_outcome, select = -c(SNP.x) )

harmonised_data_depressive_symp_pp <- harmonise_data(exposure_dat = depressive_symp_exposure, outcome_dat = merged_pp_outcome, action=2)
depressive_symp_pp_results <- mr(dat = harmonised_data_depressive_symp_pp)

# Send results to excel
# Need to add back in SWB and depression if we can get some output

results_summary_PP_out <- rbind(swb_pp_results,
                         neuroticism_pp_results,
                         depressive_symp_pp_results)

results_summary_PP_out$id.outcome = "PP"

write.xlsx(results_summary_PP_out, file = "Results/PP_outcome.xlsx")

################################################################################################
# Example lines of code

#alcohol_consumption <- extract_instruments(outcomes = "ieu-b-73")
#alcohol_consumption <- clump_data(alcohol_consumption)
#lung_cancer <- extract_outcome_data(snps = alcohol_consumption$SNP, outcomes = "ieu-a-966")
#harmonised_data <- harmonise_data(exposure_dat = alcohol_consumption, outcome_dat = lung_cancer)
#alcohol_lung <- mr(dat = harmonised_data)

# 5 lines of code to get MR results, 6 if you want them in odds ratios! Not bad