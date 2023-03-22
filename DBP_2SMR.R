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
#install.packages("readxl")

# Load Packages
library(TwoSampleMR)
library(ieugwasr)
library(MendelianRandomization)
library(MRInstruments)
library(MRlap)
library(dplyr)
library(readxl)

ieugwasr::check_access_token()

################################################################################
                      # Diastolic blood pressure exposure #
################################################################################
# Extract DBP GWAS

gwasinfo("ieu-b-39")
dbp_tophits <- tophits(id="ieu-b-39", clump=0)

dbp_gwas <- rename(dbp_tophits, c(
  "SNP"="rsid",
  "effect_allele.exposure"="ea",
  "other_allele.exposure"="nea",
  "beta.exposure"="beta",
  "se.exposure"="se",
  "eaf.exposure"="eaf",
  "pval.exposure"="p",
  "N"="n"))

dbp_gwas$exposure = "DBP"

dbp_liberal <- clump_data(dbp_gwas, clump_kb = 10000, clump_r2 = 0.05,)
dbp_stringent <- clump_data(dbp_gwas, clump_kb = 10000, clump_r2 = 0.001,)

dbp_liberal$id.exposure = "DBP - liberal clump"
dbp_stringent$id.exposure = "DBP - stringent/standard clump"

# Calculate R2 and F stat for exposure data
# Liberal DBP F stat
dbp_liberal$r2 <- (2 * (dbp_liberal$beta.exposure^2) * dbp_liberal$eaf.exposure * (1 - dbp_liberal$eaf.exposure)) /
  (2 * (dbp_liberal$beta.exposure^2) * dbp_liberal$eaf.exposure * (1 - dbp_liberal$eaf.exposure) +
     2 * dbp_liberal$N * dbp_liberal$eaf.exposure * (1 - dbp_liberal$eaf.exposure) * dbp_liberal$se.exposure^2)

dbp_liberal$F <- dbp_liberal$r2 * (dbp_liberal$N - 2) / (1 - dbp_liberal$r2)
dbp_liberal_meanF <- mean(dbp_liberal$F)

# Stringent DBP F stat
dbp_stringent$r2 <- (2 * (dbp_stringent$beta.exposure^2) * dbp_stringent$eaf.exposure * (1 - dbp_stringent$eaf.exposure)) /
  (2 * (dbp_stringent$beta.exposure^2) * dbp_stringent$eaf.exposure * (1 - dbp_stringent$eaf.exposure) +
     2 * dbp_stringent$N * dbp_stringent$eaf.exposure * (1 - dbp_stringent$eaf.exposure) * dbp_stringent$se.exposure^2)

dbp_stringent$F <- dbp_stringent$r2 * (dbp_stringent$N - 2) / (1 - dbp_stringent$r2)
dbp_stringent_meanF <- mean(dbp_stringent$F)

# Find SNPs in outcome data
ao <- available_outcomes()

anxiety_dbp_liberal <- extract_outcome_data(snps = dbp_liberal$SNP, outcomes = "ukb-b-11311")
anxiety_dbp_stringent <- extract_outcome_data(snps = dbp_stringent$SNP, outcomes = "ukb-b-11311")

swb_dbp_liberal <- extract_outcome_data(snps = dbp_liberal$SNP, outcomes = "ieu-a-1009")
swb_dbp_stringent <- extract_outcome_data(snps = dbp_stringent$SNP, outcomes = "ieu-a-1009")

neuroticism_dbp_liberal <- extract_outcome_data(snps = dbp_liberal$SNP, outcomes = "ieu-a-1007")
neuroticism_dbp_stringent <- extract_outcome_data(snps = dbp_stringent$SNP, outcomes = "ieu-a-1007")

depressive_symp_dbp_liberal <- extract_outcome_data(snps = dbp_liberal$SNP, outcomes = "ieu-a-1000")
depressive_symp_dbp_stringent <- extract_outcome_data(snps = dbp_stringent$SNP, outcomes = "ieu-a-1000")

# Harmonise datasets
harm_anxiety_dbp_liberal <- harmonise_data(exposure_dat = dbp_liberal, outcome_dat = anxiety_dbp_liberal)
harm_anxiety_dbp_stringent <- harmonise_data(exposure_dat = dbp_stringent, outcome_dat = anxiety_dbp_stringent)

harm_swb_dbp_liberal <- harmonise_data(exposure_dat = dbp_liberal, outcome_dat = swb_dbp_liberal)
harm_swb_dbp_stringent <- harmonise_data(exposure_dat = dbp_stringent, outcome_dat = swb_dbp_stringent)

harm_neuroticism_dbp_liberal <- harmonise_data(exposure_dat = dbp_liberal, outcome_dat = neuroticism_dbp_liberal)
harm_neuroticism_dbp_stringent <- harmonise_data(exposure_dat = dbp_stringent, outcome_dat = neuroticism_dbp_stringent)

harm_depressive_symp_dbp_liberal <- harmonise_data(exposure_dat = dbp_liberal, outcome_dat = depressive_symp_dbp_liberal)
harm_depressive_symp_dbp_stringent <- harmonise_data(exposure_dat = dbp_stringent, outcome_dat = depressive_symp_dbp_stringent)

# Run and view MR
# Consider setting random seed so analyses/results are completely reproducible

anxiety_dbp_liberal_results <- mr(dat = harm_anxiety_dbp_liberal)
anxiety_dbp_stringent_results <- mr(dat = harm_anxiety_dbp_stringent)

swb_dbp_liberal_results <- mr(dat = harm_swb_dbp_liberal)
swb_dbp_stringent_results <- mr(dat = harm_swb_dbp_stringent)

neuroticism_dbp_liberal_results <- mr(dat = harm_neuroticism_dbp_liberal)
neuroticism_dbp_stringent_results <- mr(dat = harm_neuroticism_dbp_stringent)

depressive_symp_dbp_liberal_results <- mr(dat = harm_depressive_symp_dbp_liberal)
depressive_symp_dbp_stringent_results <- mr(dat = harm_depressive_symp_dbp_stringent)

# Send results to excel

results_summary <- rbind(anxiety_dbp_liberal_results , 
                         anxiety_dbp_stringent_results,
                         swb_dbp_liberal_results,
                         swb_dbp_stringent_results,
                         neuroticism_dbp_liberal_results,
                         neuroticism_dbp_stringent_results, 
                         depressive_symp_dbp_liberal_results,
                         depressive_symp_dbp_stringent_results)

results_summary <- results_summary %>%
  mutate(id.outcome = recode(id.outcome, 'ukb-b-11311' = 'Anxiety', 'ieu-a-1009' = 'SWB', 'ieu-a-1007' = 'Neuroticism', 'ieu-a-1000' ='Depressive symptoms' ))

write.xlsx(results_summary, file = "Results/DBP_exposure.xlsx")


# Calculate F stats for each individual analysis
anxiety_dbp_liberal_meanF <- mean(harm_anxiety_dbp_liberal$F)
anxiety_dbp_stringent_meanF <- mean(harm_anxiety_dbp_stringent$F)

depressive_symp_dbp_liberal_meanF <- mean(harm_depressive_symp_dbp_liberal$F)
depressive_symp_dbp_stringent_meanF <- mean(harm_depressive_symp_dbp_stringent$F)

neuroticism_dbp_liberal_meanF <- mean(harm_neuroticism_dbp_liberal$F)
neuroticism_dbp_stringent_meanF <- mean(harm_neuroticism_dbp_stringent$F)

swb_dbp_liberal_meanF <- mean(harm_swb_dbp_liberal$F)
swb_dbp_stringent_meanF <- mean(harm_swb_dbp_stringent$F)

################################################################################
                        # Psychiatric traits exposure #
################################################################################

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

DBP_anxiety_outcome <- extract_outcome_data(snps = anxiety_exposure$SNP, outcomes = "ieu-b-39")
harmonised_data_anxiety_dbp <- harmonise_data(exposure_dat = anxiety_exposure, outcome_dat = DBP_anxiety_outcome)
anxiety_dbp_results <- mr(dat = harmonised_data_anxiety_dbp)

# Subjective wellbeing exposure - 
swb_exposure <- extract_instruments(outcomes = "ieu-a-1009")
swb_exposure$exposure = "SWB"
swb_exposure$id.exposure = "SWB"

DBP_swb_outcome <- extract_outcome_data(snps = swb_exposure$SNP, outcomes = "ieu-b-39")
harmonised_data_swb_dbp <- harmonise_data(exposure_dat = swb_exposure, outcome_dat = DBP_swb_outcome)
swb_dbp_results <- mr(dat = harmonised_data_swb_dbp)

# Neuroticism exposure - 
neuroticism_exposure <- extract_instruments(outcomes = "ieu-a-1007")
neuroticism_exposure$exposure = "Neuroticism"
neuroticism_exposure$id.exposure = "Neuroticism"

DBP_neuroticism_outcome <- extract_outcome_data(snps = neuroticism_exposure$SNP, outcomes = "ieu-b-39")
harmonised_data_neuroticism_dbp <- harmonise_data(exposure_dat = neuroticism_exposure, outcome_dat = DBP_neuroticism_outcome)
neuroticism_dbp_results <- mr(dat = harmonised_data_neuroticism_dbp)

# Depressive symtpoms exposure - 
depressive_symp_exposure <- extract_instruments(outcomes = "ieu-a-1000")
depressive_symp_exposure$exposure = "Depressive symptoms"
depressive_symp_exposure$id.exposure = "Depressive symptoms"

DBP_depressive_symp_outcome <- extract_outcome_data(snps = depressive_symp_exposure$SNP, outcomes = "ieu-b-39")
harmonised_data_depressive_symp_dbp <- harmonise_data(exposure_dat = depressive_symp_exposure, outcome_dat = DBP_depressive_symp_outcome)
depressive_symp_dbp_results <- mr(dat = harmonised_data_depressive_symp_dbp)

# Send results to excel

results_summary_dbp_out <- rbind(anxiety_dbp_results,
                         swb_dbp_results,
                         neuroticism_dbp_results,
                         depressive_symp_dbp_results)

results_summary_dbp_out$id.outcome = "DBP"

write.xlsx(results_summary_anxiety_exp, file = "Results/DBP_outcome.xlsx")

################################################################################################
# Example lines of code

#alcohol_consumption <- extract_instruments(outcomes = "ieu-b-73")
#alcohol_consumption <- clump_data(alcohol_consumption)
#lung_cancer <- extract_outcome_data(snps = alcohol_consumption$SNP, outcomes = "ieu-a-966")
#harmonised_data <- harmonise_data(exposure_dat = alcohol_consumption, outcome_dat = lung_cancer)
#alcohol_lung <- mr(dat = harmonised_data)

# 5 lines of code to get MR results, 6 if you want them in odds ratios! Not bad
