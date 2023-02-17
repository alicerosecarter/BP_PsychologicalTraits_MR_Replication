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

# Load Packages
library(TwoSampleMR)
library(ieugwasr)
library(MendelianRandomization)
library(MRInstruments)
library(MRlap)
library(dplyr)
library(xlsx)

ieugwasr::check_access_token()

################################################################################
                      # Systolic blood pressure exposure #
################################################################################
# Extract hypertension GWAS

gwasinfo("ukb-b-12493")
hyperten_tophits <- tophits(id="ukb-b-12493", clump=0)

hyperten_gwas <- rename(hyperten_tophits, c(
  "SNP"="rsid",
  "effect_allele.exposure"="ea",
  "other_allele.exposure"="nea",
  "beta.exposure"="beta",
  "se.exposure"="se",
  "eaf.exposure"="eaf",
  "pval.exposure"="p",
  "N"="n"))

hyperten_gwas$exposure = "hyperten"

hyperten_liberal <- clump_data(hyperten_gwas, clump_kb = 10000, clump_r2 = 0.05,)
hyperten_stringent <- clump_data(hyperten_gwas, clump_kb = 10000, clump_r2 = 0.001,)

hyperten_liberal$id.exposure = "hyperten - liberal clump"
hyperten_stringent$id.exposure = "hyperten - stringent/standard clump"

# Calculate R2 and F stat for exposure data
# Liberal hypertension F stat
hyperten_liberal$r2 <- (2 * (hyperten_liberal$beta.exposure^2) * hyperten_liberal$eaf.exposure * (1 - hyperten_liberal$eaf.exposure)) /
  (2 * (hyperten_liberal$beta.exposure^2) * hyperten_liberal$eaf.exposure * (1 - hyperten_liberal$eaf.exposure) +
     2 * hyperten_liberal$N * hyperten_liberal$eaf.exposure * (1 - hyperten_liberal$eaf.exposure) * hyperten_liberal$se.exposure^2)

hyperten_liberal$F <- hyperten_liberal$r2 * (hyperten_liberal$N - 2) / (1 - hyperten_liberal$r2)
hyperten_liberal_meanF <- mean(hyperten_liberal$F)

# Stringent hypertension F stat
hyperten_stringent$r2 <- (2 * (hyperten_stringent$beta.exposure^2) * hyperten_stringent$eaf.exposure * (1 - hyperten_stringent$eaf.exposure)) /
  (2 * (hyperten_stringent$beta.exposure^2) * hyperten_stringent$eaf.exposure * (1 - hyperten_stringent$eaf.exposure) +
     2 * hyperten_stringent$N * hyperten_stringent$eaf.exposure * (1 - hyperten_stringent$eaf.exposure) * hyperten_stringent$se.exposure^2)

hyperten_stringent$F <- hyperten_stringent$r2 * (hyperten_stringent$N - 2) / (1 - hyperten_stringent$r2)
hyperten_stringent_meanF <- mean(hyperten_stringent$F)

# Find SNPs in outcome data
ao <- available_outcomes()

anxiety_hyperten_liberal <- extract_outcome_data(snps = hyperten_liberal$SNP, outcomes = "ukb-b-11311")
anxiety_hyperten_stringent <- extract_outcome_data(snps = hyperten_stringent$SNP, outcomes = "ukb-b-11311")

swb_hyperten_liberal <- extract_outcome_data(snps = hyperten_liberal$SNP, outcomes = "ieu-a-1009")
swb_hyperten_stringent <- extract_outcome_data(snps = hyperten_stringent$SNP, outcomes = "ieu-a-1009")

neuroticism_hyperten_liberal <- extract_outcome_data(snps = hyperten_liberal$SNP, outcomes = "ieu-a-1007")
neuroticism_hyperten_stringent <- extract_outcome_data(snps = hyperten_stringent$SNP, outcomes = "ieu-a-1007")

depressive_symp_hyperten_liberal <- extract_outcome_data(snps = hyperten_liberal$SNP, outcomes = "ieu-a-1000")
depressive_symp_hyperten_stringent <- extract_outcome_data(snps = hyperten_stringent$SNP, outcomes = "ieu-a-1000")

# Harmonise datasets
harm_anxiety_hyperten_liberal <- harmonise_data(exposure_dat = hyperten_liberal, outcome_dat = anxiety_hyperten_liberal)
harm_anxiety_hyperten_stringent <- harmonise_data(exposure_dat = hyperten_stringent, outcome_dat = anxiety_hyperten_stringent)

harm_swb_hyperten_liberal <- harmonise_data(exposure_dat = hyperten_liberal, outcome_dat = swb_hyperten_liberal)
harm_swb_hyperten_stringent <- harmonise_data(exposure_dat = hyperten_stringent, outcome_dat = swb_hyperten_stringent)

harm_neuroticism_hyperten_liberal <- harmonise_data(exposure_dat = hyperten_liberal, outcome_dat = neuroticism_hyperten_liberal)
harm_neuroticism_hyperten_stringent <- harmonise_data(exposure_dat = hyperten_stringent, outcome_dat = neuroticism_hyperten_stringent)

harm_depressive_symp_hyperten_liberal <- harmonise_data(exposure_dat = hyperten_liberal, outcome_dat = depressive_symp_hyperten_liberal)
harm_depressive_symp_hyperten_stringent <- harmonise_data(exposure_dat = hyperten_stringent, outcome_dat = depressive_symp_hyperten_stringent)


# F stats for each analysis
anxiety_hyperten_liberal_meanF <- mean(harm_anxiety_hyperten_liberal$F)
anxiety_hyperten_stringent_meanF <- mean(harm_anxiety_hyperten_stringent$F)

depressive_symp_hyperten_liberal_meanF <- mean(harm_depressive_symp_hyperten_liberal$F)
depressive_symp_hyperten_stringent_meanF <- mean(harm_depressive_symp_hyperten_stringent$F)

neuroticism_hyperten_liberal_meanF <- mean(harm_neuroticism_hyperten_liberal$F)
neuroticism_hyperten_stringent_meanF <- mean(harm_neuroticism_hyperten_stringent$F)

swb_hyperten_liberal_meanF <- mean(harm_swb_hyperten_liberal$F)
swb_hyperten_stringent_meanF <- mean(harm_swb_hyperten_stringent$F)

# Run MR
anxiety_hyperten_liberal_results <- mr(dat = harm_anxiety_hyperten_liberal)
anxiety_hyperten_stringent_results <- mr(dat = harm_anxiety_hyperten_stringent)

swb_hyperten_liberal_results <- mr(dat = harm_swb_hyperten_liberal)
swb_hyperten_stringent_results <- mr(dat = harm_swb_hyperten_stringent)

neuroticism_hyperten_liberal_results <- mr(dat = harm_neuroticism_hyperten_liberal)
neuroticism_hyperten_stringent_results <- mr(dat = harm_neuroticism_hyperten_stringent)

depressive_symp_hyperten_liberal_results <- mr(dat = harm_depressive_symp_hyperten_liberal)
depressive_symp_hyperten_stringent_results <- mr(dat = harm_depressive_symp_hyperten_stringent)

# Send results to excel

results_summary <- rbind(anxiety_hyperten_liberal_results , 
                         anxiety_hyperten_stringent_results,
                         swb_hyperten_liberal_results,
                         swb_hyperten_stringent_results,
                         neuroticism_hyperten_liberal_results,
                         neuroticism_hyperten_stringent_results, 
                         depressive_symp_hyperten_liberal_results,
                         depressive_symp_hyperten_stringent_results)

results_summary <- results_summary %>%
  mutate(id.outcome = recode(id.outcome, 'ukb-b-11311' = 'Anxiety', 'ieu-a-1009' = 'SWB', 'ieu-a-1007' = 'Neuroticism', 'ieu-a-1000' ='Depressive symptoms' ))

write.xlsx(results_summary, file = "Results/hyperten_exposure.xlsx")

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

hyperten_anxiety_outcome <- extract_outcome_data(snps = anxiety_exposure$SNP, outcomes = "ukb-b-12493")
harmonised_data_anxiety_hyperten <- harmonise_data(exposure_dat = anxiety_exposure, outcome_dat = hyperten_anxiety_outcome)
anxiety_hyperten_results <- mr(dat = harmonised_data_anxiety_hyperten)

# Subjective wellbeing exposure - 
swb_exposure <- extract_instruments(outcomes = "ieu-a-1009")
swb_exposure$exposure = "SWB"
swb_exposure$id.exposure = "SWB"

hyperten_swb_outcome <- extract_outcome_data(snps = swb_exposure$SNP, outcomes = "ukb-b-12493")
harmonised_data_swb_hyperten <- harmonise_data(exposure_dat = swb_exposure, outcome_dat = hyperten_swb_outcome)
swb_hyperten_results <- mr(dat = harmonised_data_swb_hyperten)

# Neuroticism exposure - 
neuroticism_exposure <- extract_instruments(outcomes = "ieu-a-1007")
neuroticism_exposure$exposure = "Neuroticism"
neuroticism_exposure$id.exposure = "Neuroticism"

hyperten_neuroticism_outcome <- extract_outcome_data(snps = neuroticism_exposure$SNP, outcomes = "ukb-b-12493")
harmonised_data_neuroticism_hyperten <- harmonise_data(exposure_dat = neuroticism_exposure, outcome_dat = hyperten_neuroticism_outcome)
neuroticism_hyperten_results <- mr(dat = harmonised_data_neuroticism_hyperten)

# Depressive symtpoms exposure - 
depressive_symp_exposure <- extract_instruments(outcomes = "ieu-a-1000")
depressive_symp_exposure$exposure = "Depressive symptoms"
depressive_symp_exposure$id.exposure = "Depressive symptoms"

hyperten_depressive_symp_outcome <- extract_outcome_data(snps = depressive_symp_exposure$SNP, outcomes = "ukb-b-12493")
harmonised_data_depressive_symp_hyperten <- harmonise_data(exposure_dat = depressive_symp_exposure, outcome_dat = hyperten_depressive_symp_outcome)
depressive_symp_hyperten_results <- mr(dat = harmonised_data_depressive_symp_hyperten)

# Send results to excel

results_summary_hyperten_out <- rbind(anxiety_hyperten_results,
                         swb_hyperten_results,
                         neuroticism_hyperten_results,
                         depressive_symp_hyperten_results)

results_summary_hyperten_out$id.outcome = "hyperten"

write.xlsx(results_summary_hyperten_out, file = "Results/hyperten_outcome.xlsx")

################################################################################################
# Example lines of code

#alcohol_consumption <- extract_instruments(outcomes = "ieu-b-73")
#alcohol_consumption <- clump_data(alcohol_consumption)
#lung_cancer <- extract_outcome_data(snps = alcohol_consumption$SNP, outcomes = "ieu-a-966")
#harmonised_data <- harmonise_data(exposure_dat = alcohol_consumption, outcome_dat = lung_cancer)
#alcohol_lung <- mr(dat = harmonised_data)

# 5 lines of code to get MR results, 6 if you want them in odds ratios! Not bad