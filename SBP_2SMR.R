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
# Extract sbp GWAS

gwasinfo("ieu-b-38")
sbp_tophits <- tophits(id="ieu-b-38", clump=0)

sbp_gwas <- rename(sbp_tophits, c(
  "SNP"="rsid",
  "effect_allele.exposure"="ea",
  "other_allele.exposure"="nea",
  "beta.exposure"="beta",
  "se.exposure"="se",
  "eaf.exposure"="eaf",
  "pval.exposure"="p",
  "N"="n"))

sbp_gwas$exposure = "SBP"

sbp_liberal <- clump_data(sbp_gwas, clump_kb = 10000, clump_r2 = 0.05,)
sbp_stringent <- clump_data(sbp_gwas, clump_kb = 10000, clump_r2 = 0.001,)

sbp_liberal$id.exposure = "SBP - liberal clump"
sbp_stringent$id.exposure = "SBP - stringent/standard clump"

################################################################################
              # Calculate R2 and F stat for SBP exposure data #
################################################################################

# Liberal SBP F stat
sbp_liberal$r2 <- (2 * (sbp_liberal$beta.exposure^2) * sbp_liberal$eaf.exposure * (1 - sbp_liberal$eaf.exposure)) /
  (2 * (sbp_liberal$beta.exposure^2) * sbp_liberal$eaf.exposure * (1 - sbp_liberal$eaf.exposure) +
     2 * sbp_liberal$N * sbp_liberal$eaf.exposure * (1 - sbp_liberal$eaf.exposure) * sbp_liberal$se.exposure^2)

sbp_liberal$F <- sbp_liberal$r2 * (sbp_liberal$N - 2) / (1 - sbp_liberal$r2)
sbp_liberal_meanF <- mean(sbp_liberal$F)

# Stringent SBP F stat
sbp_stringent$r2 <- (2 * (sbp_stringent$beta.exposure^2) * sbp_stringent$eaf.exposure * (1 - sbp_stringent$eaf.exposure)) /
  (2 * (sbp_stringent$beta.exposure^2) * sbp_stringent$eaf.exposure * (1 - sbp_stringent$eaf.exposure) +
     2 * sbp_stringent$N * sbp_stringent$eaf.exposure * (1 - sbp_stringent$eaf.exposure) * sbp_stringent$se.exposure^2)

sbp_stringent$F <- sbp_stringent$r2 * (sbp_stringent$N - 2) / (1 - sbp_stringent$r2)
sbp_stringent_meanF <- mean(sbp_stringent$F)

################################################################################
                      # Find SNPs in outcome data #
################################################################################

ao <- available_outcomes()

anxiety_sbp_liberal <- extract_outcome_data(snps = sbp_liberal$SNP, outcomes = "ukb-b-11311")
anxiety_sbp_stringent <- extract_outcome_data(snps = sbp_stringent$SNP, outcomes = "ukb-b-11311")

swb_sbp_liberal <- extract_outcome_data(snps = sbp_liberal$SNP, outcomes = "ieu-a-1009")
swb_sbp_stringent <- extract_outcome_data(snps = sbp_stringent$SNP, outcomes = "ieu-a-1009")

neuroticism_sbp_liberal <- extract_outcome_data(snps = sbp_liberal$SNP, outcomes = "ieu-a-1007")
neuroticism_sbp_stringent <- extract_outcome_data(snps = sbp_stringent$SNP, outcomes = "ieu-a-1007")

depressive_symp_sbp_liberal <- extract_outcome_data(snps = sbp_liberal$SNP, outcomes = "ieu-a-1000")
depressive_symp_sbp_stringent <- extract_outcome_data(snps = sbp_stringent$SNP, outcomes = "ieu-a-1000")

# Harmonise datasets
harm_anxiety_sbp_liberal <- harmonise_data(exposure_dat = sbp_liberal, outcome_dat = anxiety_sbp_liberal)
harm_anxiety_sbp_stringent <- harmonise_data(exposure_dat = sbp_stringent, outcome_dat = anxiety_sbp_stringent)

harm_swb_sbp_liberal <- harmonise_data(exposure_dat = sbp_liberal, outcome_dat = swb_sbp_liberal)
harm_swb_sbp_stringent <- harmonise_data(exposure_dat = sbp_stringent, outcome_dat = swb_sbp_stringent)

harm_neuroticism_sbp_liberal <- harmonise_data(exposure_dat = sbp_liberal, outcome_dat = neuroticism_sbp_liberal)
harm_neuroticism_sbp_stringent <- harmonise_data(exposure_dat = sbp_stringent, outcome_dat = neuroticism_sbp_stringent)

harm_depressive_symp_sbp_liberal <- harmonise_data(exposure_dat = sbp_liberal, outcome_dat = depressive_symp_sbp_liberal)
harm_depressive_symp_sbp_stringent <- harmonise_data(exposure_dat = sbp_stringent, outcome_dat = depressive_symp_sbp_stringent)


# F stats for each analysis
anxiety_sbp_liberal_meanF <- mean(harm_anxiety_sbp_liberal$F)
anxiety_sbp_stringent_meanF <- mean(harm_anxiety_sbp_stringent$F)

depressive_symp_sbp_liberal_meanF <- mean(harm_depressive_symp_sbp_liberal$F)
depressive_symp_sbp_stringent_meanF <- mean(harm_depressive_symp_sbp_stringent$F)

neuroticism_sbp_liberal_meanF <- mean(harm_neuroticism_sbp_liberal$F)
neuroticism_sbp_stringent_meanF <- mean(harm_neuroticism_sbp_stringent$F)

swb_sbp_liberal_meanF <- mean(harm_swb_sbp_liberal$F)
swb_sbp_stringent_meanF <- mean(harm_swb_sbp_stringent$F)

################################################################################
                                    # Run MR #
################################################################################
# Consider setting random seed so analyses/results are completely reproducible

anxiety_sbp_liberal_results <- mr(dat = harm_anxiety_sbp_liberal)
anxiety_sbp_stringent_results <- mr(dat = harm_anxiety_sbp_stringent)

swb_sbp_liberal_results <- mr(dat = harm_swb_sbp_liberal)
swb_sbp_stringent_results <- mr(dat = harm_swb_sbp_stringent)

neuroticism_sbp_liberal_results <- mr(dat = harm_neuroticism_sbp_liberal)
neuroticism_sbp_stringent_results <- mr(dat = harm_neuroticism_sbp_stringent)

depressive_symp_sbp_liberal_results <- mr(dat = harm_depressive_symp_sbp_liberal)
depressive_symp_sbp_stringent_results <- mr(dat = harm_depressive_symp_sbp_stringent)

################################################################################
                            # Send results to excel #
################################################################################

results_summary <- rbind(anxiety_sbp_liberal_results , 
                         anxiety_sbp_stringent_results,
                         swb_sbp_liberal_results,
                         swb_sbp_stringent_results,
                         neuroticism_sbp_liberal_results,
                         neuroticism_sbp_stringent_results, 
                         depressive_symp_sbp_liberal_results,
                         depressive_symp_sbp_stringent_results)

results_summary <- results_summary %>%
  mutate(id.outcome = recode(id.outcome, 'ukb-b-11311' = 'Anxiety', 'ieu-a-1009' = 'SWB', 'ieu-a-1007' = 'Neuroticism', 'ieu-a-1000' ='Depressive symptoms' ))

write.xlsx(results_summary, file = "Results/SBP_exposure.xlsx")

################################################################################
                        # Psychiatric traits exposure #
################################################################################

available_outcomes(access_token = ieugwasr::check_access_token())

# Anxiety exposure
################################################################################

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

sbp_anxiety_outcome <- extract_outcome_data(snps = anxiety_exposure$SNP, outcomes = "ieu-b-38")
harmonised_data_anxiety_sbp <- harmonise_data(exposure_dat = anxiety_exposure, outcome_dat = sbp_anxiety_outcome)
anxiety_sbp_results <- mr(dat = harmonised_data_anxiety_sbp)

# Subjective wellbeing exposure - 
################################################################################

swb_exposure <- extract_instruments(outcomes = "ieu-a-1009")
swb_exposure$exposure = "SWB"
swb_exposure$id.exposure = "SWB"

sbp_swb_outcome <- extract_outcome_data(snps = swb_exposure$SNP, outcomes = "ieu-b-38")
harmonised_data_swb_sbp <- harmonise_data(exposure_dat = swb_exposure, outcome_dat = sbp_swb_outcome)
swb_sbp_results <- mr(dat = harmonised_data_swb_sbp)

# Neuroticism exposure - 
################################################################################

neuroticism_exposure <- extract_instruments(outcomes = "ieu-a-1007")
neuroticism_exposure$exposure = "Neuroticism"
neuroticism_exposure$id.exposure = "Neuroticism"

sbp_neuroticism_outcome <- extract_outcome_data(snps = neuroticism_exposure$SNP, outcomes = "ieu-b-38")
harmonised_data_neuroticism_sbp <- harmonise_data(exposure_dat = neuroticism_exposure, outcome_dat = sbp_neuroticism_outcome)
neuroticism_sbp_results <- mr(dat = harmonised_data_neuroticism_sbp)

# Depressive symtpoms exposure - 
################################################################################

depressive_symp_exposure <- extract_instruments(outcomes = "ieu-a-1000")
depressive_symp_exposure$exposure = "Depressive symptoms"
depressive_symp_exposure$id.exposure = "Depressive symptoms"

sbp_depressive_symp_outcome <- extract_outcome_data(snps = depressive_symp_exposure$SNP, outcomes = "ieu-b-38")
harmonised_data_depressive_symp_sbp <- harmonise_data(exposure_dat = depressive_symp_exposure, outcome_dat = sbp_depressive_symp_outcome)
depressive_symp_sbp_results <- mr(dat = harmonised_data_depressive_symp_sbp)

################################################################################
                          # Send results to excel #
################################################################################

results_summary_sbp_out <- rbind(anxiety_sbp_results,
                         swb_sbp_results,
                         neuroticism_sbp_results,
                         depressive_symp_sbp_results)

results_summary_sbp_out$id.outcome = "SBP"

write.xlsx(results_summary_sbp_out, file = "Results/SBP_outcome.xlsx")


################################################################################
        # Calculate R2 and F stat for psychiatric traits exposure data #
################################################################################

# Anxiety SBP F stat
anxiety_exposure$r2 <- (2 * (anxiety_exposure$beta.exposure^2) * anxiety_exposure$eaf.exposure * (1 - anxiety_exposure$eaf.exposure)) /
  (2 * (anxiety_exposure$beta.exposure^2) * anxiety_exposure$eaf.exposure * (1 - anxiety_exposure$eaf.exposure) +
     2 * anxiety_exposure$N * anxiety_exposure$eaf.exposure * (1 - anxiety_exposure$eaf.exposure) * anxiety_exposure$se.exposure^2)

anxiety_exposure$F <- anxiety_exposure$r2 * (anxiety_exposure$N - 2) / (1 - anxiety_exposure$r2)
anxiety_meanF <- mean(anxiety_exposure$F)

# SWB SBP F stat
swb_exposure$r2 <- (2 * (swb_exposure$beta.exposure^2) * swb_exposure$eaf.exposure * (1 - swb_exposure$eaf.exposure)) /
  (2 * (swb_exposure$beta.exposure^2) * swb_exposure$eaf.exposure * (1 - swb_exposure$eaf.exposure) +
     2 * swb_exposure$samplesize.exposure * swb_exposure$eaf.exposure * (1 - swb_exposure$eaf.exposure) * swb_exposure$se.exposure^2)

swb_exposure$F <- swb_exposure$r2 * (swb_exposure$samplesize.exposure - 2) / (1 - swb_exposure$r2)
swb_meanF <- mean(swb_exposure$F)

# Neuroticism SBP F stat
neuroticism_exposure$r2 <- (2 * (neuroticism_exposure$beta.exposure^2) * neuroticism_exposure$eaf.exposure * (1 - neuroticism_exposure$eaf.exposure)) /
  (2 * (neuroticism_exposure$beta.exposure^2) * neuroticism_exposure$eaf.exposure * (1 - neuroticism_exposure$eaf.exposure) +
     2 * neuroticism_exposure$samplesize.exposure * neuroticism_exposure$eaf.exposure * (1 - neuroticism_exposure$eaf.exposure) * neuroticism_exposure$se.exposure^2)

neuroticism_exposure$F <- neuroticism_exposure$r2 * (neuroticism_exposure$samplesize.exposure - 2) / (1 - neuroticism_exposure$r2)
neuroticism_meanF <- mean(neuroticism_exposure$F)

# Depressive symptoms SBP F stat
depressive_symp_exposure$r2 <- (2 * (depressive_symp_exposure$beta.exposure^2) * depressive_symp_exposure$eaf.exposure * (1 - depressive_symp_exposure$eaf.exposure)) /
  (2 * (depressive_symp_exposure$beta.exposure^2) * depressive_symp_exposure$eaf.exposure * (1 - depressive_symp_exposure$eaf.exposure) +
     2 * depressive_symp_exposure$samplesize.exposure * depressive_symp_exposure$eaf.exposure * (1 - depressive_symp_exposure$eaf.exposure) * depressive_symp_exposure$se.exposure^2)

depressive_symp_exposure$F <- depressive_symp_exposure$r2 * (depressive_symp_exposure$samplesize.exposure - 2) / (1 - depressive_symp_exposure$r2)
depressive_symp_meanF <- mean(depressive_symp_exposure$F)


################################################################################################
# Example lines of code

#alcohol_consumption <- extract_instruments(outcomes = "ieu-b-73")
#alcohol_consumption <- clump_data(alcohol_consumption)
#lung_cancer <- extract_outcome_data(snps = alcohol_consumption$SNP, outcomes = "ieu-a-966")
#harmonised_data <- harmonise_data(exposure_dat = alcohol_consumption, outcome_dat = lung_cancer)
#alcohol_lung <- mr(dat = harmonised_data)

# 5 lines of code to get MR results, 6 if you want them in odds ratios! Not bad
