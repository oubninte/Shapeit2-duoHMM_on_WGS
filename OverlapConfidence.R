
## Objective
# This script is designed to calculate the overlap confidence interval for concordance. 
# The script takes input data of overlaps detected

# Clear all variables from the workspace
rm(list=ls())

# ============================================================================
# SECTION 1: Load and merge recombination data with pedigree information
# ============================================================================

# Define the base path 
rpath="~/links/projects/rrg-girardsi/genealogy_sims/results/Samir/"

# Load pedigree data
peds=read.table(paste0(rpath,"/P1/Candiate_recombinations/WGS.ped"))

# Load recombination events from Shapeit2 and Merlin 
Recom.SC=read.csv2(paste0(rpath,"/P1/paper3_Analyzes/Recom.SC.csv"), header = T)

# Merge recombination data with pedigree to add family IDs
Recom.SC2= merge(Recom.SC,peds[,1:2], by.x="CHILD", by.y="V2", all.x=TRUE )
# Rename the family ID column 
colnames(Recom.SC2)[ncol(Recom.SC2)]="FID"

# Group recombination events by family and calculate overlap statistics
OverlapByfamily <- Recom.SC2 %>%
  group_by(FID) %>%
  summarise(
    Npair = n(), ## Total number of recombination events detected in WGS
    NOvScMr = sum(ov_Mr > 0), # Count of WGS recombinations detected by Merlin
    NOvScOE = sum(ov_OE>0), # Count of recombinations detected by Merlin on SNP chip
    NOvScMrOE = sum(ov_OE>0 & ov_Mr>0) # Count of recombinations detected by all three Method-data
  )


# ============================================================================
# ANALYSIS 1: Confidence interval for Shapeit2 WGS recombinations detected by Merlin
# ============================================================================

# Recombination events  count
sum_NOvScMr <- sum(OverlapByfamily$NOvScMr)

# Calculate the proportion of WGS recombinations detected by Merlin
sum_Npair <- sum(OverlapByfamily$Npair)
proportion <- sum_NOvScMr / sum_Npair

# Calculate the standard error for the binomial proportion
se <- sqrt(proportion * (1 - proportion) / sum_Npair)

# Set the confidence level to 95%
confidence_level <- 0.95
# Calculate the z-score corresponding to the confidence level
z_score <- qnorm(1 - (1 - confidence_level) / 2)

# Calculate the margin of error using z-score and standard error
margin_of_error <- z_score * se

# Calculate the lower bound of the confidence interval
ci_lower <- proportion - margin_of_error
# Calculate the upper bound of the confidence interval
ci_upper <- proportion + margin_of_error

# Print the proportion as percentage with 95% confidence interval
cat("\n ",round(proportion*100,1), "% (95% CI: ", round(ci_lower*100,1), "% to  ", round(ci_upper*100,1)," %)\n")



# ============================================================================
# ANALYSIS 2: Confidence interval for recombinations detected in all three methods
# ============================================================================

# Recombination events detected by all three methods
sum_NOvScMrOE <- sum(OverlapByfamily$NOvScMrOE)
# Recombination events detected by Shapeit2 on both WGS and SNP chip
sum_NOvScOE <- sum(OverlapByfamily$NOvScOE)

# Calculate the proportion of triple-detection events 
proportion <- sum_NOvScMrOE / sum_NOvScOE

# Calculate the standard error for this proportion
se <- sqrt(proportion * (1 - proportion) / sum_NOvScOE)

# Set the confidence level to 95%
confidence_level <- 0.95
# Calculate the z-score for the confidence level
z_score <- qnorm(1 - (1 - confidence_level) / 2)

# Calculate the margin of error
margin_of_error <- z_score * se

# Calculate the confidence interval bounds
ci_lower <- proportion - margin_of_error
ci_upper <- proportion + margin_of_error

# Print the proportion as percentage with 95% confidence interval
cat("\n ",round(proportion*100,1), "% (95% CI: ", round(ci_lower*100,1), "% to  ", round(ci_upper*100,1)," %)\n")


# ============================================================================
# ANALYSIS 3: Confidence interval for Shapeit2 recombinations on SNP chip detected by Merlin
# ============================================================================


# Load recombination events on SNP chip (OE – OmniExpress chip)
Recom.OE=read.csv2(paste0(rpath,"/P1/paper3_Analyzes/Recom.OE.csv"), header = T)
Recom.OE= merge(Recom.OE,peds[,1:2], by.x="CHILD", by.y="V2", all.x=TRUE )
colnames(Recom.OE)[ncol(Recom.OE)]="FID"

OverlapByfamily <- Recom.OE %>%
  group_by(FID) %>%
  summarise(
    Npair = n(), 
    NOvOEMr = sum(ov_Mr > 0), 
  )



# Overlapping recombination events across all families
sum_NOvOEMr <- sum(OverlapByfamily$NOvOEMr)
sum_Npair <- sum(OverlapByfamily$Npair)

# Calculate the proportion of SNP chip recombinations inferred by Shapeit2 detected by Merlin
proportion <- sum_NOvOEMr / sum_Npair

# Calculate the standard error for the binomial proportion
se <- sqrt(proportion * (1 - proportion) / sum_Npair)

# Set the confidence level to 95%
confidence_level <- 0.95
# Calculate the z-score for the confidence level
z_score <- qnorm(1 - (1 - confidence_level) / 2)

# Calculate the margin of error
margin_of_error <- z_score * se

# Calculate the confidence interval bounds
ci_lower <- proportion - margin_of_error
ci_upper <- proportion + margin_of_error

# Print the proportion as percentage with 95% confidence interval
cat("\n ",round(proportion*100,1), "% (95% CI: ", round(ci_lower*100,1), "% to  ", round(ci_upper*100,1)," %)\n")
