rm(list=ls())

#This version erasing the old that not include chrom 9.
#This version is adapted for rorqual paths


# get recombination events per family
rpath="~/links/projects/rrg-girardsi/genealogy_sims/results/Samir/"
peds=read.table(paste0(rpath,"/P1/Candiate_recombinations/WGS.ped"))
Recom.SC=read.csv2(paste0(rpath,"/P1/paper3_Analyzes/Recom.SC.csv"), header = T)

Recom.SC2= merge(Recom.SC,peds[,1:2], by.x="CHILD", by.y="V2", all.x=TRUE )
colnames(Recom.SC2)[ncol(Recom.SC2)]="FID"

OverlapByfamily <- Recom.SC2 %>%
  group_by(FID) %>%
  summarise(
    Npair = n(), ## Number recombination WGS
    NOvScMr = sum(ov_Mr > 0), # Number of overlap WGS and Merlin in family
    NOvScOE = sum(ov_OE>0), # Number of overlap WGS Merlin and OE in family
    NOvScMrOE = sum(ov_OE>0 & ov_Mr>0) # Number of overlap WGS Merlin and OE in family
  )


#****************
#Méthode 1: IC for "Shapeit2 recombination events on WGS were detected by Merlin" 
#****************

# Calculate sums
sum_NOvScMr <- sum(OverlapByfamily$NOvScMr)
sum_Npair <- sum(OverlapByfamily$Npair)

# Calculate proportion
proportion <- sum_NOvScMr / sum_Npair

# Calculate standard error
se <- sqrt(proportion * (1 - proportion) / sum_Npair)

# Determine confidence level
confidence_level <- 0.95
z_score <- qnorm(1 - (1 - confidence_level) / 2)

# Calculate margin of error
margin_of_error <- z_score * se

# Calculate confidence interval
ci_lower <- proportion - margin_of_error
ci_upper <- proportion + margin_of_error

# Print results
cat("\n ",round(proportion*100,1), "% (95% CI: ", round(ci_lower*100,1), "% to  ", round(ci_upper*100,1)," %)\n")



#*****************
#Méthode 1: IC for "Proportion of Shapeit2 recombination (in both SNP chip & WGS) detected by Merlin" 
#*****************
# Calculate sums
sum_NOvScMrOE <- sum(OverlapByfamily$NOvScMrOE)
sum_NOvScOE <- sum(OverlapByfamily$NOvScOE)

# Calculate proportion
proportion <- sum_NOvScMrOE / sum_NOvScOE

# Calculate standard error
se <- sqrt(proportion * (1 - proportion) / sum_NOvScOE)

# Determine confidence level
confidence_level <- 0.95
z_score <- qnorm(1 - (1 - confidence_level) / 2)

# Calculate margin of error
margin_of_error <- z_score * se

# Calculate confidence interval
ci_lower <- proportion - margin_of_error
ci_upper <- proportion + margin_of_error

# Print results
cat("\n ",round(proportion*100,1), "% (95% CI: ", round(ci_lower*100,1), "% to  ", round(ci_upper*100,1)," %)\n")




# #*******Method 2***********
# #*
# data=OverlapByfamily
# # Calculate sums
# sum_Npair <- sum(data$Npair)
# sum_NOvScMr <- sum(data$NOvScMr)
# 
# # Calculate the proportion
# p_hat <- sum_NOvScMr / sum_Npair
# # Calculate standard error using the delta method
# var_NOvScMr <- sum((data$NOvScMr - mean(data$NOvScMr))^2) / (nrow(data) - 1)
# var_Npair <- sum((data$Npair - mean(data$Npair))^2) / (nrow(data) - 1)
# cov_NOvScMr_Npair <- cov(data$NOvScMr, data$Npair)
# 
# se <- sqrt((var_NOvScMr / sum_NOvScMr^2) + (var_Npair / sum_Npair^2) - (2 * cov_NOvScMr_Npair / (sum_NOvScMr * sum_Npair)))
# 
# # Confidence interval
# z <- qnorm(0.975) # 95% confidence
# ci_lower <- p_hat - z * se
# ci_upper <- p_hat + z * se
# 
# # Display results
# cat("95% Confidence Interval: [", ci_lower, ", ", ci_upper, "]\n")



#****************
#Méthode 1: IC for "Shapeit2 recombination events on WGS were detected by Merlin" 
#****************


# get recombination events per family
Recom.OE=read.csv2(paste0(rpath,"/P1/paper3_Analyzes/Recom.OE.csv"), header = T)
Recom.OE= merge(Recom.OE,peds[,1:2], by.x="CHILD", by.y="V2", all.x=TRUE )
colnames(Recom.OE)[ncol(Recom.OE)]="FID"

OverlapByfamily <- Recom.OE %>%
  group_by(FID) %>%
  summarise(
    Npair = n(),
    NOvOEMr = sum(ov_Mr > 0), # Number of overlap OE and Merlin in family
  )



# Calculate sums
sum_NOvOEMr <- sum(OverlapByfamily$NOvOEMr)
sum_Npair <- sum(OverlapByfamily$Npair)

# Calculate proportion
proportion <- sum_NOvOEMr / sum_Npair

# Calculate standard error
se <- sqrt(proportion * (1 - proportion) / sum_Npair)

# Determine confidence level
confidence_level <- 0.95
z_score <- qnorm(1 - (1 - confidence_level) / 2)

# Calculate margin of error
margin_of_error <- z_score * se

# Calculate confidence interval
ci_lower <- proportion - margin_of_error
ci_upper <- proportion + margin_of_error

# Print results
cat("\n ",round(proportion*100,1), "% (95% CI: ", round(ci_lower*100,1), "% to  ", round(ci_upper*100,1)," %)\n")

