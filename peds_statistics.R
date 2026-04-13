
# Clear all variables from the environment
rm(list=ls())
# Load required libraries for pedigree analysis and data manipulation
library(kinship2)
library(dplyr)
library(xtable)

#This script analyzes pedigree (family tree) data to generate descriptive statistics for research 
#Overall goal: Provide comprehensive statistical summaries and visualizations of the pedigree structure

# Load family structure information
ped=read.table("/home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/P1/paper3_Analyzes/WGS_subset.ped", header = T)
pedigree_data=ped[c(1:5)]

# Function to compute the number of generations for each individual
# by traversing up the paternal and maternal lineages
compute_generations <- function(data, individual) {
  # Initialize generation counters for paternal and maternal lines
  generations1=generations2=0
  # Start with the given individual
  current_id <- individual
  
  # Traverse paternal lineage (through fathers) until reaching a founder (ID=0)
  while (TRUE) {
    if (current_id == 0) {
      break
    }
    generations1 <- generations1 + 1
    # Move to the father of the current individual
    current_id <- data[data$IID == current_id, "dadid"]
  }
  
  # Reset to given individual and traverse maternal lineage (through mothers)
  current_id <- individual
  while (TRUE) {
    if (current_id == 0) {
      break
    }
    generations2 <- generations2 + 1
    # Move to the mother of the current individual
    current_id <- data[data$IID == current_id, "momid"]
  }
  
  # Return the maximum generation depth from either lineage
  return(max(generations1, generations2))
}

# Compute the number of generations for each individual
pedigree_data$generations <- sapply(pedigree_data$IID, function(x) compute_generations(pedigree_data, x))

# Calculate descriptive statistics for each family
stats <- pedigree_data %>%
  # Group data by family ID
  group_by(FID) %>%
  summarise(
    # Total number of subjects in family
    n_subjects = n(),
    # Count of females (sex=2)
    n_female = sum(sex == 2, na.rm = TRUE),
    # Count of males (sex=1)
    n_male = sum(sex == 1, na.rm = TRUE),
    # Maximum generation depth in family
    max_generation = max(generations)
  )


# Identify founders: individuals with no recorded parents (dadid=0 or NA and momid=0 or NA)
founders <- filter(pedigree_data, (dadid == "0" | is.na(dadid)) & (momid == "0" | is.na(momid)))

# Identify nonfounders: individuals with at least one parent recorded in pedigree
nonfounders <- filter(pedigree_data, !(dadid == "0" | is.na(dadid)) | !(momid == "0" | is.na(momid)))

# Count the number of founders and nonfounders
f <- nrow(founders)
n <- nrow(nonfounders)

# Print summary statistics for the pedigree
cat("Number of founders:", f, "\n")
cat("Number of nonfounders:", n, "\n")
cat("Number of individuals:", length(unique(pedigree_data$IID)), "\n")
cat("Number of families:", length(unique(pedigree_data$FID)), "\n")
# Display the distribution of families by generation
table(stats$max_generation)


# Calculate summary statistics (min, quartiles, median, mean, max) for each variable
stats=stats[,2:ncol(stats)]
summary_stats <- data.frame(
  # Variable names
  variable = c("n_subjects", "n_female", "n_male", "max_generation"),
  # Sum of all values across families
  N = sapply(stats, sum, na.rm = TRUE),
  min = sapply(stats, min, na.rm = TRUE),
  Q1 = sapply(stats, function(x) quantile(x, 0.25, na.rm = TRUE)),
  median = sapply(stats, median, na.rm = TRUE),
  mean = sapply(stats, mean, na.rm = TRUE),
  Q3 = sapply(stats, function(x) quantile(x, 0.75, na.rm = TRUE)),
  max = sapply(stats, max, na.rm = TRUE)
)

# Convert the summary statistics to LaTeX format for publication
latex_table <- xtable(summary_stats)
print(latex_table, type = "latex")



# Load a subset of large families that were split into smaller sub-families
ped_cuted=read.table("/home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/SFTP hub/PHD -P1/Data/sous_familles.ped", header = FALSE)
# Assign standard pedigree column names
colnames(ped_cuted)= c("FID","IID", "dadid", "momid", "sex", "affected")
#ped2=read.table("/home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/P1/paper3_Analyzes/Sous_Familles_tri.ped", header = FALSE)
#file in /PHD -P1/Data/DataMerlin/sous_familles.ped is NOT the latest version of sous_familles.ped, the latest contain IID= GS60704629 of FID= GSF9249_1
# Summary statistics: number of subjects, females, males per family, and maximum generation depth

# Extract the original family ID from the modified FID column 
# This identifies which large families were split into sub-families
bigF <- sapply(strsplit(ped_cuted$FID, "_"), `[`, 1)
cat("Number of big families:", length(unique(bigF)), "\n")

# Descriptive statistics for large families before they were split
# Extract relevant columns from original pedigree data
pedigree_data=ped[c(1:5)]
# Filter to include only individuals from the original large families
pedigree_data =pedigree_data[pedigree_data$FID  %in% bigF,] 

# Compute the number of generations for each individual in large families
pedigree_data$generations <- sapply(pedigree_data$IID, function(x) compute_generations(pedigree_data, x))

# Calculate statistics for each large family (pre-split)
stats <- pedigree_data %>%
  # Group by family ID
  group_by(FID) %>%
  summarise(
    n_subjects = n(),
    n_female = sum(sex == 2, na.rm = TRUE),
    n_male = sum(sex == 1, na.rm = TRUE),
    max_generation = max(generations)
  )

# Display the distribution of families by generation
table(stats$max_generation)

# Calculate summary statistics for large families before splitting
stats=stats[,2:ncol(stats)]
summary_stats <- data.frame(
  variable = c("Number of Individuals", "Number of Females", "Number of Males", "Maximum Generation"),
  N =sapply(stats, sum, na.rm = TRUE),
  min = sapply(stats, min, na.rm = TRUE),
  Q1 = sapply(stats, function(x) quantile(x, 0.25, na.rm = TRUE)),
  median = sapply(stats, median, na.rm = TRUE),
  mean = sapply(stats, mean, na.rm = TRUE),
  Q3 = sapply(stats, function(x) quantile(x, 0.75, na.rm = TRUE)),
  max = sapply(stats, max, na.rm = TRUE)
)

# Create LaTeX table object with 1 decimal place precision
latex_table <- xtable(summary_stats, digits = 1)

# Print LaTeX table without row names
print(latex_table, include.rownames = FALSE, digits = 1)


# Descriptive statistics for large families after being split into smaller sub-families
# Extract relevant columns from the split pedigree data
pedigree_data=ped_cuted[c(1:5)]

# Function to compute generations for individuals in split families
# This version accounts for the modified family structure after splitting
compute_generations2 <- function(data, i) {
  generations1=generations2=0
  individual=data$IID[i]
  current_id <- individual
  data=data[data$FID %in% data$FID[i],]
  
  # Traverse paternal lineage until reaching a founder (ID=0)
  while (TRUE) {
    if (current_id == 0) {
      break
    }
    generations1 <- generations1 + 1
    # Move to the father of the current individual
    current_id <- data[data$IID == current_id, "dadid"]
  }
  
  # Reset to given individual and traverse maternal lineage
  current_id <- individual
  while (TRUE) {
    if (current_id == 0) {
      break
    }
    generations2 <- generations2 + 1
    # Move to the mother of the current individual
    current_id <- data[data$IID == current_id, "momid"]
  }
  
  # Return the maximum generation depth from either lineage
  return(max(generations1, generations2))
}


# Compute the number of generations for each individual in split families
pedigree_data$generations <- sapply(1:nrow(pedigree_data), function(x) compute_generations2(pedigree_data, x))

# Calculate statistics for each sub-family (post-split)
stats <- pedigree_data %>%
  # Group by family ID
  group_by(FID) %>%
  summarise(
    # Total number of subjects
    n_subjects = n(),
    # Count of females (coded as "F" in this dataset)
    n_female = sum(sex == "F", na.rm = TRUE),
    # Count of males (coded as "M" in this dataset)
    n_male = sum(sex == "M", na.rm = TRUE),
    # Maximum generation depth
    max_generation = max(generations)
  )

# Calculate summary statistics for split families
stats=stats[,2:ncol(stats)]
summary_stats <- data.frame(
  variable = c("Number of Individuals", "Number of Females", "Number of Males", "Maximum Generation"),
  # Sum of all values
  N =sapply(stats, sum, na.rm = TRUE),
  min = sapply(stats, min, na.rm = TRUE),
  Q1 = sapply(stats, function(x) quantile(x, 0.25, na.rm = TRUE)),
  median = sapply(stats, median, na.rm = TRUE),
  mean = sapply(stats, mean, na.rm = TRUE),
  Q3 = sapply(stats, function(x) quantile(x, 0.75, na.rm = TRUE)),
  max = sapply(stats, max, na.rm = TRUE)
)

# Display the distribution of sub-families by generation
table(stats$max_generation)

# Create LaTeX table object with 1 decimal place precision
latex_table <- xtable(summary_stats, digits = 1)

# Print LaTeX table without row names
print(latex_table, include.rownames = FALSE, digits = 1)

# Load kinship2 library for pedigree visualization
library(kinship2)
subfamily1_ped=ped[ped$FID%in%c("GSF0543"),]
subfamily1_ped[subfamily1_ped==0]=NA
# Create a pedigree object with individual IDs, father IDs, mother IDs, and sex information
f11=pedigree(id=subfamily1_ped$IID, dadid=subfamily1_ped$dadid, momid=subfamily1_ped$momid, sex=subfamily1_ped$sex)
# Plot the pedigree structure
plot(f11)
