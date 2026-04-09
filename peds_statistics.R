# ped=read.table("/home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/P1/Candiate_recombinations/WGS.ped")
# ped_61f=ped[ped$V1 %in% unique(ped$V1[ped$V2 %in% Recom.SC$CHILD]), ]
# colnames(ped_61f)= c("FID","IID", "dadid", "momid", "sex", "affected")
# write.table(ped_61f, file = "/home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/P1/paper3_Analyzes/WGS_subset.ped", sep = " ", quote=F, row.names=F, col.names=T)


rm(list=ls())
library(kinship2)
library(dplyr)
library(xtable)


ped=read.table("/home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/P1/paper3_Analyzes/WGS_subset.ped", header = T)

#compute number of generation
# Create pedigree object
pedigree_data=ped[c(1:5)]
compute_generations <- function(data, individual) {
  generations1=generations2=0
  current_id <- individual
  while (TRUE) {
    if (current_id == 0) {
      break
    }
    generations1 <- generations1 + 1
    current_id <- data[data$IID == current_id, "dadid"]
  }
  
  current_id <- individual
  while (TRUE) {
    if (current_id == 0) {
      break
    }
    generations2 <- generations2 + 1
    current_id <- data[data$IID == current_id, "momid"]
  }
  
  return(max(generations1, generations2))
}

# Compute the number of generations for each individual
pedigree_data$generations <- sapply(pedigree_data$IID, function(x) compute_generations(pedigree_data, x))

# Calculate statistics
stats <- pedigree_data %>%
  group_by(FID) %>%
  summarise(
    n_subjects = n(),
    n_female = sum(sex == 2, na.rm = TRUE),
    n_male = sum(sex == 1, na.rm = TRUE),
    max_generation = max(generations)
  )


# Identify founders
founders <- filter(pedigree_data, (dadid == "0" | is.na(dadid)) & (momid == "0" | is.na(momid)))

# Identify nonfounders
nonfounders <- filter(pedigree_data, !(dadid == "0" | is.na(dadid)) | !(momid == "0" | is.na(momid)))

# Count the number of founders and nonfounders
f <- nrow(founders)
n <- nrow(nonfounders)

# Print the results
cat("Number of founders:", f, "\n")
cat("Number of nonfounders:", n, "\n")
cat("Number of individuals:", length(unique(pedigree_data$IID)), "\n")
cat("Number of families:", length(unique(pedigree_data$FID)), "\n")
# number de family by generation
table(stats$max_generation)


# Calculate summary statistics for each variable
stats=stats[,2:ncol(stats)]
summary_stats <- data.frame(
  variable = c("n_subjects", "n_female", "n_male", "max_generation"),
  N = sapply(stats, sum, na.rm = TRUE),
  min = sapply(stats, min, na.rm = TRUE),
  Q1 = sapply(stats, function(x) quantile(x, 0.25, na.rm = TRUE)),
  median = sapply(stats, median, na.rm = TRUE),
  mean = sapply(stats, mean, na.rm = TRUE),
  Q3 = sapply(stats, function(x) quantile(x, 0.75, na.rm = TRUE)),
  max = sapply(stats, max, na.rm = TRUE)
)

# Convert to LaTeX table
latex_table <- xtable(summary_stats)
print(latex_table, type = "latex")



#Big family ... cutted
ped_cuted=read.table("/home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/SFTP hub/PHD -P1/Data/sous_familles.ped", header = FALSE)
colnames(ped_cuted)= c("FID","IID", "dadid", "momid", "sex", "affected")
#ped2=read.table("/home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/P1/paper3_Analyzes/Sous_Familles_tri.ped", header = FALSE)
#file in /PHD -P1/Data/DataMerlin/sous_familles.ped is NOT the latest version of sous_familles.ped, the latest contain IID= GS60704629 of FID= GSF9249_1
#Number of subjects, number female, male,
#min , Q1, median, mean, Q3, Maxof number of subjects per family, number female per family, number of male per family, number of generation

# "Big Families"
bigF <- sapply(strsplit(ped_cuted$FID, "_"), `[`, 1)
cat("Number of big families:", length(unique(bigF)), "\n")

#Descriptive of pedigree of big familles before cuting them in 2 smaller famillies
pedigree_data=ped[c(1:5)]
pedigree_data =pedigree_data[pedigree_data$FID  %in% bigF,] 

# Compute the number of generations for each individual
pedigree_data$generations <- sapply(pedigree_data$IID, function(x) compute_generations(pedigree_data, x))

# Calculate statistics
stats <- pedigree_data %>%
  group_by(FID) %>%
  summarise(
    n_subjects = n(),
    n_female = sum(sex == 2, na.rm = TRUE),
    n_male = sum(sex == 1, na.rm = TRUE),
    max_generation = max(generations)
  )

# number de family by generation
table(stats$max_generation)


# Calculate summary statistics for each variable
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

# Create the xtable object
latex_table <- xtable(summary_stats, digits = 1)

# Print the xtable object, excluding row names
print(latex_table, include.rownames = FALSE, digits = 1)



#Descriptive of pedigree of big families after cuting them in 2 smaller families
pedigree_data=ped_cuted[c(1:5)]

compute_generations2 <- function(data, i) {
  generations1=generations2=0
  individual=data$IID[i]
  current_id <- individual
  data=data[data$FID %in% data$FID[i],]
  while (TRUE) {
    if (current_id == 0) {
      break
    }
    generations1 <- generations1 + 1
    current_id <- data[data$IID == current_id, "dadid"]
  }
  
  current_id <- individual
  while (TRUE) {
    if (current_id == 0) {
      break
    }
    generations2 <- generations2 + 1
    current_id <- data[data$IID == current_id, "momid"]
  }
  
  return(max(generations1, generations2))
}


# Compute the number of generations for each individual
pedigree_data$generations <- sapply(1:nrow(pedigree_data), function(x) compute_generations2(pedigree_data, x))

# Calculate statistics
stats <- pedigree_data %>%
  group_by(FID) %>%
  summarise(
    n_subjects = n(),
    n_female = sum(sex == "F", na.rm = TRUE),
    n_male = sum(sex == "M", na.rm = TRUE),
    max_generation = max(generations)
  )

# Calculate summary statistics for each variable
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

#generation per sub-family
table(stats$max_generation)

# Create the xtable object
latex_table <- xtable(summary_stats, digits = 1)

# Print the xtable object, excluding row names
print(latex_table, include.rownames = FALSE, digits = 1)






library(kinship2)
subfamily1_ped=ped[ped$FID%in%c("GSF0543"),]
subfamily1_ped[subfamily1_ped==0]=NA
f11=pedigree(id=subfamily1_ped$IID, dadid=subfamily1_ped$dadid, momid=subfamily1_ped$momid, sex=subfamily1_ped$sex)
plot(f11)



