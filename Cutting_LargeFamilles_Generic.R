
# Load necessary library
library(FamAgg)

# Read the pedigree data
structure.famille <- read.table("/home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/P1/Candiate_recombinations/WGS.ped", header = FALSE)

# Rename columns for clarity
colnames(structure.famille) <- c("famille", "Individu", "father", "mother", "sex")

# Filter the pedigree for the specified family
family_id <- "GSF8540"
mbped <- structure.famille[structure.famille$famille == family_id, 1:5]
colnames(mbped) <- c("family", "id", "father", "mother", "sex")

# Convert family and id columns to character type
mbped$family <- as.character(mbped$family)
mbped$id <- as.character(mbped$id)

# Create FAData object
fadi <- FAData(pedigree = mbped, header = FALSE, sep = "\t", id.col = "id",
               family.col = "family", father.col = "father", mother.col = "mother",
               sex.col = "sex")

# Find founders and their children
founders <- findFounders(fadi, family = family_id)
children <- getChildren(fadi, id = founders, max.generations = 1)

# Calculate kinship matrix
kinship_matrix <- kinship(fadi)
kinship_df <- as.data.frame(as.matrix(kinship_matrix))
kinship_df <- kinship_df[names(kinship_df) %in% children]
kinship_df <- kinship_df[row.names(kinship_df) %in% children,]

# Mark upper triangular matrix to "."
for (i in seq_len(nrow(kinship_df))) {
  for (j in seq_len(ncol(kinship_df))) {
    if (i <= j) {
      kinship_df[i, j] <- "."
    }
  }
}

# Identify siblings (kinship 0.250)
siblings_list <- list()
for (i in seq_len(nrow(kinship_df))) {
  tmp <- kinship_df[kinship_df[i, ] == 0.250, i, drop = FALSE]
  if (nrow(tmp) > 0) {
    siblings_list[[i]] <- c(colnames(tmp), rownames(tmp))
  }
}
siblings_list[[length(siblings_list) + 1]] <- children[!children %in% unlist(siblings_list)]

# Filter out empty lists
siblings_list <- Filter(length, siblings_list)

# Determine the splitting point for subfamilies
split_index <- 0
for (i in seq_len(length(siblings_list))) {
  if (sum(lengths(siblings_list[1:i])) > (length(children) / 2)) {
    split_index <- i - 1
    break
  }
}

# Build pedigree for subfamilies
subfamily1_ids <- c(unlist(siblings_list[1:split_index]), getChildren(fadi, id = unlist(siblings_list[1:split_index]), max.generations = 1))
subfamily2_ids <- c(unlist(siblings_list[(split_index + 1):length(siblings_list)]), getChildren(fadi, id = unlist(siblings_list[(split_index + 1):length(siblings_list)]), max.generations = 1))

subfamily1_ped <- buildPed(fadi, id = subfamily1_ids, prune = TRUE, cex = 0.5)
subfamily2_ped <- buildPed(fadi, id = subfamily2_ids, prune = TRUE, cex = 0.5)

# Annotate subfamilies
subfamily1_ped$familyAD <- paste0(subfamily1_ped$family, "_1")
subfamily2_ped$familyAD <- paste0(subfamily2_ped$family, "_2")

# Combine subfamilies into one data frame
fAD <- rbind(subfamily1_ped, subfamily2_ped)
fAD <- cbind(familyAD = fAD[, "familyAD"], fAD[, 1:5])
row.names(fAD) <- NULL

# Output the final annotated pedigree data
fAD



# # Optional: Plot the pedigrees
# f11=pedigree(id=subfamily1_ped$id, dadid=subfamily1_ped$father, momid=subfamily1_ped$mother, sex=subfamily1_ped$sex)
# plot(f11)
# f21=pedigree(id=subfamily2_ped$id, dadid=subfamily2_ped$father, momid=subfamily2_ped$mother, sex=subfamily2_ped$sex)
# plot(f21)
# f=mbped
# f[f==0]=NA
# f0=pedigree(id=f$id, dadid=f$father, momid=f$mother, sex=f$sex)
# plot(f0)




