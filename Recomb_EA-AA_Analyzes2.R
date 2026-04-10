# Clear all variables and load required libraries for data manipulation and visualization
rm(list=ls())
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggforce)

# This script analyzes recombination events from Shapeit2-duoHMM and Merlin methods
# across African-American (AA) and European-American (EA) populations on WGS data

# Set the base path to project data on rorqual server
rpath="~/links/projects/rrg-girardsi/genealogy_sims/results/Samir/"

# Load pedigree file containing family structure information
ped=read.table(paste0(rpath,"/P1/Candiate_recombinations/WGS.ped"))

# ==================================================
# LOAD AND COMBINE AFRICAN-AMERICAN (AA) RECOMBINATION DATA
# ==================================================
# Define chromosomes to process (1-22 excluding sex chromosomes)
chr=c(1:9,10:22)
# Construct file paths for AA recombination data (Shapeit2-duoHMM on OmniExpress chip)
chemin2Recomb=paste0(rpath,"/SFTP hub/PHD -P1/Data/to delete/data21122020/")
cheminAA = paste0(chemin2Recomb, "chr",chr,".AA.OmniExpress.clean.recombinations")

# Load first chromosome's AA recombination data and add chromosome identifier
Recom.AA= read.delim2(cheminAA[1], header = TRUE, sep = "")
Recom.AA=data.frame("chr"=chr[1], Recom.AA )

# Iteratively load and combine recombination data from remaining chromosomes
for (j in 2:length(chr)) {
  tmp=read.delim2(cheminAA[j], header = TRUE, sep = "")
  Recom.AA = rbind(Recom.AA,data.frame("chr"=chr[j], tmp))
}

# Calculate recombination interval length (END - START + 1)
Recom.AA$"longueur"=Recom.AA$END-Recom.AA$START+1
# Extract unique parent-child pairs in AA population
AA_pair=(unique(paste(Recom.AA$PARENT,Recom.AA$CHILD,sep="-")))
# Filter recombination events by confidence threshold (probability >= 0.5)
Recom.AA=Recom.AA[Recom.AA$PROB_RECOMBINATION >=  0.5,]


# ==================================================
# LOAD AND COMBINE EUROPEAN-AMERICAN (EA) RECOMBINATION DATA
# ==================================================
# Construct file paths for EA recombination data (Shapeit2-duoHMM on OmniExpress chip)
cheminEA = paste0(chemin2Recomb, "chr",chr,".EA.OmniExpress.clean.recombinations")
# Load first chromosome's EA recombination data and add chromosome identifier
Recom.EA= read.delim2(cheminEA[1], header = TRUE, sep = "")
Recom.EA=data.frame("chr"=chr[1], Recom.EA )

# Iteratively load and combine recombination data from remaining chromosomes
for (j in 2:length(chr)) {
  tmp=read.delim2(cheminEA[j], header = TRUE, sep = "")
  Recom.EA = rbind(Recom.EA,data.frame("chr"=chr[j], tmp))
}

# Calculate recombination interval length
Recom.EA$"longueur"=Recom.EA$END-Recom.EA$START+1
# Extract unique parent-child pairs in EA population
EA_pair=(unique(paste(Recom.EA$PARENT,Recom.EA$CHILD,sep="-")))
# Filter recombination events by confidence threshold
Recom.EA=Recom.EA[Recom.EA$PROB_RECOMBINATION >=  0.5,]


# ==================================================
# LOAD AND COMBINE MERLIN RECOMBINATION DATA
# ==================================================
# Create chromosome name labels for file paths
chrNum = paste0("chr", chr)

# Construct file paths for Merlin recombination data
CheminRM=paste0(rpath,"/P1/paper3_Analyzes/Merlin_Recomb/", 
                "chr", chr, ".merlin.recombinaisons")
CheminRM.3innuc=paste0(rpath,"/P1/paper3_Analyzes/Merlin_Recomb/", 
                       "chr", chr, ".merlin.3innuc.recombinaisons")

# Load first chromosome's Merlin recombination 
Recom.Merlin = read.delim2(CheminRM[1], header = TRUE, sep = "") 
Recom.Merlin=rbind(Recom.Merlin, read.table(CheminRM.3innuc[1], header = T) )
Recom.Merlin$"chr"=chr[1]

# Iteratively load and combine Merlin data from remaining chromosomes
for (j in 2:length(chr)) {
  tmp=read.delim2(CheminRM[j], header = TRUE, sep = "")
  tmp=rbind(tmp, read.table(CheminRM.3innuc[j], header = T) )
  tmp$"chr"=chr[j]
  Recom.Merlin=rbind(Recom.Merlin, tmp )
}

# Reset row names for clean indexing
row.names(Recom.Merlin) = NULL
# Remove records with missing START or END positions
Recom.Merlin=Recom.Merlin[!is.na(Recom.Merlin$START),]
Recom.Merlin=Recom.Merlin[!is.na(Recom.Merlin$END),]
# Remove duplicate records
Recom.Merlin=Recom.Merlin[!duplicated(Recom.Merlin),]


# ==================================================
# FILTER TO COMMON PARENT-CHILD PAIRS ACROSS METHODS
# ==================================================
# Identify parent-child pairs present in both Merlin and Shapeit2 
common_indiv <- unique(intersect(c(paste(Recom.AA$PARENT,Recom.AA$CHILD,sep="-"),
                                   paste(Recom.EA$PARENT,Recom.EA$CHILD,sep="-"))
                                 , paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-")))

# Filter all three datasets to retain only common parent-child pairs for fair comparison
Recom.Merlin=Recom.Merlin[paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-") %in% common_indiv , ]
Recom.AA=Recom.AA[paste(Recom.AA$PARENT,Recom.AA$CHILD,sep="-") %in% common_indiv , ]
Recom.EA=Recom.EA[paste(Recom.EA$PARENT,Recom.EA$CHILD,sep="-") %in% common_indiv , ]

# Expected number of recombination events per chromosome based on genetic map
NR_expected_common= data.frame(chr=chr,
                               Expected=round(c(284.8,274.6,233.0,216.4,213.0, 197.5, 190.9,
                                                171.5, 168.5, 175.5, 165.9, 168.1,124.1, 124.9, 127.0,
                                                133.5, 138.8, 121.3, 102.6, 102.6, 63.7, 69.5)/100 , 2))


# ==================================================
# ANALYZE RECOMBINATION COUNTS BY CHROMOSOME AND METHOD
# ==================================================
# Count total recombination events by chromosome for each method
df=merge(as.data.frame(table(Recom.Merlin$chr)), as.data.frame(table(Recom.AA$chr)), by="Var1")
df=merge(df, as.data.frame(table(Recom.EA$chr)), by="Var1")

# Count Merlin recombinations specific to AA parent-child pairs
df=merge(df, as.data.frame(table(Recom.Merlin$chr[paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-") %in% AA_pair])), by="Var1")
# Count Merlin recombinations specific to EA parent-child pairs
df=merge(df, as.data.frame(table(Recom.Merlin$chr[paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-") %in% EA_pair])), by="Var1")

# Standardize column names for clarity
colnames(df)=c("chr", "Merlin_OE", "Shapeit2_AA", "Shapeit2_EA", "Merlin_AA", "Merlin_EA")
# Merge with expected recombination counts
df=merge(df, NR_expected_common, by="chr")

# Normalize recombination counts to per-meiosis rates 
df$Merlin_OE=round(df$Merlin_OE/length(unique(paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-"))),2)
df$Shapeit2_AA=round(df$Shapeit2_AA/length(unique(paste(Recom.AA$PARENT,Recom.AA$CHILD,sep="-"))),2)
df$Shapeit2_EA=round(df$Shapeit2_EA/length(unique(paste(Recom.EA$PARENT,Recom.EA$CHILD,sep="-"))),2)
df$Merlin_AA=round(df$Merlin_AA/sum(AA_pair  %in% paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-")),2)
df$Merlin_EA=round(df$Merlin_EA/sum(EA_pair  %in% paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-")),2)
df

# ==================================================
# VISUALIZE RECOMBINATION RATES BY CHROMOSOME
# ==================================================
# Convert dataframe from wide to long format for ggplot visualization
df_melt <- df %>% pivot_longer(cols = -chr)
df_melt=df_melt[df_melt$name !="Merlin_OE",]
# Create line plot comparing recombination rates across methods by chromosome
p=ggplot(df_melt, aes(chr, value, group=name, colour=name)) +
  geom_line() +
  labs(x="Chromosome number", y="Number of recombination per meiosis", colour=" ")
p


# ==================================================
# STATISTICAL COMPARISON OF RECOMBINATION COUNTS BETWEEN AA AND EA
# ==================================================
# Load pedigree file for family assignment
peds=read.table(paste0(rpath,"/P1/Candiate_recombinations/WGS.ped"))


# Calculate recombination counts per family for AA population
# Merge AA recombinations with pedigree to get family IDs
Recom.AA2= merge(Recom.AA,peds[,1:2], by.x="CHILD", by.y="V2", all.x=TRUE )
colnames(Recom.AA2)[ncol(Recom.AA2)]="FID"

# Summarize recombination statistics by family for AA population
LenByfamily <- Recom.AA2 %>%
  group_by(FID) %>%
  summarise(
    NRecom = n(),  # Total number of recombinations
    num_pairs = n_distinct(paste(PARENT, CHILD, sep = "-"))  # Number of unique parent-child pairs
  )


# Calculate recombination rate (events per meiosis) for each AA family
group1<- LenByfamily$NRecom/LenByfamily$num_pairs


# Calculate recombination counts per family for EA population
# Merge EA recombinations with pedigree to get family IDs
Recom.EA2= merge(Recom.EA,peds[,1:2], by.x="CHILD", by.y="V2", all.x=TRUE )
colnames(Recom.EA2)[ncol(Recom.EA2)]="FID"

# Summarize recombination statistics by family for EA population
LenByfamily <- Recom.EA2 %>%
  group_by(FID) %>%
  summarise(
    NRecom = n(),  # Total number of recombinations
    num_pairs = n_distinct(paste(PARENT, CHILD, sep = "-"))  # Number of unique parent-child pairs
  )


# Calculate recombination rate (events per meiosis) for each EA family
group2<- LenByfamily$NRecom/LenByfamily$num_pairs


# Perform Mann-Whitney-Wilcoxon non-parametric test to compare recombination counts between populations
test_result <- wilcox.test(group1, group2, alternative = "two.sided")

# Display test results
cat("Statistique de test:", test_result$statistic, "\n")
cat("P-value:", test_result$p.value, "\n")


# ==================================================
# ANALYZE RECOMBINATION INTERVAL LENGTH DISTRIBUTIONS
# ==================================================
# Create dataframe with recombination interval lengths for all methods
df_length_Recom=data.frame(type="SHapeit2_AA",Len=Recom.AA$END-Recom.AA$START+1, chr=Recom.AA$chr  )
df_length_Recom=rbind(df_length_Recom, data.frame(type="SHapeit2_EA",Len=Recom.EA$END-Recom.EA$START+1 , chr=Recom.EA$chr))
df_length_Recom=rbind(df_length_Recom, data.frame(type="Merlin_OE",Len=Recom.Merlin$END-Recom.Merlin$START+1 , chr=Recom.Merlin$chr ))

# Add Merlin AA-specific recombination lengths
tmp=Recom.Merlin[paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-") %in% AA_pair,]
df_length_Recom=rbind(df_length_Recom, data.frame(type="Merlin_AA",Len=tmp$END-tmp$START+1 , chr=tmp$chr ))

# Add Merlin EA-specific recombination lengths
tmp=Recom.Merlin[paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-") %in% EA_pair,]
df_length_Recom=rbind(df_length_Recom, data.frame(type="Merlin_EA",Len=tmp$END-tmp$START+1 , chr=tmp$chr ))

# Load OmniExpress (OE) Shapeit2 recombination data
Recom.OE=read.csv2(paste0(rpath,"/P1/paper3_Analyzes/Recom.OE.csv"))
# Filter to common parent-child pairs only
Recom.OE=Recom.OE[paste(Recom.OE$PARENT,Recom.OE$CHILD,sep="-") %in% common_indiv, ]

# Add OE recombination lengths to combined dataframe
df_length_Recom=rbind(df_length_Recom, data.frame(type="SHapeit2_OE",Len=Recom.OE$longueur , chr=Recom.OE$chr ))

# Create boxplot comparing recombination interval length distributions (log-transformed)
ggplot(df_length_Recom, aes(type, log10(Len), colour = type)) +
  geom_boxplot() +
  labs(x = "Method_data", y = "log-length of recombination intervals") +
  theme(legend.position = "none")


# ==================================================
# STATISTICAL COMPARISON OF RECOMBINATION INTERVAL LENGTHS BETWEEN AA AND EA
# ==================================================
# Load required library for Levene's test 
library(car)

# Merge AA recombinations with pedigree to assign family IDs
Recom.AA2= merge(Recom.AA,peds[,1:2], by.x="CHILD", by.y="V2", all.x=TRUE )
colnames(Recom.AA2)[ncol(Recom.AA2)]="FID"

# Calculate mean recombination interval length per family for AA population
LenByfamily <- Recom.AA2 %>%
  group_by(FID) %>%
  summarise(
    Mean_len = mean((longueur))  # Average interval length for family
  )
group1 <- LenByfamily$Mean_len

# Merge EA recombinations with pedigree to assign family IDs
Recom.EA2= merge(Recom.EA,peds[,1:2], by.x="CHILD", by.y="V2", all.x=TRUE )
colnames(Recom.EA2)[ncol(Recom.EA2)]="FID"

# Calculate mean recombination interval length per family for EA population
LenByfamily <- Recom.EA2 %>%
  group_by(FID) %>%
  summarise(
    Mean_len = mean((longueur))  # Average interval length for family
  )
group2 <- LenByfamily$Mean_len


# Test normality assumption for AA population (Shapiro-Wilk test)
# p-value > 0.05 indicates data is normally distributed
shapiro_test_group1 <- shapiro.test(group1)
# Test normality assumption for EA population (Shapiro-Wilk test)
shapiro_test_group2 <- shapiro.test(group2)

# Test equality of variances between populations (Levene's test)
# p-value > 0.05 indicates variances are equal
levene_test <- leveneTest(c(group1, group2) ~ factor(c(rep(1, length(group1)), rep(2, length(group2)))))

# Select appropriate statistical test based on normality and variance assumptions
if (shapiro_test_group1$p.value > 0.05 && shapiro_test_group2$p.value > 0.05) {
  # Data is normally distributed: use t-test
  if (levene_test$`Pr(>F)`[1] > 0.05) {
    # Equal variances: use standard t-test
    test_result <- t.test(group1, group2, var.equal = TRUE)
  } else {
    # Unequal variances: use Welch's t-test
    test_result <- t.test(group1, group2, var.equal = FALSE)
  }
} else {
  # Data is not normally distributed: use non-parametric Wilcoxon test
  test_result <- wilcox.test(group1, group2)
}

# Display statistical test results
list(
  shapiro_test_group1 = shapiro_test_group1,
  shapiro_test_group2 = shapiro_test_group2,
  levene_test = levene_test,
  test_result = test_result
)


# ==================================================
# COMPARE MEAN LOG-TRANSFORMED INTERVAL LENGTHS
# ==================================================
# Calculate mean log-transformed recombination interval length per AA family
LenByfamily <- Recom.AA2 %>%
  group_by(FID) %>%
  summarise(
    Mean_len = mean(log(longueur))  # Average log-transformed interval length
  )
group1 <- LenByfamily$Mean_len

# Merge EA recombinations with pedigree for family assignment
Recom.EA2= merge(Recom.EA,peds[,1:2], by.x="CHILD", by.y="V2", all.x=TRUE )
colnames(Recom.EA2)[ncol(Recom.EA2)]="FID"

# Calculate mean log-transformed recombination interval length per EA family
LenByfamily <- Recom.EA2 %>%
  group_by(FID) %>%
  summarise(
    Mean_len = mean(log(longueur))  # Average log-transformed interval length
  )
group2 <- LenByfamily$Mean_len


# Calculate mean values for both populations
mean_group1 <- mean((group1))
mean_group2 <- mean((group2))

# Calculate absolute difference between population means
absolute_difference <- abs(mean_group1 - mean_group2)
absolute_difference


# ==================================================
# ANALYZE RECOMBINATION OVERLAP BETWEEN METHODS
# ==================================================
# Load GenomicRanges library for efficient interval overlap calculations
library(GenomicRanges)

# Function to calculate percentage overlap of intervals between two datasets
# Parameters: df1 (query intervals), df2 (subject intervals)
# Returns: vector of overlap percentages for each interval in df1
prop_overlap <- function(df1, df2) {
  # Convert intervals to GRanges objects, using chr-CHILD as seqname for matching child-specific events
  gr1 <- GRanges(seqnames = paste0(df1$chr,"-",df1$CHILD), ranges = IRanges(start = df1$START, end = df1$END))
  gr2 <- GRanges(seqnames = paste0(df2$chr,"-",df2$CHILD), ranges = IRanges(start = df2$START, end = df2$END))
  
  # Find overlapping intervals between the two GRanges objects
  ov <- findOverlaps(gr1, gr2)
  # Initialize overlap percentage vector (0 for no overlap)
  overlap =rep(0, nrow(df1))
  # Calculate overlap percentage as (intersection length / interval1 length) * 100
  overlap[queryHits(ov)] <- width(pintersect(gr1[queryHits(ov)], gr2[subjectHits(ov)])) / width(gr1[queryHits(ov)])
  
  # Convert to percentage
  return(overlap*100)
}

# Calculate overlap percentages between Merlin and AA Shapeit2 recombinations
Recom.Merlin$"ov_AA" <- prop_overlap(Recom.Merlin, Recom.AA)
# Calculate overlap percentages between Merlin and EA Shapeit2 recombinations
Recom.Merlin$"ov_EA" <- prop_overlap(Recom.Merlin, Recom.AA)

# Calculate overlap percentages between AA Shapeit2 and Merlin recombinations
Recom.AA$"ov_Mr" <- prop_overlap(Recom.AA, Recom.Merlin)

# Calculate overlap percentages between EA Shapeit2 and Merlin recombinations
Recom.EA$"ov_Mr" <- prop_overlap(Recom.EA, Recom.Merlin)


# Calculate concordance rate: percentage of Shapeit2-AA events with any overlap to Merlin events
table(Recom.AA$ov_Mr>0)[2]*100/nrow(Recom.AA)

# Calculate concordance rate: percentage of Shapeit2-EA events with any overlap to Merlin events
table(Recom.EA$ov_Mr>0)[2]*100/nrow(Recom.EA)


# ==================================================
# CALCULATE CONCORDANCE RATES BY CHROMOSOME
# ==================================================
# Filter AA Shapeit2 recombinations with overlap to Merlin events
subset_data <- Recom.AA[Recom.AA$ov_Mr > 0, ]
# Count concordant events per chromosome for AA
result <- as.data.frame(table(subset_data$chr))
colnames(result)=c("chr", "shapeit2AA-MerlinOE")
Concordance_df=result

# Filter EA Shapeit2 recombinations with overlap to Merlin events
subset_data <- Recom.EA[Recom.EA$ov_Mr > 0, ]
# Count concordant events per chromosome for EA
result <- as.data.frame(table(subset_data$chr))
colnames(result)=c("chr", "shapeit2EA-MerlinOE")
Concordance_df=merge(Concordance_df, result, by="chr")

# Normalize AA concordance to per-meiosis rate for AA families
Concordance_df$`shapeit2AA-Merlin`=Concordance_df$`shapeit2AA-MerlinOE`/sum(AA_pair  %in% paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-"))
# Normalize EA concordance to per-meiosis rate for EA families
Concordance_df$`shapeit2EA-Merlin`=Concordance_df$`shapeit2EA-MerlinOE`/sum(EA_pair  %in% paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-"))

# Round concordance rates to 2 decimal places for readability
Concordance_df$`shapeit2AA-Merlin`=round(Concordance_df$`shapeit2AA-Merlin`, 2)
Concordance_df$`shapeit2EA-Merlin`=round(Concordance_df$`shapeit2EA-Merlin`, 2)
