# Clear all objects from the workspace
rm(list=ls())
# Load required libraries 
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggforce)

#Primary Goal : Compare recombination detection methods across the genome by analyzing recombination events identified using the three different approaches

# Path to the data directory in rorqual storage
rpath="~/links/projects/rrg-girardsi/genealogy_sims/results/Samir/"

# Load the pedigree file containing family information
ped=read.table(paste0(rpath,"/P1/Candiate_recombinations/WGS.ped"))

# Load recombination data from Shapeit2-duoHMM on WGS 
chemin2Recomb=paste0(rpath,"/SFTP hub/PHD -P1/Data/")
chr=c(1:9,10:22)
cheminSC = paste0(chemin2Recomb, "chr",chr,".pass.clean.recombinations")
cheminSC.3innuc = paste0(chemin2Recomb, "chr",chr,".pass.clean.3innuc.recombinations")

Recom.SC= read.delim2(cheminSC[1], header = TRUE, sep = "")
Recom.SC=rbind(Recom.SC,read.delim2(cheminSC.3innuc[1], header = TRUE, sep = ""))
# Add chromosome column to the data
Recom.SC=data.frame("chr"=chr[1], Recom.SC )

# Loop through remaining chromosomes and combine all recombination data
for (j in 2:length(chr)) {
  tmp=read.delim2(cheminSC[j], header = TRUE, sep = "")
  tmp=rbind(tmp,read.delim2(cheminSC.3innuc[j], header = TRUE, sep = ""))
  #tmp=tmp[tmp$PROB_RECOMBINATION >=  0.5,]
  Recom.SC = rbind(Recom.SC,data.frame("chr"=chr[j], tmp))
  #print(nrow(Recom.SC))
}

# Calculate recombination segment length (END - START + 1)
Recom.SC$"longueur"=Recom.SC$END-Recom.SC$START+1
# Extract unique parent-child pairs from Shapeit2-duoHMM on WGS data
CS_pair=(unique(paste(Recom.SC$PARENT,Recom.SC$CHILD,sep="-")))
# Filter recombinations to keep only those with probability >= 0.5
Recom.SC=Recom.SC[Recom.SC$PROB_RECOMBINATION >=  0.5,]

# Define paths for OmniExpress microarray recombination data
path2RecombOE=paste0(rpath,"/SFTP hub/PHD -P1/Data/")
cheminOE = paste0(path2RecombOE, "chr",chr,".OmniExpress.clean.recombinations")
cheminOE.3innuc = paste0(path2RecombOE,"chr",chr,".OmniExpress.clean.3innuc.recombinations")

# Load recombination data from OmniExpress microarray 
Recom.OE= read.delim2(cheminOE[1], header = TRUE, sep = "")
Recom.OE=rbind(Recom.OE,read.delim2(cheminOE.3innuc[1], header = TRUE, sep = "") )
Recom.OE=data.frame("chr"=chr[1], Recom.OE )

# Loop through remaining chromosomes and combine all OmniExpress recombination data
for (j in 2:length(chr)) {
  tmp=read.delim2(cheminOE[j], header = TRUE, sep = "")
  tmp=rbind(tmp, read.delim2(cheminOE.3innuc[j], header = TRUE, sep = "") )
  #tmp=tmp[tmp$PROB_RECOMBINATION >=  0.5,]
  
  Recom.OE = rbind(Recom.OE,data.frame("chr"=chr[j], tmp ))
  #print(nrow(Recom.OE))
}
# Calculate recombination segment length for OmniExpress data
Recom.OE$"longueur"=Recom.OE$END-Recom.OE$START+1
# Extract unique parent-child pairs from OmniExpress data
OE_pair=(unique(paste(Recom.OE$PARENT,Recom.OE$CHILD,sep="-")))
# Filter recombinations to keep only those with probability >= 0.5
Recom.OE=Recom.OE[Recom.OE$PROB_RECOMBINATION >= 0.5,]

# Prepare chromosome identifiers for Merlin data files
chrNum = paste0("chr", chr)

# Paths for Merlin recombination data (affects subset of families)
CheminRM=paste0(rpath,"/P1/paper3_Analyzes/Merlin_Recomb/", 
                "chr", chr, ".merlin.recombinaisons")
CheminRM.3innuc=paste0(rpath,"/P1/paper3_Analyzes/Merlin_Recomb/", 
                       "chr", chr, ".merlin.3innuc.recombinaisons")

# Load Merlin recombination data for the first chromosome
Recom.Merlin = read.delim2(CheminRM[1], header = TRUE, sep = "") 
Recom.Merlin=rbind(Recom.Merlin, read.table(CheminRM.3innuc[1], header = T) )
Recom.Merlin$"chr"=chr[1]

# Loop through remaining chromosomes and combine all Merlin recombination data
for (j in 2:length(chr)) {
  tmp=read.delim2(CheminRM[j], header = TRUE, sep = "")
  tmp=rbind(tmp, read.table(CheminRM.3innuc[j], header = T) )
  tmp$"chr"=chr[j]
  Recom.Merlin=rbind(Recom.Merlin, tmp )
}
#length(unique(paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-")))

# Extract families that have children with recombinations in Shapeit2-duoHMM on WGS data
famsSC=(unique(ped$V1[ped$V2 %in% Recom.SC$CHILD]))
# Filter Merlin data to only include children from families present in Shapeit2-duoHMM on WGS data
Recom.Merlin=  Recom.Merlin[Recom.Merlin$CHILD %in% ped$V2[ped$V1 %in% famsSC], ]
# Reset row names
row.names(Recom.Merlin) = NULL
# Remove rows with missing START positions
Recom.Merlin=Recom.Merlin[!is.na(Recom.Merlin$START),]
# Remove rows with missing END positions
Recom.Merlin=Recom.Merlin[!is.na(Recom.Merlin$END),]
# Remove duplicate recombination records
Recom.Merlin=Recom.Merlin[!duplicated(Recom.Merlin),]



#########################################
# Analyze recombination overlap across methods
#########################################

# Identify parent-child pairs present in all three recombination datasets
common_indiv <- unique(Reduce(intersect, list(paste(Recom.OE$PARENT,Recom.OE$CHILD,sep="-"),
                                               paste(Recom.SC$PARENT,Recom.SC$CHILD,sep="-"),
                                               paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-"))))

# Filter all three datasets to keep only parent-child pairs common across all methods
Recom.Merlin=Recom.Merlin[paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-") %in% common_indiv , ]
Recom.OE=Recom.OE[paste(Recom.OE$PARENT,Recom.OE$CHILD,sep="-") %in% common_indiv , ]
Recom.SC=Recom.SC[paste(Recom.SC$PARENT,Recom.SC$CHILD,sep="-") %in% common_indiv , ]

# Expected number of recombinations per chromosome scaled by the number of common individuals
NR_expected_common= data.frame(chr=chr,
                               Expected=round(c(284.8,274.6,233.0,216.4,213.0, 197.5, 190.9,
                                                171.5, 168.5, 175.5, 165.9, 168.1,124.1, 124.9, 127.0,
                                                133.5, 138.8, 121.3, 102.6, 102.6, 63.7, 69.5)*length(common_indiv)/100 , 0))


# Calculate recombination overlap between datasets using genomic ranges

# Load GenomicRanges library for handling genomic intervals
library(GenomicRanges)

# Function to calculate the overlap percentage between recombination intervals in two datasets
# df1: primary dataset, df2: secondary dataset for comparison
# Returns: vector of overlap percentages (0-100) for each record in df1
prop_overlap <- function(df1, df2) {
  # Convert df1 recombinations to GRanges objects with unique identifiers per child
  gr1 <- GRanges(seqnames = paste0(df1$chr,"-",df1$CHILD), ranges = IRanges(start = df1$START, end = df1$END))
  # Convert df2 recombinations to GRanges objects with unique identifiers per child
  gr2 <- GRanges(seqnames = paste0(df2$chr,"-",df2$CHILD), ranges = IRanges(start = df2$START, end = df2$END))
  
  # Find overlapping intervals between gr1 and gr2
  ov <- findOverlaps(gr1, gr2)
  # Initialize overlap vector with zeros (no overlap)
  overlap =rep(0, nrow(df1))
  # Calculate overlap percentage as intersection width divided by gr1 width
  overlap[queryHits(ov)] <- width(pintersect(gr1[queryHits(ov)], gr2[subjectHits(ov)])) / width(gr1[queryHits(ov)])
  
  # Convert overlap fraction to percentage and return
  return(overlap*100)
}

# Calculate overlap between Merlin and Shapeit2-duoHMM on OmniExpress recombinations
Recom.Merlin$"ov_OE" <- prop_overlap(Recom.Merlin, Recom.OE)
# Calculate overlap between Merlin recombinations and Shapeit2-duoHMM on WGS recombinations
Recom.Merlin$"ov_SC" <- prop_overlap(Recom.Merlin, Recom.SC)

# Calculate overlap between OmniExpress and Merlin recombinations
Recom.OE$"ov_Mr" <- prop_overlap(Recom.OE, Recom.Merlin)
# Calculate overlap between Shapeit2-duoHMM on  OmniExpress and WGS recombinations
Recom.OE$"ov_SC" <- prop_overlap(Recom.OE, Recom.SC)

# Calculate overlap between Shapeit2-duoHMM  recombinations on WGS and OmniExpress 
Recom.SC$"ov_OE" <- prop_overlap(Recom.SC, Recom.OE)
# Calculate overlap between Shapeit2-duoHMM on WGS recombinations  and Merlin recombinations
Recom.SC$"ov_Mr" <- prop_overlap(Recom.SC, Recom.Merlin)


# Calculate Shapeit2-duoHMM on WGS concordance 
#table(Recom.Merlin$ov_SC>0)[2]*100/nrow(Recom.SC)
table(Recom.SC$ov_Mr>0)[2]*100/nrow(Recom.SC)

# Calculate OmniExpress concordance 
#table(Recom.Merlin$ov_SC>0)[2]*100/nrow(Recom.OE)
table(Recom.OE$ov_Mr>0)[2]*100/nrow(Recom.OE)
# Save OmniExpress recombination data with overlap values to CSV
write.csv2(Recom.OE, file=paste0(rpath,"/P1/paper3_Analyzes/Recom.OE.csv"), row.names = FALSE)

# Calculate  Shapeit2-duoHMM on WGS & OE  concordance 
nrow(Recom.SC[Recom.SC$ov_OE>0 & Recom.SC$ov_Mr>0,])*100/table(Recom.SC$ov_OE>0)[2]
write.csv2(Recom.SC, file=paste0(rpath,"/P1/paper3_Analyzes/Recom.SC.csv"), row.names = FALSE)


# Generate concordance table by chromosome
subset_data <- Recom.SC[Recom.SC$ov_Mr > 0, ]
result <- as.data.frame(table(subset_data$chr))
colnames(result)=c("chr", "shapeit2WGS-MerlinOE")
Concordance_df=result

subset_data <- Recom.OE[Recom.OE$ov_Mr > 0, ]
result <- as.data.frame(table(subset_data$chr))
colnames(result)=c("chr", "shapeit2OE-MerlinOE")
Concordance_df=merge(Concordance_df, result, by="chr")

# Shapeit2-duoHMM on WGS recombinations with OmniExpress overlap per chromosome
subset_data <- Recom.SC[Recom.SC$ov_OE > 0, ]
result <- as.data.frame(table(subset_data$chr))
colnames(result)=c("chr", "shapeit2SC-shapeit2OE")
# Merge with the concordance table
Concordance_df=merge(Concordance_df, result, by="chr")

# Count Shapeit2-duoHMM on WGS and OmniExpress recombinations overlapping Merlin per chromosome
subset_data <- Recom.SC[Recom.SC$ov_OE>0 & Recom.SC$ov_Mr>0,]
result <- as.data.frame(table(subset_data$chr))
colnames(result)=c("chr", "shapeit2WGS-shapeit2OE-MerlinOE")
# Merge with the concordance table
Concordance_df=merge(Concordance_df, result, by="chr")
