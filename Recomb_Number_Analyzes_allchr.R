# ============================================================================
# RECOMBINATION NUMBER ANALYSIS ACROSS ALL CHROMOSOMES
# ============================================================================
# This script analyzes recombination events detected by multiple methods:
# - Shapeit2 on WGS (Whole Genome Sequencing) data
# - Shapeit2 on OmniExpress data
# - Merlin on OmniExpress data
# The analysis includes comparative statistics and visualizations of recombination rates

# Clear workspace and load required libraries
rm(list=ls())
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggforce)

# ============================================================================
# SECTION 1: DATA INPUT SETUP
# ============================================================================
# Define base path to project results directory
rpath="~/links/projects/rrg-girardsi/genealogy_sims/results/Samir/"

# Load pedigree file (family information)
ped=read.table(paste0(rpath,"/P1/Candiate_recombinations/WGS.ped"))


# Construct file paths for all 22 chromosomes
chemin2Recomb=paste0(rpath,"/SFTP hub/PHD -P1/Data/")
chr=c(1:22)
# Paths to Shapeit2 WGS recombination files 
cheminSC = paste0(chemin2Recomb, "chr",chr,".pass.clean.recombinations")
cheminSC.3innuc = paste0(chemin2Recomb, "chr",chr,".pass.clean.3innuc.recombinations")

# Chromosome 9 requires special handling because data is split into 9p and 9q arms
# The following code shows how to concatenate these files (9p (29) et 9q (9)) if needed
#Recom.SC.chr9= read.delim2(paste0(chemin2Recomb, "chr9p",".pass.clean.recombinations"),
#                           header = TRUE, sep = "")
#Recom.SC.chr9=rbind(Recom.SC.chr9,
#                    read.delim2(paste0(chemin2Recomb, "chr9q",".pass.clean.recombinations"), header = TRUE, sep = ""))
#write.table(Recom.SC.chr9, paste0(chemin2Recomb, "chr9",".pass.clean.recombinations"), quote = F , row.names = F)

#Recom.SC.chr9= read.delim2(paste0(chemin2Recomb, "chr9p",".pass.clean.3innuc.recombinations"),
#                           header = TRUE, sep = "")
#Recom.SC.chr9=rbind(Recom.SC.chr9,
#                    read.delim2(paste0(chemin2Recomb, "chr9q",".pass.clean.3innuc.recombinations"), header = TRUE, sep = ""))
#write.table(Recom.SC.chr9, paste0(chemin2Recomb, "chr9",".pass.clean.3innuc.recombinations"), quote = F , row.names = F)


# ============================================================================
# SECTION 2: LOAD SHAPEIT2 WGS DATA (Recom.SC)
# ============================================================================
# Load chromosome 1 data from both pass.clean and 3innuc files
Recom.SC= read.delim2(cheminSC[1], header = TRUE, sep = "")
Recom.SC=rbind(Recom.SC,read.delim2(cheminSC.3innuc[1], header = TRUE, sep = ""))
Recom.SC=data.frame("chr"=chr[1], Recom.SC )

# Load recombination data from remaining chromosomes (2-22)
for (j in 2:length(chr)) {
  tmp=read.delim2(cheminSC[j], header = TRUE, sep = "")
  tmp=rbind(tmp,read.delim2(cheminSC.3innuc[j], header = TRUE, sep = ""))
  #tmp=tmp[tmp$PROB_RECOMBINATION >=  0.5,]
  Recom.SC = rbind(Recom.SC,data.frame("chr"=chr[j], tmp))
  #print(nrow(Recom.SC))
}

# ============================================================================
# SECTION 2a: POST-PROCESS SHAPEIT2 WGS DATA
# ============================================================================
# Calculate recombination interval length 
Recom.SC$"longueur"=Recom.SC$END-Recom.SC$START+1
# Identify unique parent-child pairs in the dataset
CS_pair=(unique(paste(Recom.SC$PARENT,Recom.SC$CHILD,sep="-")))

# Filter recombination events: keep only those with probability >= 0.5
Recom.SC=Recom.SC[Recom.SC$PROB_RECOMBINATION >=  0.5,]

# ============================================================================
# SECTION 3: LOAD SHAPEIT2 OMNIEXPRESS DATA (Recom.OE)
# ============================================================================
# Define paths to OmniExpress recombination data files
path2RecombOE=paste0(rpath,"/SFTP hub/PHD -P1/Data/")
cheminOE = paste0(path2RecombOE, "chr",chr,".OmniExpress.clean.recombinations")
cheminOE.3innuc = paste0(path2RecombOE,"chr",chr,".OmniExpress.clean.3innuc.recombinations")

# Load chromosome 1 data from both files
Recom.OE= read.delim2(cheminOE[1], header = TRUE, sep = "")
Recom.OE=rbind(Recom.OE,read.delim2(cheminOE.3innuc[1], header = TRUE, sep = "") )
Recom.OE=data.frame("chr"=chr[1], Recom.OE )

# Load remaining chromosomes (2-22)
for (j in 2:length(chr)) {
  tmp=read.delim2(cheminOE[j], header = TRUE, sep = "")
  tmp=rbind(tmp, read.delim2(cheminOE.3innuc[j], header = TRUE, sep = "") )
  # Optional: filter by probability (commented out)
  #tmp=tmp[tmp$PROB_RECOMBINATION >=  0.5,]
  
  # Append to main dataset
  Recom.OE = rbind(Recom.OE,data.frame("chr"=chr[j], tmp ))
  #print(nrow(Recom.OE))
}

# ============================================================================
# SECTION 3a: POST-PROCESS SHAPEIT2 OMNIEXPRESS DATA
# ============================================================================
# Calculate recombination interval length for OE data
Recom.OE$"longueur"=Recom.OE$END-Recom.OE$START+1
# Identify unique parent-child pairs in OE dataset
OE_pair=(unique(paste(Recom.OE$PARENT,Recom.OE$CHILD,sep="-")))
# Filter by recombination probability threshold >= 0.5
Recom.OE=Recom.OE[Recom.OE$PROB_RECOMBINATION >= 0.5,]

# ============================================================================
# SECTION 4: LOAD MERLIN DATA (Recom.Merlin)
# ============================================================================
# Create chromosome labels with "chr" prefix for Merlin files
chrNum = paste0("chr", chr)

# NOTE: Merlin recombination files are produced by:
# - "/lustre09/project/6033529/genealogy_sims/results/Samir/P1/paper3_Analyzes/Recomb_Merlin_UncutFam_allchr.R"
# - "/lustre09/project/6033529/genealogy_sims/results/Samir/P1/paper3_Analyzes/Recomb_Merlin_allchr.R"
# Files include all chromosomes including chr9


# Define paths to Merlin recombination output files
CheminRM=paste0(rpath, "/P1/paper3_Analyzes/Merlin_Recomb/", "chr", chr, ".merlin.recombinaisons")
CheminRM.3innuc=paste0(rpath, "/P1/paper3_Analyzes/Merlin_Recomb/", "chr", chr, ".merlin.3innuc.recombinaisons")

# Load chromosome 1 Merlin data
Recom.Merlin = read.delim2(CheminRM[1], header = TRUE, sep = "") 
Recom.Merlin=rbind(Recom.Merlin, read.table(CheminRM.3innuc[1], header = T) )
# Add chromosome identifier
Recom.Merlin$"chr"=chr[1]

# Load remaining chromosomes (2-22)
for (j in 2:length(chr)) {
  tmp=read.delim2(CheminRM[j], header = TRUE, sep = "")
  tmp=rbind(tmp, read.table(CheminRM.3innuc[j], header = T) )
  tmp$"chr"=chr[j]
  # Append to main dataset
  Recom.Merlin=rbind(Recom.Merlin, tmp )
}

# ============================================================================
# SECTION 4a: FILTER AND CLEAN MERLIN DATA
# ============================================================================
# Count unique parent-child pairs in Merlin dataset (expected: ~61 families)
length(unique(paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-")))

# Identify families that have recombination events in the WGS data
famsSC=(unique(ped$V1[ped$V2 %in% Recom.SC$CHILD]))
# Restrict Merlin data to families with WGS recombination events
Recom.Merlin=  Recom.Merlin[Recom.Merlin$CHILD %in% ped$V2[ped$V1 %in% famsSC], ]
row.names(Recom.Merlin) = NULL

# Remove records with missing 
Recom.Merlin=Recom.Merlin[!is.na(Recom.Merlin$START),]
Recom.Merlin=Recom.Merlin[!is.na(Recom.Merlin$END),]
# Remove duplicate recombination records
Recom.Merlin=Recom.Merlin[!duplicated(Recom.Merlin),]

# ============================================================================
# SECTION 5: IDENTIFY COMMON PARENT-CHILD PAIRS ACROSS METHODS
# ============================================================================
# Find parent-child pairs that are present in all three datasets
# This ensures fair comparison between methods
common_indiv <- unique(Reduce(intersect, list(paste(Recom.OE$PARENT,Recom.OE$CHILD,sep="-"),
                                               paste(Recom.SC$PARENT,Recom.SC$CHILD,sep="-"),
                                               paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-"))))

# ============================================================================
# SECTION 5a: FILTER ALL DATASETS TO COMMON PAIRS
# ============================================================================
# Restrict Merlin data to common parent-child pairs
Recom.Merlin=Recom.Merlin[paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-") %in% common_indiv , ]
# Restrict OmniExpress data to common parent-child pairs
Recom.OE=Recom.OE[paste(Recom.OE$PARENT,Recom.OE$CHILD,sep="-") %in% common_indiv , ]
# Restrict WGS data to common parent-child pairs
Recom.SC=Recom.SC[paste(Recom.SC$PARENT,Recom.SC$CHILD,sep="-") %in% common_indiv , ]

# ============================================================================
# SECTION 6: EXPECTED RECOMBINATION NUMBERS
# ============================================================================
# Create reference table with expected recombination counts per chromosome
# Based on known recombination rates scaled to the number of common parent-child pairs
# Values represent expected recombination events per 100 meioses, scaled accordingly
NR_expected_common= data.frame(chr=chr,
                               Expected=round(c(284.8,274.6,233.0,216.4,213.0, 197.5, 190.9,
                                                171.5, 168.5, 175.5, 165.9, 168.1,124.1, 124.9, 127.0,
                                                133.5, 138.8, 121.3, 102.6, 102.6, 63.7, 69.5)*length(common_indiv)/100 , 0))

# ============================================================================
# SECTION 7: VISUALIZATION 1 - RECOMBINATION COUNT BY CHROMOSOME
# ============================================================================
# Create contingency table: number of recombination events per chromosome per method
df=merge(as.data.frame(table(Recom.Merlin$chr)), as.data.frame(table(Recom.OE$chr)), by="Var1")
df=merge(df, as.data.frame(table(Recom.SC$chr)), by="Var1")
colnames(df)[1:4]=c("chr", "Merlin_OE", "Shapeit2_OE", "Shapeit2_WGS")
df=merge(df, NR_expected_common, by="chr")
# Convert dataframe from wide to long format for ggplot2
df_melt <- df %>% pivot_longer(cols = -chr)

# Plot 1: Line plot of recombination event counts per chromosome across methods
ggplot(df_melt, aes(chr, value, group=name, colour=name)) +
  geom_line() +
  labs(x="Chromosome number", y="Number of recombination events", colour=" ") 
  # Optional: add faceted zoom on specific range
  #facet_zoom(y=(value <= 900), zoom.size = 1.2, show.area =F)


# Plot 2: Recombination rate per meiosis 
ggplot(df_melt, aes(chr, value/length(common_indiv), group = name, colour = name)) +
  geom_line() +
  # Optional: add dotted line for missing data
  #geom_line(data = subset(df_melt, !is.na(value)), linetype = "dotted", size = 0.5) +
  geom_point(size = 0.5) +
  labs(x = "Chromosome number",
       y = "Number of recombination events per meiosis",
       colour = " ") +
  # Add secondary y-axis showing total counts (right side)
  scale_y_continuous(    sec.axis = sec_axis(~ . * 295, name = "Number of recombination events")  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_minimal()


# ============================================================================
# SECTION 8: VISUALIZATION 2 - RECOMBINATION INTERVAL LENGTH DISTRIBUTION
# ============================================================================
# Create dataframe with recombination interval lengths for all methods
df_length_Recom=data.frame(type="SHapeit2_OE",Len=Recom.OE$END-Recom.OE$START+1, chr=Recom.OE$chr  )
df_length_Recom=rbind(df_length_Recom, data.frame(type="SHapeit2_WGS",Len=Recom.SC$END-Recom.SC$START+1 , chr=Recom.SC$chr))
df_length_Recom=rbind(df_length_Recom, data.frame(type="Merlin_OE",Len=Recom.Merlin$END-Recom.Merlin$START+1 , chr=Recom.Merlin$chr ))


# Plot 3: Box plot of log-transformed recombination interval lengths across methods
ggplot(df_length_Recom, aes(type, log10(Len), colour = type)) +
  geom_boxplot() +
  # Alternative: use ggplot2's built-in log scale
  #scale_y_log10()+
  labs(x = "Method_data", y = "log-length of recombination intervals") +
  theme(legend.position = "none")

# Plot 4: Box plot of log-transformed lengths faceted by chromosome
ggplot(df_length_Recom, aes(type, log10(Len), colour = type)) +
  geom_boxplot() +
  #scale_y_log10() +
  labs(x = "", y = "log-length of recombination intervals", colour = "Method-Data") +
  facet_wrap(~ chr, ncol = 6)+
  theme(axis.text.x = element_blank())


# ============================================================================
# SECTION 8a: SUMMARY STATISTICS FOR RECOMBINATION INTERVALS
# ============================================================================
# Calculate descriptive statistics (count, min, quartiles, median, mean, max) for each method
summary_stats <- df_length_Recom %>%
  group_by(type) %>%
  summarise(N = n(),
            min = min(Len),
            Q1 = quantile(Len, 0.25),
            median = median(Len),
            mean = mean(Len),
            Q3 = quantile(Len, 0.75),
            max = max(Len))

# Print summary statistics table
print(summary_stats)

# ============================================================================
# SECTION 9: SUMMARY TABLE - TOTAL COUNTS AND RATES
# ============================================================================
# Create summary table with total recombination counts and normalized rates
recap1=data.frame("Number of recombination"=colSums(df[2:5]), 
                  "Number of recombination per meiosis"= round(colSums(df[2:5])/length(common_indiv),2) )


# ============================================================================
# SECTION 10: VISUALIZATION 3 - RECOMBINATION SEGMENT DISPLAY
# ============================================================================
# Sort each dataset by recombination start position (for plotting)
Recom.Merlin=Recom.Merlin[order(Recom.Merlin$START),]
Recom.SC=Recom.SC[order(Recom.SC$START),]
Recom.OE=Recom.OE[order(Recom.OE$START),]

# Combine selected columns from all three datasets
Recom.All=rbind(Recom.SC[c(1:5)], Recom.OE[c(1:5)], Recom.Merlin[c(6,1:4)])
# Add method identifier column: SC=Shapeit2 WGS, OE=Shapeit2 OE, MR=Merlin
Recom.All$Method= c(rep("SC", nrow(Recom.SC)), rep("OE", nrow(Recom.OE)), rep("MR", nrow(Recom.Merlin)))
Recom.All=Recom.All[order(Recom.All$START),]
# Initialize Y-axis position column (for vertical stacking visualization)
Recom.All$Y=NA

# Assign random Y-coordinates for each recombination event within method and chromosome
# This prevents overlap when plotting segments on the same chromosome
for (chrom in unique(Recom.All$chr)) {
  for (method in unique(Recom.All$Method)) {
    # Count number of recombination events for this chromosome-method combination
    nr=nrow(Recom.All[Recom.All$chr==chrom & Recom.All$Method==method,])
    # Assign random Y positions from 1 to nr
    Recom.All$Y[Recom.All$chr==chrom & Recom.All$Method==method]=sample(1:nr, nr)
  }
}

# ============================================================================
# SECTION 10a: PLOT RECOMBINATION SEGMENTS BY METHOD
# ============================================================================
# Create plot showing recombination segments (Example  : for chromosome 1 )
# Each segment represents a recombination interval detected by each method
p=ggplot(Recom.All[Recom.All$chr %in% c(1) ,], aes(x = 1, y = 1, color = factor(Method))) +
  # Draw horizontal line segments for each recombination event
  geom_segment(aes(x = START, xend = END, 
                   y = Y, 
                   yend = Y), size = 1)+
  # Separate plots for each method for clarity
  facet_wrap(~ Method, ncol = 1)+
  theme(legend.position = "none")+
  ylab("Y Axis Label")+
  xlab(NULL)+
  labs(title="Recombination distribution by Method")

# Display the plot
p

# ============================================================================
# SECTION 11: VISUALIZATION 4 - RECOMBINATION RATE BY GENOMIC POSITION
# ============================================================================
# Calculate recombination frequency in 1 Mb windows across the genome
# Create the output dataframe with counts per Mb window
result <- Recom.All %>%
  # Convert genomic position to megabases (1e6 base pairs)
  mutate(position = START / 1e6) %>%
  # Group by chromosome and method
  group_by(chr, Method, position = floor(position)+1) %>%
  # Count number of recombination events in each window
  summarise(rate = n(), .groups = 'drop') %>%
  # Rename position column 
  rename(position_Mbp = position)

# Plot 5: Line plot of recombination rate per Mb for chromosome 1 by method
ggplot(result[result$chr %in% c(1),], aes(x = position_Mbp  , y = rate, color = Method)) +
  geom_line() +
  facet_wrap(~ chr, scales = "free_x", ncol = 3) +
  labs(title = "",
       x = "Position (Mb)",
       y = "Number of recombination by Mb") +
  theme_minimal() +
  theme(strip.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  # Add horizontal line for expected recombination rate (red dashed line)
  geom_hline(yintercept = length(common_indiv)/100, linetype = "dashed", color = "red") +
  # Add text label showing expected value
  geom_text(aes(x = -05, y = 05, label = paste0("Expected=",length(common_indiv)/100)), vjust = -0.050, hjust = -0.050, color = "black", size = 3)



# ============================================================================
# SECTION 12: VISUALIZATION 5 - SMOOTHED RECOMBINATION RATE (LOESS METHOD)
# ============================================================================
# Plot 6: Smoothed recombination rate for chromosome 9
# Uses LOESS smoothing to reveal trends in recombination rate across genomic position
ggplot(result[result$chr == 9,], aes(x = position_Mbp, y = rate, color = Method)) +
  # Draw raw data line 
  geom_line(alpha = 0.3) +
  # Overlay LOESS smoothed line (span=0.1 controls smoothing intensity)
  stat_smooth(method = "loess", span = 0.1, se = FALSE, size = 1) +
  # Add horizontal line for expected recombination rate
  geom_hline(yintercept = length(common_indiv)/100, linetype = "dashed", color = "red") +
  # Add text annotation for expected value
  annotate("text",
           x = -5, y = 5,
           label = paste0("Expected=", length(common_indiv)/100),
           vjust = -Inf, hjust = -Inf,
           color = "black", size = 3) +
  # Facet by chromosome 
  facet_wrap(~ chr, scales = "free_x", ncol = 6) +
  scale_y_continuous(name = "Number of recombination by Mb") +
  theme_minimal()

# ============================================================================
# SECTION 13: VISUALIZATION 6 - SMOOTHED RECOMBINATION RATE (ALL CHROMOSOMES)
# ============================================================================
# Plot 7: Smoothed recombination rates for all chromosomes
# Shows recombination hotspots and coldspots across the entire genome
ggplot(result, aes(x = position_Mbp, y = rate, color = Method)) +
  geom_line(alpha = 0.3) +
  # Overlay LOESS smoothed line
  stat_smooth(method = "loess", span = 0.1, se = FALSE, size = 1) +
  # Add horizontal line for expected recombination rate
  geom_hline(yintercept = length(common_indiv)/100, linetype = "dashed", color = "red") +
  # Add text annotation for expected value
  annotate("text",
           x = -5, y = 5,
           label = paste0("Expected=", length(common_indiv)/100),
           vjust = -Inf, hjust = -Inf,
           color = "black", size = 3) +
  # Separate panels for each chromosome (6 columns of plots)
  facet_wrap(~ chr, scales = "free_x", ncol = 6) +
  scale_y_continuous(name = "Number of recombination by Mb") +
  theme_minimal()
