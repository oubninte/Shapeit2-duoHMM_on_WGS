# Clear workspace
rm(list=ls())

# Load required libraries for data manipulation and genomic analysis
library(dplyr)
library(stringr)
library(GenomicRanges)
library(BRGenomics)


# This program identifies recombinations in uncut families using Merlin outputs
# It locates recombinations with intervals based on parental heterozygote markers


# Define paths to input data directories
path_rep="/lustre09/project/6033529/genealogy_sims/results/Samir/SFTP hub/PHD -P1/Data/DataMerlin/hapAD/"
# Construct paths to genetic map files for all chromosomes (1-22)
chemin_map = c( paste0(path_rep, "/chr", c(1:9),"/merlin_map.0",c(1:9)),
            paste0(path_rep, "/chr", c(10:22), "/merlin_map.", c(10:22)))

# Path to Merlin output flow files containing haplotype information
cheminMF = paste0("/lustre09/project/6033529/genealogy_sims/results/Samir/SFTP hub/PHD -P1/Data/DataMerlin/HapAvantDecoupe/Merlin_hap_31wiped/chr", 
c(1:9,10:22),"/merlin.flow")
chrNum = paste0("chr", c(1:9,10:22))
# Paths to heterozygosity status files for each chromosome (stored as RData)
cheminHet=paste0(path_rep, "het.",c(1:9,10:22),".RData")

# Load pedigree information 
ped=read.table("/lustre09/project/6033529/genealogy_sims/results/Samir/P1/Candiate_recombinations/WGS.ped")

# Example  : Process chromosome 9 (loop structured to allow processing of multiple chromosomes c(1:22))
for (j in c(9)) {
  
  # Load heterozygosity status data for the current chromosome
  heterozygous=get(load(cheminHet[j]))
  # Rename the first four columns to standard family-related identifiers
  colnames(heterozygous)[1:4]=c("Family", "Indiv","Father","Mother")
  
  # Load genetic map file containing marker information and genetic positions
  merlin_map = read.delim2(chemin_map[j], header = T, sep = "")
  
  # Validate that the number of markers in the map matches the number of markers in heterozygosity data
  cat("Est ce que le nombre de marqueurs correspond avec le fichier map :chr", j)
  if (nrow(merlin_map)==ncol(heterozygous)-6) {
    print("Oui")
  } else {print("Non")
    break
    }
  
  # Extract and add physical position 
  merlin_map = cbind(merlin_map, "PHY.POS"=as.numeric(str_split_fixed(merlin_map$MARKER, ":", 3)[,2]))
  
  # Read and process Merlin flow file to prepare for recombination inference
  chr = read.delim2(cheminMF[j], header = FALSE, sep = "")
  chr = cbind("Family" = chr$V2, chr)
  chr$Family[chr$Family %in% c("(FOUNDER)", "(MATERNAL)", "(PATERNAL)")] = NA
  chr$Family = zoo::na.locf(chr$Family)
  colnames(chr)[2] = "ID"
  chr = chr[chr$V5 != "",]
  
  # Initialize list to store recombination events
  Recombinaison <- list()
  chr=chr[chr$V2!="(FOUNDER)",]
  
  # Loop through each individual to identify recombination events
  for (i in 1:nrow(chr)) {
    # Perform run-length encoding to identify haplotype switches (recombinations)
    rl <- rle(chr[i, ])
    # Check if there are recombinations 
    if (length(rl$values) - 4 != 0) {
      tmp <- data.frame("CHILD" = chr[i, 2], 
                        "PARENT"=chr[i, 3],
                        "Parent_ID"=if_else(chr[i, 3]=="(PATERNAL)", ped$V3[ped$V2==chr[i, 2]] , ped$V4[ped$V2==chr[i, 2]]),
                        "PosiRecombStart" = cumsum(rl$lengths[4:length(rl$values)]))
      # Store recombination information in list
      Recombinaison[[i]] <- tmp
    }
  }
  
  Recombinaison <- do.call(rbind, Recombinaison)
  Recombinaison = Recombinaison[Recombinaison$PosiRecombStart != ncol(chr)-3,] 
  Recombinaison = Recombinaison[order(Recombinaison$PosiRecombStart),]
  
  
  # Expand recombination intervals using heterozygous SNP information from parents
  
  # Identify individuals with no heterozygous markers 
  heterozygous$invariant_rows <- apply(heterozygous[,7:ncol(heterozygous)], 1, function(row) length(unique(row)) == 1)
  
  # Process recombinations where parent is heterozygous
  Recombinaison_hetero=Recombinaison[!Recombinaison$Parent_ID %in% heterozygous$Indiv[heterozygous$invariant_rows],]
  nc=ncol(Recombinaison_hetero)
  # Add SNP index column to map
  merlin_map$"SNP_ID"=1:nrow(merlin_map)
  
  # For each recombination in heterozygous parents, find flanking heterozygous markers
  for (i in 1:nrow(Recombinaison_hetero)) {
    hetero_status=as.logical(heterozygous[heterozygous$Indiv==Recombinaison_hetero$Parent_ID[i],][7:(ncol(heterozygous)-1)])
    tmp=merlin_map[hetero_status,]
    tmp$"end"=c(tmp$PHY.POS[-1], max(merlin_map$PHY.POS))
    tmp=tmp[tmp$SNP_ID==max(tmp$SNP_ID[tmp$SNP_ID<=Recombinaison_hetero$PosiRecombStart[i]]),]
    # Add start and end positions to recombination record
    Recombinaison_hetero[i,(nc+1):(nc+2)] <- c(tmp$PHY.POS, tmp$end)
  }
  # Rename new columns to indicate interval boundaries
  colnames(Recombinaison_hetero)[(nc+1):(nc+2)]=c("START", "END")
  # Mark that this recombination is in a heterozygous parent
  Recombinaison_hetero$"isParentHeteroz"=TRUE
  Recombinaison_hetero=Recombinaison_hetero[c(1,3,5:7)]

  
  # Process recombinations where parent is homozygous
  # Filter to keep only recombinations in homozygous parents
  Recombinaison_homo=Recombinaison[Recombinaison$Parent_ID %in% heterozygous$Indiv[heterozygous$invariant_rows],]
  # Prepare map data with end positions for homozygous case
  tmp=merlin_map
  tmp$"end"=c(tmp$PHY.POS[-1], NA)
  # Merge recombinations with map to get physical positions
  Recombinaison_homo=merge(Recombinaison_homo, tmp, by.x="PosiRecombStart", by.y="SNP_ID", all.x=T, all.y=F, sort=F)
  # Mark that this recombination is in a homozygous parent
  Recombinaison_homo$"isParentHeteroz"=FALSE
  # Select and reorder relevant columns
  Recombinaison_homo=Recombinaison_homo[c(2,4,8:10)]

  # Standardize column names across both heterozygous and homozygous recombination data frames
  colnames(Recombinaison_homo)=colnames(Recombinaison_hetero)= c("CHILD", "PARENT", "START", "END", "isParentHeteroz" )

  # Combine homozygous and heterozygous recombination results
  recombinaison=rbind(Recombinaison_homo,Recombinaison_hetero)
  # Write final recombination results to output file
  write.table(recombinaison,paste0("/lustre09/project/6033529/genealogy_sims/results/Samir/P1/paper3_Analyzes/Merlin_Recomb/",chrNum[j], ".merlin.3innuc.recombinaisons"),quote=F,row.names=F)
  
}
