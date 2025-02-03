
#This code is designed to infer recombination events from genetic data.
#It first loads necessary libraries and data,
#then processes the data to identify recombination locations from Merlin output.
#The recombination intervals are then broadened to include heterozygote SNPs of the parent.
#The results are finally written to a file.


# Clear the workspace
rm(list=ls())

# Load necessary libraries
library(dplyr)
library(stringr)
library(GenomicRanges)
library(BRGenomics)

# Define paths to data
data_path="/home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/SFTP hub/PHD -P1/Data/DataMerlin/hapAD/"
map_paths = c( paste0(data_path, "/chr", c(1:8),"/merlin_map.0",c(1:8)),
            paste0(data_path, "/chr", c(10:22), "/merlin_map.", c(10:22)))
flow_paths = paste0(data_path, "/chr", c(1:8,10:22),"/merlin.flow")
chromosome_numbers = paste0("chr", c(1:8,10:22))
heterozygous_paths=paste0(data_path, "het.",c(1:8,10:22),".RData")

# Load pedigree
pedigree=read.table("/home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/P1/Candiate_recombinations/WGS.ped")

# Loop over chromosomes
for (j in c(3,5,21)) { #others are already runned
  
  # Load heterozygous data
  heterozygous_markers=get(load(heterozygous_paths[j]))
  colnames(heterozygous_markers)[1:4]=c("Family", "Indiv","Father","Mother")
  
  # Load Merlin map
  merlin_map = read.delim2(map_paths[j], header = T, sep = "")
  
  # Check if the number of markers matches the map file
  cat("le nombre de marqueurs correspond avec le fichier map :chr", j)
  if (nrow(merlin_map)==ncol(heterozygous_markers)-6) {
    print("Oui")
  } else {
    print("Non")
    break
  }
  
  # Add physical position to Merlin map
  merlin_map = cbind(merlin_map, "PHY.POS"=as.numeric(str_split_fixed(merlin_map$MARKER, ":", 3)[,2]))
  
  # Load and process chromosome data
  chromosome_data = read.delim2(flow_paths[j], header = FALSE, sep = "")
  chromosome_data = cbind("Family" = chromosome_data$V2, chromosome_data)
  chromosome_data$Family[chromosome_data$Family %in% c("(FOUNDER)", "(MATERNAL)", "(PATERNAL)")] = NA
  chromosome_data$Family = zoo::na.locf(chromosome_data$Family)
  colnames(chromosome_data)[2] = "ID"
  chromosome_data = chromosome_data[chromosome_data$V5 != "",]
  
  # Identify recombination locations
  recombination_locations <- list()
  chromosome_data=chromosome_data[chromosome_data$V2!="(FOUNDER)",]
  for (i in 1:nrow(chromosome_data)) {
    rl <- rle(chromosome_data[i, ])
    if (length(rl$values) - 4 != 0) {
      tmp <- data.frame("CHILD" = chromosome_data[i, 2], 
                        "PARENT"=chromosome_data[i, 3],
                        "Parent_ID"=if_else(chromosome_data[i, 3]=="(PATERNAL)", pedigree$V3[pedigree$V2==chromosome_data[i, 2]] , pedigree$V4[pedigree$V2==chromosome_data[i, 2]]),
                        "PosiRecombStart" = cumsum(rl$lengths[4:length(rl$values)]))
      recombination_locations[[i]] <- tmp
    }
  }
  
  recombination_locations <- do.call(rbind, recombination_locations)
  recombination_locations = recombination_locations[recombination_locations$PosiRecombStart != ncol(chromosome_data)-3,] 
  recombination_locations = recombination_locations[order(recombination_locations$PosiRecombStart),]
  
  # Broaden recombination interval to heterozygote SNP of the parent
  heterozygous_markers$invariant_rows <- apply(heterozygous_markers[,7:ncol(heterozygous_markers)], 1, function(row) length(unique(row)) == 1)
  recombination_heterozygous=recombination_locations[!recombination_locations$Parent_ID %in% heterozygous_markers$Indiv[heterozygous_markers$invariant_rows],]
  nc=ncol(recombination_heterozygous)
  merlin_map$"SNP_ID"=1:nrow(merlin_map)
  for (i in 1:nrow(recombination_heterozygous)) {
    hetero_status=as.logical(heterozygous_markers[heterozygous_markers$Indiv==recombination_heterozygous$Parent_ID[i],][7:(ncol(heterozygous_markers)-1)])
    tmp=merlin_map[hetero_status,]
    tmp$"end"=c(tmp$PHY.POS[-1], NA)

    if (sum(tmp$SNP_ID<=recombination_heterozygous$PosiRecombStart[i])==0) {
      tmp=tmp[tmp$SNP_ID==min(tmp$SNP_ID[tmp$SNP_ID>=recombination_heterozygous$PosiRecombStart[i]]),]
      recombination_heterozygous[i,(nc+1):(nc+2)] <- c(merlin_map$PHY.POS[merlin_map$SNP_ID==recombination_heterozygous$PosiRecombStart[i]], tmp$PHY.POS) 
    }else{
      tmp=tmp[tmp$SNP_ID==max(tmp$SNP_ID[tmp$SNP_ID<=recombination_heterozygous$PosiRecombStart[i]]),]
      recombination_heterozygous[i,(nc+1):(nc+2)] <- c(tmp$PHY.POS, tmp$end)
    }
  }
  colnames(recombination_heterozygous)[(nc+1):(nc+2)]=c("START", "END")
  recombination_heterozygous$"isParentHeteroz"=TRUE
  recombination_heterozygous=recombination_heterozygous[c(1,3,5:7)]
  
  # Process homozygote parents
  recombination_homozygous=recombination_locations[recombination_locations$Parent_ID %in% heterozygous_markers$Indiv[heterozygous_markers$invariant_rows],]
  tmp=merlin_map
  tmp$"end"=c(tmp$PHY.POS[-1], NA)
  recombination_homozygous=merge(recombination_homozygous, tmp, by.x="PosiRecombStart", by.y="SNP_ID", all.x=T, all.y=F, sort=F)
  recombination_homozygous$"isParentHeteroz"=FALSE
  recombination_homozygous=recombination_homozygous[c(2,4,8:10)]

  colnames(recombination_homozygous)=colnames(recombination_heterozygous)= c("CHILD", "PARENT", "START", "END", "isParentHeteroz" )

  # Combine results and write to file
  recombination=rbind(recombination_homozygous,recombination_heterozygous)
  #write.table(recombination,paste0("/home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/P1/paper3_Analyzes/Merlin_Recomb/",chromosome_numbers[j], ".merlin.recombinaisons"),quote=F,row.names=F)
  
}
