rm(list=ls())
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggforce)

#Path of data in rorqual
rpath="~/links/projects/rrg-girardsi/genealogy_sims/results/Samir/"

#Fichier de familles
ped=read.table(paste0(rpath,"/P1/Candiate_recombinations/WGS.ped"))

#La sequence complete
chemin2Recomb=paste0(rpath,"/SFTP hub/PHD -P1/Data/")
chr=c(1:9,10:22)
cheminSC = paste0(chemin2Recomb, "chr",chr,".pass.clean.recombinations")
cheminSC.3innuc = paste0(chemin2Recomb, "chr",chr,".pass.clean.3innuc.recombinations")
Recom.SC= read.delim2(cheminSC[1], header = TRUE, sep = "")
Recom.SC=rbind(Recom.SC,read.delim2(cheminSC.3innuc[1], header = TRUE, sep = ""))
Recom.SC=data.frame("chr"=chr[1], Recom.SC )

for (j in 2:length(chr)) {
  tmp=read.delim2(cheminSC[j], header = TRUE, sep = "")
  tmp=rbind(tmp,read.delim2(cheminSC.3innuc[j], header = TRUE, sep = ""))
  #tmp=tmp[tmp$PROB_RECOMBINATION >=  0.5,]
  Recom.SC = rbind(Recom.SC,data.frame("chr"=chr[j], tmp))
  #print(nrow(Recom.SC))
}

Recom.SC$"longueur"=Recom.SC$END-Recom.SC$START+1
CS_pair=(unique(paste(Recom.SC$PARENT,Recom.SC$CHILD,sep="-")))
Recom.SC=Recom.SC[Recom.SC$PROB_RECOMBINATION >=  0.5,]

#Recombinaison de la puce omniExpress
path2RecombOE=paste0(rpath,"/SFTP hub/PHD -P1/Data/")
cheminOE = paste0(path2RecombOE, "chr",chr,".OmniExpress.clean.recombinations")
cheminOE.3innuc = paste0(path2RecombOE,"chr",chr,".OmniExpress.clean.3innuc.recombinations")

Recom.OE= read.delim2(cheminOE[1], header = TRUE, sep = "")
Recom.OE=rbind(Recom.OE,read.delim2(cheminOE.3innuc[1], header = TRUE, sep = "") )
Recom.OE=data.frame("chr"=chr[1], Recom.OE )

for (j in 2:length(chr)) {
  tmp=read.delim2(cheminOE[j], header = TRUE, sep = "")
  tmp=rbind(tmp, read.delim2(cheminOE.3innuc[j], header = TRUE, sep = "") )
  #tmp=tmp[tmp$PROB_RECOMBINATION >=  0.5,]
  
  Recom.OE = rbind(Recom.OE,data.frame("chr"=chr[j], tmp ))
  #print(nrow(Recom.OE))
}
Recom.OE$"longueur"=Recom.OE$END-Recom.OE$START+1
OE_pair=(unique(paste(Recom.OE$PARENT,Recom.OE$CHILD,sep="-")))
Recom.OE=Recom.OE[Recom.OE$PROB_RECOMBINATION >= 0.5,]

#Recombinaison Merlin
chrNum = paste0("chr", chr)

#les recombinaisons Merlin dans les fichiers recombinaison* == chr*.merlin.recombinaisons (concerne une partie des familles)
CheminRM=paste0(rpath,"/P1/paper3_Analyzes/Merlin_Recomb/", 
                "chr", chr, ".merlin.recombinaisons")
CheminRM.3innuc=paste0(rpath,"/P1/paper3_Analyzes/Merlin_Recomb/", 
                       "chr", chr, ".merlin.3innuc.recombinaisons")
Recom.Merlin = read.delim2(CheminRM[1], header = TRUE, sep = "") 
Recom.Merlin=rbind(Recom.Merlin, read.table(CheminRM.3innuc[1], header = T) )
Recom.Merlin$"chr"=chr[1]

for (j in 2:length(chr)) {
  tmp=read.delim2(CheminRM[j], header = TRUE, sep = "")
  tmp=rbind(tmp, read.table(CheminRM.3innuc[j], header = T) )
  tmp$"chr"=chr[j]
  Recom.Merlin=rbind(Recom.Merlin, tmp )
}
#length(unique(paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-")))

famsSC=(unique(ped$V1[ped$V2 %in% Recom.SC$CHILD]))
Recom.Merlin=  Recom.Merlin[Recom.Merlin$CHILD %in% ped$V2[ped$V1 %in% famsSC], ]
row.names(Recom.Merlin) = NULL
Recom.Merlin=Recom.Merlin[!is.na(Recom.Merlin$START),]
Recom.Merlin=Recom.Merlin[!is.na(Recom.Merlin$END),]
Recom.Merlin=Recom.Merlin[!duplicated(Recom.Merlin),]



#########################################
#Studing overlap....
#########################################
common_indiv <- unique(Reduce(intersect, list(paste(Recom.OE$PARENT,Recom.OE$CHILD,sep="-"),
                                              paste(Recom.SC$PARENT,Recom.SC$CHILD,sep="-"),
                                              paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-"))))

#Filter by probabily for SC and OE
Recom.Merlin=Recom.Merlin[paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-") %in% common_indiv , ]
Recom.OE=Recom.OE[paste(Recom.OE$PARENT,Recom.OE$CHILD,sep="-") %in% common_indiv , ]
Recom.SC=Recom.SC[paste(Recom.SC$PARENT,Recom.SC$CHILD,sep="-") %in% common_indiv , ]

NR_expected_common= data.frame(chr=chr,
                               Expected=round(c(284.8,274.6,233.0,216.4,213.0, 197.5, 190.9,
                                                171.5, 168.5, 175.5, 165.9, 168.1,124.1, 124.9, 127.0,
                                                133.5, 138.8, 121.3, 102.6, 102.6, 63.7, 69.5)*length(common_indiv)/100 , 0))


#Recombination length distribution when overlaped

# Load the required library
library(GenomicRanges)

# Function to calculate overlap percentage
prop_overlap <- function(df1, df2) {
  gr1 <- GRanges(seqnames = paste0(df1$chr,"-",df1$CHILD), ranges = IRanges(start = df1$START, end = df1$END))
  gr2 <- GRanges(seqnames = paste0(df2$chr,"-",df2$CHILD), ranges = IRanges(start = df2$START, end = df2$END))
  
  ov <- findOverlaps(gr1, gr2)
  #df1$overlap <- 0
  overlap =rep(0, nrow(df1))
  #df1$overlap[queryHits(ov)] <- width(pintersect(gr1[queryHits(ov)], gr2[subjectHits(ov)])) / width(gr1[queryHits(ov)])
  overlap[queryHits(ov)] <- width(pintersect(gr1[queryHits(ov)], gr2[subjectHits(ov)])) / width(gr1[queryHits(ov)])
  
  # Convert to percentage and return
  return(overlap*100)
}

Recom.Merlin$"ov_OE" <- prop_overlap(Recom.Merlin, Recom.OE)
Recom.Merlin$"ov_SC" <- prop_overlap(Recom.Merlin, Recom.SC)

# Apply the function to each pair of dataframes
#Recom.OE <- calculate_overlap(Recom.OE, Recom.Merlin)
Recom.OE$"ov_Mr" <- prop_overlap(Recom.OE, Recom.Merlin)
Recom.OE$"ov_SC" <- prop_overlap(Recom.OE, Recom.SC)

Recom.SC$"ov_OE" <- prop_overlap(Recom.SC, Recom.OE)
Recom.SC$"ov_Mr" <- prop_overlap(Recom.SC, Recom.Merlin)


# Shapeit2-duoHMM on WGS concordance in %
#table(Recom.Merlin$ov_SC>0)[2]*100/nrow(Recom.SC)
table(Recom.SC$ov_Mr>0)[2]*100/nrow(Recom.SC)

# Shapeit2-duoHMM on OE concordance in %
#table(Recom.Merlin$ov_SC>0)[2]*100/nrow(Recom.OE)
table(Recom.OE$ov_Mr>0)[2]*100/nrow(Recom.OE)
write.csv2(Recom.OE, file=paste0(rpath,"/P1/paper3_Analyzes/Recom.OE.csv"), row.names = FALSE)

# Shapeit2-duoHMM on WGS & OE  concordance in %
nrow(Recom.SC[Recom.SC$ov_OE>0 & Recom.SC$ov_Mr>0,])*100/table(Recom.SC$ov_OE>0)[2]
write.csv2(Recom.SC, file=paste0(rpath,"/P1/paper3_Analyzes/Recom.SC.csv"), row.names = FALSE)


# Concordance per chromosome
# Recomb. shapeit2WGS-MerlinOE
subset_data <- Recom.SC[Recom.SC$ov_Mr > 0, ]
result <- as.data.frame(table(subset_data$chr))
colnames(result)=c("chr", "shapeit2WGS-MerlinOE")
Concordance_df=result

# Recomb. shapeit2OE-MerlinOE
subset_data <- Recom.OE[Recom.OE$ov_Mr > 0, ]
result <- as.data.frame(table(subset_data$chr))
colnames(result)=c("chr", "shapeit2OE-MerlinOE")
Concordance_df=merge(Concordance_df, result, by="chr")

# Recomb. shapeit2SC-shapeit2OE
subset_data <- Recom.SC[Recom.SC$ov_OE > 0, ]
result <- as.data.frame(table(subset_data$chr))
colnames(result)=c("chr", "shapeit2SC-shapeit2OE")
Concordance_df=merge(Concordance_df, result, by="chr")

# Recomb. "shapeit2WGS-shapeit2OE-MerlinOE"
subset_data <- Recom.SC[Recom.SC$ov_OE>0 & Recom.SC$ov_Mr>0,]
result <- as.data.frame(table(subset_data$chr))
colnames(result)=c("chr", "shapeit2WGS-shapeit2OE-MerlinOE")
Concordance_df=merge(Concordance_df, result, by="chr")


