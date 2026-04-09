#/lustre03/project/6033529/genealogy_sims/results/Samir/SFTP hub/PHD -P1/Data/to delete/data21122020/chr17.AA.OmniExpress.clean.recombinations


rm(list=ls())
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggforce)

#This version erasing the old that not include chrom 9.
#This version is adapted for rorqual paths

#Path of data in rorqual
rpath="~/links/projects/rrg-girardsi/genealogy_sims/results/Samir/"

#Fichier de familles
ped=read.table(paste0(rpath,"/P1/Candiate_recombinations/WGS.ped"))

#AA African-American  recombination  shapeit2-duoHMM on OE chip data
chemin2Recomb=paste0(rpath,"/SFTP hub/PHD -P1/Data/to delete/data21122020/")
chr=c(1:9,10:22)
cheminAA = paste0(chemin2Recomb, "chr",chr,".AA.OmniExpress.clean.recombinations")
Recom.AA= read.delim2(cheminAA[1], header = TRUE, sep = "")
Recom.AA=data.frame("chr"=chr[1], Recom.AA )

for (j in 2:length(chr)) {
  tmp=read.delim2(cheminAA[j], header = TRUE, sep = "")
  Recom.AA = rbind(Recom.AA,data.frame("chr"=chr[j], tmp))
}

Recom.AA$"longueur"=Recom.AA$END-Recom.AA$START+1
AA_pair=(unique(paste(Recom.AA$PARENT,Recom.AA$CHILD,sep="-")))
Recom.AA=Recom.AA[Recom.AA$PROB_RECOMBINATION >=  0.5,]


#AA African-American  recombination  shapeit2-duoHMM on OE chip data
#chemin2Recomb="/lustre03/project/6033529/genealogy_sims/results/Samir/SFTP hub/PHD -P1/Data/to delete/data21122020/"
chr=c(1:9,10:22)
cheminEA = paste0(chemin2Recomb, "chr",chr,".EA.OmniExpress.clean.recombinations")
Recom.EA= read.delim2(cheminEA[1], header = TRUE, sep = "")
Recom.EA=data.frame("chr"=chr[1], Recom.EA )

for (j in 2:length(chr)) {
  tmp=read.delim2(cheminEA[j], header = TRUE, sep = "")
  Recom.EA = rbind(Recom.EA,data.frame("chr"=chr[j], tmp))
}

Recom.EA$"longueur"=Recom.EA$END-Recom.EA$START+1
EA_pair=(unique(paste(Recom.EA$PARENT,Recom.EA$CHILD,sep="-")))
Recom.EA=Recom.EA[Recom.EA$PROB_RECOMBINATION >=  0.5,]


#Recombinaison Merlin
chrNum = paste0("chr", chr)

#les recombinaisons Merlin dans les fichiers recombinaison* == chr*.merlin.recombinaisons (concerne une partie des familles)
#CheminRM = paste0("/home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/SFTP hub/PHD -P1/Data/Recombinaison_Merlin/recombinaison_",chrNum)
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
row.names(Recom.Merlin) = NULL
Recom.Merlin=Recom.Merlin[!is.na(Recom.Merlin$START),]
Recom.Merlin=Recom.Merlin[!is.na(Recom.Merlin$END),]
Recom.Merlin=Recom.Merlin[!duplicated(Recom.Merlin),]


#########################################
#Restreintre au pair en commun entre merlin  et union(EA,AA)
#########################################
common_indiv <- unique(intersect(c(paste(Recom.AA$PARENT,Recom.AA$CHILD,sep="-"),
                                   paste(Recom.EA$PARENT,Recom.EA$CHILD,sep="-"))
                                 , paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-")))

#Filter by probabily for SC and OE
Recom.Merlin=Recom.Merlin[paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-") %in% common_indiv , ]
Recom.AA=Recom.AA[paste(Recom.AA$PARENT,Recom.AA$CHILD,sep="-") %in% common_indiv , ]
Recom.EA=Recom.EA[paste(Recom.EA$PARENT,Recom.EA$CHILD,sep="-") %in% common_indiv , ]

NR_expected_common= data.frame(chr=chr,
                               Expected=round(c(284.8,274.6,233.0,216.4,213.0, 197.5, 190.9,
                                                171.5, 168.5, 175.5, 165.9, 168.1,124.1, 124.9, 127.0,
                                                133.5, 138.8, 121.3, 102.6, 102.6, 63.7, 69.5)/100 , 2))


#Plot the number of recombination event by chromosome
# Before filtre
df=merge(as.data.frame(table(Recom.Merlin$chr)), as.data.frame(table(Recom.AA$chr)), by="Var1")
df=merge(df, as.data.frame(table(Recom.EA$chr)), by="Var1")
#colnames(df)[1:4]=c("chr", "Merlin_OE", "Shapeit2_AA", "Shapeit2_EA")
df=merge(df, as.data.frame(table(Recom.Merlin$chr[paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-") %in% AA_pair])), by="Var1")
df=merge(df, as.data.frame(table(Recom.Merlin$chr[paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-") %in% EA_pair])), by="Var1")
colnames(df)=c("chr", "Merlin_OE", "Shapeit2_AA", "Shapeit2_EA", "Merlin_AA", "Merlin_EA")
df=merge(df, NR_expected_common, by="chr")

#per meiosis
df$Merlin_OE=round(df$Merlin_OE/length(unique(paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-"))),2)
df$Shapeit2_AA=round(df$Shapeit2_AA/length(unique(paste(Recom.AA$PARENT,Recom.AA$CHILD,sep="-"))),2)
df$Shapeit2_EA=round(df$Shapeit2_EA/length(unique(paste(Recom.EA$PARENT,Recom.EA$CHILD,sep="-"))),2)
df$Merlin_AA=round(df$Merlin_AA/sum(AA_pair  %in% paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-")),2)
df$Merlin_EA=round(df$Merlin_EA/sum(EA_pair  %in% paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-")),2)
df

#comparer les % de chevauchement dans un groupe versus l'autre groupe

# Assuming df is your dataframe
# Melt the dataframe to long format for ggplot
df_melt <- df %>% pivot_longer(cols = -chr)
df_melt=df_melt[df_melt$name !="Merlin_OE",]
p=ggplot(df_melt, aes(chr, value, group=name, colour=name)) +
  geom_line() +
  labs(x="Chromosome number", y="Number of recombination per meiosis", colour=" ")
#facet_zoom(y=(value <= 900), zoom.size = 1.2, show.area =F)
p

#########################################
#Comparaison des deux groupes statistiquement: nombre de recombination
#########################################

# Appliquer la fonction et afficher les résultats
peds=read.table(paste0(rpath,"/P1/Candiate_recombinations/WGS.ped"))


# Sample data
# get recombination events per family
Recom.AA2= merge(Recom.AA,peds[,1:2], by.x="CHILD", by.y="V2", all.x=TRUE )
colnames(Recom.AA2)[ncol(Recom.AA2)]="FID"

LenByfamily <- Recom.AA2 %>%
  group_by(FID) %>%
  summarise(
    NRecom = n(),
    num_pairs = n_distinct(paste(PARENT, CHILD, sep = "-"))
  )

group1<- LenByfamily$NRecom/LenByfamily$num_pairs

Recom.EA2= merge(Recom.EA,peds[,1:2], by.x="CHILD", by.y="V2", all.x=TRUE )
colnames(Recom.EA2)[ncol(Recom.EA2)]="FID"

LenByfamily <- Recom.EA2 %>%
  group_by(FID) %>%
  summarise(
    NRecom = n(),
    num_pairs = n_distinct(paste(PARENT, CHILD, sep = "-"))
  )

group2<- LenByfamily$NRecom/LenByfamily$num_pairs


# Appliquer le test de Mann-Whitney-Wilcoxon
test_result <- wilcox.test(group1, group2, alternative = "two.sided")

# Afficher les résultats
cat("Statistique de test:", test_result$statistic, "\n")
cat("P-value:", test_result$p.value, "\n")



#Recombination length distribution in general
df_length_Recom=data.frame(type="SHapeit2_AA",Len=Recom.AA$END-Recom.AA$START+1, chr=Recom.AA$chr  )
df_length_Recom=rbind(df_length_Recom, data.frame(type="SHapeit2_EA",Len=Recom.EA$END-Recom.EA$START+1 , chr=Recom.EA$chr))
df_length_Recom=rbind(df_length_Recom, data.frame(type="Merlin_OE",Len=Recom.Merlin$END-Recom.Merlin$START+1 , chr=Recom.Merlin$chr ))
tmp=Recom.Merlin[paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-") %in% AA_pair,]
df_length_Recom=rbind(df_length_Recom, data.frame(type="Merlin_AA",Len=tmp$END-tmp$START+1 , chr=tmp$chr ))
tmp=Recom.Merlin[paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-") %in% EA_pair,]
df_length_Recom=rbind(df_length_Recom, data.frame(type="Merlin_EA",Len=tmp$END-tmp$START+1 , chr=tmp$chr ))

Recom.OE=read.csv2(paste0(rpath,"/P1/paper3_Analyzes/Recom.OE.csv"))
Recom.OE=Recom.OE[paste(Recom.OE$PARENT,Recom.OE$CHILD,sep="-") %in% common_indiv, ]

df_length_Recom=rbind(df_length_Recom, data.frame(type="SHapeit2_OE",Len=Recom.OE$longueur , chr=Recom.OE$chr ))

# Create the plot with log scale on the y-axis
ggplot(df_length_Recom, aes(type, log10(Len), colour = type)) +
  geom_boxplot() +
  #scale_y_log10()+
  labs(x = "Method_data", y = "log-length of recombination intervals") +
  theme(legend.position = "none")


#########################################
#Comparaison des deux groupes statistiquement : la longueur des intervalles
#########################################

# Load libraries
library(car) # For Levene's test

# Sample data
# get recombination events per family
peds=read.table(paste0(rpath,"/P1/Candiate_recombinations/WGS.ped"))
Recom.AA2= merge(Recom.AA,peds[,1:2], by.x="CHILD", by.y="V2", all.x=TRUE )
colnames(Recom.AA2)[ncol(Recom.AA2)]="FID"

LenByfamily <- Recom.AA2 %>%
  group_by(FID) %>%
  summarise(
    Mean_len = mean((longueur))
  )
group1 <- LenByfamily$Mean_len

Recom.EA2= merge(Recom.EA,peds[,1:2], by.x="CHILD", by.y="V2", all.x=TRUE )
colnames(Recom.EA2)[ncol(Recom.EA2)]="FID"

LenByfamily <- Recom.EA2 %>%
  group_by(FID) %>%
  summarise(
    Mean_len = mean((longueur))
  )
group2 <- LenByfamily$Mean_len


# Shapiro-Wilk test for normality : Si le p-value > 0.05, les données sont normalement distribuées.
shapiro_test_group1 <- shapiro.test(group1)
shapiro_test_group2 <- shapiro.test(group2)

# Levene's test for equality of variances :  Si le p-value > 0.05, les variances peuvent être considérées comme égales.
levene_test <- leveneTest(c(group1, group2) ~ factor(c(rep(1, length(group1)), rep(2, length(group2)))))

# Perform appropriate test
if (shapiro_test_group1$p.value > 0.05 && shapiro_test_group2$p.value > 0.05) {
  # Data is normal
  if (levene_test$`Pr(>F)`[1] > 0.05) {
    # Variances are equal
    test_result <- t.test(group1, group2, var.equal = TRUE)
  } else {
    # Variances are not equal
    test_result <- t.test(group1, group2, var.equal = FALSE)
  }
} else {
  # Data is not normal
  test_result <- wilcox.test(group1, group2)
}

# Output results
list(
  shapiro_test_group1 = shapiro_test_group1,
  shapiro_test_group2 = shapiro_test_group2,
  levene_test = levene_test,
  test_result = test_result
)



################################
#Mean of log-transformed value

LenByfamily <- Recom.AA2 %>%
  group_by(FID) %>%
  summarise(
    Mean_len = mean(log(longueur))
  )
group1 <- LenByfamily$Mean_len

Recom.EA2= merge(Recom.EA,peds[,1:2], by.x="CHILD", by.y="V2", all.x=TRUE )
colnames(Recom.EA2)[ncol(Recom.EA2)]="FID"

LenByfamily <- Recom.EA2 %>%
  group_by(FID) %>%
  summarise(
    Mean_len = mean(log(longueur))
  )
group2 <- LenByfamily$Mean_len


# Calculer les moyennes des deux groupes
mean_group1 <- mean((group1))
mean_group2 <- mean((group2))

# Calculer la différence absolue
absolute_difference <- abs(mean_group1 - mean_group2)
absolute_difference
# Calculer le seuil (5% de la moyenne des deux groupes)
#threshold <- 0.05 * mean(c(mean_group1, mean_group2))

# Comparer la différence au seuil
# if (absolute_difference < threshold) {
#   cat("La différence absolue est inférieure au seuil.\n")
#   cat("Différence absolue:", absolute_difference, "\n")
#   cat("Seuil:", threshold, "\n")
#   } else {
#   cat("La différence absolue dépasse le seuil.\n")
#   cat("Différence absolue:", absolute_difference, "\n")
#   cat("Seuil:", threshold, "\n")
# 
# }






#########################################
#Studing overlap....
#########################################


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
  # Convert to percentage
  #df1$overlap <- df1$overlap * 100
  
  #return(df1)
  return(overlap*100)
}

# Apply the function to each pair of dataframes
Recom.Merlin$"ov_AA" <- prop_overlap(Recom.Merlin, Recom.AA)
Recom.Merlin$"ov_EA" <- prop_overlap(Recom.Merlin, Recom.AA)

Recom.AA$"ov_Mr" <- prop_overlap(Recom.AA, Recom.Merlin)

Recom.EA$"ov_Mr" <- prop_overlap(Recom.EA, Recom.Merlin)


# Shapeit2-duoHMM on OE-AA concordance in %
table(Recom.AA$ov_Mr>0)[2]*100/nrow(Recom.AA)

# Shapeit2-duoHMM on OE-EA concordance in %
table(Recom.EA$ov_Mr>0)[2]*100/nrow(Recom.EA)


# Concordance per chromosome
# Recomb. shapeit2AA-MerlinOE
subset_data <- Recom.AA[Recom.AA$ov_Mr > 0, ]
result <- as.data.frame(table(subset_data$chr))
colnames(result)=c("chr", "shapeit2AA-MerlinOE")
Concordance_df=result

# Recomb. shapeit2EA-MerlinOE
subset_data <- Recom.EA[Recom.EA$ov_Mr > 0, ]
result <- as.data.frame(table(subset_data$chr))
colnames(result)=c("chr", "shapeit2EA-MerlinOE")
Concordance_df=merge(Concordance_df, result, by="chr")

Concordance_df$`shapeit2AA-Merlin`=Concordance_df$`shapeit2AA-MerlinOE`/sum(AA_pair  %in% paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-"))
Concordance_df$`shapeit2EA-Merlin`=Concordance_df$`shapeit2EA-MerlinOE`/sum(EA_pair  %in% paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-"))


Concordance_df$`shapeit2AA-Merlin`=round(Concordance_df$`shapeit2AA-Merlin`, 2)
Concordance_df$`shapeit2EA-Merlin`=round(Concordance_df$`shapeit2EA-Merlin`, 2)







