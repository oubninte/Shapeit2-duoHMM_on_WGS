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
#path to the files
chemin2Recomb=paste0(rpath,"/SFTP hub/PHD -P1/Data/")
chr=c(1:22)
cheminSC = paste0(chemin2Recomb, "chr",chr,".pass.clean.recombinations")
cheminSC.3innuc = paste0(chemin2Recomb, "chr",chr,".pass.clean.3innuc.recombinations")

#Concatenate of chr9 : concatenate 9q and 9p data and save files (9p (29) et 9q (9))
#Recom.SC.chr9= read.delim2(paste0(chemin2Recomb, "chr9p",".pass.clean.recombinations"),
#                           header = TRUE, sep = "")
#Recom.SC.chr9=rbind(Recom.SC.chr9,
#                    read.delim2(paste0(chemin2Recomb, "chr9q",".pass.clean.recombinations"), header = TRUE, sep = ""))
#write.table(Recom.SC.chr9, paste0(chemin2Recomb, "chr9",".pass.clean.recombinations"), quote = F , row.names = F)

#chr9 .3innuc
#Recom.SC.chr9= read.delim2(paste0(chemin2Recomb, "chr9p",".pass.clean.3innuc.recombinations"),
#                           header = TRUE, sep = "")
#Recom.SC.chr9=rbind(Recom.SC.chr9,
#                    read.delim2(paste0(chemin2Recomb, "chr9q",".pass.clean.3innuc.recombinations"), header = TRUE, sep = ""))
#write.table(Recom.SC.chr9, paste0(chemin2Recomb, "chr9",".pass.clean.3innuc.recombinations"), quote = F , row.names = F)


#load the two files of chr1 data : we should have 61 famlies 
Recom.SC= read.delim2(cheminSC[1], header = TRUE, sep = "")
Recom.SC=rbind(Recom.SC,read.delim2(cheminSC.3innuc[1], header = TRUE, sep = ""))
Recom.SC=data.frame("chr"=chr[1], Recom.SC )

#load the two files of each other chroms
for (j in 2:length(chr)) {
  tmp=read.delim2(cheminSC[j], header = TRUE, sep = "")
  tmp=rbind(tmp,read.delim2(cheminSC.3innuc[j], header = TRUE, sep = ""))
  #tmp=tmp[tmp$PROB_RECOMBINATION >=  0.5,]
  Recom.SC = rbind(Recom.SC,data.frame("chr"=chr[j], tmp))
  #print(nrow(Recom.SC))
}

#Computing length recomb interval
Recom.SC$"longueur"=Recom.SC$END-Recom.SC$START+1
CS_pair=(unique(paste(Recom.SC$PARENT,Recom.SC$CHILD,sep="-")))

#keeping recomb with probab > 0.5
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

#les recombinaisons Merlin dans les fichiers recombinaison* == chr*.merlin.recombinaisons (concerne une partie des familles) produit par ce script
#"/lustre09/project/6033529/genealogy_sims/results/Samir/P1/paper3_Analyzes/Recomb_Merlin_UncutFam_allchr.R" 
#and partie sont produit par ce script "/lustre09/project/6033529/genealogy_sims/results/Samir/P1/paper3_Analyzes/Recomb_Merlin_allchr.R"
#allchr inclus chr9 et adapter pour rorqual

CheminRM=paste0(rpath, "/P1/paper3_Analyzes/Merlin_Recomb/", "chr", chr, ".merlin.recombinaisons")
CheminRM.3innuc=paste0(rpath, "/P1/paper3_Analyzes/Merlin_Recomb/", "chr", chr, ".merlin.3innuc.recombinaisons")

Recom.Merlin = read.delim2(CheminRM[1], header = TRUE, sep = "") 
Recom.Merlin=rbind(Recom.Merlin, read.table(CheminRM.3innuc[1], header = T) )
Recom.Merlin$"chr"=chr[1]

for (j in 2:length(chr)) {
  tmp=read.delim2(CheminRM[j], header = TRUE, sep = "")
  tmp=rbind(tmp, read.table(CheminRM.3innuc[j], header = T) )
  tmp$"chr"=chr[j]
  Recom.Merlin=rbind(Recom.Merlin, tmp )
}

length(unique(paste(Recom.Merlin$PARENT,Recom.Merlin$CHILD,sep="-")))
famsSC=(unique(ped$V1[ped$V2 %in% Recom.SC$CHILD]))
Recom.Merlin=  Recom.Merlin[Recom.Merlin$CHILD %in% ped$V2[ped$V1 %in% famsSC], ]
row.names(Recom.Merlin) = NULL
Recom.Merlin=Recom.Merlin[!is.na(Recom.Merlin$START),]
Recom.Merlin=Recom.Merlin[!is.na(Recom.Merlin$END),]
Recom.Merlin=Recom.Merlin[!duplicated(Recom.Merlin),]

#########################################
#Restreintre au pair en commun
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

#Plot the number of recombination event by chromosome
# Before filtre
df=merge(as.data.frame(table(Recom.Merlin$chr)), as.data.frame(table(Recom.OE$chr)), by="Var1")
df=merge(df, as.data.frame(table(Recom.SC$chr)), by="Var1")
colnames(df)[1:4]=c("chr", "Merlin_OE", "Shapeit2_OE", "Shapeit2_WGS")
df=merge(df, NR_expected_common, by="chr")
#df=merge(df, NR_expected_common, by="chr")
#colnames(df)[5:6]=c("E.Commom", "E.Global")
#df=merge(df, NR_expected, by="chr")

# Assuming df is your dataframe
# Melt the dataframe to long format for ggplot
df_melt <- df %>% pivot_longer(cols = -chr)

ggplot(df_melt, aes(chr, value, group=name, colour=name)) +
  geom_line() +
  labs(x="Chromosome number", y="Number of recombination events", colour=" ") 
  #facet_zoom(y=(value <= 900), zoom.size = 1.2, show.area =F)

# #Add Ingo recommandations in figure2 of paper:
# df_melt <- df %>% pivot_longer(cols = -chr)
# df_melt=rbind(df_melt, data.frame(
#   chr = rep("9", 4),
#   name = c("Merlin_OE", "Shapeit2_OE", "Shapeit2_WGS", "Expected"),
#   value = c(rep(NA, 3), 168.5*length(common_indiv)/100 ) 
# ))
# df_melt$chr <- factor(df_melt$chr, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"))


ggplot(df_melt, aes(chr, value/length(common_indiv), group = name, colour = name)) +
  geom_line() +
  #geom_line(data = subset(df_melt, !is.na(value)), linetype = "dotted", size = 0.5) +
  geom_point(size = 0.5) +
  labs(x = "Chromosome number",
       y = "Number of recombination events per meiosis",
       colour = " ") +
  #scale_x_discrete(drop = FALSE) +  # Ensure missing levels (like chr9) are shown on the axis
  scale_y_continuous(    sec.axis = sec_axis(~ . * 295, name = "Number of recombination events")  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels for better readability
  theme_minimal()


#Recombination length distribution in general
df_length_Recom=data.frame(type="SHapeit2_OE",Len=Recom.OE$END-Recom.OE$START+1, chr=Recom.OE$chr  )
df_length_Recom=rbind(df_length_Recom, data.frame(type="SHapeit2_WGS",Len=Recom.SC$END-Recom.SC$START+1 , chr=Recom.SC$chr))
df_length_Recom=rbind(df_length_Recom, data.frame(type="Merlin_OE",Len=Recom.Merlin$END-Recom.Merlin$START+1 , chr=Recom.Merlin$chr ))


# Create the plot with log scale on the y-axis
ggplot(df_length_Recom, aes(type, log10(Len), colour = type)) +
  geom_boxplot() +
  #scale_y_log10()+
  labs(x = "Method_data", y = "log-length of recombination intervals") +
  theme(legend.position = "none")

# Create the plot with log scale on the y-axis and facet by chromosome
ggplot(df_length_Recom, aes(type, log10(Len), colour = type)) +
  geom_boxplot() +
  #scale_y_log10() +
  labs(x = "", y = "log-length of recombination intervals", colour = "Method-Data") +
  facet_wrap(~ chr, ncol = 6)+
  theme(axis.text.x = element_blank())


# Calculate descriptive statistics for each 'type'
summary_stats <- df_length_Recom %>%
  group_by(type) %>%
  summarise(N = n(),
            min = min(Len),
            Q1 = quantile(Len, 0.25),
            median = median(Len),
            mean = mean(Len),
            Q3 = quantile(Len, 0.75),
            max = max(Len))

print(summary_stats)

#Recap :
#colSums(df[2:5]) #Number of recombination by method
#round(colSums(df[2:5])/length(common_indiv),2) #Number of recombination per meiosis by method
recap1=data.frame("Number of recombination"=colSums(df[2:5]), 
                  "Number of recombination per meiosis"= round(colSums(df[2:5])/length(common_indiv),2) )


# fusion .... plot recombination segments by method
Recom.Merlin=Recom.Merlin[order(Recom.Merlin$START),]
Recom.SC=Recom.SC[order(Recom.SC$START),]
Recom.OE=Recom.OE[order(Recom.OE$START),]

Recom.All=rbind(Recom.SC[c(1:5)], Recom.OE[c(1:5)], Recom.Merlin[c(6,1:4)])
Recom.All$Method= c(rep("SC", nrow(Recom.SC)), rep("OE", nrow(Recom.OE)), rep("MR", nrow(Recom.Merlin)))
Recom.All=Recom.All[order(Recom.All$START),]
Recom.All$Y=NA
for (chrom in unique(Recom.All$chr)) {
  for (method in unique(Recom.All$Method)) {
    nr=nrow(Recom.All[Recom.All$chr==chrom & Recom.All$Method==method,])
    Recom.All$Y[Recom.All$chr==chrom & Recom.All$Method==method]=sample(1:nr, nr)
  }
}

#plot
p=ggplot(Recom.All[Recom.All$chr %in% c(1) ,], aes(x = 1, y = 1, color = factor(Method))) +
  geom_segment(aes(x = START, xend = END, 
                   y = Y, 
                   yend = Y), size = 1)+
  facet_wrap(~ Method, ncol = 1)+
  theme(legend.position = "none")+           # Hide the legend
  ylab("Y Axis Label")+    # Add label for Y axis
  xlab(NULL)+
  labs(title="Recombination distribution by Method")

p

# Number of recombination per Mbp 
# Create the output dataframe
result <- Recom.All %>%
  mutate(position = START / 1e6) %>%
  group_by(chr, Method, position = floor(position)+1) %>%
  summarise(rate = n(), .groups = 'drop') %>%
  rename(position_Mbp = position)

#plot
ggplot(result[result$chr %in% c(1),], aes(x = position_Mbp  , y = rate, color = Method)) +
  geom_line() +
  facet_wrap(~ chr, scales = "free_x", ncol = 3) +
  labs(title = "",
       x = "Position (Mb)",
       y = "Number of recombination by Mb") +
  theme_minimal() +
  theme(strip.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept = length(common_indiv)/100, linetype = "dashed", color = "red") +
  geom_text(aes(x = -05, y = 05, label = paste0("Expected=",length(common_indiv)/100)), vjust = -0.050, hjust = -0.050, color = "black", size = 3)



#Lissage: methode de LOESS (span contrôle le lissage)
#chromosome 9
ggplot(result[result$chr == 9,], aes(x = position_Mbp, y = rate, color = Method)) +
  geom_line(alpha = 0.3) +
  stat_smooth(method = "loess", span = 0.1, se = FALSE, size = 1) +
  geom_hline(yintercept = length(common_indiv)/100, linetype = "dashed", color = "red") +
  annotate("text",
           x = -5, y = 5,
           label = paste0("Expected=", length(common_indiv)/100),
           vjust = -Inf, hjust = -Inf,
           color = "black", size = 3) +
  facet_wrap(~ chr, scales = "free_x", ncol = 6) +
  scale_y_continuous(name = "Number of recombination by Mb") +
  theme_minimal()

# All chromosome
ggplot(result, aes(x = position_Mbp, y = rate, color = Method)) +
  geom_line(alpha = 0.3) +
  stat_smooth(method = "loess", span = 0.1, se = FALSE, size = 1) +
  geom_hline(yintercept = length(common_indiv)/100, linetype = "dashed", color = "red") +
  annotate("text",
           x = -5, y = 5,
           label = paste0("Expected=", length(common_indiv)/100),
           vjust = -Inf, hjust = -Inf,
           color = "black", size = 3) +
  facet_wrap(~ chr, scales = "free_x", ncol = 6) +
  scale_y_continuous(name = "Number of recombination by Mb") +
  theme_minimal()



