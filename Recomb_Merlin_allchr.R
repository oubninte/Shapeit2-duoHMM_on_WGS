rm(list=ls())

# Merlin recombination inference
library(dplyr)
library(stringr)
library(GenomicRanges)
library(BRGenomics)


#L'objectifs �largir les intervalles de merlin basant sur l'h�t�rozygocit� du parent
#update to include chr9 and path of data in rorqual
#ped1=read.table("/home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/SFTP hub/PHD -P1/Data/sous_familles.ped", header = FALSE)
#Dans ce programme  : Localisation des recombinaison


# paths to data
path_rep="/lustre09/project/6033529/genealogy_sims/results/Samir/SFTP hub/PHD -P1/Data/DataMerlin/hapAD/"
chemin_map = c( paste0(path_rep, "/chr", c(1:9),"/merlin_map.0",c(1:9)),
            paste0(path_rep, "/chr", c(10:22), "/merlin_map.", c(10:22)))
cheminMF = paste0(path_rep, "/chr", c(1:9,10:22),"/merlin.flow")
chrNum = paste0("chr", c(1:9,10:22))
cheminHet=paste0(path_rep, "het.",c(1:9,10:22),".RData")

#j=21
#loading pedigree
ped=read.table("/lustre09/project/6033529/genealogy_sims/results/Samir/P1/Candiate_recombinations/WGS.ped")

for (j in c(9)) { #others are already runned
  
  #check if the number of markers matches the map file => yes
  
  heterozygous=get(load(cheminHet[j]))
  colnames(heterozygous)[1:4]=c("Family", "Indiv","Father","Mother")
  
  #ncol(heterozygous)-6
  merlin_map = read.delim2(chemin_map[j], header = T, sep = "")
  
  #nrow(merlin_map)
  
  cat("le nombre de marqueurs correspond avec le fichier map :chr", j)
  if (nrow(merlin_map)==ncol(heterozygous)-6) {
    print("Oui")
  } else {print("Non")
    break
    }
  
  #Marker physical position in genetic map
  merlin_map = cbind(merlin_map, "PHY.POS"=as.numeric(str_split_fixed(merlin_map$MARKER, ":", 3)[,2]))
  
  # data processing : preparting to infer recombinaiton
  chr = read.delim2(cheminMF[j], header = FALSE, sep = "")
  chr = cbind("Family" = chr$V2, chr)
  chr$Family[chr$Family %in% c("(FOUNDER)", "(MATERNAL)", "(PATERNAL)")] = NA
  chr$Family = zoo::na.locf(chr$Family)
  colnames(chr)[2] = "ID"
  chr = chr[chr$V5 != "",]
  
  #Locatinzation of recombination (identify the index of marker)
  Recombinaison <- list()
  chr=chr[chr$V2!="(FOUNDER)",]
  for (i in 1:nrow(chr)) {
    rl <- rle(chr[i, ])
    if (length(rl$values) - 4 != 0) {
      tmp <- data.frame("CHILD" = chr[i, 2], 
                        "PARENT"=chr[i, 3],
                        "Parent_ID"=if_else(chr[i, 3]=="(PATERNAL)", ped$V3[ped$V2==chr[i, 2]] , ped$V4[ped$V2==chr[i, 2]]),
                        "PosiRecombStart" = cumsum(rl$lengths[4:length(rl$values)]))
      Recombinaison[[i]] <- tmp
    }
  }
  
  Recombinaison <- do.call(rbind, Recombinaison)
  Recombinaison = Recombinaison[Recombinaison$PosiRecombStart != ncol(chr)-3,] 
  Recombinaison = Recombinaison[order(Recombinaison$PosiRecombStart),]
  
  
  #broading recombination interval to heterozygote SNP of the parent
  
  # if the parent if is herozygote
  heterozygous$invariant_rows <- apply(heterozygous[,7:ncol(heterozygous)], 1, function(row) length(unique(row)) == 1)
  Recombinaison_hetero=Recombinaison[!Recombinaison$Parent_ID %in% heterozygous$Indiv[heterozygous$invariant_rows],]
  nc=ncol(Recombinaison_hetero)
  merlin_map$"SNP_ID"=1:nrow(merlin_map)
  for (i in 1:nrow(Recombinaison_hetero)) {
    hetero_status=as.logical(heterozygous[heterozygous$Indiv==Recombinaison_hetero$Parent_ID[i],][7:(ncol(heterozygous)-1)])
    tmp=merlin_map[hetero_status,]
    tmp$"end"=c(tmp$PHY.POS[-1], NA)

    if (sum(tmp$SNP_ID<=Recombinaison_hetero$PosiRecombStart[i])==0) {
      tmp=tmp[tmp$SNP_ID==min(tmp$SNP_ID[tmp$SNP_ID>=Recombinaison_hetero$PosiRecombStart[i]]),]
      #The start will be the same given by merlin, the end will be the first heterozygot SNP after the start (it could be the next one)
      Recombinaison_hetero[i,(nc+1):(nc+2)] <- c(merlin_map$PHY.POS[merlin_map$SNP_ID==Recombinaison_hetero$PosiRecombStart[i]], tmp$PHY.POS) 
    }else{
      tmp=tmp[tmp$SNP_ID==max(tmp$SNP_ID[tmp$SNP_ID<=Recombinaison_hetero$PosiRecombStart[i]]),]
      Recombinaison_hetero[i,(nc+1):(nc+2)] <- c(tmp$PHY.POS, tmp$end)
    }
  }
  colnames(Recombinaison_hetero)[(nc+1):(nc+2)]=c("START", "END")
  Recombinaison_hetero$"isParentHeteroz"=TRUE
  Recombinaison_hetero=Recombinaison_hetero[c(1,3,5:7)]

  
  #Recombinaison_hetero_newMethod=Recombinaison_hetero
  # if the parent if is NOT herozygote
  Recombinaison_homo=Recombinaison[Recombinaison$Parent_ID %in% heterozygous$Indiv[heterozygous$invariant_rows],]
  tmp=merlin_map
  tmp$"end"=c(tmp$PHY.POS[-1], NA)
  Recombinaison_homo=merge(Recombinaison_homo, tmp, by.x="PosiRecombStart", by.y="SNP_ID", all.x=T, all.y=F, sort=F)
  Recombinaison_homo$"isParentHeteroz"=FALSE
  Recombinaison_homo=Recombinaison_homo[c(2,4,8:10)]

  colnames(Recombinaison_homo)=colnames(Recombinaison_hetero)= c("CHILD", "PARENT", "START", "END", "isParentHeteroz" )

  recombinaison=rbind(Recombinaison_homo,Recombinaison_hetero)
  write.table(recombinaison,paste0("/lustre09/project/6033529/genealogy_sims/results/Samir/P1/paper3_Analyzes/Merlin_Recomb/",chrNum[j], ".merlin.recombinaisons"),quote=F,row.names=F)
  
}
