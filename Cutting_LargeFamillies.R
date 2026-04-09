#Cut large families to supportable one by Merlin
#Characterized cutted ones (many number many genotyped one  statistics summary) nombre de generation


#"D:\copie PhD\PHD\Automne 2020\Activité de recherche\ARA20_3.7b.*"

library(tidyr)
library(FamAgg)

structure.famille=read.delim2("D:/DATA/WGS.ped", header = FALSE, sep="")
colnames(structure.famille)[1]="famille"
colnames(structure.famille)[2]="Individu"

f="GSF8540"

mbped= structure.famille[structure.famille$famille==f,]
mbped <- mbped[, 1:5]
colnames(mbped) <- c("family", "id", "father", "mother", "sex")

mbped$family=as.character(mbped$family)
mbped$id=as.character(mbped$id)

fadi=FAData(pedigree=mbped, header=FALSE, sep="\t", id.col="id",
            family.col="family", father.col="father", mother.col="mother",
            sex.col="sex")
#*******************
fond=findFounders(fadi, family = f)
chil2=getChildren(fadi, id=fond, max.generations=1)
kin2 <- kinship(fadi)
kin=as.matrix(kin2)
kin=as.data.frame(kin)
kin=kin[names(kin)%in%chil2]
kin=kin[row.names(kin)%in%chil2,]
for (i in 1:length(kin)) {
  
  for (j in 1:length(kin)) {
    if (i<=j) {
      kin[i,j]="."
    }
  }
}
listfrere=vector(mode="list")
for (i in 1:length(kin)) {
  tmp=kin[kin[i]==0.250,][i]
  if (length(row.names(tmp))!=0) {
    listfrere[[i]]= c(colnames(tmp),rownames(tmp))
  }
}
listfrere[[i+1]]=chil2[!chil2 %in% unlist(listfrere)]

listfrere=Filter(length, listfrere)
m=0
for (i in 1:length(listfrere)) {
  if (sum(lengths(listfrere[1:i]))>(length(chil2)/2)  ) {
    m=i-1
    break
  }
}

f1=buildPed(fadi, id=c(unlist(listfrere[1:m]),getChildren(fadi, id=unlist(listfrere[1:m]), max.generations=1)) ,prune=TRUE, cex=0.5)

f2=buildPed(fadi, id=c(unlist(listfrere[m+1:length(listfrere)]),getChildren(fadi, id=unlist(listfrere[m+1:length(listfrere)]), max.generations=1)) ,prune=TRUE, cex=0.5)
f1$"familyAD"=paste0(f1$family,"_1") 
f2$"familyAD"=paste0(f2$family,"_2")
fAD=rbind(f1,f2)
fAD=cbind("familyAD"=fAD[,6],fAD[,1:5])
row.names(fAD)=NULL
fAD

library(tidyr)
library(FamAgg)

structure.famille=read.delim2("D:/DATA/WGS.ped", header = FALSE, sep="")
colnames(structure.famille)[1]="famille"
colnames(structure.famille)[2]="Individu"


chr.shapeit2=read.delim2(paste0("D:/DATA/Data14092020/chr1.OmniExpress.clean.recombinations"), header = TRUE, sep="")
famille.shapeit2.chr=structure.famille[structure.famille$Individu %in% as.character(unique(chr.shapeit2$CHILD)),]
famille.shapeit2.chr=as.character(unique(famille.shapeit2.chr$famille))

famille_decouper=NULL
for (f in famille.shapeit2.chr[! famille.shapeit2.chr%in%c("GSF2010", "GSF1046","GSF8540","GSF0536","GSF3511","GSF6397", "GSF7670","GSF5640","GSF8978","GSF3160","GSF9249","GSF2741","GSF3028","GSF5519","GSF9696","GSF7966")]) {
  mbped= structure.famille[structure.famille$famille==f,]
  mbped <- mbped[, 1:5]
  colnames(mbped) <- c("family", "id", "father", "mother", "sex")
  
  mbped$family=as.character(mbped$family)
  mbped$id=as.character(mbped$id)
  fadi=FAData(pedigree=mbped, header=FALSE, sep="\t", id.col="id",
              family.col="family", father.col="father", mother.col="mother",
              sex.col="sex")
  fond=findFounders(fadi, family = f)
  chil2=getChildren(fadi, id=fond, max.generations=1)
  kin2 <- kinship(fadi)
  kin=as.matrix(kin2)
  kin=as.data.frame(kin)
  kin=kin[names(kin)%in%chil2]
  kin=kin[row.names(kin)%in%chil2,]
  for (i in 1:length(kin)) {
    for (j in 1:length(kin)) {
      if (i<=j) {
        kin[i,j]="."
      }
    }
  }
  listfrere=vector(mode="list")
  tmp=NULL
  for (i in 1:length(kin)) {
    tmp=kin[kin[i]==0.250,][i]
    if (length(row.names(tmp))!=0) {
      listfrere[[i]]= c(colnames(tmp),rownames(tmp))
    }
  }
  listfrere[[i+1]]=chil2[!chil2 %in% unlist(listfrere)]
  
  listfrere=Filter(length, listfrere)
  listfrere=listfrere[order(sapply(listfrere, length))]
  m=0
  for (i in 1:length(listfrere)) {
    if (length(unique(unlist(listfrere[1:i])))>=(length(chil2)/2)  ) {
      m=length(unique(unlist(listfrere[1:i])))-1
      break
    }
  }
  f1=buildPed(fadi, id=c(unlist(listfrere[1:m]),getChildren(fadi, id=unlist(listfrere[1:m]), max.generations=1)) ,prune=TRUE, cex=0.5)
  f2=buildPed(fadi, id=c(chil2[!chil2 %in% unlist(listfrere[1:m])],getChildren(fadi, id=chil2[!chil2 %in% unlist(listfrere[1:m])], max.generations=1)) ,prune=TRUE, cex=0.5)
  
  f1$"familyAD"=paste0(f1$family,"_1") 
  f2$"familyAD"=paste0(f2$family,"_2")
  fAD=rbind(f1,f2)
  fAD=cbind("familyAD"=fAD[,6],fAD[,1:5])
  famille_decouper=rbind(famille_decouper,fAD)
}
length(famille.shapeit2.chr[! famille.shapeit2.chr%in%c("GSF2010", "GSF1046","GSF8540","GSF0536","GSF3511","GSF6397", "GSF7670","GSF5640","GSF8978","GSF3160","GSF9249","GSF2741","GSF3028","GSF5519","GSF9696","GSF7966")])
length(unique(famille_decouper$family))
print("**************************************************************************************")

# "GSF5640","GSF8978" non exécuté
# ,"GSF3160","GSF9249"
#,"GSF2741","GSF3028","GSF5519","GSF9696","GSF7966" speciale

for (f in c("GSF1046","GSF8540","GSF0536","GSF3511","GSF6397","GSF7670")) {
  mbped= structure.famille[structure.famille$famille==f,]
  mbped <- mbped[, 1:5]
  colnames(mbped) <- c("family", "id", "father", "mother", "sex")
  
  mbped$family=as.character(mbped$family)
  mbped$id=as.character(mbped$id)
  fadi=FAData(pedigree=mbped, header=FALSE, sep="\t", id.col="id",
              family.col="family", father.col="father", mother.col="mother",
              sex.col="sex")
  fond=findFounders(fadi, family = f)
  chil2=getChildren(fadi, id=fond, max.generations=1)
  kin2 <- kinship(fadi)
  kin=as.matrix(kin2)
  kin=as.data.frame(kin)
  kin=kin[names(kin)%in%chil2]
  kin=kin[row.names(kin)%in%chil2,]
  for (i in 1:length(kin)) {
    
    for (j in 1:length(kin)) {
      if (i<=j) {
        kin[i,j]="."
      }
    }
  }
  listfrere=vector(mode="list")
  tmp=NULL
  for (i in 1:length(kin)) {
    tmp=kin[kin[i]==0.250,][i]
    if (length(row.names(tmp))!=0) {
      listfrere[[i]]= c(colnames(tmp),rownames(tmp))
    }
  }
  listfrere[[i+1]]=chil2[!chil2 %in% unlist(listfrere)]
  
  listfrere=Filter(length, listfrere)
  listfrere=listfrere[order(sapply(listfrere, length))]
  m=0
  for (i in 1:length(listfrere)) {
    if (length(unique(unlist(listfrere[1:i])))>=(length(chil2)/2)  ) {
      m=i-1
      break
    }
  }
  f1=buildPed(fadi, id=c(unlist(listfrere[1:m]),getChildren(fadi, id=unlist(listfrere[1:m]), max.generations=1)) ,prune=TRUE, cex=0.5)
  f2=buildPed(fadi, id=c(chil2[!chil2 %in% unlist(listfrere[1:m])],getChildren(fadi, id=chil2[!chil2 %in% unlist(listfrere[1:m])], max.generations=1)) ,prune=TRUE, cex=0.5)
  f1$"familyAD"=paste0(f1$family,"_1") 
  f2$"familyAD"=paste0(f2$family,"_2")
  fAD=rbind(f1,f2)
  fAD=cbind("familyAD"=fAD[,6],fAD[,1:5])
  famille_decouper=rbind(famille_decouper,fAD)
}
length(unique(famille_decouper$family))
print("**************************************************************************************")

# #,"GSF2741","GSF3028","GSF5519","GSF9696" ,"GSF7966" speciale
#c("GSF5640","GSF8978","GSF3160","GSF9249")

for (f in c("GSF8978","GSF3160")) {
  mbped= structure.famille[structure.famille$famille==f,]
  mbped <- mbped[, 1:5]
  colnames(mbped) <- c("family", "id", "father", "mother", "sex")
  
  mbped$family=as.character(mbped$family)
  mbped$id=as.character(mbped$id)
  fadi=FAData(pedigree=mbped, header=FALSE, sep="\t", id.col="id",
              family.col="family", father.col="father", mother.col="mother",
              sex.col="sex")
  fond=findFounders(fadi, family = f)
  chil2=getChildren(fadi, id=fond, max.generations=1)
  
  f1=buildPed(fadi, id=c(chil2[1:as.integer(length(chil2)/2)],getChildren(fadi, id=chil2[1:as.integer(length(chil2)/2)], max.generations=1), fond) ,prune=TRUE, cex=0.5) 
  f2=buildPed(fadi, id=c(chil2[(as.integer(length(chil2)/2)+1):length(chil2)],getChildren(fadi, id=chil2[(as.integer(length(chil2)/2)+1):length(chil2)], max.generations=1)) ,prune=TRUE, cex=0.5)
  
  f1$"familyAD"=paste0(f1$family,"_1") 
  f2$"familyAD"=paste0(f2$family,"_2")
  fAD=rbind(f1,f2)
  fAD=cbind("familyAD"=fAD[,6],fAD[,1:5])
  famille_decouper=rbind(famille_decouper,fAD)
}
length(unique(famille_decouper$family))
print("**************************************************************************************")


# #,"GSF2741","GSF3028","GSF5519","GSF9696" ,"GSF7966" speciale
#c("GSF5640","GSF8978","GSF3160","GSF9249")
#,"GSF9249"

f="GSF5640"
mbped= structure.famille[structure.famille$famille==f,]
mbped <- mbped[, 1:5]
colnames(mbped) <- c("family", "id", "father", "mother", "sex")

mbped$family=as.character(mbped$family)
mbped$id=as.character(mbped$id)
fadi=FAData(pedigree=mbped, header=FALSE, sep="\t", id.col="id",
            family.col="family", father.col="father", mother.col="mother",
            sex.col="sex")
fond=findFounders(fadi, family = f)
chil2=getChildren(fadi, id=fond, max.generations=1)
fond2=mbped[(mbped$father==0)&(mbped$mother==0),]$id
fond2=fond2[!fond2 %in% fond]
f1=buildPed(fadi, id=c(fond2, getChildren(fadi, id=fond2, max.generations=1)) ,prune=TRUE, cex=0.5) 
f2=buildPed(fadi, id=c(chil2, getChildren(fadi, id=chil2, max.generations=1)) ,prune=TRUE, cex=0.5) 

f1$"familyAD"=paste0(f1$family,"_1") 
f2$"familyAD"=paste0(f2$family,"_2")
fAD=rbind(f1,f2)
fAD=cbind("familyAD"=fAD[,6],fAD[,1:5])
famille_decouper=rbind(famille_decouper,fAD)
length(unique(famille_decouper$family))

print("**************************************************************************************")



# #,"GSF2741","GSF3028","GSF5519","GSF9696" ,"GSF7966" speciale

f="GSF9249"
mbped= structure.famille[structure.famille$famille==f,]
mbped <- mbped[, 1:5]
colnames(mbped) <- c("family", "id", "father", "mother", "sex")

mbped$family=as.character(mbped$family)
mbped$id=as.character(mbped$id)
fadi=FAData(pedigree=mbped, header=FALSE, sep="\t", id.col="id",
            family.col="family", father.col="father", mother.col="mother",
            sex.col="sex")
fond=findFounders(fadi, family = f)
chil2=getChildren(fadi, id=fond, max.generations=1)
fond2=mbped[(mbped$father==0)&(mbped$mother==0),]$id
fond2=fond2[!fond2 %in% fond]

f1=buildPed(fadi, id=c(fond2[c(4,5)], getChildren(fadi, id=fond2[c(4,5)], max.generations=1),getChildren(fadi, id=getChildren(fadi, id=fond2[c(4,5)], max.generations=1), max.generations=1)) ,prune=TRUE, cex=0.5)
f2=buildPed(fadi, id=c(chil2[-1], getChildren(fadi, id=chil2[-1], max.generations=1)) ,prune=TRUE, cex=0.5)

#plotPed(fadi, family = f, cex=0.5)

f1$"familyAD"=paste0(f1$family,"_1") 
f2$"familyAD"=paste0(f2$family,"_2")
fAD=rbind(f1,f2)
fAD=cbind("familyAD"=fAD[,6],fAD[,1:5])
famille_decouper=rbind(famille_decouper,fAD)
length(unique(famille_decouper$family))


# #,"GSF3028","GSF5519","GSF9696" ,"GSF7966" speciale

f="GSF2741"
mbped= structure.famille[structure.famille$famille==f,]
mbped <- mbped[, 1:5]
colnames(mbped) <- c("family", "id", "father", "mother", "sex")

mbped$family=as.character(mbped$family)
mbped$id=as.character(mbped$id)
fadi=FAData(pedigree=mbped, header=FALSE, sep="\t", id.col="id",
            family.col="family", father.col="father", mother.col="mother",
            sex.col="sex")
fond=findFounders(fadi, family = f)
chil2=getChildren(fadi, id=fond, max.generations=1)
fond2=mbped[(mbped$father==0)&(mbped$mother==0),]$id
fond2=fond2[!fond2 %in% fond]

f1=buildPed(fadi, id=c(chil2[3:4], getChildren(fadi, id=chil2[3:4], max.generations=1)) ,prune=TRUE, cex=0.5) 
f2=buildPed(fadi, id=c(chil2[1:2],"GS47871020", getChildren(fadi, id=c(chil2[1:2],"GS47871020"), max.generations=1)) ,prune=TRUE, cex=0.5)  
f1$"familyAD"=paste0(f1$family,"_1") 
f2$"familyAD"=paste0(f2$family,"_2")
fAD=rbind(f1,f2)
fAD=cbind("familyAD"=fAD[,6],fAD[,1:5])
famille_decouper=rbind(famille_decouper,fAD)
length(unique(famille_decouper$family))

print("**************************************************************************************")



# #"GSF5519","GSF9696" ,"GSF7966" speciale

f="GSF3028"
mbped= structure.famille[structure.famille$famille==f,]
mbped <- mbped[, 1:5]
colnames(mbped) <- c("family", "id", "father", "mother", "sex")

mbped$family=as.character(mbped$family)
mbped$id=as.character(mbped$id)
fadi=FAData(pedigree=mbped, header=FALSE, sep="\t", id.col="id",
            family.col="family", father.col="father", mother.col="mother",
            sex.col="sex")
fond=findFounders(fadi, family = f)
chil2=getChildren(fadi, id=fond, max.generations=1)
fond2=mbped[(mbped$father==0)&(mbped$mother==0),]$id
fond2=fond2[!fond2 %in% fond]
f1=buildPed(fadi, id=c(chil2[2:length(chil2)], getChildren(fadi, id=chil2[2:length(chil2)], max.generations=1)) ,prune=TRUE, cex=0.5)
f2=buildPed(fadi, id=c(fond2[7:length(fond2)], getChildren(fadi, id=fond2[7:length(fond2)], max.generations=1)) ,prune=TRUE, cex=0.5)

f1$"familyAD"=paste0(f1$family,"_1") 
f2$"familyAD"=paste0(f2$family,"_2")
fAD=rbind(f1,f2)
fAD=cbind("familyAD"=fAD[,6],fAD[,1:5])
famille_decouper=rbind(famille_decouper,fAD)
length(unique(famille_decouper$family))

print("**************************************************************************************")


# #,"GSF9696" ,"GSF7966" speciale

f="GSF5519"
mbped= structure.famille[structure.famille$famille==f,]
mbped <- mbped[, 1:5]
colnames(mbped) <- c("family", "id", "father", "mother", "sex")

mbped$family=as.character(mbped$family)
mbped$id=as.character(mbped$id)
fadi=FAData(pedigree=mbped, header=FALSE, sep="\t", id.col="id",
            family.col="family", father.col="father", mother.col="mother",
            sex.col="sex")

fond=findFounders(fadi, family = f)
chil2=getChildren(fadi, id=fond, max.generations=1)
fond2=mbped[(mbped$father==0)&(mbped$mother==0),]$id
fond2=fond2[!fond2 %in% fond]

f1=buildPed(fadi, id=c(chil2[1:2], getChildren(fadi, id=chil2[1:2], max.generations=1)) ,prune=TRUE, cex=0.5)
f2=buildPed(fadi, id=c(fond2[c(1,3:5)], getChildren(fadi, id=fond2[c(1,3:5)], max.generations=1)) ,prune=TRUE, cex=0.5)

f1$"familyAD"=paste0(f1$family,"_1") 
f2$"familyAD"=paste0(f2$family,"_2")
fAD=rbind(f1,f2)
fAD=cbind("familyAD"=fAD[,6],fAD[,1:5])
famille_decouper=rbind(famille_decouper,fAD)
length(unique(famille_decouper$family))

print("**************************************************************************************")


# #"GSF7966" speciale

f="GSF9696"
mbped= structure.famille[structure.famille$famille==f,]
mbped <- mbped[, 1:5]
colnames(mbped) <- c("family", "id", "father", "mother", "sex")

mbped$family=as.character(mbped$family)
mbped$id=as.character(mbped$id)
fadi=FAData(pedigree=mbped, header=FALSE, sep="\t", id.col="id",
            family.col="family", father.col="father", mother.col="mother",
            sex.col="sex")
fond=findFounders(fadi, family = f)
chil2=getChildren(fadi, id=fond, max.generations=1)
fond2=mbped[(mbped$father==0)&(mbped$mother==0),]$id
fond2=fond2[!fond2 %in% fond]

f1=buildPed(fadi, id=c(fond2[1:3], getChildren(fadi, id=fond2[1:3], max.generations=1)) ,prune=TRUE, cex=0.5)
f2=buildPed(fadi, id=c(chil2[-c(1,5)], getChildren(fadi, id=chil2[-c(1,5)], max.generations=1)) ,prune=TRUE, cex=0.5)

f1$"familyAD"=paste0(f1$family,"_1") 
f2$"familyAD"=paste0(f2$family,"_2")
fAD=rbind(f1,f2)
fAD=cbind("familyAD"=fAD[,6],fAD[,1:5])
famille_decouper=rbind(famille_decouper,fAD)
length(unique(famille_decouper$family))

print("**************************************************************************************")

f="GSF7966"
mbped= structure.famille[structure.famille$famille==f,]
mbped <- mbped[, 1:5]
colnames(mbped) <- c("family", "id", "father", "mother", "sex")

mbped$family=as.character(mbped$family)
mbped$id=as.character(mbped$id)
fadi=FAData(pedigree=mbped, header=FALSE, sep="\t", id.col="id",
            family.col="family", father.col="father", mother.col="mother",
            sex.col="sex")

fond=findFounders(fadi, family = f)
chil2=getChildren(fadi, id=fond, max.generations=1)
fond2=mbped[(mbped$father==0)&(mbped$mother==0),]$id
fond2=fond2[!fond2 %in% fond]

f1=buildPed(fadi, id=c(chil2[1],fond, getChildren(fadi, id=chil2[1], max.generations=1)) ,prune=TRUE, cex=0.5)
f2=buildPed(fadi, id=c(fond2[-c(1,3)],chil2[c(2,3)], getChildren(fadi, id=c(fond2[-c(1,3)],chil2[c(2,3)]), max.generations=1)) ,prune=TRUE, cex=0.5)

f1$"familyAD"=paste0(f1$family,"_1") 
f2$"familyAD"=paste0(f2$family,"_2")
fAD=rbind(f1,f2)
fAD=cbind("familyAD"=fAD[,6],fAD[,1:5])
famille_decouper=rbind(famille_decouper,fAD)
length(unique(famille_decouper$family))

print("**************************************************************************************")


f="GSF2010"
mbped= structure.famille[structure.famille$famille==f,]
mbped <- mbped[, 1:5]
colnames(mbped) <- c("family", "id", "father", "mother", "sex")

mbped$family=as.character(mbped$family)
mbped$id=as.character(mbped$id)
fadi=FAData(pedigree=mbped, header=FALSE, sep="\t", id.col="id",
            family.col="family", father.col="father", mother.col="mother",
            sex.col="sex")

fond=findFounders(fadi, family = f)
chil2=getChildren(fadi, id=fond, max.generations=1)
fond2=mbped[(mbped$father==0)&(mbped$mother==0),]$id
fond2=fond2[!fond2 %in% fond]

f1=buildPed(fadi, id=c(chil2[-3], getChildren(fadi, id=chil2[-3], max.generations=1)) ,prune=TRUE, cex=0.5)
f2=buildPed(fadi, id=c(fond2[3:4], getChildren(fadi, id=fond2[3:4], max.generations=1)) ,prune=TRUE, cex=0.5)

f1$"familyAD"=paste0(f1$family,"_1") 
f2$"familyAD"=paste0(f2$family,"_2")
fAD=rbind(f1,f2)
fAD=cbind("familyAD"=fAD[,6],fAD[,1:5])
famille_decouper=rbind(famille_decouper,fAD)
length(unique(famille_decouper$family))

print("**************************************************************************************")
famille.shapeit2.chr %in% unique(famille_decouper$family)
row.names(famille_decouper)=NULL
famille_decouper = subset(famille_decouper, select = -c(2) )

#Vérifier si merlin fait l'haplotypage de sous-famille f1 (la plus grande) 



t=as.data.frame(table(famille_decouper$familyAD))
t[order(-t$Freq),]
PGSF=as.character(t[t$Freq==max(t$Freq),]$Var1)

f1.ped= famille_decouper[famille_decouper$familyAD==PGSF,]
h=rbind (c("2/ 1",   "2/ 1",   "2/ 1"), c("1/ 1",   "1/ 1",   "1/ 1" ),c("2/ 1",   "2/ 1",   "2/ 1"  ),c("2/ 1",   "2/ 1",   "2/ 1"  ))
f1.ped=cbind(f1.ped,rbind(h,h,h,h,h,h[1,]))

# fichier map
f.dat=rbind(c("M",  "1:841166:A:G"), c("M" , "1:911428:C:T"),c("M" , "1:968785:A:G"))

write.table(f1.ped, file = "D:/DATA/f1.ped", sep = " ", quote=F, row.names=F, col.names=F)
write.table(f.dat, file = "D:/DATA/f.dat", sep = " ", quote=F, row.names=F, col.names=F)
#Commande Merlin :

  #MERLIN -p f1.ped -d f.dat -m merlin_map.01  --best --horizontal --smallSwap
  #Analyse merlin et bien faite


f.ped= structure.famille[structure.famille$famille=="GSF3028",]
f.ped=f.ped[,1:5]
f.ped=cbind(f.ped,rbind(h,h,h,h,h,h,h))
write.table(f.ped, file = "D:/DATA/f.ped", sep = " ", quote=F, row.names=F, col.names=F)

#elaboration de fichier.ped des sousfamilles
sous_familles=cbind(famille_decouper, c6=rep(-9,length(famille_decouper$id)),c7=rep(0,length(famille_decouper$id)),c8=rep(0,length(famille_decouper$id)))
sous_familles[is.na(sous_familles)]=0

write.table(sous_familles, file = "D:/DATA/sous_familles.ped", sep = " ", quote=F, row.names=F, col.names=F)



f="GSF9249"
mbped=structure.famille[structure.famille$famille==f,]
mbped <- mbped[, 1:5]
colnames(mbped) <- c("family", "id", "father", "mother", "sex")

mbped$family=as.character(mbped$family)
mbped$id=as.character(mbped$id)
fadi=FAData(pedigree=mbped, header=FALSE, sep="\t", id.col="id",
       family.col="family", father.col="father", mother.col="mother",
       sex.col="sex")
plotPed(fadi,family=f, cex=0.5)

f="GSF9249_1"
mbped= sous_familles[sous_familles$family==f,]
mbped <- mbped[, 1:5]
colnames(mbped) <- c("family", "id", "father", "mother", "sex")

mbped$family=as.character(mbped$family)
mbped$id=as.character(mbped$id)
fadi=FAData(pedigree=mbped, header=FALSE, sep="\t", id.col="id",
       family.col="family", father.col="father", mother.col="mother",
       sex.col="sex")
plotPed(fadi,family=f, cex=0.5)

f="GSF9249_2"
mbped= sous_familles[sous_familles$family==f,]
mbped <- mbped[, 1:5]
colnames(mbped) <- c("family", "id", "father", "mother", "sex")

mbped$family=as.character(mbped$family)
mbped$id=as.character(mbped$id)
fadi=FAData(pedigree=mbped, header=FALSE, sep="\t", id.col="id",
       family.col="family", father.col="father", mother.col="mother",
       sex.col="sex")
plotPed(fadi,family=f, cex=0.5)

f1.ped=cbind(mbped,rbind(h,h,h,h,h[1:2,]))
write.table(f1.ped, file = "D:/DATA/f1.ped", sep = " ", quote=F, row.names=F, col.names=F)

