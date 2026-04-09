# Cut large families into smaller sub-families that are manageable for Merlin linkage analysis

# Load required libraries for pedigree analysis
library(tidyr)
library(FamAgg)

# Load pedigree data 
structure.famille=read.delim2("D:/DATA/WGS.ped", header = FALSE, sep="")

# Rename columns to standard pedigree format
colnames(structure.famille)[1]="famille"
colnames(structure.famille)[2]="Individu"

# ========== INITIAL EXAMPLE: Split family GSF8540 ==========
f="GSF8540"

# Extract pedigree data for the target family
mbped= structure.famille[structure.famille$famille==f,]

# Keep only first 5 columns (family, id, father, mother, sex)
mbped <- mbped[, 1:5]
colnames(mbped) <- c("family", "id", "father", "mother", "sex")

# Convert family and individual IDs to character type
mbped$family=as.character(mbped$family)
mbped$id=as.character(mbped$id)

# Create FAData object for pedigree analysis using FamAgg package
fadi=FAData(pedigree=mbped, header=FALSE, sep="\t", id.col="id",
            family.col="family", father.col="father", mother.col="mother",
            sex.col="sex")

# Find founder individuals 
fond=findFounders(fadi, family = f)

# Get direct children of founders 
chil2=getChildren(fadi, id=fond, max.generations=1)

# Calculate kinship coefficients between all individuals
kin2 <- kinship(fadi)
kin=as.matrix(kin2)
kin=as.data.frame(kin)

# Filter kinship matrix to only include children of founders
kin=kin[names(kin)%in%chil2]
kin=kin[row.names(kin)%in%chil2,]

# Mask upper triangle of kinship matrix 
for (i in 1:length(kin)) {
  for (j in 1:length(kin)) {
    if (i<=j) {
      kin[i,j]="."
    }
  }
}

# Create list to group siblings based on kinship coefficient of 0.250 
listfrere=vector(mode="list")
for (i in 1:length(kin)) {
  tmp=kin[kin[i]==0.250,][i]
  if (length(row.names(tmp))!=0) {
    listfrere[[i]]= c(colnames(tmp),rownames(tmp))
  }
}
listfrere[[i+1]]=chil2[!chil2 %in% unlist(listfrere)]

# Remove empty list elements
listfrere=Filter(length, listfrere)
# Find split point
m=0
for (i in 1:length(listfrere)) {
  if (sum(lengths(listfrere[1:i]))>(length(chil2)/2)  ) {
    m=i-1
    break
  }
}

# Build first sub-family with founders and their children
f1=buildPed(fadi, id=c(unlist(listfrere[1:m]),getChildren(fadi, id=unlist(listfrere[1:m]), max.generations=1)) ,prune=TRUE, cex=0.5)

# Build second sub-family with remaining children
f2=buildPed(fadi, id=c(unlist(listfrere[m+1:length(listfrere)]),getChildren(fadi, id=unlist(listfrere[m+1:length(listfrere)]), max.generations=1)) ,prune=TRUE, cex=0.5)

# Add subfamily identifiers
f1$"familyAD"=paste0(f1$family,"_1") 
f2$"familyAD"=paste0(f2$family,"_2")

# Combine both subfamilies
fAD=rbind(f1,f2)

# Reorder columns 
fAD=cbind("familyAD"=fAD[,6],fAD[,1:5])
row.names(fAD)=NULL
fAD

# ========== MAIN PROCESSING: Automated family splitting (case 1) ==========
# Recall : Load libraries and data
library(tidyr)
library(FamAgg)

structure.famille=read.delim2("D:/DATA/WGS.ped", header = FALSE, sep="")
colnames(structure.famille)[1]="famille"
colnames(structure.famille)[2]="Individu"

# Families that have individuals in shapeit2 data
chr.shapeit2=read.delim2(paste0("D:/DATA/Data14092020/chr1.OmniExpress.clean.recombinations"), header = TRUE, sep="")
famille.shapeit2.chr=structure.famille[structure.famille$Individu %in% as.character(unique(chr.shapeit2$CHILD)),]
famille.shapeit2.chr=as.character(unique(famille.shapeit2.chr$famille))

# Initialize dataframe to accumulate split families
famille_decouper=NULL

for (f in famille.shapeit2.chr[! famille.shapeit2.chr%in%c("GSF2010", "GSF1046","GSF8540","GSF0536","GSF3511","GSF6397", "GSF7670","GSF5640","GSF8978","GSF3160","GSF9249","GSF2741","GSF3028","GSF5519","GSF9696","GSF7966")]) {
  mbped= structure.famille[structure.famille$famille==f,]
  mbped <- mbped[, 1:5]
  colnames(mbped) <- c("family", "id", "father", "mother", "sex")
  
  mbped$family=as.character(mbped$family)
  mbped$id=as.character(mbped$id)
  # Create pedigree object
  fadi=FAData(pedigree=mbped, header=FALSE, sep="\t", id.col="id",
              family.col="family", father.col="father", mother.col="mother",
              sex.col="sex")
  # Find founders and their children
  fond=findFounders(fadi, family = f)
  chil2=getChildren(fadi, id=fond, max.generations=1)
  # Calculate kinship matrix
  kin2 <- kinship(fadi)
  kin=as.matrix(kin2)
  kin=as.data.frame(kin)
  kin=kin[names(kin)%in%chil2]
  kin=kin[row.names(kin)%in%chil2,]
  # Mask upper triangle
  for (i in 1:length(kin)) {
    for (j in 1:length(kin)) {
      if (i<=j) {
        kin[i,j]="."
      }
    }
  }
  # Group siblings
  listfrere=vector(mode="list")
  tmp=NULL
  for (i in 1:length(kin)) {
    tmp=kin[kin[i]==0.250,][i]
    if (length(row.names(tmp))!=0) {
      listfrere[[i]]= c(colnames(tmp),rownames(tmp))
    }
  }
  listfrere[[i+1]]=chil2[!chil2 %in% unlist(listfrere)]
  
  # Remove empty lists and sort by size
  listfrere=Filter(length, listfrere)
  listfrere=listfrere[order(sapply(listfrere, length))]
  # Find split point to balance subfamily sizes
  m=0
  for (i in 1:length(listfrere)) {
    if (length(unique(unlist(listfrere[1:i])))>=(length(chil2)/2)  ) {
      m=length(unique(unlist(listfrere[1:i])))-1
      break
    }
  }
  # Build and label subfamilies
  f1=buildPed(fadi, id=c(unlist(listfrere[1:m]),getChildren(fadi, id=unlist(listfrere[1:m]), max.generations=1)) ,prune=TRUE, cex=0.5)
  f2=buildPed(fadi, id=c(chil2[!chil2 %in% unlist(listfrere[1:m])],getChildren(fadi, id=chil2[!chil2 %in% unlist(listfrere[1:m])], max.generations=1)) ,prune=TRUE, cex=0.5)
  
  f1$"familyAD"=paste0(f1$family,"_1") 
  f2$"familyAD"=paste0(f2$family,"_2")
  fAD=rbind(f1,f2)
  fAD=cbind("familyAD"=fAD[,6],fAD[,1:5])
  # Accumulate results
  famille_decouper=rbind(famille_decouper,fAD)
}
length(famille.shapeit2.chr[! famille.shapeit2.chr%in%c("GSF2010", "GSF1046","GSF8540","GSF0536","GSF3511","GSF6397", "GSF7670","GSF5640","GSF8978","GSF3160","GSF9249","GSF2741","GSF3028","GSF5519","GSF9696","GSF7966")])G[...]
# Report number of unique families processed
length(unique(famille_decouper$family))
print("**************************************************************************************")

# ========== MAIN PROCESSING: Automated family splitting (case 2) ==========
# Process families with specific splitting criteria
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
  # Mask upper triangle
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


# ========== MAIN PROCESSING: Automated family splitting (case 3) ==========

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

# ========== MANUAL HANDLING: Families requiring custom splitting ==========

# ========== SPECIAL CASE 1: GSF5640 ==========
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
# Identify additional founders not found by primary method
fond2=mbped[(mbped$father==0)&(mbped$mother==0),]$id
fond2=fond2[!fond2 %in% fond]
# First subfamily with secondary founders
f1=buildPed(fadi, id=c(fond2, getChildren(fadi, id=fond2, max.generations=1)) ,prune=TRUE, cex=0.5) 
# Second subfamily with primary children
f2=buildPed(fadi, id=c(chil2, getChildren(fadi, id=chil2, max.generations=1)) ,prune=TRUE, cex=0.5) 

f1$"familyAD"=paste0(f1$family,"_1") 
f2$"familyAD"=paste0(f2$family,"_2")
fAD=rbind(f1,f2)
fAD=cbind("familyAD"=fAD[,6],fAD[,1:5])
famille_decouper=rbind(famille_decouper,fAD)
length(unique(famille_decouper$family))

print("**************************************************************************************")

# ========== SPECIAL CASE 2: GSF9249 ) ==========
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

# Split using specific indices: selected founders and their descendants
f1=buildPed(fadi, id=c(fond2[c(4,5)], getChildren(fadi, id=fond2[c(4,5)], max.generations=1),getChildren(fadi, id=getChildren(fadi, id=fond2[c(4,5)], max.generations=1), max.generations=1)) ,prune=TRUE, cex=0.5)
f2=buildPed(fadi, id=c(chil2[-1], getChildren(fadi, id=chil2[-1], max.generations=1)) ,prune=TRUE, cex=0.5)

#plotPed(fadi, family = f, cex=0.5)
f1$"familyAD"=paste0(f1$family,"_1") 
f2$"familyAD"=paste0(f2$family,"_2")
fAD=rbind(f1,f2)
fAD=cbind("familyAD"=fAD[,6],fAD[,1:5])
famille_decouper=rbind(famille_decouper,fAD)
length(unique(famille_decouper$family))


# ========== SPECIAL CASE 3: GSF2741 ==========
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

# Split using specific child indices and a manual individual ID
f1=buildPed(fadi, id=c(chil2[3:4], getChildren(fadi, id=chil2[3:4], max.generations=1)) ,prune=TRUE, cex=0.5) 
f2=buildPed(fadi, id=c(chil2[1:2],"GS47871020", getChildren(fadi, id=c(chil2[1:2],"GS47871020"), max.generations=1)) ,prune=TRUE, cex=0.5)  
f1$"familyAD"=paste0(f1$family,"_1") 
f2$"familyAD"=paste0(f2$family,"_2")
fAD=rbind(f1,f2)
fAD=cbind("familyAD"=fAD[,6],fAD[,1:5])
famille_decouper=rbind(famille_decouper,fAD)
length(unique(famille_decouper$family))

print("**************************************************************************************")

# ========== SPECIAL CASE 4: GSF3028 ==========
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
# First subfamily
f1=buildPed(fadi, id=c(chil2[2:length(chil2)], getChildren(fadi, id=chil2[2:length(chil2)], max.generations=1)) ,prune=TRUE, cex=0.5)
# Second subfamily
f2=buildPed(fadi, id=c(fond2[7:length(fond2)], getChildren(fadi, id=fond2[7:length(fond2)], max.generations=1)) ,prune=TRUE, cex=0.5)

f1$"familyAD"=paste0(f1$family,"_1") 
f2$"familyAD"=paste0(f2$family,"_2")
fAD=rbind(f1,f2)
fAD=cbind("familyAD"=fAD[,6],fAD[,1:5])
famille_decouper=rbind(famille_decouper,fAD)
length(unique(famille_decouper$family))

print("**************************************************************************************")

# ========== SPECIAL CASE 5: GSF5519 ==========
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

# Split using specific child and founder indices
f1=buildPed(fadi, id=c(chil2[1:2], getChildren(fadi, id=chil2[1:2], max.generations=1)) ,prune=TRUE, cex=0.5)
f2=buildPed(fadi, id=c(fond2[c(1,3:5)], getChildren(fadi, id=fond2[c(1,3:5)], max.generations=1)) ,prune=TRUE, cex=0.5)

f1$"familyAD"=paste0(f1$family,"_1") 
f2$"familyAD"=paste0(f2$family,"_2")
fAD=rbind(f1,f2)
fAD=cbind("familyAD"=fAD[,6],fAD[,1:5])
famille_decouper=rbind(famille_decouper,fAD)
length(unique(famille_decouper$family))

print("**************************************************************************************")

# ========== SPECIAL CASE 6: GSF9696 ==========
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

# First subfamily
f1=buildPed(fadi, id=c(fond2[1:3], getChildren(fadi, id=fond2[1:3], max.generations=1)) ,prune=TRUE, cex=0.5)
# Second subfamily
f2=buildPed(fadi, id=c(chil2[-c(1,5)], getChildren(fadi, id=chil2[-c(1,5)], max.generations=1)) ,prune=TRUE, cex=0.5)

f1$"familyAD"=paste0(f1$family,"_1") 
f2$"familyAD"=paste0(f2$family,"_2")
fAD=rbind(f1,f2)
fAD=cbind("familyAD"=fAD[,6],fAD[,1:5])
famille_decouper=rbind(famille_decouper,fAD)
length(unique(famille_decouper$family))

print("**************************************************************************************")

# ========== SPECIAL CASE 7: GSF7966 ==========
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

# Split with specific indices for mixed founder and child combinations
f1=buildPed(fadi, id=c(chil2[1],fond, getChildren(fadi, id=chil2[1], max.generations=1)) ,prune=TRUE, cex=0.5)
f2=buildPed(fadi, id=c(fond2[-c(1,3)],chil2[c(2,3)], getChildren(fadi, id=c(fond2[-c(1,3)],chil2[c(2,3)]), max.generations=1)) ,prune=TRUE, cex=0.5)

f1$"familyAD"=paste0(f1$family,"_1") 
f2$"familyAD"=paste0(f2$family,"_2")
fAD=rbind(f1,f2)
fAD=cbind("familyAD"=fAD[,6],fAD[,1:5])
famille_decouper=rbind(famille_decouper,fAD)
length(unique(famille_decouper$family))

print("**************************************************************************************")

# ========== SPECIAL CASE 8: GSF2010 ==========
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

# Split on children and founder subsets
f1=buildPed(fadi, id=c(chil2[-3], getChildren(fadi, id=chil2[-3], max.generations=1)) ,prune=TRUE, cex=0.5)
f2=buildPed(fadi, id=c(fond2[3:4], getChildren(fadi, id=fond2[3:4], max.generations=1)) ,prune=TRUE, cex=0.5)

f1$"familyAD"=paste0(f1$family,"_1") 
f2$"familyAD"=paste0(f2$family,"_2")
fAD=rbind(f1,f2)
fAD=cbind("familyAD"=fAD[,6],fAD[,1:5])
famille_decouper=rbind(famille_decouper,fAD)
length(unique(famille_decouper$family))

print("**************************************************************************************")

# ========== FINALIZATION AND OUTPUT ==========
# Verify processing results
famille.shapeit2.chr %in% unique(famille_decouper$family)
# Reset row names for clean output
row.names(famille_decouper)=NULL
# Remove the individual ID column (column 2) from output
famille_decouper = subset(famille_decouper, select = -c(2) )

#Verify if Merlin performs haplotyping of first subfamily (largest) 
# Create frequency table of subfamily sizes
t=as.data.frame(table(famille_decouper$familyAD))
t[order(-t$Freq),]
# Identify the largest subfamily
PGSF=as.character(t[t$Freq==max(t$Freq),]$Var1)

# Extract the largest subfamily for Merlin analysis
f1.ped= famille_decouper[famille_decouper$familyAD==PGSF,]
# Create genotype matrix: rows are marker combinations, columns are 3 markers with specific genotypes
h=rbind (c("2/ 1",   "2/ 1",   "2/ 1"), c("1/ 1",   "1/ 1",   "1/ 1" ),c("2/ 1",   "2/ 1",   "2/ 1"  ),c("2/ 1",   "2/ 1",   "2/ 1"  ))
# Attach genotypes to pedigree (repeating pattern for all individuals)
f1.ped=cbind(f1.ped,rbind(h,h,h,h,h,h[1,]))

# Define marker map file: marker type and position/variant information
f.dat=rbind(c("M",  "1:841166:A:G"), c("M" , "1:911428:C:T"),c("M" , "1:968785:A:G"))

# Write pedigree and marker files for Merlin input
write.table(f1.ped, file = "D:/DATA/f1.ped", sep = " ", quote=F, row.names=F, col.names=F)
write.table(f.dat, file = "D:/DATA/f.dat", sep = " ", quote=F, row.names=F, col.names=F)
# Merlin command for linkage analysis:
  #MERLIN -p f1.ped -d f.dat -m merlin_map.01  --best --horizontal --smallSwap
  #Merlin analysis performed successfully

# Extract GSF3028 family from original data for analysis
f.ped= structure.famille[structure.famille$famille=="GSF3028",]
f.ped=f.ped[,1:5]
# Attach genotypes to GSF3028 pedigree
f.ped=cbind(f.ped,rbind(h,h,h,h,h,h,h))
write.table(f.ped, file = "D:/DATA/f.ped", sep = " ", quote=F, row.names=F, col.names=F)

# ========== CREATE OUTPUT PEDIGREE FILE FOR ALL SUBFAMILIES ==========
# Create pedigree file for all subfamilies with genotype and phenotype columns
sous_familles=cbind(famille_decouper, c6=rep(-9,length(famille_decouper$id)),c7=rep(0,length(famille_decouper$id)),c8=rep(0,length(famille_decouper$id)))
# Replace missing values with 0 (standard PED file format)
sous_familles[is.na(sous_familles)]=0

# Write complete subfamilies pedigree file
write.table(sous_familles, file = "D:/DATA/sous_familles.ped", sep = " ", quote=F, row.names=F, col.names=F)

# ========== VISUALIZATION OF SUBFAMILIES ==========
# Plot original GSF9249 family
f="GSF9249"
mbped=structure.famille[structure.famille$famille==f,]
mbped <- mbped[, 1:5]
colnames(mbped) <- c("family", "id", "father", "mother", "sex")

mbped$family=as.character(mbped$family)
mbped$id=as.character(mbped$id)
fadi=FAData(pedigree=mbped, header=FALSE, sep="\t", id.col="id",
       family.col="family", father.col="father", mother.col="mother",
       sex.col="sex")
# Display original family pedigree
plotPed(fadi,family=f, cex=0.5)

# Plot first subfamily of GSF9249
f="GSF9249_1"
mbped= sous_familles[sous_familles$family==f,]
mbped <- mbped[, 1:5]
colnames(mbped) <- c("family", "id", "father", "mother", "sex")

mbped$family=as.character(mbped$family)
mbped$id=as.character(mbped$id)
fadi=FAData(pedigree=mbped, header=FALSE, sep="\t", id.col="id",
       family.col="family", father.col="father", mother.col="mother",
       sex.col="sex")
# Display first split subfamily pedigree
plotPed(fadi,family=f, cex=0.5)

# Plot second subfamily of GSF9249
f="GSF9249_2"
mbped= sous_familles[sous_familles$family==f,]
mbped <- mbped[, 1:5]
colnames(mbped) <- c("family", "id", "father", "mother", "sex")

mbped$family=as.character(mbped$family)
mbped$id=as.character(mbped$id)
fadi=FAData(pedigree=mbped, header=FALSE, sep="\t", id.col="id",
       family.col="family", father.col="father", mother.col="mother",
       sex.col="sex")
# Display second split subfamily pedigree
plotPed(fadi,family=f, cex=0.5)


