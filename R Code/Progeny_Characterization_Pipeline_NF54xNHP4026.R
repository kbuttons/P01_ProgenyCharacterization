library("qtl") #load the package
library("vcfR")
library("igraph")
library("LinkageMapView")
library(qtlcharts)

### Setwd
setwd("/Users/kbuttons/Documents/Research/P01/Data/Core C/NF54xNHP4026/")

### read in file 
crossSNPs <- read.vcfR("NF54_NHP4026_combo.snps.recal.sel.vcf")
growthCross1 <- read.delim("Progeny_cross1_growth.csv", na.strings = c("", "NA"), sep=",",header=TRUE,as.is=TRUE)

### Setwd for saving files
setwd("/Users/kbuttons/Documents/Research/P01/Data/Core C/NF54xNHP4026/characterization/")

### Save Parent IDs
parent1 <- "Sample_NF54_GFP_LUC"
parent2 <- "Sample_ART-1"

### Set cross progeny IDs
cross1 <- colnames(crossSNPs@gt)[which(!colnames(crossSNPs@gt) %in% c("FORMAT",growthCross1[,1]))]
cross2 <- colnames(crossSNPs@gt)[which(colnames(crossSNPs@gt) %in% c(growthCross1[,1]))]


### Put these two columns into vectors 
parent1GT <- rep(-1, length(crossSNPs@gt[,parent1]))
parent2GT <- rep(-1, length(crossSNPs@gt[,parent2]))
parent1Coverage <- rep(-1, length(crossSNPs@gt[,parent1]))
parent2Coverage <- rep(-1, length(crossSNPs@gt[,parent2]))
parent1GQ <- rep(-1, length(crossSNPs@gt[,parent1]))
parent2GQ <- rep(-1, length(crossSNPs@gt[,parent2]))

### Define progeny for each cross and parents
parents1 <- c(parent1,parent2)
cross1 <- colnames(crossSNPs@gt)[which(!colnames(crossSNPs@gt) %in% c(parents1,"FORMAT"))]
filename1 <- "NF54xNHP4026_SNPs_6-20-19.csv"
crossName1 <- "NF54xNHP4026"

### loop through the entire vector and insert the 0/0 or 0/1 into vectors
for(i in 1:length(crossSNPs@gt[,parent1])){
  parent1GT[i] <- unlist(strsplit(crossSNPs@gt[i,parent1],":"))[1]
  parent2GT[i] <- unlist(strsplit(crossSNPs@gt[i,parent2],":"))[1]
  parent1Coverage[i] <- unlist(strsplit(crossSNPs@gt[i,parent1],":"))[3]
  parent2Coverage[i] <- unlist(strsplit(crossSNPs@gt[i,parent2],":"))[3]
  parent1GQ[i] <- unlist(strsplit(crossSNPs@gt[i,parent1],":"))[4]
  parent2GQ[i] <- unlist(strsplit(crossSNPs@gt[i,parent2],":"))[4]
}


hist(as.numeric(parent1GQ))
hist(as.numeric(parent2GQ))
hist(as.numeric(parent1Coverage),breaks=200)
hist(as.numeric(parent2Coverage),breaks=200)

plot(as.numeric(parent2Coverage),as.numeric(parent2GQ))

### SNPs where Parents 1 and 2 are homozygotes and parents are bi-allelic
biallelicParentalSNPs<- which((parent1GT=="0/0"&parent2GT=="1/1")|(parent1GT=="1/1"&parent2GT=="0/0"))
### Filter Parent SNPs to exclude anything with low coverage (<10) and low quality (GQ < 99)
coverageParentalSNPs <- which(as.numeric(parent1Coverage)>=10&as.numeric(parent2Coverage)>=10)
GQParentalSNPs <- which(as.numeric(parent1GQ)>=99&as.numeric(parent2GQ)>=99)

### Keep SNPs that are biparental and decent coverage and high quality
parentalSNPs <- intersect(intersect(biallelicParentalSNPs,coverageParentalSNPs),GQParentalSNPs)

### Store parental genotypes at SNPs where parents meet this above condition
SNPMatrix <- matrix(-1,ncol=(length(colnames(crossSNPs@gt))+6),nrow=length(crossSNPs@gt[,parent1]))
for(i in 1:length(crossSNPs@gt[,parent1])){
  SNPMatrix[i,1:7] <- crossSNPs@fix[i,1:7]
  for(j in 2:length(colnames(crossSNPs@gt))){
    SNPMatrix[i,j+6] <- unlist(strsplit(crossSNPs@gt[i,j],"[:]"))[1]
  }
}

colnames(SNPMatrix) <- c(colnames(crossSNPs@fix[,1:7]),colnames(crossSNPs@gt[,-1]))


### Store coverage at SNPs where parents meet this above condition
SNPCoverageMatrix <- data.frame(crossSNPs@fix[,1:7])
SNPCoverage <- rep(-1,length(crossSNPs@gt[,parent1]))

### For each sample which is a column in crossSNPs@gt (except the first column which gives format info) fill in a vector SNPCoverage and append to SNPCoverageMatrix
for(j in 2:length(colnames(crossSNPs@gt))){
  ### Fill in SNPCoverage Matrix for each row i in column j
  for(i in 1:length(crossSNPs@gt[,parent1])){
    SNPCoverage[i] <- as.numeric(unlist(strsplit(crossSNPs@gt[i,j],"[:]"))[3])
  }
  SNPCoverageMatrix <- cbind(SNPCoverageMatrix,SNPCoverage)
}

colnames(SNPCoverageMatrix) <- c(colnames(crossSNPs@fix[,1:7]),colnames(crossSNPs@gt[,-1]))
#SNPCoverageMatrix[,1] <- gsub("_v3","",gsub("Pf3D7_","",SNPCoverageMatrix[,1]))

### Store GQ score at SNPs where parents meet this above condition
SNPGQMatrix <- data.frame(crossSNPs@fix[,1:7])
SNPGQ <- rep(-1,length(crossSNPs@gt[,parent1]))

### For each sample which is a column in crossSNPs@gt (except the first column which gives format info) fill in a vector SNPCoverage and append to SNPCoverageMatrix
for(j in 2:length(colnames(crossSNPs@gt))){
  ### Fill in SNPCoverage Matrix for each row i in column j
  for(i in 1:length(crossSNPs@gt[,parent1])){
    SNPGQ[i] <- as.numeric(unlist(strsplit(crossSNPs@gt[i,j],"[:]"))[4])
  }
  SNPGQMatrix <- cbind(SNPGQMatrix,SNPGQ)
}

colnames(SNPGQMatrix) <- c(colnames(crossSNPs@fix[,1:7]),colnames(crossSNPs@gt[,-1]))
#SNPGQMatrix[,1] <- gsub("_v3","",gsub("Pf3D7_","",SNPGQMatrix[,1]))


hist(as.numeric(as.character(SNPGQMatrix[,parent1])),xlab="GQ Score",main=paste0(parent1," GQ Score Distribution"),breaks=50)
hist(as.numeric(as.character(SNPGQMatrix[,parent2])),xlab="GQ Score",main=paste0(parent1," GQ Score Distribution"),breaks=50)

### subset SNPs based on homozygosity of parents
parentalHomozySNPs <- SNPMatrix[parentalSNPs,]
parentalHomozySNPsGQ <- SNPGQMatrix[parentalSNPs,]
parentalHomozySNPsCoverage <- SNPCoverageMatrix[parentalSNPs,]

### Calculate mean coverage and GQ score for each SNP
MeanSNPGQ <- rowMeans(parentalHomozySNPsGQ[,8:dim(parentalHomozySNPsGQ)[2]],na.rm=TRUE)
MeanSNPCoverage <- rowMeans(parentalHomozySNPsCoverage[,8:dim(parentalHomozySNPsGQ)[2]],na.rm=TRUE)

hist(MeanSNPGQ,breaks=50)
hist(MeanSNPCoverage,breaks=20)

### Filter SNPs based on mean coverage and GQ across samples
goodSNPCoverage <- which(MeanSNPCoverage>3)
goodSNPGQ <- which(MeanSNPGQ>20)

### subset SNPs based on mean coverage and GQ score
filteredSNPs <- parentalHomozySNPs[intersect(goodSNPCoverage,goodSNPGQ),]
filteredSNPsGQ <- parentalHomozySNPsGQ[intersect(goodSNPCoverage,goodSNPGQ),]
filteredSNPsCoverage <- parentalHomozySNPsCoverage[intersect(goodSNPCoverage,goodSNPGQ),]

### filter out bad quality SNPs at an individual level
filteredSNPs[which(filteredSNPsCoverage<3|filteredSNPsGQ<10)] <- "./."


### Calculate mean coverage and GQ score for each sample
MeanGQ <- colMeans(filteredSNPsGQ[,8:dim(parentalHomozySNPsGQ)[2]],na.rm=TRUE)
MeanCoverage <- colMeans(filteredSNPsCoverage[,8:dim(parentalHomozySNPsGQ)[2]],na.rm=TRUE)

### Display sample GQ and Coverage distributions and set cut-offs
hist(MeanGQ,breaks=20)
hist(MeanCoverage,breaks=50)

### Filter samples based on mean coverage and GQ across samples
goodCoverage <- which(MeanCoverage>2)
goodGQ <- which(MeanGQ>10)

### Sort SNPs based on location
filteredSNPs <- filteredSNPs[order(filteredSNPs[,"CHROM"],as.numeric(filteredSNPs[,"POS"])),]

### Count heterozygous calls in matrix
hetCounts <- matrix(0,nrow=length(colnames(crossSNPs@gt[,-1])),ncol=1)
hetProp <- matrix(0,nrow=length(colnames(crossSNPs@gt[,-1])),ncol=1)
missingCounts <- matrix(0,nrow=length(colnames(crossSNPs@gt[,-1])),ncol=1)

for(j in 2:length(colnames(crossSNPs@gt))){
  heterozygousCount <- length(which(filteredSNPs[,j+6]=="0/1"|filteredSNPs[,j+6]=="1/0"))
  hetCounts[j-1,1] <- heterozygousCount
  hetProp[j-1,1] <- heterozygousCount/length(which(filteredSNPs[,j+6]=="0/1"|filteredSNPs[,j+6]=="1/0"|filteredSNPs[,j+6]=="1/1"|filteredSNPs[,j+6]=="0/0"))
  missingCounts[j-1,1] <- length(which(filteredSNPs[,j+6]=="./."))
}

rownames(hetCounts) <- colnames(crossSNPs@gt[,-1])
rownames(hetProp) <- colnames(crossSNPs@gt[,-1])
rownames(missingCounts) <- colnames(crossSNPs@gt[,-1])

### histogram of heterozygous call counts
hist(hetCounts,breaks=100,xlab="Number of heterozygous SNPs per progeny",ylab="Counts")
hist(hetProp,breaks=50)
hist(missingCounts/length(filteredSNPs[,1]),breaks=50,xlab="Proportion of missing data per progeny")
rownames(hetProp)[which(hetProp<=0.7)]
missingKeep <- which(missingCounts/length(filteredSNPs[,1])<=0.8)
missingDrop <- which(missingCounts/length(filteredSNPs[,1])>0.8)

### Take intersection of samples with little missing data, good coverage and good GQ scores
samplesKeep <- intersect(intersect(goodCoverage,goodGQ),missingKeep)

### lambda for Poisson dist
meanErrorRate <- mean(hetProp,na.rm=TRUE)

### subset individual progeny with less than 50% missing data
filteredSNPsNM <- filteredSNPs[,c(1:7,7+samplesKeep)]

hetCountsNM <- matrix(0,nrow=length(samplesKeep),ncol=1)
hetPropNM <- matrix(0,nrow=length(samplesKeep),ncol=1)
missingCountsNM <- matrix(0,nrow=length(samplesKeep),ncol=1)

for(j in 2:(length(samplesKeep)+1)){
  heterozygousCount <- length(which(filteredSNPsNM[,j+6]=="0/1"|filteredSNPsNM[,j+6]=="1/0"))
  hetCountsNM[j-1,1] <- heterozygousCount
  hetPropNM[j-1,1] <- heterozygousCount/length(which(filteredSNPsNM[,j+6]=="0/1"|filteredSNPsNM[,j+6]=="1/0"|filteredSNPsNM[,j+6]=="1/1"|filteredSNPsNM[,j+6]=="0/0"))
  missingCountsNM[j-1,1] <- length(which(filteredSNPsNM[,j+6]=="./."))
}

rownames(hetCountsNM) <- colnames(crossSNPs@gt[,-1])[samplesKeep]
rownames(hetPropNM) <- colnames(crossSNPs@gt[,-1])[samplesKeep]
rownames(missingCountsNM) <- colnames(crossSNPs@gt[,-1])[samplesKeep]

hist(hetCountsNM,breaks=200,xlab="Number of heterozygous SNPs per progeny",ylab="Counts")
hist(hetPropNM,breaks=50)
hist(missingCountsNM/length(filteredSNPsNM[,1]),breaks=50,xlab="Proportion of missing data per progeny")

### Check for SNPs with high heterozygosity across samples
heterozygousSNPCount <- cbind(filteredSNPsNM[,1:7],rep(0,length(filteredSNPsNM[,1])))

for(i in 1:length(filteredSNPsNM[,1])){
  heterozygousSNPCount[i,8] <- length(which(filteredSNPsNM[i,8:length(filteredSNPsNM[1,])]=="0/1"|filteredSNPsNM[i,8:length(filteredSNPsNM[1,])]=="1/0"))
}

### Sliding window size of average cross over to look for above average number of mutations
### From histogram, take 0.02 as our average heterozygous error rate for this which gives us an expectation of 50 SNPs between errors
windowSize = 40
meanErrorRate <- mean(hetPropNM)
slidingWindowHetz <- matrix(0,nrow=dim(filteredSNPsNM)[1]/windowSize,ncol=length(colnames(filteredSNPsNM[,c(-1,-2,-3,-4,-5,-6,-7)])))
slidingWindowProb <- matrix(0,nrow=dim(filteredSNPsNM)[1]/windowSize,ncol=length(colnames(filteredSNPsNM[,c(-1,-2,-3,-4,-5,-6,-7)])))

for(i in 1:dim(slidingWindowHetz)[1]){
  for(j in 1:(length(colnames(filteredSNPsNM))-7)){
    slidingWindowHetz[i,j] <- length(which(filteredSNPsNM[(((i-1)*windowSize)+1):(i*windowSize),j+7]=="0/1"|filteredSNPsNM[(((i-1)*windowSize)+1):(i*windowSize),j+7]=="1/0"))
    slidingWindowProb[i,j] <- dbinom(slidingWindowHetz[i,j],windowSize,meanErrorRate)
  }
}

slidingWindowHetzOverlap <- matrix(0,nrow=dim(filteredSNPs)[1]-windowSize,ncol=length(colnames(filteredSNPsNM[,c(-1,-2,-3,-4,-5,-6,-7)])))

for(i in 1:(dim(filteredSNPs)[1]-windowSize)){
  for(j in 1:(length(colnames(filteredSNPsNM))-7)){
    slidingWindowHetzOverlap[i,j] <- length(which(filteredSNPs[(i):(i+windowSize),j+6]=="0/1"|filteredSNPs[(i):(i+windowSize),j+6]=="1/0"))
  }
}

colnames(slidingWindowHetz) <- colnames(filteredSNPsNM[,c(-1,-2,-3,-4,-5,-6,-7)])
colnames(slidingWindowProb) <- colnames(filteredSNPsNM[,c(-1,-2,-3,-4,-5,-6,-7)])
colnames(slidingWindowHetzOverlap) <- colnames(filteredSNPsNM[,c(-1,-2,-3,-4,-5,-6,-7)])

hist(slidingWindowHetz,breaks=100,main=bquote("Number of SNPs in" ~ .(windowSize) ~ "SNP sliding window"))
boxplot(slidingWindowHetz,ylab=bquote("Number of SNPs in" ~ .(windowSize) ~ "SNP sliding window"),xlab="Samples")
abline(h=10)

hist(slidingWindowProb,breaks=30,main=bquote("Number of SNPs in" ~ .(windowSize) ~ "SNP sliding window"))
boxplot(slidingWindowProb,ylab=bquote("Prob " ~ .(windowSize) ~ "SNP sliding window"),xlab="Samples")


### Count how many windows per progeny have greater than 10 mixed SNPs

# Choose cutoff for determining windows with significant number of mixed SNPs based on Poisson dist p-values corrected for number of tests
numTests <- (length(parent2GT)/windowSize)*(length(colnames(crossSNPs@gt))-1)
poissonPval <- cbind(seq(1,windowSize),rep(0,windowSize),rep(0,windowSize))
for(i in 1:windowSize){
  poissonPval[i,2] <- dpois(poissonPval[i,1],windowSize*meanErrorRate)
  poissonPval[i,3] <- abs(poissonPval[i,2]-(0.05/numTests))
}
pvalCutoff <- poissonPval[which.min(poissonPval[,3])]

highHetWindowCount <- matrix(0,ncol=1,nrow=(length(slidingWindowHetz[1,])))
for(j in 1:length(colnames(slidingWindowHetz))){
  highHetWindowCount[j,] <- length(which(slidingWindowHetz[,j]> pvalCutoff))
}
rownames(highHetWindowCount)  <- colnames(slidingWindowHetz)

### Check which columns have the maximum heterozygous counts for a window above expected
maxHet <- matrix(0,ncol=1,nrow=length(colnames(slidingWindowHetz)))
for(i in 1:length(colnames(slidingWindowHetz))){
  maxHet[i,1] <- max(slidingWindowHetz[,i])
}
rownames(maxHet) <- colnames(slidingWindowHetz)
colnames(slidingWindowHetz)[which(maxHet>15)]

### Create table of possibly mixed samples

possibleMixed <- merge(hetCountsNM,highHetWindowCount,by="row.names",all=TRUE)
write.csv(possibleMixed,"mixed_info_6-20-19.csv")

### HeatMap of genome scan

library(RColorBrewer)
library(gplots)
library(vegan)
library(grDevices)

## Make vector of colors for values below threshold
rc1 <- colorpanel(8,"yellow", "yellow")
## Make vector of colors for values above threshold
rc2 <- colorpanel(8,"lightblue", "darkblue")
rampcols <- c(rc1, rc2)
## In your example, this line sets the color for values between 49 and 51. 
#rampcols[c(-nHalf, nHalf+1)] <- rgb(t(col2rgb("green")), maxColorValue=256) 
rb1 <- seq(0, (pvalCutoff), length.out=8)
rb2 <- seq((pvalCutoff), windowSize, length.out=8)[-1]
rampbreaks <- c(-4,rb1,rb2,windowSize+10)


mat <- slidingWindowHetz[,colnames(slidingWindowHetz) %in% cross1]

png(paste0(crossName1,"_hetwindow_6-20-19.png"), height=20, width=24, units="in", res=220)
par(oma=c(0,0,0,0))

heatmap.2(mat, Rowv=FALSE, Colv=FALSE, dendrogram="none", 
           trace="none", density="none", scale="none",col = rampcols, 
           breaks = rampbreaks, cexRow=1,cexCol=1.5,na.rm=TRUE,na.color="grey90",
           key=FALSE,margins=c(25, 15), srtCol=90, labCol=NULL, labRow=NULL)
 
legend(0,1,c("background error","above background heterzygousity","high heterozygosity"),col=c("yellow","lightblue","darkblue"),pch=c(15,15,15),bty="n",y.intersp=1,cex=1.5)

dev.off()

### exclude mixed progeny

clonal <- possibleMixed[which(possibleMixed[,3]==0),1]

clonalNMHomozySNPs <- cbind(filteredSNPsNM[,c(1:7)],filteredSNPsNM[,clonal])

### Filter alleles to remove alleles with missing data and heterozygous SNP calls
# Find alleles with no heterozygous SNP calls)
hetAlleles <- rep(0,length(clonalNMHomozySNPs[,1]))
for (i in 1:length(clonalNMHomozySNPs[,1])){
  hetAlleles[i] <- length(which(clonalNMHomozySNPs[i,]=="0/1"|clonalNMHomozySNPs[i,]=="1/0"))
}

homozygousSNPAlleles <- which(hetAlleles<round(dim(clonalNMHomozySNPs)[2]*.1))

# Find alleles with no missing data
missingAlleles <- rep(0,length(clonalNMHomozySNPs[,1]))
for (i in 1:length(clonalNMHomozySNPs[,1])){
  missingAlleles[i] <- length(which(clonalNMHomozySNPs[i,]=="./."))
}

nomissingSNPAlleles <- which(missingAlleles<40)

#Exclude SNPs with missing data and heterzygous calls
filteredClonalNMHomozySNPs <- clonalNMHomozySNPs[intersect(nomissingSNPAlleles,homozygousSNPAlleles),]

filteredClonalNMSNPs <- filteredClonalNMHomozySNPs

### Phase SNPs with MKK2835 as 0 and NHP1337 as 1
parent1FilteredSNPs <- filteredClonalNMHomozySNPs[,parent1]
parent2FilteredSNPs <- filteredClonalNMHomozySNPs[,parent2]
for(i in colnames(filteredClonalNMHomozySNPs)[8:length(filteredClonalNMHomozySNPs[1,])]){
  for(j in 1:length(filteredClonalNMSNPs[,parent1])){
    if(filteredClonalNMHomozySNPs[j,i]==parent1FilteredSNPs[j]){
      filteredClonalNMSNPs[j,i]<- 0
    }
    if(filteredClonalNMHomozySNPs[j,i]==parent2FilteredSNPs[j]){
      filteredClonalNMSNPs[j,i]<- 1
    }
  }
}
filteredClonalNMSNPs[which(filteredClonalNMSNPs=="./.")] <- "NA"
filteredClonalNMSNPs[which(filteredClonalNMSNPs=="0/.")] <- "NA"
filteredClonalNMSNPs[which(filteredClonalNMSNPs=="./0")] <- "NA"
filteredClonalNMSNPs[which(filteredClonalNMSNPs=="1/.")] <- "NA"
filteredClonalNMSNPs[which(filteredClonalNMSNPs=="./1")] <- "NA"
filteredClonalNMSNPs[which(filteredClonalNMSNPs=="0/1")] <- "NA"
filteredClonalNMSNPs[which(filteredClonalNMSNPs=="1/0")] <- "NA"
filteredClonalNMSNPs[which(filteredClonalNMSNPs=="NA")] <- "NA"

### Create matrix for heatmap
filteredClonalNMSNPs_map <- matrix(0,ncol=length(filteredClonalNMHomozySNPs[1,]),nrow=length(filteredClonalNMHomozySNPs[,1]))
filteredClonalNMSNPs_map[which(filteredClonalNMSNPs=="1")] <- 1
filteredClonalNMSNPs_map[which(is.na(filteredClonalNMSNPs)==TRUE)] <- -1

### Inheritance Heat Map 
## Make vector of colors for values below threshold
rc1 <- colorpanel(8,"white", "red")    
## Make vector of colors for values above threshold
rc2 <- colorpanel(8,"black", hsv(0.175,1,0.95))
rampcols <- c(rc1, rc2)
## In your example, this line sets the color for values between 49 and 51. 
#rampcols[c(-nHalf, nHalf+1)] <- rgb(t(col2rgb("green")), maxColorValue=256) 
rb1 <- seq(-1, 0, length.out=8)
rb2 <- seq(1, 3, length.out=8)[-1]
rampbreaks <- c(-4,rb1,rb2,4)

mat <- t(filteredClonalNMSNPs_map[,-c(1,2,3,4,5,6,7)])

png("SNPmap_MKK2835xNHP1337_6-20-19.png", height=20, width=24, units="in", res=220)
par(oma=c(0,0,0,0))

heatmap.2(mat, Rowv=FALSE, Colv=FALSE, dendrogram="none", 
           trace="none", density="none", scale="none",col = rampcols, 
           breaks = rampbreaks, cexRow=1,cexCol=2,colRow="black",na.rm=TRUE,na.color="grey90",
           key=FALSE,margins=c(5, 15), srtCol=360,colCol="black")
 
legend(0,1,c("MKK2835","NHP1337"),col=c("black","red"),pch=c(15,15),bty="n",y.intersp=1,cex=1.5)

dev.off()

### Write cross files
cross1SNPs <- cbind(paste(gsub("_v3","",gsub("Pf3D7_","",filteredClonalNMSNPs[,1]))), filteredClonalNMSNPs[,2],filteredClonalNMSNPs[,colnames(filteredClonalNMSNPs) %in% parents1],filteredClonalNMSNPs[,colnames(filteredClonalNMSNPs) %in% cross1])
cross1SNPs_sorted <- cross1SNPs[order(cross1SNPs[,1],as.numeric(cross1SNPs[,2]),decreasing=FALSE),]
rownames(cross1SNPs_sorted) <- paste("M",seq(1:length(cross1SNPs_sorted[,1])),sep="",collapse=NULL)
write.csv(t(cross1SNPs_sorted),file=filename1)


### Characterize Unique Recombinants

# Read in cross_data
#unique_recombinants <- function(filename, crossName, cross, parents, growthCross = NULL){
filename <- filename1
crossName <- crossName1
cross <- cross1
parentsAll <- parents1
growthCross <- NULL
OldMapFile <- "NF54xNHP4026_progenyOldMap.csv"
 SNPdata_full <- read.cross("csv", file=filename, genotypes=c(0,1), alleles=c(0,1), na.strings="NA", estimate.map=FALSE)

  SNPDataCSV <- read.csv(file=filename,header=FALSE, as.is=TRUE, na.strings="NA")

  ### Take into account previous cross and subset by those included
  progenyInOldMap <- read.csv(file=OldMapFile, header=FALSE, as.is = TRUE)
  
  oldCrossExcluded <- NULL
  
  clonenames_full <- SNPDataCSV[c(-1,-2,-3),1]
  
  SNPdata <- subset(SNPdata_full, ind=(!clonenames_full %in% oldCrossExcluded))
  clonenamesAll <- clonenames_full[which(!clonenames_full %in% oldCrossExcluded)]
  
  ### If no growth data use all progeny
  SNPdata <- SNPdata_full
  clonenamesAll <- SNPDataCSV[c(-1,-2,-3),1]
  
  ### Compare Genotypes
  cg <- comparegeno(SNPdata)
  hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="Proportion matching genotypes", main="Histogram of relatedness among genotypes")

  ### Identify clones with lots of matching alleles
  wh <- rbind(which(cg > 0.90, arr=TRUE),which(is.na(cg) == TRUE, arr=TRUE))
  
  ### Visualize progeny clusters in igraph
  progenyGraph <- graph_from_edgelist(wh,directed=FALSE)
  progenyGraphSim <- simplify(progenyGraph, remove.multiple = TRUE, remove.loops = TRUE)
  graphClusters <- cbind(seq(1,length(clonenamesAll)),unlist(components(progenyGraph)$membership))
  graphClusterNumbers <- which(components(progenyGraph)$csize > 1)
  
  progenyClustersList <- list()
  progenyDrop <- c()
  
  set.seed(1389191)
  
  #Default to random if no growth data is entered, except for clusters that contain parents
  if(is.null(growthCross)){
    i <- 1
    while (i <= length(graphClusterNumbers)){
      progenyClustersList <- c(progenyClustersList,list(which(graphClusters[,2]==graphClusterNumbers[i])))
      clusterMembers <- clonenamesAll[which(graphClusters[,2]==graphClusterNumbers[i])]
      if( any(clusterMembers %in% parentsAll)){
        progenyKeep <- intersect(clusterMembers,parentsAll) 
      } else if( any(clusterMembers %in% progenyInOldMap[,1])){
        progenyKeep <- intersect(clusterMembers,progenyInOldMap[,1])
      }  else{
        progenyKeep <- sample(clusterMembers,1)
      }
      progenyDrop <- c(progenyDrop,setdiff(clusterMembers,progenyKeep))
      i = i + 1
    }
  }
  if(!(is.null(growthCross))){
    # Find any progeny with missing growth data and set growth data to average
    if(any(is.na(growthCross[,"Growth"]))){
      missingGrowth <- which(is.na(growthCross[,"Growth"]))
      noMissingGrowth <- growthCross[,"Growth"][!is.na(growthCross[,"Growth"])]
      averageGrowth <- mean(noMissingGrowth)
      growthCross[missingGrowth,2] <- averageGrowth
    }
    # Loop through genotype clusters and keep cluster member with lowest growth rate or parents if parents are in cluster
    i=1
    while (i <= length(graphClusterNumbers)){
      progenyClustersList <- c(progenyClustersList,list(which(graphClusters[,2]==graphClusterNumbers[i])))
      clusterMembers <- clonenamesAll[which(graphClusters[,2]==graphClusterNumbers[i])]
      # If parents are in the cluster keep the parent, else sample 1 from cluster members with minimum growth values
      if( any(clusterMembers %in% parentsAll)){
        progenyKeep <- intersect(clusterMembers,parentsAll) 
      }  else if( any(clusterMembers %in% progenyInOldMap[,1])){
        progenyKeep <- intersect(clusterMembers,progenyInOldMap[,1])
      }  else{
        progenyKeep <- clusterMembers[which.min(growthCross[growthCross[,"Progeny_ID"] %in% clusterMembers,2])]
      }
      # Append individuals in the cluster that are not selected for keeping to the drop list
      progenyDrop <- c(progenyDrop,setdiff(clusterMembers,progenyKeep))
      i = i + 1
    }
  } 
  
  cloningInfo <- read.delim("cloningRounds.csv",sep=",",as.is=TRUE,header=FALSE)
  
  coords <- layout_with_fr(progenyGraphSim)
  rownames(coords) <- clonenamesAll
  layoutCoords <- cbind(coords,unlist(components(progenyGraph)$membership))
  TestCoords <- coords
  center <- colMeans(coords)
  radius <- sqrt((coords[which(layoutCoords[,3]==2),1]-center[1])^2+(coords[which(layoutCoords[,3]==2),2]-center[2])^2)
  theta <- atan((coords[which(layoutCoords[,3]==2),2]-center[2])/(coords[which(layoutCoords[,3]==2),1]-center[1]))
  TestCoords[which(layoutCoords[,3]==2),1] <- radius*5*cos(theta*pi)+center[1]
  TestCoords[which(layoutCoords[,3]==2),2] <- radius*5*sin(theta*pi)+center[2]
  
  for(j in graphClusterNumbers){
    NewCoords <- TestCoords[which(layoutCoords[,3]==j),1:2]
    if(length(NewCoords[,1])>=4){
    numPoints <- length(TestCoords[which(layoutCoords[,3]==j),1])
    center <- colMeans(NewCoords)
    numVer <- 1
    thetaVals <- 0
    radiusVals <- seq(0,7,by=.9)
    verCount <- 1
    NewCoords[verCount,1] <- radiusVals[1]*cos(0)+center[1]
    NewCoords[verCount,2] <- radiusVals[1]*sin(0)+center[2]
    verLayer <- c(7,14,21,28,35,38)
    verCount <- verCount + 1
    i=2
    while(verCount <= numPoints){
      if((verCount+verLayer[i-1])<=numPoints){
        thetaVals <- seq(0,(2*pi),by=(2*pi/(verLayer[i-1])))
        NewCoords[verCount:(verCount+verLayer[i-1]),1] <- radiusVals[i]*cos(thetaVals)+center[1]
        NewCoords[verCount:(verCount+verLayer[i-1]),2] <- radiusVals[i]*sin(thetaVals)+center[2]
        verCount <- verCount+verLayer[i-1]
        i=i+1
      } else if((verCount+verLayer[i-1])>numPoints){
        thetaVals <- seq(0,(2*pi),by=(2*pi/(numPoints-(verCount-1))))
        NewCoords[verCount:(verCount+length(thetaVals)-2),1] <- radiusVals[i]*cos(thetaVals[-1])+center[1]
        NewCoords[verCount:(verCount+length(thetaVals)-2),2] <- radiusVals[i]*sin(thetaVals[-1])+center[2]
        verCount <- verCount+verLayer[i-1]
        i=i+1
      }
    }
    
    TestCoords[which(layoutCoords[,3]==j),1] <- NewCoords[,1]
    TestCoords[which(layoutCoords[,3]==j),2] <- NewCoords[,2]
    }
  }
  
  plot(progenyGraphSim,layout=TestCoords)
  
  png(paste0(crossName,"_genotype_clusters_7-21-20_colorcoded_newlegend.png"), height=24, width=40, units="in", res=220)
  #plot.igraph(progenyGraph,vertex.size=7,vertex.label.cex=2,mark.groups=list(c(7,11),c(8,13),c(23,26,41)))
  plot.igraph(progenyGraphSim,vertex.label=rep("",length(clonenamesAll)),edge.color="black",vertex.color=cloningInfo[,3],vertex.size=7,vertex.label.cex=2.5,mark.groups=progenyClustersList,mark.shape=1,mark.border=c("darkorange","#66CC66","#66CC66","#66CC66","#66CC66","#66CC66","black","#66CC66","#66CC66","#66CC66","#66CC66","#66CC66","#66CC66","#66CC66","#66CC66","#66CC66","#66CC66","#66CC66"),mark.col=c("darkorange","#66CC66","#66CC66","#66CC66","#66CC66","#66CC66","gray20","#66CC66","#66CC66","#66CC66","#66CC66","#66CC66","#66CC66","#66CC66","#66CC66","#66CC66","#66CC66","#66CC66"))
  legend("topright",c("Parent","NF54GFPLucXNHP4026 - Round 1 (0 days)","NF54GFPLucXNHP4026 - Round 2 (14 days)","NF54GFPLucXNHP4026 - Round 3 (19 days)","NF54WTxNHP4026 (5 days)"),col=c("gold","#1A237E","#039BE5","#B3E5FC","red"),pch=19,cex=3)
  dev.off()
  
  ### Drop non-unique progeny from each cluster
  SNPDataSubset <- subset(SNPdata, ind=(!clonenamesAll %in% progenyDrop))
  
  ### Find markers with exact same pattern to thin out duplicate markers
  markerClusters <- findDupMarkers(SNPDataSubset, exact.only=FALSE)
  
  ### Look for skewed alleles
  gt <- geno.table(SNPDataSubset)
  gt[gt$P.value < 0.05/totmar(SNPDataSubset),]
  
  beginSNP <- vector()
  beginLocation <- vector()
  endSNP <- vector()
  endLocation <- vector()
  middleSNP <- vector()
  middleLocation <- vector()
  snpDrop <- vector()
  snpKeep <- vector()
  
  ## Write csv of matching clusters
  filenameClusters <- paste("clusers_", crossName, ".csv", sep = "")
  for(i in 1:length(names(markerClusters))){
    beginSNP <- c(beginSNP, names(markerClusters)[i])
    beginIndex <- gsub("M","", names(markerClusters)[i])
    beginIndex <- as.numeric(beginIndex)
    beginLocation <- c(beginLocation, SNPDataCSV[[beginIndex+1]][2] )
    if(length(markerClusters[[i]]) ==  1){
      middleSNP <- c(middleSNP, names(markerClusters)[i])
      middleLocation <- c(middleLocation, SNPDataCSV[[beginIndex+1]][2] )
    }
    else{
      middleSNP <- c(middleSNP, markerClusters[[names(markerClusters[i])]][length(markerClusters[[i]])/2])
      middleIndex <- gsub("M","", markerClusters[[names(markerClusters[i])]][length(markerClusters[[i]])/2])
      middleIndex <- as.numeric(middleIndex)
      middleLocation <- c(middleLocation, SNPDataCSV[[middleIndex+1]][2])
    }
    endSNP <- c(endSNP, markerClusters[[names(markerClusters[i])]][length(markerClusters[[i]])] )
    endIndex <- gsub("M", "", markerClusters[[names(markerClusters[i])]][length(markerClusters[[i]])])
    endIndex <- as.numeric(endIndex)
    endLocation <- c(endLocation, SNPDataCSV[[endIndex+1]][2])
    
    dropIt <- names(markerClusters[i])
    dropIt <- c(dropIt, markerClusters[[names(markerClusters[i])]][beginIndex:endIndex-1])
    dropIt <- c(names(markerClusters[i]), markerClusters[[names(markerClusters[i])]][-length(markerClusters[[i]])/2])
    snpDrop <- c(snpDrop, dropIt)
    snpKeep <- c(snpKeep, markerClusters[[names(markerClusters[i])]][length(markerClusters[[i]])/2])
  } # End of For loop
  
  allClusters <- cbind( beginSNP, beginLocation, middleSNP, middleLocation, endSNP, endLocation)
  #row.names(allClusters)<- cbind("Middle", "Location Middle", "Begin", "Location Begin", "End", "Location End")
  write.csv(allClusters, file=filenameClusters)

  ## Drop the SNPs to Drop
  SNPDataSubset <- drop.markers(SNPDataSubset, snpDrop)

  # Create CSV File for JoinMap
  filenameJoinMap <- paste("JoinMap_", crossName, ".csv", sep = "")
  editedParents <- gsub("1","a",filteredClonalNMSNPs[,colnames(filteredClonalNMSNPs) %in% parents])
  editedParents <- gsub("0","b",editedParents)
  
  
  # Progeny Names from SNPDataSubset
  progenyNames <- as.vector(pull.pheno(SNPDataSubset))
  
  
  ### Append between chromomsome breaks to SNPDataCSV
 # chrEndMarkers <- as.character(c("MC1","02",rep(3,length(progenyNames))))
  chrEndMarkers <- data.frame(stringsAsFactors = FALSE)[1:(length(progenyNames)+3), ]
  count <- 1
  chromosome <- c("01","02","03","04","05","06","07","08","09","10","11","12","13","14")
  for(i in 2:14){
    for(j in 1:4){
      chrEndMarkers <- cbind(chrEndMarkers,data.frame(c(paste("MC",count,sep=""),chromosome[i],as.character(j),rep("3",length(progenyNames))),stringsAsFactors = FALSE))
      names(chrEndMarkers)[count] <- paste("MC",count,sep="")
      count <- count + 1
    }
  }
  
  #test <- data.frame(c(paste("MC",count,sep=""),chromosome[i],as.character(j),rep("3",length(progenyNames))),stringsAsFactors = FALSE)
  
  ### Subset SNPDataCSV to include only clonal unique recombinant progeny
  SNPDataCSVSubset <- SNPDataCSV[c(1,2,3,which(SNPDataCSV[,1] %in% progenyNames)),]
  SNPDataCSVSubset <- SNPDataCSVSubset[,c(1,which(!SNPDataCSV[1,] %in% snpDrop))]
  
  
  SNPDataCSVCleanedSubset <- SNPDataCSV[c(1,2,3,which(SNPDataCSV[,1] %in% progenyNames)),]
  write.csv(SNPDataCSVCleanedSubset,paste("AllFinal_filtered_SNPS_",crossName1,"071219.csv"))
  
  write.csv(SNPDataCSVSubset,paste("MapFinal_filtered_SNPS_",crossName1,"071219.csv"))
  
  SNPDataCSVCleanedSubset <- read.cross("csv", file=paste("MapFinal_filtered_SNPS_",crossName1,"071219.csv"), genotypes=c(0,1), alleles=c(0,1), na.strings="NA", estimate.map=FALSE)

  
  ### For each progeny go across the genome and add a marker for every 5 kb
  genomeLength <- 23292622
  windowSize <- 5000
  numIntervals <- round((genomeLength/windowSize))
  SNPwindowMap <- matrix(-1000,ncol=(dim(SNPDataCSVSubset)[1]-3),nrow=(numIntervals))
  row.names(SNPwindowMap) <- seq(1,genomeLength,by=windowSize)
  colnames(SNPwindowMap) <- SNPDataCSVSubset[c(-1,-2,-3),1]
    
  chromosomeLength <- cbind(seq(1,15),c(0,640851,947102,1067971,1200490,1343557,1418242,1445207,1472805,1541735,1687656,2038340,2271494,2925236,3291936))
  SNPDistance <- rep(NA,length(SNPDataCSVSubset[1,]))
  for (i in 2:length(SNPDataCSVSubset[1,])){
    SNPDistance[i] <- (as.numeric(as.character(SNPDataCSVSubset[3,i]))+sum(chromosomeLength[1:which(chromosomeLength[,1]==as.numeric(as.character(SNPDataCSVSubset[2,i]))),2]))
  }
  
  # Read in core, non-core definition from 
  GenomeDef <- read.delim("PFal_CoreGenome_Milesetal.csv",sep=",",as.is=TRUE)
  DefineCore <- GenomeDef[which(GenomeDef[,"Name"]=="Core"),c("Chrom","Start","End")]
  DefineCore$StartCum <- rep(NA,length(DefineCore[,1]))
  DefineCore$EndCum <- rep(NA,length(DefineCore[,1]))
  for (i in 1:14){
    DefineCore[which(DefineCore[,"Chrom"]==i),"StartCum"] <- DefineCore[which(DefineCore[,"Chrom"]==i),"Start"] + sum(chromosomeLength[1:i,2])
    DefineCore[which(DefineCore[,"Chrom"]==i),"EndCum"] <- DefineCore[which(DefineCore[,"Chrom"]==i),"End"] + sum(chromosomeLength[1:i,2])
  }
  
  CoreWindows <- c()
  for(i in 1:length(DefineCore[,1])){
    CurrentCore <- rownames(SNPwindowMap)[intersect(which(as.numeric(rownames(SNPwindowMap)) <= DefineCore[i,"EndCum"]),which(as.numeric(rownames(SNPwindowMap)) >= DefineCore[i,"StartCum"]))]
    CoreWindows <- c(CoreWindows,CurrentCore)
  }
  
  NonCoreWindow <- rownames(SNPwindowMap)[which(!rownames(SNPwindowMap) %in% CoreWindows)]
  
  totalDistance <- rep(0,15)
  for(i in chromosomeLength[2:15,1]){
    totalDistance[i] <- sum(chromosomeLength[which(chromosomeLength[,1]<=chromosomeLength[i,1]),2])
  }
  
  # Create matrix to append and then sort
  ChrEndLocMarkers <- matrix(3,ncol=(dim(SNPDataCSVSubset)[1]-3),nrow=130)
  row.names(ChrEndLocMarkers) <- c(totalDistance[2:14],(totalDistance[2:14]+1),(totalDistance[2:14]+2),(totalDistance[2:14]+3),(totalDistance[2:14]+4),(totalDistance[2:14]+5),(totalDistance[2:14]+6),(totalDistance[2:14]+7),(totalDistance[2:14]+8),(totalDistance[2:14]+9))
  
  #loop through progeny
  for(i in 1:(dim(SNPDataCSVSubset)[1]-3)){
    # Inititalize chrCounter
    chrCounter <- 2
    # start with first maker as current marker
    currentMarker <- as.numeric(SNPDataCSVSubset[3+i,(which(is.na(SNPDataCSVSubset[3+i,-1])==FALSE)[1]+1)])
    #currentProgenySNPs <- rep(NA,dim(SNPwindowMap)[1])
    #loop through 5kb windows and pull dominant SNP
    for(j in 1:(genomeLength/windowSize)){
      currentSNPs <- SNPDataCSVSubset[3+i,intersect(which(SNPDistance>=as.numeric(row.names(SNPwindowMap)[j])),which(SNPDistance<as.numeric(row.names(SNPwindowMap)[j+1])))]
      #If no SNPs in window set to NA
      if(length(currentSNPs)==0){
        SNPwindowMap[j,i] <- currentMarker
      } else if(length(currentSNPs)>0){ #if SNPs in window make a table to count NA, 0 and 1
        SNPtable <- c(length(which(is.na(currentSNPs)==TRUE)),length(which(currentSNPs==0)),length(which(currentSNPs==1)))
        names(SNPtable) <- c(NA,0,1)
        #SNPtable <- sort(SNPtable,decreasing=TRUE)
        # If there are only NAs then set to NA
        if(SNPtable[1]==length(currentSNPs)){
          SNPwindowMap[j,i] <- currentMarker
        } else{
          SNPwindowMap[j,i] <- as.numeric(names(which.max(SNPtable[2:3])))
          currentMarker <- as.numeric(names(which.max(SNPtable[2:3])))
        }
      }
      # If at the end of the chromosome, increment chromosome and reset currentMarker to first non NA marker on chromosome
      if(as.numeric(rownames(SNPwindowMap)[j])>=totalDistance[chrCounter]){
        chrCounter <- chrCounter + 1
        SNPsInChr <- intersect(which(SNPDistance>=totalDistance[chrCounter-1]),which(SNPDistance<totalDistance[chrCounter]))
        currentMarker <- as.numeric(SNPDataCSVSubset[3+i,(which(is.na(SNPDataCSVSubset[3+i,SNPsInChr])==FALSE)[1])+SNPsInChr[1]-1])
      }
    }
  }
  
  # Need to set the non-core region to -1
  SNPwindowCore <- SNPwindowMap[,c(2,1,3:length(SNPwindowMap[1,]))]
  for (i in NonCoreWindow){
    SNPwindowCore[i,] <- rep(-1,length(SNPwindowCore[i,]))
  }
  
  
  # Bind together chromosome breaks and phased 5kb data
  SNPWindowMapChrEndMarkers <- rbind(SNPwindowCore,ChrEndLocMarkers)
  SNPDataChrEndMarkersSort <- SNPWindowMapChrEndMarkers[order(as.numeric(rownames(SNPWindowMapChrEndMarkers)),decreasing=FALSE),]
  
  
  ### Make labels for heatMap
  chrLabels <- rep(" ",length(rownames(SNPDataChrEndMarkersSort)))
  chrLengthCount <- 1
  totalDistance <- rep(0,15)
  for(i in chromosomeLength[2:15,1]){
    totalDistance[i] <- sum(chromosomeLength[which(chromosomeLength[,1]<=chromosomeLength[i,1]),2])
    chrMarkers <- intersect(which(as.numeric(rownames(SNPDataChrEndMarkersSort))<=totalDistance[i]),which(as.numeric(rownames(SNPDataChrEndMarkersSort))>=totalDistance[i-1]))
    chrLabels[chrLengthCount + (round(length(chrMarkers)/3)*2)] <- i-1
    chrLengthCount <- chrLengthCount + length(chrMarkers)
  }
  
  ### Sort by cloning round
  CrossProg <- colnames(SNPDataChrEndMarkersSort)
  CRSort <- cloningInfo[which(cloningInfo[,1]==CrossProg[1]),c(1,2)]
  for(i in 2:length(CrossProg)){
    CRSort <- rbind(CRSort,cloningInfo[which(cloningInfo[,1]==CrossProg[i]),c(1,2)])
  }
  CRSort$Type <- 2
  BreakPos <- rbind(c("",1,1),c("",2,1),c("",3,1),c("",4,1))
  colnames(BreakPos) <- c("V1","V2","Type")
  CRSortWBreaks <- rbind(CRSort,BreakPos)
  CRSortWBreaks[which(CRSortWBreaks[,1] %in% c("Sample_NF54_GFP_LUC","Sample_ART-1")),2] <- 0
  
  ### Add breaks for cloning round to SNPDataChrEndMarkersSort
  CRBreaks <- cbind(rep(-4,length(SNPDataChrEndMarkersSort[,1])),rep(-4,length(SNPDataChrEndMarkersSort[,1])),rep(-4,length(SNPDataChrEndMarkersSort[,1])),rep(-4,length(SNPDataChrEndMarkersSort[,1])))
  SNPDataChrEndMarkersSortCRBreaks <- cbind(SNPDataChrEndMarkersSort,CRBreaks)
  
  
  ### Inheritance Heat Map 
  ## Make vector of colors for values below threshold
  rc1 <- colorpanel(8,"grey60", "black")    
  ## Make vector of colors for values above threshold
  rc2 <- colorpanel(8,"red", "gold")
  rampcols <- c("white",rc1, rc2)
  ## In your example, this line sets the color for values between 49 and 51. 
  #rampcols[c(-nHalf, nHalf+1)] <- rgb(t(col2rgb("green")), maxColorValue=256) 
  rb1 <- seq(-1, 0, length.out=8)
  rb2 <- seq(1, 3, length.out=8)[-1]
  rampbreaks <- c(-5,-4,rb1,rb2,4)
  
  #mat <- joinMapData_map[c(seq(3,dim(joinMapData)[1]),2,1),]
  mat <- t(SNPDataChrEndMarkersSortCRBreaks)[order(CRSortWBreaks[,2],CRSortWBreaks[,3]),]
  
  png(paste0("SNPmap_",crossName1,"_filtered_coverage_mapmarkers_physloc_071219_NF54_sorted_breaks_labels.png"), height=20, width=24, units="in", res=220)
  par(oma=c(0,0,0,0))
  
  heatmap.2(mat, Rowv=FALSE, Colv=FALSE, dendrogram="none", 
            trace="none", density="none", scale="none", col = rampcols, 
            breaks = rampbreaks, cexRow=1, cexCol=3, colRow=NULL, na.rm=TRUE, na.color="gray80",
            key=FALSE, margins=c(6, 15), srtCol=360, colCol=NULL, labCol=chrLabels)
  
  legend(0,1,c(paste0(parent1," Genotype"),paste0(parent2," Genotype")),col=c("black","red"),pch=c(15,15),bty="n",y.intersp=1,cex=1.5)
  
  dev.off()
  
  ### Subset SNPWindowCore by BioRep
  SNPwindowCoreCross1 <- SNPwindowCore[,which(colnames(SNPwindowCore) %in% cross1)]
  SNPwindowCoreCross2 <- SNPwindowCore[,which(colnames(SNPwindowCore) %in% cross2)]
  
  C1CR1Prog <- cloningInfo[which(cloningInfo[,2] == 1),1]
  C1CR2Prog <- cloningInfo[which(cloningInfo[,2] == 2),1]
  C1CR3Prog <- cloningInfo[which(cloningInfo[,2] == 3),1]
  C2CR1Prog <- cloningInfo[which(cloningInfo[,2] == 4),1]
  
  ### Find marker frequency in SNPWindowCore
  markerFreq5K <- data.frame(as.numeric(rownames(SNPwindowCore)))
  markerFreq5K <- cbind(markerFreq5K,(rowSums(SNPwindowCore,na.rm=TRUE)/dim(SNPwindowCore)[2]))
  markerFreq5KCross1 <- data.frame(as.numeric(rownames(SNPwindowCore)))
  markerFreq5KCross1 <- cbind(markerFreq5KCross1,(rowSums(SNPwindowCoreCross1,na.rm=TRUE)/dim(SNPwindowCoreCross1)[2]))
  markerFreq5KCross2 <- data.frame(as.numeric(rownames(SNPwindowCore)))
  markerFreq5KCross2 <- cbind(markerFreq5KCross2,(rowSums(SNPwindowCoreCross2,na.rm=TRUE)/dim(SNPwindowCoreCross2)[2]))
  markerFreq5KC1CR1 <- data.frame(as.numeric(rownames(SNPwindowCore)))
  markerFreq5KC1CR1 <- cbind(markerFreq5KC1CR1,(rowSums(SNPwindowCore[,which(colnames(SNPwindowCore) %in% C1CR1Prog)],na.rm=TRUE)/dim(SNPwindowCore[,which(colnames(SNPwindowCore) %in% C1CR1Prog)])[2]))
  markerFreq5KC1CR2 <- data.frame(as.numeric(rownames(SNPwindowCore)))
  markerFreq5KC1CR2 <- cbind(markerFreq5KC1CR2,(rowSums(SNPwindowCore[,which(colnames(SNPwindowCore) %in% C1CR2Prog)],na.rm=TRUE)/dim(SNPwindowCore[,which(colnames(SNPwindowCore) %in% C1CR2Prog)])[2]))
  markerFreq5KC1CR3 <- data.frame(as.numeric(rownames(SNPwindowCore)))
  markerFreq5KC1CR3 <- cbind(markerFreq5KC1CR3,(rowSums(SNPwindowCore[,which(colnames(SNPwindowCore) %in% C1CR3Prog)],na.rm=TRUE)/dim(SNPwindowCore[,which(colnames(SNPwindowCore) %in% C1CR3Prog)])[2]))
  markerFreq5KC2CR1 <- data.frame(as.numeric(rownames(SNPwindowCore)))
  markerFreq5KC2CR1 <- cbind(markerFreq5KC2CR1,(rowSums(SNPwindowCore[,which(colnames(SNPwindowCore) %in% C2CR1Prog)],na.rm=TRUE)/dim(SNPwindowCore[,which(colnames(SNPwindowCore) %in% C2CR1Prog)])[2]))
  
  AlleleFrequency <- cbind(markerFreq5K,markerFreq5KCross1[,2],markerFreq5KCross2[,2])
  colnames(AlleleFrequency) <- c("GenomeLoc","AllProgeny","BioRep1","BioRep2")
  AlleleFrequency$Chr <- AlleleFrequency[,1]
  AlleleFrequency$Loc <- AlleleFrequency[,1]
  for(i in 2:length(totalDistance)){
    ChrMarkers <- intersect(which(as.numeric(AlleleFrequency[,"GenomeLoc"])<=totalDistance[i]),which(as.numeric(AlleleFrequency[,"GenomeLoc"])>=totalDistance[i-1]))
    AlleleFrequency[ChrMarkers,"Chr"] <- i-1
    AlleleFrequency[ChrMarkers,"Loc"] <- AlleleFrequency[ChrMarkers,"GenomeLoc"]-totalDistance[i-1]
  }
  write.csv(AlleleFrequency,"NF54xNHP4026_AlleleFrequency.csv")
  
  chrLabels <- rep(" ",length(rownames(SNPDataChrEndMarkersSort)))
  chrLengthCount <- 1
  totalDistance <- rep(0,15)
  for(i in chromosomeLength[2:15,1]){
    totalDistance[i] <- sum(chromosomeLength[which(chromosomeLength[,1]<=chromosomeLength[i,1]),2])
    chrMarkers <- intersect(which(as.numeric(markerFreq5K[,1])<=totalDistance[i]),which(as.numeric(markerFreq5K[,1])>=totalDistance[i-1]))
    chrLabels[chrLengthCount + (round(length(chrMarkers)/3)*2)] <- i-1
    chrLengthCount <- chrLengthCount + length(chrMarkers)
  }
  
  skewedLoci <- markerFreq5K[c(intersect(which(markerFreq5K[,2]<=0.36),which(markerFreq5K[,2]>0)),which(markerFreq5K[,2]>=0.64)),]
  skewedLociCross1 <- markerFreq5KCross1[c(intersect(which(markerFreq5KCross1[,2]<=0.33),which(markerFreq5KCross1[,2]>0)),which(markerFreq5KCross1[,2]>=0.67)),]
  skewedLociCross2 <- markerFreq5KCross2[c(intersect(which(markerFreq5KCross2[,2]<=0.25),which(markerFreq5KCross2[,2]>0)),which(markerFreq5KCross2[,2]>=0.75)),]
  skewedLociBoth <- intersect(skewedLociCross1[,1],skewedLociCross2[,1])
  write.csv(skewedLociBoth,"NF54xNHP4026 skewed loci in BioRep 1 and 2.csv")
  
  
  
  png(paste0("FrequencySkewMap_",crossName1,"_filtered_All_goodcoverage_100119_p0.001.png"), height=10, width=28, units="in", res=220)
  par(mar=c(7,8,2,0))
  plot(markerFreq5K[,1],markerFreq5K[,2],ylim=c(0,1.0),xlim=c(800000,22500000),ylab="",xlab="",xaxt="none",yaxt="none",pch=19,cex=1)
  abline(h=c(0.69,0.31),lty=1,lwd=3,col="black")
  axis(1, at=c(320425.5,1114402,2121938,3256169,4528192,5909092,7340816,8799822,10307092,11921788,13784786,15939703,18538068,21646654),tick=FALSE,seq(1,14),cex.axis=3.5,line=1)
  axis(2, at=c(0.0,0.2,0.4,0.6,0.8,1.0),tick=TRUE,c(0.0,0.2,0.4,0.6,0.8,1.0),cex.axis=3.5,line=1)
  abline(v=c(640851,1587953,2655924,3856414,5199971,6618213,8063420,9536225,11077960,12765616,14803956,17075450,20000686))
  abline(h=0.5,lty=1,col="grey70",lwd=3)
  title(xlab="Genome Location (Chr)",ylab="NHP4026 Allele Frequency", line=5, cex.lab=3.5)
  dev.off()  
    
  png(paste0("FrequencySkewMap_",crossName1,"_filtered_All_goodcoverage_Bioreps_07212020_p0.001.png"), height=10, width=28, units="in", res=220)
  par(mar=c(7,8,2,0))
  plot(markerFreq5K[,1],markerFreq5K[,2],ylim=c(0,1.0),xlim=c(800000,22500000),ylab="",xlab="",xaxt="none",yaxt="none",pch=19,cex=1)
  points(markerFreq5KCross1[,1],markerFreq5KCross1[,2],pch=19,cex=1,col="red")
  points(markerFreq5KCross2[,1],markerFreq5KCross2[,2],pch=19,cex=1,col="blue")
  abline(h=c(0.69,0.31),lty=1,col="black",lwd=3)
  abline(h=c(0.81,0.19),lty=1,col="blue",lwd=3)
  abline(h=c(0.71,0.29),lty=1,col="red",lwd=3)
  abline(h=0.5,lty=1,col="grey70",lwd=3)
  axis(1, at=c(320425.5,1114402,2121938,3256169,4528192,5909092,7340816,8799822,10307092,11921788,13784786,15939703,18538068,21646654),tick=FALSE,seq(1,14),cex.axis=3.5,line=1)
  axis(2, at=c(0.0,0.2,0.4,0.6,0.8,1.0),tick=TRUE,c(0.0,0.2,0.4,0.6,0.8,1.0),cex.axis=3.5,line=1)
  abline(v=c(640851,1587953,2655924,3856414,5199971,6618213,8063420,9536225,11077960,12765616,14803956,17075450,20000686))
  title(xlab="Genome Location (Chr)",ylab="NHP4026 Allele Frequency", line=5, cex.lab=3.5)
  legend("topleft",c("All Progeny (84)","NF54GFPLuc x NHP4026 (57)","NF54 WT x NHP4026 2 (27)"),pch=19,col=c("black","red","blue"),bty="o",box.col="white",cex=2)
  dev.off()
  
  png(paste0("FrequencySkewMap_",crossName1,"_filtered_All_goodcoverage_cloningrounds.png"), height=10, width=28, units="in", res=220)
  par(mar=c(7,8,2,0))
  plot(markerFreq5KC1CR1[,1],markerFreq5KC1CR1[,2],ylim=c(0,1.0),xlim=c(800000,22500000),ylab="",xlab="",xaxt="none",yaxt="none",pch=15,cex=1,col="#1A237E")
  points(markerFreq5KC1CR2[,1],markerFreq5KC1CR2[,2],pch=16,cex=1,col="#039BE5")
  points(markerFreq5KC1CR3[,1],markerFreq5KC1CR3[,2],pch=17,cex=1,col="#B3E5FC")
  points(markerFreq5KC2CR1[,1],markerFreq5KC2CR1[,2],pch=18,cex=1,col="red")
  abline(h=c(0.69,0.31),lty=1,col="black",lwd=3)
  #abline(h=c(0.81,0.19),lty=1,col="blue",lwd=3)
  #abline(h=c(0.71,0.29),lty=1,col="red",lwd=3)
  abline(h=0.5,lty=1,col="grey70",lwd=3)
  axis(1, at=c(320425.5,1114402,2121938,3256169,4528192,5909092,7340816,8799822,10307092,11921788,13784786,15939703,18538068,21646654),tick=FALSE,seq(1,14),cex.axis=3.5,line=1)
  axis(2, at=c(0.0,0.2,0.4,0.6,0.8,1.0),tick=TRUE,c(0.0,0.2,0.4,0.6,0.8,1.0),cex.axis=3.5,line=1)
  abline(v=c(640851,1587953,2655924,3856414,5199971,6618213,8063420,9536225,11077960,12765616,14803956,17075450,20000686))
  title(xlab="Genome Location (Chr)",ylab="NHP4026 Allele Frequency", line=5, cex.lab=3.5)
  legend("topleft",c("NF54GFPLuc x NHP4026 CR1 (20) - Day 0","NF54GFPLuc x NHP4026 CR2 (12) - Day 14","NF54GFPLuc x NHP4026 CR3 (25) - Day 19","NF54 x NHP4026 CR1 (28) - Day 5"),pch=c(15,16,17,18),col=c("#1A237E","#039BE5","#B3E5FC","red"),bty="o",box.col="white",cex=1.5)
  dev.off()
  
  ### Count recombinations at each locus
  # Loop i through progeny
  #for(i in 1:(length(joinMapDataNumeric[,1]))){
  # Loop j through markers
  #for(j in 1:(length(joinMapDataNumeric[1,])-1)){
  # For all transitions between 0 and 1, ignoring NAs
  #if((joinMapDataNumeric[i,j]==1 | joinMapDataNumeric[i,j]==0) & (joinMapDataNumeric[i,j+1]==1 | joinMapDataNumeric[i,j+1]==0)){
  #  if(joinMapDataNumeric[i,j] != joinMapDataNumeric[i,j+1]){
  #    recombHotSpot[j,5] = as.numeric(recombHotSpot[j,5]) + 1
  #  }
  #}
  #}
  #}
  
  ### Count recombinations at each locus
  # Loop i through progeny
  for(i in 1:(length(joinMapDataNumeric[,1]))){
    # Loop j through markers
    markerType <- joinMapDataNumeric[i,1]
    for(j in 2:length(joinMapDataNumeric[1,])){
      # For all transitions between 0 and 1, ignoring NAs
      if(is.na(joinMapDataNumeric[i,j])==FALSE){
        if(markerType != joinMapDataNumeric[i,j]){
          recombHotSpot[j,5] = as.numeric(recombHotSpot[j,5]) + 1
        }
        markerType <- joinMapDataNumeric[i,j]
      }
    }
  }
  
  #recomb_coverage <- merge(recombHotSpot,SNPCoverageMatrix[,c(1,2,4,5,6,7,8,275)],by=c("CHROM","POS"))
  #recomb_coverage <- recomb_coverage[order(as.numeric(as.character(recomb_coverage[,"CHROM"])),as.numeric(as.character(recomb_coverage[,"POS"])),decreasing=FALSE),]
  #joinMapData_coverage <- joinMapData_map[,which(as.numeric(as.character(recomb_coverage[,"MKK2835_PL1_2G_S11_L005"]))>3)]
  #rownames(joinMapData_coverage) <- rownames(joinMapData)
  #colnames(joinMapData_coverage) <- colnames(joinMapData)[which(as.numeric(as.character(recomb_coverage[,"MKK2835_PL1_2G_S11_L005"]))>3)]
  #write.csv(joinMapData_coverage,"MKK2835xNHP1337Markers_UniqueRecombinantProgeny_sorted.csv")
  
  
  ### Make labels for plots
  chrLabels <- rep(" ",length(recombHotSpot[1,]))
  chrLengthCount <- 1
  for(i in c("01","02","03","04","05","06","07","08","09","10","11","12","13","14")){
    chrMarkers <- which(recombHotSpot[,2]==i)
    chrLabels[chrLengthCount + (round(length(chrMarkers)/3)*2)] <- as.numeric(i)
    chrLengthCount <- chrLengthCount + length(chrMarkers)
  }
  
  
  png(paste0("FrequencySkewMap_",crossName1,"_filtered_All_goodcoverage_071219.png"), height=10, width=24, units="in", res=220)
  plot(recombHotSpot[,4],recombHotSpot[,6],ylim=c(0,1.0),ylab="Allele Frequency",xlab="Genome Location (Chr)",xaxt="none",pch=19,cex=1,cex.lab=2.5,cex.axis=2)
  axis(1, at=c(320425.5,1114402,2121938,3256169,4528192,5909092,7340816,8799822,10307092,11921788,13784786,15939703,18538068,21646654),tick=FALSE,seq(1,14),cex.axis=2)
  abline(v=c(640851,1587953,2655924,3856414,5199971,6618213,8063420,9536225,11077960,12765616,14803956,17075450,20000686,23292622))
  dev.off()
  
  
  
  ### Single SNP analysis

  ### bind together markers and between chromosome bands
  SNPDataChrEndMarkers <- cbind(SNPDataCSVSubset,chrEndMarkers)
  
  SNPDataChrEndMarkersSort <- SNPDataChrEndMarkers[,order(SNPDataChrEndMarkers[2,],as.numeric(SNPDataChrEndMarkers[3,]),decreasing=FALSE)]
  
  # Take a subset based on progeny name & remove excess columns
  #joinMapData <- subset(filteredClonalNMSNPs, colnames(filteredClonalNMSNPs) %in% progenyNames)
  joinMapData <- SNPDataCSVSubset[c(-1,-2,-3),c(-1)]
  #joinMapData <- joinMapData[,which(!colnames(joinMapData) %in% c("POS", "ID", "REF", "ALT", "QUAL", "FILTER"))]
  # Add Marker Names
  #joinMapData[,1] <- paste("M",seq(1:length(joinMapData[,1])),sep="",collapse=NULL)
  rownames(joinMapData) <- SNPDataCSVSubset[c(-1,-2,-3),1]
  colnames(joinMapData) <- SNPDataCSVSubset[1,c(-1)]
  
  ### Make labels for heatMap
  chrLabels <- rep(" ",length(SNPDataChrEndMarkersSort[1,]))
  chrLengthCount <- 1
  for(i in c("01","02","03","04","05","06","07","08","09","10","11","12","13","14")){
    chrMarkers <- which(SNPDataChrEndMarkersSort[2,]==i)
    chrLabels[chrLengthCount + (round(length(chrMarkers)/3)*2)] <- as.numeric(i)
    chrLengthCount <- chrLengthCount + length(chrMarkers)
  }
  
  #joinMapData <- joinMapData[which(rownames(joinMapData) %in% progenyNames),]
  #joinMapData <- joinMapData[,which(!colnames(joinMapData) %in% snpDrop)]

  ### Create map matrix
  joinMapData_map <- matrix(-1,ncol=length(joinMapData[1,]),nrow=length(joinMapData[,1]))
  joinMapData_map[which(joinMapData=="1")] <- 1
  joinMapData_map[which(is.na(joinMapData)==TRUE)] <- -1
  joinMapData_map[which(joinMapData=="0")] <- 0
  joinMapData_map[which(joinMapData=="3")] <- 3
  
  ### Inheritance Heat Map 
  ## Make vector of colors for values below threshold
  rc1 <- colorpanel(8,"grey70", "black")    
  ## Make vector of colors for values above threshold
  rc2 <- colorpanel(8,"red", "white")
  rampcols <- c(rc1, rc2)
  ## In your example, this line sets the color for values between 49 and 51. 
  #rampcols[c(-nHalf, nHalf+1)] <- rgb(t(col2rgb("green")), maxColorValue=256) 
  rb1 <- seq(-1, 0, length.out=8)
  rb2 <- seq(1, 3, length.out=8)[-1]
  rampbreaks <- c(-4,rb1,rb2,4)
  
  #mat <- joinMapData_map[c(seq(3,dim(joinMapData)[1]),2,1),]
  mat <- t(SNPwindowMap)
  
  png(paste0("SNPmap_",crossName1,"_filtered_coverage_mapmarkers_physloc.png"), height=20, width=24, units="in", res=220)
  par(oma=c(0,0,0,0))
  
  heatmap.2(mat, Rowv=FALSE, Colv=FALSE, dendrogram="none", 
            trace="none", density="none", scale="none", col = rampcols, 
            breaks = rampbreaks, cexRow=1, cexCol=3, colRow=NULL, na.rm=TRUE, na.color="gray80",
            key=FALSE, margins=c(5, 15), srtCol=360, colCol=NULL, labCol=chrLabels)
  
  legend(0,1,c(paste0(parent1," Genotype"),paste0(parent2," Genotype")),col=c("black","red"),pch=c(15,15),bty="n",y.intersp=1,cex=1.5)
  
  dev.off()
  
  
  
  
  ### Create numeric matrix without columns between chromosomes
  joinMapDataNumeric <- joinMapData_map
  colnames(joinMapDataNumeric) <- colnames(joinMapData)
  rownames(joinMapDataNumeric) <- rownames(joinMapData)
  
  ### Find Recombination Hotspots
  recombHotSpot <- data.frame(colnames(joinMapData)[-1])
  recombHotSpot <- cbind(recombHotSpot,cross1SNPs_sorted[which(rownames(cross1SNPs_sorted) %in% colnames(joinMapData)[-1]),1])
  recombHotSpot <- cbind(recombHotSpot,cross1SNPs_sorted[which(rownames(cross1SNPs_sorted) %in% colnames(joinMapData)[-1]),2])
  
  # Find location based on genome size of each marker, chromosome length from https://www.genome.jp/kegg-bin/show_organism?org=pfa
  chromosomeLength <- cbind(seq(1,15),c(0,640851,947102,1067971,1200490,1343557,1418242,1445207,1472805,1541735,1687656,2038340,2271494,2925236,3291936))
  markerDistance <- rep(NA,length(recombHotSpot[,1]))
  for (i in 1:length(recombHotSpot[,1])){
    markerDistance[i] <- (as.numeric(as.character(recombHotSpot[i,3]))+sum(chromosomeLength[1:which(chromosomeLength[,1]==as.numeric(as.character(recombHotSpot[i,2]))),2]))
  }
  recombHotSpot <- cbind(recombHotSpot,markerDistance)
  recombHotSpot <- cbind(recombHotSpot,as.numeric(rep(0,length(recombHotSpot[,1]))))
  recombHotSpot <- cbind(recombHotSpot,rowSums(joinMapData_map[,which(colnames(joinMapData) %in% growthCross1[,1])])/length(which(colnames(joinMapData) %in% cross1)))
  recombHotSpot <- cbind(recombHotSpot,rowSums(joinMapData_map[,which(colnames(joinMapData) %in% cross2)])/length(which(colnames(joinMapData) %in% cross2)))
  recombHotSpot <- cbind(recombHotSpot,(colSums(joinMapDataNumeric,na.rm=TRUE)/(length(progenyNames)-sapply(joinMapData[,which(!colnames(joinMapData) %in% colnames(chrEndMarkers))], function(x) sum(is.na(x))))))
  
  colnames(recombHotSpot) <- c("Marker","CHROM","POS","GenomeLoc","CrossOverCount","Frequency")
  
  ### Count recombinations at each locus
  # Loop i through progeny
  #for(i in 1:(length(joinMapDataNumeric[,1]))){
    # Loop j through markers
    #for(j in 1:(length(joinMapDataNumeric[1,])-1)){
      # For all transitions between 0 and 1, ignoring NAs
      #if((joinMapDataNumeric[i,j]==1 | joinMapDataNumeric[i,j]==0) & (joinMapDataNumeric[i,j+1]==1 | joinMapDataNumeric[i,j+1]==0)){
      #  if(joinMapDataNumeric[i,j] != joinMapDataNumeric[i,j+1]){
      #    recombHotSpot[j,5] = as.numeric(recombHotSpot[j,5]) + 1
      #  }
      #}
    #}
  #}
  
  ### Count recombinations at each locus
  # Loop i through progeny
  for(i in 1:(length(joinMapDataNumeric[,1]))){
    # Loop j through markers
    markerType <- joinMapDataNumeric[i,1]
    for(j in 2:length(joinMapDataNumeric[1,])){
      # For all transitions between 0 and 1, ignoring NAs
      if(is.na(joinMapDataNumeric[i,j])==FALSE){
        if(markerType != joinMapDataNumeric[i,j]){
          recombHotSpot[j,5] = as.numeric(recombHotSpot[j,5]) + 1
        }
        markerType <- joinMapDataNumeric[i,j]
      }
    }
  }
  
  #recomb_coverage <- merge(recombHotSpot,SNPCoverageMatrix[,c(1,2,4,5,6,7,8,275)],by=c("CHROM","POS"))
  #recomb_coverage <- recomb_coverage[order(as.numeric(as.character(recomb_coverage[,"CHROM"])),as.numeric(as.character(recomb_coverage[,"POS"])),decreasing=FALSE),]
  #joinMapData_coverage <- joinMapData_map[,which(as.numeric(as.character(recomb_coverage[,"MKK2835_PL1_2G_S11_L005"]))>3)]
  #rownames(joinMapData_coverage) <- rownames(joinMapData)
  #colnames(joinMapData_coverage) <- colnames(joinMapData)[which(as.numeric(as.character(recomb_coverage[,"MKK2835_PL1_2G_S11_L005"]))>3)]
  #write.csv(joinMapData_coverage,"MKK2835xNHP1337Markers_UniqueRecombinantProgeny_sorted.csv")
  
  
  ### Make labels for plots
  chrLabels <- rep(" ",length(recombHotSpot[1,]))
  chrLengthCount <- 1
  for(i in c("01","02","03","04","05","06","07","08","09","10","11","12","13","14")){
    chrMarkers <- which(recombHotSpot[,2]==i)
    chrLabels[chrLengthCount + (round(length(chrMarkers)/3)*2)] <- as.numeric(i)
    chrLengthCount <- chrLengthCount + length(chrMarkers)
  }

  
  png(paste0("FrequencySkewMap_",crossName1,"_filtered_All_goodcoverage_071219.png"), height=10, width=24, units="in", res=220)
  plot(recombHotSpot[,4],recombHotSpot[,6],ylim=c(0,1.0),ylab="Allele Frequency",xlab="Genome Location (Chr)",xaxt="none",pch=19,cex=1,cex.lab=2.5,cex.axis=2)
  axis(1, at=c(320425.5,1114402,2121938,3256169,4528192,5909092,7340816,8799822,10307092,11921788,13784786,15939703,18538068,21646654),tick=FALSE,seq(1,14),cex.axis=2)
  abline(v=c(640851,1587953,2655924,3856414,5199971,6618213,8063420,9536225,11077960,12765616,14803956,17075450,20000686,23292622))
  dev.off()
  
  png("FrequencySkewMap_NF54xNHP4026_filtered_gfp.png", height=10, width=24, units="in", res=220)
  
  plot(recombHotSpot[,4],recombHotSpot[,8],ylab="Allele Frequency",xlab="Genome Location (BP)",pch=19,cex=1)
  points(recombHotSpot[,4],recombHotSpot[,7],pch=1,col="red",cex=1)
  dev.off()
  
  png("FrequencySkewMap_NF54xNHP4026_filtered.png", height=10, width=24, units="in", res=220)
  
  plot(recombHotSpot[,4],recombHotSpot[,8],ylab="Allele Frequency",xlab="Genome Location (BP)",pch=19,cex=1)
  points(recombHotSpot[,4],recombHotSpot[,7],pch=1,col="red",cex=1)
  points(recombHotSpot[,4],recombHotSpot[,6],pch=1,col="blue",cex=1)
  dev.off()
  
  skewedAlleles <- recombHotSpot
  skewedAlleles$AlleleFrequency <- rowSums(joinMapData_map)/(length(joinMapData_map[1,]))
  
  png("RecombMap_NF54xNHP4026_filtered.png", height=10, width=24, units="in", res=220)
  plot(skewedAlleles[,4],skewedAlleles[,7],ylab="Allele Frequency",xlab="Genome Location (BP)",pch=19,cex=1)
  dev.off()
  
  ### Want recombinations in 1 kb bins across genome not counting chromosome ends
  chromInfo <- cbind(seq(1,14),rep(NA,14),rep(NA,14))
  colnames(chromInfo) <- c("CHROM","Start","End")
  count <- 1
  for(i in c("01","02","03","04","05","06","07","08","09","10","11","12","13","14")){
    chrMarkers <- which(recombHotSpot[,2]==i)
    chromInfo[count,2] <- min(as.numeric(recombHotSpot[chrMarkers,4]))
    chromInfo[count,3] <- max(as.numeric(recombHotSpot[chrMarkers,4]))
    count <- count + 1
  }
  
  chrCounter <- 1
  binSize <- 1000
  recombFreq <- c(99000,0)
  
  for (i in 1:14){
    markerBin <- chromInfo[i,2]
    while((markerBin + binSize) <= chromInfo[i,3]){
      recombFreq <- rbind(recombFreq,c(markerBin,sum(as.numeric(recombHotSpot[intersect(which(as.numeric(recombHotSpot[,4])>=markerBin),which(as.numeric(recombHotSpot[,4])<(markerBin+binSize))),5]))))
      markerBin <- markerBin + binSize
    }
  }
  
  png(paste0("RecombMap_",crossName1,"_filtered_5-22-19.png"), height=10, width=24, units="in", res=220)
  plot(recombFreq[,1],recombFreq[,2],ylab="Recombination Frequency",xlab="Genome Location (Chr)",xaxt="none",pch=19,cex=1,cex.lab=2.5,cex.axis=2)
  axis(1, at=c(320425.5,1114402,2121938,3256169,4528192,5909092,7340816,8799822,10307092,11921788,13784786,15939703,18538068,21646654),tick=FALSE,seq(1,14),cex.axis=2)
  abline(v=c(640851,1587953,2655924,3856414,5199971,6618213,8063420,9536225,11077960,12765616,14803956,17075450,20000686,23292622))
  dev.off()
  
  
  recombFreq <- recombHotSpot[which(recombHotSpot[,5]>0),]
  lengthRecombBlock <- recombFreq[,5]
  for(i in 2:length(recombFreq[,1])){
    if(recombFreq[i,2]==recombFreq[i-1,2]){
      lengthRecombBlock[i] <- as.numeric(recombFreq[i,4])-as.numeric(recombFreq[i-1,4])
    }
  }
 
  
  recombFreq <- cbind(recombFreq,lengthRecombBlock)
  
  
  png("Histogram of recombination event length 81.png", height=10, width=24, units="in", res=220)
  par(mar=c(8,10,2,0))
  hist(as.numeric(lengthRecombBlock)/1000,breaks=200,xlim=c(0,250),ylim=c(0,100),cex.lab=3.5,cex.axis=3.5,col="gray70",main="",axes=FALSE,xlab="",ylab="")
  abline(v=median(as.numeric(lengthRecombBlock)/1000))
  text(x=median(as.numeric(lengthRecombBlock)/1000),y=50,labels=round(median(as.numeric(lengthRecombBlock)/1000),digits=1),pos=4,cex=2.5)
  axis(1, at=seq(0,250,by=50),tick=TRUE,c("","","","","",""),cex.axis=3,line=1)
  axis(2, at=c(0,20,40,60,80,100),tick=TRUE,c(0,20,40,60,80,100),cex.axis=3,line=1)
  mtext(c(0,50,100,150,200,250),1,line=3,at=c(0,50,100,150,200,250),cex=3)
  title(xlab="Distance between markers (kb)",ylab="Frequency", line=6, cex.lab=3.5)
  
  #abline(v=median(as.numeric(lengthRecombBlock)))
  #text(x=median(as.numeric(lengthRecombBlock)),y=20,labels=round(median(as.numeric(lengthRecombBlock)),digits=0),pos=3,cex=2)
  
  dev.off()
  
  
  recombFreq[which(recombFreq[,9]>150000),]
  
  
  ### Resample N progeny
  B=1000
  lengthBetweenMarkers <- rep(0,B)
  for(k in 1:B){
    N <- 35
    progenySample <- sample(seq(2,length(joinMapDataNumeric[,1])),N,replace=FALSE)
    recombHotSpotSampled36 <- recombHotSpot[,c(1:5)]
    recombHotSpotSampled36[,5] <- rep(0,length(recombHotSpotSampled36[,1]))
  
    for(i in 1:(length(joinMapDataNumeric[1,])-1)){
      for(j in progenySample){
        if((joinMapDataNumeric[j,i] == 0 & joinMapDataNumeric[j,i+1] == 1)|(joinMapDataNumeric[j,i] == 1 & joinMapDataNumeric[j,i+1] == 0)){
          recombHotSpotSampled36[i,5] = as.numeric(recombHotSpotSampled36[i,5]) + 1
        }
      }
    }
  
    recombFreqSampled36 <- recombHotSpotSampled36[which(recombHotSpotSampled36[,5]>0),]
    lengthRecombBlockSampled36 <- recombFreqSampled36[,5]
    for(i in 2:length(recombFreqSampled36[,1])){
      if(recombFreqSampled36[i,2]==recombFreqSampled36[i-1,2]){
        lengthRecombBlockSampled36[i] <- as.numeric(recombFreqSampled36[i,4])-as.numeric(recombFreqSampled36[i-1,4])
      }
    }
  
  
    recombFreqSampled36 <- cbind(recombFreqSampled36,lengthRecombBlockSampled36)
    lengthBetweenMarkers[k] <- median(as.numeric(lengthRecombBlockSampled36))
  }
  
  lengthBetweenMarkersS35 <- lengthBetweenMarkers
  
  png("Histogram of recombination event length sampled35.png", height=10, width=24, units="in", res=220)
  par(mar=c(8,10,2,0))
  hist(as.numeric(lengthRecombBlockSampled36)/1000,breaks=200,xlim=c(0,250),ylim=c(0,100),cex.lab=3.5,cex.axis=3.5,col="gray70",main="",axes=FALSE,xlab="",ylab="")
  abline(v=median(as.numeric(lengthRecombBlockSampled36)/1000))
  text(x=median(as.numeric(lengthRecombBlockSampled36)/1000),y=50,labels=round(median(as.numeric(lengthRecombBlockSampled36)/1000),digits=1),pos=4,cex=2.5)
  axis(1, at=seq(0,250,by=50),tick=TRUE,c("","","","","",""),cex.axis=3,line=1)
  axis(2, at=c(0,20,40,60,80,100),tick=TRUE,c(0,20,40,60,80,100),cex.axis=3,line=1)
  mtext(c(0,50,100,150,200,250),1,line=3,at=c(0,50,100,150,200,250),cex=3)
  title(xlab="Distance between markers (kb)",ylab="Frequency", line=6, cex.lab=3.5)
  dev.off()
  
  
  png("Histogram of average length between markers sampled35.png", height=10, width=24, units="in", res=220)
  par(mai=c(1,1.5,1,0))
  hist(lengthBetweenMarkers/1000,breaks=20,xlab="Mean length between markers in 1000 resampled datasets",cex.lab=2,cex.axis=1.5,col="gray70",xlim=c(0,40))
  dev.off()
  
  write.csv(lengthBetweenMarkers,paste0(crossName1,"_",N,"crossoverdist35median.csv"))
  
  ### Comparison histogram of medians
  NF54xNHP4026_N35 <- read.delim("NF54xNHP4026_35crossoverdist35median.csv",sep=",")
  NF54xNHP4026_N57 <- read.delim("NF54xNHP4026_57crossoverdist57median.csv",sep=",")
  MKK2835xNHP1337_N35 <- read.delim("/Users/kbuttons/Documents/Research/P01/Data/Core C/MKK2835xNHP1337/characterization/MKK2835xNHP1337_35crossoverdistmedian.csv",sep=",")
  
  png("Histogram of average length between markers sampled35.png", height=10, width=24, units="in", res=220)
  par(mar=c(8,10,2,0))
  MKK2835xNHP1337_hist <- hist(MKK2835xNHP1337_N35[,2]/1000,breaks=20,plot=FALSE)
  NF54xNHP4026_N35_hist <- hist(NF54xNHP4026_N35[,2]/1000,breaks=20,plot=FALSE)
  NF54xNHP4026_N57_hist <- hist(NF54xNHP4026_N57[,2]/1000,breaks=10,plot=FALSE)
  plot(rep(12.305,250),seq(1,250),bty="n",type="l",lwd=3,col="black",xlim=c(10,25),ylim=c(0,250),xaxt="none",yaxt="none",xlab="",ylab="")
  plot(NF54xNHP4026_N57_hist,add=TRUE,cex.lab=2,cex.axis=1.5,col=rgb(1,0,0,0.4),xlim=c(10,25),ylim=c(0,250))
  plot(MKK2835xNHP1337_hist,add=TRUE,col=rgb(0.7,0.7,0.7,0.4))
  plot(NF54xNHP4026_N35_hist,add=TRUE,col=rgb(0,0,1,0.4))
  axis(1, at=seq(10,25,by=5),tick=TRUE,c("","","",""),cex.axis=3,line=1)
  axis(2, at=c(0,50,100,150,200,250),tick=TRUE,c(0,50,100,150,200,250),cex.axis=3,line=1)
  mtext(c(10,15,20,25),1,line=3,at=c(10,15,20,25),cex=3)
  title(xlab="Median length between markers in 1000 resampled progeny sets",ylab="Frequency", line=6, cex.lab=3.5)
  text(x=11.5,y=150,labels=round(12.305,digits=1),pos=4,cex=2.5)
  legend("topright",legend=c("NF54xNHP4026 N=35","MKK2835xNHP1337 N=35","NF54xNHP4026 N=57","MKK2835xNHP1337 N=57"),pch=c(15,15,15,NA),lty=c(NA,NA,NA,1),lwd=c(NA,NA,NA,3),col=c(rgb(0,0,1,0.4),rgb(0.7,0.7,0.7,0.4),rgb(1,0,0,0.4),"black"),cex=3)
  dev.off()
  

  
  crossOverInfo <- data.frame(t(as.numeric(c(0,0,0,1,0,0,0))))
  crossOverCounter <- 1
  
  ### Record a vector of Xover points for each individual
  # i loops over clones
  for(i in 2:length(joinMapDataNumeric[,1])){ 
    # j loops over markers
    for(j in 1:(length(joinMapDataNumeric[1,])-1)){
      if(joinMapDataNumeric[i,] != joinMapDataNumeric[i,j+1]){
        crossOverCounter <- crossOverCounter + 1
        marker <- joinMapData[j,"CHROM"]
        crossOverInfoNew <- data.frame(i,colnames(joinMapDataNumeric)[i],marker,as.numeric(as.character(recombHotSpot[which(recombHotSpot[,"Marker"]==marker),2])),as.numeric(as.character(recombHotSpot[which(recombHotSpot[,"Marker"]==marker),3])),recombHotSpot[which(recombHotSpot[,"Marker"]==marker),4],0)
        colnames(crossOverInfoNew) <- colnames(crossOverInfo)
        # If in the same chromosome, subtract distance to last cross over
        if(crossOverInfoNew[4]==crossOverInfo[crossOverCounter-1,4]){
          crossOverInfoNew[7] <- crossOverInfoNew[6]-crossOverInfo[crossOverCounter-1,6]
        }
        #If on a new chromosome, subtract distance to start of the chromosome and add a distance between last cross over and end of the previous chromosome
        if(crossOverInfoNew[4]!=crossOverInfo[crossOverCounter-1,4]){
          crossOverInfoPrev <- data.frame(c(i,colnames(joinMapData)[i],"EC",crossOverInfo[crossOverCounter-1,4:6],0))
          colnames(crossOverInfoPrev) <- colnames(crossOverInfo)
          crossOverInfoPrev[7] <- sum(chromosomeLength[which(chromosomeLength[,1]==(crossOverInfo[crossOverCounter-1,4])+1):1,2])-crossOverInfoPrev[6]
          crossOverInfo <- rbind(crossOverInfo,crossOverInfoPrev)
          crossOverInfoNew[7] <- crossOverInfoNew[6]-sum(chromosomeLength[which(chromosomeLength[,1]==crossOverInfoNew[1,4]):1,2])
          crossOverCounter <- crossOverCounter + 1
        }   
        crossOverInfo <- rbind(crossOverInfo,crossOverInfoNew)
      }
    }
  }
  
  
  
  
  ### Record a vector of Xover points for each individual
  # i loops over clones
  for(i in 2:length(joinMapData[1,])){ 
    # j loops over markers
    for(j in 1:(length(joinMapData[,1])-1)){
      if((joinMapDataNumeric[j,i] == 0 & joinMapDataNumeric[j,i+1] == 1)|(joinMapDataNumeric[j,i] == 1 & joinMapDataNumeric[j,i+1] == 0)){
        crossOverCounter <- crossOverCounter + 1
        marker <- joinMapData[j,"CHROM"]
        crossOverInfoNew <- t(data.frame(c(i,colnames(joinMapData)[i],marker,recombHotSpot[which(recombHotSpot[,"Marker"]==marker),2],as.numeric(as.character(recombHotSpot[which(recombHotSpot[,"Marker"]==marker),3])),recombHotSpot[which(recombHotSpot[,"Marker"]==marker),4],0)))
        colnames(crossOverInfoNew) <- colnames(crossOverInfo)
        # If in the same chromosome, subtract distance to last cross over
        if(crossOverInfoNew[4]==as.numeric(as.character(crossOverInfo[crossOverCounter-1,4]))){
          crossOverInfoNew[7] <- crossOverInfoNew[6]-as.numeric(as.character(crossOverInfo[crossOverCounter-1,6]))
        }
        #If on a new chromosome, subtract distance to start of the chromosome and add a distance between last cross over and end of the previous chromosome
        if(crossOverInfoNew[4]!=as.numeric(as.character(crossOverInfo[crossOverCounter-1,4]))){
          crossOverInfoPrev <- data.frame(c(i,colnames(joinMapData)[i],"EC",crossOverInfo[crossOverCounter-1,4:6],0))
          crossOverInfoPrev[7] <- chromosomeLength[which(chromosomeLength[,1]==(as.numeric(as.character(crossOverInfo[crossOverCounter,4])))+1),2]-crossOverInfoPrev[6]
          colnames(crossOverInfoPrev) <- colnames(crossOverInfo)
          crossOverInfo <- rbind(crossOverInfo,crossOverInfoPrev)
          crossOverInfoNew[7] <- crossOverInfoNew[6]-chromosomeLength[which(chromosomeLength[,1]==as.numeric(as.character(crossOverInfoNew[4]))),2]
          crossOverCounter <- crossOverCounter + 1
        }   
        crossOverInfo <- rbind(crossOverInfo,crossOverInfoNew)
      }
    }
  }
  
  
  
  markerIndex <- which(recombHotSpot[,"Marker"]==marker)
  crossOverInfo <- c(marker,recombHotSpot[which(recombHotSpot[,"Marker"]==marker),2:4])
  
  # Make substiution 0 = b and 1=a
  joinMapData <-gsub("0","b", joinMapData)
  joinMapData <- gsub("1","a", joinMapData)
  

  
  ### Markers from SNPDataSubset
  allMarkers <- vector()
  for(i in 1:14){
    allMarkers <- c(allMarkers, names(SNPDataSubset$geno[[i]]$map))
  }
  # Take subset of markers
  joinMapData <- subset(joinMapData, joinMapData[,1] %in% allMarkers)
  # Write to csv
  write.csv((joinMapData),file=filenameJoinMap)
  
  ### Create SNP Map and Save in Local Directory
  # outfile = file.path(getwd(), "linkage1.pdf")
  # lmv.linkage.plot(SNPDataSubset,outfile,mapthese=c("01","02","03","04","05","06","07"), dupnbr = TRUE)
  # outfile2 = file.path(getwd(), "linkage2.pdf")
  # lmv.linkage.plot(SNPDataSubset,outfile2,mapthese=c("08","09","10","11","12","13","14"), dupnbr = TRUE)
  
  return(components(progenyGraph)$no)
}

crossProgeny <- matrix(0,nrow=3,ncol=1)
crossProgeny[1,1] <- unique_recombinants(filename1,crossName1,cross1, parents1,growthCross1)
crossProgeny[2,1] <- unique_recombinants(filename2,crossName2, cross2, parents2)
crossProgeny[3,1] <- unique_recombinants(filenameAll,crossNameAll, crossAll, parentsAll) 
