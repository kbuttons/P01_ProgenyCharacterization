filteredClonalNMSNPsHet <- filteredSNPs[,c(1:7,which(colnames(filteredSNPs) %in% c("Sample_ART-1","Sample_NF54_GFP_LUC","AB041C0_D2_S37_L002","P01_ND6C11_LPGen_S33","P01_ND6C8_LPGen_S20","P01_ND6H1_LPGen_S44")))]

### Phase SNPs with MKK2835 as 0 and NHP1337 as 1
parent1FilteredSNPs <- filteredSNPs[,parent1]
parent2FilteredSNPs <- filteredSNPs[,parent2]
for(i in colnames(filteredSNPs)[which(colnames(filteredSNPs) %in% c("Sample_ART-1","Sample_NF54_GFP_LUC","AB041C0_D2_S37_L002","P01_ND6C11_LPGen_S33","P01_ND6C8_LPGen_S20","P01_ND6H1_LPGen_S44"))]
){
  for(j in 1:length(filteredClonalNMSNPsHet[,parent1])){
    if(filteredSNPs[j,i]==parent1FilteredSNPs[j]){
      filteredClonalNMSNPsHet[j,i]<- 0
    }
    if(filteredSNPs[j,i]==parent2FilteredSNPs[j]){
      filteredClonalNMSNPsHet[j,i]<- 1
    }
  }
}
filteredClonalNMSNPsHet[which(filteredClonalNMSNPsHet=="./.")] <- "NA"
filteredClonalNMSNPsHet[which(filteredClonalNMSNPsHet=="0/.")] <- "NA"
filteredClonalNMSNPsHet[which(filteredClonalNMSNPsHet=="./0")] <- "NA"
filteredClonalNMSNPsHet[which(filteredClonalNMSNPsHet=="1/.")] <- "NA"
filteredClonalNMSNPsHet[which(filteredClonalNMSNPsHet=="./1")] <- "NA"
filteredClonalNMSNPsHet[which(filteredClonalNMSNPsHet=="2/.")] <- "NA"
filteredClonalNMSNPsHet[which(filteredClonalNMSNPsHet=="./2")] <- "NA"
filteredClonalNMSNPsHet[which(filteredClonalNMSNPsHet=="0/1")] <- 3
filteredClonalNMSNPsHet[which(filteredClonalNMSNPsHet=="1/0")] <- 3
filteredClonalNMSNPsHet[which(filteredClonalNMSNPsHet=="NA")] <- "NA"
filteredClonalNMSNPsHet[which(filteredClonalNMSNPsHet=="2/2")] <- "NA"
filteredClonalNMSNPsHet[which(filteredClonalNMSNPsHet=="2/0")] <- 3
filteredClonalNMSNPsHet[which(filteredClonalNMSNPsHet=="0/2")] <- 3
filteredClonalNMSNPsHet[which(filteredClonalNMSNPsHet=="2/1")] <- 3
filteredClonalNMSNPsHet[which(filteredClonalNMSNPsHet=="1/2")] <- 3

### Create matrix for heatmap
filteredClonalNMSNPs_map <- matrix(NA,ncol=length(filteredClonalNMSNPsHet[1,]),nrow=length(filteredClonalNMSNPsHet[,1]))
filteredClonalNMSNPs_map[which(filteredClonalNMSNPsHet=="0")] <- 0
filteredClonalNMSNPs_map[which(filteredClonalNMSNPsHet=="1")] <- 1
filteredClonalNMSNPs_map[which(filteredClonalNMSNPsHet=="3")] <- 3
filteredClonalNMSNPs_map[which(is.na(filteredClonalNMSNPsHet)==TRUE)] <- -1

colnames(filteredClonalNMSNPs_map) <- colnames(filteredClonalNMSNPsHet)

filteredClonalNMSNPs_map <- filteredClonalNMSNPs_map[order(filteredClonalNMSNPs_map[,1],as.numeric(filteredClonalNMSNPs_map[,2]),decreasing=FALSE),]

### Inheritance Heat Map 
## Make vector of colors for values below threshold
rc1 <- colorpanel(8,"white", "red")    
## Make vector of colors for values above threshold
rc2 <- colorpanel(8,"black", "blue")
rampcols <- c(rc1, rc2)
## In your example, this line sets the color for values between 49 and 51. 
#rampcols[c(-nHalf, nHalf+1)] <- rgb(t(col2rgb("green")), maxColorValue=256) 
rb1 <- seq(-1, 0, length.out=8)
rb2 <- seq(1, 3, length.out=8)[-1]
rampbreaks <- c(-4,rb1,rb2,4)


mat <- t(filteredClonalNMSNPs_map[which(filteredClonalNMSNPsHet[,1] %in% c("Pf3D7_06_v3","Pf3D7_07_v3")),-c(1,2,3,4,5,6,7)])

library(RColorBrewer)
library(gplots)
library(vegan)
library(grDevices)

png("SNPmap_NF54xNHP4026_AB041cluster.png", height=20, width=24, units="in", res=220)
par(oma=c(0,0,0,0))

heatmap.2(mat, Rowv=FALSE, Colv=FALSE, dendrogram="none", 
           trace="none", density="none", scale="none",col = rampcols, 
           breaks = rampbreaks, cexRow=1,cexCol=2,colRow="black",na.rm=TRUE,na.color="grey90",
           key=FALSE,margins=c(5, 15), srtCol=360,colCol="black",labCol=filteredClonalNMSNPsHet[which(filteredClonalNMSNPsHet[,1] %in% c("Pf3D7_06_v3","Pf3D7_07_v3")),1])

legend(0,1,c("NF54","NHP4026","mixed"),col=c("black","red","blue"),pch=c(15,15,15),bty="n",y.intersp=1,cex=1.5)

dev.off()

### Write cross files
cross1SNPs <- cbind(paste(gsub("_v3","",gsub("Pf3D7_","",filteredClonalNMSNPsHet[,1]))), filteredClonalNMSNPsHet[,2],filteredClonalNMSNPsHet[,8:13])
cross1SNPs_sorted <- cross1SNPs[order(cross1SNPs[,1],as.numeric(cross1SNPs[,2]),decreasing=FALSE),]
rownames(cross1SNPs_sorted) <- paste("M",seq(1:length(cross1SNPs_sorted[,1])),sep="",collapse=NULL)
write.csv(t(cross1SNPs_sorted),file=filename1)


filename <- filename1
crossName <- crossName1
cross <- cross1
parentsAll <- parents1
growthCross <- NULL
#OldMapFile <- "NF54xNHP4026_progenyOldMap.csv"
SNPdata_full <- read.cross("csv", file=filename, genotypes=c(0,1), alleles=c(0,1), na.strings="NA", estimate.map=FALSE)

SNPDataCSV <- read.csv(file=filename,header=FALSE, as.is=TRUE, na.strings="NA")

### Take into account previous cross and subset by those included
progenyInOldMap <- NULL

oldCrossExcluded <- NULL

clonenames_full <- SNPDataCSV[c(-1,-2,-3),1]

SNPdata <- subset(SNPdata_full, ind=(!clonenames_full %in% oldCrossExcluded))
clonenamesAll <- clonenames_full[which(!clonenames_full %in% oldCrossExcluded)]

SNPdata <- subset(SNPdata_full, ind=(!clonenames_full %in% c(parent1,parent2)))
clonenamesAll <- SNPDataCSV[c(-1,-2,-3,-4,-5),1]

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

clusterInfo <- data.frame(components(progenyGraph)$membership)
rownames(clusterInfo) <- clonenamesAll
write.csv(clusterInfo,"GenotypeClusters.csv")

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

cloningInfo <- read.delim("cloningRounds.csv",sep=",",as.is=TRUE,header=TRUE)
cloneColor <- c()
CRTSNP <- c()
for(i in 1:length(clonenamesAll)){
  cloneColor <- c(cloneColor,cloningInfo[which(cloningInfo[,1]==clonenamesAll[i]),3])
  #CRTSNP <- c(CRTSNP,cloningInfo[which(cloningInfo[,"ID"]==clonenamesAll[i]),"CRT.PCH"])
}
cloneColorRim <- cloneColor
cloneColor[which(clonenamesAll %in% c("AC051","AC052","AC031","AC035","AC036","AC037","AC001","AC002"))] <- "white"

#coords <- layout_with_fr(progenyGraphSim)
#rownames(coords) <- clonenamesAll
#layoutCoords <- cbind(coords,unlist(components(progenyGraph)$membership))
#TestCoords <- coords
#TestCoords[which(layoutCoords[,3]==42),1] <- coords[which(layoutCoords[,3]==42),1]+rnorm(15,0,1.5)
#TestCoords[which(layoutCoords[,3]==42),2] <- coords[which(layoutCoords[,3]==42),2]+rnorm(15,0,1.5)
#plot(progenyGraphSim,layout=TestCoords)

# Get coordinates from a simplified graph of the progeny and add progeny labels and cluster information
coords <- layout_with_fr(progenyGraphSim)
rownames(coords) <- clonenamesAll
layoutCoords <- cbind(coords,unlist(components(progenyGraph)$membership))


# Pull clusters that need spreading apart
coordNeeded <- which(components(progenyGraph)$csize > 5)

# Get center of clusters in coordinates


TestCoords <- coords

# Set layers by total number of clonal
verLayerSum <- c(1,8)
verLayer <- c(1,7)
i=2
while(sum(verLayer) < length(clonal)){
  verLayer <- c(verLayer,c(verLayer[i]+7))
  verLayerSum <- c(verLayerSum, sum(verLayer))
  i=i+1
}

for(i in coordNeeded) {
  NewCoords <- TestCoords[which(layoutCoords[,3]==i),1:2]
  center <- colMeans(TestCoords[which(layoutCoords[,3]==i),1:2])
  numPoints <- length(TestCoords[which(layoutCoords[,3]==i),1])
  outLayer <- which.min(abs(verLayerSum-numPoints))
  if(verLayerSum[outLayer] < numPoints){
    outLayer <- outLayer+1
  }
  numVer <- 1
  thetaVals <- 0
  radiusVals <- seq(0,outLayer,by=.9)
  verCount <- 1
  NewCoords[verCount,1] <- radiusVals[1]*cos(0)+center[1]
  NewCoords[verCount,2] <- radiusVals[1]*sin(0)+center[2]
  verCount <- verCount + 1
  for(j in 2:(outLayer)){
    if(verLayerSum[j] <= numPoints){
      thetaVals <- seq(0,(2*pi),by=(2*pi/((verLayer[j]-1))))
      NewCoords[(verLayerSum[j-1]+1):verLayerSum[j],1] <- radiusVals[j]*cos(thetaVals)+center[1]
      NewCoords[(verLayerSum[j-1]+1):verLayerSum[j],2] <- radiusVals[j]*sin(thetaVals)+center[2]
    }
    else {
      thetaVals <- seq(0,(2*pi),by=(2*pi/(verLayer[outLayer])))
      NewCoords[(verLayerSum[j-1]+1):numPoints,1] <- radiusVals[j]*cos(thetaVals[1:length((verLayerSum[j-1]+1):numPoints)])+center[1]
      NewCoords[(verLayerSum[j-1]+1):numPoints,2] <- radiusVals[j]*sin(thetaVals[1:length((verLayerSum[j-1]+1):numPoints)])+center[2]
    }
  }
  
  TestCoords[which(layoutCoords[,3]==i),1] <- NewCoords[,1]
  TestCoords[which(layoutCoords[,3]==i),2] <- NewCoords[,2]
}

plot(progenyGraphSim,layout=TestCoords)


test <- progenyGraphSim
test2 <- test + vertices(c("AC051_1","AC052_1","AC031_1","AC035_1","AC036_1","AC037_1","AC001_1","AC002_1","AC058_1","AC059_1","AC060_1"),color="white")
coords2 <- rbind(TestCoords,TestCoords[which(rownames(TestCoords) %in% c("AC051","AC052","AC031","AC035","AC036","AC037","AC001","AC002","AC058","AC059","AC060")),])
rownames(coords2)[146:156] <- c("AC051_1","AC052_1","AC031_1","AC035_1","AC036_1","AC037_1","AC001_1","AC002_1","AC058_1","AC059_1","AC060_1")


borderCol <- c(rep("#66CC66",8),"black",rep("#66CC66",11),"darkorange")
markCol <- c(rep("#66CC66",8),"gray40",rep("#66CC66",11),"darkorange")

edgeCol <- get.edgelist(progenyGraphSim)
edgeCol <- cbind(edgeCol,rep(NA,length(edgeCol[,1])))
for (i in 1:length(progenyClustersList)){
  edgeCol[which(edgeCol[,1] %in% unlist(progenyClustersList[i])),3] <- markCol[i]
}

png(paste0(crossName,"_genotype_clusters_colorcoded_newlegend_84Prog_MSadded_solid.png"), height=24, width=40, units="in", res=220)
#plot.igraph(progenyGraph,vertex.size=7,vertex.label.cex=2,mark.groups=progenyClustersList,vertex.color=cloneColor,vertex.shape=CRTSNP)
plot.igraph(progenyGraphSim,layout=TestCoords,vertex.label=rep("",length(clonenamesAll)),edge.color=edgeCol[,3],vertex.color=cloneColor,vertex.frame.color=cloneColor,vertex.size=5.5,vertex.label.cex=2.5,mark.groups=progenyClustersList,mark.shape=1,mark.border=borderCol,mark.col=markCol,mark.expand=7.5)
#plot.igraph(test2,layout=coords2,vertex.label=rep("",length(clonenamesAll)),edge.color=edgeCol[,3],vertex.color=c(cloneColorRim,"white","white","white","white","white","white","white","white","white","white","white"),vertex.frame.color=c(cloneColorRim,"white","white","white","white","white","white","white","white","white","white","white"),vertex.size=c(rep(5.5,145),3,3,3,3,3,3,3,3,3,3,3),vertex.label.cex=2.5,mark.groups=progenyClustersList,mark.shape=1,mark.border=borderCol,mark.col=markCol,mark.expand=7.5)

#plot.igraph(progenyGraphSim,vertex.label=rep("",length(clonenamesAll)),vertex.color=cloningInfo[,3],vertex.size=7,vertex.label.cex=2.5,mark.groups=progenyClustersList)
#legend("topright",c("Parental Progeny","NF54HT-GFP–luc X NHP4026 - Round 1 (0 days)","NF54HT-GFP–luc X NHP4026 - Round 2 (14 days)","NF54HT-GFP–luc X NHP4026 - Round 3 (19 days)","NF54WT X NHP4026 (5 days)"),col=c("gold","#1A237E","#039BE5","#B3E5FC","red"),pch=19,cex=3)
legend("topright",c("Cloning Round","NF54HT-GFP–luc X NHP4026","Round 1 (0 days)","NF54HT-GFP–luc X NHP4026","Round 2 (14 days)","NF54HT-GFP–luc X NHP4026","Round 3 (19 days)","NF54WT X NHP4026 (5 days)"),col=c("white","#1A237E","white","#039BE5","white","#B3E5FC","white","red"),pch=19,cex=3)
legend("bottomright",c("Genotype Clusters","Parental","Recombinant - within","cloning round","Recombinant - between","cloning rounds"),col=c("white","dark orange","#66CC66","white","gray40","white"),pch=19,cex=3,pt.cex=7)

dev.off()

plot.igraph(test2,layout=coords2,vertex.label=rep("",length(clonenamesAll)),edge.color=edgeCol[,3],vertex.color=c(cloneColorRim,"white","white","white","white","white","white","white","white","white","white","white"),vertex.frame.color=c(cloneColorRim,"white","white","white","white","white","white","white","white","white","white","white"),vertex.size=c(rep(5.5,145),3,3,3,3,3,3,3,3,3,3,3),vertex.label.cex=2.5,mark.groups=progenyClustersList,mark.shape=1,mark.border=borderCol,mark.col=markCol,mark.expand=7.5)

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
SNPDataSubset <- SNPdata

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
SNPDataCSVSubset <- SNPDataCSV
#SNPDataCSVSubset <- SNPDataCSVSubset[,c(which(!SNPDataCSV[1,] %in% snpDrop))]


SNPDataCSVCleanedSubset <- SNPDataCSV[c(1,2,3,which(SNPDataCSV[,1] %in% progenyNames)),]
write.csv(SNPDataCSVCleanedSubset,paste("AllFinal_filtered_SNPS_",crossName1,"071219.csv"))

write.csv(SNPDataCSVSubset,paste("MapFinal_filtered_SNPS_",crossName1,"071219.csv"))

#SNPDataCSVCleanedSubset <- read.cross("csv", file=paste("MapFinal_filtered_SNPS_",crossName1,"071219.csv"), genotypes=c(0,1), alleles=c(0,1), na.strings="NA", estimate.map=FALSE)
SNPDataCSVCleanedSubset <- SNPDataCSVCleanedSubset[,-1]

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
ChrEndLocMarkers <- matrix(5,ncol=(dim(SNPDataCSVSubset)[1]-3),nrow=130)
row.names(ChrEndLocMarkers) <- c(totalDistance[2:14],(totalDistance[2:14]+1),(totalDistance[2:14]+2),(totalDistance[2:14]+3),(totalDistance[2:14]+4),(totalDistance[2:14]+5),(totalDistance[2:14]+6),(totalDistance[2:14]+7),(totalDistance[2:14]+8),(totalDistance[2:14]+9))

#loop through progeny
for(i in 1:(dim(SNPDataCSVSubset)[1]-3)){
  # Inititalize chrCounter
  chrCounter <- 2
  # start with first maker as current marker
  currentMarker <- as.numeric(SNPDataCSVSubset[3+i,(which(is.na(SNPDataCSVSubset[3+i,])==FALSE)[1]+1)])
  #currentProgenySNPs <- rep(NA,dim(SNPwindowMap)[1])
  #loop through 5kb windows and pull dominant SNP
  for(j in 1:(genomeLength/windowSize)){
    currentSNPs <- SNPDataCSVSubset[3+i,intersect(which(SNPDistance>=as.numeric(row.names(SNPwindowMap)[j])),which(SNPDistance<as.numeric(row.names(SNPwindowMap)[j+1])))]
    #If no SNPs in window set to NA
    if(length(currentSNPs)==0){
      SNPwindowMap[j,i] <- currentMarker
    } else if(length(currentSNPs)>0){ #if SNPs in window make a table to count NA, 0 and 1
      SNPtable <- c(length(which(is.na(currentSNPs)==TRUE)),length(which(currentSNPs==0)),length(which(currentSNPs==1)),length(which(currentSNPs==3)))
      names(SNPtable) <- c(NA,0,1,3)
      #SNPtable <- sort(SNPtable,decreasing=TRUE)
      # If there are only NAs then set to NA
      if(SNPtable[1]==length(currentSNPs)){
        SNPwindowMap[j,i] <- currentMarker
      } else{
        SNPwindowMap[j,i] <- as.numeric(names(which.max(SNPtable[2:4])))
        currentMarker <- as.numeric(names(which.max(SNPtable[2:4])))
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

### Inheritance Heat Map 
## Make vector of colors for values below threshold
rc1 <- colorpanel(2,"grey60", "black")    
## Make vector of colors for values above threshold
rc2 <- colorpanel(2,"red", "blue")
rampcols <- c(rc1, rc2, "yellow","yellow")
## In your example, this line sets the color for values between 49 and 51. 
#rampcols[c(-nHalf, nHalf+1)] <- rgb(t(col2rgb("green")), maxColorValue=256) 
rb1 <- seq(-1, 0, length.out=2)
rb2 <- seq(1, 3, length.out=2)
rampbreaks <- c(-4,rb1,rb2,5,6)

#mat <- joinMapData_map[c(seq(3,dim(joinMapData)[1]),2,1),]
mat <- t(SNPDataChrEndMarkersSort)

png(paste0("SNPmap_",crossName1,"_filtered_coverage_mapmarkers_physloc_071219_NF54_new.png"), height=10, width=24, units="in", res=220)
par(oma=c(0,0,0,0))

heatmap.2(mat, Rowv=FALSE, Colv=FALSE, dendrogram="none", 
          trace="none", density="none", scale="none", col = rampcols, 
          breaks = rampbreaks, cexRow=3, cexCol=3, colRow=NULL, na.rm=TRUE, na.color="gray80",
          key=FALSE, margins=c(6, 15), srtCol=360, colCol=NULL, labCol=chrLabels)

legend(0,1,c(paste0(parent1," Genotype"),paste0(parent2," Genotype"),"both genotypes"),col=c("black","red","blue"),pch=c(15,15),bty="n",y.intersp=1,cex=1.5)

dev.off()

