library("VariantAnnotation") #load the package
library("vcfR")
library("memgene")

### Setw
setwd("/Users/kbuttons/Documents/Research/P01/Data/Core C/")


### read in file 
Cross_SNPs <- read.vcfR("Cross.HQ.filt.recal.vcf")
colnames(Cross_SNPs@gt)
Cross_SNPs@gt[,c("Sample_ART-1","Sample_NF54_GFP_LUC")]

### Put these two columns into vectors 
Parent1_GT <- rep(-1, length(Cross_SNPs@gt[,"Sample_ART-1"]))
Parent2_GT <- rep(-1, length(Cross_SNPs@gt[,"Sample_NF54_GFP_LUC"]))

### loop through the entire vector and insert the 0/0 or 0/1 into vectors
for(i in 1:length(Cross_SNPs@gt[,"Sample_ART-1"])){
  Parent1_GT[i] <- unlist(strsplit(Cross_SNPs@gt[i,"Sample_ART-1"],":"))[1]
  Parent2_GT[i] <- unlist(strsplit(Cross_SNPs@gt[i,"Sample_NF54_GFP_LUC"],":"))[1]
}

### SNPs where Parents 1 and 2 are homozygotes	
Parental_SNPs <- which((Parent1_GT=="0/0"|Parent1_GT=="1/1")&(Parent2_GT=="0/0"|Parent2_GT=="1/1"))


### Store parental genotypes at SNPs where parents meet this above condition
SNPmatrix <- matrix(-1,ncol=(length(colnames(Cross_SNPs@gt))+6),nrow=length(Cross_SNPs@gt[,"Sample_ART-1"]))
for(i in 1:length(Cross_SNPs@gt[,"Sample_ART-1"])){
  SNPmatrix[i,1:7] <- Cross_SNPs@fix[i,1:7]
  for(j in 2:length(colnames(Cross_SNPs@gt))){
    SNPmatrix[i,j+6] <- unlist(strsplit(Cross_SNPs@gt[i,j],"[:]"))[1]
  }
}

colnames(SNPmatrix) <- c(colnames(Cross_SNPs@fix[,1:7]),colnames(Cross_SNPs@gt[,-1]))

### subset SNPs based on homozygosity of parents
Parental_homozySNPs <- SNPmatrix[Parental_SNPs,]

### Count heterozygous calls in matrix
het_counts <- matrix(0,nrow=length(colnames(Cross_SNPs@gt[,-1])),ncol=1)
for(j in 2:length(colnames(Cross_SNPs@gt))){
  het_counts[j-1,1] <- length(which(Parental_homozySNPs[,j+6]=="0/1"|Parental_homozySNPs[,j+6]=="1/0"))
}

rownames(het_counts) <- colnames(Cross_SNPs@gt[,-1])

### histogram of heterozygous call counts
hist(het_counts,breaks=100)
rownames(het_counts)[which(het_counts>700)]

### Sliding window size of average cross over to look for above average number of mutations
### From histogram, take 130 as our average heterozygous error rate for this 
window_size = 30
sliding_window_hetz <- matrix(0,nrow=dim(Parental_homozySNPs)[1]/window_size,ncol=length(colnames(Cross_SNPs@gt[,-1])))

for(i in 1:dim(sliding_window_hetz)[1]){
  for(j in 2:length(colnames(Cross_SNPs@gt))){
    sliding_window_hetz[i,j-1] <- length(which(Parental_homozySNPs[(((i-1)*window_size)+1):(i*window_size),j+6]=="0/1"|Parental_homozySNPs[(((i-1)*window_size)+1):(i*window_size),j+6]=="1/0"))
  }
}

colnames(sliding_window_hetz) <- colnames(Cross_SNPs@gt[,-1])

hist(sliding_window_hetz,breaks=30)
boxplot(sliding_window_hetz)


### Check which columns have the maximum heterozygous counts for a window above expected
max_het <- rep(0,length(colnames(Cross_SNPs@gt[,-1])))
for(i in 1:length(colnames(Cross_SNPs@gt[,-1]))){
  max_het[i] <- max(sliding_window_hetz[,i])
}

colnames(sliding_window_hetz)[which(max_het>15)]

