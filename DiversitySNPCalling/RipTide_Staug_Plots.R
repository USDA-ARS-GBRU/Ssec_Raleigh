# library
library(ggplot2)
library(ggExtra)
library("scales")

setwd('/Users/jcg233/Documents/GDF/WGS/RipTide_Staug')
propGenome1XF='RipTide_Staug_propGenome1X.txt'
propGenome1X = read.delim(propGenome1XF)
dim(propGenome1X)
names(propGenome1X)
head(propGenome1X)
unique(propGenome1X$Plate)
colnames(propGenome1X)[2]="PlateX"
colnames(propGenome1X)[5]="ReactVol"
colnames(propGenome1X)[8]="PropGenomeCov"
propGenome1X$Plate=substr(propGenome1X$PlateX,1,nchar(propGenome1X$PlateX)-3)
unique(propGenome1X$Plate)
trimUnderscore = substr(propGenome1X$Plate,nchar(propGenome1X$Plate),nchar(propGenome1X$Plate))=="_"
unique(propGenome1X[trimUnderscore,]$Plate)
propGenome1X[trimUnderscore,]$Plate = substr(propGenome1X[trimUnderscore,]$Plate,1,nchar(propGenome1X[trimUnderscore,]$Plate)-1)
unique(propGenome1X$Plate)
names(propGenome1X)
unique(propGenome1X$ReactVol)

# make a short plate name (Plate_ReactVol) for the box plot
propGenome1X$Plate_ReactVol = "filler"
propGenome1X[grep("potato_05X",propGenome1X$PlateX),]$Plate_ReactVol = "SwPotat_05X"
propGenome1X[grep("potato_1X",propGenome1X$PlateX),]$Plate_ReactVol = "SwPotat_1X"
propGenome1X[grep("D03_05X",propGenome1X$PlateX),]$Plate_ReactVol = "Vitis1_05X"
propGenome1X[grep("D03_1X",propGenome1X$PlateX),]$Plate_ReactVol = "Vitis1_1X"
propGenome1X[grep("D04_05X",propGenome1X$PlateX),]$Plate_ReactVol = "Vitis2_05X"
propGenome1X[grep("D04_1X",propGenome1X$PlateX),]$Plate_ReactVol = "Vitis2_1X"
propGenome1X[grep("D05_05X",propGenome1X$PlateX),]$Plate_ReactVol = "Vitis3_05X"
unique(propGenome1X$Plate_ReactVol)


# grouped boxplot
bp = ggplot(propGenome1X, aes(x=Plate_ReactVol, y=PropGenomeCov, fill=ReactVol)) + 
  scale_fill_manual(values = c('steelblue2','steelblue4')) +
  geom_boxplot()
bp
ggsave('GenomeCovByPlateBoxPlot.pdf',bp,height=6,width=8)


halfReactPropGenome =  propGenome1X[propGenome1X$ReactVol=="05X",c("Sample","Plate","PropGenomeCov")]
fullReactPropGenome =  propGenome1X[propGenome1X$ReactVol=="1X",c("Sample","Plate","PropGenomeCov")]
dim(halfReactPropGenome)
dim(fullReactPropGenome)
names(halfReactPropGenome)
names(fullReactPropGenome)
names(halfReactPropGenome)[3] = "PropGenomeHalfReact"
names(fullReactPropGenome)[3] = "PropGenomeFullReact"
PropGenome = merge(halfReactPropGenome,fullReactPropGenome)
names(PropGenome)
dim(PropGenome)
grep("BLANK",PropGenome$Sample) # none (these were skipped during bash parsing of *.thresholds.bed.gz)
PropGenome$Sample # no BLANKs
(myCorTest=cor.test(PropGenome$PropGenomeFullReact,PropGenome$PropGenomeHalfReact,"two.sided"))
myLim=c(0,0.3)
plot(PropGenome$PropGenomeFullReact,PropGenome$PropGenomeHalfReact,xlab="Proportion Genome Covered - Full Reaction",ylab="Proportion Genome Covered - Half Reaction",type="p",xlim=myLim,ylim=myLim)
abline(0,1,lty=3)
text(0.05,0.28,paste0("corr = ",round(myCorTest$estimate,3)," (p ",format.pval(myCorTest$p.value,digits=6),")"),cex=1)
dev.copy(pdf,'GenomeCovHalfVsFullRxn.pdf',8,8.5)
dev.off()

names(propGenome1X)
names(propGenome1X)[4] = "ShortSampleName"
propGenome1X$Sample = paste(propGenome1X$Plate,propGenome1X$ShortSampleName,sep="_")
propGenome1X$Sample



# total reads per sample
barcodeMetricsF='CovPlots/HKHLTBGXK.sample_barcode_metrics.txt'
barcodeMetrics = read.delim(barcodeMetricsF)
colnames(barcodeMetrics)
dim(barcodeMetrics)
head(barcodeMetrics)

barcodeMetrics$Sample = substr(barcodeMetrics$barcode_name,1,nchar(barcodeMetrics$barcode_name)-3)
barcodeMetrics$ReactVol = substr(barcodeMetrics$barcode_name,nchar(barcodeMetrics$barcode_name)-2,nchar(barcodeMetrics$barcode_name))
barcodeMetrics[barcodeMetrics$ReactVol=="05X",]$Sample = substr(barcodeMetrics[barcodeMetrics$ReactVol=="05X",]$Sample,1,nchar(barcodeMetrics[barcodeMetrics$ReactVol=="05X",]$Sample)-1)
barcodeMetrics[barcodeMetrics$ReactVol=="_1X",]$ReactVol = "1X"
head(barcodeMetrics$barcode_name)
head(barcodeMetrics$Sample)
head(barcodeMetrics$ReactVol)
tail(barcodeMetrics$barcode_name)
tail(barcodeMetrics$Sample)
tail(barcodeMetrics$ReactVol)
dim(barcodeMetrics)
barcodeMetrics=barcodeMetrics[barcodeMetrics$Sample!="unmatc",]
dim(barcodeMetrics)

halfReaction = barcodeMetrics[barcodeMetrics$ReactVol=="05X",c("Sample","templates")]
fullReaction = barcodeMetrics[barcodeMetrics$ReactVol=="1X",c("Sample","templates")]
colnames(halfReaction)[2] = "nReadPairs05X"
colnames(fullReaction)[2] = "nReadPairs1X"

head(halfReaction)
head(fullReaction)
tail(halfReaction)
tail(fullReaction)

nReadPairs = merge(halfReaction,fullReaction,by="Sample")
colnames(nReadPairs)
head(nReadPairs)
dim(nReadPairs)
grep("BLANK",nReadPairs$Sample,value=TRUE) # 3 BLANKs
nReadPairs[grep("BLANK",nReadPairs$Sample),] # see the numbers for the BLANKs
(myCorTest=cor.test(nReadPairs$nReadPairs1X,nReadPairs$nReadPairs05X,"two.sided"))
myLim=c(0,1.3)
plot(nReadPairs$nReadPairs1X/1000000,nReadPairs$nReadPairs05X/1000000,xlab="Full reaction volume (millions of read pairs)",ylab="Half reaction volume (millions of read pairs)",type="p",xlim=myLim,ylim=myLim)
abline(0,1,lty=3)
text(0.2,1.1,paste0("corr = ",round(myCorTest$estimate,3)," (p ",format.pval(myCorTest$p.value,digits=4),")"),cex=1)
dev.copy(pdf,'ReadPairsHalfVsFullRxn.pdf',8,8.5)
dev.off()



names(barcodeMetrics)
names(propGenome1X)
unique(barcodeMetrics$ReactVol)
unique(propGenome1X$ReactVol)
unique(propGenome1X$PlateX)
nPairsPropGenomeCov = merge(barcodeMetrics,propGenome1X[propGenome1X$PlateX!='vDNAcad1113D05_05X',])
dim(nPairsPropGenomeCov)
names(nPairsPropGenomeCov)
nPairsPropGenomeCov = nPairsPropGenomeCov[,c('Sample','ShortSampleName','Plate','ReactVol','templates','PropGenomeCov')]
names(nPairsPropGenomeCov)
names(nPairsPropGenomeCov)[5]='ReadPairs'
head(nPairsPropGenomeCov)
tail(nPairsPropGenomeCov)
grep("BLANK",nPairsPropGenomeCov$Sample,value=TRUE) # none
plot(nPairsPropGenomeCov$ReadPairs/1000000,nPairsPropGenomeCov$PropGenomeCov,type="p",xlab="Read Pairs (millions)",ylab="Proportion Genome Covered (minMapQ30)")

plot(nPairsPropGenomeCov[nPairsPropGenomeCov$ReactVol=='05X',]$ReadPairs/1000000,nPairsPropGenomeCov[nPairsPropGenomeCov$ReactVol=='05X',]$PropGenomeCov,type="p",xlab="Read Pairs (millions)",ylab="Proportion Genome Covered (minMapQ30)",pch=3,col='red')

points(nPairsPropGenomeCov[nPairsPropGenomeCov$ReactVol=='1X',]$ReadPairs/1000000,nPairsPropGenomeCov[nPairsPropGenomeCov$ReactVol=='1X',]$PropGenomeCov,pch=1,col=alpha('blue',0.5))

p = ggplot(nPairsPropGenomeCov, aes(ReadPairs,PropGenomeCov,color=ReactVol,shape=ReactVol,size=ReactVol,alpha=ReactVol)) +
  geom_point() +
  scale_size_manual(values = c(0.5,1.5)) +
  scale_shape_manual(values = c('bullet', 'triangle open')) +
  scale_color_manual(values = c('firebrick4','steelblue')) +
  scale_alpha_manual(values = c(1,0.5))
p1 <- ggMarginal(p,type="histogram",groupColour=TRUE)
p1
ggsave('GenomeCovVsReadPairs.png',p1,dpi=600)
# ggsave('GenomeCovVsReadPairsGgsave.pdf',p1,height=6,width=8) # no difference vs dev.copy(pdf)
dev.copy(pdf,'GenomeCovVsReadPairs.pdf',8,6)  # same resolution vs png with small file, slightly tighter fonts
dev.off()


# include plate vDNAcad1113D05_05X
nPairsPropGenomeCov7Plates = merge(barcodeMetrics,propGenome1X)
dim(nPairsPropGenomeCov7Plates)
names(nPairsPropGenomeCov7Plates)
unique(nPairsPropGenomeCov7Plates$Plate_ReactVol)
nPairsPropGenomeCov7Plates = nPairsPropGenomeCov7Plates[,c('Sample','ShortSampleName','Plate','ReactVol','templates','PropGenomeCov','Plate_ReactVol')]
names(nPairsPropGenomeCov7Plates)
names(nPairsPropGenomeCov7Plates)[5]='ReadPairs'

p7 = ggplot(nPairsPropGenomeCov7Plates, aes(ReadPairs,PropGenomeCov,color=ReactVol,shape=ReactVol,size=ReactVol,alpha=ReactVol)) +
  geom_point() +
  scale_size_manual(values = c(0.5,1.5)) +
  scale_shape_manual(values = c('bullet', 'triangle open')) +
  scale_color_manual(values = c('firebrick4','steelblue')) +
  scale_alpha_manual(values = c(1,0.5))
pp7 <- ggMarginal(p7,type="histogram",groupColour=TRUE)
pp7
ggsave('GenomeCovVsReadPairs7Plates.pdf',pp7,height=6,width=8) # no difference vs dev.copy(pdf)


# grouped boxplot of nReads per sample by plate
bp = ggplot(nPairsPropGenomeCov7Plates, aes(x=Plate_ReactVol, y=ReadPairs, fill=ReactVol)) + 
  scale_fill_manual(values = c('steelblue2','steelblue4')) +
  geom_boxplot()
bp
ggsave('nReadsByPlateBoxPlot.pdf',bp,height=6,width=8)








markDuplMetricsF='MarkDuplicatesMetrics.txt'
markDuplMetrics = read.delim(markDuplMetricsF)
names(markDuplMetrics)

markDuplMetricsMolBarF='MarkDuplicatesMetrics_molBar.txt'
markDuplMetricsMolBar = read.delim(markDuplMetricsMolBarF)
names(markDuplMetricsMolBar)

keepCols = c('LIBRARY','PERCENT_DUPLICATION')
pctDup = merge(
  markDuplMetricsMolBar[,c('LIBRARY','PERCENT_DUPLICATION')],
  markDuplMetrics[,c('LIBRARY','PERCENT_DUPLICATION')],
  by.x='LIBRARY', by.y='LIBRARY'
)
names(pctDup)
names(pctDup)[2] = 'withMolBarcode'
names(pctDup)[3] = 'without'
names(pctDup)
(myCorTest=cor.test(pctDup$withMolBarcode,pctDup$without,"two.sided"))
myLim=c(0.08,0.22)
plot(pctDup$without,pctDup$withMolBarcode,xlab="without",ylab="with MolBarcode",type="p",xlim=myLim,ylim=myLim,main='Duplication Rates (with and without pseudo molecular barcodes)')
abline(0,1,lty=3)
text(0.1,0.21,paste0("corr = ",round(myCorTest$estimate,3)," (p = ",format.pval(myCorTest$p.value,digits=6),")"),cex=0.75)
dev.copy(pdf,'dupRates_withVsWithoutMolBarCodes.pdf',8,8.5)
dev.off()

# grouped boxplot
bp = ggplot(pctDup, aes(x=Plate_ReactVol, y=PropGenomeCov, fill=ReactVol)) + 
  scale_fill_manual(values = c('steelblue2','steelblue4')) +
  geom_boxplot()
bp
ggsave('GenomeCovByPlateBoxPlot.pdf',bp,height=6,width=8)










markDuplMetrics$RxnVol = "filler"
markDuplMetrics[grep("05X",markDuplMetrics$LIBRARY),]$RxnVol = "05X"
markDuplMetrics[grep("1X",markDuplMetrics$LIBRARY),]$RxnVol = "1X"
markDuplMetrics$Sample = "filler"
markDuplMetrics[grep("05X_",markDuplMetrics$LIBRARY),]$Sample = substr(markDuplMetrics[grep("05X_",markDuplMetrics$LIBRARY),]$LIBRARY,12,500)
markDuplMetrics[grep("1X_",markDuplMetrics$LIBRARY),]$Sample = substr(markDuplMetrics[grep("1X_",markDuplMetrics$LIBRARY),]$LIBRARY,11,500)
head(markDuplMetrics)
tail(markDuplMetrics)

names(markDuplMetrics)
names(barcodeMetrics)[21] = "RxnVol"
names(barcodeMetrics)

nPairsPctDupl = merge(barcodeMetrics,markDuplMetrics)
names(nPairsPctDupl)
nPairsPctDupl = nPairsPctDupl[,c("Sample","RxnVol","templates","READ_PAIR_DUPLICATES","PERCENT_DUPLICATION")]
names(nPairsPctDupl)[3] = "ReadPairs"
names(nPairsPctDupl)[4] = "DupReadPairs"
names(nPairsPctDupl)[5] = "DuplRate"
head(nPairsPctDupl)
tail(nPairsPctDupl)

markDuplMetrics[markDuplMetrics$PERCENT_DUPLICATION<0.014,] # almost all 1X vDNAcad1112D04

DvR = ggplot(nPairsPctDupl, aes(ReadPairs,DuplRate,color=RxnVol,shape=RxnVol,size=RxnVol,alpha=RxnVol)) +
  geom_point() +
  scale_size_manual(values = c(0.5,1.5)) +
  scale_shape_manual(values = c('bullet', 'triangle open')) +
  scale_color_manual(values = c('firebrick4','steelblue')) +
  scale_alpha_manual(values = c(1,0.5))
DvRh <- ggMarginal(DvR,type="histogram",groupColour=TRUE)
DvRh
ggsave('DuplRateVsReadPairs.pdf',DvRh,height=6,width=8)

DRPvR = ggplot(nPairsPctDupl, aes(ReadPairs,DupReadPairs,color=RxnVol,shape=RxnVol,size=RxnVol,alpha=RxnVol)) +
  geom_point() +
  scale_size_manual(values = c(0.5,1.5)) +
  scale_shape_manual(values = c('bullet', 'triangle open')) +
  scale_color_manual(values = c('firebrick4','steelblue')) +
  scale_alpha_manual(values = c(1,0.5))
DRPvRh <- ggMarginal(DRPvR,type="histogram",groupColour=TRUE)
DRPvRh
ggsave('DupReadPairsVsReadPairs.pdf',DRPvRh,height=6,width=8)

# plate vDNAcad1112D04 is weird. Group by plate_RxnVol
nPairsPctDupl$Plate_RxnVol = "filler"
nPairsPctDupl[grep("potato",nPairsPctDupl$Sample),]$Plate_RxnVol = "SwPotat"
nPairsPctDupl[grep("vDNAcad1111D03_",nPairsPctDupl$Sample),]$Plate_RxnVol = "Vitis1"
nPairsPctDupl[grep("vDNAcad1112D04_",nPairsPctDupl$Sample),]$Plate_RxnVol = "Vitis2"
nPairsPctDupl[grep("vDNAcad1113D05_",nPairsPctDupl$Sample),]$Plate_RxnVol = "Vitis3"
names(nPairsPctDupl)
nPairsPctDupl[nPairsPctDupl$Plate_RxnVol=="SwPotat"&nPairsPctDupl$RxnVol=="05X",]$Plate_RxnVol = "SwPotat_05X"
nPairsPctDupl[nPairsPctDupl$Plate_RxnVol=="SwPotat"&nPairsPctDupl$RxnVol=="1X",]$Plate_RxnVol = "SwPotat_1X"
nPairsPctDupl[nPairsPctDupl$Plate_RxnVol=="Vitis1"&nPairsPctDupl$RxnVol=="05X",]$Plate_RxnVol = "Vitis1_05X"
nPairsPctDupl[nPairsPctDupl$Plate_RxnVol=="Vitis1"&nPairsPctDupl$RxnVol=="1X",]$Plate_RxnVol = "Vitis1_1X"
nPairsPctDupl[nPairsPctDupl$Plate_RxnVol=="Vitis2"&nPairsPctDupl$RxnVol=="05X",]$Plate_RxnVol = "Vitis2_05X"
nPairsPctDupl[nPairsPctDupl$Plate_RxnVol=="Vitis2"&nPairsPctDupl$RxnVol=="1X",]$Plate_RxnVol = "Vitis2_1X"
nPairsPctDupl[nPairsPctDupl$Plate_RxnVol=="Vitis3"&nPairsPctDupl$RxnVol=="05X",]$Plate_RxnVol = "Vitis3_05X"
unique(nPairsPctDupl$Plate_RxnVol)
nPairsPctDupl[nPairsPctDupl$Plate_RxnVol=="Vitis3_05X",]

DRPvRByPlate = ggplot(nPairsPctDupl, aes(ReadPairs,DupReadPairs,color=Plate_RxnVol,shape=Plate_RxnVol)) +
  geom_point() +
  scale_shape_manual(values = c(0,1,2,3,5,8,6))
DRPvRByPlate
ggsave('DupReadPairsVsReadPairsByPlate.pdf',DRPvRByPlate,height=6,width=8)

DvRByPlate = ggplot(nPairsPctDupl, aes(ReadPairs,DuplRate,color=Plate_RxnVol,shape=Plate_RxnVol)) +
  geom_point() +
  scale_shape_manual(values = c(0,1,2,3,5,8,6))
DvRByPlate
ggsave('DuplRateVsReadPairsByPlate.pdf',DvRByPlate,height=6,width=8)


