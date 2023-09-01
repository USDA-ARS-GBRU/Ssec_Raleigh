
################
# RipTide_Staug
################


# Dependencies (in order of usage):

# demux:
fastqc
fgbio.jar
makeMetadataCSVFor_fgbio_Demux.sh # custom script by Jeff G 
fgbio_DemuxFastqs_RipTide.sh  # custom script by Jeff G

# alignment and SNP calling:
bwa
picard-tools-2.26.1  # version not critcal
parallel  # GNU parallel
samtools
bcftools
R         # I used version 4.0.5, but version not critical
callLowDepthDiploidSNPs.awk # custom awk script by Jeff G 
beagle4.1
VCFtools 0.1.17
plink v1.07
WGSChrCvgBinsFromBamR1ProperPairNoDup.sh # custom script by Jeff G
SNPChrDensityFromVCF.sh # custom script by Jeff G
RipTide_Staug_Plots.R   # run locally on my Mac in RStudio


# Custom scripts from above:
makeMetadataCSVFor_fgbio_Demux.sh
fgbio_DemuxFastqs_RipTide.sh
callLowDepthDiploidSNPs.awk 
WGSChrCvgBinsFromBamR1ProperPairNoDup.sh
SNPChrDensityFromVCF.sh
RipTide_Staug_Plots.R


#######################
# RipTide_Staug demux
#######################

# on screen 4 on cbsudesktop08:
cd /data/WGS/RipTide_Staug/AAAKYL5HV
./download.sh
mkdir fastqc
fastqc -o fastqc/ 13101_1465_155956_AAAKYL5HV_Plate1_Staug_ACAGTG_R1.fastq.gz &
fastqc -o fastqc/ 13101_1465_155956_AAAKYL5HV_Plate1_Staug_ACAGTG_R2.fastq.gz &

# 1.3 billion reads!
R1: 1308227895 reads
R2: 1308227895 reads

# verify file integrity...
md5sum 13101_1465_155956_AAAKYL5HV_Plate1_Staug_ACAGTG_R1.fastq.gz &
md5sum 13101_1465_155956_AAAKYL5HV_Plate1_Staug_ACAGTG_R2.fastq.gz &
cfac145e3b7b1804781d18632c2c321e  13101_1465_155956_AAAKYL5HV_Plate1_Staug_ACAGTG_R1.fastq.gz
4cd26bf2ab6accfdbdb07be97e84dda1  13101_1465_155956_AAAKYL5HV_Plate1_Staug_ACAGTG_R2.fastq.gz

ll *.fastq.gz
-rw-r----- 1 jcg233 IGD 105787417084 Mar 29 17:13 13101_1465_155956_AAAKYL5HV_Plate1_Staug_ACAGTG_R1.fastq.gz
-rw-r----- 1 jcg233 IGD 109947002661 Mar 29 17:33 13101_1465_155956_AAAKYL5HV_Plate1_Staug_ACAGTG_R2.fastq.gz

# ...vs distrib email
Sample: Plate1_Staug_ACAGTG
File: 13101_1465_155956_AAAKYL5HV_Plate1_Staug_ACAGTG_R1.fastq.gz 
Size 105787417084 bytes, MD5: cfac145e3b7b1804781d18632c2c321e
Link: http://cbsuapps.biohpc.cornell.edu/Sequencing/showseqfile.aspx?mode=http&cntrl=337954327&refid=933038

File: 13101_1465_155956_AAAKYL5HV_Plate1_Staug_ACAGTG_R2.fastq.gz 
Size 109947002661 bytes, MD5: 4cd26bf2ab6accfdbdb07be97e84dda1
Link: http://cbsuapps.biohpc.cornell.edu/Sequencing/showseqfile.aspx?mode=http&cntrl=2104595972&refid=933039



### Made these changes to my sample demux script:
    # 4 threads by default (more than 4 threads doesn't speed it up any)
    # +T at end of R1Struct & R2Struct to match any number of remaining template bases after B and M
    # use /workdir/${USER}/tmp as a temp dir (partial bam files written there for later sorting)

# http://fulcrumgenomics.github.io/fgbio/tools/latest/DemuxFastqs.html
# threads	t	Int	The number of threads to use while de-multiplexing. The performance does not increase linearly with the # of threads and seems not to improve beyond 2-4 threads.	Optional	1	1

# in this case (NextSeq 2K), R1 and R2 are 160bp
Read structures [https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures] are made up of <number><operator> pairs much like the CIGAR string in BAM files. Four kinds of operators are recognized:
    T identifies a template read
    B identifies a sample barcode read
    M identifies a unique molecular index read
    S identifies a set of bases that should be skipped or ignored
    The last <number><operator> pair may be specified using a + sign instead of number to denote "all remaining bases". This is useful if, e.g., fastqs have been trimmed and contain reads of varying length. Both reads must have template bases. Any molecular identifiers will be concatenated using the - delimiter and placed in the given SAM record tag (RX by default). Similarly, the sample barcode bases from the given read will be placed in the BC tag.


nano ${HOME}/scripts/fgbio_DemuxFastqs_RipTide.sh
# relevant edited sections:
###
runName=$1
nThreads=${2:-4}
R1Struct=${3:-"8B12M+T"}  #  +T will read any number of remaining template bases
R2Struct=${4:-"8M+T"}
...
[ -d "/workdir/${USER}/tmp" ] || mkdir -p /workdir/${USER}/tmp || { echo "\n\nERROR: Could not make temp dir /workdir/$USER/tmp\n\n" && exit 1; }
TMPDIR=/workdir/${USER}/tmp
...
java -Xmx16G -Djava.io.tmpdir=$TMPDIR -jar /programs/fgbio.jar DemuxFastqs \
   --threads $nThreads \
   --inputs $R1 $R2 \
   --metadata ${runName}_fgbioDemuxMetadata.csv \
   --read-structures $R1Struct $R2Struct \
   --output demux/${runName} \
   --metrics demux/${runName}.sample_barcode_metrics.txt
###




# FileZilla Plate1_Staug_plateMap.txt from the SLIMS project RipTide_Staug to this folder (Windows line endings will be dealt with by dos2unix)
cd /data/WGS/RipTide_Staug/AAAKYL5HV

# for BI on cbsubi: scriptDir=/workdir/scripts/jcg233Scripts/RipTide
scriptDir=${HOME}/scripts
plate=Plate1_Staug

# make a csv barcode key metadata file for the sample level demux
bash ${scriptDir}/makeMetadataCSVFor_fgbio_Demux.sh ${plate}

# run the sample level demux
bash ${scriptDir}/fgbio_DemuxFastqs_RipTide.sh ${plate} 2>&1 | tee ${plate}_fgbio_DemuxFastqs_RipTide_$(date +%Y%m%d-%Hh%Mm%Ss).log

Running /home/jcg233/scripts/fgbio_DemuxFastqs_RipTide.sh ...
[2022/03/30 16:27:22 | DemuxFastqs | Info] Assuming input metadata file is simple CSV file.
[2022/03/30 16:27:22 | FgBioMain | Info] Executing DemuxFastqs from fgbio version 1.5.0-b01fc04-SNAPSHOT as jcg233@cbsudesktop08.biohpc.cornell.edu on JRE 13.0.2+8 with snappy, IntelInflater, and IntelDeflater
[2022/03/30 16:27:22 | DemuxFastqs | Info] Auto-detected quality format as: Standard
[2022/03/30 16:27:42 | DemuxFastqs | Info] processed     1,000,000 records.  Elapsed time: 00:00:19s.  Time for last 1,000,000:   16s.  Last read position: */*
...
[2022/03/31 21:34:49 | SamWriter | Info] Wrote     44,477,392 records.  Elapsed time: 29:07:26s.  Time for last 477,392:    7s.  Last read position: */*
[2022/03/31 21:34:51 | FgBioMain | Info] DemuxFastqs completed. Elapsed time: 1,747.54 minutes.
/home/jcg233/scripts/fgbio_DemuxFastqs_RipTide.sh is finished!
# It took 29 hours to run!

# copy from cbsudesktop08 to cbsubi.biohpc.cornell
cd /data/WGS/RipTide_Staug/
rsync -rv AAAKYL5HV cbsubi:/workdir/data/RipTide_Staug/

# make the two scripts availabe to BI (on cbsubi.biohpc.cornell):
cp ${HOME}/scripts/makeMetadataCSVFor_fgbio_Demux.sh /workdir/scripts/jcg233Scripts/RipTide/
cp ${HOME}/scripts/fgbio_DemuxFastqs_RipTide.sh /workdir/scripts/jcg233Scripts/RipTide/









#####################
# Align with BWA mem
#####################


# on cbsulm28 (screen 3)
cd /workdir/jcg233
mkdir -p RipTide_Staug/AAAKYL5HV/demux

# on cbsudesktop08 (screen 3)
cd /data/WGS/RipTide_Staug/AAAKYL5HV/demux
scp *.bam cbsulm28:/workdir/jcg233/RipTide_Staug/AAAKYL5HV/demux/


# index the genome
cd /local/workdir/jcg233/refGenomes/Stenotaphrum_secundatum
bwa index Stenotaphrum_secundatum_RaleighCultivar_genome_v1.fasta
java -jar /programs/picard-tools-2.26.1/picard.jar CreateSequenceDictionary \
      -R Stenotaphrum_secundatum_RaleighCultivar_genome_v1.fasta \
      -O Stenotaphrum_secundatum_RaleighCultivar_genome_v1.dict


outBamDir=/local/workdir/jcg233/RipTide_Staug/AAAKYL5HV/bam
[ -d $outBamDir ] || mkdir $outBamDir
tempDir=/local/workdir/jcg233/RipTide_Staug/AAAKYL5HV/tmp
[ -d $tempDir ] || mkdir $tempDir
logDir=/local/workdir/jcg233/RipTide_Staug/AAAKYL5HV/logs
[ -d $logDir ] || mkdir $logDir
cd /local/workdir/jcg233/RipTide_Staug/AAAKYL5HV/demux
myRef=/local/workdir/jcg233/refGenomes/Stenotaphrum_secundatum/Stenotaphrum_secundatum_RaleighCultivar_genome_v1.fasta
for myBam in $(ls *.bam | grep -v unmatched | grep -v BLANK); do
    outBam=$(echo $myBam | cut -d "-" -f 2,3)
    outBam=${outBam/%.bam/.aligned.bam}
    java -Xmx64G -jar /programs/picard-tools-2.26.1/picard.jar SamToFastq \
        -I $myBam \
        --FASTQ /dev/stdout \
        --INTERLEAVE true \
        --NON_PF true \
        --TMP_DIR $tempDir | \
    bwa mem -M -t 16 -p $myRef /dev/stdin | \
    java -Xmx64G -jar /programs/picard-tools-2.26.1/picard.jar MergeBamAlignment \
        --ALIGNED_BAM /dev/stdin \
        --UNMAPPED_BAM $myBam \
        --OUTPUT $outBamDir/$outBam \
        -R $myRef \
        --CREATE_INDEX true \
        --ADD_MATE_CIGAR true \
        --CLIP_ADAPTERS false \
        --CLIP_OVERLAPPING_READS true \
        --INCLUDE_SECONDARY_ALIGNMENTS true \
        --MAX_INSERTIONS_OR_DELETIONS -1 \
        --PRIMARY_ALIGNMENT_STRATEGY BestMapq \
        --ATTRIBUTES_TO_RETAIN XS \
        --TMP_DIR $tempDir
done 2>&1 | tee ${logDir}/SamToFastq_bwa_mem_MergeBamAlignment_$(date +%Y%m%d-%Hh%Mm%Ss).log





#############################################################
# dedup with and without using the pseudo molecular barcodes
#############################################################


# with pseudo molecular barcodes (screen 4)
bamDir=/local/workdir/jcg233/RipTide_Staug/AAAKYL5HV/bam
dedupBamDir=/local/workdir/jcg233/RipTide_Staug/AAAKYL5HV/dedupBam_molBar
[ -d $dedupBamDir ] || mkdir $dedupBamDir
cd $bamDir
[ -f dedupCommands_molBar.txt ] && rm dedupCommands_molBar.txt
for myBam in *.aligned.bam; do
    dedupBam=${myBam/%.aligned.bam/.dedup.bam}
    echo "java -Xmx24G -jar /programs/picard-tools-2.26.1/picard.jar MarkDuplicates --BARCODE_TAG RX --CREATE_INDEX true --INPUT $myBam --OUTPUT $dedupBamDir/$dedupBam --METRICS_FILE $dedupBamDir/${dedupBam}.metrics.txt" >> dedupCommands_molBar.txt
done

parallel -j 12 <dedupCommands_molBar.txt  # uses all cores (at peaks) and ~230GB of peak RAM


# without pseudo molecular barcodes (screen 4)
bamDir=/local/workdir/jcg233/RipTide_Staug/AAAKYL5HV/bam
dedupBamDir=/local/workdir/jcg233/RipTide_Staug/AAAKYL5HV/dedupBam
[ -d $dedupBamDir ] || mkdir $dedupBamDir
cd $bamDir
[ -f dedupCommands.txt ] && rm dedupCommands.txt
for myBam in *.aligned.bam; do
    dedupBam=${myBam/%.aligned.bam/.dedup.bam}
    echo "java -Xmx24G -jar /programs/picard-tools-2.26.1/picard.jar MarkDuplicates --CREATE_INDEX true --INPUT $myBam --OUTPUT $dedupBamDir/$dedupBam --METRICS_FILE $dedupBamDir/${dedupBam}.metrics.txt" >> dedupCommands.txt
done

parallel -j 12 <dedupCommands.txt  




dedupBamDir=/local/workdir/jcg233/RipTide_Staug/AAAKYL5HV/dedupBam_molBar
metricsF=/local/workdir/jcg233/RipTide_Staug/AAAKYL5HV/MarkDuplicatesMetrics_molBar.txt
cd $dedupBamDir
noHeader=1
for metrics in *.dedup.bam.metrics.txt; do
    [ $noHeader -eq 1 ] && noHeader=0 && grep -A 1 "^## METRICS CLASS" $metrics | tail -1 
    grep -A 2 "^## METRICS CLASS" $metrics | tail -1
done > $metricsF


dedupBamDir=/local/workdir/jcg233/RipTide_Staug/AAAKYL5HV/dedupBam
metricsF=/local/workdir/jcg233/RipTide_Staug/AAAKYL5HV/MarkDuplicatesMetrics.txt
cd $dedupBamDir
noHeader=1
for metrics in *.dedup.bam.metrics.txt; do
    [ $noHeader -eq 1 ] && noHeader=0 && grep -A 1 "^## METRICS CLASS" $metrics | tail -1 
    grep -A 2 "^## METRICS CLASS" $metrics | tail -1
done > $metricsF



cd /local/workdir/jcg233/RipTide_Staug/AAAKYL5HV/dedupBam_molBar
samtools view --require-flags DUP Plate1_Staug_A01_1-CGTACGTA.dedup.bam | wc -l
# 898238 
samtools view --require-flags DUP --exclude-flags SECONDARY,SUPPLEMENTARY Plate1_Staug_A01_1-CGTACGTA.dedup.bam | wc -l
# 898238


dedupBamDir=/local/workdir/jcg233/RipTide_Staug/AAAKYL5HV/dedupBam_molBar
cd $dedupBamDir
samtools view --exclude-flags SECONDARY,SUPPLEMENTARY Plate1_Staug_A01_1-CGTACGTA.dedup.bam | cut -f5 | sort -nr | uniq -c | less
5685661 60
  14701 59
  13079 58
  34172 57
   7206 56
   7384 55
   6218 54
   8355 53
  12932 52
  10261 51
  13178 50
  10920 49
 118890 48
  43785 47
  31372 46
  39106 45
  19315 44
  24331 43
  13416 42
  13872 41
 228766 40
  22531 39
   2691 38
  13918 37
   9384 36
   5463 35
   5073 34
  13694 33
   4178 32
   8804 31
   6662 30
   4357 29
   4733 28
 130922 27
   3269 26
  32138 25
  25783 24
  16571 23
  20361 22
  19240 21
  14176 20
  15197 19
  12058 18
   8507 17
   7906 16
  14330 15
   5052 14
  10465 13
   9498 12
   7267 11
   8026 10
  19223 9
  11925 8
  19551 7
  17235 6
  13535 5
  25372 4
  21984 3
  11076 2
  20091 1
2399134 0

samtools view --exclude-flags SECONDARY,SUPPLEMENTARY Plate1_Staug_B12_90-CTACCATC.dedup.bam | cut -f5 | sort -nr | uniq -c | less
20492103 60
  40738 59
  34095 58
  99647 57
  20294 56
  20873 55
  17438 54
  22389 53
  34896 52
  27634 51
  37041 50
  29695 49
 385117 48
 115468 47
  84820 46
  88668 45
  46483 44
  54416 43
  29954 42
  33535 41
 561750 40
  46313 39
   6122 38
  22312 37
  17279 36
  10023 35
  10301 34
  33404 33
   9167 32
  17079 31
  13556 30
   8851 29
   9128 28
 358354 27
   7250 26
  67163 25
  45837 24
  29369 23
  39115 22
  39526 21
  26520 20
  28442 19
  23187 18
  16994 17
  15502 16
  29414 15
  10655 14
  19807 13
  19297 12
  14088 11
  16312 10
  40306 9
  24054 8
  37036 7
  32945 6
  24946 5
  45118 4
  43502 3
  20661 2
  37782 1
4312243 0

samtools view --exclude-flags SECONDARY,SUPPLEMENTARY Plate1_Staug_F11_86-GAAGCTTC.dedup.bam | cut -f5 | sort -nr | uniq -c | less
38711779 60
  83975 59
  67952 58
 187929 57
  41724 56
  42990 55
  35187 54
  43409 53
  64702 52
  54893 51
  73662 50
  59290 49
 716123 48
 230939 47
 165514 46
 167818 45
  89678 44
 102247 43
  56818 42
  64284 41
1054871 40
  80181 39
  10917 38
  40768 37
  30850 36
  17241 35
  17766 34
  61950 33
  16312 32
  31393 31
  23918 30
  15753 29
  15517 28
 626044 27
  12118 26
 123944 25
  83967 24
  51249 23
  70962 22
  70734 21
  44295 20
  47204 19
  36730 18
  27280 17
  24092 16
  49852 15
  15303 14
  32954 13
  30821 12
  22247 11
  26075 10
  70277 9
  40439 8
  62341 7
  52931 6
  40048 5
  73676 4
  72584 3
  31128 2
  60859 1
7125450 0




samtools view -H Plate1_Staug_A01_1-CGTACGTA.dedup.bam | awk '$1=="@SQ" && $2~/chr/' | cut -f 2,3
SN:Stsec_v1_chr01       LN:42065675
SN:Stsec_v1_chr02       LN:53133371
SN:Stsec_v1_chr03       LN:53686315
SN:Stsec_v1_chr04       LN:38127507
SN:Stsec_v1_chr05       LN:49512229
SN:Stsec_v1_chr06       LN:38705929
SN:Stsec_v1_chr07       LN:40499537
SN:Stsec_v1_chr08       LN:70895137
SN:Stsec_v1_chr09       LN:56762739

samtools view -H Plate1_Staug_A01_1-CGTACGTA.dedup.bam | awk '$1=="@SQ" && $2~/scaf/' | wc -l
# 894
samtools view -H Plate1_Staug_A01_1-CGTACGTA.dedup.bam | awk '$1=="@SQ" && $2~/scaf/' | cut -f3 | cut -d: -f2 | awk '{sum+=$1}END{print sum}'
# 22,035,071  # do all the 894 scaffolds together


samtools view -H Plate1_Staug_A01_1-CGTACGTA.dedup.bam | awk '$1=="@SQ" && $2~/chr/' | cut -f 2,3 | awk '{split($1,chr,":");split($2,ln,":");print chr[2] "\t" "1" "\t" ln[2] > chr[2]".regions.txt"}'

samtools view -H Plate1_Staug_A01_1-CGTACGTA.dedup.bam | awk '$1=="@SQ" && $2~/scaf/' | cut -f 2,3 | awk '{split($1,chr,":");split($2,ln,":");print chr[2] "\t" "1" "\t" ln[2]}' > Stsec_v1_scaffolds.regions.txt


# on screen 3    ### next time use --output-type b ???
genome=/local/workdir/jcg233/refGenomes/Stenotaphrum_secundatum/Stenotaphrum_secundatum_RaleighCultivar_genome_v1.fasta
dedupBamDir=/local/workdir/jcg233/RipTide_Staug/AAAKYL5HV/dedupBam_molBar
vcfDir=/local/workdir/jcg233/RipTide_Staug/AAAKYL5HV/mpileup
mkdir $vcfDir
cd $dedupBamDir
bamfiles=$(ls *.bam)

for regions in *.regions.txt; do
    mpileupVCF=${regions%.regions.txt}.mpileup.vcf.gz
    echo $mpileupVCF
    bcftools mpileup --threads 4 -a AD,DP,INFO/AD,INFO/ADF,INFO/ADR --incl-flags PROPER_PAIR -q 60 --min-BQ 20 --regions-file $regions --output-type z --output $vcfDir/$mpileupVCF -f $genome $bamfiles 2>&1 | tee $vcfDir/${regions%.regions.txt}.mpileup.$(date +%Y%m%d-%Hh%Mm%Ss).log &
done


cd $vcfDir
for vcf in *.mpileup.vcf.gz; do 
    bcftools index $vcf &
done


# depth histogram for SNPs on chrs
for vcf in Stsec_v1_chr*.mpileup.vcf.gz; do
    bcftools view --types snps --output-type u $vcf | bcftools query -f '%INFO/DP\n' | awk '{if($1<10000){depth[$1]++}else{depth[10000]++}}END{for(dp=1;dp<=10000;dp++){print dp "\t" depth[dp]}}' > ${vcf%.mpileup.vcf.gz}.snp.DP.hist.txt &
done

awk '{depth[$1]+=$2}END{for(dp=1;dp<=10000;dp++){print dp "\t" depth[dp]}}' Stsec_v1_chr*.snp.DP.hist.txt > RipTide_Staug.snp.DP.hist.txt

awk '{sum+=$2}END{print sum}' RipTide_Staug.snp.DP.hist.txt
# 91711293


# in R
SNPTotDepthHist=read.delim('RipTide_Staug.snp.DP.hist.txt',header=FALSE)
colnames(SNPTotDepthHist) = c('totalDepth','nSNPs')
myPlotTitle='Total depth for 91,711,293 raw turfgrass RipTide SNPs'

pdf(file='RipTide_Staug.SNPTotDepthHist.log.pdf', title='RipTide_Staug.SNPTotDepthHist.log.pdf', height=8, width=11)
plot(SNPTotDepthHist$totalDepth,SNPTotDepthHist$nSNPs,log='x',xaxt='n',type='h',col='gray30',lend='square',xlab='Total Depth (n=94)',ylab='SNPs',xaxs="i",yaxs="i",xlim=c(10,10001),ylim=c(0,100001),main=myPlotTitle)
axis(1,at=c(10,100,1000,10000))
abline(v=c(564,1128), col=c("red", "red"), lty=c(2,2), lwd=c(1.5,1.5))
dev.off()

png(file='RipTide_Staug.SNPTotDepthHist.log.png', title='RipTide_Staug.SNPTotDepthHist.log.png',height=8,width=11,units="in",res=400)
plot(SNPTotDepthHist$totalDepth,SNPTotDepthHist$nSNPs,log='x',xaxt='n',type='h',col='gray30',lend='square',xlab='Total Depth (n=94)',ylab='SNPs',xaxs="i",yaxs="i",xlim=c(10,10001),ylim=c(0,100001),main=myPlotTitle)
axis(1,at=c(10,100,1000,10000))
abline(v=c(564,1128), col=c("red", "red"), lty=c(2,2), lwd=c(1.5,1.5))
dev.off()

pdf(file='RipTide_Staug.SNPTotDepthHist.pdf', title='RipTide_Staug.SNPTotDepthHist.pdf', height=8, width=11)
plot(SNPTotDepthHist$totalDepth,SNPTotDepthHist$nSNPs,type='h',col='gray30',lend='square',xlab='Total Depth (n=94)',ylab='SNPs',xaxs="i",yaxs="i",xlim=c(0,5001),xaxp=c(0,5000,10),ylim=c(0,100001),main=myPlotTitle)
abline(v=c(564,1128), col=c("red", "red"), lty=c(2,2), lwd=c(1.5,1.5))
text(650,95000,'6x',cex=1.5,col='red')
text(1015,95000,'12x',cex=1.5,col='red')
dev.off()




bcftools filter --SnpGap 5:indel --include 'TYPE="snp"' --output-type u Stsec_v1_chr03.mpileup.vcf.gz | bcftools view --types snps --include '(INFO/AD[0]+INFO/AD[1])>=564 & (INFO/AD[0]+INFO/AD[1])<=1128 & (INFO/AD[1]/(INFO/AD[0]+INFO/AD[1]))>0.02 & (INFO/AD[2]/(INFO/AD[0]+INFO/AD[1]))<0.02' --output-type u | bcftools query -f '[ %AD{0}/%AD{1}]\n' | less -S

bcftools filter --SnpGap 5:indel --include 'TYPE="snp"' --output-type u Stsec_v1_chr03.mpileup.vcf.gz | bcftools view --types snps --include '(INFO/AD[0]+INFO/AD[1])>=564 & (INFO/AD[0]+INFO/AD[1])<=1128 & (INFO/AD[1]/(INFO/AD[0]+INFO/AD[1]))>0.02 & (INFO/AD[2]/(INFO/AD[0]+INFO/AD[1]))<0.02'  --output-type u | bcftools query -f '[ %AD{0}/%AD{1}]\n' | bcftools query -f '[ %DP]\n' | less -S

bcftools filter --SnpGap 5:indel --include 'TYPE="snp"' --output-type u Stsec_v1_chr03.mpileup.vcf.gz | bcftools view -H --types snps --include '(INFO/AD[0]+INFO/AD[1])>=564 & (INFO/AD[0]+INFO/AD[1])<=1128 & (INFO/AD[1]/(INFO/AD[0]+INFO/AD[1]))>0.02 & (INFO/AD[2]/(INFO/AD[0]+INFO/AD[1]))<0.02' --output-type v | wc -l
# 227,765

# Rare 3rd allele must be <1% of depth of first two alleles
bcftools filter --SnpGap 5:indel --include 'TYPE="snp" & (INFO/AD[0]+INFO/AD[1])>=564 & (INFO/AD[0]+INFO/AD[1])<=1128 & (INFO/AD[1]/(INFO/AD[0]+INFO/AD[1]))>0.02 & (INFO/AD[2]/(INFO/AD[0]+INFO/AD[1]))<0.01' --output-type u Stsec_v1_chr03.mpileup.vcf.gz | bcftools view -H --output-type v | wc -l
# 225,080


# filt1:  SnpGap 5:indel; polymorphic SNPs based on INFO/AD; allow rare third alleles (i.e., from sequencing errors)
vcfDir=/local/workdir/jcg233/RipTide_Staug/AAAKYL5HV/mpileup
cd $vcfDir
mkdir filt1
for vcf in Stsec_v1_*.mpileup.vcf.gz; do
    bcftools filter --SnpGap 5:indel --include 'TYPE="snp" & (INFO/AD[0]+INFO/AD[1])>=564 & (INFO/AD[0]+INFO/AD[1])<=1128 & (INFO/AD[1]/(INFO/AD[0]+INFO/AD[1]))>0.02 & (INFO/AD[2]/(INFO/AD[0]+INFO/AD[1]))<0.01' --output-type b -o filt1/${vcf/%.mpileup.vcf.gz/.mpileup.filt1.bcf.gz} $vcf &
done


genoDir=/local/workdir/jcg233/RipTide_Staug/AAAKYL5HV/genos
mkdir $genoDir
cd /local/workdir/jcg233/RipTide_Staug/AAAKYL5HV/mpileup/filt1
for bcf in Stsec_v1_*.mpileup.filt1.bcf.gz; do
    awk -v keepINFO=0 -v genoMinDp=7 -v minHetProp=0.1 -v today="$(date)" -f ~/awk/callLowDepthDiploidSNPs.awk <(bcftools view --output-type v $bcf) 2> $genoDir/${bcf%.mpileup.filt1.bcf.gz}.callLowDepthDiploidSNPs.awk.$(date +%Y%m%d-%Hh%Mm%Ss).log | bcftools view --output-type b -o $genoDir/${bcf%.mpileup.filt1.bcf.gz}.bcf &
done

cd $genoDir
grep SNPs Stsec_v1_*.callLowDepthDiploidSNPs.awk.*.log
Stsec_v1_chr01.callLowDepthDiploidSNPs.awk.20220526-16h06m08s.log:187410 SNPs
Stsec_v1_chr02.callLowDepthDiploidSNPs.awk.20220526-16h06m08s.log:243750 SNPs
Stsec_v1_chr03.callLowDepthDiploidSNPs.awk.20220526-16h06m08s.log:225080 SNPs
Stsec_v1_chr04.callLowDepthDiploidSNPs.awk.20220526-16h06m08s.log:211905 SNPs
Stsec_v1_chr05.callLowDepthDiploidSNPs.awk.20220526-16h06m08s.log:219437 SNPs
Stsec_v1_chr06.callLowDepthDiploidSNPs.awk.20220526-16h06m08s.log:171045 SNPs
Stsec_v1_chr07.callLowDepthDiploidSNPs.awk.20220526-16h06m08s.log:187432 SNPs
Stsec_v1_chr08.callLowDepthDiploidSNPs.awk.20220526-16h06m08s.log:221507 SNPs
Stsec_v1_chr09.callLowDepthDiploidSNPs.awk.20220526-16h06m08s.log:266606 SNPs
Stsec_v1_scaffolds.callLowDepthDiploidSNPs.awk.20220526-16h06m08s.log:4670 SNPs
grep SNPs Stsec_v1_*.callLowDepthDiploidSNPs.awk.*.log | cut -d: -f2 | cut -d " " -f 1
187410
243750
225080
211905
219437
171045
187432
221507
266606
4670
grep SNPs Stsec_v1_*.callLowDepthDiploidSNPs.awk.*.log | cut -d: -f2 | cut -d " " -f 1 | awk '{sum+=$1}END{print sum}'
1938842

# missing, total genos, missing rate, call rate
grep -E 'missing|expected' Stsec_v1_*.callLowDepthDiploidSNPs.awk.*.log | cut -d: -f2 | cut -d " " -f 1 | awk '{if(NR%2==1){miss+=$1}else{sum+=$1}}END{print miss,sum,miss/sum,1-miss/sum}'
90353711 182251148 0.495765 0.504235

echo $((1938842*94))
182251148

# Nb: MAF = AC/AN;  AN=2*NS
for genos in *.bcf; do
    bcftools +fill-tags $genos --output-type b -o ${genos/%.bcf/.fill-tags.bcf} -- -t AN,NS,AC,MAF,HWE,ExcHet &
done

# F<it> = (He-Ho)/He = 1-Ho/He
for genos in *.fill-tags.bcf; do
    bcftools query -f '%CHROM\t%POS\t%INFO/DP_ar\t%INFO/MAF\t%INFO/HWE\t%INFO/ExcHet[\t%GT]\n' $genos | awk 'BEGIN{FS=OFS="\t";print "chr","pos","dp_ar","maf","hwe","excHet","n","h","ho","he","fit"}{chr=$1;pos=$2;dp_ar=$3;maf=$4;hwe=$5;excHet=$6;n=h=0;for(i=7;i<=NF;i++){if($i!="./."){n++;if($i=="0/1"||$i=="1/0"){h++}}}ho=h/n;he=2*maf*(1-maf);fit=(he>0?(1-(ho/he)):"NA");print chr,pos,dp_ar,maf,hwe,excHet,n,h,ho,he,fit}' > ${genos%.fill-tags.bcf}.Fit.txt &
done


R
fit=read.delim('Stsec_v1_chr03.Fit.txt')
pdf(file='Stsec_v1_chr03.FitHist.pdf', title='Stsec_v1_chr03.FitHist.pdf', height=8, width=11)
hist(fit$fit)
dev.off()
q()




# There are too many hets (Fis is too negatively skewed), so increase the cutoff for het calls (minHetProp) to 20%
genoDir=/local/workdir/jcg233/RipTide_Staug/AAAKYL5HV/genosMinHetProp20
mkdir $genoDir
cd /local/workdir/jcg233/RipTide_Staug/AAAKYL5HV/mpileup/filt1
for bcf in Stsec_v1_*.mpileup.filt1.bcf.gz; do
    awk -v keepINFO=0 -v genoMinDp=7 -v minHetProp=0.2 -v today="$(date)" -f ~/awk/callLowDepthDiploidSNPs.awk <(bcftools view --output-type v $bcf) 2> $genoDir/${bcf%.mpileup.filt1.bcf.gz}.callLowDepthDiploidSNPs.awk.minHetProp20.$(date +%Y%m%d-%Hh%Mm%Ss).log | bcftools view --output-type b -o $genoDir/${bcf%.mpileup.filt1.bcf.gz}.minHetProp20.bcf &
done

cd /local/workdir/jcg233/RipTide_Staug/AAAKYL5HV/genosMinHetProp20

# missing, total genos, missing rate, call rate
grep -E 'missing|expected' Stsec_v1_*.callLowDepthDiploidSNPs.awk.minHetProp20.*.log | cut -d: -f2 | cut -d " " -f 1 | awk '{if(NR%2==1){miss+=$1}else{sum+=$1}}END{print miss,sum,miss/sum,1-miss/sum}'
90353711 182251148 0.495765 0.504235

# Nb: MAF = AC/AN;  AN=2*NS
for genos in *.bcf; do
    bcftools +fill-tags $genos --output-type b -o ${genos/%.bcf/.fill-tags.bcf} -- -t AN,NS,AC,MAF,HWE,ExcHet &
done

# F<it> = (He-Ho)/He = 1-Ho/He
for genos in *.fill-tags.bcf; do
    bcftools query -f '%CHROM\t%POS\t%INFO/DP_ar\t%INFO/MAF\t%INFO/HWE\t%INFO/ExcHet[\t%GT]\n' $genos | awk 'BEGIN{FS=OFS="\t";print "chr","pos","dp_ar","maf","hwe","excHet","n","h","ho","he","fit"}{chr=$1;pos=$2;dp_ar=$3;maf=$4;hwe=$5;excHet=$6;n=h=0;for(i=7;i<=NF;i++){if($i!="./."){n++;if($i=="0/1"||$i=="1/0"){h++}}}ho=h/n;he=2*maf*(1-maf);fit=(he>0?(1-(ho/he)):"NA");print chr,pos,dp_ar,maf,hwe,excHet,n,h,ho,he,fit}' > ${genos%.fill-tags.bcf}.Fit.txt &
done

# the distribution shifted, but the mode is still -0.1 to 0.  Seems better though.
R
fit=read.delim('Stsec_v1_chr03.minHetProp20.Fit.txt')
pdf(file='Stsec_v1_chr03.minHetProp20.FitHist.pdf', title='Stsec_v1_chr03.minHetProp20.FitHist.pdf', height=8, width=11)
hist(fit$fit)
dev.off()
q()




# Filter for minimum MAF of 5%
for genos in *.fill-tags.bcf; do
    bcftools view --min-af 0.05:minor --output-type b -o ${genos/%.fill-tags.bcf/.fill-tags.maf5.bcf} $genos &
done

for genos in *.maf5.bcf; do
    chr=$(echo $genos | cut -d"." -f1 | cut -d_ -f3)
    echo $chr $(bcftools view -H $genos | wc -l)
done
cchr01 111329
chr02 148004
chr03 131339
chr04 127303
chr05 126949
chr06 102601
chr07 110697
chr08 133051
chr09 154641
scaffolds 2603

for genos in *.maf5.bcf; do
    chr=$(echo $genos | cut -d"." -f1 | cut -d_ -f3)
    echo $chr $(bcftools view -H $genos | wc -l)
done | awk '{sum+=$2}END{print "total " sum}'
total 1148517



# Filter for F<it> between -0.2 and 0.2

# F<it> = (He-Ho)/He = 1-Ho/He
for genos in *.fill-tags.maf5.bcf; do
    bcftools query -f '%CHROM\t%POS\t%INFO/DP_ar\t%INFO/MAF\t%INFO/HWE\t%INFO/ExcHet[\t%GT]\n' $genos | awk 'BEGIN{FS=OFS="\t";print "chr","pos","dp_ar","maf","hwe","excHet","n","h","ho","he","fit"}{chr=$1;pos=$2;dp_ar=$3;maf=$4;hwe=$5;excHet=$6;n=h=0;for(i=7;i<=NF;i++){if($i!="./."){n++;if($i=="0/1"||$i=="1/0"){h++}}}ho=h/n;he=2*maf*(1-maf);fit=(he>0?(1-(ho/he)):"NA");print chr,pos,dp_ar,maf,hwe,excHet,n,h,ho,he,fit}' > ${genos%.fill-tags.maf5.bcf}.maf5.Fit.txt &
done

# the distribution shifted, but the mode is still -0.1 to 0.  Seems better though.
R
fit=read.delim('Stsec_v1_chr03.minHetProp20.maf5.Fit.txt')
pdf(file='Stsec_v1_chr03.minHetProp20.maf5.FitHist.pdf', title='Stsec_v1_chr03.minHetProp20.maf5.FitHist.pdf', height=8, width=11)
hist(fit$fit,main='Histogram of F<IT> for 5% minMAF SNPs on chr03',xlab='F<IT>')
dev.off()
q()


awk '$11>-0.1 && $11<=0.0' Stsec_v1_chr03.minHetProp20.Fit.txt | wc -l
108022
awk '$11>0.0 && $11<=0.1' Stsec_v1_chr03.minHetProp20.Fit.txt | wc -l
15256

for fit in *.maf5.Fit.txt; do
    ll $fit
    awk 'NR>1 && $11>=-0.2 && $11<=0.2{print $1"\t"$2}' $fit > ${fit%.Fit.txt}.keep.txt
done

for genos in *.fill-tags.maf5.bcf; do
    bcftools view --targets-file ${genos%.fill-tags.maf5.bcf}.maf5.keep.txt --output-type z -o ${genos%.fill-tags.maf5.bcf}.maf5.fis20.vcf.gz $genos &
done

wc -l *.keep.txt
   74678 Stsec_v1_chr01.minHetProp20.maf5.keep.txt
   98559 Stsec_v1_chr02.minHetProp20.maf5.keep.txt
   89327 Stsec_v1_chr03.minHetProp20.maf5.keep.txt
   77874 Stsec_v1_chr04.minHetProp20.maf5.keep.txt
   87158 Stsec_v1_chr05.minHetProp20.maf5.keep.txt
   69532 Stsec_v1_chr06.minHetProp20.maf5.keep.txt
   75841 Stsec_v1_chr07.minHetProp20.maf5.keep.txt
   77864 Stsec_v1_chr08.minHetProp20.maf5.keep.txt
  107881 Stsec_v1_chr09.minHetProp20.maf5.keep.txt
    1445 Stsec_v1_scaffolds.minHetProp20.maf5.keep.txt
  760159 total

for genos in *.maf5.fis20.vcf.gz; do
    bcftools view -H $genos | wc -l
done
74678
98559
89327
77874
87158
69532
75841
77864
107881
1445





# Impute with Beagle (assumes diploid)


# Nb: don't need the log file as beagle makes a *.beag.log file
mkdir beagle41
for genos in *.maf5.fis20.vcf.gz; do
    java -Xmx30G -jar /programs/beagle41/beagle41.jar nthreads=8 window=2000 overlap=500 err=0.001 ne=10000 gt=$genos out=beagle41/${genos%.vcf.gz}.beag
done 

# merge genos
cd beagle41
bcftools view Stsec_v1_chr01.minHetProp20.maf5.fis20.beag.vcf.gz --output-type v > RipTide_Staug.minHetProp20.maf5.fis20.beag.vcf
for chr in chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 scaffolds; do
    bcftools view -H Stsec_v1_${chr}.minHetProp20.maf5.fis20.beag.vcf.gz >> RipTide_Staug.minHetProp20.maf5.fis20.beag.vcf
done





# Convert filt genos to plink, get distance matrix, perform plink mds, & make MDS plot in R
vcftools --version
# VCFtools (0.1.17)
vcftools --vcf RipTide_Staug.minHetProp20.maf5.fis20.beag.vcf --plink-tped --out RipTide_Staug.minHetProp20.maf5.fis20.beag

plink --noweb --tfile RipTide_Staug.minHetProp20.maf5.fis20.beag --cluster --mds-plot 4 --out RipTide_Staug.minHetProp20.maf5.fis20.beag.mds
# |        PLINK!       |     v1.07      |   10/Aug/2009

awk '{if(NR==1){print "IID"}else{split($2,indiv,"_");print indiv[4]}}' RipTide_Staug.minHetProp20.maf5.fis20.beag.mds.mds > RipTide_Staug.minHetProp20.maf5.fis20.beag.mds.mds.iid



# do the same for the raw genos
cd /local/workdir/jcg233/RipTide_Staug/AAAKYL5HV/genosMinHetProp20
bcftools view Stsec_v1_chr01.minHetProp20.maf5.fis20.vcf.gz --output-type v > RipTide_Staug.minHetProp20.maf5.fis20.vcf
for chr in chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 scaffolds; do
    bcftools view -H Stsec_v1_${chr}.minHetProp20.maf5.fis20.vcf.gz >> RipTide_Staug.minHetProp20.maf5.fis20.vcf
done 

vcftools --vcf RipTide_Staug.minHetProp20.maf5.fis20.vcf --plink-tped --out RipTide_Staug.minHetProp20.maf5.fis20

plink --noweb --tfile RipTide_Staug.minHetProp20.maf5.fis20 --cluster --mds-plot 4 --out RipTide_Staug.minHetProp20.maf5.fis20.mds
# |        PLINK!       |     v1.07      |   10/Aug/2009
awk '{if(NR==1){print "IID"}else{split($2,indiv,"_");print indiv[4]}}' RipTide_Staug.minHetProp20.maf5.fis20.mds.mds > RipTide_Staug.minHetProp20.maf5.fis20.mds.mds.iid


mkdir stats
awk 'BEGIN{for(af=1;af<=100;af++){print af/100}}' > stats/afBinsOnePct.txt
grep -m1 "^#CHR" RipTide_Staug.minHetProp20.maf5.fis20.vcf | cut -f 10- | tr "\t" "\n" > stats/samples.txt
bcftools stats --af-bins stats/afBinsOnePct.txt --1st-allele-only --depth 1,2000,1 --samples-file stats/samples.txt RipTide_Staug.minHetProp20.maf5.fis20.vcf > stats/RipTide_Staug.minHetProp20.maf5.fis20.stats.vchk
cd stats
awk 'BEGIN{FS=OFS="\t";print "depth","nGenos","pctGenos","nSites","pctSites"}{if($1=="DP"){if($3=="<1"){$3=0}else if($3~/^>/){$3=substr($3,2);$3=$3+1};print $3,$4,$5,$6,$7}}' RipTide_Staug.minHetProp20.maf5.fis20.stats.vchk > minHetProp20.maf5.fis20.siteDP.txt

# sorted individual depth
grep "^PSC" RipTide_Staug.minHetProp20.maf5.fis20.stats.vchk  | cut -f 3,10 | sort -k 2,2n | awk '{split($1,samp,"_"); print (samp[4]<10?" ":"") samp[4],$0}'
58 Plate1_Staug_B08_58  0.7
19 Plate1_Staug_C03_19  0.8
 8 Plate1_Staug_H01_8   1.5
12 Plate1_Staug_D02_12  2.0
39 Plate1_Staug_G05_39  2.1
 1 Plate1_Staug_A01_1   2.4
40 Plate1_Staug_H05_40  2.8
 7 Plate1_Staug_G01_7   3.0
66 Plate1_Staug_B09_66  3.1
64 Plate1_Staug_H08_64  3.2
26 Plate1_Staug_B04_26  3.3
27 Plate1_Staug_C04_27  3.3
91 Plate1_Staug_C12_91  3.3
62 Plate1_Staug_F08_62  3.3
82 Plate1_Staug_B11_82  3.4
 9 Plate1_Staug_A02_9   3.6
65 Plate1_Staug_A09_65  3.6
14 Plate1_Staug_F02_14  3.6
41 Plate1_Staug_A06_41  3.8
51 Plate1_Staug_C07_51  4.0
50 Plate1_Staug_B07_50  4.5
67 Plate1_Staug_C09_67  4.5
75 Plate1_Staug_C10_75  4.5
11 Plate1_Staug_C02_11  4.6
49 Plate1_Staug_A07_49  4.7
34 Plate1_Staug_B05_34  4.7
78 Plate1_Staug_F10_78  4.8
55 Plate1_Staug_G07_55  4.8
20 Plate1_Staug_D03_20  4.9
61 Plate1_Staug_E08_61  4.9
42 Plate1_Staug_B06_42  5.1
 5 Plate1_Staug_E01_5   5.2
60 Plate1_Staug_D08_60  5.3
22 Plate1_Staug_F03_22  5.3
63 Plate1_Staug_G08_63  5.8
17 Plate1_Staug_A03_17  6.0
38 Plate1_Staug_F05_38  6.1
46 Plate1_Staug_F06_46  6.3
 3 Plate1_Staug_C01_3   6.5
83 Plate1_Staug_C11_83  6.6
80 Plate1_Staug_H10_80  6.8
94 Plate1_Staug_F12_94  7.5
 2 Plate1_Staug_B01_2   7.6
89 Plate1_Staug_A12_89  7.7
35 Plate1_Staug_C05_35  7.8
43 Plate1_Staug_C06_43  7.9
21 Plate1_Staug_E03_21  7.9
68 Plate1_Staug_D09_68  8.0
69 Plate1_Staug_E09_69  8.4
32 Plate1_Staug_H04_32  8.6
44 Plate1_Staug_D06_44  8.8
92 Plate1_Staug_D12_92  8.8
45 Plate1_Staug_E06_45  8.8
93 Plate1_Staug_E12_93  8.8
57 Plate1_Staug_A08_57  8.9
25 Plate1_Staug_A04_25  9.0
90 Plate1_Staug_B12_90  9.0
81 Plate1_Staug_A11_81  9.2
73 Plate1_Staug_A10_73  9.3
59 Plate1_Staug_C08_59  9.4
85 Plate1_Staug_E11_85  9.4
10 Plate1_Staug_B02_10  9.6
48 Plate1_Staug_H06_48  10.1
18 Plate1_Staug_B03_18  10.4
76 Plate1_Staug_D10_76  10.6
 4 Plate1_Staug_D01_4   10.7
37 Plate1_Staug_E05_37  10.7
16 Plate1_Staug_H02_16  11.2
88 Plate1_Staug_H11_88  11.2
54 Plate1_Staug_F07_54  11.3
36 Plate1_Staug_D05_36  11.4
23 Plate1_Staug_G03_23  11.9
31 Plate1_Staug_G04_31  12.1
74 Plate1_Staug_B10_74  12.3
33 Plate1_Staug_A05_33  13.0
84 Plate1_Staug_D11_84  13.1
72 Plate1_Staug_H09_72  13.2
47 Plate1_Staug_G06_47  13.3
53 Plate1_Staug_E07_53  13.4
87 Plate1_Staug_G11_87  13.5
56 Plate1_Staug_H07_56  14.0
77 Plate1_Staug_E10_77  14.6
28 Plate1_Staug_D04_28  14.7
79 Plate1_Staug_G10_79  15.3
70 Plate1_Staug_F09_70  15.4
71 Plate1_Staug_G09_71  15.6
52 Plate1_Staug_D07_52  16.0
15 Plate1_Staug_G02_15  16.3
29 Plate1_Staug_E04_29  16.7
24 Plate1_Staug_H03_24  17.0
 6 Plate1_Staug_F01_6   17.8
30 Plate1_Staug_F04_30  18.2
13 Plate1_Staug_E02_13  18.4
86 Plate1_Staug_F11_86  18.5

grep "^PSC" RipTide_Staug.minHetProp20.maf5.fis20.stats.vchk  | cut -f 3,10 | sort -k 2,2nr | awk 'BEGIN{print "sample\tdepth"}{split($1,samp,"_"); print samp[4] "\t" $2}' > minHetProp20.maf5.fis20.sampDP.sorted.txt



# bargraph of number of reads per sample
cd /workdir/jcg233/RipTide_Staug/AAAKYL5HV/demux
R
barcode_metrics = read.delim('Plate1_Staug.sample_barcode_metrics.txt')
names(barcode_metrics)
 [1] "barcode_name"                         
 [2] "library_name"                         
 [3] "barcode"                              
 [4] "templates"                            
 [5] "pf_templates"                         
 [6] "perfect_matches"                      
 [7] "pf_perfect_matches"                   
 [8] "one_mismatch_matches"                 
 [9] "pf_one_mismatch_matches"              
[10] "q20_bases"                            
[11] "q30_bases"                            
[12] "total_number_of_bases"                
[13] "fraction_matches"                     
[14] "ratio_this_barcode_to_best_barcode"   
[15] "pf_fraction_matches"                  
[16] "pf_ratio_this_barcode_to_best_barcode"
[17] "pf_normalized_matches"                
[18] "frac_q20_bases"                       
[19] "frac_q30_bases"                       


barcodeNameSplit = as.data.frame(do.call(rbind, strsplit(barcode_metrics$barcode_name,'_')))
colnames(barcodeNameSplit)=c('plate1','who','well','sample')
barcodeNameSplit$plate = paste(barcodeNameSplit$plate1,barcodeNameSplit$who,sep="_")
barcodeNameSplit = barcodeNameSplit[,c("plate","well","sample")]
barcode_metrics = cbind(barcodeNameSplit,barcode_metrics)

pdf(file='RipTide_Staug_readsPerSampHist.pdf', height=8, width=11)
barplot(barcode_metrics$pf_templates,names.arg=barcode_metrics$sample,main='Total Reads per Sample',xlab='sample',ylab='read pairs',las=2,cex.names=0.5)
dev.off()

pdf(file='RipTide_Staug_readsPerSampHist.pdf', height=8, width=11)
x <- barplot(barcode_metrics$pf_templates,main='Total Reads per Sample',ylab='read pairs',xaxs='i')
axis(1,x,labels=barcode_metrics$sample,tick=FALSE,xlab='sample',las=2,cex.axis=0.5,mgp=c(3,0,2))
dev.off()

barcode_metrics_sorted = barcode_metrics[order(barcode_metrics$pf_templates,decreasing=TRUE),]
pdf(file='RipTide_Staug_readsPerSampHistSorted.pdf', height=8, width=11)
x <- barplot(barcode_metrics_sorted$pf_templates,main='Total Reads per Sample (sorted)',ylab='read pairs',xaxs='i')
axis(1,x,labels=barcode_metrics_sorted$sample,tick=FALSE,xlab='sample',las=2,cex.axis=0.5,mgp=c(3,0,2))
dev.off()

dim(barcode_metrics_sorted)
barcode_metrics_sorted_samplesOnly = barcode_metrics_sorted[-c(1,96,97),]
pdf(file='RipTide_Staug_readsPerSampHistSamplesOnlySorted.pdf', height=8, width=11)
x <- barplot(barcode_metrics_sorted_samplesOnly$pf_templates,main='Total Reads per Sample (sorted)',ylab='read pairs',xaxs='i')
axis(1,x,labels=barcode_metrics_sorted_samplesOnly$sample,tick=FALSE,xlab='sample',las=2,cex.axis=0.5,mgp=c(3,0,2))
dev.off()






# bargraph of site depth (screen 8)
cd /workdir/jcg233/RipTide_Staug/AAAKYL5HV/genosMinHetProp20/stats
R
sampDP = read.delim('minHetProp20.maf5.fis20.sampDP.sorted.txt')
pdf(file='RipTide_Staug_sampDPSorted.pdf', height=8, width=11)
x <- barplot(sampDP$depth,main='Average depth per sample in 760,159 final SNPs (sorted)',ylab='depth',xaxs='i')
axis(1,x,labels=sampDP$sample,tick=FALSE,xlab='sample',las=2,cex.axis=0.5,mgp=c(3,0,2))
dev.off()



# back to the other R session on screen 7 with wd = cd /workdir/jcg233/RipTide_Staug/AAAKYL5HV/demux
sampDP = read.delim('/workdir/jcg233/RipTide_Staug/AAAKYL5HV/genosMinHetProp20/stats/minHetProp20.maf5.fis20.sampDP.sorted.txt')
readsSampDP = merge(sampDP,barcode_metrics_sorted_samplesOnly[,c('sample','pf_templates')])
pdf(file='RipTide_Staug_sampDPVsTotalReads.pdf', height=8, width=10)
plot(readsSampDP$pf_templates,
     readsSampDP$depth,
     main="Average depth per sample in 760,159 final SNPs vs. num read pairs per sample",
     xlab="read pairs",
     ylab="average depth"
)
text(2.75e+07,7.45,adj=c(0,0),'2',col='red',cex=1)
dev.off()

# the sample with the most depth has only average depth!
readsSampDP[readsSampDP$pf_templates > 2.5e+07,]
   sample depth pf_templates
2       2   7.6     28131275   ****
6       6  17.8     27821346
24     24  17.0     26028309
29     29  16.7     26135742
30     30  18.2     27856848
86     86  18.5     25028475




# read pair density by chr histograms
bamDir=/workdir/jcg233/RipTide_Staug/AAAKYL5HV/dedupBam_molBar
cd $bamDir
mkdir CovPlots
cd $bamDir/CovPlots
species=Stenotaphrum_secundatum
bam=Plate1_Staug_A04_25-ATCCGGTA.dedup.bam
bam=Plate1_Staug_B12_90-CTACCATC.dedup.bam
$HOME/scripts/WGSChrCvgBinsFromBamR1ProperPairNoDup.sh $species 100000 30 $bamDir $bam





# SNP density by chr (before and after filt1, etc)
cp $HOME/scripts/WGSChrCvgBinsFromBamR1ProperPairNoDup.sh $HOME/scripts/SNPChrDensityFromVCF.sh
nano $HOME/scripts/SNPChrDensityFromVCF.sh  # modify appropriately

cd /local/workdir/jcg233/RipTide_Staug/AAAKYL5HV/mpileup/filt1
chrLenDir=/workdir/jcg233/RipTide_Staug/AAAKYL5HV/dedupBam_molBar/CovPlots
chrLengths=Stenotaphrum_secundatum.chrLengths.txt
chrNumToLength=Stenotaphrum_secundatum.chrNumToLength.txt
for chr in {1..9}; do
    echo ${chr}
    $HOME/scripts/SNPChrDensityFromVCF.sh $chrLenDir/$chrLengths $chrLenDir/$chrNumToLength 100000 Stsec_v1_chr0${chr}.mpileup.filt1.bcf.gz &
done

[ -f Stsec_v1.mpileup.filt1_coverage100KBins.txt ] && rm Stsec_v1.mpileup.filt1_coverage100KBins.txt
for chr in {1..9}; do
    awk -v chr=$chr '$1==chr' Stsec_v1_chr0${chr}.mpileup.filt1_coverage100KBins.txt >> Stsec_v1.mpileup.filt1_coverage100KBins.txt
done






# compare number of reads per sample to DNA concs (Jing and Ashley's)
cd /Users/jcg233/Documents/GDF/WGS/RipTide_Staug
awk '
BEGIN{
    FS=OFS="\t"
    for(r=1;r<=8;r++){
        rowToRowLetter[r]=sprintf("%c",r+64)
    }
    print "well", "concJing"
}
{
    if ($1=="conc") {conc = 1; row = 0}
    if (conc) {
        if (row > 0) {
            for (col=1;col<=12;col++) {
                print rowToRowLetter[row] (col<10?"0":"") col , $(col+1)
            }
        }
        row++
    }
}
' riptide_DNAQuant_Jing.txt > riptide_DNAQuant_Jing_byWell.txt





cd /workdir/jcg233/RipTide_Staug/AAAKYL5HV/demux
R
barcode_metrics = read.delim('Plate1_Staug.sample_barcode_metrics.txt')
names(barcode_metrics)
 [1] "barcode_name"                         
 [2] "library_name"                         
 [3] "barcode"                              
 [4] "templates"                            
 [5] "pf_templates"                         
 [6] "perfect_matches"                      
 [7] "pf_perfect_matches"                   
 [8] "one_mismatch_matches"                 
 [9] "pf_one_mismatch_matches"              
[10] "q20_bases"                            
[11] "q30_bases"                            
[12] "total_number_of_bases"                
[13] "fraction_matches"                     
[14] "ratio_this_barcode_to_best_barcode"   
[15] "pf_fraction_matches"                  
[16] "pf_ratio_this_barcode_to_best_barcode"
[17] "pf_normalized_matches"                
[18] "frac_q20_bases"                       
[19] "frac_q30_bases"                       


barcode_metrics$barcode_name = as.character(barcode_metrics$barcode_name)
barcodeNameSplit = as.data.frame(do.call(rbind, strsplit(barcode_metrics$barcode_name,'_')))
colnames(barcodeNameSplit)=c('plate1','who','well','sample')
barcodeNameSplit$plate = paste(barcodeNameSplit$plate1,barcodeNameSplit$who,sep="_")
barcodeNameSplit = barcodeNameSplit[,c("plate","well","sample")]
barcode_metrics = cbind(barcodeNameSplit,barcode_metrics)
names(barcode_metrics)
 [1] "plate"                                
 [2] "well"                                 
 [3] "sample"                               
 [4] "barcode_name"                         
 [5] "library_name"                         
 [6] "barcode"                              
 [7] "templates"                            
 [8] "pf_templates"                         
 [9] "perfect_matches"                      
[10] "pf_perfect_matches"                   
[11] "one_mismatch_matches"                 
[12] "pf_one_mismatch_matches"              
[13] "q20_bases"                            
[14] "q30_bases"                            
[15] "total_number_of_bases"                
[16] "fraction_matches"                     
[17] "ratio_this_barcode_to_best_barcode"   
[18] "pf_fraction_matches"                  
[19] "pf_ratio_this_barcode_to_best_barcode"
[20] "pf_normalized_matches"                
[21] "frac_q20_bases"                       
[22] "frac_q30_bases"                       

jing_concs = read.delim('../riptide_DNAQuant_Jing_byWell.txt')
names(jing_concs)
[1] "well"     "concJing"

ashley_concs = read.delim('../GrowthBreedingCycle_StAug_PicogreenAndNotes.txt')
names(ashley_concs)
[1] "sample"   "Genotype" "Type"     "Ploidy"   "conc"   
colnames(ashley_concs)[5]="concAshley"
barcode_metrics_concs = merge(barcode_metrics[1:94,],jing_concs,"well")
barcode_metrics_concs = merge(barcode_metrics_concs,ashley_concs,"sample")



regConcs <- lm(concAshley~concJing,data=barcode_metrics_concs) 
summary(regConcs)
(coefConcs = round(coef(regConcs),2))
(myCorTest=cor.test(barcode_metrics_concs$concJing,barcode_metrics_concs$concAshley,"two.sided"))
xLim=c(0,180)
yLim=c(0,180)
plot(barcode_metrics_concs$concJing,
     barcode_metrics_concs$concAshley,
     xlim=xLim,
     ylim=yLim,
     main="RipTide_Staug: Ashley vs. Jing DNA concentrations",
     xlab="Jing concentration (ng/ul)",
     ylab="Ashley concentration (ng/ul)",
     col='blue',
     pch=1
)
abline(regConcs,col='blue',lty=3,lwd=2)
text(100,100,paste0("corr = ",round(myCorTest$estimate,3)," (p ",format.pval(myCorTest$p.value,digits=4),")"),cex=1)
dev.copy(pdf,'AshleyVsJingConc.pdf',width=7,height=7.7)
dev.off()


regJing <- lm(templates~concJing,data=barcode_metrics_concs) 
summary(regJing)
(coefConcs = round(coef(regJing),2))
(myCorTest=cor.test(barcode_metrics_concs$concJing,barcode_metrics_concs$templates,"two.sided"))
xLim=c(0,80)
# yLim=c(0,180)
plot(barcode_metrics_concs$concJing,
     barcode_metrics_concs$templates,
     xlim=xLim,
    #  ylim=yLim,
     main="RipTide_Staug: nReadPairs vs. Jing DNA concentration",
     xlab="Jing concentration (ng/ul)",
     ylab="Read Pairs",
     col='blue',
     pch=1
)
abline(regJing,col='blue',lty=2,lwd=2)
text(60,3e+06,paste0("corr = ",round(myCorTest$estimate,3)," (p ",format.pval(myCorTest$p.value,digits=4),")"),cex=0.8)
dev.copy(pdf,'ReadsVsJingConc.pdf',width=7,height=7.6)
dev.off()

regAshley <- lm(templates~concAshley,data=barcode_metrics_concs) 
summary(regAshley)
(coefConcs = round(coef(regAshley),2))
(myCorTest=cor.test(barcode_metrics_concs$concAshley,barcode_metrics_concs$templates,"two.sided"))
xLim=c(0,200)
# yLim=c(0,180)
plot(barcode_metrics_concs$concAshley,
     barcode_metrics_concs$templates,
     xlim=xLim,
    #  ylim=yLim,
     main="RipTide_Staug: nReadPairs vs. Ashley DNA concentration",
     xlab="Ashley concentration (ng/ul)",
     ylab="Read Pairs",
     col='purple4',
     pch=1
)
abline(regAshley,col='purple4',lty=2,lwd=2)
text(175,2e+06,paste0("corr = ",round(myCorTest$estimate,3)," (p ",format.pval(myCorTest$p.value,digits=4),")"),cex=0.8)
dev.copy(pdf,'ReadsVsAshleyConc.pdf',width=7,height=7.6)
dev.off()




