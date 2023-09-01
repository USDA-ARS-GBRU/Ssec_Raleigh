# Assumes that the SNPs have been prefiltered using bcftools filter or bcftools view
# e.g.:  SnpGap 5:indel; polymorphic SNPs based on INFO/AD; allow rare third alleles (i.e., from sequencing errors)

# Assumes that the genos from mpileup are in the format PL:DP:AD   (PL gets dropped from the output)


# USAGE example:
# awk -v keepINFO=1 -v genoMinDp=7 -v minHetProp=0.1 -v today="$(date)" -f /path/to/callLowDepthDiploidSNPs.awk \
#     <(bcftools view --output-type v $bcf) \
#     2> $genoDir/${bcf%.mpileup.filt1.bcf.gz}.callLowDepthDiploidSNPs.awk.$(date +%Y%m%d-%Hh%Mm%Ss).log | \
#     bcftools view --output-type b -o $genoDir/${bcf%.mpileup.filt1.bcf.gz}.bcf.gz

BEGIN{
    FS=OFS="\t"
    if (genoMinDp=="") {genoMinDp=7}
    if (minHetProp=="") {minHetProp=0.1}
    if (keepINFO=="") {keepINFO=1}
}
{
    if ($1 ~ /^#/) {
        if ($1 ~ /^##FORMAT=<ID=DP/) {
            print  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        }
        if ($1 ~ /^#CHROM/) {
            print  "##awkCommand=awk -v keepINFO=" keepINFO " -v genoMinDp=" genoMinDp " -v minHetProp=" minHetProp " -f callLowDepthDiploidSNPs.awk; Date=" today
            nindiv = NF - 9
        }
        if ($1 !~ /^##INFO=/ && $1 !~ /^##ALT=/ && $1 !~ /^##FORMAT=<ID=PL/) {
            print
        }
        if ($1 ~ /^##INFO=/ && keepINFO) {
            print
        }
        if ($1 ~ /^##INFO=<ID=DP/) {
            if (!keepINFO) {print}  # always print the original DP tag
            print "##INFO=<ID=DP_ar,Number=1,Type=Integer,Description=\"Read depth in alt and ref for all samples\">"
            print "##INFO=<ID=DP_ar" genoMinDp ",Number=1,Type=Integer,Description=\"Read depth in alt and ref for samples with minDP " genoMinDp "\">"
        }
    } else {
        # print $1,$2,$4,$5,$6,annos[anno],$11,$12,$38  # testing filtering
        outline = $1
        for (i=2;i<=4;i++) {
            outline = outline "\t" $i
        }
        nAlt=split($5,alt,",")
        # I only care about the first two alleles (ref and alt1)
        outline = outline "\t" alt[1] "\t" "." "\t" "."
        
        # assumes that the genos from mpileup are in the format PL:DP:AD  # drop PL
        genoLine = "GT:DP:AD"
        dp = dpMinDp = 0;
        for (i=10;i<=NF;i++) {
            split($i,format,":")
            split(format[3],ad,",")
            # I only care about the first two alleles (ref and alt1)
            dp += ad[1]+ad[2]
            if (ad[1]+ad[2] < genoMinDp) {
                geno = "./."
                nmiss++
            } else {
                dpMinDp += ad[1]+ad[2]
                if (ad[1] >= ad[2]) {
                    if (ad[2]/(ad[1]+ad[2]) >= minHetProp) {
                        geno = "0/1"
                        nhet++
                    } else {
                        geno = "0/0"
                        nref++
                    }
                } else {
                    if (ad[1]/(ad[1]+ad[2]) >= minHetProp) {
                        geno = "0/1"
                        nhet++
                    } else {
                        geno = "1/1"
                        nalt++
                    }
                }
            }
            genoLine = genoLine "\t" geno ":" ad[1]+ad[2] ":" ad[1] "," ad[2]
        }

        if (keepINFO) {
            outline = outline "\t" $8 ";DP_ar=" dp ";DP_ar" genoMinDp "=" dpMinDp "\t" genoLine
        } else {
            nAnno=split($8,annos,";")
            for (a=1; a<= nAnno; a++) {
                if (annos[a] ~ /^DP=/) {
                    origDP = annos[a]
                    break
                }
            }
            outline = outline "\t" origDP ";DP_ar=" dp ";DP_ar" genoMinDp "=" dpMinDp "\t" genoLine
        }
        print outline
        nSNPs++
    }
}
END {
    print "" > "/dev/stderr"
    print nindiv " samples" > "/dev/stderr"
    print nSNPs " SNPs" > "/dev/stderr"
    print "" > "/dev/stderr"
    print nmiss " missing" > "/dev/stderr"
    print nhet " hets" > "/dev/stderr"
    print nref " homozygous ref" > "/dev/stderr"
    print nalt " homozygous alt" > "/dev/stderr"
    print "" > "/dev/stderr"
    print nmiss+nhet+nref+nalt " total genos" > "/dev/stderr"
    print nindiv * nSNPs " expected total genos" > "/dev/stderr"
    print "" > "/dev/stderr"
}
