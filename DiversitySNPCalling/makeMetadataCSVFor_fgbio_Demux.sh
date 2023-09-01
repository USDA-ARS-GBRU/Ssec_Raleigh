# cd into the directory containing the fastq files before running this script
# copy the plate maps from SLIMs to this directory, renaming as ${PlateName}_plateMap.txt

barcodeF=${HOME}/scripts/RipTideWellToInlineBarcode.txt


plate=$1
suffix=$2
[ -f ${plate}_plateMap.txt ] || { echo "\n\nERROR: Plate map file ${plate}_plateMap.txt not found\n\n" && exit 1; }
dos2unix ${plate}_plateMap.txt


if [ "$suffix" == "" ]; then
    fgbioDemuxMetadataF=${plate}_fgbioDemuxMetadata.csv
else
    fgbioDemuxMetadataF=${plate}_${suffix}_fgbioDemuxMetadata.csv
fi

awk -v suffix=$suffix '
BEGIN {
    FS="\t"
    OFS=","
    print "Sample_ID","Sample_Name","Library_ID","Sample_Plate","Sample_Well"," Sample_Barcode","Sample_Project","Description"
    row[1] = "A"
    row[2] = "B"
    row[3] = "C"
    row[4] = "D"
    row[5] = "E"
    row[6] = "F"
    row[7] = "G"
    row[8] = "H"
}
{
    if (NR==FNR) {
        wellToBarcode[$1] = $2
    } else if (FNR>2) {
        project=$1
        plate=$3
        well=$4
        sample=$6
        outline[well] = well "," plate "_" well "_" sample (suffix?"_" suffix:"") "," "RipTide" suffix "_" plate "_" well "_" sample "," plate "," well "," wellToBarcode[well] "," project "," ""
    }
}
END {
  for (col=1; col<=12; col++) {
      for (rowNum=1; rowNum<=8; rowNum++) {
          well = row[rowNum] (col<10?"0":"") col
          if (well in outline) {
              print outline[well]
          }
      }
  }
}
' $barcodeF ${plate}_plateMap.txt > $fgbioDemuxMetadataF



