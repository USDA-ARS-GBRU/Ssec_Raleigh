# cd into the directory containing the fastq files and demux metadata csv file before running this script

echo
echo
echo "Running $0 ..."
echo
echo

runName=$1
nThreads=${2:-4}
R1Struct=${3:-"8B12M+T"}  #  +T will read any number of remaining template bases
R2Struct=${4:-"8M+T"}

[ -f "${runName}_fgbioDemuxMetadata.csv" ] || { echo "\n\nERROR: Demux metadata file ${runName}_fgbioDemuxMetadata.csv not found\n\n" && exit 1; }

[ -d "/workdir/${USER}/tmp" ] || mkdir -p /workdir/${USER}/tmp || { echo "\n\nERROR: Could not make temp dir /workdir/$USER/tmp\n\n" && exit 1; }
TMPDIR=/workdir/${USER}/tmp

nR1Matches=$(ls *_${runName}_*_R1.fastq.gz | wc -l)
if [[ nR1Matches -eq 0 ]]; then
    echo "\n\nERROR: No matches found for *_${runName}_*_R1.fastq.gz\n\n" && exit 1;
elif [[ nR1Matches -gt 1 ]]; then
    echo "\n\nERROR: Too many matches found for *_${runName}_*_R1.fastq.gz:\n" && ls *_${runName}_*_R1.fastq.gz && exit 1;
else
    R1=$(ls *_${runName}_*_R1.fastq.gz)
fi

nR2Matches=$(ls *_${runName}_*_R2.fastq.gz | wc -l)
if [[ nR2Matches -eq 0 ]]; then
    echo "\n\nERROR: No matches found for *_${runName}_*_R2.fastq.gz\n\n" && exit 1;
elif [[ nR2Matches -gt 1 ]]; then
    echo "\n\nERROR: Too many matches found for *_${runName}_*_R2.fastq.gz:\n" && ls *_${runName}_*_R2.fastq.gz && exit 1;
else
    R2=$(ls *_${runName}_*_R2.fastq.gz)
fi

[ -d demux ] || mkdir demux
[ -d demux/${runName} ] || mkdir demux/${runName}


java -Xmx16G -Djava.io.tmpdir=$TMPDIR -jar /programs/fgbio.jar DemuxFastqs \
   --threads $nThreads \
   --inputs $R1 $R2 \
   --metadata ${runName}_fgbioDemuxMetadata.csv \
   --read-structures $R1Struct $R2Struct \
   --output demux/${runName} \
   --metrics demux/${runName}.sample_barcode_metrics.txt


echo
echo
echo "$0 is finished!"
echo

exit 0
