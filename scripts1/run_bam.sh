#!/bin/bash

#SBATCH --job-name=fspls
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32000 # mb
#SBATCH --time=100:00:00
#SBATCH --output=fspls.stdout
#SBATCH --error=fspls.stderr
#SBATCH --cpus-per-task=8

# USAGE:
# sbatch --job-name=R  runR.sh
#sbatch --array=1-10 runR.sh  covar_file  dataType options
#options can include fspls.numSamplesTrain=100   for example
##use fspls.install=TRUE for installing
## Load required modules

##need todo.txt
#sh run_bam.sh 1  who
#sh run_bam.sh 1 sparsely

shopt -s expand_aliases

while read line; do export $line; done <~/.sfx
alias curlget='curl -X 'GET' --user ${USER}:${PASS}'
alias curlpost='curl -X 'POST' --user ${USER}:${PASS}'


module load SAMtools/1.16.1
d1=$(date +%Y%m%d%H%M%S)

n=$SLURM_ARRAY_TASK_ID
if [ ! $n ];then 
	n=$1
fi
mode=$2
if [ ! $mode ] ;then
	mode="who"
fi
#mode="who"
##WHO
#cho "WHO"
sampleID=$(tail -n +$n todo.txt | head -n 1)
bamf="./${sampleID}/bam/${sampleID}.bam"
header="header.${n}"
echo $header
outputf="${sampleID}.bam.results"
maxrecords=900000
sampleID=$(echo $bamf | rev | cut -f 1 -d '/' | rev)
echo "HERE DEETS " $sampleID $bamf $mode

if [ $mode == "who" ] ; then
echo "running WHO ${sampleID}"

curlget ${URL}/who/clear/Coin/WHO/2023?sampleID=${sampleID}
samtools view -H ${bamf} > ${header}; curlget ${URL}/who/snps/Coin/WHO/2023 | sed 's/NC_000962.3/Chromosome/g' | samtools view -F 4 -F 256 ${bamf} -L - | head -n ${maxrecords} |  split - prefix -l 4000 --filter="cat $header - | samtools view -b - | curl -k -X 'POST' --user ${USER}:${PASS} -F upload=@- ${URL}/who/resistance/bam/Coin/WHO/2023?sampleID=${sampleID}"
cmd="curlpost  ${URL}/who/resistance/sample/Coin/WHO/2023?sampleID=${sampleID} > ${outputf}.who.tsv"
echo $cmd
curlpost  ${URL}/who/resistance/sample/Coin/WHO/2023?sampleID=${sampleID} > ${outputf}.who.tsv

cmd="curlget ${URL}/who/samples/Coin/WHO/2023"
echo $cmd
curlget ${URL}/who/samples/Coin/WHO/2023
fi

#exit 1
if [ $mode == "sparsely" ]; then
echo "running sparsely"

cmd="curlpost  ${URL}/sparsely/load/Coin/TB"
echo $cmd 
curlpost  ${URL}/sparsely/load/Coin/TB 
cmd="curlget ${URL}/sparsely/samples/Coin/TB"
echo $cmd
curlget ${URL}/sparsely/samples/Coin/TB



cmd="curlget ${URL}/sparsely/clear/Coin/TB?sampleID=${sampleID}"
echo $cmd
curlget ${URL}/sparsely/clear/Coin/TB?sampleID=${sampleID}

#POST BIGGER
echo "running samtools"
samtools view -H ${bamf} > ${header};  samtools view -F 4 -F 256 ${bamf} | head -n $maxrecords | split - prefix -l 4000 --filter="cat ${header} - | samtools view -b - | curl -k -X 'POST' --user ${USER}:${PASS} -F upload=@- ${URL}/sparsely/bam/Coin/TB?sampleID=${sampleID}"
echo "done"

cmd="curlpost  -F flags='{"type":"geno","tool":"mykrobe"}'  ${URL}/sparsely/sample/Coin/TB?sampleID=${sampleID}"
echo $cmd
curlpost  -F flags='{"type":"geno","tool":"mykrobe"}'  ${URL}/sparsely/sample/Coin/TB?sampleID=${sampleID}   > ${outputf}.mykrobe.tsv
curlpost  -F flags='{"type":"pheno","tool":"mykrobe"}'  ${URL}/sparsely/sample/Coin/TB?sampleID=${sampleID}   > ${outputf}.pheno.tsv

fi



d2=$(date +%Y%m%d%H%M%S)
 
echo "$(($d2-$d1))"
