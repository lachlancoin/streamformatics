#!/bin/bash

#SBATCH --job-name=fspls_%A_%a
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32000 # mb
#SBATCH --time=100:00:00
#SBATCH --output=fspls._%A_%a.stdout
#SBATCH --error=fspls._%A_%a.stderr
#SBATCH --cpus-per-task=8

# USAGE:
# sbatch --job-name=R  runR.sh
#sbatch --array=1-10 runR.sh  covar_file  dataType options
#options can include fspls.numSamplesTrain=100   for example
##use fspls.install=TRUE for installing
## Load required modules
#module load gcc/8.3.0 
#module load openmpi/3.1.4
#module load r-bundle-bioconductor/3.9-r-3.6.0
shopt -s expand_aliases
while read line; do export $line; done </home/lcoin/.sfx
alias curlget='curl -X 'GET' --user ${USER}:${PASS}'
alias curlpost='curl -X 'POST' --user ${USER}:${PASS}'

echo $USER
echo $PASS

#module load SAMtools/1.16.1
d1=$(date +%Y%m%d%H%M%S)

n=$SLURM_ARRAY_TASK_ID
if [ ! $n ];then
  n=$1
fi

sampleID=$(tail -n +$n todo.txt | head -n 1)
outp=${sampleID}.out
cov=${sampleID}/bam/${sampleID}.cov.gz
vcf=${sampleID}/vcf/${sampleID}.targets.vcf.gz

sampleID=$(echo $vcf | rev | cut -f 1 -d '/' | rev)
mode=$2
if [ ! $mode ]; then
 mode="who"
fi

if [ $mode == "who" ]; then
echo "WHO"


cmd="curlget ${URL}/who/clear/Coin/WHO/2023?sampleID=${sampleID}"
echo $cmd
curlget ${URL}/who/clear/Coin/WHO/2023?sampleID=${sampleID}

#get resistance for VCF file
cmd="curlpost -F upload=@${vcf}  ${URL}/who/resistance/vcf/Coin/WHO/2023"
echo $cmd
curlpost -F upload=@${vcf}  ${URL}/who/resistance/vcf/Coin/WHO/2023 > ${outp}.who.tsv

##ACCESS RESULTS BY SAMPLE ID
cmd="curlpost  ${URL}/who/resistance/sample/Coin/WHO/2023?sampleID=${sampleID} >${outp}.who.tsv"
echo $cmd
curlpost  ${URL}/who/resistance/sample/Coin/WHO/2023?sampleID=${sampleID} >${outp}.who.1.tsv


fi
#exit 1;
if [ $mode == "sparsely" ]; then
##SPARSELY
echo "sparsely"
cmd="curlpost  ${URL}/sparsely/load/Coin/TB"
echo $cmd 
curlpost  ${URL}/sparsely/load/Coin/TB 
cmd="curlget ${URL}/sparsely/samples/Coin/TB"
echo $cmd
curlget ${URL}/sparsely/samples/Coin/TB


curlget ${URL}/who/clear/Coin/WHO/2023?sampleID=${sampleID}

#get resistance for VCF file
cmd="curlpost -F upload=@${vcf}  ${URL}/who/resistance/vcf/Coin/WHO/2023"
echo $cmd
#curlpost -F upload=@${vcf}  ${URL}/who/resistance/vcf/Coin/WHO/2023 > ${outp}.who.tsv


cmd="curlget ${URL}/sparsely/clear/Coin/TB?sampleID=${sampleID}"
echo $cmd
curlget ${URL}/sparsely/clear/Coin/TB?sampleID=${sampleID}

##POST VCF
cmd="curlpost -F upload=@${vcf} -F upload=@${cov}  ${URL}/sparsely/vcf/Coin/TB" 
echo $cmd
curlpost -F upload=@${vcf} -F upload=@${cov}  ${URL}/sparsely/vcf/Coin/TB  > ${outp}.pheno.tsv

##ACCESS RESULTS BY SAMPLE ID
#cmd="curlget ${URL}/sparsely/samples/Coin/TB"
echo $cmd
#curlget ${URL}/sparsely/samples/Coin/TB


cmd="curlpost  -F flags='{"type":"geno","tool":"mykrobe"}'  ${URL}/sparsely/sample/Coin/TB?sampleID=${sampleID}"
echo $cmd
curlpost  -F flags='{"type":"geno","tool":"mykrobe"}'  ${URL}/sparsely/sample/Coin/TB?sampleID=${sampleID}   > ${outp}.mykrobe.tsv


fi


#POST BIGGER





d2=$(date +%Y%m%d%H%M%S)
 
echo "$(($d2-$d1))"
