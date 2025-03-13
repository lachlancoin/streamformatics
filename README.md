#SETUP
Client for accessing streamformatics

You need a file called  "~/.sfx" with following format

URL=https://api.streamformatics.org
USER={username}
PASS={pass}
API_KEY={key}


Note that you choose username and password and you need to ask for the key from the admin
##CONVENIENT ALIASES
##CAN ADD FOLLOWING TO .bashrc
#The -k option is only required for localhost
while read line; do export $line; done <~/.sfx

#alias curlget='curl -k -X 'GET' --user ${USER}:${PASS}'
#alias curlpost='curl -k -X 'POST' --user ${USER}:${PASS}'



##TEST THE ENDPONT
curlget ${URL}/admin/echo?msg=test

##CURL EXAMPLES

##REGISTER USER FOR ORG
##NEED THE KEY

curlget ${URL}/admin/register/Coin?key=${API_KEY}

##WHO ENDPOINTS

#clear results by sample ID
vcf=ERR10225529.vcf.gz

curlget ${URL}/who/clear/Coin/WHO/2023?sampleID=${vcf}

#get resistance for VCF file
curlpost -F upload=@${vcf}  ${URL}/who/resistance/vcf/Coin/WHO/2023

##ACCESS RESULTS BY SAMPLE ID
curlpost  ${URL}/who/resistance/sample/Coin/WHO/2023?sampleID=${vcf}


##QUERY BAM FOR RESISTANCE
##need to define USER, PASS and bamf
##NODE THAT THIS ONE REQUIRES SAMTOOLS INSTALLED LOCALLY
bamf=tb-s1.bam
#bamf=ERR10225529.bam
samtools view -H ${bamf} > header; curlget ${URL}/who/snps/Coin/WHO/2023 | samtools view ${bamf} -L - | split - prefix -l 4000 --filter="cat header - | samtools view -b - | curl -k -X 'POST' --user ${USER}:${PASS} -F upload=@- ${URL}/who/resistance/bam/Coin/WHO/2023?sampleID=${bamf}"

##ACCESS RESULTS BY SAMPLE ID
curlpost  ${URL}/who/resistance/sample/Coin/WHO/2023?sampleID=${bamf}



##SPARSELY ENDPOINTS
##NEED VCF AND COV FILE
vcf=ERR10225529.vcf.gz
cov=ERR10225529.bam.cov.gz
#OPTIONAL - PRELOAD DATABASE INTO MEMORY
curlpost  ${URL}/sparsely/load/Coin/TB 

##POST VCF
curlpost -F upload=@${vcf} -F upload=@${cov}  ${URL}/sparsely/vcf/Coin/TB 

##ACCESS RESULTS BY SAMPLE ID
curlpost  -F flags='{"type":"geno","tool":"mykrobe"}'  ${URL}/sparsely/sample/Coin/TB?sampleID=${vcf}

#POST SMALL BAM
bamf=tb-s1.bam
curlpost -F upload=@${bamf} ${URL}/sparsely/bam/Coin/TB?sampleID=${bamf}

#CLEAR SAMPLE ID
curlget ${URL}/sparsely/clear/Coin/TB?sampleID=${vcf}


#POST BIGGER
bamf=ERR10225529.bam
samtools view -H ${bamf} > header;  samtools view ${bamf} | split - prefix -l 1000 --filter="cat header - | samtools view -b - | curl -k -X 'POST' --user ${USER}:${PASS} -F upload=@- ${URL}/sparsely/bam/Coin/TB?sampleID=${bamf}"


