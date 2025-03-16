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
alias curlget='curl -k -X 'GET' --user ${USER}:${PASS}'
alias curlpost='curl -k -X 'POST' --user ${USER}:${PASS}'



##TEST THE ENDPONT
curlget ${URL}/admin/echo?msg=test

##CURL EXAMPLES

##REGISTER USER FOR ORG
##NEED THE KEY

curlget ${URL}/admin/register/Coin?key=${API_KEY}

##MAKE NEW DATABASE
curlget ${URL}/admin/make/Coin/RAPIDS5/Discovery

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
samtools view -H ${bamf} > header; curlget ${URL}/who/snps/Coin/WHO/2023 | samtools view ${bamf} -L - | split - prefix -l 1000 --filter="cat header - | samtools view -b - | curl -k -X 'POST' --user ${USER}:${PASS} -F upload=@- ${URL}/who/resistance/bam/Coin/WHO/2023?sampleID=${bamf}"

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


##FSPLS ENDPOINTS
##OPTIONAL -PRELOAD DATABASE
curlpost ${URL}/fspls/data/load/Coin/golub
curlpost ${URL}/fspls/data/pheno/Coin/golub

#flags='{"pthresh":0.05,"topn":10,"beam":1,"train":"golub_data","test":"golub_data"}'
#phens='{"x.x"}'
curlpost -F phens='["x.x"]' -F flags='{"pthresh":0.05,"topn":10,"beam":1,"train":"golub_data","test":"golub_data"}' ${URL}/fspls/data/train/Coin/golub
curlpost -F phens='["x.x"]' -F flags='{"pthresh":0.05,"topn":10,"beam":1,"train":"golub_data","test":"golub_data"}' ${URL}/fspls/data/train/Coin/golub > vars.json
curlpost -F phens='["x.x"]' -F flags='{"pthresh":0.05,"topn":10,"beam":1,"train":"golub_data","test":"golub_data"}' -F upload=@vars.json ${URL}/fspls/data/makeModels/Coin/golub > models.json
curlpost -F phens='["x.x"]' -F flags='{"pthresh":0.05,"topn":10,"beam":1,"train":"golub_data","test":"golub_data"}' -F upload=@models.json ${URL}/fspls/data/evaluate/Coin/golub > eval.csv

curlpost  -F upload=@eval.csv ${URL}/fspls/data/plot/Coin > eval.png



###UPLOADING TRANSCRIPTOMICS ENDPONTS
##stream kallisto file
file="/home/unimelb.edu.au/lcoin/Data/sAPI/Coin/LCAH/6202_12_S49_R1_001.tsv.gz"
curlpost -F flags='{"format":"kallisto"}' -F upload=@${file}  ${URL}/dist/stream/Coin/RAPIDS4/Discovery


pheno_f = "/home/unimelb.edu.au/lcoin/Data/sAPI/Coin/LCAH/Discovery_meta.csv"
curlpost -F flags='{"slug_phen":"function(x)x","slug_sample":"function(x)toupper(x)","id":"Sequencing_Sample_ID","sep":","}' -F upload=@${pheno_f}  ${URL}/dist/upload_pheno/Coin/RAPIDS4/Discovery


curlget ${URL}/dist/samples/Coin/RAPIDS4/Discovery
