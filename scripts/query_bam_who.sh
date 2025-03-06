#conda activate japsa

if [ ! -f ~/.sfx ];then
echo "need sparsely file with credentials"
fi 


if [ ! $1 ]; then
echo ' need bam'
exit 1
fi

bf=$1
split=$2
if [ !$2 ]; then
split=4000
fi

while read line; do export $line; done <~/.sfx
user=$(whoami)


WWW="who/snps/Coin/WHO"

>&2 echo "URL ${URL}"
>&2 echo "PASS ${PASS}"
>&2 echo "USER ${user}"
>&2 echo "WWW ${WWW}"



url2=${URL}/${WWW}
WWW="who/resistance/bam/Coin/WHO/2023"
url3=${URL}/${WWW}
WWW="who/clear/Coin/WHO/2023"
url4=${URL}/${WWW}
WWW="who/resistance/sample/Coin/WHO/2023"
url5=${URL}/${WWW}



>&2 echo "url2 ${url2}"
flags='{}'
flags1='{\"sampleID\":\"'$bf'\"}'
flags2='{"sampleID":"'$bf'"}'

>&2 echo $flags | jq
#cmd="curl -k  --user \"${user}:${PASS}\"  -F flags='${flags}' ${url2}  | samtools view $1 -b -L - | curl -k  --user \"${user}:${PASS}\" -F upload=@-   $url3 | jq"
#>&2 echo $cmd
cmd="curl -X 'GET' -k  --user \"${user}:${PASS}\" -F flags='${flags2}'  $url4"
>&2 echo $cmd
curl -X 'GET' -k  --user "${user}:${PASS}" -F flags="${flags2}"  $url4
#curl -k  --user "${user}:${PASS}" -F flags="${flags}"  $url2  | samtools view $1 -b -L - | curl -k  --user "${user}:${PASS}" -F upload=@- -F flags="${flags}"  $url3
#curl -k  --user "${user}:${PASS}" -F flags="${flags}"  $url2  | samtools view $1  -L -  | split - prefix -l ${split} --filter='cat header - | samtools view -b - > $FILE.bam'
cmd="samtools view -H $bf > header; curl -X 'GET' -k  --user \"${user}:${PASS}\" -F flags=\"${flags}\"  $url2  | tail -n +2 | samtools view $bf  -L -  | split - prefix -l ${split} --filter=\"cat header - | samtools view -b -  | curl -X 'POST' -k  --user \"${user}:${PASS}\"  -F flags='${flags1}' -F upload=@- ${url3}\""
>&2 echo $cmd

#curl -k  --user "${user}:${PASS}" -F flags="${flags}"  $url2  | samtools view $bf  -L  > tmp.1.bam

#flags1='{}'
samtools view -H $bf > header; curl -X 'GET' -k  --user "${user}:${PASS}" -F flags="${flags}"  $url2 | tail -n +2 | samtools view $bf  -L -  | split - prefix -l ${split} --filter="cat header - | samtools view -b -  | curl -X 'POST' -k  --user \"${user}:${PASS}\"  -F flags=${flags1} -F upload=@- ${url3}"


cmd="curl -k  --user \"${user}:${PASS}\" -F flags='${flags2}'  $url5"
echo $cmd

curl -X 'GET' -k  --user "${user}:${PASS}" -F flags="${flags2}"  $url5


