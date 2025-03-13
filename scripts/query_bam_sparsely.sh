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



>&2 echo "URL ${URL}"
>&2 echo "PASS ${PASS}"
>&2 echo "USER ${user}"



url2="${URL}/sparsely/clear/Coin/TB?sampleID=${bf}"
url3="${URL}/sparsely/bam/Coin/TB?sampleID=${bf}"
url4="${URL}/sparsely/resistance/Coin/TB?sampleID=${bf}"



#flags1='{\"sampleID\":\"'$bf'\"}'
#flags2='{"sampleID":"'$bf'"}'
#flags1=$(echo $flags2 | sed 's/"/\\"/g')

#>&2 echo $flags | jq
>&2 echo "simple command without splitting not executed " 
cmd="curl -X 'POST'  -k  --user \"${user}:${PASS}\"  -F flags='${flags}' ${url2}  | samtools view $1 -b -L - | curl -X 'POST' -k  --user \"${user}:${PASS}\" -F upload=@-   $url3 "
>&2 echo $cmd
#curl -k  --user "${user}:${PASS}" -F flags="${flags}"  $url2  | samtools view $1 -b -L - | curl -k  --user "${user}:${PASS}" -F upload=@- -F flags="${flags}"  $url3

##first clear optional
url2="${URL}/sparsely/clear/Coin/TB"
cmd="curl -X 'POST' -k  --user \"${user}:${PASS}\"   $url2"
>&2 echo $cmd
curl -X 'POST' -k  --user "${user}:${PASS}"   $url2


#curl -k  --user "${user}:${PASS}" -F flags="${flags}"  $url2  | samtools view $1  -L -  | split - prefix -l ${split} --filter='cat header - | samtools view -b - > $FILE.bam'
cmd="samtools view -H $bf > header; samtools view $bf | split - prefix -l ${split} --filter=\"cat header - | samtools view -b -  | curl -X 'POST' -k  --user \"${user}:${PASS}\"  -F flags='${flags1}' -F upload=@- ${url3}\""
>&2 echo $cmd
samtools view -H $bf > header;  samtools view $bf  | split - prefix -l ${split} --filter="cat header - | samtools view -b -  | curl -X 'POST' -k  --user \"${user}:${PASS}\"   -F upload=@- ${url3}"


cmd="curl -X 'GET' -k  --user \"${user}:${PASS}\" -F flags='${flags2}'  $url4"
echo $cmd

curl -X 'GET' -k  --user "${user}:${PASS}"  $url4


