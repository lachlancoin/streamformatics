##THIS GETS ABX PROFILE AND LINEAGE FOR A VCF
##USAGE
#sh query_vcf.sh  $PATH_TO_VCF
if [ ! -f ~/.sparsely ];then
echo "need sparsely file with credentials"
fi 

while read line; do export $line; done <~/.sparsely
user=$(whoami)

if [ ! $1 ] ; then 
echo "need vcf"
exit
fi

WWW="who/resistance/vcf/Coin/WHO/2023"

>&2 echo "URL ${URL}"
>&2 echo "PASS ${PASS}"
>&2 echo "USER ${user}"
>&2 echo "WWW ${WWW}"



url2=${URL}/${WWW}
>&2 echo "url2 ${url2}"
#curl -F upload=@ERR036186.masked.vcf.gz  $url | jq
#--data '{\"flags\":{\"cov\" : 5, \"norm\" : TRUE}}' 
#cmd="curl -k --user \"${user}:${PASS}\" -F upload=@${1} -F upload=@${2}  -F flags='{\"cov\" : 5, \"norm\" : true\" }' $url2 "
#xy
flags='{"include":["1) Assoc w R","2) Assoc w R - Interim"]}'
flags1='"exclude":["3) Uncertain significance","4) Not assoc w R - Interim","5) Not assoc w R"]} '
#flags='{"cov":20,"norm":true}'

>&2 echo $flags | jq
cmd="curl -X 'POST'  --user \"${user}:${PASS}\" -F upload=@${1} -F flags=${flags} $url2"
>&2 echo $cmd
curl -X 'POST' -k  --user "${user}:${PASS}" -F upload=@${1} -F flags="${flags}"  $url2  


#jq -r '(map(keys) | add | unique) as $cols | map(. as $row | $cols | map($row[.])) as $rows | $cols, $rows[] | @csv'


