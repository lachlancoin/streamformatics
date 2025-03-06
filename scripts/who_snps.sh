##THIS GETS ABX PROFILE AND LINEAGE FOR A VCF
##USAGE
#sh query_vcf.sh  $PATH_TO_VCF
if [ ! -f ~/.sparsely ];then
echo "need sparsely file with credentials"
fi 

while read line; do export $line; done <~/.sparsely
user=$(whoami)


WWW="who/TB/snps"

>&2 echo "URL ${URL}"
>&2 echo "PASS ${PASS}"
>&2 echo "USER ${user}"
>&2 echo "WWW ${WWW}"



url2=${URL}/${WWW}
>&2 echo "url2 ${url2}"
flags='{"chrom":"NC"}'

>&2 echo $flags | jq
cmd="curl  --user \"${user}:${PASS}\"  -F flags=${flags} $url2"
>&2 echo $cmd
curl -k  --user "${user}:${PASS}" -F flags="${flags}"  $url2  


#jq -r '(map(keys) | add | unique) as $cols | map(. as $row | $cols | map($row[.])) as $rows | $cols, $rows[] | @csv'


