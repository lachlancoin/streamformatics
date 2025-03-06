##USAGE
#sh upload_pheno.sh  $PATH_TO_CSV

while read line; do export $line; done </.sparsely
user=$(whoami)

if [ ! $1 ] ; then 
echo "need pheno"
exit
fi


export pheno=$1

#path="mtb"
export WWW="/Coin/TB/dist/upload_pheno"

echo "URL ${URL}"
echo "PASS ${PASS}"
echo "USER ${user}"
echo "WWW ${WWW}"



url2=${URL}/${WWW}
echo "url2 ${url2}"

curl -k --user "${USER}:${PASS}" -F upload=@${1}  $url2 | jq


