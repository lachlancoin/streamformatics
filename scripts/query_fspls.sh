#USAGE
#sh query.sh  example.json


while read line; do export $line; done <~/.sfx
user=$(whoami)

WWW="fspls/data/train/Coin/LCAH"

>&2 echo "URL ${URL}"
>&2 echo "PASS ${PASS}"
>&2 echo "USER ${user}"
>&2 echo "WWW ${WWW}"



url2=${URL}/${WWW}
>&2 echo "url2 ${url2}"
flags='{"include":["1) Assoc w R","2) Assoc w R - Interim"]}'
flags='{"beam":1,"pthresh":1e-10,"return":"eval","train":"Discovery","max":2,"topn":1}'
phens='["disease_class.DBDV","gn_od_00.gn","gn_od_24.gn"]'
#flags='{"cov":20,"norm":true}'
>&2 echo $flags | jq


curl -X 'GET' -k  --user "${user}:${PASS}" -F flags="${flags}" -F phens="${phens}"  $url2
cmd="curl -X 'GET' -k  --user \"${user}:${PASS}\" -F flags='${flags}' -F phens='${phens}'  $url2"
echo $cmd

