#USAGE
#sh query.sh  example.json


#set env vars
while read line; do export $line; done < ~/.sfx
user=$(whoami)

if [ ! $1 ] ; then 
echo "need json"
exit
fi


export JSON=$1

#path="all"

export WWW="sparsely/score"

echo "URL ${URL}"
echo "PASS ${PASS}"
echo "USER ${user}"
echo "WWW ${WWW}"


ab="curl -k --user \"${USER}:${PASS}\" -H \"Content-Type: application/json\" -d @${1} ${URL}/${WWW})"
echo $ab
#exit ;

d1=$(date +%s)
ab=$(curl -k --user "${USER}:${PASS}" -H "Content-Type: application/json" -d @${1} ${URL}/${WWW})
d2=$(date +%s)
d3=$(($d2-$d1))

echo $ab | jq
echo $d3

