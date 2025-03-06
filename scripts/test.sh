##a simple program to test if url is working
##USAGE 
##sh  test.sh

#set env vars
while read line; do export $line; done <~/.sparsely
user=$(whoami)
export WWW="admin/echo?msg=testing"

echo "URL ${URL}"
echo "PASS ${PASS}"
echo "USER ${user}"
echo "WWW ${WWW}"


cmd="curl -k ${URL}/${WWW}"
echo $cmd

