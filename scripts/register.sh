#!/bin/bash

#USAGE
##sh register.sh
#THIS REGISTERS THE USER
#set env vars
while read line; do export $line; done < ~/.sparsely
user=$(whoami)

export WWW="admin/register/Coin"
#export WWW="/dist/register"

echo "URL ${URL}"
echo "PASS ${PASS}"
echo "API_KEY ${API_KEY}"
echo "USER ${user}"
echo "WWW ${WWW}"


cmd="curl -k --user ${user}:${PASS}  ${URL}/${WWW}?key=${API_KEY}"
echo $cmd

ab=$(curl -k --user ${user}:${PASS}  ${URL}/${WWW}?key=${API_KEY} | jq)
echo $ab

##NOW EXPORT THE PASSWORD
#THE FOLLOWING DOES NOT WORK FOR SOME REASON

#export SPASS=$pass
echo $pass
