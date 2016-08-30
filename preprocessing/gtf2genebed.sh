#!/bin/bash

# script that takes a gtf file and extracts gene start/end position
# outputs as bed format with gene id in last column

# Prints a header

# can provide optional argument -t if you want the transcript start/end positions instead with transript ids

set -o nounset -o errexit -o pipefail

# checking for valid command line input
usage="usage: gtf2TSS [-t] input.gtf\nOptions:\n\t-t for transcript information instead of gene information"
if [ $# -gt 2 ] || [ $# -lt 1 ]; then
   echo -e $usage
   exit
fi

if [ $# -eq 1 ]; then
    gtf=$1
    kind="gene"
elif [ "$1" == "-t" ]; then
    gtf=$2
    kind="transcript"
else
    echo -e $usage
    exit
fi

# skip header
awk -v type=$kind 'BEGIN{OFS="\t"}{
if (substr($1,1,1)=="#") {
  next    
}
if ($3!=type) {
  next  
}

if (type=="gene"){
  name=substr($10,2,length($10)-3) 
} else {
  name=substr($12,2,length($12)-3)
}
print $1,$4-1,$5,name

}' $gtf
