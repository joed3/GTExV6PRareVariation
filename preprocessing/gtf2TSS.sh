#!/bin/bash

# script that takes a gtf file and extracts the TSS information
# outputs to stdout the following information (bed format, i.e. 0-based and excludes end position):
# chr start end transcripy_id(ENSG/T)

# NOTE: This script finds the TSS information for transcripts not genes.

# The script adds chr if it is missing
# if the file lacks strand information, the script prints a message and exits

set -o nounset -o errexit -o pipefail

if [ $# -ne 1 ]; then
   echo "usage: gtf2TSS input.gtf"
   exit
fi

gtf=$1

# skip header
awk 'BEGIN{OFS="\t"}{
if (substr($1,1,1)=="#") {
  next    
}
if ($3!="transcript") {
  next
}

if (substr($1,1,1)!="c") {
  chr="chr"$1
} else {
  chr=$1
}

if ($7=="+") {
  start=$4-1
  end=$4
} else if ($7=="-") {
  start=$5-1
  end=$5
} else {
  print "no strand information: exiting"
  exit 1
}

name=substr($12,2,length($12)-3)

print chr,start,end,name

}' $gtf