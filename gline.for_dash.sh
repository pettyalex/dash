#!/bin/sh

PIPE=`dirname $0`
GERMLINE="$PIPE/gline.sh"

for f in $GERMLINE; do
  if [ ! -f $f ]; then
    echo "The file $f does not exist, please check that it is referenced properly"
    exit
  fi
done

$GERMLINE $1 $2 $3 "$4 -haploid -bin_out -min_m 1 -bits 32 -err_hom 1 -err_het 1"
