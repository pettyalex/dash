#!/bin/sh

# --- SOFTWARE SOURCES --- #
PIPE=`dirname $0`
GERMLINE="$PIPE/bin/germline"
# --- SOFTWARE SOURCES --- #

if [ $# -lt 3 ]; then
  echo "---"
  echo " GERMLINE IBD detection script"
  echo " $PIPE gline.sh"
  echo "---"
  echo -e "Usage:\t$0 [phased ped file] [map file] [output name] [optional parameters]"
  echo -e "Output:\tGenerates [output] match file with IBD segments"
  exit
fi

PED=$1
MAP=$2
OUT=$3

for f in $GERMLINE $PED $MAP; do
  if [ ! -f $f ]; then
    echo "The file $f does not exist, please check that it is referenced properly"
    exit
  fi
done

echo -e "1\tRunning GERMLINE"
echo '1' > $OUT.run
echo $MAP >> $OUT.run
echo $PED >> $OUT.run
echo $OUT >> $OUT.run

$GERMLINE $4 < $OUT.run

rm $OUT.run 
