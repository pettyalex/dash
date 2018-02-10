#!/bin/sh

# --- SOFTWARE SOURCES --- #
PIPE=`dirname $0`
DASH="$PIPE/bin/dash_cc"
PLINK="$PIPE/bin/plink"
BPARSE="$PIPE/bin/parse_bmatch"
# --- SOFTWARE SOURCES --- #

if [ $# -lt 3 ]; then
  echo "---"
  echo " DASH-CC Haplotype association script"
  echo " $PIPE"
  echo "---"
  echo -e "Usage:\t$0\t[bmatch prefix] [fam file] [output name]"
  echo -e "Output:\tGenerates [output] bim/bed/fam file with coded haploypes, and [output] cmap file with haplotype boundaries"
  exit
fi

MATCH=$1
FAM=$2
OUT="$3.dash_cc"

for f in $DASH $PLINK $BPARSE $MATCH.bmatch $MATCH.bsid $MATCH.bmid $FAM; do
  if [ ! -f $f ]; then
    echo "The file $f does not exist, please check that it is referenced properly"
    exit
  fi
done

CHR=`$BPARSE $MATCH.bmatch $MATCH.bsid $MATCH.bmid | head -n1 | awk '{ print $5 }'`

echo "Using chromosome $CHR for DASH output"
$BPARSE $MATCH.bmatch $MATCH.bsid $MATCH.bmid | cut -f 1,2,4 | $DASH $FAM $OUT

echo "DASH run completed, converting to binary PLINK format"
cut -f 1-5 $OUT.clst | awk -v c=$CHR '{ print c,c"_"$1,$2,$3,$4,$5 }' > $OUT.cmap
awk -v c=$CHR '{ print c,$2,0,int(($3+$4)/2) }' $OUT.cmap > $OUT.map

$PLINK --maf 0.001 --file $OUT --make-bed --recode --out $OUT

rm -f $OUT.ped $OUT.map $OUT.log
