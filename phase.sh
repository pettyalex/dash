#!/bin/sh

# BEAGLE Heap Size
HEAP=6500

# --- SOFTWARE SOURCES --- #

# Full path to the pipeline
PIPE=`dirname $0`
# Full path to temporary directory (optional)
TMPDIR=''
# Self-setting files from pipeline
BEAGLE="$PIPE/bin/beagle.jar"
PED_TO_BGL="$PIPE/bin/ped_to_bgl"
BGL_TO_PED="$PIPE/bin/bgl_to_ped"
SEARCH="$PIPE/bin/search"

# --- SOFTWARE SOURCES --- #

if [ $# -lt 3 ]; then
  echo "---"
  echo " BEAGLE Phasing script"
  echo " $PIPE phase.sh"
  echo "---"
  echo -e "Usage:\t\t$0 [ped file] [map file] [output name]"
  echo -e "Output:\t\tGenerates [output name].phased ped/map files for input to GERMLINE"
  exit
fi

PED=$1
MAP=$2
OUT=$3

for f in $BEAGLE $PED_TO_BGL $BGL_TO_PED $SEARCH $PED $MAP; do
  if [ ! -f $f ]; then
    echo "The file $f does not exist, please check that it is referenced properly"
    exit
  fi
done

echo "Converting to BEAGLE"
$PED_TO_BGL \
$PED \
$MAP \
> $OUT.pre_phase.bgl

nr_ind=`cat $PED | wc -l`
if [ $nr_ind -lt 250 ]; then
  samples=20;
elif [ $nr_ind -lt 500 ]; then
  samples=10;
elif [ $nr_ind -lt 1000 ]; then
  samples=4;
elif [ $nr_ind -lt 2000 ]; then
  samples=2;
else
  samples=1;
fi

echo "Run BEAGLE"
if [ -d $TMPDIR ]; then
 JAVA="java -Djava.io.tmpdir=$TMPDIR -Xmx${HEAP}m"
else
 JAVA="java -Xmx${HEAP}m"
fi

$JAVA -jar $BEAGLE \
unphased=$OUT.pre_phase.bgl log=$OUT.bgl \
niterations=10 nsamples=$samples missing=0

tail -n+2 $OUT.pre_phase.bgl.phased > $OUT.phased.bgl
rm $OUT.pre_phase.bgl.phased $OUT.pre_phase.bgl
if [ -s $OUT.pre_phase.bgl.gprobs ]; then
  cut -f 1 -d ' ' $OUT.pre_phase.bgl.gprobs | $SEARCH $MAP 2 > $OUT.phased.map
else
  cp $MAP $OUT.phased.map
fi

awk '{ print $1,$2,$3,$4,$5,$6; }' $PED > $OUT.fam
$BGL_TO_PED \
$OUT.phased.bgl \
$OUT.fam 0 \
> $OUT.phased.ped

rm -f $OUT.bgl.log $OUT.pre_phase.bgl.r2 $OUT.pre_phase.bgl.gprobs
rm -f $OUT.fam $OUT.phased.bgl
