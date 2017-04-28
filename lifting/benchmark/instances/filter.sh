#!/bin/bash
if [ "$#" -ne 1 ]; then
    echo "usage: $0 input/input.txt.bz2"
    exit 1
fi

IN=$1
OUTDIR=filtered/

bzcat $IN | sed -f filter.sed | sort | uniq | shuf > $OUTDIR/all.txt
cat $OUTDIR/all.txt | awk '$3==4 {print}' > $OUTDIR/quads.txt
cat $OUTDIR/all.txt | awk '$3==3 {print}' > $OUTDIR/triangles.txt

