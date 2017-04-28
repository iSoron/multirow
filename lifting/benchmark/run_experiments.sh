#!/bin/bash
function title()
{
    COLS=`tput cols`    
    printf '\n'
    ruby -e 'print "'"$*"'".center('$COLS')'
    printf '\n%*s\n' "$COLS" '' | tr ' ' -    
}

RUN=../../build/lifting/benchmark/lifting-benchmark.run
if [ ! -f $RUN ]; then
    echo "not found: $RUN"
    echo "please build the project before running this script"
    exit 1
fi

INSTANCES="instances/filtered/all.txt"
# SAMPLES_SLOW=10
# SAMPLES_MEDIUM=100
# SAMPLES_FAST=1000
SAMPLES_SLOW=1
SAMPLES_MEDIUM=1
SAMPLES_FAST=1
SEED=1240

# ORIGINAL
# ------------------------------------------------------------------------------

ANSWERS="answers/orig-$SEED.txt"
COMMON_OPTS="--seed $SEED --sets $INSTANCES"

title Generating answers
$RUN $COMMON_OPTS --bound --write-answers $ANSWERS || exit
COMMON_OPTS="$COMMON_OPTS --check-answers $ANSWERS"

DIR=orig-100
mkdir -p $DIR; rm -f $DIR/*log $DIR/*yaml

title Bound Original
$RUN $COMMON_OPTS --samples $SAMPLES_FAST --bound --log $DIR/bound-nopre.log --stats $DIR/bound-nopre.yaml || exit

title Bound Pre-processing
$RUN $COMMON_OPTS --samples $SAMPLES_FAST --bound --preprocess --log $DIR/bound-pre.log --stats $DIR/bound-pre.yaml || exit

title Naive Bounding-Box
$RUN $COMMON_OPTS --samples $SAMPLES_MEDIUM --naive --log $DIR/naive-bbox.log --stats $DIR/naive-bbox.yaml || exit

title Naive Fixed-M
M=50
$RUN $COMMON_OPTS --samples $SAMPLES_MEDIUM --naive --fixed-bounds $M --log $DIR/naive-fixed-$M.log --stats $DIR/naive-fixed-$M.yaml || exit

title MIP
$RUN $COMMON_OPTS --samples $SAMPLES_SLOW --mip --log $DIR/mip.log --stats $DIR/mip.yaml || exit

# SHEAR
# ------------------------------------------------------------------------------

ANSWERS=answers/shear-$SEED.txt
COMMON_OPTS="--shear --seed $SEED --sets $INSTANCES"

title Generating answers
$RUN $COMMON_OPTS --bound --write-answers $ANSWERS || exit
COMMON_OPTS="$COMMON_OPTS --check-answers $ANSWERS"

DIR=shear-100
mkdir -p $DIR; rm -f $DIR/*log $DIR/*yaml

title Bound Pre-processing + Shear
$RUN $COMMON_OPTS --samples $SAMPLES_FAST --bound --preprocess --log $DIR/bound-pre.log --stats $DIR/bound-pre.yaml || exit

title Bound Original + Shear
$RUN $COMMON_OPTS --samples $SAMPLES_MEDIUM --bound --log $DIR/bound-nopre.log --stats $DIR/bound-nopre.yaml || exit

title Naive Bounding-Box + Shear
$RUN $COMMON_OPTS --samples $SAMPLES_SLOW --naive --log $DIR/naive-bbox.log --stats $DIR/naive-bbox.yaml || exit

title Naive Fixed-M + Shear
M=50
$RUN $COMMON_OPTS --samples $SAMPLES_MEDIUM --naive --fixed-bounds $M --log $DIR/naive-fixed-$M.log --stats $DIR/naive-fixed-$M.yaml || exit

title MIP
$RUN $COMMON_OPTS --samples $SAMPLES_SLOW --mip --log $DIR/mip.log --stats $DIR/mip.yaml || exit

# TABLES
# ------------------------------------------------------------------------------

for TABLE in orig-100 shear-100; do
    echo Writing file tables/$TABLE.csv
    scripts/table.rb $TABLE/*yaml > tables/$TABLE.csv || exit
done

