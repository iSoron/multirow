#!/bin/bash
function title()
{
    COLS=`tput cols`    
    printf '\n'
    ruby -e 'print "'"$*"'".center('$COLS')'
    printf '\n%*s\n' "$COLS" '' | tr ' ' -    
}

RUN=../../build/lifting/benchmark/lifting-benchmark.run
make -C ../../build lifting-benchmark.run || exit 1

if [ ! -f $RUN ]; then
    echo "not found: $RUN"
    echo "please build the project before running this script"
    exit 1
fi

INSTANCES="instances/filtered/all.txt"
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

title Heuristic
$RUN $COMMON_OPTS --samples 100000 --heuristic --log $DIR/heur.log --stats $DIR/heur.yaml || exit

title Bound Original
$RUN $COMMON_OPTS --samples 1000 --bound --log $DIR/bound.log --stats $DIR/bound.yaml || exit

title Bound Pre-processing
$RUN $COMMON_OPTS --samples 1000 --bound --preprocess --log $DIR/bound-pre.log --stats $DIR/bound-pre.yaml || exit

title Naive Bounding-Box
$RUN $COMMON_OPTS --samples 100 --naive --log $DIR/naive-bbox.log --stats $DIR/naive-bbox.yaml || exit

title Naive Fixed-M
M=50
$RUN $COMMON_OPTS --samples 100 --naive --fixed-bounds $M --log $DIR/naive-fixed-$M.log --stats $DIR/naive-fixed-$M.yaml || exit

title Naive Bounding-Box Pre-processing
$RUN $COMMON_OPTS --samples 100 --naive --preprocess --log $DIR/naive-bbox-pre.log --stats $DIR/naive-bbox-pre.yaml || exit

title MIP
$RUN $COMMON_OPTS --samples 10 --mip --log $DIR/mip.log --stats $DIR/mip.yaml || exit

title MIP Pre-processing
$RUN $COMMON_OPTS --samples 10 --mip --preprocess --log $DIR/mip-pre.log --stats $DIR/mip-pre.yaml || exit

# SHEAR
# ------------------------------------------------------------------------------

ANSWERS=answers/shear-$SEED.txt
COMMON_OPTS="--shear --seed $SEED --sets $INSTANCES"

title Generating answers
$RUN $COMMON_OPTS --bound --write-answers $ANSWERS || exit
COMMON_OPTS="$COMMON_OPTS --check-answers $ANSWERS"

DIR=shear-100
mkdir -p $DIR; rm -f $DIR/*log $DIR/*yaml

title Heuristic + Shear
$RUN $COMMON_OPTS --samples 100000 --heuristic --log $DIR/heur.log --stats $DIR/heur.yaml || exit

title Bound Original + Shear
$RUN $COMMON_OPTS --samples 1000 --bound --log $DIR/bound.log --stats $DIR/bound.yaml || exit

title Bound Pre-processing + Shear
$RUN $COMMON_OPTS --samples 1000 --bound --preprocess --log $DIR/bound-pre.log --stats $DIR/bound-pre.yaml || exit

title Naive Bounding-Box + Shear
$RUN $COMMON_OPTS --samples 10 --naive --log $DIR/naive-bbox.log --stats $DIR/naive-bbox.yaml || exit

title Naive Fixed-M + Shear
M=50
$RUN $COMMON_OPTS --samples 100 --naive --fixed-bounds $M --log $DIR/naive-fixed-$M.log --stats $DIR/naive-fixed-$M.yaml || exit

title MIP + Shear
$RUN $COMMON_OPTS --samples 10 --mip --log $DIR/mip.log --stats $DIR/mip.yaml || exit

title MIP Pre-processing + Shear
$RUN $COMMON_OPTS --samples 10 --mip --preprocess --log $DIR/mip-pre.log --stats $DIR/mip-pre.yaml || exit


# TABLES
# ------------------------------------------------------------------------------

for TABLE in orig-100 shear-100; do
    echo Writing file tables/$TABLE.csv
    scripts/table.rb $TABLE/*yaml > tables/$TABLE.csv || exit
done

