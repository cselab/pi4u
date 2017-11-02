#! /bin/bash

BASE_DIR="runs"

CURRENT_DIR=$(pwd)

export TORC_WORKERS=24


Ilist=(1 2 3 4 5)



if [ ! -d "$BASE_DIR" ]; then
        echo "${BASE_DIR} does not exist"
        exit
fi


for I in "${Ilist[@]}"
do

        II="$(printf "%03d" ${I})"

        RUN_DIR="${BASE_DIR}/run_${II}"

        if [ ! -d "$RUN_DIR" ]; then
                echo "${RUN_DIR} does not exist"
                exit
        fi

        cd $RUN_DIR

        mpirun -np 1   ./sample_theta

        cd $CURRENT_DIR

done
