#! /bin/bash

BASE_DIR="runs"

CURRENT_DIR=$(pwd)

SAVE_DIR="${CURRENT_DIR}/data/theta"

export TORC_WORKERS=24


Ilist=(1)
# Ilist=(1 2 3 4 5)



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

        mpirun -np 1	./sample_theta

		FNAME="${SAVE_DIR}/theta_${II}.txt"
		cp final.txt $FNAME

		FNAME="${SAVE_DIR}//evidence_${II}.txt"
		cp fitness.txt $FNAME

        cd $CURRENT_DIR

done
