#! /bin/bash

BASE_DIR="runs"

CURRENT_DIR=$(pwd)

SAVE_DIR="${CURRENT_DIR}/data/theta"
mkdir -p $SAVE_DIR


export TORC_WORKERS=4


# Ilist=(1)
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

        mpirun -np 1	./sample_theta_fast
       	#./sample_theta_fast

		FNAME="${SAVE_DIR}/theta_${II}.txt"
		cp final.txt $FNAME

		FNAME="${SAVE_DIR}/evidence_${II}.txt"
		cp log_evidence.txt $FNAME

        cd $CURRENT_DIR

done
