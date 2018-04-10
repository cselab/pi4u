#! /bin/bash

BASE_DIR="runs"
DATA_DIR="./data"
DATA_FILE_PREFIX=

EXEC_FILE="../source/sample_theta_fast"
PAR_FILE="./theta.par"

Ilist=(1)
#Ilist=(1 2 3 4 5)



if [ ! -d "$BASE_DIR" ]; then
	mkdir $BASE_DIR
fi


for I in "${Ilist[@]}"
do

	II="$(printf "%03d" ${I})"

	RUN_DIR="${BASE_DIR}/run_${II}"

	if [ -d "$RUN_DIR" ]; then
		#echo "The directory '${RUN_DIR}' exists. Delete the directory and rerun the script."
		rm -rf $RUN_DIR
	fi

	mkdir $RUN_DIR


	DATA_FILE="${DATA_DIR}/${DATA_PREFIX}${II}.dat"



	cp ${DATA_FILE} "$RUN_DIR/data.txt"

	cp ${EXEC_FILE}    ${RUN_DIR}
	cp ${PAR_FILE}     ${RUN_DIR}/tmcmc.par

done
