#! /bin/bash

BASE_DIR="runs"
DATA_DIR="./data"
DATA_FILE_PREFIX=

EXEC_FILE="../source/sample_posterior_theta_fast"
PAR_FILE="./posterior_theta.par"

THETA_DIR="./data/theta"
THETA_FILE_PREFIX='theta_'
EV_FILE_PREFIX='evidence_'

PSI_DIR="./data/psi"

# Ilist=(1)
Ilist=(1 2 3 4 5)



if [ ! -d "$BASE_DIR" ]; then
	mkdir $BASE_DIR
fi


for I in "${Ilist[@]}"
do

	II="$(printf "%03d" ${I})"

	RUN_DIR="${BASE_DIR}/posterior_run_${II}"

	if [ -d "$RUN_DIR" ]; then
		echo "The directory '${RUN_DIR}' exists. Delete the directory and rerun the script."; exit
		# rm -rf $RUN_DIR
	fi

	mkdir $RUN_DIR


	DATA_FILE="${DATA_DIR}/${DATA_FILE_PREFIX}${II}.dat"
	THETA_FILE="${THETA_DIR}/${THETA_FILE_PREFIX}${II}.txt"
	EV_FILE="${THETA_DIR}/${EV_FILE_PREFIX}${II}.txt"
	PSI_FILE="${PSI_DIR}/psi.txt"


	cp ${DATA_FILE}  "$RUN_DIR/data.txt"
	cp ${THETA_FILE} "$RUN_DIR/theta.txt"
	cp ${THETA_FILE} "$RUN_DIR/evidence.txt"
	cp ${PSI_FILE} 	 "$RUN_DIR/psi.txt"


	cp ${EXEC_FILE}    ${RUN_DIR}
	cp ${PAR_FILE}     ${RUN_DIR}/tmcmc.par


done
