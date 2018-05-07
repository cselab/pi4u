#! /bin/bash

BASE_DIR="runs"

DATA_DIR="./data"
DATA_FILE_PREFIX=
DATA_FILE_POSTFIX=".dat"


THETA_DIR="./data/theta"
THETA_FILE_PREFIX="theta_"
THETA_FILE_POSTFIX=".txt"
EV_FILE_PREFIX="evidence_"
EV_FILE_POSTFIX=".txt"

PSI_DIR="./data/psi"
PSI_FILE="psi.txt"


EXEC_FILE="../../../build/sample_posterior_theta_fast"

PAR_FILE="./tmcmc_theta.par"
PRIOR_THETA_FILE="./priors_theta.par"
PRIOR_FILE="./priors_aux_3.par"

DB_FILE="./db_theta.par"



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


	DATA_FILE="${DATA_DIR}/${DATA_FILE_PREFIX}${II}${DATA_FILE_POSTFIX}"
	THETA_FILE="${THETA_DIR}/${THETA_FILE_PREFIX}${II}${THETA_FILE_POSTFIX}"
	EV_FILE="${THETA_DIR}/${EV_FILE_PREFIX}${II}${EV_FILE_POSTFIX}"
	PSI_FILE_LOC="${PSI_DIR}/${PSI_FILE}"


	cp ${DATA_FILE}             "$RUN_DIR/data.txt"
	cp ${THETA_FILE}            "$RUN_DIR/theta.txt"
	cp ${EV_FILE}               "$RUN_DIR/evidence.txt"
	cp ${PSI_FILE_LOC}          "$RUN_DIR/psi.txt"

	cp ${EXEC_FILE}             ${RUN_DIR}
	cp ${PAR_FILE}              ${RUN_DIR}/tmcmc.par
	cp ${PRIOR_THETA_FILE}      ${RUN_DIR}/priors_theta.par
	cp ${PRIOR_FILE}            ${RUN_DIR}/priors.par
	cp ${DB_FILE}               ${RUN_DIR}/db_theta.par 


    echo "------------>  $II"


done
