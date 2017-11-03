#! /bin/bash

BASE_DIR="runs"

EXEC_FILE="../source/sample_psi"
PAR_FILE="./psi.par"



RUN_DIR="${BASE_DIR}/psi"

if [ ! -d "$RUN_DIR" ]; then
	mkdir $RUN_DIR
fi


cp ${EXEC_FILE}    ${RUN_DIR}

cp ${PAR_FILE}     ${RUN_DIR}/tmcmc.par



