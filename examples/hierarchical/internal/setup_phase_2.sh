#! /bin/bash

BASE_DIR="runs"

EXEC_FILE="../../../build/sample_psi"
PAR_FILE="./tmcmc_psi.par"
PRIOR_FILE="./priors_psi.par"
PRIOR_THETA_FILE="./priors_theta.par"

DB_FILE="./db_psi.par"

RUN_DIR="${BASE_DIR}/psi"


if [ ! -d "$RUN_DIR" ]; then
	mkdir -p $RUN_DIR
fi


cp ${EXEC_FILE}             ${RUN_DIR}
cp ${PAR_FILE}              ${RUN_DIR}/tmcmc.par
cp ${PRIOR_FILE}            ${RUN_DIR}/priors.par
cp ${PRIOR_THETA_FILE}      ${RUN_DIR}/priors_theta.par

cp ${DB_FILE}               ${RUN_DIR}/db_psi.par

