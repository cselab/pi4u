#! /bin/bash

BASE_DIR="runs"
DATA_DIR="./data"
DATA_FILE_PREFIX=

EXEC_FILE="../source/sample_theta"
PAR_FILE="./theta.par"
MODEL_FILE="model/my_model.py"
MODEL_SCRIPT="model/doall.sh"
LL_FILE="model/log_like.py"



Ilist=(1)



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

	MODEL_DIR="$RUN_DIR/model"
	mkdir $MODEL_DIR

	DATA_FILE="${DATA_DIR}/${DATA_PREFIX}${II}.dat"



	cp ${DATA_FILE} "${MODEL_DIR}/data.txt"

	cp ${EXEC_FILE}    ${RUN_DIR}
	cp ${PAR_FILE}     ${RUN_DIR}/tmcmc.par
	cp ${MODEL_FILE}   ${MODEL_DIR}
	cp ${MODEL_SCRIPT} ${MODEL_DIR}
	cp ${LL_FILE}      ${MODEL_DIR}



done
