export CLUS_DIR=$(pwd)
export MACHINE_NAME="slurmsimcont"
export RUN_NAME="test2"
export dtstart=59
export replica=1

slurmsim -v run_sim  -d \
	-e ${CLUS_DIR}/etc \
	-a ${CLUS_DIR}/etc/sacctmgr.script \
	-w ${CLUS_DIR}/workload/jobs26607.events \
	-r ${CLUS_DIR}/results/${MACHINE_NAME}/${RUN_NAME}/dtstart_${dtstart}_${replica} \
	-dtstart $dtstart
