#!/bin/bash

# === Konfiguration ===
SERVER="jab49wd@julia2.hpc.uni-wuerzburg.de"         
LOCAL_DIR="C:/Forschung/SEHTMT/HTMT3/arraysimulation/"  # WSL-Pfad
REMOTE_DIR="/home/jab49wd/simulationhtmt/"             
SLURM_SCRIPT="submit.sh"

ssh "$SERVER" "mkdir -p $REMOTE_DIR/code $REMOTE_DIR/vines"

# === 1. Ordner zum HPC schieben ===
echo "sending files to hpc."
scp -r "$LOCAL_DIR/conditions.rds" "$SERVER:$REMOTE_DIR/code/" &&
scp -r "$LOCAL_DIR/setup.R" "$SERVER:$REMOTE_DIR/code/" &&
scp -r "$LOCAL_DIR/sim.R" "$SERVER:$REMOTE_DIR/code/" && 
scp -r "$LOCAL_DIR/vines/" "$SERVER:$REMOTE_DIR/vines/"
if [ $? -ne 0 ]; then
    echo "upload failed."
    exit 1
fi
echo "upload succeded."

echo "sending sh file to hpc."
scp -r "$LOCAL_DIR/$SLURM_SCRIPT" "$SERVER:$REMOTE_DIR/"
if [ $? -ne 0 ]; then 
	echo "upload of sh file failed."
	exit 1
fi
echo "upload sh file succeded."

# === 2. SLURM-Job submitten ===
echo "submit SLURM-job."
JOB_OUTPUT=$(ssh "$SERVER" "cd $REMOTE_DIR && sbatch $SLURM_SCRIPT")
if [ $? -ne 0 ]; then
    echo "job-submit failed."
    exit 1
fi

# Job-ID extrahieren
JOB_ID=$(echo "$JOB_OUTPUT" | grep -oP '\d+')
echo "Job submitted with ID: $JOB_ID"

read -p "Drücke Enter zum Schließen..."

