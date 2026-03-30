#!/bin/bash

# === Konfiguration ===
SERVER="jab49wd@julia2.hpc.uni-wuerzburg.de"         
LOCAL_DIR="C:/Forschung/SEHTMT/HTMT3/claude_fix_all/"  # WSL-Pfad
REMOTE_DIR="/home/jab49wd/revision1sim"             
SLURM_SCRIPT_PATH="C:/Forschung/SEHTMT/HTMT3/shfiles/HTMTsim_all.sh"
SLURM_SCRIPT="HTMTsim_all.sh"

# === 1. Ordner zum HPC schieben ===
echo "sending files to hpc."
scp -r "$LOCAL_DIR"* "$SERVER:$REMOTE_DIR/"
if [ $? -ne 0 ]; then
    echo "upload failed."
    exit 1
fi
echo "upload succeded."

echo "sending sh file to hpc."
scp -r "$SLURM_SCRIPT_PATH" "$SERVER:$REMOTE_DIR/"
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

