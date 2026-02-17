#!/bin/bash
#BSUB -J pleio_a1[1-3]                 # sim = 1..3
#BSUB -q long
#BSUB -n 4
#BSUB -M 64
#BSUB -R rusage[mem=64]
#BSUB -W 25:00
#BSUB -o /rsrch8/home/biostatistics/swang25/other_summary_statistics/CRC_IBD/errDir/pleio_a1_%I.out
#BSUB -e /rsrch8/home/biostatistics/swang25/other_summary_statistics/CRC_IBD/errDir/pleio_a1_%I.err
#BSUB -cwd /rsrch8/home/biostatistics/swang25/other_summary_statistics/CRC_IBD

module load R/4.4.1
export R_LIBS_USER="/home/swang25/R/ubuntu/4.4.1"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

AID=1
SIM=${LSB_JOBINDEX}
SEED=$((1000 + 100*AID + SIM))
echo "[$(date)] host=$(hostname) AID=$AID SIM=$SIM SEED=$SEED" 1>&2

/risapps/rhel8/R/4.4.1/lib64/R/bin/Rscript \
  /rsrch8/home/biostatistics/swang25/other_summary_statistics/CRC_IBD/sample_approach_clean.R \
  "$AID" "$SIM" "$SEED"
