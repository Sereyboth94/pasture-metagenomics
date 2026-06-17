#!/bin/bash
#SBATCH --job-name=kraken2_rhizo_humanfree
#SBATCH --account=lincoln04266
#SBATCH --partition=milan
#SBATCH --qos=normal
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --array=1-48%10
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail


module load Kraken2/2.1.6-GCC-12.3.0

# -------------------------
# Paths
# -------------------------
# Human-filtered nonhost reads
IN_DIR="/home/soths/00_nesi_projects/lincoln04415_nobackup/Metagenomics/Fastq/Rhizosphere/HumanFree"

KRAKEN_DB="/home/soths/00_nesi_projects/lincoln04415_nobackup/Metagenomics/Databases/Kraken2"

OUT_BASE="/home/soths/00_nesi_projects/lincoln04415_nobackup/Metagenomics/Fastq/Rhizosphere/Kraken2_Rhizosphere_HumanFree"
REP_DIR="${OUT_BASE}/reports"
KRAK_DIR="${OUT_BASE}/kraken"
LOG_DIR="${OUT_BASE}/logs"

mkdir -p "${OUT_BASE}" "${REP_DIR}" "${KRAK_DIR}" "${LOG_DIR}" logs

# -------------------------
# Kraken2 settings
# -------------------------
# For soil/rhizosphere, 0.05 is a good starting point (reduce false positives).
# If too many reads become unclassified, try 0.02; if you want stricter, try 0.1.
CONFIDENCE="0.05"

# -------------------------
# Build R1 list (humanfree reads)
# -------------------------
mapfile -t R1_LIST < <(
  find "${IN_DIR}" -maxdepth 1 -type f \( \
    -name "*_humanfree_1.fq.gz" -o \
    -name "*_humanfree_1.fastq.gz" \
  \) | sort
)

N="${#R1_LIST[@]}"
ID="${SLURM_ARRAY_TASK_ID}"

if [[ "${N}" -eq 0 ]]; then
  echo "ERROR: No humanfree R1 files found in ${IN_DIR}"
  exit 1
fi

if [[ "${ID}" -gt "${N}" ]]; then
  echo "Task ${ID} exceeds samples detected (${N}). Exiting."
  exit 0
fi

R1="${R1_LIST[$((ID-1))]}"

# -------------------------
# Derive R2 + sample name
# -------------------------
if [[ "${R1}" == *"_humanfree_1.fq.gz" ]]; then
  R2="${R1/_humanfree_1.fq.gz/_humanfree_2.fq.gz}"
  SAMPLE="$(basename "${R1}" _humanfree_1.fq.gz)"
elif [[ "${R1}" == *"_humanfree_1.fastq.gz" ]]; then
  R2="${R1/_humanfree_1.fastq.gz/_humanfree_2.fastq.gz}"
  SAMPLE="$(basename "${R1}" _humanfree_1.fastq.gz)"
else
  echo "ERROR: Unexpected R1 naming: ${R1}"
  exit 1
fi

if [[ ! -f "${R2}" ]]; then
  echo "ERROR: Missing R2 for ${SAMPLE}"
  echo "R1=${R1}"
  echo "Expected R2=${R2}"
  exit 1
fi

REPORT="${REP_DIR}/${SAMPLE}.kreport.txt"
KRAKEN_OUT="${KRAK_DIR}/${SAMPLE}.kraken.txt"

echo "[$(date)] SAMPLE=${SAMPLE}"
echo "[$(date)] R1=${R1}"
echo "[$(date)] R2=${R2}"
echo "[$(date)] DB=${KRAKEN_DB}"
echo "[$(date)] OUT=${OUT_BASE}"
echo "[$(date)] CONFIDENCE=${CONFIDENCE}"

# -------------------------
# Run Kraken2
# -------------------------
kraken2 \
  --db "${KRAKEN_DB}" \
  --paired \
  --gzip-compressed \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --confidence "${CONFIDENCE}" \
  --report "${REPORT}" \
  --output "${KRAKEN_OUT}" \
  "${R1}" "${R2}"

echo "[$(date)] Finished ${SAMPLE}"
