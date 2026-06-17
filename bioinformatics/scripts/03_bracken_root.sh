#!/bin/bash
#SBATCH --job-name=bracken_root
#SBATCH --account=lincoln04266
#SBATCH --partition=milan
#SBATCH --qos=normal
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-48%12
#SBATCH --output=logs/bracken_root_%A_%a.out
#SBATCH --error=logs/bracken_root_%A_%a.err

set -euo pipefail

module purge
module load NeSI/zen3
module load Bracken/2.7-GCC-11.3.0

# -------------------------
# Paths
# -------------------------
KRAKEN_REPORT_DIR="/home/soths/00_nesi_projects/lincoln04415_nobackup/Metagenomics/Fastq/Root/Kraken2_Root/reports"
DB="/home/soths/00_nesi_projects/lincoln04415_nobackup/Metagenomics/Databases/Kraken2"

OUT_BASE="/home/soths/00_nesi_projects/lincoln04415_nobackup/Metagenomics/Fastq/Root/Bracken_Root"
OUT_P="${OUT_BASE}/phylum"
OUT_F="${OUT_BASE}/family"
OUT_G="${OUT_BASE}/genus"
OUT_S="${OUT_BASE}/species"   # optional but keep for completeness

mkdir -p "${OUT_BASE}" "${OUT_P}" "${OUT_F}" "${OUT_G}" "${OUT_S}" logs

# -------------------------
# Build ordered list of reports
# -------------------------
mapfile -t REP_LIST < <(find "${KRAKEN_REPORT_DIR}" -maxdepth 1 -type f -name "*.kreport.txt" | sort)

N="${#REP_LIST[@]}"
ID="${SLURM_ARRAY_TASK_ID}"

if [[ "${N}" -eq 0 ]]; then
  echo "ERROR: No kreport files found in ${KRAKEN_REPORT_DIR}"
  exit 1
fi

if [[ "${ID}" -gt "${N}" ]]; then
  echo "Task ${ID} exceeds reports detected (${N}). Exiting."
  exit 0
fi

REPORT="${REP_LIST[$((ID-1))]}"
SAMPLE="$(basename "${REPORT}" .kreport.txt)"

# -------------------------
# Bracken settings
# -------------------------
READLEN=150   # change if needed
THRESH=10     # min reads threshold

echo "[$(date)] SAMPLE=${SAMPLE}"
echo "[$(date)] REPORT=${REPORT}"
echo "[$(date)] DB=${DB}"
echo "[$(date)] READLEN=${READLEN}"
echo "[$(date)] bracken=$(which bracken)"

# -------------------------
# Run Bracken ranks
# -------------------------
# Phylum
bracken -d "${DB}" -i "${REPORT}" -o "${OUT_P}/${SAMPLE}.bracken.P.txt" -r "${READLEN}" -l P -t "${THRESH}"

# Family
bracken -d "${DB}" -i "${REPORT}" -o "${OUT_F}/${SAMPLE}.bracken.F.txt" -r "${READLEN}" -l F -t "${THRESH}"

# Genus
bracken -d "${DB}" -i "${REPORT}" -o "${OUT_G}/${SAMPLE}.bracken.G.txt" -r "${READLEN}" -l G -t "${THRESH}"

# Species (optional; comment out if you don’t want it)
bracken -d "${DB}" -i "${REPORT}" -o "${OUT_S}/${SAMPLE}.bracken.S.txt" -r "${READLEN}" -l S -t "${THRESH}"

echo "[$(date)] Done ${SAMPLE}"