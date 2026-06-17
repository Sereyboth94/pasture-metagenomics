#!/bin/bash
#SBATCH --job-name=humann4_root_full
#SBATCH --account=lincoln04266
#SBATCH --partition=milan
#SBATCH --qos=normal
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G
#SBATCH --array=1-48%6
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

# -------------------------
# Clean modules + conda env
# -------------------------
module load Miniconda3/23.10.0-1

ENV="/home/soths/00_nesi_projects/lincoln04415_nobackup/Metatrans_Test/envs/biobakery4_py312_clean"
source activate "$ENV"

# Force env binaries FIRST (avoid ~/.local humann v3.9)
export PATH="$ENV/bin:$PATH"
hash -r

# -------------------------
# MetaPhlAn wrapper (HUMAnN legacy args -> MetaPhlAn4)
# also strips flags MetaPhlAn4 doesn't support
# -------------------------
WRAPDIR="$HOME/bin_humann_wrappers"
mkdir -p "$WRAPDIR"

cat > "$WRAPDIR/metaphlan" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
REAL="$HOME/00_nesi_projects/lincoln04415_nobackup/Metatrans_Test/envs/biobakery4_py312_clean/bin/metaphlan"

# HUMAnN version check: only first line (avoid DB warning noise)
if [[ "${1:-}" == "--version" ]]; then
  "$REAL" --version 2>/dev/null | head -n 1
  exit 0
fi

args=()
i=1
while [[ $i -le $# ]]; do
  a="${!i}"
  case "$a" in
    --bowtie2db)  i=$((i+1)); args+=(--db_dir "${!i}") ;;
    --bowtie2out) i=$((i+1)); args+=(--mapout "${!i}") ;;
    --add_viruses|--unclassified_estimation)
      # MetaPhlAn 4.2+ no longer supports these; drop them
      ;;
    *) args+=("$a") ;;
  esac
  i=$((i+1))
done

exec "$REAL" "${args[@]}"
EOF

chmod +x "$WRAPDIR/metaphlan"
export PATH="$WRAPDIR:$PATH"
hash -r

# -------------------------
# Paths
# -------------------------
INDIR="/home/soths/00_nesi_projects/lincoln04415_nobackup/Metagenomics/Fastq/Root/Merged_Fastq"
OUTBASE="/home/soths/00_nesi_projects/lincoln04415_nobackup/Metagenomics/Fastq/Root/Humann_Output"

DBROOT="/home/soths/00_nesi_projects/lincoln04415_nobackup/Metatrans_Test/db"
HUMANN_NUC_DB="$DBROOT/chocophlan"
HUMANN_PROT_DB="$DBROOT/uniref"
HUMANN_UTIL_DB="$DBROOT/utility_mapping"

METAPHLAN_DB="$DBROOT/mpa_vJan25_CHOCOPhlAnSGB_202503"
METAPHLAN_INDEX="mpa_vJan25_CHOCOPhlAnSGB_202503"
HUMANN_EXPECTED_DB="vOct22_CHOCOPhlAnSGB_202403"

mkdir -p "$OUTBASE" logs

# -------------------------
# Sanity: binaries
# -------------------------
echo "humann:        $(which humann)"
humann --version || true
echo "humann_config: $(which humann_config)"
echo "metaphlan:     $(which metaphlan)"
metaphlan --version || true

# -------------------------
# Ensure HUMAnN config uses FULL folders
# (writes into the env's humann cfg, not ~/.local)
# -------------------------
humann_config --update database_folders nucleotide "$HUMANN_NUC_DB"
humann_config --update database_folders protein "$HUMANN_PROT_DB"
humann_config --update database_folders utility_mapping "$HUMANN_UTIL_DB"

# -------------------------
# DB sanity
# -------------------------
[[ -d "$HUMANN_NUC_DB" ]] || { echo "ERROR: missing $HUMANN_NUC_DB"; exit 1; }
[[ -d "$HUMANN_PROT_DB" ]] || { echo "ERROR: missing $HUMANN_PROT_DB"; exit 1; }
[[ -d "$HUMANN_UTIL_DB" ]] || { echo "ERROR: missing $HUMANN_UTIL_DB"; exit 1; }
[[ -f "$HUMANN_UTIL_DB/metacyc_reactions_level4ec_only.uniref.bz2" ]] || { echo "ERROR: missing MetaCyc reactions in utility_mapping"; exit 1; }
[[ -e "$HUMANN_UTIL_DB/metacyc_pathways_structured_filtered_v24_subreactions" ]] || { echo "ERROR: missing MetaCyc pathways struct in utility_mapping"; exit 1; }

[[ -d "$METAPHLAN_DB" ]] || { echo "ERROR: missing $METAPHLAN_DB"; exit 1; }
[[ -f "$METAPHLAN_DB/${METAPHLAN_INDEX}.1.bt2l" ]] || { echo "ERROR: missing MetaPhlAn index bt2l"; exit 1; }

# -------------------------
# Sample list
# -------------------------
cd "$INDIR"
mapfile -t SAMPLES < <(ls -1 *_merged.fq.gz | sed 's/_merged\.fq\.gz$//' | sort)

N=${#SAMPLES[@]}
IDX=$((SLURM_ARRAY_TASK_ID - 1))
[[ "$IDX" -ge 0 && "$IDX" -lt "$N" ]] || { echo "ERROR: task $SLURM_ARRAY_TASK_ID out of range (N=$N)"; exit 1; }

ID="${SAMPLES[$IDX]}"
IN="${INDIR}/${ID}_merged.fq.gz"
OUTDIR="${OUTBASE}/${ID}"
THREADS="${SLURM_CPUS_PER_TASK}"

mkdir -p "$OUTDIR"

echo "[$(date)] Sample:  $ID"
echo "[$(date)] Input:   $IN"
echo "[$(date)] Outdir:  $OUTDIR"
echo "[$(date)] Threads: $THREADS"

[[ -s "$IN" ]] || { echo "ERROR: Missing $IN"; exit 1; }

# -------------------------
# Skip if final outputs look real
# -------------------------
FINAL_GF="${OUTDIR}/${ID}_merged_2_genefamilies.tsv"
FINAL_RXN="${OUTDIR}/${ID}_merged_3_reactions.tsv"
FINAL_PWY="${OUTDIR}/${ID}_merged_4_pathabundance.tsv"

if [[ -s "$FINAL_GF" && -s "$FINAL_RXN" && -s "$FINAL_PWY" ]]; then
  if [[ "$(wc -l < "$FINAL_PWY")" -gt 5 ]]; then
    echo "[$(date)] Outputs already exist and look OK. Skipping $ID."
    exit 0
  fi
  echo "[$(date)] pathabundance too small; re-running $ID."
fi

# -------------------------
# 1) MetaPhlAn profile (required: rel_ab_w_read_stats)
# -------------------------
MP_RAW="${OUTDIR}/${ID}_metaphlan_rel_ab_w_read_stats.tsv"
MP_OK="${OUTDIR}/${ID}_metaphlan_for_humann.tsv"
MP_MAPOUT="${OUTDIR}/${ID}_metaphlan_mapout.bz2"

echo "[$(date)] Running MetaPhlAn..."
metaphlan "$IN" \
  --input_type fastq \
  --bowtie2db "$METAPHLAN_DB" \
  -x "$METAPHLAN_INDEX" \
  --nproc "$THREADS" \
  --offline \
  -t rel_ab_w_read_stats \
  --bowtie2out "$MP_MAPOUT" \
  -o "$MP_RAW"

[[ -s "$MP_RAW" ]] || { echo "ERROR: MetaPhlAn output empty: $MP_RAW"; exit 1; }

# Patch header tag if HUMAnN alpha expects a specific DB label
sed "1s/${METAPHLAN_INDEX}/${HUMANN_EXPECTED_DB}/" "$MP_RAW" > "$MP_OK"
[[ -s "$MP_OK" ]] || { echo "ERROR: patched MetaPhlAn output empty: $MP_OK"; exit 1; }

# -------------------------
# 2) HUMAnN (force full utility mapping)
# -------------------------
echo "[$(date)] Running HUMAnN..."
humann \
  -i "$IN" \
  -o "$OUTDIR" \
  --threads "$THREADS" \
  --nucleotide-database "$HUMANN_NUC_DB" \
  --protein-database "$HUMANN_PROT_DB" \
  --utility-database "$HUMANN_UTIL_DB" \
  --taxonomic-profile "$MP_OK"

# -------------------------
# Validate pathabundance
# -------------------------
PWY_LINES=$(wc -l < "$FINAL_PWY" || echo 0)
echo "[$(date)] pathabundance lines: $PWY_LINES"
if [[ "$PWY_LINES" -le 5 ]]; then
  echo "ERROR: pathabundance looks wrong (<=5 lines). Check log:"
  echo "  ${OUTDIR}/${ID}_merged_0.log"
  echo "Look for utility_DEMO or pathways file paths."
  exit 2
fi

# -------------------------
# Cleanup (keep key outputs + log)
# -------------------------
echo "[$(date)] Cleaning intermediates..."
rm -f "$MP_RAW" "$MP_OK" "$MP_MAPOUT" 2>/dev/null || true

rm -f \
  "$OUTDIR"/*_bowtie2out*.bz2 \
  "$OUTDIR"/*_bowtie2out*.txt \
  "$OUTDIR"/*diamond*.daa \
  "$OUTDIR"/*tmp* \
  "$OUTDIR"/*temp* 2>/dev/null || true

echo "[$(date)] Done $ID"
