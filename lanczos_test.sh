#!/bin/bash

INPUT_CONFIG="param/cluster/sc.config"
SLURM_TEMPLATE="slurm/sc_cascade.slurm"

# Timestamped output directory
DATE_TAG=$(date +"%Y%m%d_%H%M")
OUTPUT_DIR="auto_generated_lanczos_${DATE_TAG}"

DOS_OPTIONS=("sc" "bcc" "fcc")
N_VALUES=(2000 4000 8000 12000 20000)

# Memory setup
declare -A MEM_MAP
MEM_MAP[2000]="8gb"
MEM_MAP[4000]="8gb"
MEM_MAP[8000]="24gb"
MEM_MAP[12000]="32gb"
MEM_MAP[20000]="98gb"

# CPU setup
declare -A CPU_MAP
CPU_MAP[2000]="4"
CPU_MAP[4000]="8"
CPU_MAP[8000]="12"
CPU_MAP[12000]="40"
CPU_MAP[20000]="40"

mkdir -p "$OUTPUT_DIR"

for dos in "${DOS_OPTIONS[@]}"; do
    for Nval in "${N_VALUES[@]}"; do
        
        CONFIG_FILE="${OUTPUT_DIR}/${dos}_N${Nval}.config"
        SLURM_FILE="${OUTPUT_DIR}/${dos}_N${Nval}.slurm"
        LOG_FILE="${OUTPUT_DIR}/output_${dos}_N${Nval}_${DATE_TAG}.txt"

        MEM="${MEM_MAP[$Nval]}"
        CPUS="${CPU_MAP[$Nval]}"

        # --- Generate modified config file ---
        sed \
            -e "s/^dos .*/dos $dos/" \
            -e "s/^N .*/N $Nval/" \
            "$INPUT_CONFIG" > "$CONFIG_FILE"

        echo "Generated config: $CONFIG_FILE"

        # --- Generate slurm batch script ---
        sed \
            -e "s#./build_CascadeLake/latticecut .*#./build_CascadeLake/latticecut ${CONFIG_FILE}#" \
            -e "s/^#SBATCH --mem=.*/#SBATCH --mem=${MEM}/" \
            -e "s/^#SBATCH --cpus-per-task=.*/#SBATCH --cpus-per-task=${CPUS}/" \
            -e "s#^#SBATCH --output=.*#SBATCH --output=${LOG_FILE}#" \
            -e "s/^#SBATCH --job-name=.*/#SBATCH --job-name=${dos}_N${Nval}/" \
            "$SLURM_TEMPLATE" > "$SLURM_FILE"

        echo "Generated slurm: $SLURM_FILE"
        echo " → Log will write to: $LOG_FILE"

        # Submit job
        sbatch "$SLURM_FILE"
        echo "Submitted: dos=$dos | N=$Nval | CPUs=$CPUS | Mem=$MEM"
        echo
    done
done