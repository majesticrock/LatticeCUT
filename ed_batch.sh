#!/bin/bash

# Set architecture: "IceLake" or "CascadeLake"
arch="CascadeLake"
LATTICE_TYPE="bcc"

input_file="params/for_auto.txt"
readarray -t NEW_VALUES < "${input_file}"

declare -A TOKENS
TOKENS=(
  [1]="phonon_coupling"
  [g]="phonon_coupling"
  [2]="local_interaction"
  [U]="local_interaction"
  [3]="fermi_energy"
  [E]="fermi_energy"
  [4]="omega_debye"
  [w]="omega_debye"
)


echo "Select the parameter that is to be varied:"
echo "1 / g: phonon_coupling"
echo "2 / U: local_interaction"
echo "3 / E: fermi_energy"
echo "4 / w: omega_debye"
read -p "Enter your choice: " choice
TOKEN=${TOKENS[$choice]}

if [ -z "$TOKEN" ]; then
  echo "Invalid choice. Exiting."
  exit 1
fi

CURRENT_TIME=$(date +"%Y%m%d_%H%M%S")

rm -rf auto_generated_${CURRENT_TIME}/
mkdir -p auto_generated_${CURRENT_TIME}

for NEW_VALUE in "${NEW_VALUES[@]}"; do
  NEW_NAME=$(echo "$NEW_VALUE" | sed 's/ /_/g')

  # Replace the config value
  while read line; do
    if [[ $line == \#* ]]; then
      continue
    fi

    TOKEN_NAME=$(echo "$line" | awk '{print $1}')
    TOKEN_VALUE=$(echo "$line" | cut -d' ' -f2-)

    if [[ "$TOKEN_NAME" == "$TOKEN" ]]; then
      sed "s/$TOKEN_NAME $TOKEN_VALUE/$TOKEN_NAME $NEW_VALUE/" params/cluster/${LATTICE_TYPE}.config > auto_generated_${CURRENT_TIME}/$NEW_NAME.config
      break
    fi
  done < params/cluster/${LATTICE_TYPE}.config

  # Create the slurm file with modifications
  slurm_path="auto_generated_${CURRENT_TIME}/$NEW_NAME.slurm"
  sed -e "s|#SBATCH --job-name=${LATTICE_TYPE}_ed|#SBATCH --job-name=${LATTICE_TYPE}_${NEW_NAME}_${CURRENT_TIME}|" \
      -e "s|#SBATCH --output=/home/althueser/phd/cpp/LatticeCUT/output_${LATTICE_TYPE}.txt|#SBATCH --output=/home/althueser/phd/cpp/LatticeCUT/output_${CURRENT_TIME}_$NEW_NAME.txt|" \
      -e "s|^#SBATCH --constraint=.*|#SBATCH --constraint=${arch}|" \
      -e "s|./build_.*/latticecut .*|./build_${arch}/latticecut auto_generated_${CURRENT_TIME}/$NEW_NAME.config|" \
      slurm/${LATTICE_TYPE}_ed.slurm > "$slurm_path"

  # Submit the job
  sbatch "$slurm_path"
done
