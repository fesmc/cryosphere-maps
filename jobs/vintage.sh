#!/bin/bash
# The vintage maps at full resolution (plots/vintage/: PDF, 100 dpi PNG and small PNG).
# Usage (from the repo root): sbatch jobs/vintage.sh
#SBATCH --job-name=cryomaps-vintage
#SBATCH --account=envi.p_forclima
#SBATCH --partition=smp
#SBATCH --qos=12h
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/%x-%j.out

set -euo pipefail
for region in greenland antarctica; do
    julia --project=. scripts/vintage.jl "$region" full
done
