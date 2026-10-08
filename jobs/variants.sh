#!/bin/bash
# The poster variants shown in the website gallery, at full resolution
# (plots/variants/; only their small PNGs are tracked).
# Usage (from the repo root): sbatch jobs/variants.sh
#SBATCH --job-name=cryomaps-variants
#SBATCH --account=envi.p_forclima
#SBATCH --partition=smp
#SBATCH --qos=12h
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x-%j.out

set -euo pipefail
for opt in cmap=ember dark surface bed tier=2; do
    julia --project=. scripts/greenland.jl full $opt
    julia --project=. scripts/antarctica.jl full $opt
done
