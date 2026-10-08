#!/bin/bash
# Step 1 on albedo: regrid raw data onto the poster grids.
# Usage (from the repo root): sbatch jobs/prepare.sh [greenland] [antarctica]
#SBATCH --job-name=cryomaps-prepare
#SBATCH --account=envi.p_forclima
#SBATCH --partition=smp
#SBATCH --qos=12h
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=logs/%x-%j.out

set -euo pipefail
julia --project=. scripts/prepare.jl "$@"
