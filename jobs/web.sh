#!/bin/bash
# Step 3 on albedo: web map assets (tiles, vector layers) from the full-resolution grids.
# Usage (from the repo root): sbatch jobs/web.sh
#SBATCH --job-name=cryomaps-web
#SBATCH --account=envi.p_forclima
#SBATCH --partition=smp
#SBATCH --qos=12h
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x-%j.out

set -euo pipefail
julia --project=. scripts/web.jl full
