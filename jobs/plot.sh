#!/bin/bash
# Step 2 on albedo: render the A0 posters from the prepared grids.
# Usage (from the repo root): sbatch jobs/plot.sh [velocity|surface|bed]
#SBATCH --job-name=cryomaps-plot
#SBATCH --account=envi.p_forclima
#SBATCH --partition=smp
#SBATCH --qos=12h
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x-%j.out

set -euo pipefail
julia --project=. scripts/greenland.jl "$@"
julia --project=. scripts/antarctica.jl "$@"
