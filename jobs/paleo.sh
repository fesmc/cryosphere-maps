#!/bin/bash
# Paleo maps on albedo: the grids from PaleoMIST (scripts/prepare_paleo.jl) and the
# small PNGs of both raster styles (scripts/paleo.jl).
# Usage (from the repo root): sbatch jobs/paleo.sh [time=20] [nh] [antarctica]
#SBATCH --job-name=cryomaps-paleo
#SBATCH --account=envi.p_forclima
#SBATCH --partition=smp
#SBATCH --qos=12h
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=logs/%x-%j.out

set -euo pipefail
t=time=20; regions=()
for a in "$@"; do
    case $a in time=*) t=$a ;; *) regions+=("$a") ;; esac
done
[ ${#regions[@]} -eq 0 ] && regions=(nh antarctica)
julia --project=. scripts/prepare_paleo.jl "$t" "${regions[@]}"
for r in "${regions[@]}"; do
    for style in surface bed; do
        julia --project=. scripts/paleo.jl "$r" "$t" "$style"
    done
done
