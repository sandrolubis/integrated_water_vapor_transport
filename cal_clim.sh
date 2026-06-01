#!/bin/bash
set -e

qdir="/pscratch/sd/s/slubis/ERA5_Data/SH/daily_mean"
udir="/pscratch/sd/s/slubis/ERA5_Data/U/daily_mean"
vdir="/pscratch/sd/s/slubis/ERA5_Data/V/daily_mean"

outdir="/pscratch/sd/s/slubis/ERA5_Data/daily_climatology"
mkdir -p "$outdir"

yr1=1979
yr2=2025
nyears=$((yr2 - yr1 + 1))

for var in q u v; do

    if [ "$var" = "q" ]; then
        indir="$qdir"
    elif [ "$var" = "u" ]; then
        indir="$udir"
    else
        indir="$vdir"
    fi

    filelist=""
    count=0

    for yr in $(seq $yr1 $yr2); do
        file="${indir}/${var}.${yr}.nc"

        if [ -f "$file" ]; then
            filelist="$filelist $file"
            count=$((count + 1))
        else
            echo "Missing file: $file"
            exit 1
        fi
    done

    echo "======================================"
    echo "Variable: $var"
    echo "Number of files found: $count"
    echo "Expected number of files: $nyears"
    echo "======================================"

    if [ "$count" -ne "$nyears" ]; then
        echo "ERROR: Not all years are included for $var"
        exit 1
    fi

    merged="${outdir}/${var}_merged_${yr1}_${yr2}.nc"
    clim="${outdir}/${var}_daily_clim_${yr1}_${yr2}.nc"

    echo "Merging $var from $yr1 to $yr2..."
    cdo -O mergetime $filelist "$merged"

    echo "Checking merged file time length..."
    cdo ntime "$merged"

    echo "Computing daily climatology..."
    cdo -O ydaymean "$merged" "$clim"

    echo "Checking daily climatology time length..."
    cdo ntime "$clim"

    echo "Saved: $clim"

donecal_clim.sh
