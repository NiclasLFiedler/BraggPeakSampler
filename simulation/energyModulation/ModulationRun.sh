#!/bin/bash

JSON_IN="../../analysis/config.json"
JSON_OUT="../../analysis/config_tmp.json"

DETECTOR_MATERIALS=("h2o")

ENERGY=220
DETECTOR_SELECT=0
# ── helper: apply jq patch, run simulation, analyse ──────────────────────────
run_simulation() {
    local mat="$1"
    local modulation="$2"
    local target_select="$3"
    local target_thickness="$4"

    echo "========================================"
    echo "  detectorType     = $mat"
    echo "  detectorSelect   = $DETECTOR_SELECT"
    echo "  targetSelect     = $target_select"
    echo "  targetThickness  = $target_thickness"
    echo "  modulation       = $modulation"
    echo "========================================"

    jq \
        --argjson ds  "$DETECTOR_SELECT" \
        --argjson ts  "$target_select" \
        --argjson tt  "$target_thickness" \
        --argjson mod "$modulation" \
        --arg     mat "$mat" \
        '
        .detectorSelect = $ds
        | .targetSelect = $ts
        | .pmod = $mod
        | (.detectors[] | select(.detectorID == 0) | .detectorType)    = $mat
        | (.detectors[] | select(.detectorID == 0) | .targetThickness) = $tt
        ' \
        "$JSON_IN" > "$JSON_OUT" && mv "$JSON_OUT" "$JSON_IN"

    if [ $? -ne 0 ]; then
        echo "JSON modification failed — skipping"
        return 1
    fi

    cd build
    make -j
    ./braggtheory run_p.mac
    cd ../data_analysis
    ./merge_p.sh
    output_file="data_${modulation}_${target_thickness}.root"
    mv "raw_data_p.root" $output_file
    mv $output_file data
    cd ..
    # python3 depthAnalysis.py
}

# ── main loops ────────────────────────────────────────────────────────────────
for (( m=0; m<${#DETECTOR_MATERIALS[@]}; m++ )); do

    MAT="${DETECTOR_MATERIALS[m]}"

    echo ""
    echo "###############################################"
    echo "  Detector material : $MAT   crystalSize[2] = $SIZE"
    echo "###############################################"

    # --- no-target run (targetSelect=0, targetThickness=0) ---
    echo ""
    echo "--- No-target run ---"
    run_simulation "$MAT" 0 0 0

    # --- thickness sweep: 10 → 300 in steps of 10 (targetSelect=0) ---
    for (( t=100; t<=200; t+=100 )); do
        for (( mod=100; mod<=800; mod+=100 )); do
            echo ""
            echo "--- Thickness sweep: $t mm --- mod = $mod ---"
            run_simulation "$MAT" "$mod" 2 "$t"
        done
    done
done

# ── final convolution ─────────────────────────────────────────────────────────
# python3 convolution.py