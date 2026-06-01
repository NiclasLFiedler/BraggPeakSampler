#!/bin/bash

JSON_IN="../../analysis/config.json"
JSON_OUT="../../analysis/config_tmp.json"

DETECTOR_MATERIALS=("pbwo4" "h2o")
DETECTOR_SIZES=(8 40)

RESOLUTION=300
DETECTOR_SELECT=0
# ── helper: apply jq patch, run simulation, analyse ──────────────────────────
run_simulation() {
    local mat="$1"
    local size="$2"
    local target_select="$3"
    local target_thickness="$4"

    echo "========================================"
    echo "  detectorType     = $mat"
    echo "  crystalSize[2]   = $size"
    echo "  detectorSelect   = $DETECTOR_SELECT"
    echo "  targetSelect     = $target_select"
    echo "  targetThickness  = $target_thickness"
    echo "  nLayers          = $RESOLUTION"
    echo "========================================"

    jq \
        --argjson ds  "$DETECTOR_SELECT" \
        --argjson ts  "$target_select" \
        --argjson tt  "$target_thickness" \
        --argjson sr  "$RESOLUTION" \
        --arg     mat "$mat" \
        --argjson sz  "$size" \
        '
        .detectorSelect = $ds
        | .targetSelect = $ts
        | (.detectors[] | select(.detectorID == 0) | .detectorType)    = $mat
        | (.detectors[] | select(.detectorID == 0) | .crystalSize[2])  = $sz
        | (.detectors[] | select(.detectorID == 0) | .targetThickness) = $tt
        | (.detectors[] | select(.detectorID == 0) | .nLayers)         = $sr
        ' \
        "$JSON_IN" > "$JSON_OUT" && mv "$JSON_OUT" "$JSON_IN"

    if [ $? -ne 0 ]; then
        echo "JSON modification failed — skipping"
        return 1
    fi

    cd build
    make -j8
    ./idealDD run_p.mac
    cd ..
    ./merge_p.sh
    python3 depthAnalysis.py
}

# ── main loops ────────────────────────────────────────────────────────────────
for (( m=0; m<${#DETECTOR_MATERIALS[@]}; m++ )); do

    MAT="${DETECTOR_MATERIALS[m]}"
    SIZE="${DETECTOR_SIZES[m]}"

    echo ""
    echo "###############################################"
    echo "  Detector material : $MAT   crystalSize[2] = $SIZE"
    echo "###############################################"

    # --- no-target run (targetSelect=0, targetThickness=0) ---
    echo ""
    echo "--- No-target run ---"
    run_simulation "$MAT" "$SIZE" 0 0

    # --- thickness sweep: 10 → 300 in steps of 10 (targetSelect=0) ---
    for (( t=10; t<=250; t+=10 )); do
        echo ""
        echo "--- Thickness sweep: $t mm ---"
        run_simulation "$MAT" "$SIZE" 2 "$t"
    done
done

# ── final convolution ─────────────────────────────────────────────────────────
# python3 convolution.py