#!/bin/bash

# -------- CONFIG --------
JSON_IN="../../analysis/config.json"
JSON_OUT="../../analysis/config_tmp.json"

# -------- PARAMETER LISTS --------
TARGET_SELECT_LIST=(0 2)
TARGET_THICKNESS_LIST=(0 200)
RESOLUTION=(4000)

# TARGET_SELECT_LIST=(1)
# TARGET_THICKNESS_LIST=(52)
# -------- SAFETY CHECK --------
if [ "${#TARGET_SELECT_LIST[@]}" -ne "${#TARGET_THICKNESS_LIST[@]}" ]; then
  echo "Error: TARGET_SELECT_LIST and TARGET_THICKNESS_LIST must have same length"
  exit 1
fi

# -------- LOOP (ZIPPED) --------
for (( i=0; i<${#TARGET_SELECT_LIST[@]}; i++ )); do

  TARGET_SELECT="${TARGET_SELECT_LIST[i]}"
  TARGET_THICKNESS="${TARGET_THICKNESS_LIST[i]}"
  SIMULATION_RESOLUTION="${RESOLUTION}"
  echo "========================================"
  echo "Run $((i+1)):"
  echo "  targetSelect     = $TARGET_SELECT"
  echo "  targetThickness  = $TARGET_THICKNESS"
  echo "  RESOLUTION  = $RESOLUTION"
  echo "========================================"

  jq \
    --argjson ts "$TARGET_SELECT" \
    --argjson tt "$TARGET_THICKNESS" \
    --argjson sr "$SIMULATION_RESOLUTION" \
    '
    .targetSelect = $ts
    | (.detectors[] | select(.detectorID == 0) | .targetThickness) = $tt
    | (.detectors[] | select(.detectorID == 0) | .nLayers) = $sr
    ' \
    "$JSON_IN" > "$JSON_OUT"  && mv "$JSON_OUT" "$JSON_IN"
      

    if [ $? -ne 0 ]; then
      echo "JSON modification failed — skipping"
      continue
    fi
    cd build
    make -j8
    ./idealDD run_p.mac
    cd ..
    python3 depthAnalysis.py
done

python3 convolution.py