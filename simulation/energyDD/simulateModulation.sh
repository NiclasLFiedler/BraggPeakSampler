#!/bin/bash

JSON_IN="../../analysis/config.json"
JSON_OUT="../../analysis/config_tmp.json"

TARGET_SELECT_LIST=(0 2 2)
TARGET_THICKNESS_LIST=(0 100 50)
RESOLUTION=(500)

for (( j=0; j<${#RESOLUTION[@]}; j++ )); do
  for (( i=0; i<${#TARGET_SELECT_LIST[@]}; i++ )); do

    TARGET_SELECT="${TARGET_SELECT_LIST[i]}"
    TARGET_THICKNESS="${TARGET_THICKNESS_LIST[i]}"
    SIMULATION_RESOLUTION="${RESOLUTION[j]}"
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
      ./merge_p.sh
      python3 depthAnalysis.py
  done
done

python3 convolution.py