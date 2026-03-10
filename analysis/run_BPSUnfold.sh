#!/bin/bash

# -------- CONFIG --------
JSON_IN="config.json"
JSON_OUT="config_tmp.json"

# -------- PARAMETER LISTS --------
TARGET_SELECT_LIST=(0 1 1 1 1 1 1)
TARGET_THICKNESS_LIST=(0 52 53 54 55 56 57)

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

  echo "========================================"
  echo "Run $((i+1)):"
  echo "  targetSelect     = $TARGET_SELECT"
  echo "  targetThickness  = $TARGET_THICKNESS"
  echo "========================================"

  jq \
    --argjson ts "$TARGET_SELECT" \
    --argjson tt "$TARGET_THICKNESS" \
    '
    .targetSelect = $ts
    | (.detectors[] | select(.detectorID == 0) | .targetThickness) = $tt
    ' \
    "$JSON_IN" > "$JSON_OUT"  && mv "$JSON_OUT" "$JSON_IN"
      

    if [ $? -ne 0 ]; then
      echo "JSON modification failed — skipping"
      continue
    fi
    # root -q 'analysis.cpp()'
    python3 BPSUnfolding.py
done