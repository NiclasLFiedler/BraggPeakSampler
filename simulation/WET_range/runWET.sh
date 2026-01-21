#!/bin/bash

# -------- CONFIG --------
JSON_IN="../../../analysis/config.json"
JSON_OUT="../../../analysis/config_tmp.json"

# -------- PARAMETER LISTS --------
nLayers=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23)
# 24 25 26 27 28 29 30 31 32)

cd build
make -j8

# -------- LOOP (ZIPPED) --------
for (( i=0; i<${#nLayers[@]}; i++ )); do

  nLayer="${nLayers[i]}"

  echo "========================================"
  echo "Run $((i+1)):"
  echo "  nLayers     = $nLayer"
  echo "========================================"

  jq \
    --argjson tt "$nLayer" \
    '
    (.detectors[] | select(.detectorID == 0) | .nLayers) = $tt
    ' \
    "$JSON_IN" > "$JSON_OUT"  && mv "$JSON_OUT" "$JSON_IN"
      

    output_file="raw_data_${nLayer}.root"

    # Run the simulation with the modified macro file
    ./braggtheory run_p.mac
    cd ../data_analysis
    ./merge_p.sh
    mv "raw_data_p.root" $output_file
    root -l -b -q "data_analysis_p.cpp(${nLayer})"
    cd ../build
    echo "Simulation for $nLayer completed. Output saved"

done