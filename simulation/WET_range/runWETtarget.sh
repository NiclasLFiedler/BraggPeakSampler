#!/bin/bash

# -------- CONFIG --------
JSON_IN="../../../analysis/config.json"
JSON_OUT="../../../analysis/config_tmp.json"

# -------- PARAMETER LISTS --------
nLayers=(23 24 25 26)
energies=(205 210 215 220)
nolayer=(0)

cd build
make -j8

# -------- LOOP (ZIPPED) --------
for (( i=0; i<${#nLayers[@]}; i++ )); do

  nLayer="${nLayers[i]}"
  energy="${energies[i]}"

  echo "========================================"
  echo "Run $((i+1)):"
  echo "  energies     = $energy"
  echo "  nLayers      = $nLayer"
  echo "========================================"

  jq \
    --argjson ts "$energy" \
    --argjson tt "$nLayer" \
    '
    (.detectors[] | select(.detectorID == 0) | .beamEnergy) = $ts
    | (.detectors[] | select(.detectorID == 0) | .nLayers) = $tt
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

# -------- LOOP (ZIPPED) --------
for (( i=0; i<${#nLayers[@]}; i++ )); do

  energy="${energies[i]}"
  nLayer="${nolayer[0]}"
  echo "========================================"
  echo "Run $((i+1)):"
  echo "  energies     = $energy"
  echo "  nLayers      = $nLayer"
  echo "========================================"

  jq \
    --argjson ts "$energy" \
    --argjson tt "$nLayer" \
    '
    (.detectors[] | select(.detectorID == 0) | .beamEnergy) = $ts
    | (.detectors[] | select(.detectorID == 0) | .nLayers) = $tt
    ' \
    "$JSON_IN" > "$JSON_OUT"  && mv "$JSON_OUT" "$JSON_IN"
      

    output_file="raw_data_${energy}.root"

    # Run the simulation with the modified macro file
    ./braggtheory run_p.mac
    cd ../data_analysis
    ./merge_p.sh
    mv "raw_data_p.root" $output_file
    root -l -b -q "data_analysis_p.cpp(${energy})"
    cd ../build
    echo "Simulation for $energy MeV completed. Output saved"
done