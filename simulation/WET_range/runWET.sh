#!/bin/bash

# -------- CONFIG --------
JSON_IN="../../../analysis/config.json"
JSON_OUT="../../../analysis/config_tmp.json"
cd build
make -j8

energies=(221.6 193.72 192.62 192.15 191.68 190.65 190.03)
nolayer=(0)
nLayers=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31)
thickness=(0 52 53 54 55 56 57)
# -------- LOOP (ZIPPED) --------
for (( i=0; i<${#energies[@]}; i++ )); do

  energy="${energies[i]}"
  energy_label=${energy%.*}  # Truncates "221.6" to "221"
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
      

    output_file="raw_data_${energy_label}.root"

    # Run the simulation with the modified macro file
    ./braggtheory run_p.mac
    cd ../data_analysis
    ./merge_p.sh
    mv "raw_data_p.root" $output_file
    root -l -b -q "data_analysis_p.cpp(${energy})"
    rm $output_file
    rangeoutput="range_${energy_label}.root"
    fileoutput="WET${thickness[i]}"
    mv $rangeoutput $fileoutput
    cd ../build
    echo "Simulation for $energy MeV completed. Output saved"

  for (( j=0; j<${#nLayers[@]}; j++ )); do
    nLayer="${nLayers[j]}"

    echo "========================================"
    echo "Run $((j+1)):"
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
    rm $output_file
    rangeoutput="range_${nLayer}.root"
    fileoutput="WET${thickness[i]}"
    mv $rangeoutput $fileoutput
    cd ../build
    echo "Simulation for $nLayer completed. Output saved"

  done
done