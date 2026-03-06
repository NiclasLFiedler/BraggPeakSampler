#!/bin/bash

# -------- CONFIG --------
JSON_IN="../../../analysis/config.json"
JSON_OUT="../../../analysis/config_tmp.json"

# -------- PARAMETER LISTS --------
#notarget
#nLayers=(28 29 30 31 32)
#energies=(225 230 235 240 245)

#pmma 52mm 195 MeV
nLayers=(22)
# energies=(195 200 205 210 215 220)
energies=(196)
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
    rm $output_file
    rangeoutput="range_${nLayer}.root"
    mv $rangeoutput targettestpbwo4/
    cd ../build
    echo "Simulation for $nLayer completed. Output saved"
done

#notarget
#energies=(221.6 225 230 235 240 245)

#pmma 52 194MeV
energies=(193.7 196 199 202 205 208)

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
    mv $rangeoutput targettest/
    cd ../build
    echo "Simulation for $energy MeV completed. Output saved"
done