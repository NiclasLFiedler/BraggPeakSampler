#!/bin/bash

TARGERT=0
jq --argjson size "$TARGERT" '.targetSelect = $size' ../../analysis/config.json > temp.json && mv temp.json ../../analysis/config.json

cd build
make -j8
./braggsampler run_p.mac

# Check if the Geant4 simulation ran successfully
if [ $? -eq 0 ]; then
  echo "Geant4 simulation completed successfully."
else
  echo "Error: Geant4 simulation failed."
  exit 1  # Exit the script with a non-zero status code to indicate failure
fi

cd ../data
./merge_p.sh
mv raw_data_p.root ../../../data/simulation/notarget/input/notarget.root
cd ../../../analysis
root -q 'analysis.cpp()'
