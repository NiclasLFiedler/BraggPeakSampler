#!/bin/bash

# Run hadd command to merge all data_t*.root files into output_data.root
hadd -f raw_data_muon.root data_t*.root
rm data_t*.root
# Check if hadd command was successful
if [ $? -eq 0 ]; then
  echo "Merging successful. Output file: output_data.root"
else
  echo "Error: Merging failed."
fi