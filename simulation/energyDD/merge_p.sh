#!/bin/bash

# Run hadd command to merge all data_t*.root files into output_data.root
cd data/temp
hadd -f data.root data_t*.root
rm data_t*.root
cd ../..
# Check if hadd command was successful
if [ $? -eq 0 ]; then
  echo "Merging successful. Output file: output_data.root"
else
  echo "Error: Merging failed."
fi