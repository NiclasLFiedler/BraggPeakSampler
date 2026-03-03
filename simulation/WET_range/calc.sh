#!/bin/bash

./runWET.sh
./runWETtarget.sh
cd data_analysis
python3 range_energy.py