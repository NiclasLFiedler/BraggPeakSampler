#!/bin/bash

TARGERT=0
jq --argjson size "$TARGERT" '.targetSelect = $size' ../../analysis/config.json > temp.json && mv temp.json ../../analysis/config.json

cd build
make -j8
./braggsampler
