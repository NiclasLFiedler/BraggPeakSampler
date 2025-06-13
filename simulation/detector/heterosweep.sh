#!/bin/bash

CONFIG_FILE="../../analysis/config.json"
#SCRIPT_TO_RUN="./homotarget_run.sh"  # Change to your script's name
SCRIPT_TO_RUN="./heterotarget_run.sh"  # Change to your script's name
NOTARGET_RUN="./notarget_run.sh"
# Update JSON values function
update_json() {
    key=$1
    value=$2
    jq ".$key = $value" "$CONFIG_FILE" > tmp.json && mv tmp.json "$CONFIG_FILE"
}

# Define loop ranges
source ~/venv/bin/activate

for value1 in 0; do  # Change these as needed
  for value2 in 0; do  # Change these as needed
      echo "Setting config.json: param1=$value1, param2=$value2"
      
      # Modify JSON file
      update_json "pmod" "$value1"
      update_json "heteroThickness" "$value2"

      # Execute the other script
      bash "$NOTARGET_RUN"
      cd ../../analysis
      python3 braggfit.py
      cd ../simulation/detector
      cp ../../data/simulation/notarget/output/notargetMeans.root ../../data/heterosweep/output/"$value1"um_"$value2"mmMeans.root
      cp ../../data/simulation/notarget/output/notarget_fit.root ../../data/heterosweep/output/"$value1"um_"$value2"mmFit.root
  done
done

for value1 in 100 200 300 400 500 600 700 800; do  # Change these as needed
    for value2 in 50 100 150 200; do  # Change these as needed
        echo "Setting config.json: param1=$value1, param2=$value2"
        
        # Modify JSON file
        update_json "pmod" "$value1"
        update_json "heteroThickness" "$value2"

        # Execute the other script
        bash "$SCRIPT_TO_RUN"
        cd ../../analysis
        python3 braggfit.py
        cd ../simulation/detector
        cp ../../data/simulation/heterotarget/output/heterotargetMeans.root ../../data/heterosweep/output/"$value1"um_"$value2"mmMeans.root
        cp ../../data/simulation/heterotarget/output/heterotarget_fit.root ../../data/heterosweep/output/"$value1"um_"$value2"mmFit.root
    done
done

#cd ../../analysis
#python3 heterosweep.py

deactivate
