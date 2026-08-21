#!/bin/bash

cd build
make -j8
# List of beam energies
beam_energies=("3" "5" "10" "15" "20" "30" "40" "50" "60" "70" "80" "90" "100" "125" "150" "175" "200" "225" "250")
# beam_energies=("50" "100" "150" "200")
# Loop over each energy
for energy in "${beam_energies[@]}"; do
    # Modify the beam energy in your macro file or set it directly
    sed -i "s/\/set\/beamenergy [0-9]* MeV/\/set\/beamenergy $energy MeV/" run_p.mac

    # Set a unique output filename based on the energy
    output_file="lung/range_${energy}.root"

    # Run the simulation with the modified macro file
    ./braggtheory run_p.mac
    cd ../data_analysis
    ./merge_p.sh
    mv "raw_data_p.root" $output_file
        
    cd ../build
    echo "Simulation with $energy MeV completed. Output saved"
done
cd ../data_analysis
activate_venv
python3 analyseProjected.py
python3 range_energy.py
rm lung/range_*.root

