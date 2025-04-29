import ROOT
import sys
from array import array

def lightyieldToROOT(input_file, output_file):
    # Create a new ROOT file
    root_file = ROOT.TFile(output_file, "RECREATE")
    
    # Create a ROOT tree
    tree = ROOT.TTree("tree", "Histogram data")

    # Variable to hold the channel number
    channel = array('i', [0])

    # Create a branch in the tree
    tree.Branch("channel", channel, "channel/I")

    # Open the input text file
    with open(input_file, 'r') as f:
        for line in f:
            # Split each line into channel number and count
            parts = line.split()
            if len(parts) == 2:
                channel_nb = int(parts[0])
                count_nb = int(parts[1])
                
                # Fill the branch with the channel number, repeating it `count_nb` times
                for _ in range(count_nb):
                    channel[0] = channel_nb
                    tree.Fill()

    # Write the tree to the ROOT file and close the file
    tree.Write()
    root_file.Close()
    
def main():
   for i in range(0, 15):
        filename = f'BPS_CH{i}.his.txt'
        lightyieldToROOT(filename, filename[:-8]+".root")

if __name__ == "__main__":
    main()

