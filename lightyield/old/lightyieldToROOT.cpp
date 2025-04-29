#include <TFile.h>
#include <TH1F.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <sstream>

void ToROOT(const char* inputFileName, const char* outputFileName) {
    std::ifstream inputFile(inputFileName);
    if (!inputFile.is_open()) {
        std::cerr << "Error opening input file: " << inputFileName << std::endl;
        return;
    }

    // Read the data and find the maximum channel number
    int maxChannel = 0;
    int channel;
    float count;
    while (inputFile >> channel >> count) {
        if (channel > maxChannel) {
            maxChannel = channel;
        }
    }
    inputFile.clear(); // Clear EOF flag
    inputFile.seekg(0); // Rewind file to read again

    // Create a histogram with a number of bins equal to maxChannel + 1
    TH1F* hist = new TH1F("hist", "Channel Histogram", maxChannel + 1, 0, maxChannel + 1);

    // Fill the histogram with data from the file
    while (inputFile >> channel >> count) {
        hist->SetBinContent(channel + 1, count); // ROOT histograms are 1-based indexing
    }

    // Close the input file
    inputFile.close();

    // Create the output ROOT file
    TFile* outFile = TFile::Open(outputFileName, "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error creating output file: " << outputFileName << std::endl;
        return;
    }

    // Write the histogram to the file
    hist->Write();
    outFile->Close();

    std::cout << "Histogram saved to " << outputFileName << std::endl;

    // Cleanup
    delete hist;
}
void lightyieldToROOT(){
    Char_t filename[100]; // temp store 
    Char_t output[100]; // temp store 
    for(int i = 0; i< 15; i++){
        sprintf(filename, "BPS_CH%i.his.txt", i);
        sprintf(output, "BPS_CH%i.root", i);
        ToROOT(filename, output);
    }
}