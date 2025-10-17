
TH2D* TransposeTH2D(TH2D* original) {
  // Get binning and axis info
  int nBinsX = original->GetNbinsX();
  int nBinsY = original->GetNbinsY();
  double xMin = original->GetXaxis()->GetXmin();
  double xMax = original->GetXaxis()->GetXmax();
  double yMin = original->GetYaxis()->GetXmin();
  double yMax = original->GetYaxis()->GetXmax();

  // Create transposed histogram
  TH2D* transposed = new TH2D(
    TString(original->GetName()) + "_transposed",
    TString(original->GetTitle()) + " (Transposed);" +
    original->GetYaxis()->GetTitle() + ";" +
    original->GetXaxis()->GetTitle(),
    nBinsY, yMin, yMax,  // New X-axis (old Y-axis)
    nBinsX, xMin, xMax   // New Y-axis (old X-axis)
  );

  // Copy bin contents and errors
  for (int i = 1; i <= nBinsX; i++) {
    for (int j = 1; j <= nBinsY; j++) {
      double content = original->GetBinContent(i, j);
      double error = original->GetBinError(i, j);
      transposed->SetBinContent(j, i, content); // Swap indices
      transposed->SetBinError(j, i, error);
    }
  }

  // Copy statistics (entries, etc.)
  transposed->SetEntries(original->GetEntries());
  return transposed;
}

void RootUnfold(){
    // Open the ROOT file
    TFile* file = TFile::Open("responseMatrix.root");

    // Check if file opened successfully
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file!" << std::endl;
        return;
    }

    // List contents to see what's inside
    file->ls();

    // Load measured data and response matrix
    TH2D* hResponse = (TH2D*)file->Get("h_responseMatrix");
    TH1D* hMeasData = (TH1D*)file->Get("h_measurement");

    TH2D* hResponse_T = TransposeTH2D(hResponse);

    TCanvas *c2 = new TCanvas("c2", "Unfolding", 800, 600);
    // Set up response object
    hMeasData->Draw("COLZ");

    RooUnfoldResponse response(nullptr, nullptr, hResponse);

    // Unfold with Bayesian method (4 iterations)
    RooUnfoldSvd unfold(&response, hMeasData, 4);
    TH1D* hUnfolded = (TH1D*)unfold.Hunfold();

    // Draw results
    TCanvas *c1 = new TCanvas("c1", "Unfolding", 800, 600);
    hMeasData->SetLineColor(kRed);
    hMeasData->SetLineWidth(2);
    hUnfolded->SetLineColor(kBlue);
    hUnfolded->SetLineWidth(2);
    hMeasData->Draw("HIST");
    hUnfolded->Draw("HIST SAME");
    c1->BuildLegend();

    // Optional: Save unfolded spectrum
    TFile *fout = new TFile("unfolded.root", "RECREATE");
    hUnfolded->Write("unfolded_spectrum");
    fout->Close();
}