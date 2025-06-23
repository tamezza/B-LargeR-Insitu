#include "../../Util/AtlasStyle.C"
#include <iostream>
#include <string>
#include <vector>
#include <sstream> // For std::stringstream
#include <iomanip> // For std::fixed, std::setprecision, std::setw
#include <fstream> // For ofstream
#include <cmath>   // For std::pow and std::sqrt
#include <algorithm> // For std::max

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TLine.h"
#include "TLatex.h"
#include "TBox.h"
#include "TSystem.h" // For gSystem->Exec
#include "TROOT.h" // For gROOT
#include "TAxis.h" // For TAxis::SetBinLabel

// Function to draw ATLAS label (modified to accept working point)
void DrawAtlasLabel(const char* wp_label) {
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(72);
    latex.SetTextColor(kBlack);

    // Draw "ATLAS" label
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.15, 0.925, "ATLAS"); // Adjusted position for consistency

    // Draw "Internal" label
    latex.SetTextFont(42);
    latex.DrawLatex(0.25, 0.925, "Internal"); // Adjusted position for consistency

    // Draw the energy and luminosity text
    latex.DrawLatex(0.15, 0.885, "#sqrt{s} = 13 TeV, 140 fb^{-1}"); // Adjusted position for consistency
    
    // Draw working point label
    latex.SetTextSize(0.03);
    latex.DrawLatex(0.15, 0.84, Form("Working Point: %s", wp_label)); // Adjusted position
}

// Function to save fit parameters to a file
void SaveFitParameters(const char* filename, const char* fitType, 
    double mean, double meanErr, 
    double sigma, double sigmaErr,
    double amp, double ampErr,
    double chi2, double ndf) {

    ofstream outFile(filename);

    if (!outFile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    outFile << "Fit Type: " << fitType << endl;
    outFile << "Mean: " << mean << " +- " << meanErr << endl;
    outFile << "Sigma: " << sigma << " +- " << sigmaErr << endl;
    outFile << "Amplitude: " << amp << " +- " << ampErr << endl;
    outFile << "Chi2/NDF: " << chi2/ndf << " (" << chi2 << "/" << ndf << ")" << endl;

    outFile.close();
    cout << "Fit parameters saved to " << filename << endl;
}

// Function to create directories for different mass cut scenarios (modified to accept working point)
void CreateDirectories(const char* wp_label) {
    // Base directory for the specific working point
    TString baseDirCmd = Form("mkdir -p ../Plots/MJB/FitRdbNoMcut/%s", wp_label);
    gSystem->Exec(baseDirCmd.Data());

    // Define mass cut labels
    //const int nMassCuts = 6;
    //const char* massCutLabels[nMassCuts] = {"NoCut", "50-150", "50-200", "50-250", "50-300", "50-350"};

    const int nMassCuts = 1;
    const char* massCutLabels[nMassCuts] = {"NoCut"};

    // Create a directory for each mass cut within the working point directory
    for (int i = 0; i < nMassCuts; i++) {
        TString dirCmd = Form("mkdir -p ../Plots/MJB/FitRdbNoMcut/%s/%s", wp_label, massCutLabels[i]);
        gSystem->Exec(dirCmd.Data());
    }
}

void RunRdbAnalysis(const char* wp_label) {
    // Load ATLAS style settings for consistent plotting
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    gROOT->LoadMacro("AtlasLabels.C"); // Keep this if AtlasLabels.C defines other needed functions
    
    // Create directories for different mass cuts
    CreateDirectories(wp_label);

    // Define the mass cuts we want to analyze
    const int nMassCuts = 1;
    const char* massCutLabels[nMassCuts] = {"NoCut"};//, "50-150", "50-200", "50-250", "50-300", "50-350"};
    const double massCutLow[nMassCuts] = {-1};//, 50, 50, 50, 50, 50};  // -1 means no lower cut
    const double massCutHigh[nMassCuts] = {-1};//, 150, 200, 250, 300, 350}; // -1 means no upper cut

    //-------------------------------------------
    // Open input files and get trees (modified paths)
    //-------------------------------------------
    TString mcFileName = Form("../../NtupleSlim/MJB_0p%swp_slim/ZbbJets_0p%swp_slim.root", wp_label, wp_label);
    TString dataFileName = Form("../../NtupleSlim/MJB_0p%swp_slim/data_0p%swp_slim.root", wp_label, wp_label);
    TString mulijetsFileName = Form("../../NtupleSlim/MJB_0p%swp_slim/Multijets_0p%swp_slim.root", wp_label, wp_label);

    TFile *fmc = TFile::Open(mcFileName.Data(), "READ");
    TFile *fdata = TFile::Open(dataFileName.Data(), "READ");
    TFile *fmulijets = TFile::Open(mulijetsFileName.Data(), "READ");

    if (!fmc || fmc->IsZombie()) {
        cerr << "Error: MC input file could not be opened for working point " << wp_label << endl;
        /*if (fmc) delete fmc; // Clean up if partially opened
        if (fdata) delete fdata;
        if (fmulijets) delete fmulijets;*/
        return;
    }
    if (!fdata || fdata->IsZombie()) {
        cerr << "Error: Data input file could not be opened for working point " << wp_label << endl;
        /*if (fmc) delete fmc;
        if (fdata) delete fdata;
        if (fmulijets) delete fmulijets;*/
        return;
    }
    if (!fmulijets || fmulijets->IsZombie()) {
        cerr << "Error: Multijets input file could not be opened for working point " << wp_label << endl;
        /*if (fmc) delete fmc;
        if (fdata) delete fdata;
        if (fmulijets) delete fmulijets;*/
        return;
    }

    TTree *tree_mc = (TTree*)fmc->Get("nominal");
    TTree *tree_data = (TTree*)fdata->Get("nominal");
    TTree *tree_mulijets = (TTree*)fmulijets->Get("nominal");

    if (!tree_mc) {
        cerr << "Error: 'nominal' tree not found in MC file for working point " << wp_label << endl;
        fmc->Close(); fdata->Close(); fmulijets->Close();
        //delete fmc; delete fdata; delete fmulijets;
        return;
    }
    if (!tree_data) {
        cerr << "Error: 'nominal' tree not found in Data file for working point " << wp_label << endl;
        fmc->Close(); fdata->Close(); fmulijets->Close();
        //delete fmc; delete fdata; delete fmulijets;
        return;
    }
    if (!tree_mulijets) {
        cerr << "Error: 'nominal' tree not found in Multijets file for working point " << wp_label << endl;
        fmc->Close(); fdata->Close(); fmulijets->Close();
        //delete fmc; delete fdata; delete fmulijets;
        return;
    }


    // Arrays to store histograms for each mass cut
    TH1F* h_mc_ljet_m1[nMassCuts];
    TH1F* h_data_ljet_m1[nMassCuts];
    TH1F* h_multijets_ljet_m1[nMassCuts];
    TH1F* h_mc_R_DB[nMassCuts];
    TH1F* h_data_R_DB[nMassCuts];
    TH1F* h_multijets_ljet_R_DB[nMassCuts];
    TH1F* h_combined_R_DB[nMassCuts];

    // Arrays to store fit results for comparison plots and summary table
    double mc_means[nMassCuts], mc_mean_errs[nMassCuts];
    double data_means[nMassCuts], data_mean_errs[nMassCuts];
    double bkg_means[nMassCuts], bkg_mean_errs[nMassCuts];
    double combined_means[nMassCuts], combined_mean_errs[nMassCuts];
    
    double mc_sigmas[nMassCuts], mc_sigma_errs[nMassCuts];
    double data_sigmas[nMassCuts], data_sigma_errs[nMassCuts];
    double bkg_sigmas[nMassCuts], bkg_sigma_errs[nMassCuts];
    double combined_sigmas[nMassCuts], combined_sigma_errs[nMassCuts];

    double mc_chi2_ndf[nMassCuts];
    double data_chi2_ndf[nMassCuts];
    double bkg_chi2_ndf[nMassCuts];
    double combined_chi2_ndf[nMassCuts];
    
    // Initialize histograms for each mass cut
    for (int m = 0; m < nMassCuts; m++) {
        // Jet mass histograms (used for scale factor determination)
        h_mc_ljet_m1[m] = new TH1F(Form("h_mc_ljet_m1_%s_%s", wp_label, massCutLabels[m]), "ljet_m", 20, 50, 150);
        h_data_ljet_m1[m] = new TH1F(Form("h_data_m1_%s_%s", wp_label, massCutLabels[m]), "ljet_m", 20, 50, 150);
        h_multijets_ljet_m1[m] = new TH1F(Form("h_multijets_ljet_m1_%s_%s", wp_label, massCutLabels[m]), "ljet_m", 20, 50, 150);
        h_mc_ljet_m1[m]->Sumw2(); // Enable error calculation
        h_data_ljet_m1[m]->Sumw2();
        h_multijets_ljet_m1[m]->Sumw2();

        // R_DB histograms (main analysis focus)
        h_mc_R_DB[m] = new TH1F(Form("h_mc_R_DB_%s_%s", wp_label, massCutLabels[m]), "R_DB", 60, 0, 2);
        h_data_R_DB[m] = new TH1F(Form("h_data_R_DB_%s_%s", wp_label, massCutLabels[m]), "R_DB", 60, 0, 2);
        h_multijets_ljet_R_DB[m] = new TH1F(Form("h_multijets_ljet_R_DB_%s_%s", wp_label, massCutLabels[m]), "R_DB", 60, 0, 2);
        h_mc_R_DB[m]->Sumw2(); // Enable error calculation
        h_data_R_DB[m]->Sumw2();
        h_multijets_ljet_R_DB[m]->Sumw2();


        // Set initial histogram styling
        h_data_ljet_m1[m]->SetFillColor(kBlack);
        h_multijets_ljet_m1[m]->SetFillColor(kBlue);
        h_mc_ljet_m1[m]->SetFillColor(kRed);

        h_data_R_DB[m]->SetFillColor(kBlack);
        h_multijets_ljet_R_DB[m]->SetFillColor(kBlue);
        h_mc_R_DB[m]->SetFillColor(kRed);

        // Set line colors to match fill colors for better visibility
        h_data_ljet_m1[m]->SetLineColor(kBlack);
        h_multijets_ljet_m1[m]->SetLineColor(kBlue);
        h_mc_ljet_m1[m]->SetLineColor(kRed);

        h_data_R_DB[m]->SetLineColor(kBlack);
        h_multijets_ljet_R_DB[m]->SetLineColor(kBlue);
        h_mc_R_DB[m]->SetLineColor(kRed);
    }

    //-------------------------------------------
    // Set up variables for reading from trees
    //-------------------------------------------
    std::vector<float>* ljet_m = nullptr;
    double weight, weight_multijets;
    double R_DB;

    // Set branch addresses for all trees
    tree_mc->SetBranchAddress("ljet_m", &ljet_m);
    tree_data->SetBranchAddress("ljet_m", &ljet_m);
    tree_mulijets->SetBranchAddress("ljet_m", &ljet_m);

    tree_mc->SetBranchAddress("R_DB", &R_DB);
    tree_data->SetBranchAddress("R_DB", &R_DB);
    tree_mulijets->SetBranchAddress("R_DB", &R_DB);

    tree_mc->SetBranchAddress("total_weight", &weight);
    tree_mulijets->SetBranchAddress("total_weight", &weight_multijets);

    //-------------------------------------------
    // Fill histograms from trees with different mass cuts
    //-------------------------------------------
    // MC signal Histograms
    Long64_t nentries = tree_mc->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree_mc->GetEntry(i);
        // Ensure ljet_m is not empty before accessing its first element
        if (ljet_m->empty()) continue;
        double jetMass = (*ljet_m)[0]/1000.0;  // Convert MeV to GeV
        
        // Apply different mass cuts
        for (int m = 0; m < nMassCuts; m++) {
            bool passesCut = true;
            
            // Apply lower mass cut if specified
            if (massCutLow[m] > 0 && jetMass < massCutLow[m]) {
                passesCut = false;
            }
            
            // Apply upper mass cut if specified
            if (massCutHigh[m] > 0 && jetMass > massCutHigh[m]) {
                passesCut = false;
            }
            
            // Fill histograms if passes cut
            if (passesCut) {
                h_mc_ljet_m1[m]->Fill(jetMass, weight);
                h_mc_R_DB[m]->Fill(R_DB, weight);
            }
        }
    }

    // Data Histograms
    nentries = tree_data->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree_data->GetEntry(i);
        // Ensure ljet_m is not empty before accessing its first element
        if (ljet_m->empty()) continue;
        double jetMass = (*ljet_m)[0]/1000.0;  // Convert MeV to GeV
        
        // Apply different mass cuts
        for (int m = 0; m < nMassCuts; m++) {
            bool passesCut = true;
            
            // Apply lower mass cut if specified
            if (massCutLow[m] > 0 && jetMass < massCutLow[m]) {
                passesCut = false;
            }
            
            // Apply upper mass cut if specified
            if (massCutHigh[m] > 0 && jetMass > massCutHigh[m]) {
                passesCut = false;
            }
            
            // Fill histograms if passes cut
            if (passesCut) {
                h_data_ljet_m1[m]->Fill(jetMass);
                h_data_R_DB[m]->Fill(R_DB);
            }
        }
    }

    // Multijet background Histograms
    nentries = tree_mulijets->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree_mulijets->GetEntry(i);
        // Ensure ljet_m is not empty before accessing its first element
        if (ljet_m->empty()) continue;
        double jetMass = (*ljet_m)[0]/1000.0;  // Convert MeV to GeV
        
        // Apply different mass cuts
        for (int m = 0; m < nMassCuts; m++) {
            bool passesCut = true;
            
            // Apply lower mass cut if specified
            if (massCutLow[m] > 0 && jetMass < massCutLow[m]) {
                passesCut = false;
            }
            
            // Apply upper mass cut if specified
            if (massCutHigh[m] > 0 && jetMass > massCutHigh[m]) {
                passesCut = false;
            }
            
            // Fill histograms if passes cut
            if (passesCut) {
                h_multijets_ljet_m1[m]->Fill(jetMass, weight_multijets);
                h_multijets_ljet_R_DB[m]->Fill(R_DB, weight_multijets);
            }
        }
    }

    // Now process each mass cut
    for (int m = 0; m < nMassCuts; m++) {
        cout << "==========================================" << endl;
        cout << "Processing mass cut: " << massCutLabels[m] << " for Working Point: " << wp_label << endl;
        cout << "==========================================" << endl;
        
        //-------------------------------------------
        // Calculate scale factor from sidebands
        //-------------------------------------------
        // Define sideband regions
        int binLow1 = h_multijets_ljet_m1[m]->FindBin(50);
        int binHigh1 = h_multijets_ljet_m1[m]->FindBin(65);
        int binLow2 = h_multijets_ljet_m1[m]->FindBin(110);
        int binHigh2 = h_multijets_ljet_m1[m]->FindBin(150);
        
        // Create arrays to store x and y values for the fit
        std::vector<double> x_values;
        std::vector<double> y_values;
        std::vector<double> y_errors;
        
        // Extract bin contents for first sideband region
        for (int bin = binLow1; bin <= binHigh1; ++bin) {
            double bin_center = h_multijets_ljet_m1[m]->GetBinCenter(bin);
            double data_content = h_data_ljet_m1[m]->GetBinContent(bin);
            double mc_content = h_multijets_ljet_m1[m]->GetBinContent(bin);
            
            // Avoid division by zero
            if (mc_content > 0) {
                x_values.push_back(bin_center);
                y_values.push_back(data_content / mc_content);
                
                // Calculate error on the ratio
                double data_error = h_data_ljet_m1[m]->GetBinError(bin);
                double mc_error = h_multijets_ljet_m1[m]->GetBinError(bin);
                double ratio_error = sqrt(pow(data_error/mc_content, 2) + 
                                        pow(data_content*mc_error/(mc_content*mc_content), 2));
                y_errors.push_back(ratio_error);
            }
        }
        
        // Extract bin contents for second sideband region
        for (int bin = binLow2; bin <= binHigh2; ++bin) {
            double bin_center = h_multijets_ljet_m1[m]->GetBinCenter(bin);
            double data_content = h_data_ljet_m1[m]->GetBinContent(bin);
            double mc_content = h_multijets_ljet_m1[m]->GetBinContent(bin);
            
            // Avoid division by zero
            if (mc_content > 0) {
                x_values.push_back(bin_center);
                y_values.push_back(data_content / mc_content);
                
                // Calculate error on the ratio
                double data_error = h_data_ljet_m1[m]->GetBinError(bin);
                double mc_error = h_multijets_ljet_m1[m]->GetBinError(bin);
                double ratio_error = sqrt(pow(data_error/mc_content, 2) + 
                                        pow(data_content*mc_error/(mc_content*mc_content), 2));
                y_errors.push_back(ratio_error);
            }
        }
        
        // Create a TGraphErrors for the fit
        TGraphErrors* ratio_graph = new TGraphErrors(x_values.size());
        for (size_t i = 0; i < x_values.size(); ++i) {
            ratio_graph->SetPoint(i, x_values[i], y_values[i]);
            ratio_graph->SetPointError(i, 0, y_errors[i]);
        }
        
        // Define the fit function (linear)
        TF1* fit_func = new TF1(Form("fit_func_%s_%s", wp_label, massCutLabels[m]), "[0] + [1] * x", 50, 150);
        fit_func->SetParameters(1.0, 0.0);  // Initial guess
        
        // Perform the fit
        ratio_graph->Fit(Form("fit_func_%s_%s", wp_label, massCutLabels[m]), "R");  // "R" restricts fit to the specified range
        
        // Get fit parameters
        double p0 = fit_func->GetParameter(0);  // Intercept
        double p1 = fit_func->GetParameter(1);  // Slope
        double p0_err = fit_func->GetParError(0);
        double p1_err = fit_func->GetParError(1);
        
        printf("Mass cut: %s, Working Point: %s, Linear fit results: SF(m) = %.4f (+- %.4f) + %.8f (+- %.8f) * m\n", 
            massCutLabels[m], wp_label, p0, p0_err, p1, p1_err); 
        
        // Save the fit result to a canvas
        TCanvas* fit_canvas = new TCanvas(Form("fit_canvas_%s_%s", wp_label, massCutLabels[m]), 
                                        Form("Scale Factor Fit (%s - %s)", wp_label, massCutLabels[m]), 800, 600);
        ratio_graph->SetTitle(Form("Data/MC Ratio in Sidebands (%s - %s);Large-R Jet Mass [GeV];Data/MC", wp_label, massCutLabels[m]));
        ratio_graph->SetMarkerStyle(20);
        ratio_graph->SetMarkerColor(kBlue);
        ratio_graph->SetLineColor(kBlue);
        ratio_graph->Draw("AP");
        fit_func->SetLineColor(kRed);
        fit_func->Draw("same");
        
        // Add fit parameters to the canvas
        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(42);
        latex.SetTextSize(0.035);
        latex.DrawLatex(0.2, 0.85, Form("SF(m) = %.4f + %.8f * m", p0, p1));
        latex.DrawLatex(0.2, 0.80, Form("#chi^{2}/ndf = %.2f/%d", fit_func->GetChisquare(), fit_func->GetNDF()));
        
        DrawAtlasLabel(wp_label); // Pass working point label
        
        if (m == 0) {
            latex.DrawLatex(0.2, 0.75, "No Mass Cut");
        } else {
            latex.DrawLatex(0.2, 0.75, Form("Mass Cut: %d < m < %d GeV", (int)massCutLow[m], (int)massCutHigh[m]));
        }
        
        fit_canvas->SaveAs(Form("../Plots/MJB/FitRdbNoMcut/%s/%s/scale_factor_fit.pdf", wp_label, massCutLabels[m]));
        
        //-------------------------------------------
        // Apply scale factor to background histogram
        //-------------------------------------------
        // Make a copy first to preserve the original
        TH1F* h_multijets_ljet_m1_scaled = (TH1F*)h_multijets_ljet_m1[m]->Clone(Form("h_multijets_ljet_m1_scaled_%s_%s", wp_label, massCutLabels[m]));
        
        // Apply mass-dependent scale factor bin by bin
        for (int bin = 1; bin <= h_multijets_ljet_m1_scaled->GetNbinsX(); ++bin) {
            double bin_center = h_multijets_ljet_m1_scaled->GetBinCenter(bin);
            double content = h_multijets_ljet_m1_scaled->GetBinContent(bin);
            double error = h_multijets_ljet_m1_scaled->GetBinError(bin);
            
            // Calculate bin-by-bin scale factor
            double scale_factor = p0 + p1 * bin_center;
            
            // Apply the scale factor
            h_multijets_ljet_m1_scaled->SetBinContent(bin, content * scale_factor);
            h_multijets_ljet_m1_scaled->SetBinError(bin, error * scale_factor);
        }
        
        // Calculate average scale factor for R_DB histogram
        double avgScaleFactor = 0.0;
        int nBins = 0;
        
        // Calculate average SF in the full mass range
        for (int bin = 1; bin <= h_multijets_ljet_m1[m]->GetNbinsX(); ++bin) {
            double bin_center = h_multijets_ljet_m1[m]->GetBinCenter(bin);
            double scale_factor = p0 + p1 * bin_center;
            avgScaleFactor += scale_factor;
            nBins++;
        }
        
        avgScaleFactor /= nBins;
        printf("Mass cut: %s, Working Point: %s, Average Scale Factor: %.4f\n", massCutLabels[m], wp_label, avgScaleFactor);
        
        // Apply average scale factor to R_DB histogram
        h_multijets_ljet_R_DB[m]->Scale(avgScaleFactor);
        
        // Estimate background in the signal region using scaled histogram
        double estimatedBackground = h_multijets_ljet_m1_scaled->Integral(binHigh1, binLow2);
        
        // Calculate signal events
        double dataSignal = h_data_ljet_m1[m]->Integral(binHigh1, binLow2);
        double signalEvents = dataSignal - estimatedBackground;
        
        // Print signal events
        std::cout << "Mass cut: " << massCutLabels[m] << ", Working Point: " << wp_label << ", Number of signal events: " << signalEvents << std::endl;
        
        //-------------------------------------------
        // Create mass stack plot
        //-------------------------------------------
        TCanvas* mass_canvas = new TCanvas(Form("mass_canvas_%s_%s", wp_label, massCutLabels[m]), 
                                        Form("Jet Mass Distribution (%s - %s)", wp_label, massCutLabels[m]), 800, 600);
        
        // Create a stack for background and signal MC
        THStack* mc_stack = new THStack("mc_stack", "Large-R Jet Mass;M [GeV];Events");
        
        // Add scaled background to stack
        h_multijets_ljet_m1_scaled->SetFillColor(kBlue);
        h_multijets_ljet_m1_scaled->SetLineColor(kBlue);
        mc_stack->Add(h_multijets_ljet_m1_scaled);
        
        // Add signal to stack
        h_mc_ljet_m1[m]->SetFillColor(kRed);
        h_mc_ljet_m1[m]->SetLineColor(kRed);
        mc_stack->Add(h_mc_ljet_m1[m]);
        
        // Draw data points
        h_data_ljet_m1[m]->Draw("E");
        mc_stack->Draw("HIST SAME");
        h_data_ljet_m1[m]->Draw("E SAME");
        
        // Create legend
        TLegend* mass_legend = new TLegend(0.6, 0.7, 0.89, 0.89);
        mass_legend->SetBorderSize(0);
        mass_legend->AddEntry(h_data_ljet_m1[m], "Data", "lp");
        mass_legend->AddEntry(h_mc_ljet_m1[m], "Z (#rightarrow bb) + Jets", "f");
        mass_legend->AddEntry(h_multijets_ljet_m1_scaled, "Multijets", "f");
        mass_legend->Draw();
        
        DrawAtlasLabel(wp_label); // Pass working point label
        
        TLatex mass_text;
        mass_text.SetNDC();
        mass_text.SetTextFont(42);
        mass_text.SetTextSize(0.035);
        
        if (m == 0) {
            mass_text.DrawLatex(0.2, 0.75, "No Mass Cut");
        } else {
            mass_text.DrawLatex(0.2, 0.75, Form("Mass Cut: %d < m < %d GeV", (int)massCutLow[m], (int)massCutHigh[m]));
        }
        
        mass_canvas->SaveAs(Form("../Plots/MJB/FitRdbNoMcut/%s/%s/jet_mass_comparison.pdf", wp_label, massCutLabels[m]));

        //-------------------------------------------
        // Create combined MC histogram (signal + background)
        //-------------------------------------------
        h_combined_R_DB[m] = (TH1F*)h_mc_R_DB[m]->Clone(Form("h_combined_R_DB_%s_%s", wp_label, massCutLabels[m]));
        h_combined_R_DB[m]->Add(h_multijets_ljet_R_DB[m]);
        h_combined_R_DB[m]->SetFillColor(kGreen-2);
        h_combined_R_DB[m]->SetLineColor(kGreen-2);
        h_combined_R_DB[m]->SetTitle("Combined MC Signal + Background");

        //-------------------------------------------
        // Fit MC Signal R_DB distribution
        //-------------------------------------------
        TCanvas* mc_canvas = new TCanvas(Form("mc_canvas_%s_%s", wp_label, massCutLabels[m]), 
                                      Form("MC Signal R_DB (%s - %s)", wp_label, massCutLabels[m]), 800, 600);
        mc_canvas->cd();

        TF1 *gausFit_mc = new TF1(Form("gausFit_mc_%s_%s", wp_label, massCutLabels[m]), "gaus", 0.9, 1.1);
        h_mc_R_DB[m]->Fit(gausFit_mc, "R");
        
        double chi_mc = gausFit_mc->GetChisquare();
        int ndf_mc = gausFit_mc->GetNDF();
        
        double amp_mc = gausFit_mc->GetParameter(0);
        double Mean_mc = gausFit_mc->GetParameter(1);
        double Sigma_mc = gausFit_mc->GetParameter(2);
        
        double ampErr_mc = gausFit_mc->GetParError(0);
        double meanErr_mc = gausFit_mc->GetParError(1);
        double sigmaErr_mc = gausFit_mc->GetParError(2);
        
        // Store results for comparison
        mc_means[m] = Mean_mc;
        mc_mean_errs[m] = meanErr_mc;
        mc_sigmas[m] = Sigma_mc;
        mc_sigma_errs[m] = sigmaErr_mc;
        mc_chi2_ndf[m] = (ndf_mc > 0) ? chi_mc / ndf_mc : 0; // Handle ndf = 0 case

        // Add legend and labels
        TLegend *Legend1 = new TLegend(0.6, 0.75, 0.9, 0.9);
        Legend1->SetBorderSize(0);
        Legend1->SetTextSize(0.025);
        Legend1->AddEntry(h_mc_R_DB[m], "Z (#rightarrow bb) + Jets", "l");
        Legend1->AddEntry(gausFit_mc, "Gaussian Fit", "l");
        
        TLatex text;
        text.SetTextSize(0.03);
        text.DrawLatexNDC(0.65, 0.7, Form("#mu = %.3f #pm %.3f", Mean_mc, meanErr_mc));
        text.DrawLatexNDC(0.65, 0.65, Form("#sigma = %.3f #pm %.3f", Sigma_mc, sigmaErr_mc));
        text.DrawLatexNDC(0.65, 0.6, Form("#chi^{2}/NDF = %.3f", mc_chi2_ndf[m]));
        Legend1->Draw();
        
        DrawAtlasLabel(wp_label); // Pass working point label
        
        // Save MC signal R_DB fit plot
        mc_canvas->SaveAs(Form("../Plots/MJB/FitRdbNoMcut/%s/%s/Fit_R_DB_mc.pdf", wp_label, massCutLabels[m]));
        
        // Save MC fit parameters
        SaveFitParameters(Form("../Plots/MJB/FitRdbNoMcut/%s/%s/mc_fit_parameters.txt", wp_label, massCutLabels[m]),
                          Form("MC: Z (#rightarrow bb) + Jets (Mass Cut: %s, WP: %s)", massCutLabels[m], wp_label),
                          Mean_mc, meanErr_mc, 
                          Sigma_mc, sigmaErr_mc,
                          amp_mc, ampErr_mc,
                          chi_mc, ndf_mc);
                          

        // Fit background (Multijet) R_DB distribution
        TCanvas* bkg_canvas = new TCanvas(Form("bkg_canvas_%s_%s", wp_label, massCutLabels[m]), 
                                          Form("Background R_DB Fit (%s - %s)", wp_label, massCutLabels[m]), 800, 600);
        bkg_canvas->cd();
        
        TF1* gausFit_bkg = new TF1(Form("gausFit_bkg_%s_%s", wp_label, massCutLabels[m]), "gaus", 0.9, 1.1);
        h_multijets_ljet_R_DB[m]->Fit(gausFit_bkg, "R");
        
        h_multijets_ljet_R_DB[m]->Draw("HIST");
        gausFit_bkg->SetLineColor(kViolet-4);
        gausFit_bkg->Draw("same");
        
        // Get background fit parameters
        double chi_bkg = gausFit_bkg->GetChisquare();
        int ndf_bkg = gausFit_bkg->GetNDF();
        
        double amp_bkg = gausFit_bkg->GetParameter(0);
        double Mean_bkg = gausFit_bkg->GetParameter(1);
        double Sigma_bkg = gausFit_bkg->GetParameter(2);
        
        double ampErr_bkg = gausFit_bkg->GetParError(0);
        double meanErr_bkg = gausFit_bkg->GetParError(1);
        double sigmaErr_bkg = gausFit_bkg->GetParError(2);

        // Store background results
        bkg_means[m] = Mean_bkg;
        bkg_mean_errs[m] = meanErr_bkg;
        bkg_sigmas[m] = Sigma_bkg;
        bkg_sigma_errs[m] = sigmaErr_bkg;
        bkg_chi2_ndf[m] = (ndf_bkg > 0) ? chi_bkg / ndf_bkg : 0; // Handle ndf = 0 case

        // Add legend and labels for background
        TLegend *Legend_bkg = new TLegend(0.6, 0.75, 0.9, 0.9);
        Legend_bkg->SetBorderSize(0);
        Legend_bkg->SetTextSize(0.025);
        Legend_bkg->AddEntry(h_multijets_ljet_R_DB[m], "Multijets", "l");
        Legend_bkg->AddEntry(gausFit_bkg, "Gaussian Fit", "l");
        
        TLatex text_bkg;
        text_bkg.SetTextSize(0.03);
        text_bkg.DrawLatexNDC(0.65, 0.7, Form("#mu = %.3f #pm %.3f", Mean_bkg, meanErr_bkg));
        text_bkg.DrawLatexNDC(0.65, 0.65, Form("#sigma = %.3f #pm %.3f", Sigma_bkg, sigmaErr_bkg));
        text_bkg.DrawLatexNDC(0.65, 0.6, Form("#chi^{2}/NDF = %.3f", bkg_chi2_ndf[m]));
        Legend_bkg->Draw();
        
        DrawAtlasLabel(wp_label); // Pass working point label
        
        // Save background R_DB fit plot
        bkg_canvas->SaveAs(Form("../Plots/MJB/FitRdbNoMcut/%s/%s/Fit_R_DB_bkg.pdf", wp_label, massCutLabels[m]));
        
        // Save background fit parameters
        SaveFitParameters(Form("../Plots/MJB/FitRdbNoMcut/%s/%s/bkg_fit_parameters.txt", wp_label, massCutLabels[m]),
                          Form("Background: Multijets (Mass Cut: %s, WP: %s)", massCutLabels[m], wp_label),
                          Mean_bkg, meanErr_bkg, 
                          Sigma_bkg, sigmaErr_bkg,
                          amp_bkg, ampErr_bkg,
                          chi_bkg, ndf_bkg);

        // Fit data R_DB distribution
        TCanvas* data_canvas = new TCanvas(Form("data_canvas_%s_%s", wp_label, massCutLabels[m]), 
                                           Form("Data R_DB Fit (%s - %s)", wp_label, massCutLabels[m]), 800, 600);
        data_canvas->cd();
        
        TF1* gausFit_data = new TF1(Form("gausFit_data_%s_%s", wp_label, massCutLabels[m]), "gaus", 0.9, 1.1);
        h_data_R_DB[m]->Fit(gausFit_data, "R");
        
        h_data_R_DB[m]->Draw("EP");
        gausFit_data->SetLineColor(kViolet-4);
        gausFit_data->Draw("same");
        
        // Get data fit parameters
        double chi_data = gausFit_data->GetChisquare();
        int ndf_data = gausFit_data->GetNDF();
        
        double amp_data = gausFit_data->GetParameter(0);
        double Mean_data = gausFit_data->GetParameter(1);
        double Sigma_data = gausFit_data->GetParameter(2);
        
        double ampErr_data = gausFit_data->GetParError(0);
        double meanErr_data = gausFit_data->GetParError(1);
        double sigmaErr_data = gausFit_data->GetParError(2);

        // Store data results
        data_means[m] = Mean_data;
        data_mean_errs[m] = meanErr_data;
        data_sigmas[m] = Sigma_data;
        data_sigma_errs[m] = sigmaErr_data;
        data_chi2_ndf[m] = (ndf_data > 0) ? chi_data / ndf_data : 0; // Handle ndf = 0 case

        // Add legend and labels for data
        TLegend *Legend_data = new TLegend(0.6, 0.75, 0.9, 0.9);
        Legend_data->SetBorderSize(0);
        Legend_data->SetTextSize(0.025);
        Legend_data->AddEntry(h_data_R_DB[m], "Data", "lep");
        Legend_data->AddEntry(gausFit_data, "Gaussian Fit", "l");
        
        TLatex text_data;
        text_data.SetTextSize(0.03);
        text_data.DrawLatexNDC(0.65, 0.7, Form("#mu = %.3f #pm %.3f", Mean_data, meanErr_data));
        text_data.DrawLatexNDC(0.65, 0.65, Form("#sigma = %.3f #pm %.3f", Sigma_data, sigmaErr_data));
        text_data.DrawLatexNDC(0.65, 0.6, Form("#chi^{2}/NDF = %.3f", data_chi2_ndf[m]));
        Legend_data->Draw();
        
        DrawAtlasLabel(wp_label); // Pass working point label
        
        // Save data R_DB fit plot
        data_canvas->SaveAs(Form("../Plots/MJB/FitRdbNoMcut/%s/%s/Fit_R_DB_data.pdf", wp_label, massCutLabels[m]));
        
        // Save data fit parameters
        SaveFitParameters(Form("../Plots/MJB/FitRdbNoMcut/%s/%s/data_fit_parameters.txt", wp_label, massCutLabels[m]),
                          Form("Data (Mass Cut: %s, WP: %s)", massCutLabels[m], wp_label),
                          Mean_data, meanErr_data, 
                          Sigma_data, sigmaErr_data,
                          amp_data, ampErr_data,
                          chi_data, ndf_data);

        // Fit combined MC (signal + background) R_DB distribution
        TCanvas* combined_canvas = new TCanvas(Form("combined_canvas_%s_%s", wp_label, massCutLabels[m]), 
                                               Form("Combined MC R_DB Fit (%s - %s)", wp_label, massCutLabels[m]), 800, 600);
        combined_canvas->cd();
        
        TF1* gausFit_combined = new TF1(Form("gausFit_combined_%s_%s", wp_label, massCutLabels[m]), "gaus", 0.9, 1.1);
        h_combined_R_DB[m]->Fit(gausFit_combined, "R");
        
        h_combined_R_DB[m]->Draw("HIST");
        gausFit_combined->SetLineColor(kViolet-4);
        gausFit_combined->Draw("same");
        
        // Get combined MC fit parameters
        double chi_combined = gausFit_combined->GetChisquare();
        int ndf_combined = gausFit_combined->GetNDF();
        
        double amp_combined = gausFit_combined->GetParameter(0);
        double Mean_combined = gausFit_combined->GetParameter(1);
        double Sigma_combined = gausFit_combined->GetParameter(2);
        
        double ampErr_combined = gausFit_combined->GetParError(0);
        double meanErr_combined = gausFit_combined->GetParError(1);
        double sigmaErr_combined = gausFit_combined->GetParError(2);

        // Store combined MC results
        combined_means[m] = Mean_combined;
        combined_mean_errs[m] = meanErr_combined;
        combined_sigmas[m] = Sigma_combined;
        combined_sigma_errs[m] = sigmaErr_combined;
        combined_chi2_ndf[m] = (ndf_combined > 0) ? chi_combined / ndf_combined : 0; // Handle ndf = 0 case

        // Add legend and labels for combined MC
        TLegend *Legend_combined = new TLegend(0.6, 0.75, 0.9, 0.9);
        Legend_combined->SetBorderSize(0);
        Legend_combined->SetTextSize(0.025);
        Legend_combined->AddEntry(h_combined_R_DB[m], "Combined MC", "l");
        Legend_combined->AddEntry(gausFit_combined, "Gaussian Fit", "l");
        
        TLatex text_combined;
        text_combined.SetTextSize(0.03);
        text_combined.DrawLatexNDC(0.65, 0.7, Form("#mu = %.3f #pm %.3f", Mean_combined, meanErr_combined));
        text_combined.DrawLatexNDC(0.65, 0.65, Form("#sigma = %.3f #pm %.3f", Sigma_combined, sigmaErr_combined));
        text_combined.DrawLatexNDC(0.65, 0.6, Form("#chi^{2}/NDF = %.3f", combined_chi2_ndf[m]));
        Legend_combined->Draw();
        
        DrawAtlasLabel(wp_label); // Pass working point label
        
        // Save combined MC R_DB fit plot
        combined_canvas->SaveAs(Form("../Plots/MJB/FitRdbNoMcut/%s/%s/Fit_R_DB_combined.pdf", wp_label, massCutLabels[m]));
        
        // Save combined MC fit parameters
        SaveFitParameters(Form("../Plots/MJB/FitRdbNoMcut/%s/%s/combined_fit_parameters.txt", wp_label, massCutLabels[m]),
                          Form("Combined MC (Mass Cut: %s, WP: %s)", massCutLabels[m], wp_label),
                          Mean_combined, meanErr_combined, 
                          Sigma_combined, sigmaErr_combined,
                          amp_combined, ampErr_combined,
                          chi_combined, ndf_combined);
        //-------------------------------------------
        // New: Stack plot of all three distributions with Ratio plot Data/MC
        //-------------------------------------------
        TCanvas* combinedstack_canvas = new TCanvas(Form("combinedstack_canvas_masscut%s_%s", massCutLabels[m], wp_label), 
                                                      Form("R_DB Comparison (Mass Bin: %s, WP: %s)", massCutLabels[m], wp_label), 800, 800); // Increased height for better ratio plot
        combinedstack_canvas->cd();
            
        // Configure pads
        TPad* pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
        pad1->SetBottomMargin(0.02);  // Reduced bottom margin
        pad1->SetTopMargin(0.08); // Add some top margin
        pad1->SetLeftMargin(0.12); // Keep consistent with other plots
        pad1->SetRightMargin(0.05); // Keep consistent
        pad1->Draw();
        pad1->cd();  // Switch to pad1
            
        // Clone histograms (same as original logic, adapted for current scope)
        // Note: h_combined_R_DB[m] is already created and represents h_combined_MC_R_DB_ptbin[m][bin] equivalent
        TH1F* h_data_clone_comb = (TH1F*)h_data_R_DB[m]->Clone(Form("h_data_clone_comb_masscut%s_%s", massCutLabels[m], wp_label));
        TH1F* h_mc_clone_comb = (TH1F*)h_mc_R_DB[m]->Clone(Form("h_mc_clone_comb_masscut%s_%s", massCutLabels[m], wp_label));
        TH1F* h_bkg_clone_comb = (TH1F*)h_multijets_ljet_R_DB[m]->Clone(Form("h_bkg_clone_comb_masscut%s_%s", massCutLabels[m], wp_label));
        TH1F* h_combined_clone_comb = (TH1F*)h_combined_R_DB[m]->Clone(Form("h_combined_clone_comb_masscut%s_%s", massCutLabels[m], wp_label));
            
        // Style clones
        h_data_clone_comb->SetMarkerStyle(20);
        h_data_clone_comb->SetMarkerColor(kBlack);
        //h_data_clone_comb->SetLineColor(kBlack);
            
        h_mc_clone_comb->SetFillColor(kRed);
        h_mc_clone_comb->SetLineColor(kRed);
        h_mc_clone_comb->SetFillStyle(1001); // Solid fill
            
        h_bkg_clone_comb->SetFillColor(kBlue);
        h_bkg_clone_comb->SetLineColor(kBlue);
        h_bkg_clone_comb->SetFillStyle(1001); // Solid fill
            
        h_combined_clone_comb->SetFillColor(kGreen+2); // This will be used for the error band if needed, or total MC color
        h_combined_clone_comb->SetLineColor(kGreen+2);
        h_combined_clone_comb->SetFillStyle(1001);
            
        // Create and draw stack
        THStack* stack_comb = new THStack(Form("stack_masscut%s_%s", massCutLabels[m], wp_label), 
                                        Form("R_{DB} Distribution (Mass Bin: %s, WP: %s);;Events", massCutLabels[m], wp_label));  // Note empty x-axis title
        stack_comb->Add(h_bkg_clone_comb);
        stack_comb->Add(h_mc_clone_comb);
        stack_comb->Draw("HIST");
        h_data_clone_comb->Draw("EP SAME");
            
        // Add the fitted curves for data and combined MC
        // Need to recreate the TF1 objects here as they were deleted after previous canvases
        TF1 *gausFit_data_recreate = new TF1(Form("gausFit_data_recreate_masscut%s_%s", massCutLabels[m], wp_label), "gaus", 0.9, 1.1);
        gausFit_data_recreate->SetParameters(amp_data, Mean_data, Sigma_data); // Use stored parameters from this loop iteration
        gausFit_data_recreate->SetLineColor(kViolet-4);
        gausFit_data_recreate->SetLineStyle(1); // line for distinction
        gausFit_data_recreate->SetLineWidth(2);
        gausFit_data_recreate->Draw("same");
            
        TF1 *gausFit_combined_recreate = new TF1(Form("gausFit_combined_recreate_masscut%s_%s", massCutLabels[m], wp_label), "gaus", 0.9, 1.1);
        gausFit_combined_recreate->SetParameters(amp_combined, Mean_combined, Sigma_combined); // Use stored parameters
        gausFit_combined_recreate->SetLineColor(kGreen+3); // A slightly darker green for fit line
        gausFit_combined_recreate->SetLineStyle(1); // line for distinction
        gausFit_combined_recreate->SetLineWidth(2);
        gausFit_combined_recreate->Draw("same");

        // Hide x-axis labels on top pad
        stack_comb->GetXaxis()->SetLabelSize(0);
        stack_comb->GetXaxis()->SetTitleSize(0);
            
        // Style y-axis
        stack_comb->GetYaxis()->SetTitleSize(0.05);
        stack_comb->GetYaxis()->SetTitle("Events");
        stack_comb->GetYaxis()->SetLabelSize(0.04);
            
        stack_comb->SetMaximum(stack_comb->GetMaximum() * 1.5); // Leave room for legend
            
        TH1F* h_mc_total_for_error_band = (TH1F*)stack_comb->GetStack()->Last(); // Get the total MC (signal+background)
            
        TGraphAsymmErrors* err_band = new TGraphAsymmErrors(h_mc_total_for_error_band->GetNbinsX()); // Initialize with number of bins
            
        // Fill the error band points with appropriate errors
        for (int i = 0; i < h_mc_total_for_error_band->GetNbinsX(); ++i) {
            double x_val = h_mc_total_for_error_band->GetBinCenter(i + 1); // Get bin center
            double y_val = h_mc_total_for_error_band->GetBinContent(i + 1); // Get bin content
            double error = h_mc_total_for_error_band->GetBinError(i + 1);
            err_band->SetPoint(i, x_val, y_val); // Set x and y value for the point
            err_band->SetPointEYhigh(i, error);
            err_band->SetPointEYlow(i, error);
        }
            
        // Style the error band
        err_band->SetFillColorAlpha(kBlack, 0.35); // Semi-transparent
        err_band->SetFillStyle(3004);              // Hatch pattern
        err_band->SetLineWidth(0);                 // No outline
        err_band->Draw("2 SAME");                  // "2" draw mode for error band
            
        // Legend (adjusted position)
        TLegend* stack_legend_comb = new TLegend(0.6, 0.65, 0.89, 0.89);
        stack_legend_comb->SetBorderSize(0);          // No border
        stack_legend_comb->SetFillStyle(0);           // Transparent background
        stack_legend_comb->SetTextSize(0.025);        // Readable text size
        stack_legend_comb->AddEntry(h_data_clone_comb, "Data", "lep");
        stack_legend_comb->AddEntry(h_mc_clone_comb, "Z (#rightarrow bb) + Jets", "f");
        stack_legend_comb->AddEntry(h_bkg_clone_comb, "Multijets", "f");
        stack_legend_comb->AddEntry(gausFit_data_recreate, "Data Fit", "l");
        stack_legend_comb->AddEntry(gausFit_combined_recreate, "Total MC Fit", "l");
        stack_legend_comb->AddEntry(err_band, "MC Stat. Unc.", "f");
        stack_legend_comb->Draw();
            
        // Add ATLAS label to the plot
        DrawAtlasLabel(wp_label);
        TLatex text_stack_comb;
        text_stack_comb.SetNDC();
        text_stack_comb.SetTextFont(42);
        text_stack_comb.SetTextSize(0.035);
        text_stack_comb.DrawLatex(0.18, 0.73, Form("Mass Bin: %s", massCutLabels[m]));
            
        pad1->Modified();
        pad1->Update();
            
        // Switch to ratio pad
        combinedstack_canvas->cd();
        TPad* pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
        pad2->SetTopMargin(0.05); // Adjusted for spacing from upper pad
        pad2->SetBottomMargin(0.35); // Increased for x-axis label visibility
        pad2->SetLeftMargin(0.12);
        pad2->SetRightMargin(0.05);
        //pad2->SetGridy(); // Add horizontal grid lines
        pad2->Draw();
        pad2->cd();
            
        // Create ratio histogram
        TH1F* h_ratio = (TH1F*)h_data_clone_comb->Clone(Form("h_ratio_masscut%s_%s", massCutLabels[m], wp_label));
        h_ratio->Divide(h_combined_clone_comb); // Divide by the combined MC clone
            
        // Style ratio plot
        h_ratio->SetTitle(""); // Remove title
        h_ratio->SetMarkerStyle(20); // Changed to 20 for consistency with data points
        h_ratio->SetMarkerSize(1.0);
        h_ratio->SetLineColor(kBlack);
            
        // Dynamic Y-axis for ratio
        double min_ratio = std::numeric_limits<double>::max();
        double max_ratio = std::numeric_limits<double>::min();
        bool has_valid_ratio = false;
        for (int i = 1; i <= h_ratio->GetNbinsX(); ++i) {
            double content = h_ratio->GetBinContent(i);
            if (std::isfinite(content) && content != 0) { // Check for valid content
                if (content < min_ratio) min_ratio = content;
                if (content > max_ratio) max_ratio = content;
                has_valid_ratio = true;
            }
        }
        
        double y_ratio_min = 0.8; // Default lower bound
        double y_ratio_max = 1.2; // Default upper bound

        if (has_valid_ratio) {
            double y_ratio_padding = (max_ratio - min_ratio) * 0.1; // 10% padding
            y_ratio_min = min_ratio - y_ratio_padding;
            y_ratio_max = max_ratio + y_ratio_padding;

            // Constrain dynamic range to typical ratio plot bounds if it goes too wide
            if (y_ratio_min < 0.5) y_ratio_min = 0.5;
            if (y_ratio_max > 1.5) y_ratio_max = 1.5;
        }

        h_ratio->GetYaxis()->SetRangeUser(y_ratio_min, y_ratio_max);
        h_ratio->GetYaxis()->SetTitle("Data/MC"); // MC: ("Sig MC + Multijets)")
        h_ratio->GetYaxis()->SetNdivisions(505); // Fewer divisions for cleaner look
        h_ratio->GetYaxis()->SetTitleSize(0.14);
        h_ratio->GetYaxis()->SetTitleOffset(0.4); // Adjusted for larger title size
        h_ratio->GetYaxis()->SetLabelSize(0.10);
        h_ratio->GetXaxis()->SetTitle("R_{DB}");
        h_ratio->GetXaxis()->SetTitleSize(0.14);
        h_ratio->GetXaxis()->SetTitleOffset(1.0);
        h_ratio->GetXaxis()->SetLabelSize(0.10);
        h_ratio->GetYaxis()->CenterTitle();
            
        // Draw the ratio plot
        h_ratio->Draw("EP"); // Changed from "El same" to "EP" for error bars and points
            
        // Create MC uncertainty band for ratio
        TGraphAsymmErrors* err_band_ratio = new TGraphAsymmErrors(h_combined_clone_comb->GetNbinsX());
        for (int i = 0; i < h_combined_clone_comb->GetNbinsX(); ++i) {
            double x = h_combined_clone_comb->GetBinCenter(i + 1);
            double content = h_combined_clone_comb->GetBinContent(i+1);
            double error = h_combined_clone_comb->GetBinError(i+1);
                
            // Calculate relative error for the ratio
            double rel_error = (content != 0) ? error / content : 0;
                
            err_band_ratio->SetPoint(i, x, 1.0); // Set y to 1.0 for the ratio plot band
            err_band_ratio->SetPointEYhigh(i, rel_error);
            err_band_ratio->SetPointEYlow(i, rel_error);
        }
            
        // Style the ratio error band
        err_band_ratio->SetFillColorAlpha(kBlack, 0.35);
        err_band_ratio->SetFillStyle(3002);  // Consistent style
        err_band_ratio->Draw("2 SAME");
            
        // Add a horizontal reference line at y=1
        double xMin = h_ratio->GetXaxis()->GetXmin();
        double xMax = h_ratio->GetXaxis()->GetXmax();
        TLine* line = new TLine(xMin, 1.0, xMax, 1.0);
        line->SetLineColor(kBlack);
        line->SetLineWidth(2);
        line->SetLineStyle(2);  // Dashed line
        line->Draw("SAME");
            
        // Update both pads to ensure proper display
        pad2->Update();
            
        // Save canvas
        combinedstack_canvas->Update();
        combinedstack_canvas->SaveAs(Form("../Plots/MJB/FitRdbNoMcut/%s/%s/R_DB_stack_ratio_combined.pdf", wp_label, massCutLabels[m]));


        // Clean up individual canvas and objects inside the loop
        /*delete fit_canvas;
        delete ratio_graph;
        delete fit_func;
        delete h_multijets_ljet_m1_scaled; // Delete the cloned histogram
        delete mass_canvas;
        delete mc_stack;
        delete mass_legend;
        delete mc_canvas;
        delete Legend1;
        delete bkg_canvas;
        delete Legend_bkg;
        delete data_canvas;
        delete Legend_data;
        delete combined_canvas;
        delete Legend_combined;*/
    } // End of loop over mass cuts

    // Create comparison plots
    // 1. Mean Comparison
    TCanvas* mean_canvas = new TCanvas(Form("mean_canvas_%s", wp_label), Form("R_DB Mean vs Mass Cut (%s)", wp_label), 800, 600);
    mean_canvas->cd();
    
    double massCutX[nMassCuts];
    for (int m = 0; m < nMassCuts; m++) {
        massCutX[m] = m; // Use indices as x-coordinates
    }
    
    // Create graphs for mean comparison
    TGraphErrors* g_mc_mean = new TGraphErrors(nMassCuts, massCutX, mc_means, nullptr, mc_mean_errs);
    TGraphErrors* g_data_mean = new TGraphErrors(nMassCuts, massCutX, data_means, nullptr, data_mean_errs);
    TGraphErrors* g_bkg_mean = new TGraphErrors(nMassCuts, massCutX, bkg_means, nullptr, bkg_mean_errs);
    TGraphErrors* g_combined_mean = new TGraphErrors(nMassCuts, massCutX, combined_means, nullptr, combined_mean_errs);
    
    // Style graphs
    g_mc_mean->SetMarkerStyle(20);
    g_mc_mean->SetMarkerColor(kRed);
    g_mc_mean->SetLineColor(kRed);
    
    g_data_mean->SetMarkerStyle(21);
    g_data_mean->SetMarkerColor(kBlack);
    g_data_mean->SetLineColor(kBlack);
    
    g_bkg_mean->SetMarkerStyle(22);
    g_bkg_mean->SetMarkerColor(kBlue);
    g_bkg_mean->SetLineColor(kBlue);
    
    g_combined_mean->SetMarkerStyle(23);
    g_combined_mean->SetMarkerColor(kGreen+2);
    g_combined_mean->SetLineColor(kGreen+2);
    
    // Draw graphs
    g_mc_mean->SetTitle(Form("R_{MJB} Mean vs Mass Cut (%s);Mass Cut;Mean R_{MJB}", wp_label));
    g_mc_mean->GetYaxis()->SetRangeUser(0.95, 1.05);
    g_mc_mean->Draw("APL");
    g_data_mean->Draw("PL SAME");
    g_bkg_mean->Draw("PL SAME");
    g_combined_mean->Draw("PL SAME");
    
    // Add custom x-axis labels
    TAxis *xaxis = g_mc_mean->GetXaxis();
    for (int m = 0; m < nMassCuts; m++) {
        xaxis->SetBinLabel(xaxis->FindBin(m), massCutLabels[m]);
    }
    
    // Add legend
    TLegend* mean_legend = new TLegend(0.6, 0.7, 0.89, 0.89);
    mean_legend->SetBorderSize(0);
    mean_legend->AddEntry(g_mc_mean, "MC Signal", "pl");
    mean_legend->AddEntry(g_data_mean, "Data", "pl");
    mean_legend->AddEntry(g_bkg_mean, "Multijets", "pl");
    mean_legend->AddEntry(g_combined_mean, "Combined MC", "pl");
    mean_legend->Draw();
    
    DrawAtlasLabel(wp_label); // Pass working point label
    mean_canvas->SaveAs(Form("../Plots/MJB/FitRdbNoMcut/%s/R_DB_mean_vs_masscut.pdf", wp_label));
    /*delete mean_canvas;
    delete g_mc_mean;
    delete g_data_mean;
    delete g_bkg_mean;
    delete g_combined_mean;
    delete mean_legend;*/


    // 2. Sigma Comparison
    TCanvas* sigma_canvas = new TCanvas(Form("sigma_canvas_%s", wp_label), Form("R_DB Sigma vs Mass Cut (%s)", wp_label), 800, 600);
    sigma_canvas->cd();
    
    // Create graphs for sigma comparison
    TGraphErrors* g_mc_sigma = new TGraphErrors(nMassCuts, massCutX, mc_sigmas, nullptr, mc_sigma_errs);
    TGraphErrors* g_data_sigma = new TGraphErrors(nMassCuts, massCutX, data_sigmas, nullptr, data_sigma_errs);
    TGraphErrors* g_bkg_sigma = new TGraphErrors(nMassCuts, massCutX, bkg_sigmas, nullptr, bkg_sigma_errs);
    TGraphErrors* g_combined_sigma = new TGraphErrors(nMassCuts, massCutX, combined_sigmas, nullptr, combined_sigma_errs);
    
    // Style sigma graphs
    g_mc_sigma->SetMarkerStyle(20);
    g_mc_sigma->SetMarkerColor(kRed);
    g_mc_sigma->SetLineColor(kRed);
    
    g_data_sigma->SetMarkerStyle(21);
    g_data_sigma->SetMarkerColor(kBlack);
    g_data_sigma->SetLineColor(kBlack);
    
    g_bkg_sigma->SetMarkerStyle(22);
    g_bkg_sigma->SetMarkerColor(kBlue);
    g_bkg_sigma->SetLineColor(kBlue);
    
    g_combined_sigma->SetMarkerStyle(23);
    g_combined_sigma->SetMarkerColor(kGreen+2);
    g_combined_sigma->SetLineColor(kGreen+2);
    
    // Draw sigma graphs
    g_mc_sigma->SetTitle(Form("R_{MJB} Sigma vs Mass Cut (%s);Mass Cut;#sigma R_{MJB}", wp_label));
    g_mc_sigma->GetYaxis()->SetRangeUser(0.0, 0.2);
    g_mc_sigma->Draw("APL");
    g_data_sigma->Draw("PL SAME");
    g_bkg_sigma->Draw("PL SAME");
    g_combined_sigma->Draw("PL SAME");
    
    // Add custom x-axis labels for sigma plot
    TAxis *xaxis_sigma = g_mc_sigma->GetXaxis();
    for (int m = 0; m < nMassCuts; m++) {
        xaxis_sigma->SetBinLabel(xaxis_sigma->FindBin(m), massCutLabels[m]);
    }
    
    // Add legend for sigma plot
    TLegend* sigma_legend = new TLegend(0.6, 0.7, 0.89, 0.89);
    sigma_legend->SetBorderSize(0);
    sigma_legend->AddEntry(g_mc_sigma, "MC Signal", "pl");
    sigma_legend->AddEntry(g_data_sigma, "Data", "pl");
    sigma_legend->AddEntry(g_bkg_sigma, "Multijets", "pl");
    sigma_legend->AddEntry(g_combined_sigma, "Combined MC", "pl");
    sigma_legend->Draw();
    
    DrawAtlasLabel(wp_label); // Pass working point label
    sigma_canvas->SaveAs(Form("../Plots/MJB/FitRdbNoMcut/%s/R_DB_sigma_vs_masscut.pdf", wp_label));
    /*delete sigma_canvas;
    delete g_mc_sigma;
    delete g_data_sigma;
    delete g_bkg_sigma;
    delete g_combined_sigma;
    delete sigma_legend;*/


    // 3. Ratio Comparison (Means)
    // Calculate ratios
    double mean_data_mc_ratio[nMassCuts], mean_data_mc_ratio_err[nMassCuts];
    double mean_data_combined_ratio[nMassCuts], mean_data_combined_ratio_err[nMassCuts];
    double mean_bkg_mc_ratio[nMassCuts], mean_bkg_mc_ratio_err[nMassCuts];
    
    for (int m = 0; m < nMassCuts; m++) {
        // Data/MC ratio for means
        if (mc_means[m] != 0) {
            mean_data_mc_ratio[m] = data_means[m] / mc_means[m];
            mean_data_mc_ratio_err[m] = mean_data_mc_ratio[m] * sqrt(
                pow(data_mean_errs[m]/data_means[m], 2) + 
                pow(mc_mean_errs[m]/mc_means[m], 2));
        } else {
            mean_data_mc_ratio[m] = 0.0;
            mean_data_mc_ratio_err[m] = 0.0;
        }
            
        // Data/Combined ratio for means
        if (combined_means[m] != 0) {
            mean_data_combined_ratio[m] = data_means[m] / combined_means[m];
            mean_data_combined_ratio_err[m] = mean_data_combined_ratio[m] * sqrt(
                pow(data_mean_errs[m]/data_means[m], 2) + 
                pow(combined_mean_errs[m]/combined_means[m], 2));
        } else {
            mean_data_combined_ratio[m] = 0.0;
            mean_data_combined_ratio_err[m] = 0.0;
        }
            
        // Background/MC ratio for means
        if (mc_means[m] != 0) {
            mean_bkg_mc_ratio[m] = bkg_means[m] / mc_means[m];
            mean_bkg_mc_ratio_err[m] = mean_bkg_mc_ratio[m] * sqrt(
                pow(bkg_mean_errs[m]/bkg_means[m], 2) + 
                pow(mc_mean_errs[m]/mc_means[m], 2));
        } else {
            mean_bkg_mc_ratio[m] = 0.0;
            mean_bkg_mc_ratio_err[m] = 0.0;
        }
    }
    
    // Create ratio comparison plot for Means
    TCanvas* mean_ratio_canvas = new TCanvas(Form("mean_ratio_canvas_%s", wp_label), Form("R_DB Mean Ratio vs Mass Cut (%s)", wp_label), 800, 600);
    mean_ratio_canvas->cd();
    
    // Create graphs for mean ratio comparison
    TGraphErrors* g_mean_data_mc_ratio = new TGraphErrors(nMassCuts, massCutX, mean_data_mc_ratio, nullptr, mean_data_mc_ratio_err);
    TGraphErrors* g_mean_data_combined_ratio = new TGraphErrors(nMassCuts, massCutX, mean_data_combined_ratio, nullptr, mean_data_combined_ratio_err);
    TGraphErrors* g_mean_bkg_mc_ratio = new TGraphErrors(nMassCuts, massCutX, mean_bkg_mc_ratio, nullptr, mean_bkg_mc_ratio_err);
    
    // Style mean ratio graphs
    g_mean_data_mc_ratio->SetMarkerStyle(20);
    g_mean_data_mc_ratio->SetMarkerColor(kRed);
    g_mean_data_mc_ratio->SetLineColor(kRed);
    
    g_mean_data_combined_ratio->SetMarkerStyle(21);
    g_mean_data_combined_ratio->SetMarkerColor(kBlack);
    g_mean_data_combined_ratio->SetLineColor(kBlack);
    
    g_mean_bkg_mc_ratio->SetMarkerStyle(22);
    g_mean_bkg_mc_ratio->SetMarkerColor(kBlue);
    g_mean_bkg_mc_ratio->SetLineColor(kBlue);
    
    // Draw mean ratio graphs
    g_mean_data_mc_ratio->SetTitle(Form("R_{MJB} Mean Ratio vs Mass Cut (%s);Mass Cut;Mean Ratio", wp_label));
    g_mean_data_mc_ratio->GetYaxis()->SetRangeUser(0.95, 1.05);
    g_mean_data_mc_ratio->Draw("APL");
    g_mean_data_combined_ratio->Draw("PL SAME");
    g_mean_bkg_mc_ratio->Draw("PL SAME");
    
    // Add reference line at 1.0
    TLine* mean_ratio_line = new TLine(-0.5, 1.0, nMassCuts-0.5, 1.0);
    mean_ratio_line->SetLineStyle(2);
    mean_ratio_line->SetLineColor(kRed);
    mean_ratio_line->Draw();
    
    // Add custom x-axis labels for mean ratio plot
    TAxis *xaxis_mean_ratio = g_mean_data_mc_ratio->GetXaxis();
    for (int m = 0; m < nMassCuts; m++) {
        xaxis_mean_ratio->SetBinLabel(xaxis_mean_ratio->FindBin(m), massCutLabels[m]);
    }
    
    // Add legend for mean ratio plot
    TLegend* mean_ratio_legend = new TLegend(0.6, 0.7, 0.89, 0.89);
    mean_ratio_legend->SetBorderSize(0);
    mean_ratio_legend->AddEntry(g_mean_data_mc_ratio, "Data/MC Signal", "pl");
    mean_ratio_legend->AddEntry(g_mean_data_combined_ratio, "Data/Combined MC", "pl");
    mean_ratio_legend->AddEntry(g_mean_bkg_mc_ratio, "Multijets/MC Signal", "pl");
    mean_ratio_legend->AddEntry(mean_ratio_line, "Ratio = 1", "l");
    mean_ratio_legend->Draw();
    
    DrawAtlasLabel(wp_label); // Pass working point label
    mean_ratio_canvas->SaveAs(Form("../Plots/MJB/FitRdbNoMcut/%s/R_DB_mean_ratio_vs_masscut.pdf", wp_label));
    /*delete mean_ratio_canvas;
    delete g_mean_data_mc_ratio;
    delete g_mean_data_combined_ratio;
    delete g_mean_bkg_mc_ratio;
    delete mean_ratio_line;
    delete mean_ratio_legend;*/

    // 4. Ratio Comparison (Sigmas)
    double sigma_data_mc_ratio[nMassCuts], sigma_data_mc_ratio_err[nMassCuts];
    double sigma_data_combined_ratio[nMassCuts], sigma_data_combined_ratio_err[nMassCuts];
    double sigma_bkg_mc_ratio[nMassCuts], sigma_bkg_mc_ratio_err[nMassCuts];

    for (int m = 0; m < nMassCuts; m++) {
        // Data/MC sigma ratio
        if (mc_sigmas[m] != 0) {
            sigma_data_mc_ratio[m] = data_sigmas[m] / mc_sigmas[m];
            sigma_data_mc_ratio_err[m] = sigma_data_mc_ratio[m] * sqrt(
                pow(data_sigma_errs[m]/data_sigmas[m], 2) +
                pow(mc_sigma_errs[m]/mc_sigmas[m], 2));
        } else {
            sigma_data_mc_ratio[m] = 0.0;
            sigma_data_mc_ratio_err[m] = 0.0;
        }

        // Data/Combined sigma ratio
        if (combined_sigmas[m] != 0) {
            sigma_data_combined_ratio[m] = data_sigmas[m] / combined_sigmas[m];
            sigma_data_combined_ratio_err[m] = sigma_data_combined_ratio[m] * sqrt(
                pow(data_sigma_errs[m]/data_sigmas[m], 2) +
                pow(combined_sigma_errs[m]/combined_sigmas[m], 2));
        } else {
            sigma_data_combined_ratio[m] = 0.0;
            sigma_data_combined_ratio_err[m] = 0.0;
        }

        // Background/MC sigma ratio
        if (mc_sigmas[m] != 0) {
            sigma_bkg_mc_ratio[m] = bkg_sigmas[m] / mc_sigmas[m];
            sigma_bkg_mc_ratio_err[m] = sigma_bkg_mc_ratio[m] * sqrt(
                pow(bkg_sigma_errs[m]/bkg_sigmas[m], 2) +
                pow(mc_sigma_errs[m]/mc_sigmas[m], 2));
        } else {
            sigma_bkg_mc_ratio[m] = 0.0;
            sigma_bkg_mc_ratio_err[m] = 0.0;
        }
    }

    // Create ratio comparison plot for Sigmas
    TCanvas* sigma_ratio_canvas = new TCanvas(Form("sigma_ratio_canvas_%s", wp_label), Form("R_DB Sigma Ratio vs Mass Cut (%s)", wp_label), 800, 600);
    sigma_ratio_canvas->cd();

    // Create graphs for sigma ratio comparison
    TGraphErrors* g_sigma_data_mc_ratio = new TGraphErrors(nMassCuts, massCutX, sigma_data_mc_ratio, nullptr, sigma_data_mc_ratio_err);
    TGraphErrors* g_sigma_data_combined_ratio = new TGraphErrors(nMassCuts, massCutX, sigma_data_combined_ratio, nullptr, sigma_data_combined_ratio_err);
    TGraphErrors* g_sigma_bkg_mc_ratio = new TGraphErrors(nMassCuts, massCutX, sigma_bkg_mc_ratio, nullptr, sigma_bkg_mc_ratio_err);

    // Style sigma ratio graphs
    g_sigma_data_mc_ratio->SetMarkerStyle(20);
    g_sigma_data_mc_ratio->SetMarkerColor(kRed);
    g_sigma_data_mc_ratio->SetLineColor(kRed);

    g_sigma_data_combined_ratio->SetMarkerStyle(21);
    g_sigma_data_combined_ratio->SetMarkerColor(kBlack);
    g_sigma_data_combined_ratio->SetLineColor(kBlack);

    g_sigma_bkg_mc_ratio->SetMarkerStyle(22);
    g_sigma_bkg_mc_ratio->SetMarkerColor(kBlue);
    g_sigma_bkg_mc_ratio->SetLineColor(kBlue);

    // Draw sigma ratio graphs
    g_sigma_data_mc_ratio->SetTitle(Form("R_{MJB} Sigma Ratio vs Mass Cut (%s);Mass Cut;#sigma Ratio", wp_label));
    g_sigma_data_mc_ratio->GetYaxis()->SetRangeUser(0.95, 1.05);
    g_sigma_data_mc_ratio->Draw("APL");
    g_sigma_data_combined_ratio->Draw("PL SAME");
    g_sigma_bkg_mc_ratio->Draw("PL SAME");

    // Add reference line at 1.0
    TLine* sigma_ratio_line = new TLine(-0.5, 1.0, nMassCuts-0.5, 1.0);
    sigma_ratio_line->SetLineStyle(2);
    sigma_ratio_line->SetLineColor(kRed);
    sigma_ratio_line->Draw();

    // Add custom x-axis labels for sigma ratio plot
    TAxis *xaxis_sigma_ratio = g_sigma_data_mc_ratio->GetXaxis();
    for (int m = 0; m < nMassCuts; m++) {
        xaxis_sigma_ratio->SetBinLabel(xaxis_sigma_ratio->FindBin(m), massCutLabels[m]);
    }

    // Add legend for sigma ratio plot
    TLegend* sigma_ratio_legend = new TLegend(0.6, 0.7, 0.89, 0.89);
    sigma_ratio_legend->SetBorderSize(0);
    sigma_ratio_legend->AddEntry(g_sigma_data_mc_ratio, "Data/MC Signal", "pl");
    sigma_ratio_legend->AddEntry(g_sigma_data_combined_ratio, "Data/Combined MC", "pl");
    sigma_ratio_legend->AddEntry(g_sigma_bkg_mc_ratio, "Multijets/MC Signal", "pl");
    sigma_ratio_legend->AddEntry(sigma_ratio_line, "Ratio = 1", "l");
    sigma_ratio_legend->Draw();

    DrawAtlasLabel(wp_label); // Pass working point label
    sigma_ratio_canvas->SaveAs(Form("../Plots/MJB/FitRdbNoMcut/%s/R_DB_sigma_ratio_vs_masscut.pdf", wp_label));
    /*delete sigma_ratio_canvas;
    delete g_sigma_data_mc_ratio;
    delete g_sigma_data_combined_ratio;
    delete g_sigma_bkg_mc_ratio;
    delete sigma_ratio_line;
    delete sigma_ratio_legend;
    */


    // Fit constant to Data/Combined MC Mean Ratio graph
    TF1 *fit_const_mean_ratio = new TF1(Form("fit_const_mean_ratio_%s", wp_label), "[0]", -0.5, nMassCuts-0.5);
    g_mean_data_combined_ratio->Fit(fit_const_mean_ratio, "RQ"); // "R" for range, "Q" for quiet

    double final_data_combined_mean_ratio_fit_val = fit_const_mean_ratio->GetParameter(0);
    double final_data_combined_mean_ratio_fit_err = fit_const_mean_ratio->GetParError(0);
    double final_data_combined_mean_ratio_chi2 = fit_const_mean_ratio->GetChisquare();
    int final_data_combined_mean_ratio_ndf = fit_const_mean_ratio->GetNDF();

    // Fit constant to Data/Combined MC Sigma Ratio graph
    TF1 *fit_const_sigma_ratio = new TF1(Form("fit_const_sigma_ratio_%s", wp_label), "[0]", -0.5, nMassCuts-0.5);
    g_sigma_data_combined_ratio->Fit(fit_const_sigma_ratio, "RQ"); // "R" for range, "Q" for quiet

    double final_data_combined_sigma_ratio_fit_val = fit_const_sigma_ratio->GetParameter(0);
    double final_data_combined_sigma_ratio_fit_err = fit_const_sigma_ratio->GetParError(0);
    double final_data_combined_sigma_ratio_chi2 = fit_const_sigma_ratio->GetChisquare();
    int final_data_combined_sigma_ratio_ndf = fit_const_sigma_ratio->GetNDF();


    // Create summary table
    ofstream summaryTable(Form("../Plots/MJB/FitRdbNoMcut/%s/R_DB_summary_table.txt", wp_label));
    if (summaryTable.is_open()) {
        summaryTable << "R_DB SUMMARY TABLE FOR ALL MASS CUTS (Working Point: " << wp_label << ")" << endl;
        summaryTable << "====================================" << endl << endl;
        
        for (int m = 0; m < nMassCuts; m++) {
            summaryTable << "Mass Cut: " << massCutLabels[m] << endl;
            summaryTable << "-----------------------------" << endl;
            
            summaryTable << left 
                << setw(20) << "Sample" 
                << setw(20) << "Mean" 
                << setw(20) << "Mean Error" 
                << setw(20) << "Sigma" 
                << setw(20) << "Sigma Error"
                << setw(20) << "Chi2/NDF"
                << endl;
            
            summaryTable << "-----------------------------" << endl;
            
            // MC Signal
            summaryTable << left 
                << setw(20) << "MC Signal" 
                << setw(20) << Form("%.4f", mc_means[m]) 
                << setw(20) << Form("%.4f", mc_mean_errs[m]) 
                << setw(20) << Form("%.4f", mc_sigmas[m]) 
                << setw(20) << Form("%.4f", mc_sigma_errs[m]) 
                << setw(20) << Form("%.3f", mc_chi2_ndf[m]) 
                << endl;
            
            // Data
            summaryTable << left 
                << setw(20) << "Data" 
                << setw(20) << Form("%.4f", data_means[m]) 
                << setw(20) << Form("%.4f", data_mean_errs[m]) 
                << setw(20) << Form("%.4f", data_sigmas[m]) 
                << setw(20) << Form("%.4f", data_sigma_errs[m]) 
                << setw(20) << Form("%.3f", data_chi2_ndf[m]) 
                << endl;
            
            // Multijets
            summaryTable << left 
                << setw(20) << "Multijets" 
                << setw(20) << Form("%.4f", bkg_means[m]) 
                << setw(20) << Form("%.4f", bkg_mean_errs[m]) 
                << setw(20) << Form("%.4f", bkg_sigmas[m]) 
                << setw(20) << Form("%.4f", bkg_sigma_errs[m]) 
                << setw(20) << Form("%.3f", bkg_chi2_ndf[m]) 
                << endl;
            
            // Combined MC
            summaryTable << left 
                << setw(20) << "Combined MC" 
                << setw(20) << Form("%.4f", combined_means[m]) 
                << setw(20) << Form("%.4f", combined_mean_errs[m]) 
                << setw(20) << Form("%.4f", combined_sigmas[m]) 
                << setw(20) << Form("%.4f", combined_sigma_errs[m]) 
                << setw(20) << Form("%.3f", combined_chi2_ndf[m]) 
                << endl;
            
            summaryTable << endl;

            // Ratios for Means
            summaryTable << "----------------------------------------------------------------------------------------------------------------------------------------\n";
            summaryTable << "Ratios of Means\n";
            summaryTable << "----------------------------------------------------------------------------------------------------------------------------------------\n";
            summaryTable << std::left 
                         << std::setw(25) << "Ratio Type" 
                         << std::setw(20) << "Value" 
                         << std::setw(20) << "Error" 
                         << std::setw(20) << "Relative Impr." << std::endl;
            summaryTable << "----------------------------------------------------------------------------------------------------------------------------------------\n";
            
            double relative_improvement_mean_data_mc = (mean_data_mc_ratio[m] != 0) ? (1.0 - mean_data_mc_ratio[m]) * 100.0 : 0.0;
            double relative_improvement_mean_data_combined = (mean_data_combined_ratio[m] != 0) ? (1.0 - mean_data_combined_ratio[m]) * 100.0 : 0.0;
            double relative_improvement_mean_bkg_mc = (mean_bkg_mc_ratio[m] != 0) ? (1.0 - mean_bkg_mc_ratio[m]) * 100.0 : 0.0;

            summaryTable << std::left
                         << std::setw(25) << "Data/MC Signal Mean"
                         << std::setw(20) << Form("%.4f", mean_data_mc_ratio[m])
                         << std::setw(20) << Form("%.4f", mean_data_mc_ratio_err[m])
                         << std::setw(20) << Form("%.2f %%", relative_improvement_mean_data_mc) << std::endl;

            summaryTable << std::left
                         << std::setw(25) << "Data/Combined MC Mean"
                         << std::setw(20) << Form("%.4f", mean_data_combined_ratio[m])
                         << std::setw(20) << Form("%.4f", mean_data_combined_ratio_err[m])
                         << std::setw(20) << Form("%.2f %%", relative_improvement_mean_data_combined) << std::endl;

            summaryTable << std::left
                         << std::setw(25) << "Multijets/MC Signal Mean"
                         << std::setw(20) << Form("%.4f", mean_bkg_mc_ratio[m])
                         << std::setw(20) << Form("%.4f", mean_bkg_mc_ratio_err[m])
                         << std::setw(20) << Form("%.2f %%", relative_improvement_mean_bkg_mc) << std::endl;
            summaryTable << std::endl;


            // Ratios for Sigmas
            summaryTable << "----------------------------------------------------------------------------------------------------------------------------------------\n";
            summaryTable << "Ratios of Sigmas (Resolutions)\n";
            summaryTable << "----------------------------------------------------------------------------------------------------------------------------------------\n";
            summaryTable << std::left 
                         << std::setw(25) << "Ratio Type" 
                         << std::setw(20) << "Value" 
                         << std::setw(20) << "Error" 
                         << std::setw(20) << "Relative Impr." << std::endl;
            summaryTable << "----------------------------------------------------------------------------------------------------------------------------------------\n";

            double relative_improvement_sigma_data_mc = (sigma_data_mc_ratio[m] != 0) ? (1.0 - sigma_data_mc_ratio[m]) * 100.0 : 0.0;
            double relative_improvement_sigma_data_combined = (sigma_data_combined_ratio[m] != 0) ? (1.0 - sigma_data_combined_ratio[m]) * 100.0 : 0.0;
            double relative_improvement_sigma_bkg_mc = (sigma_bkg_mc_ratio[m] != 0) ? (1.0 - sigma_bkg_mc_ratio[m]) * 100.0 : 0.0;
            
            summaryTable << std::left
                         << std::setw(25) << "Data/MC Signal Sigma"
                         << std::setw(20) << Form("%.4f", sigma_data_mc_ratio[m])
                         << std::setw(20) << Form("%.4f", sigma_data_mc_ratio_err[m])
                         << std::setw(20) << Form("%.2f %%", relative_improvement_sigma_data_mc) << std::endl;

            summaryTable << std::left
                         << std::setw(25) << "Data/Combined MC Sigma"
                         << std::setw(20) << Form("%.4f", sigma_data_combined_ratio[m])
                         << std::setw(20) << Form("%.4f", sigma_data_combined_ratio_err[m])
                         << std::setw(20) << Form("%.2f %%", relative_improvement_sigma_data_combined) << std::endl;

            summaryTable << std::left
                         << std::setw(25) << "Multijets/MC Signal Sigma"
                         << std::setw(20) << Form("%.4f", sigma_bkg_mc_ratio[m])
                         << std::setw(20) << Form("%.4f", sigma_bkg_mc_ratio_err[m])
                         << std::setw(20) << Form("%.2f %%", relative_improvement_sigma_bkg_mc) << std::endl;
            summaryTable << std::endl;

            // Overall Data/Combined MC Ratio Statistics (from constant fits)
            summaryTable << "========================================================================================================================\n";
            summaryTable << "Overall Data/Combined MC Ratio Statistics (from constant fits across mass cuts)\n";
            summaryTable << "========================================================================================================================\n";
            summaryTable << std::left 
                         << std::setw(35) << "Metric" 
                         << std::setw(20) << "Fitted Value" 
                         << std::setw(20) << "Fitted Error" 
                         << std::setw(20) << "Chi2/NDF" 
                         << std::setw(20) << "Relative Impr." << std::endl;
            summaryTable << "------------------------------------------------------------------------------------------------------------------------\n";
            
            double overall_relative_improvement_mean = (final_data_combined_mean_ratio_fit_val != 0) ? (1.0 - final_data_combined_mean_ratio_fit_val) * 100.0 : 0.0;
            double overall_relative_improvement_sigma = (final_data_combined_sigma_ratio_fit_val != 0) ? (1.0 - final_data_combined_sigma_ratio_fit_val) * 100.0 : 0.0;

            summaryTable << std::left
                         << std::setw(35) << "Mean Ratio (Data/Combined MC)"
                         << std::setw(20) << Form("%.4f", final_data_combined_mean_ratio_fit_val)
                         << std::setw(20) << Form("%.4f", final_data_combined_mean_ratio_fit_err)
                         << std::setw(20) << Form("%.3f", (final_data_combined_mean_ratio_ndf > 0) ? final_data_combined_mean_ratio_chi2 / final_data_combined_mean_ratio_ndf : 0.0)
                         << std::setw(20) << Form("%.2f %%", overall_relative_improvement_mean) << std::endl;
            
            summaryTable << std::left
                         << std::setw(35) << "Sigma Ratio (Data/Combined MC)"
                         << std::setw(20) << Form("%.4f", final_data_combined_sigma_ratio_fit_val)
                         << std::setw(20) << Form("%.4f", final_data_combined_sigma_ratio_fit_err)
                         << std::setw(20) << Form("%.3f", (final_data_combined_sigma_ratio_ndf > 0) ? final_data_combined_sigma_ratio_chi2 / final_data_combined_sigma_ratio_ndf : 0.0)
                         << std::setw(20) << Form("%.2f %%", overall_relative_improvement_sigma) << std::endl;
            
            summaryTable << std::endl << std::endl;
        }
        
        summaryTable.close();
        std::cout << "Summary table saved to ../Plots/MJB/FitRdbNoMcut/" << wp_label << "/R_DB_summary_table.txt" << std::endl;
    } else {
        std::cerr << "Error: Unable to open summary table file for writing for WP " << wp_label << "." << std::endl;
    }

    // Clean up histograms (moved outside the mass cut loop)
    /*for (int m = 0; m < nMassCuts; m++) {
        delete h_mc_ljet_m1[m];
        delete h_data_ljet_m1[m];
        delete h_multijets_ljet_m1[m];
        delete h_mc_R_DB[m];
        delete h_data_R_DB[m];
        delete h_multijets_ljet_R_DB[m];
        delete h_combined_R_DB[m];
    }*/
    
    // Clean up TGraphs from the end of the script
    //delete fit_const_mean_ratio;
    //delete fit_const_sigma_ratio;


    // Close input files
    if (fmc) { fmc->Close();}// delete fmc; }
    if (fdata) { fdata->Close();} //delete fdata; }
    if (fmulijets) { fmulijets->Close();} //delete fmulijets; }

    cout << "R_DB analysis with mass cuts completed successfully for Working Point: " << wp_label << endl;
}

// Function to run the full analysis for different working points
void RdbNoMcut() {
    const int nWorkingPoints = 3;
    const char* workingPoints[nWorkingPoints] = {"25", "46", "74"};

    // Disable ROOT global graphics for batch processing
    gROOT->SetBatch(kTRUE);

    for (int i = 0; i < nWorkingPoints; ++i) {
        RunRdbAnalysis(workingPoints[i]);
    }
}
