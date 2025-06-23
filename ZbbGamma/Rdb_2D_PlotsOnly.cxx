#include "../../Util/AtlasStyle.C"
#include <iostream>
#include <string>
#include <vector>
#include <iomanip> // For std::setw
#include <fstream> // For file output
#include <cmath>   // For std::pow and std::sqrt

// Forward declarations for ROOT classes
class TFile;
class TTree;
class TH1F;
class TH2F;
class TCanvas;
class TGraphErrors;
class TF1;
class TLegend;
class TLine;
class TLatex;
class TPaletteAxis;

// Global constants for pT and mass bins (must be defined globally or passed)
const int nPtBins = 2;
const double ptBins[nPtBins + 1] = {150, 300, 450}; // 450-600, 600-900, 900-1200
const char* ptBinLabels[nPtBins] = {"150 < p_{T} < 300 GeV", "300 < p_{T} < 450 GeV"};

/*const int nMassCuts = 6;
const char* massCutLabels[nMassCuts] = {"NoCut", "50-150", "50-200", "50-250", "50-300", "50-350"};
const double massCutLow[nMassCuts] = {-1, 50, 50, 50, 50, 50};  // -1 means no lower cut
const double massCutHigh[nMassCuts] = {-1, 150, 200, 250, 300, 350}; // -1 means no upper cut*/
// Define global constants for mass bins
const int nMassCuts = 3;
const char* massCutLabels[nMassCuts] = {"50-80GeV", "80-110GeV", "#geq110 GeV"}; // For display
const char* fileSystemMassCutLabels[nMassCuts] = {"50-80GeV", "80-110GeV", "MassGE110GeV"}; // For file system
const double massCutLow[nMassCuts] = {50, 80, 110};  // Lower bound for each bin
const double massCutHigh[nMassCuts] = {80, 110, -1}; // Upper bound for each bin (-1 means no upper bound)


// Function to draw the ATLAS label, consistent with your provided macros
void DrawAtlasLabel() {
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(72);
    latex.SetTextColor(kBlack);

    // Draw "ATLAS" label
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.2, 0.88, "ATLAS");

    // Draw "Internal" label
    latex.SetTextFont(42);
    latex.DrawLatex(0.3, 0.88, "Internal");

    // Draw the energy and luminosity text
    latex.DrawLatex(0.2, 0.84, "#sqrt{s} = 13 TeV, 140 fb^{-1}");
}

// Create directories for each working point and mass cut
void Create2DPlotDirectories(const char* wp_label) {
    // Base directory for the specific working point
    TString baseDir = Form("../Plots/DB/2DplotsRdbpTbin/%s", wp_label);
    gSystem->Exec(Form("mkdir -p %s", baseDir.Data())); // Create the working point's main directory
    gSystem->Exec(Form("mkdir -p %s/Summary", baseDir.Data())); // Create a summary subdirectory
    
    // Individual mass cut subdirectories are not strictly needed for this script
    // but can be added if you plan to save per-mass-cut 2D plots.
}

/**
 * Helper function to dynamically set the palette range for a TH2F histogram.
 * The range is set symmetrically around 1.0, includes all data points with padding,
 * and ensures a minimum visual spread.
 * @param hist The TH2F histogram to adjust.
 */
void SetDynamicRatioPaletteRange(TH2F* hist) {
    double minVal = 1e18; // Initialize with a very large number
    double maxVal = -1e18; // Initialize with a very small number
    bool hasContent = false;

    // Find the actual min and max content of the histogram
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        for (int j = 1; j <= hist->GetNbinsY(); ++j) {
            double content = hist->GetBinContent(i, j);
            if (content != 0) { // Only consider non-zero content for range determination
                minVal = std::min(minVal, content);
                maxVal = std::max(maxVal, content);
                hasContent = true;
            }
        }
    }

    if (!hasContent) {
        // If no content, set a default sensible range around 1.0
        hist->SetMinimum(0.9);
        hist->SetMaximum(1.1);
        return;
    }

    // Calculate the raw range of the data
    double rawRange = maxVal - minVal;

    // Ensure a minimum visual spread, even if all values are very close
    const double minDisplayRange = 0.02; // e.g., 0.02 means 1% above and 1% below 1.0
    if (rawRange < minDisplayRange) {
        rawRange = minDisplayRange;
    }

    // Determine the symmetric extent needed around 1.0
    double center = 1.0;
    double distToMin = center - minVal;
    double distToMax = maxVal - center;
    double maxAbsDist = std::max(std::abs(distToMin), std::abs(distToMax));

    // Add 2% padding to the max absolute distance
    double paddedExtent = maxAbsDist * 1.02;

    // Set the new symmetric range
    double newMin = center - paddedExtent;
    double newMax = center + paddedExtent;

    // Apply the range to the histogram
    hist->SetMinimum(newMin);
    hist->SetMaximum(newMax);
}


// Main processing function for a single working point
void Process2DPlotsOnly(const char* wp_label) {
    std::cout << "\n--- Processing 2D Plots for Working Point: " << wp_label << " ---" << std::endl;

    // Create directories for the current working point
    Create2DPlotDirectories(wp_label);

    // Get current username for path construction
    std::string username = std::getenv("USER") ? std::getenv("USER") : "tamezza";
    //std::string basePath = "root://eosuser.cern.ch//eos/user/t/" + username + "/Documents/EasyBjets/ZbbJets/";
    // If running locally and files are in relative path, uncomment the line below and comment the above one
    std::string basePath = "../../";

    // Construct file paths for this working point
    std::string dataPath = basePath + "NtupleSlim/DB_" + wp_label + "wp_slim/data_" + wp_label + "wp_slim.root";
    std::string zbbPath = basePath + "NtupleSlim/DB_" + wp_label + "wp_slim/Zbby_" + wp_label + "wp_slim.root";
    std::string mjPath = basePath + "NtupleSlim/DB_" + wp_label + "wp_slim/gamma_jets_" + wp_label + "wp_slim.root";
    
    // Open ROOT files
    TFile* fdata = TFile::Open(dataPath.c_str(), "READ");
    TFile* fmc = TFile::Open(zbbPath.c_str(), "READ");
    TFile* fmultijets = TFile::Open(mjPath.c_str(), "READ");

    if (!fdata || !fmc || !fmultijets) {
        std::cerr << "Error: Could not open one or more ROOT files for WP " << wp_label << ". Skipping." << std::endl;
        // Clean up any opened files before returning
        if (fdata) fdata->Close();
        if (fmc) fmc->Close();
        if (fmultijets) fmultijets->Close();
        return;
    }

    TTree *tree_mc = (TTree*)fmc->Get("nominal");
    TTree *tree_data = (TTree*)fdata->Get("nominal");
    TTree *tree_mulijets = (TTree*)fmultijets->Get("nominal");

    if (!tree_mc || !tree_data || !tree_mulijets) {
        std::cerr << "Error: Could not retrieve 'nominal' tree from one or more ROOT files for WP " << wp_label << ". Skipping." << std::endl;
        fdata->Close(); fmc->Close(); fmultijets->Close();
        return;
    }
    
    // Create histograms for each mass cut and pT bin
    // Add wp_label to histogram names to ensure uniqueness across different working points
    TH1F* h_mc_R_DB_ptbin[nMassCuts][nPtBins];
    TH1F* h_data_R_DB_ptbin[nMassCuts][nPtBins];
    TH1F* h_multijets_R_DB_ptbin[nMassCuts][nPtBins];
    TH1F* h_combined_MC_R_DB_ptbin[nMassCuts][nPtBins];
    
    // Initialize histograms for each mass cut and pT bin
    for (int m = 0; m < nMassCuts; m++) {
        for (int i = 0; i < nPtBins; i++) {
            TString histName = Form("h_mc_R_DB_masscut%s_ptbin%d_%s", massCutLabels[m], i, wp_label);
            TString histTitle = Form("R_DB for %s (Mass Cut: %s, WP: %s)", ptBinLabels[i], massCutLabels[m], wp_label);
            h_mc_R_DB_ptbin[m][i] = new TH1F(histName, histTitle, 50, 0, 2);
            // h_mc_R_DB_ptbin[m][i]->SetFillColor(kRed); // Not needed for mean extraction
            // h_mc_R_DB_ptbin[m][i]->SetLineColor(kRed);
            
            histName = Form("h_data_R_DB_masscut%s_ptbin%d_%s", massCutLabels[m], i, wp_label);
            h_data_R_DB_ptbin[m][i] = new TH1F(histName, histTitle, 50, 0, 2);
            // h_data_R_DB_ptbin[m][i]->SetFillColor(kBlack);
            // h_data_R_DB_ptbin[m][i]->SetLineColor(kBlack);
            
            histName = Form("h_multijets_R_DB_masscut%s_ptbin%d_%s", massCutLabels[m], i, wp_label);
            h_multijets_R_DB_ptbin[m][i] = new TH1F(histName, histTitle, 50, 0, 2);
            // h_multijets_R_DB_ptbin[m][i]->SetFillColor(kBlue);
            // h_multijets_R_DB_ptbin[m][i]->SetLineColor(kBlue);
            
            histName = Form("h_combined_MC_R_DB_masscut%s_ptbin%d_%s", massCutLabels[m], i, wp_label);
            h_combined_MC_R_DB_ptbin[m][i] = new TH1F(histName, histTitle, 50, 0, 2);
            // h_combined_MC_R_DB_ptbin[m][i]->SetFillColor(kGreen+2);
            // h_combined_MC_R_DB_ptbin[m][i]->SetLineColor(kGreen+2);
        }
    }

    // Declare variables for tree branches
    std::vector<float>* ljet_pt = nullptr;
    std::vector<float>* ljet_m = nullptr;
    double weight, weight_mulijets;
    double R_DB;

    // Set branch addresses
    tree_mc->SetBranchAddress("ljet_pt", &ljet_pt);
    tree_data->SetBranchAddress("ljet_pt", &ljet_pt);
    tree_mulijets->SetBranchAddress("ljet_pt", &ljet_pt);

    tree_mc->SetBranchAddress("ljet_m", &ljet_m);
    tree_data->SetBranchAddress("ljet_m", &ljet_m);
    tree_mulijets->SetBranchAddress("ljet_m", &ljet_m);

    tree_mc->SetBranchAddress("R_DB", &R_DB);
    tree_data->SetBranchAddress("R_DB", &R_DB);
    tree_mulijets->SetBranchAddress("R_DB", &R_DB);

    tree_mc->SetBranchAddress("total_weight", &weight);
    tree_mulijets->SetBranchAddress("total_weight", &weight_mulijets);

    // Fill MC signal histograms with different mass cuts
    Long64_t nentries = tree_mc->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree_mc->GetEntry(i);
        if (ljet_pt->empty() || ljet_m->empty()) continue; // Skip if no jets
        double leadingJetPt = (*ljet_pt)[0]/1000.0; // Convert MeV to GeV
        double jetMass = (*ljet_m)[0]/1000.0;      // Convert MeV to GeV
        
        // Find the appropriate pT bin
        for (int bin = 0; bin < nPtBins; bin++) {
            if (leadingJetPt >= ptBins[bin] && leadingJetPt < ptBins[bin+1]) {
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
                    
                    // Fill histogram if passes cut
                    if (passesCut) {
                        h_mc_R_DB_ptbin[m][bin]->Fill(R_DB, weight);
                    }
                }
                break;
            }
        }
    }

    // Fill data histograms with different mass cuts
    nentries = tree_data->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree_data->GetEntry(i);
        if (ljet_pt->empty() || ljet_m->empty()) continue; // Skip if no jets
        double leadingJetPt = (*ljet_pt)[0]/1000.0; // Convert MeV to GeV
        double jetMass = (*ljet_m)[0]/1000.0;      // Convert MeV to GeV
        
        // Find the appropriate pT bin
        for (int bin = 0; bin < nPtBins; bin++) {
            if (leadingJetPt >= ptBins[bin] && leadingJetPt < ptBins[bin+1]) {
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
                    
                    // Fill histogram if passes cut
                    if (passesCut) {
                        h_data_R_DB_ptbin[m][bin]->Fill(R_DB);
                    }
                }
                break;
            }
        }
    }

    // Fill background histograms with different mass cuts
    nentries = tree_mulijets->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree_mulijets->GetEntry(i);
        if (ljet_pt->empty() || ljet_m->empty()) continue; // Skip if no jets
        double leadingJetPt = (*ljet_pt)[0]/1000.0; // Convert MeV to GeV
        double jetMass = (*ljet_m)[0]/1000.0;      // Convert MeV to GeV
        
        // Find the appropriate pT bin
        for (int bin = 0; bin < nPtBins; bin++) {
            if (leadingJetPt >= ptBins[bin] && leadingJetPt < ptBins[bin+1]) {
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
                    
                    // Fill histogram if passes cut
                    if (passesCut) {
                        h_multijets_R_DB_ptbin[m][bin]->Fill(R_DB, weight_mulijets);
                    }
                }
                break;
            }
        }
    }

    // Calculate scale factors for each mass cut separately
    for (int m = 0; m < nMassCuts; m++) {
        // Create histograms for scale factor calculation from jet mass
        TH1F *h_mc_ljet_m1 = new TH1F(Form("h_mc_ljet_m1_masscut%s_%s", massCutLabels[m], wp_label), "ljet_m", 20, 50, 150);
        TH1F *h_data_ljet_m1 = new TH1F(Form("h_data_m1_masscut%s_%s", massCutLabels[m], wp_label), "ljet_m", 20, 50, 150);
        TH1F *h_multijets_ljet_m1 = new TH1F(Form("h_multijets_ljet_m1_masscut%s_%s", massCutLabels[m], wp_label), "ljet_m", 20, 50, 150);

        // Fill MC signal jet mass histogram with mass cut
        nentries = tree_mc->GetEntries();
        for (Long64_t i = 0; i < nentries; i++) {
            tree_mc->GetEntry(i);
            if (ljet_m->empty()) continue;
            double jetMass = (*ljet_m)[0]/1000.0;
            
            bool passesCut = true;
            if (massCutLow[m] > 0 && jetMass < massCutLow[m]) {
                passesCut = false;
            }
            if (massCutHigh[m] > 0 && jetMass > massCutHigh[m]) {
                passesCut = false;
            }
            
            if (passesCut) {
                h_mc_ljet_m1->Fill(jetMass, weight);
            }
        }

        // Fill data jet mass histogram with mass cut
        nentries = tree_data->GetEntries();
        for (Long64_t i = 0; i < nentries; i++) {
            tree_data->GetEntry(i);
            if (ljet_m->empty()) continue;
            double jetMass = (*ljet_m)[0]/1000.0;
            
            bool passesCut = true;
            if (massCutLow[m] > 0 && jetMass < massCutLow[m]) {
                passesCut = false;
            }
            if (massCutHigh[m] > 0 && jetMass > massCutHigh[m]) {
                passesCut = false;
            }
            
            if (passesCut) {
                h_data_ljet_m1->Fill(jetMass);
            }
        }

        // Fill multijet background jet mass histogram with mass cut
        nentries = tree_mulijets->GetEntries();
        for (Long64_t i = 0; i < nentries; i++) {
            tree_mulijets->GetEntry(i);
            if (ljet_m->empty()) continue;
            double jetMass = (*ljet_m)[0]/1000.0;
            
            bool passesCut = true;
            if (massCutLow[m] > 0 && jetMass < massCutLow[m]) {
                passesCut = false;
            }
            if (massCutHigh[m] > 0 && jetMass > massCutHigh[m]) {
                passesCut = false;
            }
            
            if (passesCut) {
                h_multijets_ljet_m1->Fill(jetMass, weight_mulijets);
            }
        }
        
        // Define sideband regions (same as original code)
        int binLow1 = h_multijets_ljet_m1->FindBin(50);
        int binHigh1 = h_multijets_ljet_m1->FindBin(65);
        int binLow2 = h_multijets_ljet_m1->FindBin(110);
        int binHigh2 = h_multijets_ljet_m1->FindBin(150);
        
        // Create vectors to store x and y values for the fit
        std::vector<double> x_values;
        std::vector<double> y_values;
        std::vector<double> y_errors;
        
        // Extract bin contents for first sideband region
        for (int bin = binLow1; bin <= binHigh1; ++bin) {
            double bin_center = h_multijets_ljet_m1->GetBinCenter(bin);
            double data_content = h_data_ljet_m1->GetBinContent(bin);
            double mc_content = h_multijets_ljet_m1->GetBinContent(bin);
            
            // Avoid division by zero
            if (mc_content > 0) {
                x_values.push_back(bin_center);
                y_values.push_back(data_content / mc_content);
                
                // Calculate error on the ratio
                double data_error = h_data_ljet_m1->GetBinError(bin);
                double mc_error = h_multijets_ljet_m1->GetBinError(bin);
                double ratio_error = std::sqrt(std::pow(data_error/mc_content, 2) + 
                                         std::pow(data_content*mc_error/(mc_content*mc_content), 2));
                y_errors.push_back(ratio_error);
            }
        }
        
        // Extract bin contents for second sideband region
        for (int bin = binLow2; bin <= binHigh2; ++bin) {
            double bin_center = h_multijets_ljet_m1->GetBinCenter(bin);
            double data_content = h_data_ljet_m1->GetBinContent(bin);
            double mc_content = h_multijets_ljet_m1->GetBinContent(bin);
            
            // Avoid division by zero
            if (mc_content > 0) {
                x_values.push_back(bin_center);
                y_values.push_back(data_content / mc_content);
                
                // Calculate error on the ratio
                double data_error = h_data_ljet_m1->GetBinError(bin);
                double mc_error = h_multijets_ljet_m1->GetBinError(bin);
                double ratio_error = std::sqrt(std::pow(data_error/mc_content, 2) + 
                                         std::pow(data_content*mc_error/(mc_content*mc_content), 2));
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
        TF1* fit_func = new TF1(Form("fit_func_masscut%s_%s", massCutLabels[m], wp_label), "[0] + [1] * x", 50, 150);
        fit_func->SetParameters(1.0, 0.0);  // Initial guess
        
        // Perform the fit
        ratio_graph->Fit(Form("fit_func_masscut%s_%s", massCutLabels[m], wp_label), "RQ");  // "R" restricts fit to the specified range, "Q" for quiet
        
        // Get fit parameters
        double p0 = fit_func->GetParameter(0);  // Intercept
        double p1 = fit_func->GetParameter(1);  // Slope
        
        // Calculate average scale factor
        double avgScaleFactor = 0.0;
        int nBinsAvg = 0;
        
        // Calculate average SF in the full mass range
        for (int bin = 1; bin <= h_multijets_ljet_m1->GetNbinsX(); ++bin) {
            double bin_center = h_multijets_ljet_m1->GetBinCenter(bin);
            double scale_factor = p0 + p1 * bin_center;
            avgScaleFactor += scale_factor;
            nBinsAvg++;
        }
        
        if (nBinsAvg > 0) {
            avgScaleFactor /= nBinsAvg;
        } else {
            avgScaleFactor = 1.0; // Default to 1 if no bins to average
        }
        
        // Apply scale factor to background histograms
        for (int bin = 0; bin < nPtBins; bin++) {
            h_multijets_R_DB_ptbin[m][bin]->Scale(avgScaleFactor);
            
            // Create the combined MC histogram by adding signal and scaled background
            h_combined_MC_R_DB_ptbin[m][bin]->Add(h_mc_R_DB_ptbin[m][bin]);
            h_combined_MC_R_DB_ptbin[m][bin]->Add(h_multijets_R_DB_ptbin[m][bin]);
        }
        
        // Clean up
        delete h_mc_ljet_m1;
        delete h_data_ljet_m1;
        delete h_multijets_ljet_m1;
        delete ratio_graph;
        delete fit_func;
    }

    // Create arrays to store fit results for use in the 2D plots
    double mc_means[nMassCuts][nPtBins], mc_mean_errs[nMassCuts][nPtBins];
    double data_means[nMassCuts][nPtBins], data_mean_errs[nMassCuts][nPtBins];
    double bkg_means[nMassCuts][nPtBins], bkg_mean_errs[nMassCuts][nPtBins];
    double combined_means[nMassCuts][nPtBins], combined_mean_errs[nMassCuts][nPtBins];
    
    // For each mass cut, perform fits
    for (int m = 0; m < nMassCuts; m++) {
        // Perform fits for each pT bin with this mass cut
        for (int bin = 0; bin < nPtBins; bin++) {
            // MC signal fit
            TF1 *gausFit_mc = new TF1(Form("gausFit_mc_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label), "gaus", 0.9, 1.1);
            h_mc_R_DB_ptbin[m][bin]->Fit(gausFit_mc, "RQ");  // 'Q' for quiet mode
            
            // Get fit parameters
            mc_means[m][bin] = gausFit_mc->GetParameter(1);
            mc_mean_errs[m][bin] = gausFit_mc->GetParError(1);
            
            // Data fit
            TF1 *gausFit_data = new TF1(Form("gausFit_data_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label), "gaus", 0.9, 1.1);
            h_data_R_DB_ptbin[m][bin]->Fit(gausFit_data, "RQ");  // 'Q' for quiet mode
            
            // Get fit parameters
            data_means[m][bin] = gausFit_data->GetParameter(1);
            data_mean_errs[m][bin] = gausFit_data->GetParError(1);
            
            // Background fit
            TF1 *gausFit_bkg = new TF1(Form("gausFit_bkg_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label), "gaus", 0.9, 1.1);
            h_multijets_R_DB_ptbin[m][bin]->Fit(gausFit_bkg, "RQ");  // 'Q' for quiet mode
            
            // Get fit parameters
            bkg_means[m][bin] = gausFit_bkg->GetParameter(1);
            bkg_mean_errs[m][bin] = gausFit_bkg->GetParError(1);
            
            // Combined MC fit
            TF1 *gausFit_combined = new TF1(Form("gausFit_combined_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label), "gaus", 0.9, 1.1);
            h_combined_MC_R_DB_ptbin[m][bin]->Fit(gausFit_combined, "RQ");  // 'Q' for quiet mode
            
            // Get fit parameters
            combined_means[m][bin] = gausFit_combined->GetParameter(1);
            combined_mean_errs[m][bin] = gausFit_combined->GetParError(1);
            
            // Clean up fit functions
            delete gausFit_mc;
            delete gausFit_data;
            delete gausFit_bkg;
            delete gausFit_combined;
        }
    }
    
    // Create a summary plot comparing all mass cuts in one figure (for data/combined MC ratio)
    TCanvas* all_ratio_canvas = new TCanvas(Form("all_ratio_canvas_%s", wp_label), Form("Data/Combined MC Ratio for All pT and Mass Cuts (WP: %s)", wp_label), 800, 600);

    // Set canvas margins to make room for the ATLAS label
    all_ratio_canvas->SetTopMargin(0.12);  // Increased top margin for ATLAS label
    all_ratio_canvas->SetRightMargin(0.15); // Increased right margin for palette
    all_ratio_canvas->SetLeftMargin(0.12);  // Left margin
    all_ratio_canvas->SetBottomMargin(0.12); // Bottom margin
    all_ratio_canvas->cd();

    // Create a 2D plot to visualize all the ratios
    TH2F* ratio_map = new TH2F(Form("ratio_map_%s", wp_label), Form("Data/Combined MC Ratio (WP: %s);Mass Cut;p_{T} Bin", wp_label), nMassCuts, -0.5, nMassCuts-0.5, nPtBins, -0.5, nPtBins-0.5);

    // Fill the 2D plot
    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < nPtBins; bin++) {
            double ratio = (combined_means[m][bin] != 0) ? data_means[m][bin] / combined_means[m][bin] : 0.0;
            ratio_map->SetBinContent(m+1, bin+1, ratio);
            
            // Calculate error
            double ratio_err = (data_means[m][bin] != 0 && combined_means[m][bin] != 0) ? ratio * std::sqrt(
                std::pow(data_mean_errs[m][bin]/data_means[m][bin], 2) + 
                std::pow(combined_mean_errs[m][bin]/combined_means[m][bin], 2)) : 0.0;
            ratio_map->SetBinError(m+1, bin+1, ratio_err);
        }
    }

    // Style the plot
    ratio_map->SetStats(0);
    ratio_map->SetMarkerSize(1.2); // Make the text larger for readability

    // Set custom bin labels
    for (int m = 0; m < nMassCuts; m++) {
        ratio_map->GetXaxis()->SetBinLabel(m+1, massCutLabels[m]);
    }

    for (int bin = 0; bin < nPtBins; bin++) {
        // Simplify pT bin labels for readability
        TString simplifiedLabel = Form("%d-%d GeV", (int)ptBins[bin], (int)ptBins[bin+1]);
        ratio_map->GetYaxis()->SetBinLabel(bin+1, simplifiedLabel.Data());
    }

    ratio_map->GetXaxis()->SetLabelSize(0.025);
    ratio_map->GetYaxis()->SetLabelSize(0.025);
    ratio_map->GetXaxis()->SetTitleSize(0.035);
    ratio_map->GetYaxis()->SetTitleSize(0.035);
    ratio_map->GetZaxis()->SetTitleSize(0.03);

    // Remove original title (we'll add this with ATLAS label)
    ratio_map->SetTitle("");

    gStyle->SetNumberContours(100);

    // Set a nice range centered at 1
    //ratio_map->SetMinimum(0.97);
    //ratio_map->SetMaximum(1.01);

    SetDynamicRatioPaletteRange(ratio_map); // Apply dynamic range

    // For the text on the plot, make it more readable
    gStyle->SetPaintTextFormat("4.3f"); // Show 3 decimal places

    // Draw the plot
    ratio_map->Draw("COLZ TEXT");

    // Update the canvas to make sure the palette is created
    gPad->Update();

    // Get a pointer to the palette
    TPaletteAxis* palette = (TPaletteAxis*)ratio_map->GetListOfFunctions()->FindObject("palette");

    // Only proceed if we found the palette
    if (palette) {
        // Adjust palette position
        palette->SetX1NDC(0.86);
        palette->SetX2NDC(0.89);
        palette->SetY1NDC(0.12); // Match the bottom margin
        palette->SetY2NDC(0.88); // Match the top margin (adjusted for ATLAS label)
        palette->SetLabelSize(0.03);
        palette->SetTitleOffset(1.3);
        
        // Create vertical text for the palette
        TLatex palette_text;
        palette_text.SetNDC();
        palette_text.SetTextSize(0.035);
        palette_text.SetTextAlign(22); // Center-center aligned
        palette_text.SetTextAngle(90); // Rotate text 90 degrees
        palette_text.DrawLatex(0.97, 0.5, "#LT R_{DB}^{Data} #GT / #LT R_{DB}^{MC} #GT"); // Position to the right of palette
    }

    // Draw ATLAS label outside the plot (in the top margin)
    DrawAtlasLabel();

    all_ratio_canvas->SaveAs(Form("../Plots/DB/2DplotsRdbpTbin/%s/Summary/R_DB_all_ratio_summary.pdf", wp_label));
    delete all_ratio_canvas;
    delete ratio_map;

    //------------------------------------------------------------
    // Create 2D plot for data means
    //------------------------------------------------------------
    TCanvas* data_means_canvas = new TCanvas(Form("data_means_canvas_%s", wp_label), Form("Data R_{DB} Mean Values (WP: %s)", wp_label), 800, 600);

    // Set canvas margins
    data_means_canvas->SetTopMargin(0.12);
    data_means_canvas->SetRightMargin(0.15);
    data_means_canvas->SetLeftMargin(0.12);
    data_means_canvas->SetBottomMargin(0.12);
    data_means_canvas->cd();

    // Create a 2D plot
    TH2F* data_means_map = new TH2F(Form("data_means_map_%s", wp_label), Form("Data R_{DB} Mean Values (WP: %s);Mass Cut;p_{T} Bin", wp_label), nMassCuts, -0.5, nMassCuts-0.5, nPtBins, -0.5, nPtBins-0.5);

    // Fill the 2D plot
    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < nPtBins; bin++) {
            data_means_map->SetBinContent(m+1, bin+1, data_means[m][bin]);
            data_means_map->SetBinError(m+1, bin+1, data_mean_errs[m][bin]);
        }
    }

    // Style the plot
    data_means_map->SetStats(0);
    data_means_map->SetMarkerSize(1.2);

    // Set custom bin labels
    for (int m = 0; m < nMassCuts; m++) {
        data_means_map->GetXaxis()->SetBinLabel(m+1, massCutLabels[m]);
    }

    for (int bin = 0; bin < nPtBins; bin++) {
        TString simplifiedLabel = Form("%d-%d GeV", (int)ptBins[bin], (int)ptBins[bin+1]);
        data_means_map->GetYaxis()->SetBinLabel(bin+1, simplifiedLabel.Data());
    }

    data_means_map->GetXaxis()->SetLabelSize(0.025);
    data_means_map->GetYaxis()->SetLabelSize(0.025);
    data_means_map->GetXaxis()->SetTitleSize(0.035);
    data_means_map->GetYaxis()->SetTitleSize(0.035);
    data_means_map->GetZaxis()->SetTitleSize(0.03);

    // Remove original title
    data_means_map->SetTitle("");

    // Set the color palette - use a different color scheme for direct values
    gStyle->SetPalette(kViridis);

    // Set the range based on actual data values
    /*
    double min_val = 0.95;
    double max_val = 1.05;
    data_means_map->SetMinimum(min_val);
    data_means_map->SetMaximum(max_val);
    */

    SetDynamicRatioPaletteRange(data_means_map); // Apply dynamic range

    // For the text on the plot
    gStyle->SetPaintTextFormat("4.3f");

    // Draw the plot
    data_means_map->Draw("COLZ TEXT");

    // Update the canvas
    gPad->Update();

    // Get a pointer to the palette
    TPaletteAxis* data_palette = (TPaletteAxis*)data_means_map->GetListOfFunctions()->FindObject("palette");

    // Only proceed if we found the palette
    if (data_palette) {
        // Adjust palette position
        data_palette->SetX1NDC(0.86);
        data_palette->SetX2NDC(0.89);
        data_palette->SetY1NDC(0.12);
        data_palette->SetY2NDC(0.88);
        data_palette->SetLabelSize(0.03);
        
        // Create vertical text for the palette
        TLatex data_palette_text;
        data_palette_text.SetNDC();
        data_palette_text.SetTextSize(0.035);
        data_palette_text.SetTextAlign(22);
        data_palette_text.SetTextAngle(90);
        data_palette_text.DrawLatex(0.97, 0.5, "#LT R_{DB}^{Data} #GT");
    }

    // Draw ATLAS label outside the plot
    DrawAtlasLabel();

    data_means_canvas->SaveAs(Form("../Plots/DB/2DplotsRdbpTbin/%s/Summary/Data_R_DB_mean_values.pdf", wp_label));
    delete data_means_canvas;
    delete data_means_map;


    //------------------------------------------------------------
    // Create 2D plot for MC means
    //------------------------------------------------------------
    TCanvas* mc_means_canvas = new TCanvas(Form("mc_means_canvas_%s", wp_label), Form("MC R_{DB} Mean Values (WP: %s)", wp_label), 800, 600);

    // Set canvas margins
    mc_means_canvas->SetTopMargin(0.12);
    mc_means_canvas->SetRightMargin(0.15);
    mc_means_canvas->SetLeftMargin(0.12);
    mc_means_canvas->SetBottomMargin(0.12);
    mc_means_canvas->cd();

    // Create a 2D plot
    TH2F* mc_means_map = new TH2F(Form("mc_means_map_%s", wp_label), Form("MC R_{DB} Mean Values (WP: %s);Mass Cut;p_{T} Bin", wp_label), nMassCuts, -0.5, nMassCuts-0.5, nPtBins, -0.5, nPtBins-0.5);

    // Fill the 2D plot
    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < nPtBins; bin++) {
            mc_means_map->SetBinContent(m+1, bin+1, mc_means[m][bin]);
            mc_means_map->SetBinError(m+1, bin+1, mc_mean_errs[m][bin]);
        }
    }

    // Style the plot
    mc_means_map->SetStats(0);
    mc_means_map->SetMarkerSize(1.2);

    // Set custom bin labels
    for (int m = 0; m < nMassCuts; m++) {
        mc_means_map->GetXaxis()->SetBinLabel(m+1, massCutLabels[m]);
    }

    for (int bin = 0; bin < nPtBins; bin++) {
        TString simplifiedLabel = Form("%d-%d GeV", (int)ptBins[bin], (int)ptBins[bin+1]);
        mc_means_map->GetYaxis()->SetBinLabel(bin+1, simplifiedLabel.Data());
    }

    mc_means_map->GetXaxis()->SetLabelSize(0.025);
    mc_means_map->GetYaxis()->SetLabelSize(0.025);
    mc_means_map->GetXaxis()->SetTitleSize(0.035);
    mc_means_map->GetYaxis()->SetTitleSize(0.035);
    mc_means_map->GetZaxis()->SetTitleSize(0.03);

    // Remove original title
    mc_means_map->SetTitle("");

    // Set the color palette
    gStyle->SetPalette(kViridis);

    // Set the range based on actual data values
    /*mc_means_map->SetMinimum(min_val);
    mc_means_map->SetMaximum(max_val);*/

    SetDynamicRatioPaletteRange(mc_means_map); // Apply dynamic range

    // For the text on the plot
    gStyle->SetPaintTextFormat("4.3f");

    // Draw the plot
    mc_means_map->Draw("COLZ TEXT");

    // Update the canvas
    gPad->Update();

    // Get a pointer to the palette
    TPaletteAxis* mc_palette = (TPaletteAxis*)mc_means_map->GetListOfFunctions()->FindObject("palette");

    // Only proceed if we found the palette
    if (mc_palette) {
        // Adjust palette position
        mc_palette->SetX1NDC(0.86);
        mc_palette->SetX2NDC(0.89);
        mc_palette->SetY1NDC(0.12);
        mc_palette->SetY2NDC(0.88);
        mc_palette->SetLabelSize(0.03);
        
        // Create vertical text for the palette
        TLatex mc_palette_text;
        mc_palette_text.SetNDC();
        mc_palette_text.SetTextSize(0.035);
        mc_palette_text.SetTextAlign(22);
        mc_palette_text.SetTextAngle(90);
        mc_palette_text.DrawLatex(0.97, 0.5, "#LT R_{DB}^{Signal MC} #GT");
    }

    // Draw ATLAS label outside the plot
    DrawAtlasLabel();

    mc_means_canvas->SaveAs(Form("../Plots/DB/2DplotsRdbpTbin/%s/Summary/MC_R_DB_mean_values.pdf", wp_label));
    delete mc_means_canvas;
    delete mc_means_map;


    //------------------------------------------------------------
    // Create 2D plot for background means
    //------------------------------------------------------------
    TCanvas* bkg_means_canvas = new TCanvas(Form("bkg_means_canvas_%s", wp_label), Form("Background R_{DB} Mean Values (WP: %s)", wp_label), 800, 600);

    // Set canvas margins
    bkg_means_canvas->SetTopMargin(0.12);
    bkg_means_canvas->SetRightMargin(0.15);
    bkg_means_canvas->SetLeftMargin(0.12);
    bkg_means_canvas->SetBottomMargin(0.12);
    bkg_means_canvas->cd();

    // Create a 2D plot
    TH2F* bkg_means_map = new TH2F(Form("bkg_means_map_%s", wp_label), Form("Background R_{DB} Mean Values (WP: %s);Mass Cut;p_{T} Bin", wp_label), nMassCuts, -0.5, nMassCuts-0.5, nPtBins, -0.5, nPtBins-0.5);

    // Fill the 2D plot
    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < nPtBins; bin++) {
            bkg_means_map->SetBinContent(m+1, bin+1, bkg_means[m][bin]);
            bkg_means_map->SetBinError(m+1, bin+1, bkg_mean_errs[m][bin]);
        }
    }

    // Style the plot
    bkg_means_map->SetStats(0);
    bkg_means_map->SetMarkerSize(1.2);

    // Set custom bin labels
    for (int m = 0; m < nMassCuts; m++) {
        bkg_means_map->GetXaxis()->SetBinLabel(m+1, massCutLabels[m]);
    }

    for (int bin = 0; bin < nPtBins; bin++) {
        TString simplifiedLabel = Form("%d-%d GeV", (int)ptBins[bin], (int)ptBins[bin+1]);
        bkg_means_map->GetYaxis()->SetBinLabel(bin+1, simplifiedLabel.Data());
    }

    bkg_means_map->GetXaxis()->SetLabelSize(0.025);
    bkg_means_map->GetYaxis()->SetLabelSize(0.025);
    bkg_means_map->GetXaxis()->SetTitleSize(0.035);
    bkg_means_map->GetYaxis()->SetTitleSize(0.035);
    bkg_means_map->GetZaxis()->SetTitleSize(0.03);

    // Remove original title
    bkg_means_map->SetTitle("");

    // Set the color palette
    gStyle->SetPalette(kViridis);

    // Set the range based on actual data values
    //bkg_means_map->SetMinimum(min_val);
    //bkg_means_map->SetMaximum(max_val);

    SetDynamicRatioPaletteRange(bkg_means_map); // Apply dynamic range

    // For the text on the plot
    gStyle->SetPaintTextFormat("4.3f");

    // Draw the plot
    bkg_means_map->Draw("COLZ TEXT");

    // Update the canvas
    gPad->Update();

    // Get a pointer to the palette
    TPaletteAxis* bkg_palette = (TPaletteAxis*)bkg_means_map->GetListOfFunctions()->FindObject("palette");

    // Only proceed if we found the palette
    if (bkg_palette) {
        // Adjust palette position
        bkg_palette->SetX1NDC(0.86);
        bkg_palette->SetX2NDC(0.89);
        bkg_palette->SetY1NDC(0.12);
        bkg_palette->SetY2NDC(0.88);
        bkg_palette->SetLabelSize(0.03);
        
        // Create vertical text for the palette
        TLatex bkg_palette_text;
        bkg_palette_text.SetNDC();
        bkg_palette_text.SetTextSize(0.035);
        bkg_palette_text.SetTextAlign(22);
        bkg_palette_text.SetTextAngle(90);
        bkg_palette_text.DrawLatex(0.97, 0.5, "#LT R_{DB}^{Multijets} #GT");
    }

    // Draw ATLAS label outside the plot
    DrawAtlasLabel();

    bkg_means_canvas->SaveAs(Form("../Plots/DB/2DplotsRdbpTbin/%s/Summary/Bkg_R_DB_mean_values.pdf", wp_label));
    delete bkg_means_canvas;
    delete bkg_means_map;

    //------------------------------------------------------------
    // Create 2D plot for combined MC means
    //------------------------------------------------------------
    TCanvas* combined_means_canvas = new TCanvas(Form("combined_means_canvas_%s", wp_label), Form("Combined MC R_{DB} Mean Values (WP: %s)", wp_label), 800, 600);

    // Set canvas margins
    combined_means_canvas->SetTopMargin(0.12);
    combined_means_canvas->SetRightMargin(0.15);
    combined_means_canvas->SetLeftMargin(0.12);
    combined_means_canvas->SetBottomMargin(0.12);
    combined_means_canvas->cd();

    // Create a 2D plot
    TH2F* combined_means_map = new TH2F(Form("combined_means_map_%s", wp_label), Form("Combined MC R_{DB} Mean Values (WP: %s);Mass Cut;p_{T} Bin", wp_label), nMassCuts, -0.5, nMassCuts-0.5, nPtBins, -0.5, nPtBins-0.5);

    // Fill the 2D plot
    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < nPtBins; bin++) {
            combined_means_map->SetBinContent(m+1, bin+1, combined_means[m][bin]);
            combined_means_map->SetBinError(m+1, bin+1, combined_mean_errs[m][bin]);
        }
    }

    // Style the plot
    combined_means_map->SetStats(0);
    combined_means_map->SetMarkerSize(1.2);

    // Set custom bin labels
    for (int m = 0; m < nMassCuts; m++) {
        combined_means_map->GetXaxis()->SetBinLabel(m+1, massCutLabels[m]);
    }

    for (int bin = 0; bin < nPtBins; bin++) {
        TString simplifiedLabel = Form("%d-%d GeV", (int)ptBins[bin], (int)ptBins[bin+1]);
        combined_means_map->GetYaxis()->SetBinLabel(bin+1, simplifiedLabel.Data());
    }

    combined_means_map->GetXaxis()->SetLabelSize(0.025);
    combined_means_map->GetYaxis()->SetLabelSize(0.025);
    combined_means_map->GetXaxis()->SetTitleSize(0.035);
    combined_means_map->GetYaxis()->SetTitleSize(0.035);
    combined_means_map->GetZaxis()->SetTitleSize(0.03);

    // Remove original title
    combined_means_map->SetTitle("");

    // Set the color palette       
    gStyle->SetPalette(kViridis);

    // Set the range based on actual data values
    //combined_means_map->SetMinimum(min_val);
    //combined_means_map->SetMaximum(max_val);

    SetDynamicRatioPaletteRange(combined_means_map); // Apply dynamic range

    // For the text on the plot
    gStyle->SetPaintTextFormat("4.3f");

    // Draw the plot
    combined_means_map->Draw("COLZ TEXT");

    // Update the canvas
    gPad->Update();

    // Get a pointer to the palette
    TPaletteAxis* combined_palette = (TPaletteAxis*)combined_means_map->GetListOfFunctions()->FindObject("palette");

    // Only proceed if we found the palette
    if (combined_palette) {
        // Adjust palette position
        combined_palette->SetX1NDC(0.86);
        combined_palette->SetX2NDC(0.89);
        combined_palette->SetY1NDC(0.12);
        combined_palette->SetY2NDC(0.88);
        combined_palette->SetLabelSize(0.03);
        
        // Create vertical text for the palette
        TLatex combined_palette_text;
        combined_palette_text.SetNDC();
        combined_palette_text.SetTextSize(0.035);
        combined_palette_text.SetTextAlign(22);
        combined_palette_text.SetTextAngle(90);
        combined_palette_text.DrawLatex(0.97, 0.5, "#LT R_{DB}^{Combined MC} #GT");
    }

    // Draw ATLAS label outside the plot
    DrawAtlasLabel();

    combined_means_canvas->SaveAs(Form("../Plots/DB/2DplotsRdbpTbin/%s/Summary/Combined_MC_R_DB_mean_values.pdf", wp_label));
    delete combined_means_canvas;
    delete combined_means_map;


    // Create a summary table file
    ofstream summaryTable(Form("../Plots/DB/2DplotsRdbpTbin/%s/Summary/R_DB_summary_table.txt", wp_label));
    if (summaryTable.is_open()) {
        summaryTable << "R_DB SUMMARY TABLE FOR ALL MASS CUTS AND PT BINS (WP: " << wp_label << ")" << std::endl;
        summaryTable << "================================================" << std::endl << std::endl;
        
        for (int m = 0; m < nMassCuts; m++) {
            if (m == 0) {
                summaryTable << "MASS CUT: No Cut" << std::endl;
            } else {
                summaryTable << "MASS CUT: " << massCutLow[m] << " < m < " << massCutHigh[m] << " GeV" << std::endl;
            }
            summaryTable << "------------------------------------------------" << std::endl;
            
            // Table header
            summaryTable << std::left << std::setw(15) << "pT Bin" 
                         << std::setw(25) << "Data Mean ± Error" 
                         << std::setw(25) << "MC Mean ± Error"
                         << std::setw(25) << "#gamma + jets Mean ± Error"
                         << std::setw(25) << "Combined Mean ± Error"
                         << std::setw(15) << "Data/MC Ratio"
                         << std::setw(15) << "Data/Comb Ratio"
                         << std::setw(15) << "#gamma + jets/MC Ratio" << std::endl;
            summaryTable << "------------------------------------------------" << std::endl;
            
            for (int bin = 0; bin < nPtBins; bin++) {
                double data_mc_ratio = (mc_means[m][bin] != 0) ? data_means[m][bin] / mc_means[m][bin] : 0.0;
                double data_comb_ratio = (combined_means[m][bin] != 0) ? data_means[m][bin] / combined_means[m][bin] : 0.0;
                double bkg_mc_ratio = (mc_means[m][bin] != 0) ? bkg_means[m][bin] / mc_means[m][bin] : 0.0;
                
                summaryTable << std::left << std::setw(15) << ptBinLabels[bin]
                             << std::setw(25) << Form("%.4f ± %.4f", data_means[m][bin], data_mean_errs[m][bin])
                             << std::setw(25) << Form("%.4f ± %.4f", mc_means[m][bin], mc_mean_errs[m][bin])
                             << std::setw(25) << Form("%.4f ± %.4f", bkg_means[m][bin], bkg_mean_errs[m][bin])
                             << std::setw(25) << Form("%.4f ± %.4f", combined_means[m][bin], combined_mean_errs[m][bin])
                             << std::setw(15) << Form("%.4f", data_mc_ratio)
                             << std::setw(15) << Form("%.4f", data_comb_ratio)
                             << std::setw(15) << Form("%.4f", bkg_mc_ratio) << std::endl;
            }
            
            summaryTable << std::endl << std::endl;
        }
        
        summaryTable.close();
        std::cout << "Summary table saved to ../Plots/DB/2DplotsRdbpTbin/" << wp_label << "/Summary/R_DB_summary_table.txt" << std::endl;
    } else {
        std::cerr << "Error: Unable to open summary table file for writing for WP " << wp_label << "." << std::endl;
    }
    
    // Clean up histograms
    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < nPtBins; bin++) {
            delete h_mc_R_DB_ptbin[m][bin];
            delete h_data_R_DB_ptbin[m][bin];
            delete h_multijets_R_DB_ptbin[m][bin];
            delete h_combined_MC_R_DB_ptbin[m][bin];
        }
    }
    
    // Close input files
    fdata->Close();
    fmc->Close();
    fmultijets->Close();
    
    std::cout << "Finished processing 2D plots for Working Point: " << wp_label << std::endl;
}

// Main function to run analysis for all specified working points
void Rdb_2D_PlotsOnly() {
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    gROOT->LoadMacro("AtlasLabels.C");

    const char* workingPoints[] = {"0p25", "0p46", "0p74"};
    for (const char* wp : workingPoints) {
        Process2DPlotsOnly(wp);
    }
    std::cout << "\n--- All working points for 2D plots processed ---" << std::endl;
}