#include "../../Util/AtlasStyle.C"
#include <iostream>
#include <string>
#include <vector>
#include <iomanip> // For std::setw
#include <fstream> // For file output

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
class THStack;
class TPad;
class TAxis;
class TGraphAsymmErrors;


// Global constants for pT and mass bins (must be defined globally or passed)
const int nPtBins = 3;
const double ptBins[nPtBins + 1] = {450, 600, 900, 1200}; // 450-600, 600-900, 900-1200
const char* ptBinLabels[nPtBins] = {"450 < p_{T} < 600 GeV", "600 < p_{T} < 900 GeV", "900 < p_{T} < 1200 GeV"};

/*
const int nMassCuts = 6;
const char* massCutLabels[nMassCuts] = {"NoCut", "50-150", "50-200", "50-250", "50-300", "50-350"};
const double massCutLow[nMassCuts] = {-1, 50, 50, 50, 50, 50};  // -1 means no lower cut
const double massCutHigh[nMassCuts] = {-1, 150, 200, 250, 300, 350}; // -1 means no upper cut
*/

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

// Function to save fit parameters to a file
void SaveFitParameters(const char* filename, const char* fitType, 
                       double mean, double meanErr, 
                       double sigma, double sigmaErr,
                       double amp, double ampErr,
                       double chi2, double ndf) {

    ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    outFile << "Fit Type: " << fitType << std::endl;
    outFile << "Mean: " << mean << " ± " << meanErr << std::endl;
    outFile << "Sigma: " << sigma << " ± " << sigmaErr << std::endl;
    outFile << "Amplitude: " << amp << " ± " << ampErr << std::endl;
    outFile << "Chi2/NDF: " << chi2/ndf << " (" << chi2 << "/" << ndf << ")" << std::endl;

    outFile.close();
    std::cout << "Fit parameters saved to " << filename << std::endl;
}

// Function to create directories for each working point and mass cut
void CreateDirectories(const char* wp_label) {
    // Base directory for the specific working point
    TString baseDir = Form("../Plots/DB/PlotsRdbpTbin/%s", wp_label);
    gSystem->Exec(Form("mkdir -p %s/Summary", baseDir.Data()));
    
    // Create a directory for each mass cut within the working point's folder
    for (int i = 0; i < nMassCuts; i++) {
        TString dirCmd = Form("mkdir -p %s/%s", baseDir.Data(), fileSystemMassCutLabels[i]);
        gSystem->Exec(dirCmd.Data());
    }
}

// Main processing function for a single working point
void ProcessWorkingPoint(const char* wp_label) {
    std::cout << "\n--- Processing Working Point: " << wp_label << " ---" << std::endl;

    // Create directories for the current working point
    CreateDirectories(wp_label);

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
    TH1F* h_mc_R_DB_ptbin[nMassCuts][nPtBins];
    TH1F* h_data_R_DB_ptbin[nMassCuts][nPtBins];
    TH1F* h_multijets_R_DB_ptbin[nMassCuts][nPtBins];
    TH1F* h_combined_MC_R_DB_ptbin[nMassCuts][nPtBins];
    
    // Initialize histograms for each mass cut and pT bin
    for (int m = 0; m < nMassCuts; m++) {
        for (int i = 0; i < nPtBins; i++) {
            TString histName = Form("h_mc_R_DB_masscut%s_ptbin%d_%s", massCutLabels[m], i, wp_label);
            TString histTitle = Form("R_DB for %s (Mass Cut: %s)", ptBinLabels[i], massCutLabels[m]);
            h_mc_R_DB_ptbin[m][i] = new TH1F(histName, histTitle, 50, 0, 2);
            h_mc_R_DB_ptbin[m][i]->SetFillColor(kRed);
            h_mc_R_DB_ptbin[m][i]->SetLineColor(kRed);
            
            histName = Form("h_data_R_DB_masscut%s_ptbin%d_%s", massCutLabels[m], i, wp_label);
            h_data_R_DB_ptbin[m][i] = new TH1F(histName, histTitle, 50, 0, 2);
            h_data_R_DB_ptbin[m][i]->SetFillColor(kBlack);
            h_data_R_DB_ptbin[m][i]->SetLineColor(kBlack);
            
            histName = Form("h_multijets_R_DB_masscut%s_ptbin%d_%s", massCutLabels[m], i, wp_label);
            h_multijets_R_DB_ptbin[m][i] = new TH1F(histName, histTitle, 50, 0, 2);
            h_multijets_R_DB_ptbin[m][i]->SetFillColor(kBlue);
            h_multijets_R_DB_ptbin[m][i]->SetLineColor(kBlue);
            
            // Initialize the combined MC histogram
            histName = Form("h_combined_MC_R_DB_masscut%s_ptbin%d_%s", massCutLabels[m], i, wp_label);
            histTitle = Form("Combined MC R_DB for %s (Mass Cut: %s)", ptBinLabels[i], massCutLabels[m]);
            h_combined_MC_R_DB_ptbin[m][i] = new TH1F(histName, histTitle, 50, 0, 2);
            h_combined_MC_R_DB_ptbin[m][i]->SetFillColor(kGreen+2);
            h_combined_MC_R_DB_ptbin[m][i]->SetLineColor(kGreen+2);
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

    // Set branch addresses for all trees
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
        
        // Create arrays to store x and y values for the fit
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
        ratio_graph->Fit(Form("fit_func_masscut%s_%s", massCutLabels[m], wp_label), "R");  // "R" restricts fit to the specified range
        
        // Get fit parameters
        double p0 = fit_func->GetParameter(0);  // Intercept
        double p1 = fit_func->GetParameter(1);  // Slope
        // double p0_err = fit_func->GetParError(0);
        // double p1_err = fit_func->GetParError(1);
        
        // std::printf("Mass cut: %s, Linear fit results: SF(m) = %.4f (± %.4f) + %.8f (± %.8f) * m\n", 
        //        massCutLabels[m], p0, p0_err, p1, p1_err);
               
        // Calculate average scale factor for R_DB histogram
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
        std::printf("Mass cut: %s, Average Scale Factor: %.4f\n", massCutLabels[m], avgScaleFactor);
        
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
        delete ratio_graph; // Also delete the graph created in this scope
        delete fit_func;    // Also delete the fit function
    }

    // Create arrays to store fit results for each mass cut
    double mc_means[nMassCuts][nPtBins], mc_mean_errs[nMassCuts][nPtBins];
    double data_means[nMassCuts][nPtBins], data_mean_errs[nMassCuts][nPtBins];
    double bkg_means[nMassCuts][nPtBins], bkg_mean_errs[nMassCuts][nPtBins];
    double combined_means[nMassCuts][nPtBins], combined_mean_errs[nMassCuts][nPtBins];
    
    double mc_sigmas[nMassCuts][nPtBins], mc_sigma_errs[nMassCuts][nPtBins];
    double data_sigmas[nMassCuts][nPtBins], data_sigma_errs[nMassCuts][nPtBins];
    double bkg_sigmas[nMassCuts][nPtBins], bkg_sigma_errs[nMassCuts][nPtBins];
    double combined_sigmas[nMassCuts][nPtBins], combined_sigma_errs[nMassCuts][nPtBins];
    
    double ptBinCenters[nPtBins];
    for (int i = 0; i < nPtBins; i++) {
        ptBinCenters[i] = (ptBins[i] + ptBins[i+1]) / 2;
        // For the last bin, just use a reasonable value for plotting if it's the end of range
        if (i == nPtBins - 1) {
            ptBinCenters[i] = ptBins[i] + (ptBins[i+1] - ptBins[i]) / 2;
        }
    }

    // For each mass cut, perform fits and create plots
    for (int m = 0; m < nMassCuts; m++) {
        // Create a summary file for this mass cut
        TString summaryFileName = Form("../Plots/DB/PlotsRdbpTbin/%s/%s/R_DB_pT_binned_fit_summary.txt", wp_label, fileSystemMassCutLabels[m]);
        ofstream summaryFile(summaryFileName.Data());
        
        if (summaryFile.is_open()) {
            summaryFile << "Summary of R_DB Gaussian Fits in pT Bins (Mass Cut: " << massCutLabels[m] << ", WP: " << wp_label << ")" << std::endl;
            summaryFile << "----------------------------------------" << std::endl << std::endl;
        }
        
        // Perform fits for each pT bin with this mass cut
        for (int bin = 0; bin < nPtBins; bin++) {
            // MC signal fit
            TCanvas* mc_canvas = new TCanvas(Form("mc_canvas_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label), 
                                           Form("MC Signal R_DB for %s (Mass Cut: %s, WP: %s)", ptBinLabels[bin], massCutLabels[m], wp_label), 800, 600);
            mc_canvas->cd();
            
            TF1 *gausFit_mc = new TF1(Form("gausFit_mc_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label), "gaus", 0.9, 1.1);
            h_mc_R_DB_ptbin[m][bin]->Fit(gausFit_mc, "RQ"); // "Q" for quiet mode, "R" for range
            
            double chi_mc = gausFit_mc->GetChisquare();
            int ndf_mc = gausFit_mc->GetNDF();
            double chi2NDF_mc = (ndf_mc > 0) ? chi_mc / ndf_mc : 0.0;
            
            h_mc_R_DB_ptbin[m][bin]->Draw("HIST");
            gausFit_mc->SetLineColor(kViolet-4);
            gausFit_mc->Draw("same");
            
            h_mc_R_DB_ptbin[m][bin]->SetMaximum(h_mc_R_DB_ptbin[m][bin]->GetMaximum() * 1.25);
            h_mc_R_DB_ptbin[m][bin]->SetMinimum(0); // This line forces the y-axis to start at zero
            h_mc_R_DB_ptbin[m][bin]->GetXaxis()->SetTitle("R_{DB}");
            h_mc_R_DB_ptbin[m][bin]->GetYaxis()->SetTitle("Events");
            
            double amp_mc = gausFit_mc->GetParameter(0);
            double Mean_mc = gausFit_mc->GetParameter(1);
            double Sigma_mc = gausFit_mc->GetParameter(2);
            
            double ampErr_mc = gausFit_mc->GetParError(0);
            double meanErr_mc = gausFit_mc->GetParError(1);
            double sigmaErr_mc = gausFit_mc->GetParError(2);
            
            // Store results for comparison plots
            mc_means[m][bin] = Mean_mc;
            mc_mean_errs[m][bin] = meanErr_mc;
            mc_sigmas[m][bin] = Sigma_mc;
            mc_sigma_errs[m][bin] = sigmaErr_mc;
            
            TLegend *Legend1 = new TLegend(0.6, 0.75, 0.9, 0.9);
            Legend1->SetBorderSize(0);
            Legend1->SetTextSize(0.025);
            Legend1->AddEntry(h_mc_R_DB_ptbin[m][bin], "Z (#rightarrow bb) + #gamma ", "l");
            Legend1->AddEntry(gausFit_mc, "Gaussian Fit", "l");
            
            TLatex text;
            text.SetTextSize(0.03);
            text.DrawLatexNDC(0.65, 0.7, Form("#mu = %.3f #pm %.3f", Mean_mc, meanErr_mc));
            text.DrawLatexNDC(0.65, 0.65, Form("#sigma = %.3f #pm %.3f", Sigma_mc, sigmaErr_mc));
            text.DrawLatexNDC(0.65, 0.6, Form("#chi^{2}/NDF = %.3f", chi2NDF_mc));
            Legend1->Draw();
            
            DrawAtlasLabel();
            text.DrawLatexNDC(0.2, 0.8, ptBinLabels[bin]);
            
            // Add mass cut information
            if (m == 0) {
                text.DrawLatexNDC(0.2, 0.75, "No Mass Cut");
            } else {
                text.DrawLatexNDC(0.2, 0.75, Form("Mass Cut: %d < m < %d GeV", (int)massCutLow[m], (int)massCutHigh[m]));
            }
            
            mc_canvas->SaveAs(Form("../Plots/DB/PlotsRdbpTbin/%s/%s/Fit_R_DB_mc_ptbin%d.pdf", wp_label, fileSystemMassCutLabels[m], bin));
            
            // Save fit parameters
            SaveFitParameters(Form("../Plots/DB/PlotsRdbpTbin/%s/%s/mc_fit_parameters_ptbin%d.txt", wp_label, fileSystemMassCutLabels[m], bin),
                            Form("MC: Z (#rightarrow bb) + #gamma  (%s, Mass Cut: %s)", ptBinLabels[bin], massCutLabels[m]),
                            Mean_mc, meanErr_mc, Sigma_mc, sigmaErr_mc,
                            amp_mc, ampErr_mc, chi_mc, ndf_mc);
                            
            // Add to summary file
            if (summaryFile.is_open()) {
                summaryFile << "MC: Z (#rightarrow bb) + #gamma  (" << ptBinLabels[bin] << ")" << std::endl;
                summaryFile << "Mean: " << Mean_mc << " ± " << meanErr_mc << std::endl;
                summaryFile << "Sigma: " << Sigma_mc << " ± " << sigmaErr_mc << std::endl;
                summaryFile << "Chi2/NDF: " << chi2NDF_mc << " (" << chi_mc << "/" << ndf_mc << ")" << std::endl << std::endl;
            }
            
            delete mc_canvas; // Clean up canvas
            delete gausFit_mc; // Clean up fit function
            delete Legend1; // Clean up legend

            // Data fit
            TCanvas* data_canvas = new TCanvas(Form("data_canvas_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label), 
                                             Form("Data R_DB for %s (Mass Cut: %s, WP: %s)", ptBinLabels[bin], massCutLabels[m], wp_label), 800, 600);
            data_canvas->cd();
            
            TF1 *gausFit_data = new TF1(Form("gausFit_data_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label), "gaus", 0.9, 1.1);
            h_data_R_DB_ptbin[m][bin]->Fit(gausFit_data, "RQ");
            
            double chi_data = gausFit_data->GetChisquare();
            int ndf_data = gausFit_data->GetNDF();
            double chi2NDF_data = (ndf_data > 0) ? chi_data / ndf_data : 0.0;
            
            h_data_R_DB_ptbin[m][bin]->Draw("EP");
            gausFit_data->SetLineColor(kViolet-4);
            gausFit_data->Draw("same");
            
            h_data_R_DB_ptbin[m][bin]->SetMaximum(h_data_R_DB_ptbin[m][bin]->GetMaximum() * 1.25);
            h_data_R_DB_ptbin[m][bin]->GetXaxis()->SetTitle("R_{DB}");
            h_data_R_DB_ptbin[m][bin]->GetYaxis()->SetTitle("Events");
            
            double amp_data = gausFit_data->GetParameter(0);
            double Mean_data = gausFit_data->GetParameter(1);
            double Sigma_data = gausFit_data->GetParameter(2);
            
            double ampErr_data = gausFit_data->GetParError(0);
            double meanErr_data = gausFit_data->GetParError(1);
            double sigmaErr_data = gausFit_data->GetParError(2);
            
            // Store results for comparison plots
            data_means[m][bin] = Mean_data;
            data_mean_errs[m][bin] = meanErr_data;
            data_sigmas[m][bin] = Sigma_data;
            data_sigma_errs[m][bin] = sigmaErr_data;
            
            TLegend *Legend2 = new TLegend(0.6, 0.75, 0.9, 0.9);
            Legend2->SetBorderSize(0);
            Legend2->SetTextSize(0.025);
            Legend2->AddEntry(h_data_R_DB_ptbin[m][bin], "Data", "lep");
            Legend2->AddEntry(gausFit_data, "Gaussian Fit", "l");
            
            TLatex text2;
            text2.SetTextSize(0.03);
            text2.DrawLatexNDC(0.65, 0.7, Form("#mu = %.3f #pm %.3f", Mean_data, meanErr_data));
            text2.DrawLatexNDC(0.65, 0.65, Form("#sigma = %.3f #pm %.3f", Sigma_data, sigmaErr_data));
            text2.DrawLatexNDC(0.65, 0.6, Form("#chi^{2}/NDF = %.3f", chi2NDF_data));
            Legend2->Draw();
            
            DrawAtlasLabel();
            text2.DrawLatexNDC(0.2, 0.8, ptBinLabels[bin]);
            
            // Add mass cut information
            if (m == 0) {
                text2.DrawLatexNDC(0.2, 0.75, "No Mass Cut");
            } else {
                text2.DrawLatexNDC(0.2, 0.75, Form("Mass Cut: %d < m < %d GeV", (int)massCutLow[m], (int)massCutHigh[m]));
            }
            
            data_canvas->SaveAs(Form("../Plots/DB/PlotsRdbpTbin/%s/%s/Fit_R_DB_data_ptbin%d.pdf", wp_label, fileSystemMassCutLabels[m], bin));
            
            // Save fit parameters
            SaveFitParameters(Form("../Plots/DB/PlotsRdbpTbin/%s/%s/data_fit_parameters_ptbin%d.txt", wp_label, fileSystemMassCutLabels[m], bin),
                            Form("Data (%s, Mass Cut: %s)", ptBinLabels[bin], massCutLabels[m]),
                            Mean_data, meanErr_data, Sigma_data, sigmaErr_data,
                            amp_data, ampErr_data, chi_data, ndf_data);
                            
            // Add to summary file
            if (summaryFile.is_open()) {
                summaryFile << "Data (" << ptBinLabels[bin] << ")" << std::endl;
                summaryFile << "Mean: " << Mean_data << " ± " << meanErr_data << std::endl;
                summaryFile << "Sigma: " << Sigma_data << " ± " << sigmaErr_data << std::endl;
                summaryFile << "Chi2/NDF: " << chi2NDF_data << " (" << chi_data << "/" << ndf_data << ")" << std::endl << std::endl;
            }

            delete data_canvas;
            delete gausFit_data;
            delete Legend2;
            
            // Background fit
            TCanvas* bkg_canvas = new TCanvas(Form("bkg_canvas_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label), 
                                            Form("Background R_DB for %s (Mass Cut: %s, WP: %s)", ptBinLabels[bin], massCutLabels[m], wp_label), 800, 600);
            bkg_canvas->cd();
            
            TF1 *gausFit_bkg = new TF1(Form("gausFit_bkg_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label), "gaus", 0.9, 1.1);
            h_multijets_R_DB_ptbin[m][bin]->Fit(gausFit_bkg, "RQ");
            
            double chi_bkg = gausFit_bkg->GetChisquare();
            int ndf_bkg = gausFit_bkg->GetNDF();
            double chi2NDF_bkg = (ndf_bkg > 0) ? chi_bkg / ndf_bkg : 0.0;
            
            h_multijets_R_DB_ptbin[m][bin]->Draw("HIST");
            gausFit_bkg->SetLineColor(kViolet-4);
            gausFit_bkg->Draw("same");
            
            h_multijets_R_DB_ptbin[m][bin]->SetMaximum(h_multijets_R_DB_ptbin[m][bin]->GetMaximum() * 1.25);
            h_multijets_R_DB_ptbin[m][bin]->GetXaxis()->SetTitle("R_{DB}");
            h_multijets_R_DB_ptbin[m][bin]->GetYaxis()->SetTitle("Events");
            
            double amp_bkg = gausFit_bkg->GetParameter(0);
            double Mean_bkg = gausFit_bkg->GetParameter(1);
            double Sigma_bkg = gausFit_bkg->GetParameter(2);
            
            double ampErr_bkg = gausFit_bkg->GetParError(0);
            double meanErr_bkg = gausFit_bkg->GetParError(1);
            double sigmaErr_bkg = gausFit_bkg->GetParError(2);
            
            // Store results for comparison plots
            bkg_means[m][bin] = Mean_bkg;
            bkg_mean_errs[m][bin] = meanErr_bkg;
            bkg_sigmas[m][bin] = Sigma_bkg;
            bkg_sigma_errs[m][bin] = sigmaErr_bkg;
            
            TLegend *Legend3 = new TLegend(0.6, 0.75, 0.9, 0.9);
            Legend3->SetBorderSize(0);
            Legend3->SetTextSize(0.025);
            Legend3->AddEntry(h_multijets_R_DB_ptbin[m][bin], "#gamma + jets", "l");
            Legend3->AddEntry(gausFit_bkg, "Gaussian Fit", "l");
            
            TLatex text3;
            text3.SetTextSize(0.03);
            text3.DrawLatexNDC(0.65, 0.7, Form("#mu = %.3f #pm %.3f", Mean_bkg, meanErr_bkg));
            text3.DrawLatexNDC(0.65, 0.65, Form("#sigma = %.3f #pm %.3f", Sigma_bkg, sigmaErr_bkg));
            text3.DrawLatexNDC(0.65, 0.6, Form("#chi^{2}/NDF = %.3f", chi2NDF_bkg));
            Legend3->Draw();
            
            DrawAtlasLabel();
            text3.DrawLatexNDC(0.2, 0.8, ptBinLabels[bin]);
            
            // Add mass cut information
            if (m == 0) {
                text3.DrawLatexNDC(0.2, 0.75, "No Mass Cut");
            } else {
                text3.DrawLatexNDC(0.2, 0.75, Form("Mass Cut: %d < m < %d GeV", (int)massCutLow[m], (int)massCutHigh[m]));
            }
            
            bkg_canvas->SaveAs(Form("../Plots/DB/PlotsRdbpTbin/%s/%s/Fit_R_DB_bkg_ptbin%d.pdf", wp_label, fileSystemMassCutLabels[m], bin));
            
            // Save fit parameters
            SaveFitParameters(Form("../Plots/DB/PlotsRdbpTbin/%s/%s/bkg_fit_parameters_ptbin%d.txt", wp_label, fileSystemMassCutLabels[m], bin),
                            Form("Background: Multijets (%s, Mass Cut: %s)", ptBinLabels[bin], massCutLabels[m]),
                            Mean_bkg, meanErr_bkg, Sigma_bkg, sigmaErr_bkg,
                            amp_bkg, ampErr_bkg, chi_bkg, ndf_bkg);
                            
            // Add to summary file
            if (summaryFile.is_open()) {
                summaryFile << "Background: Multijets (" << ptBinLabels[bin] << ")" << std::endl;
                summaryFile << "Mean: " << Mean_bkg << " ± " << meanErr_bkg << std::endl;
                summaryFile << "Sigma: " << Sigma_bkg << " ± " << sigmaErr_bkg << std::endl;
                summaryFile << "Chi2/NDF: " << chi2NDF_bkg << " (" << chi_bkg << "/" << ndf_bkg << ")" << std::endl << std::endl;
            }
            delete bkg_canvas;
            delete gausFit_bkg;
            delete Legend3;
            
            // Now fit the combined MC (signal + background)
            TCanvas* combined_canvas = new TCanvas(Form("combined_canvas_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label), 
                                                 Form("Combined MC R_DB for %s (Mass Cut: %s, WP: %s)", ptBinLabels[bin], massCutLabels[m], wp_label), 800, 600);
            combined_canvas->cd();
            
            TF1 *gausFit_combined = new TF1(Form("gausFit_combined_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label), "gaus", 0.9, 1.1); //"gaus", 0.9, 1.1
            h_combined_MC_R_DB_ptbin[m][bin]->Fit(gausFit_combined, "RQ");
            
            double chi_combined = gausFit_combined->GetChisquare();
            int ndf_combined = gausFit_combined->GetNDF();
            double chi2NDF_combined = (ndf_combined > 0) ? chi_combined / ndf_combined : 0.0;
            
            h_combined_MC_R_DB_ptbin[m][bin]->Draw("HIST");
            gausFit_combined->SetLineColor(kViolet-4);
            gausFit_combined->Draw("same");
            
            h_combined_MC_R_DB_ptbin[m][bin]->SetMaximum(h_combined_MC_R_DB_ptbin[m][bin]->GetMaximum() * 1.25);
            h_combined_MC_R_DB_ptbin[m][bin]->GetXaxis()->SetTitle("R_{DB}");
            h_combined_MC_R_DB_ptbin[m][bin]->GetYaxis()->SetTitle("Events");
            
            double amp_combined = gausFit_combined->GetParameter(0);
            double Mean_combined = gausFit_combined->GetParameter(1);
            double Sigma_combined = gausFit_combined->GetParameter(2);
            
            double ampErr_combined = gausFit_combined->GetParError(0);
            double meanErr_combined = gausFit_combined->GetParError(1);
            double sigmaErr_combined = gausFit_combined->GetParError(2);
            
            // Store results for comparison plots
            combined_means[m][bin] = Mean_combined;
            combined_mean_errs[m][bin] = meanErr_combined;
            combined_sigmas[m][bin] = Sigma_combined;
            combined_sigma_errs[m][bin] = sigmaErr_combined;
            
            TLegend *combined_legend = new TLegend(0.6, 0.75, 0.9, 0.9);
            combined_legend->SetBorderSize(0);
            combined_legend->SetTextSize(0.025);
            combined_legend->AddEntry(h_combined_MC_R_DB_ptbin[m][bin], "Combined MC (Signal+Bkg)", "l");
            combined_legend->AddEntry(gausFit_combined, "Gaussian Fit", "l");
            
            TLatex combined_text;
            combined_text.SetTextSize(0.03);
            combined_text.DrawLatexNDC(0.65, 0.7, Form("#mu = %.3f #pm %.3f", Mean_combined, meanErr_combined));
            combined_text.DrawLatexNDC(0.65, 0.65, Form("#sigma = %.3f #pm %.3f", Sigma_combined, sigmaErr_combined));
            combined_text.DrawLatexNDC(0.65, 0.6, Form("#chi^{2}/NDF = %.3f", chi2NDF_combined));
            combined_legend->Draw();
            
            DrawAtlasLabel();
            combined_text.DrawLatexNDC(0.2, 0.8, ptBinLabels[bin]);
            
            // Add mass cut information
            if (m == 0) {
                combined_text.DrawLatexNDC(0.2, 0.75, "No Mass Cut");
            } else {
                combined_text.DrawLatexNDC(0.2, 0.75, Form("Mass Cut: %d < m < %d GeV", (int)massCutLow[m], (int)massCutHigh[m]));
            }
            
            combined_canvas->SaveAs(Form("../Plots/DB/PlotsRdbpTbin/%s/%s/Fit_R_DB_combined_MC_ptbin%d.pdf", wp_label, fileSystemMassCutLabels[m], bin));
            
            // Save fit parameters
            SaveFitParameters(Form("../Plots/DB/PlotsRdbpTbin/%s/%s/combined_MC_fit_parameters_ptbin%d.txt", wp_label, fileSystemMassCutLabels[m], bin),
                            Form("Combined MC (Signal+Bkg) (%s, Mass Cut: %s)", ptBinLabels[bin], massCutLabels[m]),
                            Mean_combined, meanErr_combined, Sigma_combined, sigmaErr_combined,
                            amp_combined, ampErr_combined, chi_combined, ndf_combined);
                            
            // Add to summary file
            if (summaryFile.is_open()) {
                summaryFile << "Combined MC (Signal+Bkg) (" << ptBinLabels[bin] << ")" << std::endl;
                summaryFile << "Mean: " << Mean_combined << " ± " << meanErr_combined << std::endl;
                summaryFile << "Sigma: " << Sigma_combined << " ± " << sigmaErr_combined << std::endl;
                summaryFile << "Chi2/NDF: " << chi2NDF_combined << " (" << chi_combined << "/" << ndf_combined << ")" << std::endl << std::endl;
            }
            delete combined_canvas;
            delete gausFit_combined;
            delete combined_legend;
            
            // Stack plot comparing all three distributions for this pT bin
            TCanvas* stack_canvas = new TCanvas(Form("stack_canvas_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label), 
                                              Form("R_DB Comparison for %s (Mass Cut: %s, WP: %s)", ptBinLabels[bin], massCutLabels[m], wp_label), 800, 600);
            stack_canvas->cd();
            
            // Clone histograms for comparison
            TH1F* h_data_clone = (TH1F*)h_data_R_DB_ptbin[m][bin]->Clone(Form("h_data_clone_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label));
            TH1F* h_mc_clone = (TH1F*)h_mc_R_DB_ptbin[m][bin]->Clone(Form("h_mc_clone_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label));
            TH1F* h_bkg_clone = (TH1F*)h_multijets_R_DB_ptbin[m][bin]->Clone(Form("h_bkg_clone_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label));
            TH1F* h_combined_clone = (TH1F*)h_combined_MC_R_DB_ptbin[m][bin]->Clone(Form("h_combined_clone_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label));
            
            // Style the clones
            h_data_clone->SetMarkerStyle(20);
            h_data_clone->SetMarkerColor(kBlack);
            h_data_clone->SetLineColor(kBlack);
            
            h_mc_clone->SetFillColor(kRed);
            h_mc_clone->SetLineColor(kRed);
            h_mc_clone->SetFillStyle(1001);
            
            h_bkg_clone->SetFillColor(kBlue);
            h_bkg_clone->SetLineColor(kBlue);
            h_bkg_clone->SetFillStyle(1001);
            
            h_combined_clone->SetFillColor(kGreen+2);
            h_combined_clone->SetLineColor(kGreen+2);
            h_combined_clone->SetFillStyle(1001);
            
            // Create stack
            THStack* stack = new THStack(Form("stack_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label), 
                                       Form("R_{DB} Distribution for %s;R_{DB};Events", ptBinLabels[bin]));
            stack->Add(h_bkg_clone);
            stack->Add(h_mc_clone);
            
            stack->Draw("HIST");
            h_combined_clone->Draw("HIST SAME");
            h_data_clone->Draw("E SAME");
            
            // Legend
            TLegend* stack_legend = new TLegend(0.6, 0.7, 0.89, 0.89);
            stack_legend->SetBorderSize(0);
            stack_legend->SetFillStyle(0);
            stack_legend->SetTextSize(0.025);
            stack_legend->AddEntry(h_data_clone, "Data", "lp");
            stack_legend->AddEntry(h_mc_clone, "Z (#rightarrow bb) + #gamma ", "f");
            stack_legend->AddEntry(h_bkg_clone, "#gamma + jets", "f");
            stack_legend->AddEntry(h_combined_clone, "Combined MC (Signal+Bkg)", "l");
            stack_legend->Draw();
            
            DrawAtlasLabel();
            TLatex text_stack;
            text_stack.SetTextSize(0.035);
            text_stack.DrawLatexNDC(0.2, 0.8, ptBinLabels[bin]);
            
            // Add mass cut information
            if (m == 0) {
                text_stack.DrawLatexNDC(0.2, 0.75, "No Mass Cut");
            } else {
                text_stack.DrawLatexNDC(0.2, 0.75, Form("Mass Cut: %d < m < %d GeV", (int)massCutLow[m], (int)massCutHigh[m]));
            }
            
            stack_canvas->SaveAs(Form("../Plots/DB/PlotsRdbpTbin/%s/%s/R_DB_stack_comparison_ptbin%d.pdf", wp_label, fileSystemMassCutLabels[m], bin));
            
            delete stack_canvas;
            delete h_data_clone;
            delete h_mc_clone;
            delete h_bkg_clone;
            delete h_combined_clone;
            delete stack;
            delete stack_legend;

            // Create canvas with ratio plot
            TCanvas* combinedstack_canvas = new TCanvas(Form("combinedstack_canvas_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label), 
                                                      Form("R_DB Comparison for %s (Mass Cut: %s, WP: %s)", ptBinLabels[bin], massCutLabels[m], wp_label), 800, 600);
            combinedstack_canvas->cd();
            
            // Configure pads
            TPad* pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
            pad1->SetBottomMargin(0.02);  // Reduced bottom margin
            pad1->Draw();
            pad1->cd();  // Switch to pad1
            
            // Clone histograms (same as original)
            TH1F* h_data_clone_comb = (TH1F*)h_data_R_DB_ptbin[m][bin]->Clone(Form("h_data_clone_comb_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label));
            TH1F* h_mc_clone_comb = (TH1F*)h_mc_R_DB_ptbin[m][bin]->Clone(Form("h_mc_clone_comb_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label));
            TH1F* h_bkg_clone_comb = (TH1F*)h_multijets_R_DB_ptbin[m][bin]->Clone(Form("h_bkg_clone_comb_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label));
            TH1F* h_combined_clone_comb = (TH1F*)h_combined_MC_R_DB_ptbin[m][bin]->Clone(Form("h_combined_clone_comb_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label));
            
            // Style clones (same as original)
            h_data_clone_comb->SetMarkerStyle(20);
            h_data_clone_comb->SetMarkerColor(kBlack);
            h_data_clone_comb->SetLineColor(kBlack);
            
            h_mc_clone_comb->SetFillColor(kRed);
            h_mc_clone_comb->SetLineColor(kRed);
            h_mc_clone_comb->SetFillStyle(1001);
            
            h_bkg_clone_comb->SetFillColor(kBlue);
            h_bkg_clone_comb->SetLineColor(kBlue);
            h_bkg_clone_comb->SetFillStyle(1001);
            
            h_combined_clone_comb->SetFillColor(kGreen+2);
            h_combined_clone_comb->SetLineColor(kGreen+2);
            h_combined_clone_comb->SetFillStyle(1001);
            
            // Create and draw stack
            THStack* stack_comb = new THStack(Form("stack_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label), 
                                            Form("R_{DB} Distribution for %s;;Events", ptBinLabels[bin]));  // Note empty x-axis title
            stack_comb->Add(h_bkg_clone_comb);
            stack_comb->Add(h_mc_clone_comb);
            stack_comb->Draw("HIST");
            h_data_clone_comb->Draw("EP SAME");
            
            // Add the fitted curves for data and combined MC
            // Need to recreate the TF1 objects here as they were deleted after previous canvases
            TF1 *gausFit_data_recreate = new TF1(Form("gausFit_data_recreate_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label), "gaus", 0.9, 1.1);
            gausFit_data_recreate->SetParameters(amp_data, Mean_data, Sigma_data); // Use stored parameters
            gausFit_data_recreate->SetLineColor(kViolet-4);
            gausFit_data_recreate->Draw("same");
            
            TF1 *gausFit_combined_recreate = new TF1(Form("gausFit_combined_recreate_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label), "gaus", 0.9, 1.1);
            gausFit_combined_recreate->SetParameters(amp_combined, Mean_combined, Sigma_combined); // Use stored parameters
            gausFit_combined_recreate->SetLineColor(kGreen-2);
            gausFit_combined_recreate->Draw("same");

            // Hide x-axis labels on top pad
            stack_comb->GetXaxis()->SetLabelSize(0);
            stack_comb->GetXaxis()->SetTitleSize(0);
            
            // Style y-axis
            stack_comb->GetYaxis()->SetTitleSize(0.05);
            stack_comb->GetYaxis()->SetTitle("Events");
            stack_comb->GetYaxis()->SetLabelSize(0.04);
            
            stack_comb->SetMaximum(stack_comb->GetMaximum() * 1.5); // Leave room for legend
            
            TH1F* h_mc_total = (TH1F*)stack_comb->GetStack()->Last(); // Get the total MC (signal+background)
            
            TGraphAsymmErrors* err_band = new TGraphAsymmErrors(h_mc_total);
            
            // Fill the error band points with appropriate errors
            for (int i = 0; i < err_band->GetN(); ++i) {
                double error = h_mc_total->GetBinError(i + 1);
                err_band->SetPointEYhigh(i, error);
                err_band->SetPointEYlow(i, error);
            }
            
            // Style the error band
            err_band->SetFillColorAlpha(kBlack, 0.35); // Semi-transparent
            err_band->SetFillStyle(3254);              // Hatch pattern
            err_band->SetLineWidth(0);                 // No outline
            err_band->Draw("2 SAME");                  // "2" draw mode for error band
            
            // Legend (adjusted position)
            TLegend* stack_legend_comb = new TLegend(0.6, 0.65, 0.89, 0.89);
            stack_legend_comb->SetBorderSize(0);          // No border
            stack_legend_comb->SetFillStyle(0);           // Transparent background
            stack_legend_comb->SetTextSize(0.025);        // Readable text size
            stack_legend_comb->AddEntry(h_data_clone_comb, "Data", "lep");
            stack_legend_comb->AddEntry(h_mc_clone_comb, "Z (#rightarrow bb) + #gamma ", "f");
            stack_legend_comb->AddEntry(h_bkg_clone_comb, "#gamma + jets", "f");
            stack_legend_comb->AddEntry(gausFit_data_recreate, "Data Fit", "l");
            stack_legend_comb->AddEntry(gausFit_combined_recreate, "Combined MC (Sig + Bkg) Fit", "l");
            stack_legend_comb->AddEntry(err_band, "MC Stat. Unc.", "f");
            stack_legend_comb->Draw();
            
            // Add ATLAS label to the plot
            DrawAtlasLabel();
            TLatex text_stack_comb;
            text_stack_comb.SetTextSize(0.035);
            text_stack_comb.DrawLatexNDC(0.2, 0.8, ptBinLabels[bin]);
            
            // Add mass cut information
            if (m == 0) {
                text_stack_comb.DrawLatexNDC(0.2, 0.75, "No Mass Cut");
            } else {
                text_stack_comb.DrawLatexNDC(0.2, 0.75, Form("Mass Cut: %d < m < %d GeV", (int)massCutLow[m], (int)massCutHigh[m]));
            }
            
            pad1->Modified();
            
            // Switch to ratio pad
            combinedstack_canvas->cd();
            TPad* pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
            pad2->SetTopMargin(0.05);
            pad2->SetBottomMargin(0.25);
            pad2->Draw();
            pad2->cd();
            
            // Create ratio histogram
            TH1F* h_ratio = (TH1F*)h_data_clone_comb->Clone(Form("h_ratio_masscut%s_bin%d_%s", massCutLabels[m], bin, wp_label));
            h_ratio->Divide(h_combined_clone_comb);
            
            // Style ratio plot
            h_ratio->SetTitle("");
            h_ratio->SetMarkerStyle(21);
            h_ratio->SetMarkerSize(1.0);
            h_ratio->SetLineColor(kBlack);
            h_ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
            h_ratio->GetXaxis()->SetTitle("R_{DB}");
            h_ratio->GetYaxis()->SetTitle("Data/MC"); // Data/MC: Data/(Sig+Bkg)
            h_ratio->GetYaxis()->SetTitleSize(0.12);
            h_ratio->GetYaxis()->SetTitleOffset(0.4);
            h_ratio->GetYaxis()->SetLabelSize(0.10);
            h_ratio->GetXaxis()->SetTitleSize(0.12);
            h_ratio->GetXaxis()->SetTitleOffset(1.0);
            h_ratio->GetXaxis()->SetLabelSize(0.10);
            h_ratio->GetYaxis()->SetNdivisions(505); // Fewer divisions for cleaner look
            h_ratio->GetYaxis()->CenterTitle();
            
            // Draw the ratio plot
            h_ratio->Draw("El same");
            
            // Create MC uncertainty band
            TGraphAsymmErrors* err_band_ratio = new TGraphAsymmErrors(h_combined_clone_comb);
            for (int i = 0; i < h_combined_clone_comb->GetNbinsX(); ++i) {
                double x = h_combined_clone_comb->GetBinCenter(i + 1);
                double content = h_combined_clone_comb->GetBinContent(i+1);
                double error = h_combined_clone_comb->GetBinError(i+1);
                
                // Calculate relative error for the ratio
                double rel_error = (content != 0) ? error / content : 0;
                
                err_band_ratio->SetPoint(i, x, 1.0);
                err_band_ratio->SetPointEYhigh(i, rel_error);
                err_band_ratio->SetPointEYlow(i, rel_error);
            }
            
            // Style the ratio error band
            err_band_ratio->SetFillColorAlpha(kBlack, 0.35);
            err_band_ratio->SetFillStyle(3354);  // Different style than top panel
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
            combinedstack_canvas->SaveAs(Form("../Plots/DB/PlotsRdbpTbin/%s/%s/R_DB_stack_ratio_ptbin%d.pdf", wp_label, fileSystemMassCutLabels[m], bin));
            
            // Cleanup for combinedstack_canvas
            delete pad1;
            delete pad2;
            delete combinedstack_canvas;
            delete h_data_clone_comb;
            delete h_mc_clone_comb;
            delete h_bkg_clone_comb;
            delete h_combined_clone_comb;
            delete stack_comb;
            delete gausFit_data_recreate; // Delete recreated TF1s
            delete gausFit_combined_recreate;
            delete err_band;
            delete stack_legend_comb;
            delete h_ratio;
            delete err_band_ratio;
            delete line;
        }
        
        // Close summary file
        if (summaryFile.is_open()) {
            summaryFile.close();
            std::cout << "Summary of all fit parameters saved to " << summaryFileName.Data() << std::endl;
        }
        
        // Create comparison plots of mean and sigma vs pT for this mass cut
        // Mean comparison
        TCanvas* mean_canvas = new TCanvas(Form("mean_canvas_masscut%s_%s", massCutLabels[m], wp_label), 
                                          Form("R_DB Mean vs pT (Mass Cut: %s, WP: %s)", massCutLabels[m], wp_label), 800, 600);
        mean_canvas->cd();
        
        // Create graph for each sample
        TGraphErrors* g_mc_mean = new TGraphErrors(nPtBins, ptBinCenters, mc_means[m], nullptr, mc_mean_errs[m]);
        TGraphErrors* g_data_mean = new TGraphErrors(nPtBins, ptBinCenters, data_means[m], nullptr, data_mean_errs[m]);
        TGraphErrors* g_bkg_mean = new TGraphErrors(nPtBins, ptBinCenters, bkg_means[m], nullptr, bkg_mean_errs[m]);
        TGraphErrors* g_combined_mean = new TGraphErrors(nPtBins, ptBinCenters, combined_means[m], nullptr, combined_mean_errs[m]);
        
        // Style graphs
        g_mc_mean->SetMarkerStyle(20);
        g_mc_mean->SetMarkerColor(kRed);
        g_mc_mean->SetLineColor(kRed);
        g_mc_mean->SetMarkerSize(1.2);
        
        g_data_mean->SetMarkerStyle(21);
        g_data_mean->SetMarkerColor(kBlack);
        g_data_mean->SetLineColor(kBlack);
        g_data_mean->SetMarkerSize(1.2);
        
        g_bkg_mean->SetMarkerStyle(22);
        g_bkg_mean->SetMarkerColor(kBlue);
        g_bkg_mean->SetLineColor(kBlue);
        g_bkg_mean->SetMarkerSize(1.2);
        
        g_combined_mean->SetMarkerStyle(23);
        g_combined_mean->SetMarkerColor(kGreen+2);
        g_combined_mean->SetLineColor(kGreen+2);
        g_combined_mean->SetMarkerSize(1.2);
        
        // Draw graphs
        g_mc_mean->SetTitle("R_{DB} Mean vs p_{T};p_{T} [GeV];#LT R_{DB}#GT");
        g_mc_mean->GetYaxis()->SetRangeUser(0.95, 1.05); // Adjust range as needed
        g_mc_mean->Draw("APL");
        g_data_mean->Draw("PL SAME");
        g_bkg_mean->Draw("PL SAME");
        g_combined_mean->Draw("PL SAME");
        
        // Add legend
        TLegend* mean_legend = new TLegend(0.6, 0.7, 0.89, 0.89);
        mean_legend->SetBorderSize(0);
        mean_legend->AddEntry(g_mc_mean, "Z (#rightarrow bb) + #gamma ", "pl");
        mean_legend->AddEntry(g_data_mean, "Data", "pl");
        mean_legend->AddEntry(g_bkg_mean, "#gamma + jets", "pl");
        mean_legend->AddEntry(g_combined_mean, "Combined MC (Signal+Bkg)", "pl");
        mean_legend->Draw();
        
        DrawAtlasLabel();
        
        // Add mass cut information
        TLatex mean_text;
        mean_text.SetTextSize(0.035);
        if (m == 0) {
            mean_text.DrawLatexNDC(0.2, 0.8, "No Mass Cut");
        } else {
            mean_text.DrawLatexNDC(0.2, 0.8, Form("Mass Cut: %d < m < %d GeV", (int)massCutLow[m], (int)massCutHigh[m]));
        }
        
        mean_canvas->SaveAs(Form("../Plots/DB/PlotsRdbpTbin/%s/%s/R_DB_mean_vs_pT.pdf", wp_label, fileSystemMassCutLabels[m]));
        
        delete mean_canvas;
        delete g_mc_mean;
        delete g_data_mean;
        delete g_bkg_mean;
        delete g_combined_mean;
        delete mean_legend;

        // Sigma comparison
        TCanvas* sigma_canvas = new TCanvas(Form("sigma_canvas_masscut%s_%s", massCutLabels[m], wp_label), 
                                           Form("R_DB Sigma vs pT (Mass Cut: %s, WP: %s)", massCutLabels[m], wp_label), 800, 600);
        sigma_canvas->cd();
        
        // Create graphs
        TGraphErrors* g_mc_sigma = new TGraphErrors(nPtBins, ptBinCenters, mc_sigmas[m], nullptr, mc_sigma_errs[m]);
        TGraphErrors* g_data_sigma = new TGraphErrors(nPtBins, ptBinCenters, data_sigmas[m], nullptr, data_sigma_errs[m]);
        TGraphErrors* g_bkg_sigma = new TGraphErrors(nPtBins, ptBinCenters, bkg_sigmas[m], nullptr, bkg_sigma_errs[m]);
        TGraphErrors* g_combined_sigma = new TGraphErrors(nPtBins, ptBinCenters, combined_sigmas[m], nullptr, combined_sigma_errs[m]);
        
        // Style graphs
        g_mc_sigma->SetMarkerStyle(20);
        g_mc_sigma->SetMarkerColor(kRed);
        g_mc_sigma->SetLineColor(kRed);
        g_mc_sigma->SetMarkerSize(1.2);
        
        g_data_sigma->SetMarkerStyle(21);
        g_data_sigma->SetMarkerColor(kBlack);
        g_data_sigma->SetLineColor(kBlack);
        g_data_sigma->SetMarkerSize(1.2);
        
        g_bkg_sigma->SetMarkerStyle(22);
        g_bkg_sigma->SetMarkerColor(kBlue);
        g_bkg_sigma->SetLineColor(kBlue);
        g_bkg_sigma->SetMarkerSize(1.2);
        
        g_combined_sigma->SetMarkerStyle(23);
        g_combined_sigma->SetMarkerColor(kGreen+2);
        g_combined_sigma->SetLineColor(kGreen+2);
        g_combined_sigma->SetMarkerSize(1.2);
        
        // Draw graphs
        g_mc_sigma->SetTitle("R_{DB} Width vs p_{T};p_{T} [GeV];#sigma (R_{DB})");
        g_mc_sigma->GetYaxis()->SetRangeUser(0.0, 0.2); // Adjust range as needed
        g_mc_sigma->Draw("APL");
        g_data_sigma->Draw("PL SAME");
        g_bkg_sigma->Draw("PL SAME");
        g_combined_sigma->Draw("PL SAME");
        
        // Add legend
        TLegend* sigma_legend = new TLegend(0.6, 0.7, 0.89, 0.89);
        sigma_legend->SetBorderSize(0);
        sigma_legend->AddEntry(g_mc_sigma, "Z (#rightarrow bb) + #gamma ", "pl");
        sigma_legend->AddEntry(g_data_sigma, "Data", "pl");
        sigma_legend->AddEntry(g_bkg_sigma, "#gamma + jets", "pl");
        sigma_legend->AddEntry(g_combined_sigma, "Combined MC (Signal+Bkg)", "pl");
        sigma_legend->Draw();
        
        DrawAtlasLabel();
        
        // Add mass cut information
        TLatex sigma_text;
        sigma_text.SetTextSize(0.035);
        if (m == 0) {
            sigma_text.DrawLatexNDC(0.2, 0.85, "No Mass Cut");
        } else {
            sigma_text.DrawLatexNDC(0.2, 0.85, Form("Mass Cut: %d < m < %d GeV", (int)massCutLow[m], (int)massCutHigh[m]));
        }
        
        sigma_canvas->SaveAs(Form("../Plots/DB/PlotsRdbpTbin/%s/%s/R_DB_sigma_vs_pT.pdf", wp_label, fileSystemMassCutLabels[m]));
        
        delete sigma_canvas;
        delete g_mc_sigma;
        delete g_data_sigma;
        delete g_bkg_sigma;
        delete g_combined_sigma;
        delete sigma_legend;

        // Create ratio plots (data/MC and background/MC)
        TCanvas* ratio_canvas = new TCanvas(Form("ratio_canvas_masscut%s_%s", massCutLabels[m], wp_label), 
                                           Form("R_DB Mean Ratio vs pT (Mass Cut: %s, WP: %s)", massCutLabels[m], wp_label), 800, 600);
        ratio_canvas->cd();
        
        // Calculate ratios
        double data_mc_ratio[nPtBins], data_mc_ratio_err[nPtBins];
        double bkg_mc_ratio[nPtBins], bkg_mc_ratio_err[nPtBins];
        double data_combined_ratio[nPtBins], data_combined_ratio_err[nPtBins];
        
        for (int i = 0; i < nPtBins; i++) {
            // Data/MC ratio
            if (mc_means[m][i] != 0) {
                data_mc_ratio[i] = data_means[m][i] / mc_means[m][i];
                // Error propagation for ratio
                data_mc_ratio_err[i] = data_mc_ratio[i] * std::sqrt(
                    std::pow(data_mean_errs[m][i]/data_means[m][i], 2) + 
                    std::pow(mc_mean_errs[m][i]/mc_means[m][i], 2));
            } else {
                data_mc_ratio[i] = 0.0;
                data_mc_ratio_err[i] = 0.0;
            }
                
            // Background/MC ratio
            if (mc_means[m][i] != 0) {
                bkg_mc_ratio[i] = bkg_means[m][i] / mc_means[m][i];
                // Error propagation for ratio
                bkg_mc_ratio_err[i] = bkg_mc_ratio[i] * std::sqrt(
                    std::pow(bkg_mean_errs[m][i]/bkg_means[m][i], 2) + 
                    std::pow(mc_mean_errs[m][i]/mc_means[m][i], 2));
            } else {
                bkg_mc_ratio[i] = 0.0;
                bkg_mc_ratio_err[i] = 0.0;
            }
                
            // Data/Combined MC ratio
            if (combined_means[m][i] != 0) {
                data_combined_ratio[i] = data_means[m][i] / combined_means[m][i];
                // Error propagation for ratio
                data_combined_ratio_err[i] = data_combined_ratio[i] * std::sqrt(
                    std::pow(data_mean_errs[m][i]/data_means[m][i], 2) + 
                    std::pow(combined_mean_errs[m][i]/combined_means[m][i], 2));
            } else {
                data_combined_ratio[i] = 0.0;
                data_combined_ratio_err[i] = 0.0;
            }
        }
        
        // Create graphs
        TGraphErrors* g_data_mc_ratio = new TGraphErrors(nPtBins, ptBinCenters, data_mc_ratio, nullptr, data_mc_ratio_err);
        TGraphErrors* g_bkg_mc_ratio = new TGraphErrors(nPtBins, ptBinCenters, bkg_mc_ratio, nullptr, bkg_mc_ratio_err);
        TGraphErrors* g_data_combined_ratio = new TGraphErrors(nPtBins, ptBinCenters, data_combined_ratio, nullptr, data_combined_ratio_err);
        
        // Style
        g_data_mc_ratio->SetMarkerStyle(21);
        g_data_mc_ratio->SetMarkerColor(kBlack);
        g_data_mc_ratio->SetLineColor(kBlack);
        g_data_mc_ratio->SetMarkerSize(1.2);
        
        g_bkg_mc_ratio->SetMarkerStyle(22);
        g_bkg_mc_ratio->SetMarkerColor(kBlue);
        g_bkg_mc_ratio->SetLineColor(kBlue);
        g_bkg_mc_ratio->SetMarkerSize(1.2);
        
        g_data_combined_ratio->SetMarkerStyle(23);
        g_data_combined_ratio->SetMarkerColor(kGreen+2);
        g_data_combined_ratio->SetLineColor(kGreen+2);
        g_data_combined_ratio->SetMarkerSize(1.2);
        
        // Draw
        g_data_mc_ratio->SetTitle("Ratio of R_{DB} Mean to MC;p_{T} [GeV];Ratio to MC");
        g_data_mc_ratio->GetYaxis()->SetRangeUser(0.95, 1.05); // Adjust range
        g_data_mc_ratio->Draw("APL");
        g_bkg_mc_ratio->Draw("PL SAME");
        g_data_combined_ratio->Draw("PL SAME");

        // First make sure we get the actual axis range
        double xmin = g_data_mc_ratio->GetXaxis()->GetXmin();
        double xmax = g_data_mc_ratio->GetXaxis()->GetXmax();
        
        // Add a reference line at 1
        TLine* line_ratio = new TLine(xmin, 1.0, xmax, 1.0);
        line_ratio->SetLineStyle(2);
        line_ratio->SetLineColor(kRed);
        line_ratio->Draw();
        
        // Legend
        TLegend* ratio_legend = new TLegend(0.6, 0.7, 0.89, 0.89);
        ratio_legend->SetBorderSize(0);
        ratio_legend->AddEntry(g_data_mc_ratio, "Data/Signal MC", "pl");
        ratio_legend->AddEntry(g_bkg_mc_ratio, "Multijets/Signal MC", "pl");
        ratio_legend->AddEntry(g_data_combined_ratio, "Data/Combined MC", "pl");
        ratio_legend->AddEntry(line_ratio, "Ratio = 1", "l");
        ratio_legend->Draw();
        
        DrawAtlasLabel();
        
        // Add mass cut information
        TLatex ratio_text;
        ratio_text.SetTextSize(0.035);
        if (m == 0) {
            ratio_text.DrawLatexNDC(0.2, 0.8, "No Mass Cut");
        } else {
            ratio_text.DrawLatexNDC(0.2, 0.8, Form("Mass Cut: %d < m < %d GeV", (int)massCutLow[m], (int)massCutHigh[m]));
        }
        
        ratio_canvas->SaveAs(Form("../Plots/DB/PlotsRdbpTbin/%s/%s/R_DB_mean_ratio_vs_pT.pdf", wp_label, fileSystemMassCutLabels[m]));
        
        delete ratio_canvas;
        delete g_data_mc_ratio;
        delete g_bkg_mc_ratio;
        delete g_data_combined_ratio;
        delete line_ratio;
        delete ratio_legend;

        // Create data/combined MC comparison plot
        TCanvas* data_combined_canvas = new TCanvas(Form("data_combined_canvas_masscut%s_%s", massCutLabels[m], wp_label), 
                                                  Form("Data vs Combined MC (Mass Cut: %s, WP: %s)", massCutLabels[m], wp_label), 800, 600);
        data_combined_canvas->cd();
        
        // Draw only data/combined MC ratio
        TGraphErrors* g_data_combined_ratio_only = new TGraphErrors(nPtBins, ptBinCenters, data_combined_ratio, nullptr, data_combined_ratio_err);
        g_data_combined_ratio_only->SetTitle("Ratio of Data to Combined MC;p_{T} [GeV];Data/Combined MC");
        g_data_combined_ratio_only->GetYaxis()->SetRangeUser(0.95, 1.05); // Adjust range
        g_data_combined_ratio_only->SetMarkerStyle(23);
        g_data_combined_ratio_only->SetMarkerColor(kGreen+2);
        g_data_combined_ratio_only->SetLineColor(kGreen+2);
        g_data_combined_ratio_only->SetMarkerSize(1.2);
        g_data_combined_ratio_only->Draw("APL");
        
        // Add a reference line at 1
        TLine* line2 = new TLine(ptBins[0], 1.0, ptBins[nPtBins], 1.0);
        line2->SetLineStyle(2);
        line2->SetLineColor(kRed);
        line2->Draw();
        
        // Legend
        TLegend* data_combined_legend = new TLegend(0.6, 0.7, 0.89, 0.89);
        data_combined_legend->SetBorderSize(0);
        data_combined_legend->AddEntry(g_data_combined_ratio_only, "Data/Combined MC", "pl");
        data_combined_legend->AddEntry(line2, "MC Reference", "l");
        data_combined_legend->Draw();
        
        DrawAtlasLabel();
        
        // Add mass cut information
        TLatex data_comb_text;
        data_comb_text.SetTextSize(0.035);
        if (m == 0) {
            data_comb_text.DrawLatexNDC(0.2, 0.8, "No Mass Cut");
        } else {
            data_comb_text.DrawLatexNDC(0.2, 0.8, Form("Mass Cut: %d < m < %d GeV", (int)massCutLow[m], (int)massCutHigh[m]));
        }
        
        data_combined_canvas->SaveAs(Form("../Plots/DB/PlotsRdbpTbin/%s/%s/R_DB_data_combined_ratio_vs_pT.pdf", wp_label, fileSystemMassCutLabels[m]));

        delete data_combined_canvas;
        delete g_data_combined_ratio_only;
        delete line2;
        delete data_combined_legend;
    }
    
    // Now create comparison plots across different mass cuts
    // For each pT bin, compare results across mass cuts
    for (int bin = 0; bin < nPtBins; bin++) {
        // Create comparison of means across mass cuts for this pT bin
        TCanvas* mass_mean_canvas = new TCanvas(Form("mass_mean_canvas_ptbin%d_%s", bin, wp_label), 
                                                Form("R_DB Mean vs Mass Cut (pT bin: %s, WP: %s)", ptBinLabels[bin], wp_label), 800, 600);
        mass_mean_canvas->cd();
        
        double massCutX[nMassCuts];
        for (int m = 0; m < nMassCuts; m++) {
            massCutX[m] = m; // Use indices as x-coordinates
        }
        
        // Extract means for this pT bin across different mass cuts
        double mc_mean_values[nMassCuts], mc_mean_errors[nMassCuts];
        double data_mean_values[nMassCuts], data_mean_errors[nMassCuts];
        double bkg_mean_values[nMassCuts], bkg_mean_errors[nMassCuts];
        double combined_mean_values[nMassCuts], combined_mean_errors[nMassCuts];
        
        for (int m = 0; m < nMassCuts; m++) {
            mc_mean_values[m] = mc_means[m][bin];
            mc_mean_errors[m] = mc_mean_errs[m][bin];
            data_mean_values[m] = data_means[m][bin];
            data_mean_errors[m] = data_mean_errs[m][bin];
            bkg_mean_values[m] = bkg_means[m][bin];
            bkg_mean_errors[m] = bkg_mean_errs[m][bin];
            combined_mean_values[m] = combined_means[m][bin];
            combined_mean_errors[m] = combined_mean_errs[m][bin];
        }
        
        // Create graphs
        TGraphErrors* g_mc_mean_mass = new TGraphErrors(nMassCuts, massCutX, mc_mean_values, nullptr, mc_mean_errors);
        TGraphErrors* g_data_mean_mass = new TGraphErrors(nMassCuts, massCutX, data_mean_values, nullptr, data_mean_errors);
        TGraphErrors* g_bkg_mean_mass = new TGraphErrors(nMassCuts, massCutX, bkg_mean_values, nullptr, bkg_mean_errors); 
        TGraphErrors* g_combined_mean_mass = new TGraphErrors(nMassCuts, massCutX, combined_mean_values, nullptr, combined_mean_errors);
        
        // Style graphs
        g_mc_mean_mass->SetMarkerStyle(20);
        g_mc_mean_mass->SetMarkerColor(kRed);
        g_mc_mean_mass->SetLineColor(kRed);
        g_mc_mean_mass->SetMarkerSize(1.2);
        
        g_data_mean_mass->SetMarkerStyle(21);
        g_data_mean_mass->SetMarkerColor(kBlack);
        g_data_mean_mass->SetLineColor(kBlack);
        g_data_mean_mass->SetMarkerSize(1.2);
        
        g_bkg_mean_mass->SetMarkerStyle(22); 
        g_bkg_mean_mass->SetMarkerColor(kBlue); 
        g_bkg_mean_mass->SetLineColor(kBlue); 
        g_bkg_mean_mass->SetMarkerSize(1.2); 
        
        g_combined_mean_mass->SetMarkerStyle(23);
        g_combined_mean_mass->SetMarkerColor(kGreen+2);
        g_combined_mean_mass->SetLineColor(kGreen+2);
        g_combined_mean_mass->SetMarkerSize(1.2);
        
        // Draw graphs
        g_mc_mean_mass->SetTitle(Form("R_{DB} Mean vs Mass Cut for %s;Mass Cut;#LT R_{DB}#GT", ptBinLabels[bin]));
        g_mc_mean_mass->GetYaxis()->SetRangeUser(0.95, 1.05); // Adjust range as needed
        g_mc_mean_mass->GetXaxis()->SetLimits(-0.5, nMassCuts-0.5);
        g_mc_mean_mass->Draw("APL");
        g_data_mean_mass->Draw("PL SAME");
        g_bkg_mean_mass->Draw("PL SAME"); 
        g_combined_mean_mass->Draw("PL SAME");
        
        // Add custom x-axis labels
        TAxis *xaxis = g_mc_mean_mass->GetXaxis();
        for (int m = 0; m < nMassCuts; m++) {
            xaxis->SetBinLabel(xaxis->FindBin(m), massCutLabels[m]);
        }
        xaxis->SetLabelSize(0.05);
        
        // Legend
        TLegend* mass_mean_legend = new TLegend(0.6, 0.75, 0.89, 0.89);
        mass_mean_legend->SetBorderSize(0);
        mass_mean_legend->AddEntry(g_mc_mean_mass, "Z (#rightarrow bb) + #gamma ", "pl");
        mass_mean_legend->AddEntry(g_data_mean_mass, "Data", "pl");
        mass_mean_legend->AddEntry(g_bkg_mean_mass, "#gamma + jets", "pl");
        mass_mean_legend->AddEntry(g_combined_mean_mass, "Combined MC (Signal+Bkg)", "pl");
        mass_mean_legend->Draw();
        
        DrawAtlasLabel();
        TLatex mass_mean_text;
        mass_mean_text.SetTextSize(0.035);
        mass_mean_text.DrawLatexNDC(0.2, 0.8, ptBinLabels[bin]);
        
        mass_mean_canvas->SetGridx();
        mass_mean_canvas->SaveAs(Form("../Plots/DB/PlotsRdbpTbin/%s/R_DB_mean_vs_masscut_ptbin%d.pdf", wp_label, bin));
        
        delete mass_mean_canvas;
        delete g_mc_mean_mass;
        delete g_data_mean_mass;
        delete g_bkg_mean_mass;
        delete g_combined_mean_mass;
        delete mass_mean_legend;

        // Also create comparison of sigmas across mass cuts for this pT bin
        TCanvas* mass_sigma_canvas = new TCanvas(Form("mass_sigma_canvas_ptbin%d_%s", bin, wp_label), 
                                                Form("R_DB Sigma vs Mass Cut (pT bin: %s, WP: %s)", ptBinLabels[bin], wp_label), 800, 600);
        mass_sigma_canvas->cd();
        
        // Extract sigmas for this pT bin across different mass cuts
        double mc_sigma_values[nMassCuts], mc_sigma_errors[nMassCuts];
        double data_sigma_values[nMassCuts], data_sigma_errors[nMassCuts];
        double bkg_sigma_values[nMassCuts], bkg_sigma_errors[nMassCuts]; 
        double combined_sigma_values[nMassCuts], combined_sigma_errors[nMassCuts];
        
        for (int m = 0; m < nMassCuts; m++) {
            mc_sigma_values[m] = mc_sigmas[m][bin];
            mc_sigma_errors[m] = mc_sigma_errs[m][bin];
            data_sigma_values[m] = data_sigmas[m][bin];
            data_sigma_errors[m] = data_sigma_errs[m][bin];
            bkg_sigma_values[m] = bkg_sigmas[m][bin]; 
            bkg_sigma_errors[m] = bkg_sigma_errs[m][bin];
            combined_sigma_values[m] = combined_sigmas[m][bin];
            combined_sigma_errors[m] = combined_sigma_errs[m][bin];
        }
        
        // Create graphs
        TGraphErrors* g_mc_sigma_mass = new TGraphErrors(nMassCuts, massCutX, mc_sigma_values, nullptr, mc_sigma_errors);
        TGraphErrors* g_data_sigma_mass = new TGraphErrors(nMassCuts, massCutX, data_sigma_values, nullptr, data_sigma_errors);
        TGraphErrors* g_bkg_sigma_mass = new TGraphErrors(nMassCuts, massCutX, bkg_sigma_values, nullptr, bkg_sigma_errors);
        TGraphErrors* g_combined_sigma_mass = new TGraphErrors(nMassCuts, massCutX, combined_sigma_values, nullptr, combined_sigma_errors);
        
        // Style graphs
        g_mc_sigma_mass->SetMarkerStyle(20);
        g_mc_sigma_mass->SetMarkerColor(kRed);
        g_mc_sigma_mass->SetLineColor(kRed);
        g_mc_sigma_mass->SetMarkerSize(1.2);
        
        g_data_sigma_mass->SetMarkerStyle(21);
        g_data_sigma_mass->SetMarkerColor(kBlack);
        g_data_sigma_mass->SetLineColor(kBlack);
        g_data_sigma_mass->SetMarkerSize(1.2);
        
        g_bkg_sigma_mass->SetMarkerStyle(22); 
        g_bkg_sigma_mass->SetMarkerColor(kBlue); 
        g_bkg_sigma_mass->SetLineColor(kBlue); 
        g_bkg_sigma_mass->SetMarkerSize(1.2); 
        
        g_combined_sigma_mass->SetMarkerStyle(23);
        g_combined_sigma_mass->SetMarkerColor(kGreen+2);
        g_combined_sigma_mass->SetLineColor(kGreen+2);
        g_combined_sigma_mass->SetMarkerSize(1.2);
        
        // Draw graphs
        g_mc_sigma_mass->SetTitle(Form("R_{DB} Width vs Mass Cut for %s;Mass Cut;#sigma (R_{DB})", ptBinLabels[bin]));
        g_mc_sigma_mass->GetYaxis()->SetRangeUser(0.0, 0.2); // Adjust range as needed
        g_mc_sigma_mass->GetXaxis()->SetLimits(-0.5, nMassCuts-0.5);
        g_mc_sigma_mass->Draw("APL");
        g_data_sigma_mass->Draw("PL SAME");
        g_bkg_sigma_mass->Draw("PL SAME"); 
        g_combined_sigma_mass->Draw("PL SAME");
        
        // Add custom x-axis labels
        TAxis *xaxis_sigma = g_mc_sigma_mass->GetXaxis();
        for (int m = 0; m < nMassCuts; m++) {
            xaxis_sigma->SetBinLabel(xaxis_sigma->FindBin(m), massCutLabels[m]);
        }
        xaxis_sigma->SetLabelSize(0.05);
        
        // Legend
        TLegend* mass_sigma_legend = new TLegend(0.6, 0.7, 0.89, 0.89);
        mass_sigma_legend->SetBorderSize(0);
        mass_sigma_legend->AddEntry(g_mc_sigma_mass, "Z (#rightarrow bb) + #gamma ", "pl");
        mass_sigma_legend->AddEntry(g_data_sigma_mass, "Data", "pl");
        mass_sigma_legend->AddEntry(g_bkg_sigma_mass, "#gamma + jets", "pl"); 
        mass_sigma_legend->AddEntry(g_combined_sigma_mass, "Combined MC (Signal+Bkg)", "pl");
        mass_sigma_legend->Draw();
        
        DrawAtlasLabel();
        TLatex mass_sigma_text;
        mass_sigma_text.SetTextSize(0.035);
        mass_sigma_text.DrawLatexNDC(0.2, 0.8, ptBinLabels[bin]);
        
        mass_sigma_canvas->SetGridx();
        mass_sigma_canvas->SaveAs(Form("../Plots/DB/PlotsRdbpTbin/%s/R_DB_sigma_vs_masscut_ptbin%d.pdf", wp_label, bin));
        
        delete mass_sigma_canvas;
        delete g_mc_sigma_mass;
        delete g_data_sigma_mass;
        delete g_bkg_sigma_mass;
        delete g_combined_sigma_mass;
        delete mass_sigma_legend;

        // Create data/MC ratio comparison across mass cuts
        TCanvas* mass_ratio_canvas = new TCanvas(Form("mass_ratio_canvas_ptbin%d_%s", bin, wp_label), 
                                               Form("Data/MC Ratio vs Mass Cut (pT bin: %s, WP: %s)", ptBinLabels[bin], wp_label), 800, 600);
        mass_ratio_canvas->cd();
        
        // Calculate ratios for each mass cut
        double data_mc_ratio_values[nMassCuts], data_mc_ratio_errors[nMassCuts];
        double data_combined_ratio_values[nMassCuts], data_combined_ratio_errors[nMassCuts];
        double bkg_mc_ratio_values[nMassCuts], bkg_mc_ratio_errors[nMassCuts];
        
        for (int m = 0; m < nMassCuts; m++) {
            // Data/MC ratio
            if (mc_mean_values[m] != 0) {
                data_mc_ratio_values[m] = data_mean_values[m] / mc_mean_values[m];
                data_mc_ratio_errors[m] = data_mc_ratio_values[m] * std::sqrt(
                    std::pow(data_mean_errors[m]/data_mean_values[m], 2) + 
                    std::pow(mc_mean_errors[m]/mc_mean_values[m], 2));
            } else {
                data_mc_ratio_values[m] = 0.0;
                data_mc_ratio_errors[m] = 0.0;
            }
                
            // Data/Combined ratio
            if (combined_mean_values[m] != 0) {
                data_combined_ratio_values[m] = data_mean_values[m] / combined_mean_values[m];
                data_combined_ratio_errors[m] = data_combined_ratio_values[m] * std::sqrt(
                    std::pow(data_mean_errors[m]/data_mean_values[m], 2) + 
                    std::pow(combined_mean_errors[m]/combined_mean_values[m], 2));
            } else {
                data_combined_ratio_values[m] = 0.0;
                data_combined_ratio_errors[m] = 0.0;
            }
                
            // Background/MC ratio
            if (mc_mean_values[m] != 0) {
                bkg_mc_ratio_values[m] = bkg_mean_values[m] / mc_mean_values[m];
                bkg_mc_ratio_errors[m] = bkg_mc_ratio_values[m] * std::sqrt(
                    std::pow(bkg_mean_errors[m]/bkg_mean_values[m], 2) + 
                    std::pow(mc_mean_errors[m]/mc_mean_values[m], 2));
            } else {
                bkg_mc_ratio_values[m] = 0.0;
                bkg_mc_ratio_errors[m] = 0.0;
            }
        }
        
        // Create graphs
        TGraphErrors* g_data_mc_ratio_mass = new TGraphErrors(nMassCuts, massCutX, data_mc_ratio_values, nullptr, data_mc_ratio_errors);
        TGraphErrors* g_data_combined_ratio_mass = new TGraphErrors(nMassCuts, massCutX, data_combined_ratio_values, nullptr, data_combined_ratio_errors);
        TGraphErrors* g_bkg_mc_ratio_mass = new TGraphErrors(nMassCuts, massCutX, bkg_mc_ratio_values, nullptr, bkg_mc_ratio_errors);
        
        // Style graphs
        g_data_mc_ratio_mass->SetMarkerStyle(20);
        g_data_mc_ratio_mass->SetMarkerColor(kRed);
        g_data_mc_ratio_mass->SetLineColor(kRed);
        g_data_mc_ratio_mass->SetMarkerSize(1.2);
        
        g_data_combined_ratio_mass->SetMarkerStyle(21);
        g_data_combined_ratio_mass->SetMarkerColor(kBlack);
        g_data_combined_ratio_mass->SetLineColor(kBlack);
        g_data_combined_ratio_mass->SetMarkerSize(1.2);
        
        g_bkg_mc_ratio_mass->SetMarkerStyle(22);
        g_bkg_mc_ratio_mass->SetMarkerColor(kBlue);
        g_bkg_mc_ratio_mass->SetLineColor(kBlue);
        g_bkg_mc_ratio_mass->SetMarkerSize(1.2);
        
        // Draw graphs
        g_data_mc_ratio_mass->SetTitle(Form("Ratio vs Mass Cut for %s;Mass Cut;Ratio", ptBinLabels[bin]));
        g_data_mc_ratio_mass->GetYaxis()->SetRangeUser(0.95, 1.05); // Adjust range as needed
        g_data_mc_ratio_mass->GetXaxis()->SetLimits(-0.5, nMassCuts-0.5);
        g_data_mc_ratio_mass->Draw("APL");
        g_data_combined_ratio_mass->Draw("PL SAME");
        g_bkg_mc_ratio_mass->Draw("PL SAME");
        
        // Add reference line at 1
        TLine* ratio_line = new TLine(-0.5, 1.0, nMassCuts-0.5, 1.0);
        ratio_line->SetLineStyle(2);
        ratio_line->SetLineColor(kRed);
        ratio_line->Draw();
        
        // Add custom x-axis labels
        TAxis *xaxis_ratio = g_data_mc_ratio_mass->GetXaxis();
        for (int m = 0; m < nMassCuts; m++) {
            xaxis_ratio->SetBinLabel(xaxis_ratio->FindBin(m), massCutLabels[m]);
        }
        xaxis_ratio->SetLabelSize(0.05);
        
        // Legend
        TLegend* mass_ratio_legend = new TLegend(0.6, 0.7, 0.89, 0.89);
        mass_ratio_legend->SetBorderSize(0);
        mass_ratio_legend->AddEntry(g_data_mc_ratio_mass, "Data/Signal MC", "pl");
        mass_ratio_legend->AddEntry(g_data_combined_ratio_mass, "Data/Combined MC", "pl");
        mass_ratio_legend->AddEntry(g_bkg_mc_ratio_mass, "Multijets/Signal MC", "pl");
        mass_ratio_legend->AddEntry(ratio_line, "Ratio = 1", "l");
        mass_ratio_legend->Draw();
        
        DrawAtlasLabel();
        TLatex mass_ratio_text;
        mass_ratio_text.SetTextSize(0.035);
        mass_ratio_text.DrawLatexNDC(0.2, 0.8, ptBinLabels[bin]);
        
        mass_ratio_canvas->SetGridx();
        mass_ratio_canvas->SaveAs(Form("../Plots/DB/PlotsRdbpTbin/%s/R_DB_ratio_vs_masscut_ptbin%d.pdf", wp_label, bin));

        delete mass_ratio_canvas;
        delete g_data_mc_ratio_mass;
        delete g_data_combined_ratio_mass;
        delete g_bkg_mc_ratio_mass;
        delete ratio_line;
        delete mass_ratio_legend;
    }
    
    // Create a summary plot comparing all mass cuts in one figure (for data/combined MC ratio)
    TCanvas* all_ratio_canvas = new TCanvas("all_ratio_canvas", "Data/Combined MC Ratio for All pT and Mass Cuts", 1000, 700);
    all_ratio_canvas->cd();
    
    // We'll create a 2D plot to visualize all the ratios
    TH2F* ratio_map = new TH2F("ratio_map", Form("Data/Combined MC Ratio (WP: %s);Mass Cut;p_{T} Bin", wp_label), nMassCuts, -0.5, nMassCuts-0.5, nPtBins, -0.5, nPtBins-0.5);
    
    // Fill the 2D plot
    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < nPtBins; bin++) {
            double ratio = data_means[m][bin] / combined_means[m][bin];
            ratio_map->SetBinContent(m+1, bin+1, ratio);
            
            // Calculate error
            double ratio_err = (combined_means[m][bin] != 0) ? ratio * std::sqrt(
                std::pow(data_mean_errs[m][bin]/data_means[m][bin], 2) + 
                std::pow(combined_mean_errs[m][bin]/combined_means[m][bin], 2)) : 0.0;
            ratio_map->SetBinError(m+1, bin+1, ratio_err);
        }
    }
    
    // Style the plot
    ratio_map->SetStats(0);
    ratio_map->SetMarkerSize(1.5);
    
    // Set custom bin labels
    for (int m = 0; m < nMassCuts; m++) {
        ratio_map->GetXaxis()->SetBinLabel(m+1, massCutLabels[m]);
    }
    
    for (int bin = 0; bin < nPtBins; bin++) {
        // Simplify pT bin labels for readability
        TString simplifiedLabel = Form("%d-%d GeV", (int)ptBins[bin], (int)ptBins[bin+1]);
        ratio_map->GetYaxis()->SetBinLabel(bin+1, simplifiedLabel.Data());
    }
    
    ratio_map->GetXaxis()->SetLabelSize(0.04);
    ratio_map->GetYaxis()->SetLabelSize(0.04);
    
    gStyle->SetNumberContours(100);
   
    // Set a nice range centered at 1
    ratio_map->SetMinimum(0.95);
    ratio_map->SetMaximum(1.05);
    
    // Draw the plot
    ratio_map->Draw("COLZ TEXT");
    
    // Add ATLAS label
    DrawAtlasLabel();
    
    all_ratio_canvas->SaveAs(Form("../Plots/DB/PlotsRdbpTbin/%s/R_DB_all_ratio_summary.pdf", wp_label));
    
    delete all_ratio_canvas;
    delete ratio_map;

    // Let's also create a summary heatmap for multijet/MC ratio
    TCanvas* bkg_ratio_canvas = new TCanvas("bkg_ratio_canvas", "Multijets/Signal MC Ratio for All pT and Mass Cuts", 1000, 700);
    bkg_ratio_canvas->cd();
    
    // Create a 2D plot for background to MC ratio
    TH2F* bkg_ratio_map = new TH2F("bkg_ratio_map", Form("Multijets/Signal MC Ratio (WP: %s);Mass Cut;p_{T} Bin", wp_label), nMassCuts, -0.5, nMassCuts-0.5, nPtBins, -0.5, nPtBins-0.5);
    
    // Fill the 2D plot
    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < nPtBins; bin++) {
            double ratio = (mc_means[m][bin] != 0) ? bkg_means[m][bin] / mc_means[m][bin] : 0.0;
            bkg_ratio_map->SetBinContent(m+1, bin+1, ratio);
            
            // Calculate error
            double ratio_err = (mc_means[m][bin] != 0) ? ratio * std::sqrt(
                std::pow(bkg_mean_errs[m][bin]/bkg_means[m][bin], 2) + 
                std::pow(mc_mean_errs[m][bin]/mc_means[m][bin], 2)) : 0.0;
            bkg_ratio_map->SetBinError(m+1, bin+1, ratio_err);
        }
    }
    
    // Style the plot
    bkg_ratio_map->SetStats(0);
    bkg_ratio_map->SetMarkerSize(1.5);
    
    // Set custom bin labels (same as previous)
    for (int m = 0; m < nMassCuts; m++) {
        bkg_ratio_map->GetXaxis()->SetBinLabel(m+1, massCutLabels[m]);
    }
    
    for (int bin = 0; bin < nPtBins; bin++) {
        TString simplifiedLabel = Form("%d-%d GeV", (int)ptBins[bin], (int)ptBins[bin+1]);
        bkg_ratio_map->GetYaxis()->SetBinLabel(bin+1, simplifiedLabel.Data());
    }
    
    bkg_ratio_map->GetXaxis()->SetLabelSize(0.04);
    bkg_ratio_map->GetYaxis()->SetLabelSize(0.04);
    
    // Set the range - might need to adjust based on your values
    bkg_ratio_map->SetMinimum(0.95);
    bkg_ratio_map->SetMaximum(1.05);
    
    // Draw the plot
    bkg_ratio_map->Draw("COLZ TEXT");
    
    // Add ATLAS label
    DrawAtlasLabel();
    
    bkg_ratio_canvas->SaveAs(Form("../Plots/DB/PlotsRdbpTbin/%s/R_DB_bkg_ratio_summary.pdf", wp_label));
    
    delete bkg_ratio_canvas;
    delete bkg_ratio_map;
    
    // Clean up all dynamically allocated histograms for this working point
    for (int m = 0; m < nMassCuts; m++) {
        for (int i = 0; i < nPtBins; i++) {
            delete h_mc_R_DB_ptbin[m][i];
            delete h_data_R_DB_ptbin[m][i];
            delete h_multijets_R_DB_ptbin[m][i];
            delete h_combined_MC_R_DB_ptbin[m][i];
        }
    }

    // Close input files
    fdata->Close();
    fmc->Close();
    fmultijets->Close();
    
    std::cout << "All plots and summary for WP " << wp_label << " generated successfully!" << std::endl;
}

// Main function to run analysis for all working points
void Rdb_pT_bins_MCut() {
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    gROOT->LoadMacro("AtlasLabels.C");

    const char* workingPoints[] = {"0p25", "0p46", "0p74"};
    for (const char* wp : workingPoints) {
        ProcessWorkingPoint(wp);
    }
    std::cout << "\n--- All working points processed ---" << std::endl;
}