#include "../../Util/AtlasStyle.C" // Include ATLAS style macro
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip> // For std::fixed and std::setprecision
#include <fstream> // For ofstream
#include <iomanip> // For std::setw

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

/**
 * Draws the standard ATLAS label on plots
 */
void DrawAtlasLabel() {
    TLatex latex;
    latex.SetNDC(); // Set NDC (Normalized Device Coordinates) for positioning
    latex.SetTextFont(72); // ATLAS font
    latex.SetTextColor(kBlack);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.2, 0.88, "ATLAS"); // Draw "ATLAS" label
    latex.SetTextFont(42); // Regular font
    latex.DrawLatex(0.3, 0.88, "Internal"); // Draw "Internal" label
    latex.DrawLatex(0.2, 0.84, "#sqrt{s} = 13 TeV, 140 fb^{-1}"); // Draw luminosity and energy
}

/**
 * Performs mass isolation and signal extraction analysis,
 * iterating through specified pT bins for a given working point.
 * @param wp The working point string (e.g., "25", "46", "74").
 */
void RunMassIsolationPtBinsAnalysis(const std::string& wp) {
    // Determine file-friendly and display-friendly working point labels
    std::string fileWP; // Used for file paths (e.g., "0p46", "1p25")
    std::string labelWP_display; // Used for plot titles and printouts (e.g., "0.46 % QCD Eff. WP")

    if (wp == "125") {
        fileWP = "1p25";
        labelWP_display = "1.25 % QCD Eff. WP";
    } else if (wp == "155") {
        fileWP = "1p55";
        labelWP_display = "1.55 % QCD Eff. WP";
    } else {
        fileWP = "0p" + wp;
        double wpValue = std::stod(wp) / 100.0;
        std::stringstream ss;
        ss << std::fixed << std::setprecision(2) << wpValue << " % QCD Eff. WP";
        labelWP_display = ss.str();
    }

    std::cout << "Processing Mass Isolation for Working Point: " << labelWP_display << std::endl;

    // Construct base output directory for the current working point
    std::string baseOutputDirPath = "../Plots/DB/MasspTbins/" + fileWP + "/";
    gSystem->Exec(("mkdir -p " + baseOutputDirPath).c_str());
    
    // Construct input file paths dynamically
    std::string mcPath = "../../NtupleSlim/DB_" + fileWP + "wp_slim/Zbby_" + fileWP + "wp_slim.root";
    std::string dataPath = "../../NtupleSlim/DB_" + fileWP + "wp_slim/data_" + fileWP + "wp_slim.root";
    std::string multijetsPath = "../../NtupleSlim/DB_" + fileWP + "wp_slim/gamma_jets_" + fileWP + "wp_slim.root";

    // Open input files
    TFile *fmc = TFile::Open(mcPath.c_str(), "READ");
    TFile *fdata = TFile::Open(dataPath.c_str(), "READ");
    TFile *fmultijets = TFile::Open(multijetsPath.c_str(), "READ");

    // Check if files were opened successfully
    if (!fmc || fmc->IsZombie()) {
        std::cerr << "Error: Could not open MC file " << mcPath << " for WP " << labelWP_display << "." << std::endl;
        if (fmc) delete fmc;
        if (fdata) delete fdata;
        if (fmultijets) delete fmultijets;
        return; // Exit function if files cannot be opened
    }
    if (!fdata || fdata->IsZombie()) {
        std::cerr << "Error: Could not open data file " << dataPath << " for WP " << labelWP_display << "." << std::endl;
        if (fmc) delete fmc;
        if (fdata) delete fdata;
        if (fmultijets) delete fmultijets;
        return;
    }
    if (!fmultijets || fmultijets->IsZombie()) {
        std::cerr << "Error: Could not open multijets file " << multijetsPath << " for WP " << labelWP_display << "." << std::endl;
        if (fmc) delete fmc;
        if (fdata) delete fdata;
        if (fmultijets) delete fmultijets;
        return;
    }

    TTree *tree_mc = (TTree*)fmc->Get("nominal");
    TTree *tree_data = (TTree*)fdata->Get("nominal");
    TTree *tree_multijets = (TTree*)fmultijets->Get("nominal");

    // Explicitly check if TTree pointers are valid
    if (!tree_mc) {
        std::cerr << "Error: Could not retrieve 'nominal' tree from MC file for WP " << labelWP_display << "." << std::endl;
        fmc->Close(); fdata->Close(); fmultijets->Close();
        delete fmc; delete fdata; delete fmultijets;
        return;
    }
    if (!tree_data) {
        std::cerr << "Error: Could not retrieve 'nominal' tree from data file for WP " << labelWP_display << "." << std::endl;
        fmc->Close(); fdata->Close(); fmultijets->Close();
        delete fmc; delete fdata; delete fmultijets;
        return;
    }
    if (!tree_multijets) {
        std::cerr << "Error: Could not retrieve 'nominal' tree from multijets file for WP " << labelWP_display << "." << std::endl;
        fmc->Close(); fdata->Close(); fmultijets->Close();
        delete fmc; delete fdata; delete fmultijets;
        return;
    }

    // Define pT bins for analysis
    const int nPtBins = 3;
    double ptBins[nPtBins+1] = {450, 600, 900, 1200}; // 450-600, 600-900, 900-1200
    TString ptBinLabels[nPtBins] = {"450_600", "600_900", "900_1200"};
    
    // Create output directories for each pT bin within the working point directory
    for (int i = 0; i < nPtBins; i++) {
        gSystem->Exec(Form("mkdir -p %spt%s", baseOutputDirPath.c_str(), ptBinLabels[i].Data()));
    }
    
    // Create a summary file for all pT bins for the current working point
    ofstream summaryAll((baseOutputDirPath + "all_pt_bins_summary.txt").c_str());
    if (summaryAll.is_open()) {
        summaryAll << "Large-R Jet Mass Analysis Summary - All pT Bins (WP: " << labelWP_display << ")" << endl;
        summaryAll << "=============================================================" << endl << endl;
        summaryAll << left << setw(15) << "pT Bin [GeV]" 
                   << setw(20) << "MC Mean [GeV]" 
                   << setw(20) << "Data Mean [GeV]" 
                   << setw(20) << "Delta Mean [GeV]" 
                   << setw(20) << "MC Sigma [GeV]" 
                   << setw(20) << "Data Sigma [GeV]"
                   << setw(20) << "Signal Events"
                   << endl;
        summaryAll << string(115, '-') << endl;
    }

    // Arrays to store results from all pT bins for comparison
    double mc_means[nPtBins], mc_mean_errs[nPtBins];
    double data_means[nPtBins], data_mean_errs[nPtBins];
    double mc_sigmas[nPtBins], mc_sigma_errs[nPtBins];
    double data_sigmas[nPtBins], data_sigma_errs[nPtBins];
    double mean_diffs[nPtBins], mean_diff_errs[nPtBins];
    
    // Set up variables for reading from trees
    std::vector<float>* ljet_m = nullptr;
    std::vector<float>* ljet_pt = nullptr;
    double weight, weight_multijets;

    // Set branch addresses for all trees
    tree_mc->SetBranchAddress("ljet_m", &ljet_m);
    tree_data->SetBranchAddress("ljet_m", &ljet_m);
    tree_multijets->SetBranchAddress("ljet_m", &ljet_m);
    
    tree_mc->SetBranchAddress("ljet_pt", &ljet_pt);
    tree_data->SetBranchAddress("ljet_pt", &ljet_pt);
    tree_multijets->SetBranchAddress("ljet_pt", &ljet_pt);
    
    tree_mc->SetBranchAddress("total_weight", &weight);
    tree_multijets->SetBranchAddress("total_weight", &weight_multijets);

    // Loop over pT bins
    for (int ptBin = 0; ptBin < nPtBins; ptBin++) {
        
        // Create a temporary directory for this pT bin's objects within the current WP
        // This isolates objects created in this iteration, preventing name collisions
        // and simplifying cleanup.
        TDirectory* currentDir = gDirectory; // Store current directory
        gDirectory->mkdir(Form("WP_%s_PtBin_%s", fileWP.c_str(), ptBinLabels[ptBin].Data()))->cd(); // Create and cd into new directory

        cout << "\n=============================================================" << endl;
        cout << "Processing pT bin: " << ptBins[ptBin] << "-" << ptBins[ptBin+1] << " GeV (WP: " << labelWP_display << ")" << endl;
        cout << "=============================================================" << endl;
        
        // Create histograms for large-R jet mass in this pT bin
        TH1F *h_mc_ljet_m = new TH1F(Form("h_mc_ljet_m_pt%s", ptBinLabels[ptBin].Data()), 
                                     Form("Large-R Jet Mass (p_{T}: %d-%d GeV)", (int)ptBins[ptBin], (int)ptBins[ptBin+1]), 
                                     35, 50, 150);
        TH1F *h_data_ljet_m = new TH1F(Form("h_data_ljet_m_pt%s", ptBinLabels[ptBin].Data()), 
                                       Form("Large-R Jet Mass (p_{T}: %d-%d GeV)", (int)ptBins[ptBin], (int)ptBins[ptBin+1]), 
                                       35, 50, 150);
        TH1F *h_multijets_ljet_m = new TH1F(Form("h_multijets_ljet_m_pt%s", ptBinLabels[ptBin].Data()), 
                                           Form("Large-R Jet Mass (p_{T}: %d-%d GeV)", (int)ptBins[ptBin], (int)ptBins[ptBin+1]), 
                                           35, 50, 150);

        // Enable Sumw2 for histograms to ensure proper error calculation
        //h_mc_ljet_m->Sumw2();
        //h_data_ljet_m->Sumw2();
        //h_multijets_ljet_m->Sumw2();

        // Set histogram styling
        h_data_ljet_m->SetMarkerStyle(20);
        h_data_ljet_m->SetMarkerSize(0.8);
        h_data_ljet_m->SetMarkerColor(kBlack);
        h_data_ljet_m->SetLineColor(kBlack);

        h_multijets_ljet_m->SetFillColor(kBlue);
        h_multijets_ljet_m->SetLineColor(kBlue);

        h_mc_ljet_m->SetFillColor(kRed);
        h_mc_ljet_m->SetLineColor(kRed);

        // Fill histograms for MC signal
        Long64_t nentries_mc = tree_mc->GetEntries();
        for (Long64_t i = 0; i < nentries_mc; i++) {
            tree_mc->GetEntry(i);
            
            // Apply pT cut
            double jetPt = (*ljet_pt)[0]/1000.0; // Convert MeV to GeV
            if (jetPt >= ptBins[ptBin] && jetPt < ptBins[ptBin+1]) {
                h_mc_ljet_m->Fill((*ljet_m)[0]/1000.0, weight); // Convert MeV to GeV
            }
        }

        // Fill histograms for data
        Long64_t nentries_data = tree_data->GetEntries();
        for (Long64_t i = 0; i < nentries_data; i++) {
            tree_data->GetEntry(i);
            
            // Apply pT cut
            double jetPt = (*ljet_pt)[0]/1000.0; // Convert MeV to GeV
            if (jetPt >= ptBins[ptBin] && jetPt < ptBins[ptBin+1]) {
                h_data_ljet_m->Fill((*ljet_m)[0]/1000.0); // Convert MeV to GeV
            }
        }

        // Fill histograms for multijets background
        Long64_t nentries_multijets = tree_multijets->GetEntries();
        for (Long64_t i = 0; i < nentries_multijets; i++) {
            tree_multijets->GetEntry(i);
            
            // Apply pT cut
            double jetPt = (*ljet_pt)[0]/1000.0; // Convert MeV to GeV
            if (jetPt >= ptBins[ptBin] && jetPt < ptBins[ptBin+1]) {
                h_multijets_ljet_m->Fill((*ljet_m)[0]/1000.0, weight_multijets); // Convert MeV to GeV
            }
        }

        // Define sideband regions
        int binLow1 = h_multijets_ljet_m->FindBin(50);
        int binHigh1 = h_multijets_ljet_m->FindBin(75); // Changed from 65 to 75 for consistency with previous script
        int binLow2 = h_multijets_ljet_m->FindBin(110);
        int binHigh2 = h_multijets_ljet_m->FindBin(150);

        // Create arrays to store x and y values for the fit
        std::vector<double> x_values;
        std::vector<double> y_values;
        std::vector<double> y_errors;
        
        // Extract bin contents for first sideband region
        for (int bin = binLow1; bin <= binHigh1; ++bin) {
            double bin_center = h_multijets_ljet_m->GetBinCenter(bin);
            double data_content = h_data_ljet_m->GetBinContent(bin);
            double mc_content = h_multijets_ljet_m->GetBinContent(bin);
            
            // Avoid division by zero
            if (mc_content > 0) {
                x_values.push_back(bin_center);
                y_values.push_back(data_content / mc_content);
                
                // Calculate error on the ratio
                double data_error = h_data_ljet_m->GetBinError(bin);
                double mc_error = h_multijets_ljet_m->GetBinError(bin);
                double ratio_error = sqrt(pow(data_error/mc_content, 2) + 
                                         pow(data_content*mc_error/(mc_content*mc_content), 2));
                y_errors.push_back(ratio_error);
            }
        }
        
        // Extract bin contents for second sideband region
        for (int bin = binLow2; bin <= binHigh2; ++bin) {
            double bin_center = h_multijets_ljet_m->GetBinCenter(bin);
            double data_content = h_data_ljet_m->GetBinContent(bin);
            double mc_content = h_multijets_ljet_m->GetBinContent(bin);
            
            // Avoid division by zero
            if (mc_content > 0) {
                x_values.push_back(bin_center);
                y_values.push_back(data_content / mc_content);
                
                // Calculate error on the ratio
                double data_error = h_data_ljet_m->GetBinError(bin);
                double mc_error = h_multijets_ljet_m->GetBinError(bin);
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
        TF1* fit_func = new TF1("fit_func", "[0] + [1] * x", 50, 150);
        fit_func->SetParameters(1.0, 0.0);  // Initial guess
        
        // Perform the fit
        ratio_graph->Fit("fit_func", "R");  // "R" restricts fit to the specified range
        
        // Get fit parameters
        double p0 = fit_func->GetParameter(0);  // Intercept
        double p1 = fit_func->GetParameter(1);  // Slope
        double p0_err = fit_func->GetParError(0);
        double p1_err = fit_func->GetParError(1);
        
        printf("Linear fit results: SF(m) = %.4f (± %.4f) + %.8f (± %.8f) * m\n", 
               p0, p0_err, p1, p1_err);
        
        // Save the fit result to a canvas
        TCanvas* fit_canvas = new TCanvas("fit_canvas", 
                                        Form("Scale Factor Fit (p_{T}: %d-%d GeV, WP: %s)", (int)ptBins[ptBin], (int)ptBins[ptBin+1], labelWP_display.c_str()), 
                                        800, 600);
        ratio_graph->SetTitle(Form("Data/Multijets Ratio in Sidebands (p_{T}: %d-%d GeV, WP: %s);Large-R Jet Mass [GeV];Data/MC", 
                                  (int)ptBins[ptBin], (int)ptBins[ptBin+1], labelWP_display.c_str()));
        ratio_graph->SetMarkerStyle(20);
        ratio_graph->SetMarkerColor(kBlue);
        ratio_graph->SetLineColor(kBlue);
        ratio_graph->Draw("AP");
        fit_func->SetLineColor(kRed);
        fit_func->Draw("same");
        
        // Add fit parameters to the canvas
        TLatex latex_fit_pt; // Renamed to avoid conflict with global latex
        latex_fit_pt.SetNDC();
        latex_fit_pt.SetTextFont(42);
        latex_fit_pt.SetTextSize(0.035);
        latex_fit_pt.DrawLatex(0.2, 0.85, Form("SF(m) = %.4f + %.8f * m", p0, p1));
        latex_fit_pt.DrawLatex(0.2, 0.80, Form("#chi^{2}/ndf = %.2f/%d", fit_func->GetChisquare(), fit_func->GetNDF()));
        
        DrawAtlasLabel();
        fit_canvas->SaveAs(Form("%spt%s/scale_factor_fit.pdf", baseOutputDirPath.c_str(), ptBinLabels[ptBin].Data()));
        
        // Apply scale factor to the multijets background histogram
        TH1F* h_multijets_ljet_m_scaled = (TH1F*)h_multijets_ljet_m->Clone("h_multijets_ljet_m_scaled");
        //h_multijets_ljet_m_scaled->Sumw2(); // Ensure Sumw2 is enabled for scaled histogram

        // Apply mass-dependent scale factor bin by bin
        for (int bin = 1; bin <= h_multijets_ljet_m_scaled->GetNbinsX(); ++bin) {
            double bin_center = h_multijets_ljet_m_scaled->GetBinCenter(bin);
            double content = h_multijets_ljet_m_scaled->GetBinContent(bin);
            double error = h_multijets_ljet_m_scaled->GetBinError(bin);
            
            // Calculate bin-by-bin scale factor
            double scale_factor = p0 + p1 * bin_center;
            
            // Apply the scale factor
            h_multijets_ljet_m_scaled->SetBinContent(bin, content * scale_factor);
            h_multijets_ljet_m_scaled->SetBinError(bin, error * scale_factor);
        }

        // Estimate background in the signal region using scaled histogram
        double estimatedBackground = h_multijets_ljet_m_scaled->Integral(binHigh1+1, binLow2-1);
        
        // Calculate signal events
        double dataSignal = h_data_ljet_m->Integral(binHigh1+1, binLow2-1);
        double signalEvents = dataSignal - estimatedBackground;
        
        // Print signal events
        std::cout << "Number of signal events: " << signalEvents << std::endl;
        
        // Create a polynomial-exponential fit for background modeling
        TF1 *polyExpFit = new TF1("polyExpFit", "[0]*exp([1]*x + [2]*x*x + [3]*x*x*x)", 50, 150);
        
        // Clone the data histogram for fitting the sidebands
        TH1F *h_data_sidebands = (TH1F*)h_data_ljet_m->Clone("h_data_sidebands");
        //h_data_sidebands->Sumw2(); // Ensure Sumw2 is enabled for sidebands clone
        
        // Remove the bins in the signal region from the clone
        for (int i = binHigh1 + 1; i < binLow2; ++i) {
            h_data_sidebands->SetBinContent(i, 0);
            h_data_sidebands->SetBinError(i, 0);
        }
        
        // Set initial parameters for polynomial exponential fit
        polyExpFit->SetParameter(0, h_data_sidebands->GetMaximum());    // amplitude
        polyExpFit->SetParameter(1, -0.05);  // first order term (negative for falling exp)
        polyExpFit->SetParameter(2, -0.0005); // second order term
        polyExpFit->SetParameter(3, 0.000001); // third order term
        
        // Fit the sideband regions of the data
        h_data_sidebands->Fit(polyExpFit, "R");
        
        // Get the chi-square value and number of degrees of freedom
        double chiSquare = polyExpFit->GetChisquare();
        int ndf = polyExpFit->GetNDF();
        double chi2NDF = (ndf > 0) ? chiSquare / ndf : 0.0; // Avoid division by zero
        
        // Create a clone of the data histogram to store the extracted signal
        TH1F *h_signal_ljet_m = (TH1F*)h_data_ljet_m->Clone("h_signal_ljet_m");
        h_signal_ljet_m->SetTitle("Isolated Signal (Data);M^{J} [GeV];Events"); // Set title for drawing
        //h_signal_ljet_m->Sumw2(); // Ensure Sumw2 is enabled for signal histogram

        // Create a histogram representing the background fit for subtraction
        TH1F *h_background_fit = (TH1F*)h_data_ljet_m->Clone("h_background_fit_for_subtraction");
        h_background_fit->Reset(); // Clear contents
        //h_background_fit->Sumw2(); // Enable Sumw2 for proper error handling during subtraction

        // Fill the background fit histogram with values from polyExpFit
        for (int i = 1; i <= h_background_fit->GetNbinsX(); ++i) {
            double binCenter = h_background_fit->GetBinCenter(i);
            h_background_fit->SetBinContent(i, polyExpFit->Eval(binCenter));
            // For statistical error propagation, we assume the statistical uncertainty of the fit
            // itself is negligible compared to data statistical uncertainty, or handled as systematic.
            // So, we set its statistical error to 0 for this operation.
            h_background_fit->SetBinError(i, 0.0); 
        }

        // Subtract the background fit histogram from the data histogram to isolate the signal.
        // ROOT's TH1::Add with a negative scale factor will correctly propagate errors
        // if Sumw2 is enabled on both histograms.
        h_signal_ljet_m->Add(h_background_fit, -1.0);

        // Set negative values to zero after subtraction (and adjust errors accordingly)
        for (int i = 1; i <= h_signal_ljet_m->GetNbinsX(); ++i) {
            if (h_signal_ljet_m->GetBinContent(i) < 0) {
                h_signal_ljet_m->SetBinContent(i, 0);
                h_signal_ljet_m->SetBinError(i, 0); // Error also becomes 0 if content is 0
            }
        }
        
        // Create stack plot with data overlay
        TCanvas* stack_canvas = new TCanvas("stack_canvas", 
                                          Form("Mass Distribution with Scale Factor (p_{T}: %d-%d GeV, WP: %s)", (int)ptBins[ptBin], (int)ptBins[ptBin+1], labelWP_display.c_str()), 
                                          800, 800);
        stack_canvas->Divide(1, 2);
        
        TPad* pad0 = new TPad("pad0", "pad0", 0, 0.3, 1, 1);
        pad0->SetBottomMargin(0.015);
        pad0->Draw();
        pad0->cd();
        
        THStack* hs = new THStack("hs", Form("Large-R Jet Mass (p_{T}: %d-%d GeV, WP: %s);M^{J} [GeV];Events", 
                                           (int)ptBins[ptBin], (int)ptBins[ptBin+1], labelWP_display.c_str()));
        hs->Add(h_multijets_ljet_m_scaled);
        hs->Add(h_mc_ljet_m);
        
        hs->SetMaximum(hs->GetMaximum()*1.25);
        hs->Draw("HIST");
        
        // Draw the polynomial exponential fit
        polyExpFit->SetLineColor(kGreen);
        polyExpFit->Draw("same");
        
        h_data_ljet_m->Draw("EP same");
        
        hs->GetXaxis()->SetLabelSize(0);
        h_data_ljet_m->GetXaxis()->SetLabelSize(0);
        
        // MC statistical uncertainty band for stacked histogram
        TH1F* h_mc_total = (TH1F*)hs->GetStack()->Last();
        TGraphAsymmErrors* err_band = new TGraphAsymmErrors(h_mc_total);
        for (int i = 0; i < err_band->GetN(); ++i) {
            double content = h_mc_total->GetBinContent(i + 1);
            double error = h_mc_total->GetBinError(i + 1);
            err_band->SetPointEYhigh(i, error);
            err_band->SetPointEYlow(i, error);
        }
        err_band->SetFillColorAlpha(kBlack, 0.75);
        err_band->SetFillStyle(3254);
        err_band->SetLineWidth(0);
        err_band->Draw("2 same");
        
        TLegend* legend = new TLegend(0.6, 0.6, 0.9, 0.9);
        legend->SetBorderSize(0);
        legend->SetTextSize(0.03);
        legend->AddEntry(h_data_ljet_m, "Data", "lep");
        legend->AddEntry(h_mc_ljet_m, "Z (#rightarrow bb) + #gamma", "f");
        legend->AddEntry(h_multijets_ljet_m_scaled, "#gamma + jets", "f");
        legend->AddEntry(err_band, "MC Stat. Unc.", "f");
        legend->AddEntry(polyExpFit, Form("Poly-Exp Fit (#chi^{2}/NDF = %.2f)", chi2NDF), "l");
        legend->Draw();
        
        DrawAtlasLabel();
        
        // Setup the ratio pad
        stack_canvas->cd();
        TPad* pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
        pad2->SetTopMargin(0.05);
        pad2->SetBottomMargin(0.35);
        pad2->Draw();
        pad2->cd();
        
        TH1F* h_ratio = (TH1F*)h_data_ljet_m->Clone("h_ratio");
        h_ratio->SetLineColor(kBlack);
        h_ratio->SetMinimum(0.5);
        h_ratio->SetMaximum(1.5);
        h_ratio->Divide(h_mc_total);
        h_ratio->SetMarkerStyle(21);
        h_ratio->GetXaxis()->SetTitle("M^{J} [GeV]");
        h_ratio->GetXaxis()->SetTitleSize(0.12);
        h_ratio->GetXaxis()->SetLabelSize(0.12);
        h_ratio->GetXaxis()->SetTitleOffset(1.25);
        h_ratio->GetYaxis()->SetTitle("Data/MC");
        h_ratio->GetYaxis()->SetTitleSize(0.15);
        h_ratio->GetYaxis()->SetLabelSize(0.12);
        h_ratio->GetYaxis()->SetTitleOffset(0.3);
        h_ratio->GetYaxis()->SetNdivisions(505);
        h_ratio->Draw("EP");
        
        // Draw ratio uncertainty band
        TGraphAsymmErrors* ratio_err_band = new TGraphAsymmErrors(h_mc_total);
        for (int i = 0; i < ratio_err_band->GetN(); ++i) {
            ratio_err_band->SetPointY(i, 1);
            double content = h_mc_total->GetBinContent(i + 1);
            double error = h_mc_total->GetBinError(i + 1);
            double rel_error = (content != 0) ? error / content : 0;
            ratio_err_band->SetPointEYhigh(i, rel_error);
            ratio_err_band->SetPointEYlow(i, rel_error);
        }
        ratio_err_band->SetFillColorAlpha(kBlack, 0.75);
        ratio_err_band->SetFillStyle(3254);
        ratio_err_band->Draw("2 same");
        
        // Draw reference line at 1
        TLine* line = new TLine(pad2->GetUxmin(), 1.0, pad2->GetUxmax(), 1.0);
        line->SetLineColor(kBlack);
        line->SetLineWidth(2);
        line->SetLineStyle(2);
        line->Draw("same");
        
        stack_canvas->SaveAs(Form("%spt%s/Mass_stack_scaled.pdf", baseOutputDirPath.c_str(), ptBinLabels[ptBin].Data()));
        
        // Visualize the sideband regions
        TCanvas* sideband_canvas = new TCanvas("sideband_canvas", 
                                             Form("Fit Visualization (p_{T}: %d-%d GeV, WP: %s)", (int)ptBins[ptBin], (int)ptBins[ptBin+1], labelWP_display.c_str()), 
                                             800, 600);
        sideband_canvas->cd();
        
        h_data_ljet_m->Draw("EP");
        polyExpFit->Draw("same");
        
        // Highlight sideband regions
        double xlow1 = h_data_ljet_m->GetBinLowEdge(binLow1);
        double xhigh1 = h_data_ljet_m->GetBinLowEdge(binHigh1+1);
        double xlow2 = h_data_ljet_m->GetBinLowEdge(binLow2);
        double xhigh2 = h_data_ljet_m->GetBinLowEdge(binHigh2+1);
        
        TBox *sidebandBox1 = new TBox(xlow1, 0, xhigh1, h_data_ljet_m->GetMaximum()*1.75);
        sidebandBox1->SetFillColorAlpha(kYellow, 0.3);
        sidebandBox1->Draw("same");
        
        TBox *sidebandBox2 = new TBox(xlow2, 0, xhigh2, h_data_ljet_m->GetMaximum()*1.75);
        sidebandBox2->SetFillColorAlpha(kYellow, 0.3);
        sidebandBox2->Draw("same");
        
        TLegend *legendSideband = new TLegend(0.6, 0.65, 0.9, 0.85);
        legendSideband->SetBorderSize(0);
        legendSideband->SetTextSize(0.03);
        legendSideband->AddEntry(h_data_ljet_m, "Data", "lep");
        legendSideband->AddEntry(polyExpFit, Form("Poly-Exp Fit (#chi^{2}/NDF = %.2f)", chi2NDF), "l");
        legendSideband->AddEntry(sidebandBox1, "Sideband Regions", "f");
        legendSideband->Draw();
        
        DrawAtlasLabel();
        sideband_canvas->SaveAs(Form("%spt%s/SideBandFitDebug.pdf", baseOutputDirPath.c_str(), ptBinLabels[ptBin].Data()));
        
        // Plot the isolated signal with MC signal
        TCanvas* signal_canvas = new TCanvas("signal_canvas", 
                                           Form("Isolated Signal with Gaussian Fits (p_{T}: %d-%d GeV, WP: %s)", (int)ptBins[ptBin], (int)ptBins[ptBin+1], labelWP_display.c_str()), 
                                           800, 600);
        signal_canvas->cd();
        
        h_mc_ljet_m->SetFillColor(kRed);
        h_mc_ljet_m->GetYaxis()->SetRangeUser(0, h_mc_ljet_m->GetMaximum()*1.75);
        h_mc_ljet_m->GetXaxis()->SetTitle("M^{J} [GeV]"); // Added X-axis title
        h_mc_ljet_m->GetYaxis()->SetTitle("Events"); // Added Y-axis title
        h_mc_ljet_m->Draw("HIST");
        
        h_signal_ljet_m->SetMarkerStyle(20);
        h_signal_ljet_m->SetMarkerColor(kBlack);
        h_signal_ljet_m->Draw("EP same");
        
        // Fit the MC signal with a Gaussian in the Z mass region
        // Adjust fit range based on pT bin if needed
        double fit_min = 83; // Adjusted for better initial range
        double fit_max = 105; // Adjusted for better initial range
        
        // Use histogram properties for more robust initial guesses
        TF1 *mc_gauss_fit = new TF1("mc_gauss_fit", "gaus", fit_min, fit_max);
        mc_gauss_fit->SetParameter(0, h_mc_ljet_m->GetMaximum());  // Amplitude
        mc_gauss_fit->SetParameter(1, h_mc_ljet_m->GetXaxis()->GetBinCenter(h_mc_ljet_m->GetMaximumBin()));  // Mean
        mc_gauss_fit->SetParameter(2, h_mc_ljet_m->GetRMS() * 0.7); // Sigma (scaled RMS as a rough estimate)
        
        h_mc_ljet_m->Fit(mc_gauss_fit, "R");
        
        // Get the fit parameters
        double mc_amp = mc_gauss_fit->GetParameter(0);
        double mc_mean = mc_gauss_fit->GetParameter(1);
        double mc_sigma = mc_gauss_fit->GetParameter(2);
        double mc_mean_err = mc_gauss_fit->GetParError(1);
        double mc_sigma_err = mc_gauss_fit->GetParError(2);
        double mc_chi2 = mc_gauss_fit->GetChisquare();
        int mc_ndf = mc_gauss_fit->GetNDF();
        double mc_chi2_ndf = (mc_ndf > 0) ? mc_chi2 / mc_ndf : 0.0; // Avoid division by zero
        
        // Store for comparison across pT bins
        mc_means[ptBin] = mc_mean;
        mc_mean_errs[ptBin] = mc_mean_err;
        mc_sigmas[ptBin] = mc_sigma;
        mc_sigma_errs[ptBin] = mc_sigma_err;
        
        // Fit the extracted data signal with a Gaussian
        TF1 *data_gauss_fit = new TF1("data_gauss_fit", "gaus", fit_min, fit_max);
        
        // Use histogram properties for more robust initial guesses
        data_gauss_fit->SetParameter(0, h_signal_ljet_m->GetMaximum());  // Amplitude
        data_gauss_fit->SetParameter(1, h_signal_ljet_m->GetXaxis()->GetBinCenter(h_signal_ljet_m->GetMaximumBin()));  // Mean
        data_gauss_fit->SetParameter(2, h_signal_ljet_m->GetRMS() * 0.7); // Sigma (scaled RMS as a rough estimate)
        
        h_signal_ljet_m->Fit(data_gauss_fit, "R");
        
        // Get the fit parameters
        double data_amp = data_gauss_fit->GetParameter(0);
        double data_mean = data_gauss_fit->GetParameter(1);
        double data_sigma = data_gauss_fit->GetParameter(2);
        double data_mean_err = data_gauss_fit->GetParError(1);
        double data_sigma_err = data_gauss_fit->GetParError(2);
        double data_chi2 = data_gauss_fit->GetChisquare();
        int data_ndf = data_gauss_fit->GetNDF();
        double data_chi2_ndf = (data_ndf > 0) ? data_chi2 / data_ndf : 0.0; // Avoid division by zero
        
        // Store for comparison across pT bins
        data_means[ptBin] = data_mean;
        data_mean_errs[ptBin] = data_mean_err;
        data_sigmas[ptBin] = data_sigma;
        data_sigma_errs[ptBin] = data_sigma_err;
        
        // Calculate mass shift between data and MC
        double mean_diff = data_mean - mc_mean;
        double mean_diff_err = sqrt(data_mean_err*data_mean_err + mc_mean_err*mc_mean_err);
        
        // Store for comparison across pT bins
        mean_diffs[ptBin] = mean_diff;
        mean_diff_errs[ptBin] = mean_diff_err;
        
        // Style the fits
        mc_gauss_fit->SetLineColor(kCyan);
        mc_gauss_fit->SetLineWidth(2);
        
        data_gauss_fit->SetLineColor(kViolet-4);
        data_gauss_fit->SetLineWidth(2);
        
        mc_gauss_fit->Draw("same");
        data_gauss_fit->Draw("same");
        
        // Create a legend
        TLegend *legendSignal = new TLegend(0.4, 0.65, 0.9, 0.85); 
        legendSignal->SetBorderSize(0);
        legendSignal->SetTextSize(0.018); 
        legendSignal->AddEntry(h_signal_ljet_m, "Isolated Signal (Data)", "lep");
        legendSignal->AddEntry(h_mc_ljet_m, "MC Signal", "f");
        legendSignal->AddEntry(mc_gauss_fit, Form("MC Gauss Fit: #mu = %.2f #pm %.2f, #sigma = %.2f #pm %.2f, #chi^{2}/ndf = %.2f", 
                              mc_mean, mc_mean_err, mc_sigma, mc_sigma_err, mc_chi2_ndf), "l");
        legendSignal->AddEntry(data_gauss_fit, Form("Data Gauss Fit: #mu = %.2f #pm %.2f, #sigma = %.2f #pm %.2f, #chi^{2}/ndf = %.2f", 
                              data_mean, data_mean_err, data_sigma, data_sigma_err, data_chi2_ndf), "l");
        legendSignal->Draw();
        
        DrawAtlasLabel();
        signal_canvas->SaveAs(Form("%spt%s/Isolated_Signal_GaussFits.pdf", baseOutputDirPath.c_str(), ptBinLabels[ptBin].Data()));
        
        // Create a canvas for comparing the Gaussian fits
        TCanvas* gauss_canvas = new TCanvas("gauss_canvas", 
                                          Form("Gaussian Fits Comparison (p_{T}: %d-%d GeV, WP: %s)", (int)ptBins[ptBin], (int)ptBins[ptBin+1], labelWP_display.c_str()), 
                                          800, 600);
        gauss_canvas->cd();
        
        // Create a frame for drawing just the fits
        TH1F *h_frame = new TH1F("h_frame", 
                               Form("Gaussian Fits Comparison (p_{T}: %d-%d GeV);M^{J} [GeV];Arbitrary Units", 
                                    (int)ptBins[ptBin], (int)ptBins[ptBin+1]), 
                               100, 40, 150);
        h_frame->SetStats(0);
        h_frame->GetYaxis()->SetRangeUser(0, 1.2);
        h_frame->Draw();
        
        // Normalize the Gaussian fits to have the same height
        TF1 *norm_mc_fit = new TF1("norm_mc_fit", "gaus", 40, 150);
        norm_mc_fit->SetParameters(1.0, mc_mean, mc_sigma);
        norm_mc_fit->SetLineColor(kRed);
        norm_mc_fit->SetLineWidth(2);
        
        TF1 *norm_data_fit = new TF1("norm_data_fit", "gaus", 40, 150);
        norm_data_fit->SetParameters(1.0, data_mean, data_sigma);
        norm_data_fit->SetLineColor(kBlue);
        norm_data_fit->SetLineWidth(2);
        
        norm_mc_fit->Draw("same");
        norm_data_fit->Draw("same");
        
        // Add a legend
        TLegend *legendGauss = new TLegend(0.45, 0.6, 0.9, 0.85);
        legendGauss->SetBorderSize(0);
        legendGauss->SetTextSize(0.025); 
        legendGauss->AddEntry(norm_mc_fit, Form("MC Gauss: #mu = %.2f #pm %.2f, #sigma = %.2f #pm %.2f GeV", 
                              mc_mean, mc_mean_err, mc_sigma, mc_sigma_err), "l");
        legendGauss->AddEntry(norm_data_fit, Form("Data Gauss: #mu = %.2f #pm %.2f, #sigma = %.2f #pm %.2f GeV", 
                              data_mean, data_mean_err, data_sigma, data_sigma_err), "l");
        legendGauss->Draw();
        
        // Add a label showing the difference in means
        TLatex gauss_latex;
        gauss_latex.SetNDC();
        gauss_latex.SetTextFont(42);
        gauss_latex.SetTextSize(0.03);
        gauss_latex.DrawLatex(0.2, 0.8, Form("#Delta#mu = %.2f #pm %.2f GeV", mean_diff, mean_diff_err));
        gauss_latex.DrawLatex(0.2, 0.75, Form("p_{T}: %d-%d GeV", (int)ptBins[ptBin], (int)ptBins[ptBin+1]));
        
        DrawAtlasLabel();
        gauss_canvas->SaveAs(Form("%spt%s/GaussFitsComparison.pdf", baseOutputDirPath.c_str(), ptBinLabels[ptBin].Data()));
        
        // Create a summary text file for this pT bin
        ofstream summary(Form("%spt%s/fit_summary.txt", baseOutputDirPath.c_str(), ptBinLabels[ptBin].Data()));
        if (summary.is_open()) {
            summary << "Large-R Jet Mass Analysis Summary" << endl;
            summary << "=================================" << endl << endl;
            summary << "pT bin: " << ptBins[ptBin] << "-" << ptBins[ptBin+1] << " GeV" << endl << endl;
            
            summary << "Background Scale Factor:" << endl;
            summary << "----------------------" << endl;
            summary << "SF(m) = " << p0 << " #pm " << p0_err << " + (" << p1 << " #pm " << p1_err << ") * m" << endl;
            summary << "#chi^{2}/ndf = " << fit_func->GetChisquare() << "/" << fit_func->GetNDF() << " = " << fit_func->GetChisquare()/fit_func->GetNDF() << endl << endl;
            
            summary << "Sideband Fit to Data:" << endl;
            summary << "-------------------" << endl;
            summary << "Polynomial-Exponential Fit: [0]*exp([1]*x + [2]*x*x + [3]*x*x*x)" << endl;
            summary << "[0] = " << polyExpFit->GetParameter(0) << " #pm " << polyExpFit->GetParError(0) << endl;
            summary << "[1] = " << polyExpFit->GetParameter(1) << " #pm " << polyExpFit->GetParError(1) << endl;
            summary << "[2] = " << polyExpFit->GetParameter(2) << " #pm " << polyExpFit->GetParError(2) << endl;
            summary << "[3] = " << polyExpFit->GetParameter(3) << " #pm " << polyExpFit->GetParError(3) << endl;
            summary << "#chi^{2}/ndf = " << chiSquare << "/" << ndf << " = " << chi2NDF << endl << endl;
            
            summary << "Signal Extraction:" << endl;
            summary << "----------------" << endl;
            summary << "Signal Region: " << h_data_ljet_m->GetBinLowEdge(binHigh1+1) << " - " << h_data_ljet_m->GetBinLowEdge(binLow2) << " GeV" << endl;
            summary << "Data Events in Signal Region: " << dataSignal << endl;
            summary << "Estimated Background in Signal Region: " << estimatedBackground << endl;
            summary << "Extracted Signal Events: " << signalEvents << endl << endl;
            
            summary << "Gaussian Fits:" << endl;
            summary << "-------------" << endl;
            summary << "MC Signal: #mu = " << mc_mean << " #pm " << mc_mean_err << " GeV, #sigma = " << mc_sigma << " #pm " << mc_sigma_err << " GeV" << endl;
            summary << "Data Signal: #mu = " << data_mean << " #pm " << data_mean_err << " GeV, #sigma = " << data_sigma << " #pm " << data_sigma_err << " GeV" << endl;
            summary << "Mass Shift (Data-MC): #Delta#mu = " << mean_diff << " #pm " << mean_diff_err << " GeV" << endl;
            summary << "Width Ratio (Data/MC): #sigma_{data}/#sigma_{mc} = " << data_sigma/mc_sigma << " #pm " << 
                sqrt(pow(data_sigma_err/mc_sigma, 2) + pow(data_sigma*mc_sigma_err/(mc_sigma*mc_sigma), 2)) << endl;
            summary << "Number of signal events: " << signalEvents << endl;
            
            summary.close();
            cout << "Summary saved to " << baseOutputDirPath << "pt" << ptBinLabels[ptBin].Data() << "/fit_summary.txt" << endl;
        }

        // Return to the parent directory before deleting the temporary directory
        currentDir->cd();
        gDirectory->Delete(Form("WP_%s_PtBin_%s;*", fileWP.c_str(), ptBinLabels[ptBin].Data()));

        // Clean up histograms for the current pT bin to avoid memory leaks
        /*delete h_mc_ljet_m;
        delete h_data_ljet_m;
        delete h_multijets_ljet_m;
        delete ratio_graph;
        delete fit_func;
        delete h_multijets_ljet_m_scaled;
        delete polyExpFit;
        delete h_data_sidebands;
        delete h_signal_ljet_m;
        delete h_background_fit;
        delete stack_canvas;
        delete sideband_canvas;
        delete signal_canvas;
        delete gauss_canvas;
        delete hs;
        delete err_band;
        delete legend;
        delete pad0;
        delete pad2;
        delete h_ratio;
        delete ratio_err_band;
        delete line;
        delete sidebandBox1;
        delete sidebandBox2;
        delete legendSideband;
        delete legendSignal;
        delete legendGauss;
        delete h_frame;
        delete norm_mc_fit;
        delete norm_data_fit;*/
    } // End of pT bin loop

    // Close and delete TFile objects
    if (fmc) { fmc->Close(); delete fmc; }
    if (fdata) { fdata->Close(); delete fdata; }
    if (fmultijets) { fmultijets->Close(); delete fmultijets; }

    if (summaryAll.is_open()) {
        summaryAll.close();
    }
}

// Main function to run the analysis for different working points
void MassIsolationExtractSig_PtBins() {
    // Disable ROOT global graphics for batch processing
    gROOT->SetBatch(kTRUE);

    // Load ATLAS style settings for consistent plotting
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();

    // Define working points to process
    std::vector<std::string> workingPoints = {"25", "46", "74"}; // Example working points

    // Process each working point
    for (const auto& wp : workingPoints) {
        RunMassIsolationPtBinsAnalysis(wp);
    }
}
