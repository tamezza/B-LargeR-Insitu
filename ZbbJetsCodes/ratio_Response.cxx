#include "../../Util/AtlasStyle.C" // Include ATLAS style macro
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip> // For std::fixed, std::setprecision, std::setw
#include <fstream> // For ofstream
#include <algorithm> // For std::min/max

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
#include "TH2F.h" // For 2D histograms
#include "TPaletteAxis.h" // For palette customization
#include "TStyle.h" // For gStyle

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

/**
 * Draws the standard ATLAS label on plots
 */
void DrawAtlasLabel() {
    TLatex latex;
    latex.SetNDC(); // Set NDC (Normalized Device Coordinates) for positioning
    latex.SetTextFont(72); // ATLAS font
    latex.SetTextColor(kBlack);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.15, 0.925, "ATLAS"); // Draw "ATLAS" label
    latex.SetTextFont(42); // Regular font
    latex.DrawLatex(0.25, 0.925, "Internal"); // Draw "Internal" label
    latex.DrawLatex(0.15, 0.885, "#sqrt{s} = 13 TeV, 140 fb^{-1}"); // Draw luminosity and energy
}

// Global constants for pT and mass bins
const int nPtBins = 3;
const double ptBins[nPtBins + 1] = {450, 600, 900, 1200}; // 450-600, 600-900, 900-1200
const char* ptBinLabels[nPtBins] = {"450 < p_{T} < 600 GeV", "600 < p_{T} < 900 GeV", "900 < p_{T} < 1200 GeV"};

const int nMassCuts = 3; //nMassCuts = 6;
const char* massCutLabels[nMassCuts] = {"50-80 GeV", "80-110 GeV", "#geq110 GeV"}; //{"NoCut", "50-150", "50-200", "50-250", "50-300", "50-350"};
//const double massCutLow[nMassCuts] = {-1, 50, 50, 50, 50, 50};  // -1 means no lower cut
//const double massCutHigh[nMassCuts] = {-1, 150, 200, 250, 300, 350}; // -1 means no upper cut

const double massCutLow[nMassCuts] = {50, 80, 110};
const double massCutHigh[nMassCuts] = {80, 110, -1}; // -1 means no upper cut



/**
 * Creates necessary directories for output plots for a given working point.
 * @param fileWP The file-friendly working point string (e.g., "0p46").
 */
void CreateRatioDirectories(const std::string& fileWP) {
    gSystem->Exec(Form("mkdir -p ../Plots/MJB_Ratio/%s/PlotsRdbpTbin/Summary", fileWP.c_str()));
    //gSystem->Exec(Form("mkdir -p ../Plots/MJB/%s/PlotsRdbpTbin/Summary", fileWP.c_str())); // Commented out as per original code
    //gSystem->Exec(Form("mkdir -p ../Plots/MJB_Deco/%s/PlotsRdbpTbin/Summary", fileWP.c_str())); // Commented out as per original code
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

/**
 * Creates summary plots showing Post_bJR10V00/Pre_bJR10V00 ratios vs pT bins
 * with different mass bins as separate points at each pT bin
 * @param fileWP The file-friendly working point string (e.g., "0p46").
 * @param labelWP_display Display label for the current working point.
 * @param final_ratio_means Final ratio of (Data/Combined MC) means (Post_bJR10V00/Pre_bJR10V00).
 * @param final_ratio_mean_errs Errors on final_ratio_means.
 * @param final_ratio_sigmas Final ratio of (Data/Combined MC) resolutions (Post_bJR10V00/Pre_bJR10V00).
 * @param final_ratio_sigma_errs Errors on final_ratio_sigmas.
 */
void CreateSummaryPlots(
    const std::string& fileWP,
    const std::string& labelWP_display,
    double final_ratio_means[nMassCuts][nPtBins], double final_ratio_mean_errs[nMassCuts][nPtBins],
    double final_ratio_sigmas[nMassCuts][nPtBins], double final_ratio_sigma_errs[nMassCuts][nPtBins]
) {
    // Define colors for different mass bins
    Int_t colors[nMassCuts] = {kRed, kBlue, kGreen+2};
    Int_t markers[nMassCuts] = {20, 21, 22}; // Circle, square, triangle
    
    // Use sequential positions for x-axis (1, 2, 3 for the 3 pT bins)
    double ptCenters[nPtBins];
    for (int i = 0; i < nPtBins; i++) {
        ptCenters[i] = (ptBins[i] + ptBins[i+1]) / 2.0; //Calculate pT bin centers for plotting
        //ptCenters[i] = i + 1; // Position 1, 2, 3
    }

    
    // Create offset positions for mass bins at each pT point (in GeV)
    // These offsets will shift the data points for different mass bins horizontally within each pT bin
    double offsets[nMassCuts] = {-20.0, 0.0, 20.0}; // Small offsets in GeV to separate mass points

    
    // --- Summary Plot for Mean Ratios ---
    TCanvas* summary_mean_canvas = new TCanvas("summary_mean_canvas", 
        Form("Summary: Data/Combined MC Mean Ratios (Post_bJR10V00/Pre_bJR10V00) - WP: %s", labelWP_display.c_str()), 
        800, 600);
    summary_mean_canvas->SetLeftMargin(0.15);
    summary_mean_canvas->SetBottomMargin(0.15); // Increased for x-axis labels
    summary_mean_canvas->SetTopMargin(0.12);
    summary_mean_canvas->SetRightMargin(0.1);

    // Calculate dynamic y-axis range for mean ratios
    double min_mean_value = 1e6;  // Start with a large number
    double max_mean_value = -1e6; // Start with a small number
    
    // Find the actual min/max values from the data
    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < nPtBins; bin++) {
            double value = final_ratio_means[m][bin];
            double error = final_ratio_mean_errs[m][bin];
            
            // Consider error bars in the range calculation
            double lower_bound = value - error;
            double upper_bound = value + error;
            
            if (lower_bound < min_mean_value) min_mean_value = lower_bound;
            if (upper_bound > max_mean_value) max_mean_value = upper_bound;
        }
    }
    
    // Add 10% padding to the range
    double mean_range = max_mean_value - min_mean_value;
    double mean_padding = mean_range * 0.1;
    double y_min_mean = min_mean_value - mean_padding;
    double y_max_mean = max_mean_value + mean_padding;
    
    // Ensure the range includes 1.0 (reference line) with some margin
    if (y_min_mean > 0.95) y_min_mean = 0.95;
    if (y_max_mean < 1.05) y_max_mean = 1.05;
    
    // Create frame histogram for mean ratios
    TH1F* frame_mean = new TH1F("frame_mean", "", nPtBins, ptBins);
    frame_mean->SetMinimum(y_min_mean);
    frame_mean->SetMaximum(y_max_mean);
    frame_mean->GetXaxis()->SetTitle("p_{T} [GeV]");
    frame_mean->GetYaxis()->SetTitle("#frac{JES_{bJR10v00}^{Data}/JES_{bJR10v00}^{MC}}{JES^{Data}/JES^{MC}}");
    frame_mean->GetXaxis()->SetTitleSize(0.04);
    frame_mean->GetYaxis()->SetTitleSize(0.04);
    frame_mean->GetXaxis()->SetLabelSize(0.04);
    frame_mean->GetYaxis()->SetLabelSize(0.04);
    frame_mean->GetYaxis()->SetTitleOffset(1.55); 
    frame_mean->SetStats(0);

    // Set the x-axis range to cover the full pT range
    frame_mean->GetXaxis()->SetRangeUser(ptBins[0], ptBins[nPtBins]);
    
    // Remove automatic tick marks and labels - we'll add them manually
    frame_mean->GetXaxis()->SetNdivisions(0, 0, 0, kFALSE);
    frame_mean->GetXaxis()->SetLabelSize(0); // Hide automatic labels
    frame_mean->GetXaxis()->SetLabelOffset(0.01);
    
    frame_mean->Draw("AXIS");
    
    // Draw custom tick marks at bin edges
    double ymin = frame_mean->GetMinimum();
    double ymax = frame_mean->GetMaximum();
    double yrange = ymax - ymin;
    double tick_length = 0.01 * yrange;
    
    TLatex* latex_mean = new TLatex();
    latex_mean->SetTextSize(0.04);
    latex_mean->SetTextAlign(23); // Center, top alignment
    latex_mean->SetTextColor(kBlack);
    
    // Draw tick marks and labels for each bin edge
    for (int i = 0; i <= nPtBins; i++) {
        // Draw tick mark
        TLine* tick = new TLine(ptBins[i], ymin, ptBins[i], ymin + tick_length);
        tick->SetLineColor(kBlack);
        tick->SetLineWidth(1);
        tick->Draw("same");
        
        // Draw label slightly below the axis
        latex_mean->DrawLatex(ptBins[i], ymin - 0.05 * yrange, Form("%.0f", ptBins[i]));
    }
    
    // Add horizontal line at y=1
    TLine* line_mean = new TLine(ptBins[0], 1.0, ptBins[nPtBins], 1.0);
    line_mean->SetLineStyle(2);
    line_mean->SetLineColor(kBlack);
    line_mean->SetLineWidth(2);
    line_mean->Draw("same");
    
    // Create legends - one for mass bins (left) and one for process (right)
    TLegend* legend_mean = new TLegend(0.15, 0.7, 0.6, 0.87);
    legend_mean->SetBorderSize(0);
    legend_mean->SetFillStyle(0);
    legend_mean->SetTextSize(0.035);

    TLegend* legend_process_mean = new TLegend(0.6, 0.75, 0.89, 0.89);
    legend_process_mean->SetBorderSize(0);
    legend_process_mean->SetFillStyle(0);
    legend_process_mean->SetTextSize(0.035);
    legend_process_mean->AddEntry((TObject*)0, "Z (#rightarrow bb) + Jets", "");
    
    // Create TGraphErrors for each mass bin
    TGraphErrors* graphs_mean[nMassCuts];
    
    for (int m = 0; m < nMassCuts; m++) {
        double x_points[nPtBins], y_points[nPtBins], x_errors[nPtBins], y_errors[nPtBins];
        
        for (int bin = 0; bin < nPtBins; bin++) {
            x_points[bin] = ptCenters[bin] + offsets[m];
            y_points[bin] = final_ratio_means[m][bin];
            x_errors[bin] = 0; // No x-error
            y_errors[bin] = final_ratio_mean_errs[m][bin];
        }
        
        graphs_mean[m] = new TGraphErrors(nPtBins, x_points, y_points, x_errors, y_errors);
        graphs_mean[m]->SetMarkerColor(colors[m]);
        graphs_mean[m]->SetLineColor(colors[m]);
        graphs_mean[m]->SetMarkerStyle(markers[m]);
        graphs_mean[m]->SetMarkerSize(1.2);
        graphs_mean[m]->SetLineWidth(2);
        graphs_mean[m]->Draw("P same");
        
        legend_mean->AddEntry(graphs_mean[m], Form("Mass Bin: %s", massCutLabels[m]), "p");
    }
    
    legend_mean->Draw();
    legend_process_mean->Draw(); // Draw the process legend
    DrawAtlasLabel();

    summary_mean_canvas->SaveAs(Form("../Plots/MJB_Ratio/%s/PlotsRdbpTbin/Summary/Summary_Mean_Ratios_vs_pT_%s.pdf", fileWP.c_str(), fileWP.c_str()));
    std::cout << "Summary mean ratios plot generated successfully." << std::endl;
    
    // --- Summary Plot for Resolution Ratios ---
    TCanvas* summary_sigma_canvas = new TCanvas("summary_sigma_canvas", 
        Form("Summary: Data/Combined MC Resolution Ratios (Post_bJR10V00/Pre_bJR10V00) - WP: %s", labelWP_display.c_str()), 
        800, 600);
    summary_sigma_canvas->SetLeftMargin(0.15);
    summary_sigma_canvas->SetBottomMargin(0.15); // Increased for x-axis labels
    summary_sigma_canvas->SetTopMargin(0.12);
    summary_sigma_canvas->SetRightMargin(0.1);

    // Calculate dynamic y-axis range for resolution ratios
    double min_sigma_value = 1e6;  // Start with a large number
    double max_sigma_value = -1e6; // Start with a small number
    
    // Find the actual min/max values from the data
    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < nPtBins; bin++) {
            double value = final_ratio_sigmas[m][bin];
            double error = final_ratio_sigma_errs[m][bin];
            
            // Consider error bars in the range calculation
            double lower_bound = value - error;
            double upper_bound = value + error;
            
            if (lower_bound < min_sigma_value) min_sigma_value = lower_bound;
            if (upper_bound > max_sigma_value) max_sigma_value = upper_bound;
        }
    }
    
    // Add 10% padding to the range
    double sigma_range = max_sigma_value - min_sigma_value;
    double sigma_padding = sigma_range * 0.1;
    double y_min_sigma = min_sigma_value - sigma_padding;
    double y_max_sigma = max_sigma_value + sigma_padding;
    
    // Ensure the range includes 1.0 (reference line) with some margin
    if (y_min_sigma > 0.95) y_min_sigma = 0.95;
    if (y_max_sigma < 1.05) y_max_sigma = 1.05;
    
    // Create frame histogram for resolution ratios
    TH1F* frame_sigma = new TH1F("frame_sigma", "", nPtBins, ptBins);
    frame_sigma->SetMinimum(y_min_sigma);
    frame_sigma->SetMaximum(y_max_sigma);
    frame_sigma->GetXaxis()->SetTitle("p_{T} [GeV]");
    frame_sigma->GetYaxis()->SetTitle("#frac{JER_{bJR00v10}^{Data}/JER_{bJR00v10}^{MC}}{JER^{Data}/JER^{MC}}");
    frame_sigma->GetXaxis()->SetTitleSize(0.04);
    frame_sigma->GetYaxis()->SetTitleSize(0.04);
    frame_sigma->GetXaxis()->SetLabelSize(0.04);
    frame_sigma->GetYaxis()->SetLabelSize(0.04);
    frame_sigma->GetYaxis()->SetTitleOffset(1.55);
    frame_sigma->SetStats(0);

    // Set the x-axis range to cover the full pT range
    frame_sigma->GetXaxis()->SetRangeUser(ptBins[0], ptBins[nPtBins]);
    
    // Remove automatic tick marks and labels - we'll add them manually
    frame_sigma->GetXaxis()->SetNdivisions(0, 0, 0, kFALSE);
    frame_sigma->GetXaxis()->SetLabelSize(0); // Hide automatic labels
    frame_sigma->GetXaxis()->SetLabelOffset(0.01);
    
    frame_sigma->Draw("AXIS");
    
    // Draw custom tick marks at bin edges
    double ymin_sigma = frame_sigma->GetMinimum();
    double ymax_sigma = frame_sigma->GetMaximum();
    double yrange_sigma = ymax_sigma - ymin_sigma;
    double tick_length_sigma = 0.01 * yrange_sigma;
    
    TLatex* latex_sigma = new TLatex();
    latex_sigma->SetTextSize(0.04);
    latex_sigma->SetTextAlign(23); // Center, top alignment
    latex_sigma->SetTextColor(kBlack);
    
    // Draw tick marks and labels for each bin edge
    for (int i = 0; i <= nPtBins; i++) {
        // Draw tick mark
        TLine* tick_sigma = new TLine(ptBins[i], ymin_sigma, ptBins[i], ymin_sigma + tick_length_sigma);
        tick_sigma->SetLineColor(kBlack);
        tick_sigma->SetLineWidth(1);
        tick_sigma->Draw("same");
        
        // Draw label slightly below the axis
        latex_sigma->DrawLatex(ptBins[i], ymin_sigma - 0.05 * yrange_sigma, Form("%.0f", ptBins[i]));
    }

    // Add horizontal line at y=1
    TLine* line_sigma = new TLine(ptBins[0], 1.0, ptBins[nPtBins], 1.0);
    line_sigma->SetLineStyle(2);
    line_sigma->SetLineColor(kBlack);
    line_sigma->SetLineWidth(2);
    line_sigma->Draw("same");
    
    // Create legends - one for mass bins (left) and one for process (right)
    TLegend* legend_sigma = new TLegend(0.15, 0.7, 0.6, 0.87);
    legend_sigma->SetBorderSize(0);
    legend_sigma->SetFillStyle(0);
    legend_sigma->SetTextSize(0.035);
    
    TLegend* legend_process_sigma = new TLegend(0.6, 0.75, 0.89, 0.89);
    legend_process_sigma->SetBorderSize(0);
    legend_process_sigma->SetFillStyle(0);
    legend_process_sigma->SetTextSize(0.035);
    legend_process_sigma->AddEntry((TObject*)0, "Z (#rightarrow bb) + Jets", "");

    // Create TGraphErrors for each mass bin
    TGraphErrors* graphs_sigma[nMassCuts];

    for (int m = 0; m < nMassCuts; m++) {
        double x_points[nPtBins], y_points[nPtBins], x_errors[nPtBins], y_errors[nPtBins];
        
        for (int bin = 0; bin < nPtBins; bin++) {
            x_points[bin] = ptCenters[bin] + offsets[m];
            y_points[bin] = final_ratio_sigmas[m][bin];
            x_errors[bin] = 0; // No x-error
            y_errors[bin] = final_ratio_sigma_errs[m][bin];
        }
        
        graphs_sigma[m] = new TGraphErrors(nPtBins, x_points, y_points, x_errors, y_errors);
        graphs_sigma[m]->SetMarkerColor(colors[m]);
        graphs_sigma[m]->SetLineColor(colors[m]);
        graphs_sigma[m]->SetMarkerStyle(markers[m]);
        graphs_sigma[m]->SetMarkerSize(1.2);
        graphs_sigma[m]->SetLineWidth(2);
        graphs_sigma[m]->Draw("P same");
        
        legend_sigma->AddEntry(graphs_sigma[m], Form("Mass Bin: %s", massCutLabels[m]), "p");
    }
    
    legend_sigma->Draw();
    legend_process_sigma->Draw(); // Draw the process legend
    DrawAtlasLabel();
    
    summary_sigma_canvas->SaveAs(Form("../Plots/MJB_Ratio/%s/PlotsRdbpTbin/Summary/Summary_Resolution_Ratios_vs_pT_%s.pdf", fileWP.c_str(), fileWP.c_str()));
    std::cout << "Summary resolution ratios plot generated successfully." << std::endl;
    
    // Clean up
    delete summary_mean_canvas;
    delete summary_sigma_canvas;
    delete frame_mean;
    delete frame_sigma;
    delete line_mean;
    delete line_sigma;
    delete legend_mean;
    delete legend_sigma;
    delete legend_process_mean;    // Added cleanup
    delete legend_process_sigma;   // Added cleanup
    delete latex_mean;             // Added cleanup
    delete latex_sigma;            // Added cleanup
    for (int m = 0; m < nMassCuts; m++) {
        delete graphs_mean[m];
        delete graphs_sigma[m];
    }
}

/**
 * Helper function to process R_DB data, fill histograms, perform Gaussian fits,
 * and extract mean and resolution (sigma) values and their errors.
 *
 * @param dataPath Path to the data ROOT file.
 * @param zbbPath Path to the Zbb MC ROOT file.
 * @param mjPath Path to the Multijets ROOT file.
 * @param ptBranchName Name of the pT branch in the trees.
 * @param massBranchName Name of the mass branch in the trees.
 * @param data_means Output array for data R_DB means.
 * @param data_mean_errs Output array for errors on data R_DB means.
 * @param mc_means Output array for MC R_DB means.
 * @param mc_mean_errs Output array for errors on MC R_DB means.
 * @param bkg_means Output array for background R_DB means.
 * @param bkg_mean_errs Output array for errors on background R_DB means.
 * @param combined_means Output array for combined MC R_DB means.
 * @param combined_mean_errs Output array for errors on combined MC R_DB means.
 * @param data_sigmas Output array for data R_DB resolutions.
 * @param data_sigma_errs Output array for errors on data R_DB resolutions.
 * @param mc_sigmas Output array for MC R_DB resolutions.
 * @param mc_sigma_errs Output array for errors on MC R_DB resolutions.
 * @param bkg_sigmas Output array for background R_DB resolutions.
 * @param bkg_sigma_errs Output array for errors on background R_DB resolutions.
 * @param combined_sigmas Output array for combined MC R_DB resolutions.
 * @param combined_sigma_errs Output array for errors on combined MC R_DB resolutions.
 */
void ExtractRdbMeans(
    const std::string& dataPath, const std::string& zbbPath, const std::string& mjPath,
    const char* ptBranchName, const char* massBranchName,
    double data_means[nMassCuts][nPtBins], double data_mean_errs[nMassCuts][nPtBins],
    double mc_means[nMassCuts][nPtBins], double mc_mean_errs[nMassCuts][nPtBins],
    double bkg_means[nMassCuts][nPtBins], double bkg_mean_errs[nMassCuts][nPtBins],
    double combined_means[nMassCuts][nPtBins], double combined_mean_errs[nMassCuts][nPtBins],
    double data_sigmas[nMassCuts][nPtBins], double data_sigma_errs[nMassCuts][nPtBins],
    double mc_sigmas[nMassCuts][nPtBins], double mc_sigma_errs[nMassCuts][nPtBins],
    double bkg_sigmas[nMassCuts][nPtBins], double bkg_sigma_errs[nMassCuts][nPtBins],
    double combined_sigmas[nMassCuts][nPtBins], double combined_sigma_errs[nMassCuts][nPtBins]
) {
    // Open ROOT files
    TFile* fdata = TFile::Open(dataPath.c_str(), "READ");
    TFile* fmc = TFile::Open(zbbPath.c_str(), "READ");
    TFile* fmultijets = TFile::Open(mjPath.c_str(), "READ");

    if (!fdata || fdata->IsZombie()) {
        std::cerr << "Error: Could not open data file: " << dataPath << std::endl;
        if (fdata) delete fdata;
        if (fmc) delete fmc;
        if (fmultijets) delete fmultijets;
        return;
    }
    if (!fmc || fmc->IsZombie()) {
        std::cerr << "Error: Could not open MC (Zbb) file: " << zbbPath << std::endl;
        if (fdata) delete fdata;
        if (fmc) delete fmc;
        if (fmultijets) delete fmultijets;
        return;
    }
    if (!fmultijets || fmultijets->IsZombie()) {
        std::cerr << "Error: Could not open Multijets file: " << mjPath << std::endl;
        if (fdata) delete fdata;
        if (fmc) delete fmc;
        if (fmultijets) delete fmultijets;
        return;
    }

    TTree *tree_mc = (TTree*)fmc->Get("nominal");
    TTree *tree_data = (TTree*)fdata->Get("nominal");
    TTree *tree_mulijets = (TTree*)fmultijets->Get("nominal");

    if (!tree_mc) {
        std::cerr << "Error: Could not retrieve 'nominal' tree from MC file." << std::endl;
        fdata->Close(); fmc->Close(); fmultijets->Close();
        delete fdata; delete fmc; delete fmultijets;
        return;
    }
    if (!tree_data) {
        std::cerr << "Error: Could not retrieve 'nominal' tree from data file." << std::endl;
        fdata->Close(); fmc->Close(); fmultijets->Close();
        delete fdata; delete fmc; delete fmultijets;
        return;
    }
    if (!tree_mulijets) {
        std::cerr << "Error: Could not retrieve 'nominal' tree from multijets file." << std::endl;
        fdata->Close(); fmc->Close(); fmultijets->Close();
        delete fdata; delete fmc; delete fmultijets;
        return;
    }

    // Create histograms for each Mass Bin and pT bin
    TH1F* h_mc_R_DB_ptbin[nMassCuts][nPtBins];
    TH1F* h_data_R_DB_ptbin[nMassCuts][nPtBins];
    TH1F* h_multijets_R_DB_ptbin[nMassCuts][nPtBins];
    TH1F* h_combined_MC_R_DB_ptbin[nMassCuts][nPtBins];
    
    // Initialize histograms for each Mass Bin and pT bin
    for (int m = 0; m < nMassCuts; m++) {
        for (int i = 0; i < nPtBins; i++) {
            TString histName_mc = Form("h_mc_R_DB_masscut%d_ptbin%d_%s", m, i, ptBranchName);
            TString histTitle = Form("R_DB for %s (Mass Bin: %s)", ptBinLabels[i], massCutLabels[m]);
            h_mc_R_DB_ptbin[m][i] = new TH1F(histName_mc, histTitle, 50, 0, 2);
            h_mc_R_DB_ptbin[m][i]->Sumw2(); // Enable error calculation

            TString histName_data = Form("h_data_R_DB_masscut%d_ptbin%d_%s", m, i, ptBranchName);
            h_data_R_DB_ptbin[m][i] = new TH1F(histName_data, histTitle, 50, 0, 2);
            h_data_R_DB_ptbin[m][i]->Sumw2(); // Enable error calculation
            
            TString histName_multijets = Form("h_multijets_R_DB_masscut%d_ptbin%d_%s", m, i, ptBranchName);
            h_multijets_R_DB_ptbin[m][i] = new TH1F(histName_multijets, histTitle, 50, 0, 2);
            h_multijets_R_DB_ptbin[m][i]->Sumw2(); // Enable error calculation
            
            TString histName_combined = Form("h_combined_MC_R_DB_masscut%d_ptbin%d_%s", m, i, ptBranchName);
            h_combined_MC_R_DB_ptbin[m][i] = new TH1F(histName_combined, histTitle, 50, 0, 2);
            h_combined_MC_R_DB_ptbin[m][i]->Sumw2(); // Enable error calculation
        }
    }

    // Declare variables for tree branches
    std::vector<float>* ljet_pt_branch = nullptr;
    std::vector<float>* ljet_mass_branch = nullptr;
    double weight, weight_mulijets;
    double R_DB;

    // Set branch addresses dynamically based on input branch names
    tree_mc->SetBranchAddress(ptBranchName, &ljet_pt_branch);
    tree_data->SetBranchAddress(ptBranchName, &ljet_pt_branch);
    tree_mulijets->SetBranchAddress(ptBranchName, &ljet_pt_branch);

    tree_mc->SetBranchAddress(massBranchName, &ljet_mass_branch);
    tree_data->SetBranchAddress(massBranchName, &ljet_mass_branch);
    tree_mulijets->SetBranchAddress(massBranchName, &ljet_mass_branch);

    tree_mc->SetBranchAddress("R_DB", &R_DB);
    tree_data->SetBranchAddress("R_DB", &R_DB);
    tree_mulijets->SetBranchAddress("R_DB", &R_DB);

    tree_mc->SetBranchAddress("total_weight", &weight);
    tree_mulijets->SetBranchAddress("total_weight", &weight_mulijets);

    // Fill MC signal histograms
    Long64_t nentries = tree_mc->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree_mc->GetEntry(i);
        if (ljet_pt_branch->empty() || ljet_mass_branch->empty()) continue; // Skip if no jets
        double leadingJetPt = (*ljet_pt_branch)[0]/1000.0; // Convert MeV to GeV
        double jetMass = (*ljet_mass_branch)[0]/1000.0;      // Convert MeV to GeV
        
        for (int bin = 0; bin < nPtBins; bin++) {
            if (leadingJetPt >= ptBins[bin] && leadingJetPt < ptBins[bin+1]) {
                for (int m = 0; m < nMassCuts; m++) {
                    bool passesCut = true;
                    if (massCutLow[m] > 0 && jetMass < massCutLow[m]) passesCut = false;
                    if (massCutHigh[m] > 0 && jetMass > massCutHigh[m]) passesCut = false;
                    if (passesCut) h_mc_R_DB_ptbin[m][bin]->Fill(R_DB, weight);
                }
                break;
            }
        }
    }

    // Fill data histograms
    nentries = tree_data->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree_data->GetEntry(i);
        if (ljet_pt_branch->empty() || ljet_mass_branch->empty()) continue; // Skip if no jets
        double leadingJetPt = (*ljet_pt_branch)[0]/1000.0;
        double jetMass = (*ljet_mass_branch)[0]/1000.0;
        
        for (int bin = 0; bin < nPtBins; bin++) {
            if (leadingJetPt >= ptBins[bin] && leadingJetPt < ptBins[bin+1]) {
                for (int m = 0; m < nMassCuts; m++) {
                    bool passesCut = true;
                    if (massCutLow[m] > 0 && jetMass < massCutLow[m]) passesCut = false;
                    if (massCutHigh[m] > 0 && jetMass > massCutHigh[m]) passesCut = false;
                    if (passesCut) h_data_R_DB_ptbin[m][bin]->Fill(R_DB);
                }
                break;
            }
        }
    }

    // Fill background histograms
    nentries = tree_mulijets->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree_mulijets->GetEntry(i);
        if (ljet_pt_branch->empty() || ljet_mass_branch->empty()) continue; // Skip if no jets
        double leadingJetPt = (*ljet_pt_branch)[0]/1000.0;
        double jetMass = (*ljet_mass_branch)[0]/1000.0;
        
        for (int bin = 0; bin < nPtBins; bin++) {
            if (leadingJetPt >= ptBins[bin] && leadingJetPt < ptBins[bin+1]) {
                for (int m = 0; m < nMassCuts; m++) {
                    bool passesCut = true;
                    if (massCutLow[m] > 0 && jetMass < massCutLow[m]) passesCut = false;
                    if (massCutHigh[m] > 0 && jetMass > massCutHigh[m]) passesCut = false;
                    if (passesCut) h_multijets_R_DB_ptbin[m][bin]->Fill(R_DB, weight_mulijets);
                }
                break;
            }
        }
    }

    // Calculate scale factors for each Mass Bin separately
    for (int m = 0; m < nMassCuts; m++) {
        TH1F *h_mc_ljet_mass_sf = new TH1F(Form("h_mc_ljet_mass_sf_masscut%d_%s", m, ptBranchName), "ljet_mass", 20, 50, 150);
        h_mc_ljet_mass_sf->Sumw2();
        TH1F *h_data_ljet_mass_sf = new TH1F(Form("h_data_ljet_mass_sf_masscut%d_%s", m, ptBranchName), "ljet_mass", 20, 50, 150);
        h_data_ljet_mass_sf->Sumw2();
        TH1F *h_multijets_ljet_mass_sf = new TH1F(Form("h_multijets_ljet_mass_sf_masscut%d_%s", m, ptBranchName), "ljet_mass", 20, 50, 150);
        h_multijets_ljet_mass_sf->Sumw2();

        // Fill SF histograms
        nentries = tree_mc->GetEntries();
        for (Long64_t i = 0; i < nentries; i++) {
            tree_mc->GetEntry(i);
            if (ljet_mass_branch->empty()) continue;
            double jetMass = (*ljet_mass_branch)[0]/1000.0;
            bool passesCut = true;
            if (massCutLow[m] > 0 && jetMass < massCutLow[m]) passesCut = false;
            if (massCutHigh[m] > 0 && jetMass > massCutHigh[m]) passesCut = false;
            if (passesCut) h_mc_ljet_mass_sf->Fill(jetMass, weight);
        }
        nentries = tree_data->GetEntries();
        for (Long64_t i = 0; i < nentries; i++) {
            tree_data->GetEntry(i);
            if (ljet_mass_branch->empty()) continue;
            double jetMass = (*ljet_mass_branch)[0]/1000.0;
            bool passesCut = true;
            if (massCutLow[m] > 0 && jetMass < massCutLow[m]) passesCut = false;
            if (massCutHigh[m] > 0 && jetMass > massCutHigh[m]) passesCut = false;
            if (passesCut) h_data_ljet_mass_sf->Fill(jetMass);
        }
        nentries = tree_mulijets->GetEntries();
        for (Long64_t i = 0; i < nentries; i++) {
            tree_mulijets->GetEntry(i);
            if (ljet_mass_branch->empty()) continue;
            double jetMass = (*ljet_mass_branch)[0]/1000.0;
            bool passesCut = true;
            if (massCutLow[m] > 0 && jetMass < massCutLow[m]) passesCut = false;
            if (massCutHigh[m] > 0 && jetMass > massCutHigh[m]) passesCut = false;
            if (passesCut) h_multijets_ljet_mass_sf->Fill(jetMass, weight_mulijets);
        }
        
        int binLow1 = h_multijets_ljet_mass_sf->FindBin(50);
        int binHigh1 = h_multijets_ljet_mass_sf->FindBin(65);
        int binLow2 = h_multijets_ljet_mass_sf->FindBin(110);
        int binHigh2 = h_multijets_ljet_mass_sf->FindBin(150);
        
        std::vector<double> x_values, y_values, y_errors;
        for (int bin = binLow1; bin <= binHigh1; ++bin) {
            double bin_center = h_multijets_ljet_mass_sf->GetBinCenter(bin);
            double data_content = h_data_ljet_mass_sf->GetBinContent(bin);
            double mc_content = h_multijets_ljet_mass_sf->GetBinContent(bin);
            if (mc_content > 0) {
                x_values.push_back(bin_center);
                y_values.push_back(data_content / mc_content);
                double data_error = h_data_ljet_mass_sf->GetBinError(bin);
                double mc_error = h_multijets_ljet_mass_sf->GetBinError(bin);
                double ratio_error = sqrt(pow(data_error/mc_content, 2) + pow(data_content*mc_error/(mc_content*mc_content), 2));
                y_errors.push_back(ratio_error);
            }
        }
        for (int bin = binLow2; bin <= binHigh2; ++bin) {
            double bin_center = h_multijets_ljet_mass_sf->GetBinCenter(bin);
            double data_content = h_data_ljet_mass_sf->GetBinContent(bin);
            double mc_content = h_multijets_ljet_mass_sf->GetBinContent(bin);
            if (mc_content > 0) {
                x_values.push_back(bin_center);
                y_values.push_back(data_content / mc_content);
                double data_error = h_data_ljet_mass_sf->GetBinError(bin);
                double mc_error = h_multijets_ljet_mass_sf->GetBinError(bin);
                double ratio_error = sqrt(pow(data_error/mc_content, 2) + pow(data_content*mc_error/(mc_content*mc_content), 2));
                y_errors.push_back(ratio_error);
            }
        }
        
        TGraphErrors* ratio_graph = new TGraphErrors(x_values.size());
        for (size_t i = 0; i < x_values.size(); ++i) {
            ratio_graph->SetPoint(i, x_values[i], y_values[i]);
            ratio_graph->SetPointError(i, 0, y_errors[i]);
        }
        
        TF1* fit_func = new TF1(Form("fit_func_masscut%d_%s", m, ptBranchName), "[0] + [1] * x", 50, 150);
        fit_func->SetParameters(1.0, 0.0);
        ratio_graph->Fit(Form("fit_func_masscut%d_%s", m, ptBranchName), "RQ");
        
        double p0 = fit_func->GetParameter(0);
        double p1 = fit_func->GetParameter(1);
        
        double avgScaleFactor = 0.0;
        int nBinsAvg = 0;
        for (int bin = 1; bin <= h_multijets_ljet_mass_sf->GetNbinsX(); ++bin) {
            double bin_center = h_multijets_ljet_mass_sf->GetBinCenter(bin);
            avgScaleFactor += (p0 + p1 * bin_center);
            nBinsAvg++;
        }
        if (nBinsAvg > 0) {
            avgScaleFactor /= nBinsAvg;
        } else {
            avgScaleFactor = 1.0; // Default to 1 if no bins for averaging
        }
        
        for (int bin = 0; bin < nPtBins; bin++) {
            h_multijets_R_DB_ptbin[m][bin]->Scale(avgScaleFactor);
            h_combined_MC_R_DB_ptbin[m][bin]->Add(h_mc_R_DB_ptbin[m][bin]);
            h_combined_MC_R_DB_ptbin[m][bin]->Add(h_multijets_R_DB_ptbin[m][bin]);
        }
        
        delete h_mc_ljet_mass_sf;
        delete h_data_ljet_mass_sf;
        delete h_multijets_ljet_mass_sf;
        delete ratio_graph;
        delete fit_func;
    }

    // Perform fits for each pT bin with each Mass Bin
    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < nPtBins; bin++) {
            // MC signal fit
            TF1 *gausFit_mc = new TF1(Form("gausFit_mc_masscut%d_bin%d_%s", m, bin, ptBranchName), "gaus", 0.9, 1.1);
            h_mc_R_DB_ptbin[m][bin]->Fit(gausFit_mc, "RQ");
            mc_means[m][bin] = gausFit_mc->GetParameter(1);
            mc_mean_errs[m][bin] = gausFit_mc->GetParError(1);
            mc_sigmas[m][bin] = gausFit_mc->GetParameter(2); // Store sigma
            mc_sigma_errs[m][bin] = gausFit_mc->GetParError(2); // Store sigma error
            
            // Data fit
            TF1 *gausFit_data = new TF1(Form("gausFit_data_masscut%d_bin%d_%s", m, bin, ptBranchName), "gaus", 0.9, 1.1);
            h_data_R_DB_ptbin[m][bin]->Fit(gausFit_data, "RQ");
            data_means[m][bin] = gausFit_data->GetParameter(1);
            data_mean_errs[m][bin] = gausFit_data->GetParError(1);
            data_sigmas[m][bin] = gausFit_data->GetParameter(2); // Store sigma
            data_sigma_errs[m][bin] = gausFit_data->GetParError(2); // Store sigma error
            
            // Background fit
            TF1 *gausFit_bkg = new TF1(Form("gausFit_bkg_masscut%d_bin%d_%s", m, bin, ptBranchName), "gaus", 0.9, 1.1);
            h_multijets_R_DB_ptbin[m][bin]->Fit(gausFit_bkg, "RQ");
            bkg_means[m][bin] = gausFit_bkg->GetParameter(1);
            bkg_mean_errs[m][bin] = gausFit_bkg->GetParError(1);
            bkg_sigmas[m][bin] = gausFit_bkg->GetParameter(2); // Store sigma
            bkg_sigma_errs[m][bin] = gausFit_bkg->GetParError(2); // Store sigma error
            
            // Combined MC fit
            TF1 *gausFit_combined = new TF1(Form("gausFit_combined_masscut%d_bin%d_%s", m, bin, ptBranchName), "gaus", 0.9, 1.1);
            h_combined_MC_R_DB_ptbin[m][bin]->Fit(gausFit_combined, "RQ");
            combined_means[m][bin] = gausFit_combined->GetParameter(1);
            combined_mean_errs[m][bin] = gausFit_combined->GetParError(1);
            combined_sigmas[m][bin] = gausFit_combined->GetParameter(2); // Store sigma
            combined_sigma_errs[m][bin] = gausFit_combined->GetParError(2); // Store sigma error
            
            delete gausFit_mc;
            delete gausFit_data;
            delete gausFit_bkg;
            delete gausFit_combined;
        }
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
    if (fdata) { fdata->Close(); delete fdata; }
    if (fmc) { fmc->Close(); delete fmc; }
    if (fmultijets) { fmultijets->Close(); delete fmultijets; }
}

/**
 * Generates a summary table of mean and resolution ratios.
 * @param outputFilePath Path to save the summary text file.
 * @param labelWP_display Display label for the current working point.
 * @param ratio_data_means Ratio of data R_DB means (Post_bJR10V00/Pre_bJR10V00).
 * @param ratio_data_mean_errs Errors on ratio_data_means.
 * @param ratio_mc_means Ratio of MC R_DB means (Post_bJR10V00/Pre_bJR10V00).
 * @param ratio_mc_mean_errs Errors on ratio_mc_means.
 * @param ratio_bkg_means Ratio of background R_DB means (Post_bJR10V00/Pre_bJR10V00).
 * @param ratio_bkg_mean_errs Errors on ratio_bkg_means.
 * @param ratio_combined_means Ratio of combined MC R_DB means (Post_bJR10V00/Pre_bJR10V00).
 * @param ratio_combined_mean_errs Errors on ratio_combined_means.
 * @param ratio_data_sigmas Ratio of data R_DB resolutions (Post_bJR10V00/Pre_bJR10V00).
 * @param ratio_data_sigma_errs Errors on ratio_data_sigmas.
 * @param ratio_mc_sigmas Ratio of MC R_DB resolutions (Post_bJR10V00/Pre_bJR10V00).
 * @param ratio_mc_sigma_errs Errors on ratio_mc_sigmas.
 * @param ratio_bkg_sigmas Ratio of background R_DB resolutions (Post_bJR10V00/Pre_bJR10V00).
 * @param ratio_bkg_sigma_errs Errors on ratio_bkg_sigmas.
 * @param ratio_combined_sigmas Ratio of combined MC R_DB resolutions (Post_bJR10V00/Pre_bJR10V00).
 * @param ratio_combined_sigma_errs Errors on ratio_combined_sigmas.
 * @param final_ratio_means Final ratio of (Data/Combined MC) means (Post_bJR10V00/Pre_bJR10V00).
 * @param final_ratio_mean_errs Errors on final_ratio_means.
 * @param final_ratio_sigmas Final ratio of (Data/Combined MC) resolutions (Post_bJR10V00/Pre_bJR10V00).
 * @param final_ratio_sigma_errs Errors on final_ratio_sigmas.
 */
void GenerateSummaryTable(
    const std::string& outputFilePath,
    const std::string& labelWP_display,
    double ratio_data_means[nMassCuts][nPtBins], double ratio_data_mean_errs[nMassCuts][nPtBins],
    double ratio_mc_means[nMassCuts][nPtBins], double ratio_mc_mean_errs[nMassCuts][nPtBins],
    double ratio_bkg_means[nMassCuts][nPtBins], double ratio_bkg_mean_errs[nMassCuts][nPtBins],
    double ratio_combined_means[nMassCuts][nPtBins], double ratio_combined_mean_errs[nMassCuts][nPtBins],
    double ratio_data_sigmas[nMassCuts][nPtBins], double ratio_data_sigma_errs[nMassCuts][nPtBins],
    double ratio_mc_sigmas[nMassCuts][nPtBins], double ratio_mc_sigma_errs[nMassCuts][nPtBins],
    double ratio_bkg_sigmas[nMassCuts][nPtBins], double ratio_bkg_sigma_errs[nMassCuts][nPtBins],
    double ratio_combined_sigmas[nMassCuts][nPtBins], double ratio_combined_sigma_errs[nMassCuts][nPtBins],
    double final_ratio_means[nMassCuts][nPtBins], double final_ratio_mean_errs[nMassCuts][nPtBins],
    double final_ratio_sigmas[nMassCuts][nPtBins], double final_ratio_sigma_errs[nMassCuts][nPtBins]
) {
    std::ofstream summaryFile(outputFilePath.c_str());
    if (!summaryFile.is_open()) {
        std::cerr << "Error: Could not open summary file: " << outputFilePath << std::endl;
        return;
    }

    summaryFile << "========================================================================================================================\n";
    summaryFile << "Summary of R_MJB Ratios (Post_bJR10V00 / Pre_bJR10V00) for Working Point: " << labelWP_display << "\n";
    summaryFile << "========================================================================================================================\n\n";

    summaryFile << std::fixed << std::setprecision(4);

    // Header for Mean Ratios
    summaryFile << "------------------------------------------------------------------------------------------------------------------------\n";
    summaryFile << "Mean Ratios (#LT R_{MJB}^{After} #GT / #LT R_{MJB}^{Before} #GT)\n";
    summaryFile << "------------------------------------------------------------------------------------------------------------------------\n";
    summaryFile << std::setw(15) << "Mass Bin"
                << std::setw(25) << "pT Bin"
                << std::setw(15) << "Data Ratio"
                << std::setw(15) << "MC Ratio"
                << std::setw(15) << "Bkg Ratio"
                << std::setw(15) << "Comb. MC Ratio"
                << std::setw(25) << "Final Ratio (Data/Comb)" << "\n";
    summaryFile << "------------------------------------------------------------------------------------------------------------------------\n";

    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < nPtBins; bin++) {
            summaryFile << std::setw(15) << massCutLabels[m]
                        << std::setw(25) << ptBinLabels[bin]
                        << std::setw(15) << Form("%.4f #pm %.4f", ratio_data_means[m][bin], ratio_data_mean_errs[m][bin])
                        << std::setw(15) << Form("%.4f #pm %.4f", ratio_mc_means[m][bin], ratio_mc_mean_errs[m][bin])
                        << std::setw(15) << Form("%.4f #pm %.4f", ratio_bkg_means[m][bin], ratio_bkg_mean_errs[m][bin])
                        << std::setw(15) << Form("%.4f #pm %.4f", ratio_combined_means[m][bin], ratio_combined_mean_errs[m][bin])
                        << std::setw(25) << Form("%.4f #pm %.4f", final_ratio_means[m][bin], final_ratio_mean_errs[m][bin]) << "\n";
        }
    }
    summaryFile << "------------------------------------------------------------------------------------------------------------------------\n\n";

    // Header for Resolution Ratios
    summaryFile << "------------------------------------------------------------------------------------------------------------------------\n";
    summaryFile << "Resolution Ratios (#sigma_{MJB}^{After} / #sigma_{MJB}^{Before})\n";
    summaryFile << "------------------------------------------------------------------------------------------------------------------------\n";
    summaryFile << std::setw(15) << "Mass Bin"
                << std::setw(25) << "pT Bin"
                << std::setw(15) << "Data Ratio"
                << std::setw(15) << "MC Ratio"
                << std::setw(15) << "Bkg Ratio"
                << std::setw(15) << "Comb. MC Ratio"
                << std::setw(25) << "Final Ratio (Data/Comb)" << "\n";
    summaryFile << "------------------------------------------------------------------------------------------------------------------------\n";

    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < nPtBins; bin++) {
            summaryFile << std::setw(15) << massCutLabels[m]
                        << std::setw(25) << ptBinLabels[bin]
                        << std::setw(15) << Form("%.4f #pm %.4f", ratio_data_sigmas[m][bin], ratio_data_sigma_errs[m][bin])
                        << std::setw(15) << Form("%.4f #pm %.4f", ratio_mc_sigmas[m][bin], ratio_mc_sigma_errs[m][bin])
                        << std::setw(15) << Form("%.4f #pm %.4f", ratio_bkg_sigmas[m][bin], ratio_bkg_sigma_errs[m][bin])
                        << std::setw(15) << Form("%.4f #pm %.4f", ratio_combined_sigmas[m][bin], ratio_combined_sigma_errs[m][bin])
                        << std::setw(25) << Form("%.4f #pm %.4f", final_ratio_sigmas[m][bin], final_ratio_sigma_errs[m][bin]) << "\n";
        }
    }
    summaryFile << "------------------------------------------------------------------------------------------------------------------------\n";
    summaryFile.close();
    std::cout << "Summary table generated successfully: " << outputFilePath << std::endl;
}


void process_single_wp_ratio(const std::string& wp_str) { // Renamed from ratio_Response
    // Determine file-friendly and display-friendly working point labels
    std::string fileWP; // Used for file paths (e.g., "0p46", "1p25")
    std::string labelWP_display; // Used for plot titles and printouts (e.g., "0.46 % QCD Eff. WP")

    if (wp_str == "125") {
        fileWP = "1p25";
        labelWP_display = "1.25 % QCD Eff. WP";
    } else if (wp_str == "155") {
        fileWP = "1p55";
        labelWP_display = "1.55 % QCD Eff. WP";
    } else {
        fileWP = "0p" + wp_str;
        double wpValue = std::stod(wp_str) / 100.0;
        std::stringstream ss;
        ss << std::fixed << std::setprecision(2) << wpValue << " % QCD Eff. WP";
        labelWP_display = ss.str();
    }

    std::cout << "Processing R_MJB Ratios for Working Point: " << labelWP_display << std::endl;

    // Create directories for output
    CreateRatioDirectories(fileWP);

    // --- Define file paths and branch names for "Pre_bJR10V00ration" ---
    std::string dataPath_before = Form("../../NtupleSlim/MJB_%swp_slim/data_%swp_slim.root", fileWP.c_str(), fileWP.c_str());
    std::string zbbPath_before = Form("../../NtupleSlim/MJB_%swp_slim/ZbbJets_%swp_slim.root", fileWP.c_str(), fileWP.c_str());
    std::string mjPath_before = Form("../../NtupleSlim/MJB_%swp_slim/Multijets_%swp_slim.root", fileWP.c_str(), fileWP.c_str());
    const char* ptBranchName_before = "ljet_pt";
    const char* massBranchName_before = "ljet_m";

    // --- Define file paths and branch names for "Post_bJR10V00ration" ---
    std::string dataPath_after = Form("../../NtupleSlim/MJB_Deco_%swp_slim/data_deco_%swp_slim.root", fileWP.c_str(), fileWP.c_str());
    std::string zbbPath_after = Form("../../NtupleSlim/MJB_Deco_%swp_slim/ZbbJets_deco_%swp_slim.root", fileWP.c_str(), fileWP.c_str());
    std::string mjPath_after = Form("../../NtupleSlim/MJB_Deco_%swp_slim/Multijets_deco_%swp_slim.root", fileWP.c_str(), fileWP.c_str());
    const char* ptBranchName_after = "ljet_bJR10v00_pt";
    const char* massBranchName_after = "ljet_bJR10v00_mass";

    // --- Arrays to store results for "Pre_bJR10V00ration" ---
    double data_means_before[nMassCuts][nPtBins], data_mean_errs_before[nMassCuts][nPtBins];
    double mc_means_before[nMassCuts][nPtBins], mc_mean_errs_before[nMassCuts][nPtBins];
    double bkg_means_before[nMassCuts][nPtBins], bkg_mean_errs_before[nMassCuts][nPtBins];
    double combined_means_before[nMassCuts][nPtBins], combined_mean_errs_before[nMassCuts][nPtBins];
    double data_sigmas_before[nMassCuts][nPtBins], data_sigma_errs_before[nMassCuts][nPtBins];
    double mc_sigmas_before[nMassCuts][nPtBins], mc_sigma_errs_before[nMassCuts][nPtBins];
    double bkg_sigmas_before[nMassCuts][nPtBins], bkg_sigma_errs_before[nMassCuts][nPtBins];
    double combined_sigmas_before[nMassCuts][nPtBins], combined_sigma_errs_before[nMassCuts][nPtBins];

    // --- Arrays to store results for "Post_bJR10V00ration" ---
    double data_means_after[nMassCuts][nPtBins], data_mean_errs_after[nMassCuts][nPtBins];
    double mc_means_after[nMassCuts][nPtBins], mc_mean_errs_after[nMassCuts][nPtBins];
    double bkg_means_after[nMassCuts][nPtBins], bkg_mean_errs_after[nMassCuts][nPtBins];
    double combined_means_after[nMassCuts][nPtBins], combined_mean_errs_after[nMassCuts][nPtBins];
    double data_sigmas_after[nMassCuts][nPtBins], data_sigma_errs_after[nMassCuts][nPtBins];
    double mc_sigmas_after[nMassCuts][nPtBins], mc_sigma_errs_after[nMassCuts][nPtBins];
    double bkg_sigmas_after[nMassCuts][nPtBins], bkg_sigma_errs_after[nMassCuts][nPtBins];
    double combined_sigmas_after[nMassCuts][nPtBins], combined_sigma_errs_after[nMassCuts][nPtBins];

    // --- Process data for "Pre_bJR10V00ration" ---
    std::cout << "Processing data Pre_bJR10V00rations..." << std::endl;
    ExtractRdbMeans(dataPath_before, zbbPath_before, mjPath_before,
                    ptBranchName_before, massBranchName_before,
                    data_means_before, data_mean_errs_before,
                    mc_means_before, mc_mean_errs_before,
                    bkg_means_before, bkg_mean_errs_before,
                    combined_means_before, combined_mean_errs_before,
                    data_sigmas_before, data_sigma_errs_before,
                    mc_sigmas_before, mc_sigma_errs_before,
                    bkg_sigmas_before, bkg_sigma_errs_before,
                    combined_sigmas_before, combined_sigma_errs_before);

    // --- Process data for "Post_bJR10V00ration" ---
    std::cout << "Processing data Post_bJR10V00rations..." << std::endl;
    ExtractRdbMeans(dataPath_after, zbbPath_after, mjPath_after,
                    ptBranchName_after, massBranchName_after,
                    data_means_after, data_mean_errs_after,
                    mc_means_after, mc_mean_errs_after,
                    bkg_means_after, bkg_mean_errs_after,
                    combined_means_after, combined_mean_errs_after,
                    data_sigmas_after, data_sigma_errs_after,
                    mc_sigmas_after, mc_sigma_errs_after,
                    bkg_sigmas_after, bkg_sigma_errs_after,
                    combined_sigmas_after, combined_sigma_errs_after);

    // --- Calculate the ratio of R_MJB_mean_values (Post_bJR10V00/Pre_bJR10V00) for Data, MC, Bkg, and Combined MC ---
    double ratio_data_means[nMassCuts][nPtBins];
    double ratio_data_mean_errs[nMassCuts][nPtBins];
    double ratio_mc_means[nMassCuts][nPtBins];
    double ratio_mc_mean_errs[nMassCuts][nPtBins];
    double ratio_bkg_means[nMassCuts][nPtBins];
    double ratio_bkg_mean_errs[nMassCuts][nPtBins];
    double ratio_combined_means[nMassCuts][nPtBins];
    double ratio_combined_mean_errs[nMassCuts][nPtBins];

    // --- Calculate the ratio of R_MJB_sigma_values (Post_bJR10V00/Pre_bJR10V00) for Data, MC, Bkg, and Combined MC ---
    double ratio_data_sigmas[nMassCuts][nPtBins];
    double ratio_data_sigma_errs[nMassCuts][nPtBins];
    double ratio_mc_sigmas[nMassCuts][nPtBins];
    double ratio_mc_sigma_errs[nMassCuts][nPtBins];
    double ratio_bkg_sigmas[nMassCuts][nPtBins];
    double ratio_bkg_sigma_errs[nMassCuts][nPtBins];
    double ratio_combined_sigmas[nMassCuts][nPtBins];
    double ratio_combined_sigma_errs[nMassCuts][nPtBins];

    // New arrays for Data/(Combined MC) means ratios and their errors
    double data_over_combined_before_means[nMassCuts][nPtBins];
    double data_over_combined_before_mean_errs[nMassCuts][nPtBins];
    double data_over_combined_after_means[nMassCuts][nPtBins];
    double data_over_combined_after_mean_errs[nMassCuts][nPtBins];

    // New arrays for Data/(Combined MC) resolutions ratios and their errors
    double data_over_combined_before_sigmas[nMassCuts][nPtBins];
    double data_over_combined_before_sigma_errs[nMassCuts][nPtBins];
    double data_over_combined_after_sigmas[nMassCuts][nPtBins];
    double data_over_combined_after_sigma_errs[nMassCuts][nPtBins];

    // New arrays for the final ratio: (Data/Combined MC)_after / (Data/Combined MC)_before for means
    double final_ratio_means[nMassCuts][nPtBins];
    double final_ratio_mean_errs[nMassCuts][nPtBins];

    // New arrays for the final ratio: (Data/Combined MC)_after / (Data/Combined MC)_before for sigmas
    double final_ratio_sigmas[nMassCuts][nPtBins];
    double final_ratio_sigma_errs[nMassCuts][nPtBins];


    std::cout << "\nCalculating ratios (Post_bJR10V00 / Pre_bJR10V00)..." << std::endl;
    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < nPtBins; bin++) {
            // Data Mean Ratio (Post_bJR10V00 / Pre_bJR10V00)
            if (data_means_before[m][bin] != 0) {
                ratio_data_means[m][bin] = data_means_after[m][bin] / data_means_before[m][bin];
                ratio_data_mean_errs[m][bin] = ratio_data_means[m][bin] * sqrt(
                    pow(data_mean_errs_after[m][bin] / data_means_after[m][bin], 2) +
                    pow(data_mean_errs_before[m][bin] / data_means_before[m][bin], 2)
                );
            } else {
                ratio_data_means[m][bin] = 0;
                ratio_data_mean_errs[m][bin] = 0;
                std::cerr << "Warning: Data mean division by zero for Mass Bin " << massCutLabels[m] << ", pT bin " << ptBinLabels[bin] << std::endl;
            }

            // MC Mean Ratio (Post_bJR10V00 / Pre_bJR10V00)
            if (mc_means_before[m][bin] != 0) {
                ratio_mc_means[m][bin] = mc_means_after[m][bin] / mc_means_before[m][bin];
                double err_A_div_A_sq = pow(mc_mean_errs_before[m][bin] / mc_means_before[m][bin], 2);
                double err_B_div_B_sq = pow(mc_mean_errs_after[m][bin] / mc_means_after[m][bin], 2);
                ratio_mc_mean_errs[m][bin] = ratio_mc_means[m][bin] * sqrt(err_A_div_A_sq + err_B_div_B_sq);
            } else {
                ratio_mc_means[m][bin] = 0.0;
                ratio_mc_mean_errs[m][bin] = 0.0;
                std::cerr << "Warning: MC mean division by zero for Mass Bin " << massCutLabels[m] << ", pT bin " << ptBinLabels[bin] << std::endl;
            }

            // Background Mean Ratio (Post_bJR10V00 / Pre_bJR10V00)
            if (bkg_means_before[m][bin] != 0) {
                ratio_bkg_means[m][bin] = bkg_means_after[m][bin] / bkg_means_before[m][bin];
                double err_A_div_A_sq = pow(bkg_mean_errs_before[m][bin] / bkg_means_before[m][bin], 2);
                double err_B_div_B_sq = pow(bkg_mean_errs_after[m][bin] / bkg_means_after[m][bin], 2);
                ratio_bkg_mean_errs[m][bin] = ratio_bkg_means[m][bin] * sqrt(err_A_div_A_sq + err_B_div_B_sq);
            } else {
                ratio_bkg_means[m][bin] = 0.0;
                ratio_bkg_mean_errs[m][bin] = 0.0;
                std::cerr << "Warning: Background mean division by zero for Mass Bin " << massCutLabels[m] << ", pT bin " << ptBinLabels[bin] << std::endl;
            }

            // Combined MC Mean Ratio (Post_bJR10V00 / Pre_bJR10V00)
            if (combined_means_before[m][bin] != 0) {
                ratio_combined_means[m][bin] = combined_means_after[m][bin] / combined_means_before[m][bin];
                double err_A_div_A_sq = pow(combined_mean_errs_before[m][bin] / combined_means_before[m][bin], 2);
                double err_B_div_B_sq = pow(combined_mean_errs_after[m][bin] / combined_means_after[m][bin], 2);
                ratio_combined_mean_errs[m][bin] = ratio_combined_means[m][bin] * sqrt(err_A_div_A_sq + err_B_div_B_sq);
            } else {
                ratio_combined_means[m][bin] = 0.0;
                ratio_combined_mean_errs[m][bin] = 0.0;
                std::cerr << "Warning: Combined MC mean division by zero for Mass Bin " << massCutLabels[m] << ", pT bin " << ptBinLabels[bin] << std::endl;
            }

            // Data Sigma Ratio (Post_bJR10V00 / Pre_bJR10V00)
            if (data_sigmas_before[m][bin] != 0) {
                ratio_data_sigmas[m][bin] = data_sigmas_after[m][bin] / data_sigmas_before[m][bin];
                double err_A_div_A_sq = pow(data_sigma_errs_before[m][bin] / data_sigmas_before[m][bin], 2);
                double err_B_div_B_sq = pow(data_sigma_errs_after[m][bin] / data_sigmas_after[m][bin], 2);
                ratio_data_sigma_errs[m][bin] = ratio_data_sigmas[m][bin] * sqrt(err_A_div_A_sq + err_B_div_B_sq);
            } else {
                ratio_data_sigmas[m][bin] = 0.0;
                ratio_data_sigma_errs[m][bin] = 0.0;
                std::cerr << "Warning: Data sigma division by zero for Mass Bin " << massCutLabels[m] << ", pT bin " << ptBinLabels[bin] << std::endl;
            }

            // MC Sigma Ratio (Post_bJR10V00 / Pre_bJR10V00)
            if (mc_sigmas_before[m][bin] != 0) {
                ratio_mc_sigmas[m][bin] = mc_sigmas_after[m][bin] / mc_sigmas_before[m][bin];
                double err_A_div_A_sq = pow(mc_sigma_errs_before[m][bin] / mc_sigmas_before[m][bin], 2);
                double err_B_div_B_sq = pow(mc_sigma_errs_after[m][bin] / mc_sigmas_after[m][bin], 2);
                ratio_mc_sigma_errs[m][bin] = ratio_mc_sigmas[m][bin] * sqrt(err_A_div_A_sq + err_B_div_B_sq);
            } else {
                ratio_mc_sigmas[m][bin] = 0.0;
                ratio_mc_sigma_errs[m][bin] = 0.0;
                std::cerr << "Warning: MC sigma division by zero for Mass Bin " << massCutLabels[m] << ", pT bin " << ptBinLabels[bin] << std::endl;
            }

            // Background Sigma Ratio (Post_bJR10V00 / Pre_bJR10V00)
            if (bkg_sigmas_before[m][bin] != 0) {
                ratio_bkg_sigmas[m][bin] = bkg_sigmas_after[m][bin] / bkg_sigmas_before[m][bin];
                double err_A_div_A_sq = pow(bkg_sigma_errs_before[m][bin] / bkg_sigmas_before[m][bin], 2);
                double err_B_div_B_sq = pow(bkg_sigma_errs_after[m][bin] / bkg_sigmas_after[m][bin], 2);
                ratio_bkg_sigma_errs[m][bin] = ratio_bkg_sigmas[m][bin] * sqrt(err_A_div_A_sq + err_B_div_B_sq);
            } else {
                ratio_bkg_sigmas[m][bin] = 0.0;
                ratio_bkg_sigma_errs[m][bin] = 0.0;
                std::cerr << "Warning: Background sigma division by zero for Mass Bin " << massCutLabels[m] << ", pT bin " << ptBinLabels[bin] << std::endl;
            }

            // Combined MC Sigma Ratio (Post_bJR10V00 / Pre_bJR10V00)
            if (combined_sigmas_before[m][bin] != 0) {
                ratio_combined_sigmas[m][bin] = combined_sigmas_after[m][bin] / combined_sigmas_before[m][bin];
                double err_A_div_A_sq = pow(combined_sigma_errs_before[m][bin] / combined_sigmas_before[m][bin], 2);
                double err_B_div_B_sq = pow(combined_sigma_errs_after[m][bin] / combined_sigmas_after[m][bin], 2);
                ratio_combined_sigma_errs[m][bin] = ratio_combined_sigmas[m][bin] * sqrt(err_A_div_A_sq + err_B_div_B_sq);
            } else {
                ratio_combined_sigmas[m][bin] = 0.0;
                ratio_combined_sigma_errs[m][bin] = 0.0;
                std::cerr << "Warning: Combined MC sigma division by zero for Mass Bin " << massCutLabels[m] << ", pT bin " << ptBinLabels[bin] << std::endl;
            }

            // Calculate Data/(Combined MC) for Pre_bJR10V00ration (Means)
            if (combined_means_before[m][bin] != 0) {
                data_over_combined_before_means[m][bin] = data_means_before[m][bin] / combined_means_before[m][bin];
                double err_data_sq = pow(data_mean_errs_before[m][bin] / data_means_before[m][bin], 2);
                double err_combined_sq = pow(combined_mean_errs_before[m][bin] / combined_means_before[m][bin], 2);
                data_over_combined_before_mean_errs[m][bin] = data_over_combined_before_means[m][bin] * sqrt(err_data_sq + err_combined_sq);
            } else {
                data_over_combined_before_means[m][bin] = 0.0;
                data_over_combined_before_mean_errs[m][bin] = 0.0;
                std::cerr << "Warning: Data/Combined MC (Pre_bJR10V00) mean division by zero for Mass Bin " << massCutLabels[m] << ", pT bin " << ptBinLabels[bin] << std::endl;
            }

            // Calculate Data/(Combined MC) for Post_bJR10V00ration (Means)
            if (combined_means_after[m][bin] != 0) {
                data_over_combined_after_means[m][bin] = data_means_after[m][bin] / combined_means_after[m][bin];
                double err_data_sq = pow(data_mean_errs_after[m][bin] / data_means_after[m][bin], 2);
                double err_combined_sq = pow(combined_mean_errs_after[m][bin] / combined_means_after[m][bin], 2);
                data_over_combined_after_mean_errs[m][bin] = data_over_combined_after_means[m][bin] * sqrt(err_data_sq + err_combined_sq);
            } else {
                data_over_combined_after_means[m][bin] = 0.0;
                data_over_combined_after_mean_errs[m][bin] = 0.0;
                std::cerr << "Warning: Data/Combined MC (Post_bJR10V00) mean division by zero for Mass Bin " << massCutLabels[m] << ", pT bin " << ptBinLabels[bin] << std::endl;
            }

            // Calculate Data/(Combined MC) for Pre_bJR10V00ration (Sigmas)
            if (combined_sigmas_before[m][bin] != 0) {
                data_over_combined_before_sigmas[m][bin] = data_sigmas_before[m][bin] / combined_sigmas_before[m][bin];
                double err_data_sq = pow(data_sigma_errs_before[m][bin] / data_sigmas_before[m][bin], 2);
                double err_combined_sq = pow(combined_sigma_errs_before[m][bin] / combined_sigmas_before[m][bin], 2);
                data_over_combined_before_sigma_errs[m][bin] = data_over_combined_before_sigmas[m][bin] * sqrt(err_data_sq + err_combined_sq);
            } else {
                data_over_combined_before_sigmas[m][bin] = 0.0;
                data_over_combined_before_sigma_errs[m][bin] = 0.0;
                std::cerr << "Warning: Data/Combined MC (Pre_bJR10V00) sigma division by zero for Mass Bin " << massCutLabels[m] << ", pT bin " << ptBinLabels[bin] << std::endl;
            }


            // Calculate Data/(Combined MC) for Post_bJR10V00ration (Sigmas)
            if (combined_sigmas_after[m][bin] != 0) {
                data_over_combined_after_sigmas[m][bin] = data_sigmas_after[m][bin] / combined_sigmas_after[m][bin];
                double err_data_sq = pow(data_sigma_errs_after[m][bin] / data_sigmas_after[m][bin], 2);
                double err_combined_sq = pow(combined_sigma_errs_after[m][bin] / combined_sigmas_after[m][bin], 2);
                data_over_combined_after_sigma_errs[m][bin] = data_over_combined_after_sigmas[m][bin] * sqrt(err_data_sq + err_combined_sq);
            } else {
                data_over_combined_after_sigmas[m][bin] = 0.0;
                data_over_combined_after_sigma_errs[m][bin] = 0.0;
                std::cerr << "Warning: Data/Combined MC (Post_bJR10V00) sigma division by zero for Mass Bin " << massCutLabels[m] << ", pT bin " << ptBinLabels[bin] << std::endl;
            }

            // --- Calculate Final Ratios (Data Ratio / Combined MC Ratio) ---
            // where Data_ratio = Post_Data/Pre_Data and CombinedMC_ratio = Post_CombinedMC/Pre_CombinedMC
            //
            // Calculate the final ratio: (Data/Combined MC)_after / (Data/Combined MC)_before (Means)
            if (data_over_combined_before_means[m][bin] != 0) {
                final_ratio_means[m][bin] = data_over_combined_after_means[m][bin] / data_over_combined_before_means[m][bin];
                double err_before_sq = pow(data_over_combined_before_mean_errs[m][bin] / data_over_combined_before_means[m][bin], 2);
                double err_after_sq = pow(data_over_combined_after_mean_errs[m][bin] / data_over_combined_after_means[m][bin], 2);
                final_ratio_mean_errs[m][bin] = final_ratio_means[m][bin] * sqrt(err_before_sq + err_after_sq);
            } else {
                final_ratio_means[m][bin] = 0.0;
                final_ratio_mean_errs[m][bin] = 0.0;
                std::cerr << "Warning: Final mean ratio division by zero for Mass Bin " << massCutLabels[m] << ", pT bin " << ptBinLabels[bin] << std::endl;
            }

            // Calculate the final ratio: (Data/Combined MC)_after / (Data/Combined MC)_before (Sigmas)
            if (data_over_combined_before_sigmas[m][bin] != 0) {
                final_ratio_sigmas[m][bin] = data_over_combined_after_sigmas[m][bin] / data_over_combined_before_sigmas[m][bin];
                double err_before_sq = pow(data_over_combined_before_sigma_errs[m][bin] / data_over_combined_before_sigmas[m][bin], 2);
                double err_after_sq = pow(data_over_combined_after_sigma_errs[m][bin] / data_over_combined_after_sigmas[m][bin], 2);
                final_ratio_sigma_errs[m][bin] = final_ratio_sigmas[m][bin] * sqrt(err_before_sq + err_after_sq);
            } else {
                final_ratio_sigmas[m][bin] = 0.0;
                final_ratio_sigma_errs[m][bin] = 0.0;
                std::cerr << "Warning: Final sigma ratio division by zero for Mass Bin " << massCutLabels[m] << ", pT bin " << ptBinLabels[bin] << std::endl;
            }
        }
    }

    // --- Create and save the 2D plot for the Data Mean Ratio ---
    TCanvas* ratio_data_canvas = new TCanvas("ratio_data_canvas", Form("Data R_{MJB} Mean Ratio (Post_bJR10V00/Pre_bJR10V00) - WP: %s", labelWP_display.c_str()), 800, 600);
    ratio_data_canvas->SetTopMargin(0.12);
    ratio_data_canvas->SetRightMargin(0.15);
    ratio_data_canvas->SetLeftMargin(0.12);
    ratio_data_canvas->SetBottomMargin(0.12);
    ratio_data_canvas->cd();

    TH2F* ratio_data_map = new TH2F("ratio_data_map", ";Mass Bin;p_{T} Bin", nMassCuts, -0.5, nMassCuts-0.5, nPtBins, -0.5, nPtBins-0.5);

    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < nPtBins; bin++) {
            ratio_data_map->SetBinContent(m+1, bin+1, ratio_data_means[m][bin]);
            ratio_data_map->SetBinError(m+1, bin+1, ratio_data_mean_errs[m][bin]);
        }
    }

    ratio_data_map->SetStats(0);
    ratio_data_map->SetMarkerSize(1.2);

    for (int m = 0; m < nMassCuts; m++) {
        ratio_data_map->GetXaxis()->SetBinLabel(m+1, massCutLabels[m]);
    }
    for (int bin = 0; bin < nPtBins; bin++) {
        TString simplifiedLabel = Form("%d-%d GeV", (int)ptBins[bin], (int)ptBins[bin+1]);
        ratio_data_map->GetYaxis()->SetBinLabel(bin+1, simplifiedLabel.Data());
    }

    ratio_data_map->GetXaxis()->SetLabelSize(0.025);
    ratio_data_map->GetYaxis()->SetLabelSize(0.025);
    ratio_data_map->GetXaxis()->SetTitleSize(0.035);
    ratio_data_map->GetYaxis()->SetTitleSize(0.035);
    ratio_data_map->GetZaxis()->SetTitleSize(0.03);
    ratio_data_map->SetTitle("");

    gStyle->SetNumberContours(100);
    SetDynamicRatioPaletteRange(ratio_data_map); // Apply dynamic range

    gStyle->SetPaintTextFormat("4.3f");

    ratio_data_map->Draw("COLZ TEXT");
    gPad->Update();

    TPaletteAxis* palette_data_ratio = (TPaletteAxis*)ratio_data_map->GetListOfFunctions()->FindObject("palette");
    if (palette_data_ratio) {
        palette_data_ratio->SetX1NDC(0.86);
        palette_data_ratio->SetX2NDC(0.89);
        palette_data_ratio->SetY1NDC(0.12);
        palette_data_ratio->SetY2NDC(0.88);
        palette_data_ratio->SetLabelSize(0.03);
        palette_data_ratio->SetTitleOffset(1.3);
        
        TLatex palette_text;
        palette_text.SetNDC();
        palette_text.SetTextSize(0.035);
        palette_text.SetTextAlign(22);
        palette_text.SetTextAngle(90);
        palette_text.DrawLatex(0.97, 0.5, "#LT R_{MJB}^{Data, Post_bJR10V00} #GT / #LT R_{MJB}^{Data, Pre_bJR10V00} #GT");
    }

    DrawAtlasLabel();
    ratio_data_canvas->SaveAs(Form("../Plots/MJB_Ratio/%s/PlotsRdbpTbin/Summary/Data_R_MJB_mean_ratio_Post_Pre_bJR10v00.pdf", fileWP.c_str()));
    std::cout << "\nData R_MJB mean ratio plot generated successfully." << std::endl;

    // --- Create and save the 2D plot for the MC Mean Ratio ---
    TCanvas* ratio_mc_canvas = new TCanvas("ratio_mc_canvas", Form("MC R_{MJB} Mean Ratio (Post_bJR10V00/Pre_bJR10V00) - WP: %s", labelWP_display.c_str()), 800, 600);
    ratio_mc_canvas->SetTopMargin(0.12);
    ratio_mc_canvas->SetRightMargin(0.15);
    ratio_mc_canvas->SetLeftMargin(0.12);
    ratio_mc_canvas->SetBottomMargin(0.12);
    ratio_mc_canvas->cd();

    TH2F* ratio_mc_map = new TH2F("ratio_mc_map", ";Mass Bin;p_{T} Bin", nMassCuts, -0.5, nMassCuts-0.5, nPtBins, -0.5, nPtBins-0.5);

    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < nPtBins; bin++) {
            ratio_mc_map->SetBinContent(m+1, bin+1, ratio_mc_means[m][bin]);
            ratio_mc_map->SetBinError(m+1, bin+1, ratio_mc_mean_errs[m][bin]);
        }
    }

    ratio_mc_map->SetStats(0);
    ratio_mc_map->SetMarkerSize(1.2);

    for (int m = 0; m < nMassCuts; m++) {
        ratio_mc_map->GetXaxis()->SetBinLabel(m+1, massCutLabels[m]);
    }
    for (int bin = 0; bin < nPtBins; bin++) {
        TString simplifiedLabel = Form("%d-%d GeV", (int)ptBins[bin], (int)ptBins[bin+1]);
        ratio_mc_map->GetYaxis()->SetBinLabel(bin+1, simplifiedLabel.Data());
    }

    ratio_mc_map->GetXaxis()->SetLabelSize(0.025);
    ratio_mc_map->GetYaxis()->SetLabelSize(0.025);
    ratio_mc_map->GetXaxis()->SetTitleSize(0.035);
    ratio_mc_map->GetYaxis()->SetTitleSize(0.035);
    ratio_mc_map->GetZaxis()->SetTitleSize(0.03);
    ratio_mc_map->SetTitle("");

    gStyle->SetNumberContours(100);
    SetDynamicRatioPaletteRange(ratio_mc_map); // Apply dynamic range

    gStyle->SetPaintTextFormat("4.3f");

    ratio_mc_map->Draw("COLZ TEXT");
    gPad->Update();

    TPaletteAxis* palette_mc_ratio = (TPaletteAxis*)ratio_mc_map->GetListOfFunctions()->FindObject("palette");
    if (palette_mc_ratio) {
        palette_mc_ratio->SetX1NDC(0.86);
        palette_mc_ratio->SetX2NDC(0.89);
        palette_mc_ratio->SetY1NDC(0.12);
        palette_mc_ratio->SetY2NDC(0.88);
        palette_mc_ratio->SetLabelSize(0.03);
        palette_mc_ratio->SetTitleOffset(1.3);
        
        TLatex palette_text;
        palette_text.SetNDC();
        palette_text.SetTextSize(0.035);
        palette_text.SetTextAlign(22);
        palette_text.SetTextAngle(90);
        palette_text.DrawLatex(0.97, 0.5, "#LT R_{MJB}^{MC, Post_bJR10V00} #GT / #LT R_{MJB}^{MC, Pre_bJR10V00} #GT");
    }

    DrawAtlasLabel();
    ratio_mc_canvas->SaveAs(Form("../Plots/MJB_Ratio/%s/PlotsRdbpTbin/Summary/MC_R_MJB_mean_ratio_Post_Pre_bJR10v00.pdf", fileWP.c_str()));
    std::cout << "MC R_MJB mean ratio plot generated successfully." << std::endl;

    // --- Create and save the 2D plot for the Background Mean Ratio ---
    TCanvas* ratio_bkg_canvas = new TCanvas("ratio_bkg_canvas", Form("Background R_{MJB} Mean Ratio (Post_bJR10V00/Pre_bJR10V00) - WP: %s", labelWP_display.c_str()), 800, 600);
    ratio_bkg_canvas->SetTopMargin(0.12);
    ratio_bkg_canvas->SetRightMargin(0.15);
    ratio_bkg_canvas->SetLeftMargin(0.12);
    ratio_bkg_canvas->SetBottomMargin(0.12);
    ratio_bkg_canvas->cd();

    TH2F* ratio_bkg_map = new TH2F("ratio_bkg_map", ";Mass Bin;p_{T} Bin", nMassCuts, -0.5, nMassCuts-0.5, nPtBins, -0.5, nPtBins-0.5);

    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < nPtBins; bin++) {
            ratio_bkg_map->SetBinContent(m+1, bin+1, ratio_bkg_means[m][bin]);
            ratio_bkg_map->SetBinError(m+1, bin+1, ratio_bkg_mean_errs[m][bin]);
        }
    }

    ratio_bkg_map->SetStats(0);
    ratio_bkg_map->SetMarkerSize(1.2);

    for (int m = 0; m < nMassCuts; m++) {
        ratio_bkg_map->GetXaxis()->SetBinLabel(m+1, massCutLabels[m]);
    }
    for (int bin = 0; bin < nPtBins; bin++) {
        TString simplifiedLabel = Form("%d-%d GeV", (int)ptBins[bin], (int)ptBins[bin+1]);
        ratio_bkg_map->GetYaxis()->SetBinLabel(bin+1, simplifiedLabel.Data());
    }

    ratio_bkg_map->GetXaxis()->SetLabelSize(0.025);
    ratio_bkg_map->GetYaxis()->SetLabelSize(0.025);
    ratio_bkg_map->GetXaxis()->SetTitleSize(0.035);
    ratio_bkg_map->GetYaxis()->SetTitleSize(0.035);
    ratio_bkg_map->GetZaxis()->SetTitleSize(0.03);
    ratio_bkg_map->SetTitle("");

    gStyle->SetNumberContours(100);
    SetDynamicRatioPaletteRange(ratio_bkg_map); // Apply dynamic range

    gStyle->SetPaintTextFormat("4.3f");

    ratio_bkg_map->Draw("COLZ TEXT");
    gPad->Update();

    TPaletteAxis* palette_bkg_ratio = (TPaletteAxis*)ratio_bkg_map->GetListOfFunctions()->FindObject("palette");
    if (palette_bkg_ratio) {
        palette_bkg_ratio->SetX1NDC(0.86);
        palette_bkg_ratio->SetX2NDC(0.89);
        palette_bkg_ratio->SetY1NDC(0.12);
        palette_bkg_ratio->SetY2NDC(0.88);
        palette_bkg_ratio->SetLabelSize(0.03);
        palette_bkg_ratio->SetTitleOffset(1.3);
        
        TLatex palette_text;
        palette_text.SetNDC();
        palette_text.SetTextSize(0.035);
        palette_text.SetTextAlign(22);
        palette_text.SetTextAngle(90);
        palette_text.DrawLatex(0.97, 0.5, "#LT R_{MJB}^{Background, Pre_bJR10V00} #GT / #LT R_{MJB}^{Background, Post_bJR10V00} #GT");
    }

    DrawAtlasLabel();
    ratio_bkg_canvas->SaveAs(Form("../Plots/MJB_Ratio/%s/PlotsRdbpTbin/Summary/Bkg_R_MJB_mean_ratio_Post_Pre_bJR10v00.pdf", fileWP.c_str()));
    std::cout << "Background R_MJB mean ratio plot generated successfully." << std::endl;

    // --- Create and save the 2D plot for the Combined MC Mean Ratio ---
    TCanvas* ratio_combined_canvas = new TCanvas("ratio_combined_canvas", Form("Combined MC R_{MJB} Mean Ratio (Post_bJR10V00/Pre_bJR10V00) - WP: %s", labelWP_display.c_str()), 800, 600);
    ratio_combined_canvas->SetTopMargin(0.12);
    ratio_combined_canvas->SetRightMargin(0.15);
    ratio_combined_canvas->SetLeftMargin(0.12);
    ratio_combined_canvas->SetBottomMargin(0.12);
    ratio_combined_canvas->cd();

    TH2F* ratio_combined_map = new TH2F("ratio_combined_map", ";Mass Bin;p_{T} Bin", nMassCuts, -0.5, nMassCuts-0.5, nPtBins, -0.5, nPtBins-0.5);

    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < nPtBins; bin++) {
            ratio_combined_map->SetBinContent(m+1, bin+1, ratio_combined_means[m][bin]);
            ratio_combined_map->SetBinError(m+1, bin+1, ratio_combined_mean_errs[m][bin]);
        }
    }

    ratio_combined_map->SetStats(0);
    ratio_combined_map->SetMarkerSize(1.2);

    for (int m = 0; m < nMassCuts; m++) {
        ratio_combined_map->GetXaxis()->SetBinLabel(m+1, massCutLabels[m]);
    }
    for (int bin = 0; bin < nPtBins; bin++) {
        TString simplifiedLabel = Form("%d-%d GeV", (int)ptBins[bin], (int)ptBins[bin+1]);
        ratio_combined_map->GetYaxis()->SetBinLabel(bin+1, simplifiedLabel.Data());
    }

    ratio_combined_map->GetXaxis()->SetLabelSize(0.025);
    ratio_combined_map->GetYaxis()->SetLabelSize(0.025);
    ratio_combined_map->GetXaxis()->SetTitleSize(0.035);
    ratio_combined_map->GetYaxis()->SetTitleSize(0.035);
    ratio_combined_map->GetZaxis()->SetTitleSize(0.03);
    ratio_combined_map->SetTitle("");

    gStyle->SetNumberContours(100);
    SetDynamicRatioPaletteRange(ratio_combined_map); // Apply dynamic range

    gStyle->SetPaintTextFormat("4.3f");

    ratio_combined_map->Draw("COLZ TEXT");
    gPad->Update();

    TPaletteAxis* palette_combined_ratio = (TPaletteAxis*)ratio_combined_map->GetListOfFunctions()->FindObject("palette");
    if (palette_combined_ratio) {
        palette_combined_ratio->SetX1NDC(0.86);
        palette_combined_ratio->SetX2NDC(0.89);
        palette_combined_ratio->SetY1NDC(0.12);
        palette_combined_ratio->SetY2NDC(0.88);
        palette_combined_ratio->SetLabelSize(0.03);
        palette_combined_ratio->SetTitleOffset(1.3);
        
        TLatex palette_text;
        palette_text.SetNDC();
        palette_text.SetTextSize(0.035);
        palette_text.SetTextAlign(22);
        palette_text.SetTextAngle(90);
        palette_text.DrawLatex(0.97, 0.5, "#LT R_{MJB}^{Combined MC, Pre_bJR10V00} #GT / #LT R_{MJB}^{Combined MC, Post_bJR10V00} #GT");
    }

    DrawAtlasLabel();
    ratio_combined_canvas->SaveAs(Form("../Plots/MJB_Ratio/%s/PlotsRdbpTbin/Summary/Combined_MC_R_MJB_mean_ratio_Post_Pre_bJR10v00.pdf", fileWP.c_str()));
    std::cout << "Combined MC R_MJB mean ratio plot generated successfully." << std::endl;

    // --- Create and save the 2D plot for the (Data/Combined MC) Mean Ratio (Post_bJR10V00/Pre_bJR10V00) ---
    TCanvas* final_ratio_mean_canvas = new TCanvas("final_ratio_mean_canvas", Form("Data/(Combined MC) Mean Ratio (Post_bJR10V00/Pre_bJR10V00) - WP: %s", labelWP_display.c_str()), 800, 600);
    final_ratio_mean_canvas->SetTopMargin(0.12);
    final_ratio_mean_canvas->SetRightMargin(0.15);
    final_ratio_mean_canvas->SetLeftMargin(0.12);
    final_ratio_mean_canvas->SetBottomMargin(0.12);
    final_ratio_mean_canvas->cd();

    TH2F* final_ratio_mean_map = new TH2F("final_ratio_mean_map", ";Mass Bin;p_{T} Bin", nMassCuts, -0.5, nMassCuts-0.5, nPtBins, -0.5, nPtBins-0.5);

    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < nPtBins; bin++) {
            final_ratio_mean_map->SetBinContent(m+1, bin+1, final_ratio_means[m][bin]);
            final_ratio_mean_map->SetBinError(m+1, bin+1, final_ratio_mean_errs[m][bin]);
        }
    }

    final_ratio_mean_map->SetStats(0);
    final_ratio_mean_map->SetMarkerSize(1.2);

    for (int m = 0; m < nMassCuts; m++) {
        final_ratio_mean_map->GetXaxis()->SetBinLabel(m+1, massCutLabels[m]);
    }
    for (int bin = 0; bin < nPtBins; bin++) {
        TString simplifiedLabel = Form("%d-%d GeV", (int)ptBins[bin], (int)ptBins[bin+1]);
        final_ratio_mean_map->GetYaxis()->SetBinLabel(bin+1, simplifiedLabel.Data());
    }

    final_ratio_mean_map->GetXaxis()->SetLabelSize(0.025);
    final_ratio_mean_map->GetYaxis()->SetLabelSize(0.025);
    final_ratio_mean_map->GetXaxis()->SetTitleSize(0.035);
    final_ratio_mean_map->GetYaxis()->SetTitleSize(0.035);
    final_ratio_mean_map->GetZaxis()->SetTitleSize(0.03);
    final_ratio_mean_map->SetTitle("");

    gStyle->SetNumberContours(100);
    SetDynamicRatioPaletteRange(final_ratio_mean_map); // Apply dynamic range

    gStyle->SetPaintTextFormat("4.3f");

    final_ratio_mean_map->Draw("COLZ TEXT");
    gPad->Update();

    TPaletteAxis* palette_final_ratio_mean = (TPaletteAxis*)final_ratio_mean_map->GetListOfFunctions()->FindObject("palette");
    if (palette_final_ratio_mean) {
        palette_final_ratio_mean->SetX1NDC(0.86);
        palette_final_ratio_mean->SetX2NDC(0.89);
        palette_final_ratio_mean->SetY1NDC(0.12);
        palette_final_ratio_mean->SetY2NDC(0.88);
        palette_final_ratio_mean->SetLabelSize(0.03);
        palette_final_ratio_mean->SetTitleOffset(1.3);
        
        TLatex palette_text;
        palette_text.SetNDC();
        palette_text.SetTextSize(0.035);
        palette_text.SetTextAlign(22);
        palette_text.SetTextAngle(90);
        palette_text.DrawLatex(0.97, 0.5, "(JES_{bJR10v00}^{Data}/JES_{bJR10v00}^{MC})/(JES^{Data}/JES^{MC})");
    }

    DrawAtlasLabel();
    final_ratio_mean_canvas->SaveAs(Form("../Plots/MJB_Ratio/%s/PlotsRdbpTbin/Summary/Data_Over_Combined_MC_Ratio_Post_Pre_bJR10v00_Means.pdf", fileWP.c_str()));
    std::cout << "Data/(Combined MC) mean ratio (Post_bJR10V00/Pre_bJR10V00) plot generated successfully." << std::endl;


    // --- Create and save the 2D plot for the Data Resolution Ratio ---
    TCanvas* ratio_data_sigma_canvas = new TCanvas("ratio_data_sigma_canvas", Form("Data R_{MJB} Resolution Ratio (Post_bJR10V00/Pre_bJR10V00) - WP: %s", labelWP_display.c_str()), 800, 600);
    ratio_data_sigma_canvas->SetTopMargin(0.12);
    ratio_data_sigma_canvas->SetRightMargin(0.15);
    ratio_data_sigma_canvas->SetLeftMargin(0.12);
    ratio_data_sigma_canvas->SetBottomMargin(0.12);
    ratio_data_sigma_canvas->cd();

    TH2F* ratio_data_sigma_map = new TH2F("ratio_data_sigma_map", ";Mass Bin;p_{T} Bin", nMassCuts, -0.5, nMassCuts-0.5, nPtBins, -0.5, nPtBins-0.5);

    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < nPtBins; bin++) {
            ratio_data_sigma_map->SetBinContent(m+1, bin+1, ratio_data_sigmas[m][bin]);
            ratio_data_sigma_map->SetBinError(m+1, bin+1, ratio_data_sigma_errs[m][bin]);
        }
    }

    ratio_data_sigma_map->SetStats(0);
    ratio_data_sigma_map->SetMarkerSize(1.2);

    for (int m = 0; m < nMassCuts; m++) {
        ratio_data_sigma_map->GetXaxis()->SetBinLabel(m+1, massCutLabels[m]);
    }
    for (int bin = 0; bin < nPtBins; bin++) {
        TString simplifiedLabel = Form("%d-%d GeV", (int)ptBins[bin], (int)ptBins[bin+1]);
        ratio_data_sigma_map->GetYaxis()->SetBinLabel(bin+1, simplifiedLabel.Data());
    }

    ratio_data_sigma_map->GetXaxis()->SetLabelSize(0.025);
    ratio_data_sigma_map->GetYaxis()->SetLabelSize(0.025);
    ratio_data_sigma_map->GetXaxis()->SetTitleSize(0.035);
    ratio_data_sigma_map->GetYaxis()->SetTitleSize(0.035);
    ratio_data_sigma_map->GetZaxis()->SetTitleSize(0.03);
    ratio_data_sigma_map->SetTitle("");

    gStyle->SetNumberContours(100);
    SetDynamicRatioPaletteRange(ratio_data_sigma_map); // Apply dynamic range

    gStyle->SetPaintTextFormat("4.3f");

    ratio_data_sigma_map->Draw("COLZ TEXT");
    gPad->Update();

    TPaletteAxis* palette_data_sigma_ratio = (TPaletteAxis*)ratio_data_sigma_map->GetListOfFunctions()->FindObject("palette");
    if (palette_data_sigma_ratio) {
        palette_data_sigma_ratio->SetX1NDC(0.86);
        palette_data_sigma_ratio->SetX2NDC(0.89);
        palette_data_sigma_ratio->SetY1NDC(0.12);
        palette_data_sigma_ratio->SetY2NDC(0.88);
        palette_data_sigma_ratio->SetLabelSize(0.03);
        palette_data_sigma_ratio->SetTitleOffset(1.3);
        
        TLatex palette_text;
        palette_text.SetNDC();
        palette_text.SetTextSize(0.035);
        palette_text.SetTextAlign(22);
        palette_text.SetTextAngle(90);
        palette_text.DrawLatex(0.97, 0.5, "#sigma_{MJB}^{Data, Post_bJR10V00} / #sigma_{MJB}^{Data, Pre_bJR10V00}"); 
    }

    DrawAtlasLabel();
    ratio_data_sigma_canvas->SaveAs(Form("../Plots/MJB_Ratio/%s/PlotsRdbpTbin/Summary/Data_R_MJB_sigma_ratio_Post_Pre_bJR10v00.pdf", fileWP.c_str()));
    std::cout << "Data R_MJB resolution ratio plot generated successfully." << std::endl;

    // --- Create and save the 2D plot for the MC Resolution Ratio ---
    TCanvas* ratio_mc_sigma_canvas = new TCanvas("ratio_mc_sigma_canvas", Form("MC R_{MJB} Resolution Ratio (Post_bJR10V00/Pre_bJR10V00) - WP: %s", labelWP_display.c_str()), 800, 600);
    ratio_mc_sigma_canvas->SetTopMargin(0.12);
    ratio_mc_sigma_canvas->SetRightMargin(0.15);
    ratio_mc_sigma_canvas->SetLeftMargin(0.12);
    ratio_mc_sigma_canvas->SetBottomMargin(0.12);
    ratio_mc_sigma_canvas->cd();

    TH2F* ratio_mc_sigma_map = new TH2F("ratio_mc_sigma_map", ";Mass Bin;p_{T} Bin", nMassCuts, -0.5, nMassCuts-0.5, nPtBins, -0.5, nPtBins-0.5);

    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < nPtBins; bin++) {
            ratio_mc_sigma_map->SetBinContent(m+1, bin+1, ratio_mc_sigmas[m][bin]);
            ratio_mc_sigma_map->SetBinError(m+1, bin+1, ratio_mc_sigma_errs[m][bin]);
        }
    }

    ratio_mc_sigma_map->SetStats(0);
    ratio_mc_sigma_map->SetMarkerSize(1.2);

    for (int m = 0; m < nMassCuts; m++) {
        ratio_mc_sigma_map->GetXaxis()->SetBinLabel(m+1, massCutLabels[m]);
    }
    for (int bin = 0; bin < nPtBins; bin++) {
        TString simplifiedLabel = Form("%d-%d GeV", (int)ptBins[bin], (int)ptBins[bin+1]);
        ratio_mc_sigma_map->GetYaxis()->SetBinLabel(bin+1, simplifiedLabel.Data());
    }

    ratio_mc_sigma_map->GetXaxis()->SetLabelSize(0.025);
    ratio_mc_sigma_map->GetYaxis()->SetLabelSize(0.025);
    ratio_mc_sigma_map->GetXaxis()->SetTitleSize(0.035);
    ratio_mc_sigma_map->GetYaxis()->SetTitleSize(0.035);
    ratio_mc_sigma_map->GetZaxis()->SetTitleSize(0.03);
    ratio_mc_sigma_map->SetTitle("");

    gStyle->SetNumberContours(100);
    SetDynamicRatioPaletteRange(ratio_mc_sigma_map); // Apply dynamic range

    gStyle->SetPaintTextFormat("4.3f");

    ratio_mc_sigma_map->Draw("COLZ TEXT");
    gPad->Update();

    TPaletteAxis* palette_mc_sigma_ratio = (TPaletteAxis*)ratio_mc_sigma_map->GetListOfFunctions()->FindObject("palette");
    if (palette_mc_sigma_ratio) {
        palette_mc_sigma_ratio->SetX1NDC(0.86);
        palette_mc_sigma_ratio->SetX2NDC(0.89);
        palette_mc_sigma_ratio->SetY1NDC(0.12);
        palette_mc_sigma_ratio->SetY2NDC(0.88);
        palette_mc_sigma_ratio->SetLabelSize(0.03);
        palette_mc_sigma_ratio->SetTitleOffset(1.3);
        
        TLatex palette_text;
        palette_text.SetNDC();
        palette_text.SetTextSize(0.035);
        palette_text.SetTextAlign(22);
        palette_text.SetTextAngle(90);
        palette_text.DrawLatex(0.97, 0.5, "#sigma_{MJB}^{MC, Post_bJR10V00} / #sigma_{MJB}^{MC, Pre_bJR10V00}");
    }

    DrawAtlasLabel();
    ratio_mc_sigma_canvas->SaveAs(Form("../Plots/MJB_Ratio/%s/PlotsRdbpTbin/Summary/MC_R_MJB_sigma_ratio_Post_Pre_bJR10v00.pdf", fileWP.c_str()));
    std::cout << "MC R_MJB resolution ratio plot generated successfully." << std::endl;

    // --- Create and save the 2D plot for the Background Resolution Ratio ---
    TCanvas* ratio_bkg_sigma_canvas = new TCanvas("ratio_bkg_sigma_canvas", Form("Background R_{MJB} Resolution Ratio (Post_bJR10V00/Pre_bJR10V00) - WP: %s", labelWP_display.c_str()), 800, 600);
    ratio_bkg_sigma_canvas->SetTopMargin(0.12);
    ratio_bkg_sigma_canvas->SetRightMargin(0.15);
    ratio_bkg_sigma_canvas->SetLeftMargin(0.12);
    ratio_bkg_sigma_canvas->SetBottomMargin(0.12);
    ratio_bkg_sigma_canvas->cd();

    TH2F* ratio_bkg_sigma_map = new TH2F("ratio_bkg_sigma_map", ";Mass Bin;p_{T} Bin", nMassCuts, -0.5, nMassCuts-0.5, nPtBins, -0.5, nPtBins-0.5);

    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < nPtBins; bin++) {
            ratio_bkg_sigma_map->SetBinContent(m+1, bin+1, ratio_bkg_sigmas[m][bin]);
            ratio_bkg_sigma_map->SetBinError(m+1, bin+1, ratio_bkg_sigma_errs[m][bin]);
        }
    }

    ratio_bkg_sigma_map->SetStats(0);
    ratio_bkg_sigma_map->SetMarkerSize(1.2);

    for (int m = 0; m < nMassCuts; m++) {
        ratio_bkg_sigma_map->GetXaxis()->SetBinLabel(m+1, massCutLabels[m]);
    }
    for (int bin = 0; bin < nPtBins; bin++) {
        TString simplifiedLabel = Form("%d-%d GeV", (int)ptBins[bin], (int)ptBins[bin+1]);
        ratio_bkg_sigma_map->GetYaxis()->SetBinLabel(bin+1, simplifiedLabel.Data());
    }

    ratio_bkg_sigma_map->GetXaxis()->SetLabelSize(0.025);
    ratio_bkg_sigma_map->GetYaxis()->SetLabelSize(0.025);
    ratio_bkg_sigma_map->GetXaxis()->SetTitleSize(0.035);
    ratio_bkg_sigma_map->GetYaxis()->SetTitleSize(0.035);
    ratio_bkg_sigma_map->GetZaxis()->SetTitleSize(0.03);
    ratio_bkg_sigma_map->SetTitle("");

    gStyle->SetNumberContours(100);
    SetDynamicRatioPaletteRange(ratio_bkg_sigma_map); // Apply dynamic range

    gStyle->SetPaintTextFormat("4.3f");

    ratio_bkg_sigma_map->Draw("COLZ TEXT");
    gPad->Update();

    TPaletteAxis* palette_bkg_sigma_ratio = (TPaletteAxis*)ratio_bkg_sigma_map->GetListOfFunctions()->FindObject("palette");
    if (palette_bkg_sigma_ratio) {
        palette_bkg_sigma_ratio->SetX1NDC(0.86);
        palette_bkg_sigma_ratio->SetX2NDC(0.89);
        palette_bkg_sigma_ratio->SetY1NDC(0.12);
        palette_bkg_sigma_ratio->SetY2NDC(0.88);
        palette_bkg_sigma_ratio->SetLabelSize(0.03);
        palette_bkg_sigma_ratio->SetTitleOffset(1.3);
        
        TLatex palette_text;
        palette_text.SetNDC();
        palette_text.SetTextSize(0.035);
        palette_text.SetTextAlign(22);
        palette_text.SetTextAngle(90);
        palette_text.DrawLatex(0.97, 0.5, "#sigma_{MJB}^{Background, Post_bJR10V00} / #sigma_{MJB}^{Background, Pre_bJR10V00}");
    }

    DrawAtlasLabel();
    ratio_bkg_sigma_canvas->SaveAs(Form("../Plots/MJB_Ratio/%s/PlotsRdbpTbin/Summary/Bkg_R_MJB_sigma_ratio_Post_Pre_bJR10v00.pdf", fileWP.c_str()));
    std::cout << "Background R_MJB resolution ratio plot generated successfully." << std::endl;

    // --- Create and save the 2D plot for the Combined MC Resolution Ratio ---
    TCanvas* ratio_combined_sigma_canvas = new TCanvas("ratio_combined_sigma_canvas", Form("Combined MC R_{MJB} Resolution Ratio (Post_bJR10V00/Pre_bJR10V00) - WP: %s", labelWP_display.c_str()), 800, 600);
    ratio_combined_sigma_canvas->SetTopMargin(0.12);
    ratio_combined_sigma_canvas->SetRightMargin(0.15);
    ratio_combined_sigma_canvas->SetLeftMargin(0.12);
    ratio_combined_sigma_canvas->SetBottomMargin(0.12);
    ratio_combined_sigma_canvas->cd();

    TH2F* ratio_combined_sigma_map = new TH2F("ratio_combined_sigma_map", ";Mass Bin;p_{T} Bin", nMassCuts, -0.5, nMassCuts-0.5, nPtBins, -0.5, nPtBins-0.5);

    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < nPtBins; bin++) {
            ratio_combined_sigma_map->SetBinContent(m+1, bin+1, ratio_combined_sigmas[m][bin]);
            ratio_combined_sigma_map->SetBinError(m+1, bin+1, ratio_combined_sigma_errs[m][bin]);
        }
    }

    ratio_combined_sigma_map->SetStats(0);
    ratio_combined_sigma_map->SetMarkerSize(1.2);

    for (int m = 0; m < nMassCuts; m++) {
        ratio_combined_sigma_map->GetXaxis()->SetBinLabel(m+1, massCutLabels[m]);
    }
    for (int bin = 0; bin < nPtBins; bin++) {
        TString simplifiedLabel = Form("%d-%d GeV", (int)ptBins[bin], (int)ptBins[bin+1]);
        ratio_combined_sigma_map->GetYaxis()->SetBinLabel(bin+1, simplifiedLabel.Data());
    }

    ratio_combined_sigma_map->GetXaxis()->SetLabelSize(0.025);
    ratio_combined_sigma_map->GetYaxis()->SetLabelSize(0.025);
    ratio_combined_sigma_map->GetXaxis()->SetTitleSize(0.035);
    ratio_combined_sigma_map->GetYaxis()->SetTitleSize(0.035);
    ratio_combined_sigma_map->GetZaxis()->SetTitleSize(0.03);
    ratio_combined_sigma_map->SetTitle("");

    gStyle->SetNumberContours(100);
    SetDynamicRatioPaletteRange(ratio_combined_sigma_map); // Apply dynamic range

    gStyle->SetPaintTextFormat("4.3f");

    ratio_combined_sigma_map->Draw("COLZ TEXT");
    gPad->Update();

    TPaletteAxis* palette_combined_sigma_ratio = (TPaletteAxis*)ratio_combined_sigma_map->GetListOfFunctions()->FindObject("palette");
    if (palette_combined_sigma_ratio) {
        palette_combined_sigma_ratio->SetX1NDC(0.86);
        palette_combined_sigma_ratio->SetX2NDC(0.89);
        palette_combined_sigma_ratio->SetY1NDC(0.12);
        palette_combined_sigma_ratio->SetY2NDC(0.88);
        palette_combined_sigma_ratio->SetLabelSize(0.03);
        palette_combined_sigma_ratio->SetTitleOffset(1.3);
        
        TLatex palette_text;
        palette_text.SetNDC();
        palette_text.SetTextSize(0.035);
        palette_text.SetTextAlign(22);
        palette_text.SetTextAngle(90);
        palette_text.DrawLatex(0.97, 0.5, "#sigma_{MJB}^{Combined MC, Post_bJR10V00} / #sigma_{MJB}^{Combined MC, Pre_bJR10V00}");
    }

    DrawAtlasLabel();
    ratio_combined_sigma_canvas->SaveAs(Form("../Plots/MJB_Ratio/%s/PlotsRdbpTbin/Summary/Combined_MC_R_MJB_sigma_ratio_Post_Pre_bJR10v00.pdf", fileWP.c_str()));
    std::cout << "Combined MC R_MJB resolution ratio plot generated successfully." << std::endl;


    // --- Create and save the 2D plot for the (Data/Combined MC) Resolution Ratio (Post_bJR10V00/Pre_bJR10V00) ---
    TCanvas* final_ratio_sigma_canvas = new TCanvas("final_ratio_sigma_canvas", Form("Data/(Combined MC) Resolution Ratio (Post_bJR10V00/Pre_bJR10V00) - WP: %s", labelWP_display.c_str()), 800, 600);
    final_ratio_sigma_canvas->SetTopMargin(0.12);
    final_ratio_sigma_canvas->SetRightMargin(0.15);
    final_ratio_sigma_canvas->SetLeftMargin(0.12);
    final_ratio_sigma_canvas->SetBottomMargin(0.12);
    final_ratio_sigma_canvas->cd();

    TH2F* final_ratio_sigma_map = new TH2F("final_ratio_sigma_map", ";Mass Bin;p_{T} Bin", nMassCuts, -0.5, nMassCuts-0.5, nPtBins, -0.5, nPtBins-0.5);

    for (int m = 0; m < nMassCuts; m++) {
        for (int bin = 0; bin < nPtBins; bin++) {
            final_ratio_sigma_map->SetBinContent(m+1, bin+1, final_ratio_sigmas[m][bin]);
            final_ratio_sigma_map->SetBinError(m+1, bin+1, final_ratio_sigma_errs[m][bin]);
        }
    }

    final_ratio_sigma_map->SetStats(0);
    final_ratio_sigma_map->SetMarkerSize(1.2);

    for (int m = 0; m < nMassCuts; m++) {
        final_ratio_sigma_map->GetXaxis()->SetBinLabel(m+1, massCutLabels[m]);
    }
    for (int bin = 0; bin < nPtBins; bin++) {
        TString simplifiedLabel = Form("%d-%d GeV", (int)ptBins[bin], (int)ptBins[bin+1]);
        final_ratio_sigma_map->GetYaxis()->SetBinLabel(bin+1, simplifiedLabel.Data());
    }

    final_ratio_sigma_map->GetXaxis()->SetLabelSize(0.025);
    final_ratio_sigma_map->GetYaxis()->SetLabelSize(0.025);
    final_ratio_sigma_map->GetXaxis()->SetTitleSize(0.035);
    final_ratio_sigma_map->GetYaxis()->SetTitleSize(0.035);
    final_ratio_sigma_map->GetZaxis()->SetTitleSize(0.03);
    final_ratio_sigma_map->SetTitle("");

    gStyle->SetNumberContours(100);
    SetDynamicRatioPaletteRange(final_ratio_sigma_map); // Apply dynamic range

    gStyle->SetPaintTextFormat("4.3f");

    final_ratio_sigma_map->Draw("COLZ TEXT");
    gPad->Update();

    TPaletteAxis* palette_final_ratio_sigma = (TPaletteAxis*)final_ratio_sigma_map->GetListOfFunctions()->FindObject("palette");
    if (palette_final_ratio_sigma) {
        palette_final_ratio_sigma->SetX1NDC(0.86);
        palette_final_ratio_sigma->SetX2NDC(0.89);
        palette_final_ratio_sigma->SetY1NDC(0.12);
        palette_final_ratio_sigma->SetY2NDC(0.88);
        palette_final_ratio_sigma->SetLabelSize(0.03);
        palette_final_ratio_sigma->SetTitleOffset(1.3);
        
        TLatex palette_text;
        palette_text.SetNDC();
        palette_text.SetTextSize(0.035);
        palette_text.SetTextAlign(22);
        palette_text.SetTextAngle(90);
        palette_text.DrawLatex(0.97, 0.5, "(JER_{bJR00v10}^{Data}/JER_{bJR00v10}^{MC})/(JER^{Data}/JER^{MC})");
    } 

    DrawAtlasLabel();
    final_ratio_sigma_canvas->SaveAs(Form("../Plots/MJB_Ratio/%s/PlotsRdbpTbin/Summary/Data_Over_Combined_MC_Ratio_Post_Pre_bJR10v00_Sigmas.pdf", fileWP.c_str()));
    std::cout << "Data/(Combined MC) resolution ratio (Post_bJR10V00/Pre_bJR10V00) plot generated successfully." << std::endl;

    // Generate summary table
    std::string summaryFilePath = Form("../Plots/MJB_Ratio/%s/PlotsRdbpTbin/Summary/Summary_Ratios_%s.txt", fileWP.c_str(), fileWP.c_str());
    GenerateSummaryTable(
        summaryFilePath,
        labelWP_display,
        ratio_data_means, ratio_data_mean_errs,
        ratio_mc_means, ratio_mc_mean_errs,
        ratio_bkg_means, ratio_bkg_mean_errs,
        ratio_combined_means, ratio_combined_mean_errs,
        ratio_data_sigmas, ratio_data_sigma_errs,
        ratio_mc_sigmas, ratio_mc_sigma_errs,
        ratio_bkg_sigmas, ratio_bkg_sigma_errs,
        ratio_combined_sigmas, ratio_combined_sigma_errs,
        final_ratio_means, final_ratio_mean_errs,
        final_ratio_sigmas, final_ratio_sigma_errs
    );

    // Create summary plots showing ratios vs pT with mass bins as separate points
    CreateSummaryPlots(fileWP, labelWP_display, final_ratio_means, final_ratio_mean_errs,
                    final_ratio_sigmas, final_ratio_sigma_errs);
    // Clean up canvases
    delete ratio_data_canvas;
    delete ratio_data_map;
    delete ratio_mc_canvas;
    delete ratio_mc_map;
    delete ratio_bkg_canvas;
    delete ratio_bkg_map;
    delete ratio_combined_canvas;
    delete ratio_combined_map;
    delete final_ratio_mean_canvas;
    delete final_ratio_mean_map;
    delete ratio_data_sigma_canvas;
    delete ratio_data_sigma_map;
    delete ratio_mc_sigma_canvas;
    delete ratio_mc_sigma_map;
    delete ratio_bkg_sigma_canvas;
    delete ratio_bkg_sigma_map;
    delete ratio_combined_sigma_canvas;
    delete ratio_combined_sigma_map;
    delete final_ratio_sigma_canvas;
    delete final_ratio_sigma_map;
}


// Main function to run the analysis for different working points
void ratio_Response() { // Main entry point for ROOT
    // Disable ROOT global graphics for batch processing
    gROOT->SetBatch(kTRUE);

    // Load ATLAS style settings for consistent plotting
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();

    // Define working points to process
    std::vector<std::string> workingPoints = {"25", "46", "74"}; // Example working points

    // Process each working point
    for (const auto& wp : workingPoints) {
        process_single_wp_ratio(wp); // Call the renamed function to process each WP
    }
}


