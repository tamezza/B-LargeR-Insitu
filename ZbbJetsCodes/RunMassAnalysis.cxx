#include "../../Util/AtlasStyle.C" // Include ATLAS style macro

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
 * Performs the mass isolation and signal extraction analysis for a given working point.
 * @param wp The working point string (e.g., "25", "46", "74").
 */
void MassIsolationExtractSig_MJB(const std::string& wp) {
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

    // Create output directory if it doesn't exist for the current working point
    std::string outputDirPath = "../Plots/MJB/MassIsolation/" + fileWP + "/";
    gSystem->Exec(("mkdir -p " + outputDirPath).c_str());
    
    // Construct input file paths dynamically
    std::string mcPath = "../../NtupleSlim/MJB_" + fileWP + "wp_slim/ZbbJets_" + fileWP + "wp_slim.root";
    std::string dataPath = "../../NtupleSlim/MJB_" + fileWP + "wp_slim/data_" + fileWP + "wp_slim.root";
    std::string multijetsPath = "../../NtupleSlim/MJB_" + fileWP + "wp_slim/Multijets_" + fileWP + "wp_slim.root";

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

    // Create a temporary directory for this working point's objects
    // This isolates objects created in this iteration, preventing name collisions
    // and simplifying cleanup.
    TDirectory* currentDir = gDirectory; // Store current directory
    gDirectory->mkdir(("WP_" + fileWP).c_str())->cd(); // Create and cd into new directory

    // Create histograms for large-R jet mass
    TH1F *h_mc_ljet_m = new TH1F("h_mc_ljet_m", "Large-R Jet Mass", 35, 50, 150);
    TH1F *h_data_ljet_m = new TH1F("h_data_ljet_m", "Large-R Jet Mass", 35, 50, 150);
    TH1F *h_multijets_ljet_m = new TH1F("h_multijets_ljet_m", "Large-R Jet Mass", 35, 50, 150);

    // Enable Sumw2 for data histogram to ensure proper error calculation
    //h_data_ljet_m->Sumw2();
    // Enable Sumw2 for MC histogram to ensure proper error calculation for MC signal band
    //h_mc_ljet_m->Sumw2();

    // Set histogram styling
    h_data_ljet_m->SetMarkerStyle(20);
    h_data_ljet_m->SetMarkerSize(0.8);
    h_data_ljet_m->SetMarkerColor(kBlack);
    h_data_ljet_m->SetLineColor(kBlack);

    h_multijets_ljet_m->SetFillColor(kBlue);
    h_multijets_ljet_m->SetLineColor(kBlue);

    h_mc_ljet_m->SetFillColor(kRed);
    h_mc_ljet_m->SetLineColor(kRed);

    // Set up variables for reading from trees
    std::vector<float>* ljet_m = nullptr;
    double weight, weight_multijets;

    // Set branch addresses
    tree_mc->SetBranchAddress("ljet_m", &ljet_m);
    tree_data->SetBranchAddress("ljet_m", &ljet_m);
    tree_multijets->SetBranchAddress("ljet_m", &ljet_m);
    tree_mc->SetBranchAddress("total_weight", &weight);
    tree_multijets->SetBranchAddress("total_weight", &weight_multijets);

    // Fill histograms for MC signal
    Long64_t nentries = tree_mc->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree_mc->GetEntry(i);
        h_mc_ljet_m->Fill((*ljet_m)[0]/1000.0, weight); // Convert MeV to GeV
    }

    // Fill histograms for data
    nentries = tree_data->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree_data->GetEntry(i);
        h_data_ljet_m->Fill((*ljet_m)[0]/1000.0); // Convert MeV to GeV
    }

    // Fill histograms for multijets background
    nentries = tree_multijets->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree_multijets->GetEntry(i);
        h_multijets_ljet_m->Fill((*ljet_m)[0]/1000.0, weight_multijets); // Convert MeV to GeV
    }

    // Define sideband regions
    int binLow1 = h_multijets_ljet_m->FindBin(50);
    int binHigh1 = h_multijets_ljet_m->FindBin(75);
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
    
    printf("Linear fit results for WP %s: SF(m) = %.4f (± %.4f) + %.8f (± %.8f) * m\n", 
           wp.c_str(), p0, p0_err, p1, p1_err);
    
    // Save the fit result to a canvas
    TCanvas* fit_canvas = new TCanvas("fit_canvas", ("Scale Factor Fit - WP " + labelWP_display).c_str(), 800, 600);
    ratio_graph->SetTitle(("Data/Multijets Ratio in Sidebands (WP: " + labelWP_display + ");Large-R Jet Mass [GeV];Data/MC").c_str());
    ratio_graph->SetMarkerStyle(20);
    ratio_graph->SetMarkerColor(kBlue);
    ratio_graph->SetLineColor(kBlue);
    ratio_graph->Draw("AP");
    fit_func->SetLineColor(kRed);
    fit_func->Draw("same");
    
    // Add fit parameters to the canvas
    TLatex latex_fit; // Renamed to avoid conflict
    latex_fit.SetNDC();
    latex_fit.SetTextFont(42);
    latex_fit.SetTextSize(0.035);
    latex_fit.DrawLatex(0.2, 0.81, Form("SF(m) = %.4f + %.8f * m", p0, p1));
    latex_fit.DrawLatex(0.2, 0.78, Form("#chi^{2}/ndf = %.2f/%d", fit_func->GetChisquare(), fit_func->GetNDF()));
    
    DrawAtlasLabel();
    fit_canvas->SaveAs((outputDirPath + "scale_factor_fit.pdf").c_str());
    
    // Apply scale factor to the multijets background histogram
    TH1F* h_multijets_ljet_m_scaled = (TH1F*)h_multijets_ljet_m->Clone("h_multijets_ljet_m_scaled");
    
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
    std::cout << "Number of signal events for WP " << labelWP_display << ": " << signalEvents << std::endl;
    
    // Create a polynomial-exponential fit for background modeling
    TF1 *polyExpFit = new TF1("polyExpFit", "[0]*exp([1]*x + [2]*x*x + [3]*x*x*x)", 50, 150);
    
    // Clone the data histogram for fitting the sidebands
    TH1F *h_data_sidebands = (TH1F*)h_data_ljet_m->Clone("h_data_sidebands");
    
    // Remove the bins in the signal region from the clone
    for (int i = binHigh1 + 1; i < binLow2; ++i) {
        h_data_sidebands->SetBinContent(i, 0);
        h_data_sidebands->SetBinError(i, 0);
    }
    
    // Set initial parameters for polynomial exponential fit
    polyExpFit->SetParameter(0, 80);    // amplitude
    polyExpFit->SetParameter(1, -0.05);  // first order term (negative for falling exp)
    polyExpFit->SetParameter(2, -0.0005); // second order term
    polyExpFit->SetParameter(3, 0.000001); // third order term
    
    // Fit the sideband regions of the data
    h_data_sidebands->Fit(polyExpFit, "R");
    
    // Get the chi-square value and number of degrees of freedom
    double chiSquare = polyExpFit->GetChisquare();
    int ndf = polyExpFit->GetNDF();
    //double chi2NDF = chiSquare / ndf;
    double chi2NDF = (ndf > 0) ? chiSquare / ndf : 0.0; // Avoid division by zero

    // Create a clone of the data histogram to store the signal
    TH1F *h_signal_ljet_m = (TH1F*)h_data_ljet_m->Clone("h_signal_ljet_m");
    h_signal_ljet_m->SetTitle("Isolated Signal (Data);M^{J} [GeV];Events"); // Set title for drawing
    
    // Subtract the fit from the data to isolate the signal
    
    for (int i = 1; i <= h_signal_ljet_m->GetNbinsX(); ++i) {
        double binCenter = h_signal_ljet_m->GetBinCenter(i);
        double binContent = h_signal_ljet_m->GetBinContent(i);
        double fitValue = polyExpFit->Eval(binCenter);
        h_signal_ljet_m->SetBinContent(i, binContent - fitValue);
        
        // Set negative values to zero
        if (h_signal_ljet_m->GetBinContent(i) < 0) {
            h_signal_ljet_m->SetBinContent(i, 0);
        }
    }
    

    /* Maybe this: -----------------*/
    /*
    // Create a histogram representing the background fit for subtraction
    TH1F *h_background_fit = (TH1F*)h_data_ljet_m->Clone("h_background_fit_for_subtraction");
    h_background_fit->Reset(); // Clear contents
    h_background_fit->Sumw2(); // Enable Sumw2 for proper error handling during subtraction

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
    
    ------------------- */
    
    // Create stack plot with data overlay
    TCanvas* stack_canvas = new TCanvas("stack_canvas", ("Mass Distribution with Scale Factor - WP " + labelWP_display).c_str(), 800, 800);
    stack_canvas->Divide(1, 2);
    
    TPad* pad0 = new TPad("pad0", "pad0", 0, 0.3, 1, 1);
    pad0->SetBottomMargin(0.015);
    pad0->Draw();
    pad0->cd();
    
    THStack* hs = new THStack("hs", ("Large-R Jet Mass (WP: " + labelWP_display + ");M^{J} [GeV];Events").c_str());
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
    
    // MC statistical uncertainty band
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
    legend->AddEntry(h_mc_ljet_m, "Z (#rightarrow bb) + jets", "f");
    legend->AddEntry(h_multijets_ljet_m_scaled, "Multijets", "f");
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
    
    stack_canvas->SaveAs((outputDirPath + "Mass_stack_scaled.pdf").c_str());
    
    // Visualize the sideband regions
    TCanvas* sideband_canvas = new TCanvas("sideband_canvas", ("Fit Visualization - WP " + labelWP_display).c_str(), 800, 600);
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
    legendSideband->AddEntry(polyExpFit, Form("Poly-Exp Fit (#chi^{2}/NDF = %.2f)", chi2NDF), "f");
    legendSideband->AddEntry(sidebandBox1, "Sideband Regions", "f");
    legendSideband->Draw();
    
    DrawAtlasLabel();
    sideband_canvas->SaveAs((outputDirPath + "SideBandFitDebug.pdf").c_str());


    /* ---------------------------- ISOLATED SIGNAL ---------------------------- */
    
    // Plot the isolated signal with MC signal
    TCanvas* signal_canvas = new TCanvas("signal_canvas", ("Isolated Signal with Gaussian Fits - WP " + labelWP_display).c_str(), 800, 600);
    signal_canvas->cd();
    
    h_mc_ljet_m->SetFillColor(kRed);
    h_mc_ljet_m->GetYaxis()->SetRangeUser(0, h_mc_ljet_m->GetMaximum()*1.85);
    h_mc_ljet_m->GetXaxis()->SetTitle("M^{J} [GeV]");
    h_mc_ljet_m->GetYaxis()->SetTitle("Events");
    h_mc_ljet_m->Draw("HIST");

    // Create and draw the statistical uncertainty band for the MC signal
    TGraphAsymmErrors* mc_signal_err_band = new TGraphAsymmErrors(h_mc_ljet_m);
    for (int i = 0; i < mc_signal_err_band->GetN(); ++i) {
        double content = h_mc_ljet_m->GetBinContent(i + 1);
        double error = h_mc_ljet_m->GetBinError(i + 1);
        mc_signal_err_band->SetPointEYhigh(i, error);
        mc_signal_err_band->SetPointEYlow(i, error);
    }
    /*
    mc_signal_err_band->SetFillColorAlpha(kRed-7, 0.5); // A lighter shade of red, semi-transparent
    mc_signal_err_band->SetFillStyle(1001); // Solid fill
    mc_signal_err_band->SetLineWidth(0); // No border line
    mc_signal_err_band->Draw("2 same"); // Draw as a band
    */

    
    
    h_signal_ljet_m->SetMarkerStyle(20);
    h_signal_ljet_m->SetMarkerColor(kBlack);
    h_signal_ljet_m->Draw("EP same");
    
    // Create and draw the statistical uncertainty band for the isolated data signal
    TGraphAsymmErrors* signal_err_band = new TGraphAsymmErrors(h_signal_ljet_m);
    for (int i = 0; i < signal_err_band->GetN(); ++i) {
        double content = h_signal_ljet_m->GetBinContent(i + 1);
        double error = h_signal_ljet_m->GetBinError(i + 1);
        signal_err_band->SetPointEYhigh(i, error);
        signal_err_band->SetPointEYlow(i, error);
    }
    /*
    signal_err_band->SetFillColorAlpha(kGray, 0.5); // Light grey, semi-transparent
    signal_err_band->SetFillStyle(1001); // Solid fill
    signal_err_band->SetLineWidth(0); // No border line
    signal_err_band->Draw("2 same"); // Draw as a band
    */

    // Fit the MC signal with a Gaussian in the Z mass region
    TF1 *mc_gauss_fit = new TF1("mc_gauss_fit", "gaus", 83, 105);
    
    // Set initial parameters for MC Gaussian fit
    mc_gauss_fit->SetParameter(0, h_mc_ljet_m->GetMaximum());  // Amplitude
    //mc_gauss_fit->SetParameter(1, 90);  // Mean - approximately the Z mass
    //mc_gauss_fit->SetParameter(2, 10);  // Sigma

    // Use the histogram's properties to get more robust initial guesses
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
    double mc_chi2_ndf = mc_chi2 / mc_ndf;
    
    // Fit the extracted data signal with a Gaussian
    TF1 *data_gauss_fit = new TF1("data_gauss_fit", "gaus", 83, 105);
    
    // Set initial parameters for Data Gaussian fit
    data_gauss_fit->SetParameter(0, h_signal_ljet_m->GetMaximum());  // Amplitude
    //data_gauss_fit->SetParameter(1, 90);  // Mean
    //data_gauss_fit->SetParameter(2, 10);  // Sigma

    // Use the histogram's properties to get more robust initial guesses
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
    double data_chi2_ndf = data_chi2 / data_ndf;
    
    // Style the fits
    mc_gauss_fit->SetLineColor(kCyan);
    mc_gauss_fit->SetLineWidth(2);
    
    data_gauss_fit->SetLineColor(kViolet-4);
    data_gauss_fit->SetLineWidth(2);
    
    mc_gauss_fit->Draw("same");
    data_gauss_fit->Draw("same");
    


    // Create a legend
    TLegend *legendSignal = new TLegend(0.45, 0.65, 0.9, 0.9); 
    legendSignal->SetBorderSize(0);
    legendSignal->SetTextSize(0.018); 
    legendSignal->AddEntry(h_signal_ljet_m, "Isolated Signal (Data)", "lep");
    //legendSignal->AddEntry(signal_err_band, "Isolated Signal (Stat. Unc.)", "f"); // Added entry for data error band
    legendSignal->AddEntry(h_mc_ljet_m, "MC Signal", "f");
    //legendSignal->AddEntry(mc_signal_err_band, "MC Signal (Stat. Unc.)", "f"); // Added entry for MC error band
    legendSignal->AddEntry(mc_gauss_fit, Form("MC Gauss Fit: #mu = %.2f #pm %.2f, #sigma = %.2f #pm %.2f, #chi^{2}/ndf = %.2f", 
                          mc_mean, mc_mean_err, mc_sigma, mc_sigma_err, mc_chi2_ndf), "l");
    legendSignal->AddEntry(data_gauss_fit, Form("Data Gauss Fit: #mu = %.2f #pm %.2f, #sigma = %.2f #pm %.2f, #chi^{2}/ndf = %.2f", 
                          data_mean, data_mean_err, data_sigma, data_sigma_err, data_chi2_ndf), "l");
    legendSignal->Draw();
    
    DrawAtlasLabel();
    signal_canvas->SaveAs((outputDirPath + "Isolated_Signal_GaussFits.pdf").c_str());
    
    // Create a canvas for comparing the Gaussian fits
    TCanvas* gauss_canvas = new TCanvas("gauss_canvas", ("Gaussian Fits Comparison - WP " + labelWP_display).c_str(), 800, 600);
    gauss_canvas->cd();
    
    // Create a frame for drawing just the fits
    TH1F *h_frame = new TH1F("h_frame", ("Gaussian Fits Comparison (WP: " + labelWP_display + ");M^{J} [GeV];Arbitrary Units").c_str(), 100, 50, 150);
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
    double mean_diff = data_mean - mc_mean;
    double mean_diff_err = sqrt(data_mean_err*data_mean_err + mc_mean_err*mc_mean_err);
    
    TLatex gauss_latex;
    gauss_latex.SetNDC();
    gauss_latex.SetTextFont(42);
    gauss_latex.SetTextSize(0.03);
    gauss_latex.DrawLatex(0.2, 0.8, Form("#Delta#mu = %.2f #pm %.2f GeV", mean_diff, mean_diff_err));
    
    DrawAtlasLabel();
    gauss_canvas->SaveAs((outputDirPath + "GaussFitsComparison.pdf").c_str());
    
    // Print a summary of the fit results
    cout << "Summary of fit results for WP " << labelWP_display << ":" << endl;
    cout << "======================" << endl;
    cout << "Linear scale factor: SF(m) = " << p0 << " + " << p1 << " * m" << endl;
    cout << "MC Gaussian fit: #mu = " << mc_mean << " #pm " << mc_mean_err << " GeV, #sigma = " << mc_sigma << " #pm " << mc_sigma_err << " GeV" << endl;
    cout << "Data Gaussian fit: #mu = " << data_mean << " #pm " << data_mean_err << " GeV, #sigma = " << data_sigma << " #pm " << data_sigma_err << " GeV" << endl;
    cout << "Mass shift (Data-MC): #Delta#mu = " << mean_diff << " #pm " << mean_diff_err << " GeV" << endl;
    cout << "Width ratio (Data/MC): #sigma_{data}/#sigma_{mc} = " << data_sigma/mc_sigma << " #pm " << sqrt(pow(data_sigma_err/mc_sigma, 2) + pow(data_sigma*mc_sigma_err/(mc_sigma*mc_sigma), 2)) << endl;
    cout << "Number of signal events: " << signalEvents << endl;
    
    // Create a summary text file
    ofstream summary((outputDirPath + "fit_summary.txt").c_str());
    if (summary.is_open()) {
        summary << "Large-R Jet Mass Analysis Summary for WP " << labelWP_display << endl;
        summary << "=================================" << endl << endl;
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
        
        summary.close();
        cout << "Summary saved to " << outputDirPath << "fit_summary.txt" << endl;
    }

    // Close and delete TFile objects
    if (fmc) { fmc->Close(); delete fmc; }
    if (fdata) { fdata->Close(); delete fdata; }
    if (fmultijets) { fmultijets->Close(); delete fmultijets; }

    // Go back to the original directory (e.g., root directory) before deleting the temporary directory
    gROOT->cd(); // Change to the ROOT global directory
    gDirectory->Delete(("WP_" + fileWP + ";*").c_str()); // Delete the temporary directory and its contents
}

// Main function to run the analysis
void RunMassAnalysis() {
    // Disable ROOT global graphics for batch processing
    gROOT->SetBatch(kTRUE);

    // Load ATLAS style settings for consistent plotting
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();

    // Define working points to process
    std::vector<std::string> workingPoints = {"25", "46", "74"}; // Example working points

    // Process each working point
    for (const auto& wp : workingPoints) {
        MassIsolationExtractSig_MJB(wp);
    }
    // No global cleanup calls here, as objects are managed by temporary directories.
}