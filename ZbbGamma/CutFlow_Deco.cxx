#include "../../Util/AtlasStyle.C"

void DrawAtlasLabel() {
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(72);
    latex.SetTextColor(kBlack);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.2, 0.88, "ATLAS");

    latex.SetTextFont(42);
    latex.DrawLatex(0.27, 0.88, "Internal");

    latex.DrawLatex(0.2, 0.83, "#sqrt{s} = 13 TeV, 140 fb^{-1}");
}

void CompareEfficiencies() {
    // Try to load AtlasStyle.C, but handle the case if it's not found
    bool atlasStyleLoaded = true;
    try {
        gROOT->LoadMacro("../../Util/AtlasStyle.C");
        SetAtlasStyle();
    } catch (...) {
        std::cerr << "Warning: AtlasStyle.C not found. Using default ROOT style." << std::endl;
        atlasStyleLoaded = false;
        // Set a reasonable default style
        gStyle->SetOptStat(0);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gStyle->SetLegendBorderSize(0);
    }
    
    // Define working points to compare
    //std::vector<std::string> workingPoints = {"25", "30", "37", "46", "58", "74", "94", "125", "155"};
    std::vector<std::string> workingPoints = {"25", "46", "74"};
    
    // Get username for EOS path
    std::string username = std::getenv("USER") ? std::getenv("USER") : "tamezza";
    
    // Define EOS path
    //std::string eosPath = "root://eosuser.cern.ch//eos/user/t/" + username + "/Documents/EasyBjets/ZbbJets/";
    std::string eosPath = "../../";

    // Define colors for different samples
    int dataColor = kBlack;
    int signalColor = kRed;
    int bkgColor = kBlue;
    
    // Define line styles
    int dataStyle = 1;    // Solid line
    int signalStyle = 2;  // Dashed line
    int bkgStyle = 9;     // Long dash line
    
    // Vectors to store efficiency histograms for each sample type
    std::vector<TH1F*> dataEfficiencyHistograms;
    std::vector<TH1F*> signalEfficiencyHistograms;
    std::vector<TH1F*> bkgEfficiencyHistograms;
    
    // Create efficiency histograms for each working point
    for (size_t i = 0; i < workingPoints.size(); ++i) {
        const auto& wp = workingPoints[i];

        std::string wpLabel;
        if (wp == "125") {
            wpLabel = "1p25";
        } else if (wp == "155") {
            wpLabel = "1p55";
        } else {
            wpLabel = "0p" + wp;
        }
        
        // Data histogram
        TH1F* histData = new TH1F(("h_eff_data_" + wpLabel).c_str(), 
                                  ("Data Efficiency (WP: " + wpLabel + ");Cut;Efficiency").c_str(), 
                                  8, 0, 8);
        histData->SetLineColor(dataColor);
        histData->SetLineWidth(3);
        histData->SetLineStyle(dataStyle);
        histData->SetMarkerStyle(20);
        histData->SetMarkerColor(dataColor);
        histData->SetMarkerSize(1.5);
        dataEfficiencyHistograms.push_back(histData);
        
        // Signal histogram
        TH1F* histSignal = new TH1F(("h_eff_signal_" + wpLabel).c_str(), 
                                    ("Signal Efficiency (WP: " + wpLabel + ");Cut;Efficiency").c_str(), 
                                    8, 0, 8);
        histSignal->SetLineColor(signalColor);
        histSignal->SetLineWidth(3);
        histSignal->SetLineStyle(signalStyle);
        histSignal->SetMarkerStyle(21);
        histSignal->SetMarkerColor(signalColor);
        histSignal->SetMarkerSize(1.3);
        signalEfficiencyHistograms.push_back(histSignal);
        
        // Background histogram
        TH1F* histBkg = new TH1F(("h_eff_bkg_" + wpLabel).c_str(), 
                                 ("Background Efficiency (WP: " + wpLabel + ");Cut;Efficiency").c_str(), 
                                 8, 0, 8);
        histBkg->SetLineColor(bkgColor);
        histBkg->SetLineWidth(3);
        histBkg->SetLineStyle(bkgStyle);
        histBkg->SetMarkerStyle(22);
        histBkg->SetMarkerColor(bkgColor);
        histBkg->SetMarkerSize(1.3);
        bkgEfficiencyHistograms.push_back(histBkg);
    }
    
    // Define cut names
    std::vector<std::string> cutNames = {
        "Initial Events",              // Cut 0
        "p_{T}^{#gamma}",           // Cut 1
        "#eta^{#gamma}",            // Cut 2
        "p_{T}^{J1}",               // Cut 3
        "#Delta R (#gamma, J1)",    // Cut 4: Gamma Overlap Removal
        "#Delta #Phi (#gamma , J1)",// Cut 5
        "(JVT^{J2} , p_{T}^{J2} , #eta^{J2})", //Cut6 Pileup Removal
        "#Delta R (#gamma , J2)",   // Cut 7
        "#Delta R (J1 , J2)",       // Cut 8
        "Veto",                     // Cut 9 Extra-radiation 
        "D_{xbb} Cuts",                // Cut 10
        "M^{J1}"                    // Cut 11
   };
    
    // Process each working point
    for (size_t i = 0; i < workingPoints.size(); ++i) {
        const auto& wp = workingPoints[i];

        std::string fileWP, labelWP;
        if (wp == "125") {
            fileWP = "1p25";
            labelWP = "1p25";
        } else if (wp == "155") {
            fileWP = "1p55";
            labelWP = "1p55";
        } else {
            fileWP = "0p" + wp;
            labelWP = "0p" + wp;
        }
        
        // Open all three files for this working point
        std::string dataPath = eosPath + "NtupleSlim/DB_Deco_" + fileWP + "wp_slim/data_deco_" + fileWP + "wp_slim.root";
        std::string signalPath = eosPath + "NtupleSlim/DB_Deco_" + fileWP + "wp_slim/Zbby_deco_" + fileWP + "wp_slim.root";
        std::string bkgPath = eosPath + "NtupleSlim/DB_Deco_" + fileWP + "wp_slim/gamma_jets_deco_" + fileWP + "wp_slim.root";
        
        TFile* dataFile = TFile::Open(dataPath.c_str(), "READ");
        TFile* signalFile = TFile::Open(signalPath.c_str(), "READ");
        TFile* bkgFile = TFile::Open(bkgPath.c_str(), "READ");
        
        if (!dataFile || dataFile->IsZombie()) {
            std::cerr << "Error: Could not open data file " << dataPath << std::endl;
            if (dataFile) delete dataFile;
            if (signalFile) delete signalFile;
            if (bkgFile) delete bkgFile;
            continue;
        }
        
        if (!signalFile || signalFile->IsZombie()) {
            std::cerr << "Error: Could not open signal file " << signalPath << std::endl;
            if (dataFile) delete dataFile;
            if (signalFile) delete signalFile;
            if (bkgFile) delete bkgFile;
            continue;
        }
        
        if (!bkgFile || bkgFile->IsZombie()) {
            std::cerr << "Error: Could not open background file " << bkgPath << std::endl;
            if (dataFile) delete dataFile;
            if (signalFile) delete signalFile;
            if (bkgFile) delete bkgFile;
            continue;
        }
        
        // Get the cutflow histograms
        TH1F* dataCutflow = static_cast<TH1F*>(dataFile->Get("h_cutflow"));
        TH1F* signalCutflow = static_cast<TH1F*>(signalFile->Get("h_cutflow"));
        TH1F* bkgCutflow = static_cast<TH1F*>(bkgFile->Get("h_cutflow"));
        
        if (!dataCutflow || !signalCutflow || !bkgCutflow) {
            std::cerr << "Error: Could not find cutflow histograms for WP " << labelWP << std::endl;
            dataFile->Close(); signalFile->Close(); bkgFile->Close();
            delete dataFile; delete signalFile; delete bkgFile;
            continue;
        }
        
        // Calculate efficiencies for data
        double dataInitialCount = dataCutflow->GetBinContent(1);
        for (int bin = 1; bin <= 12 && bin <= dataCutflow->GetNbinsX(); ++bin) {
            double count = dataCutflow->GetBinContent(bin);
            double efficiency = (dataInitialCount > 0) ? count / dataInitialCount : 0;
            dataEfficiencyHistograms[i]->SetBinContent(bin, efficiency);

            if (bin == 2) {
                dataEfficiencyHistograms[i]->GetXaxis()->SetBinLabel(bin, ("D_{xbb} " + labelWP + " WP").c_str());
            } else {
                dataEfficiencyHistograms[i]->GetXaxis()->SetBinLabel(bin, cutNames[bin-1].c_str());
            }
        }
        
        // Calculate efficiencies for signal
        double signalInitialCount = signalCutflow->GetBinContent(1);
        for (int bin = 1; bin <= 12 && bin <= signalCutflow->GetNbinsX(); ++bin) {
            double count = signalCutflow->GetBinContent(bin);
            double efficiency = (signalInitialCount > 0) ? count / signalInitialCount : 0;
            signalEfficiencyHistograms[i]->SetBinContent(bin, efficiency);

            if (bin == 2) {
                signalEfficiencyHistograms[i]->GetXaxis()->SetBinLabel(bin, ("D_{xbb} " + labelWP + " WP").c_str());
            } else {
                signalEfficiencyHistograms[i]->GetXaxis()->SetBinLabel(bin, cutNames[bin-1].c_str());
            }
        }
        
        // Calculate efficiencies for background
        double bkgInitialCount = bkgCutflow->GetBinContent(1);
        for (int bin = 1; bin <= 12 && bin <= bkgCutflow->GetNbinsX(); ++bin) {
            double count = bkgCutflow->GetBinContent(bin);
            double efficiency = (bkgInitialCount > 0) ? count / bkgInitialCount : 0;
            bkgEfficiencyHistograms[i]->SetBinContent(bin, efficiency);

            if (bin == 2) {
                bkgEfficiencyHistograms[i]->GetXaxis()->SetBinLabel(bin, ("D_{xbb} " + labelWP + " WP").c_str());
            } else {
                bkgEfficiencyHistograms[i]->GetXaxis()->SetBinLabel(bin, cutNames[bin-1].c_str());
            }
        }
        
        dataFile->Close(); signalFile->Close(); bkgFile->Close();
        delete dataFile; delete signalFile; delete bkgFile;
    }
    
    // Now create comparison plots for specific working points
    // Let's create a plot comparing data, signal, and background for a specific WP
    //int wpIndex = 4; // Use WP 58 (index 4) as an example - you can change this
    for (int wpIndex = 0; wpIndex < workingPoints.size(); ++wpIndex) {
        
        // Skip if histograms don't have entries
        if (dataEfficiencyHistograms[wpIndex]->GetEntries() == 0) continue;
        
        auto canvasComp = new TCanvas(("canvasComp_WP" + workingPoints[wpIndex]).c_str(), 
                                      "Efficiency Comparison - Data/Signal/Background", 800, 600);
        canvasComp->SetGridy();
        canvasComp->SetLogy();
        canvasComp->SetTopMargin(0.06);
        canvasComp->SetBottomMargin(0.12);
        canvasComp->SetLeftMargin(0.12);
        canvasComp->SetRightMargin(0.04);
        
        // Create legend
        auto legendComp = new TLegend(0.65, 0.65, 0.95, 0.90);
        legendComp->SetBorderSize(0);
        legendComp->SetFillStyle(0);
        legendComp->SetTextSize(0.03);
        
        // Get WP value for legend and title
        std::string wpValue;
        if (workingPoints[wpIndex] == "125") {
            wpValue = "1.25";
        } else if (workingPoints[wpIndex] == "155") {
            wpValue = "1.55";
        } else {
            wpValue = "0." + workingPoints[wpIndex];
        }
        
        // Set y-axis range
        dataEfficiencyHistograms[wpIndex]->GetYaxis()->SetRangeUser(1e-4, 3.0);
        dataEfficiencyHistograms[wpIndex]->GetYaxis()->SetTitle("Data Selection Efficiency");
        dataEfficiencyHistograms[wpIndex]->GetYaxis()->SetTitleOffset(1.3);
        dataEfficiencyHistograms[wpIndex]->GetXaxis()->SetTitle("Cuts");
        dataEfficiencyHistograms[wpIndex]->SetTitle(("Efficiency Comparison (" + wpValue + " % QCD Eff. WP)").c_str());
        
        // Draw histograms
        dataEfficiencyHistograms[wpIndex]->Draw("LP");
        signalEfficiencyHistograms[wpIndex]->Draw("LP SAME");
        bkgEfficiencyHistograms[wpIndex]->Draw("LP SAME");
        
        // Add to legend
        legendComp->AddEntry(dataEfficiencyHistograms[wpIndex], "Data", "lp");
        legendComp->AddEntry(signalEfficiencyHistograms[wpIndex], "Z+bb jets", "lp");
        legendComp->AddEntry(bkgEfficiencyHistograms[wpIndex], "#gamma + jets", "lp");
        
        legendComp->Draw();
        DrawAtlasLabel();
        
        // Add working point label on the plot
        TLatex wpLabel;
        wpLabel.SetNDC();
        wpLabel.SetTextFont(42);
        wpLabel.SetTextSize(0.04);
        wpLabel.DrawLatex(0.2, 0.85, (wpValue + " % QCD Eff. WP").c_str());
        
        // Save the comparison plot with WP-specific filename
        canvasComp->SaveAs(("../Plots/DB_Deco/CutFlow/EfficiencyComparison_DataSignalBkg_WP" + workingPoints[wpIndex] + ".pdf").c_str());
        //canvasComp->SaveAs(("../Plots/DB_Deco/CutFlow/EfficiencyComparison_DataSignalBkg_WP" + workingPoints[wpIndex] + ".png").c_str());
        
        // Clean up this specific canvas and legend
        delete canvasComp;
        delete legendComp;
    }
    
    // Print a summary table
    std::cout << "\n+----------+------------------+------------------+------------------+\n";
    std::cout << "| WP       | Data Final Eff.  | Signal Final Eff.| Bkg Final Eff.   |\n";
    std::cout << "+----------+------------------+------------------+------------------+\n";
    
    for (size_t i = 0; i < workingPoints.size(); ++i) {
        if (dataEfficiencyHistograms[i]->GetEntries() > 0) {
            double dataFinalEff = dataEfficiencyHistograms[i]->GetBinContent(8);
            double signalFinalEff = signalEfficiencyHistograms[i]->GetBinContent(8);
            double bkgFinalEff = bkgEfficiencyHistograms[i]->GetBinContent(8);
            
            std::string wpDisplay;
            if (workingPoints[i] == "125") {
                wpDisplay = "1p25";
            } else if (workingPoints[i] == "155") {
                wpDisplay = "1p55";
            } else {
                wpDisplay = "0p" + workingPoints[i];
            }
            
            std::cout << "| " << std::left << std::setw(8) << wpDisplay << " | " 
                     << std::right << std::setw(16) << std::scientific << std::setprecision(4) << dataFinalEff << " | " 
                     << std::setw(16) << signalFinalEff << " | "
                     << std::setw(16) << bkgFinalEff << " |\n";
        }
    }
    
    std::cout << "+----------+------------------+------------------+------------------+\n";
    
    // Optional: Create a single canvas with subplots for all working points
    auto canvasAll = new TCanvas("canvasAll", "All Working Points Comparison", 800, 600);
    canvasAll->Divide(3, 3); // 3x3 grid for 9 working points
    
    for (int wpIndex = 0; wpIndex < workingPoints.size() && wpIndex < 9; ++wpIndex) {
        if (dataEfficiencyHistograms[wpIndex]->GetEntries() == 0) continue;
        
        canvasAll->cd(wpIndex + 1);
        //gPad->SetGridy();
        gPad->SetLogy();
        gPad->SetTopMargin(0.08);
        gPad->SetBottomMargin(0.15);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        
        // Clone histograms for subplot
        TH1F* dataClone = (TH1F*)dataEfficiencyHistograms[wpIndex]->Clone(("data_clone_" + std::to_string(wpIndex)).c_str());
        TH1F* signalClone = (TH1F*)signalEfficiencyHistograms[wpIndex]->Clone(("signal_clone_" + std::to_string(wpIndex)).c_str());
        TH1F* bkgClone = (TH1F*)bkgEfficiencyHistograms[wpIndex]->Clone(("bkg_clone_" + std::to_string(wpIndex)).c_str());
        
        // Adjust text sizes for subplots
        dataClone->GetXaxis()->SetLabelSize(0.06);
        dataClone->GetYaxis()->SetLabelSize(0.06);
        dataClone->GetXaxis()->SetTitleSize(0.06);
        dataClone->GetYaxis()->SetTitleSize(0.06);
        dataClone->GetYaxis()->SetTitleOffset(1.0);
        
        // Get WP value for title
        std::string wpValue;
        if (workingPoints[wpIndex] == "125") {
            wpValue = "1.25 % WP";
        } else if (workingPoints[wpIndex] == "155") {
            wpValue = "1.55 % WP";
        } else {
            wpValue = "0." + workingPoints[wpIndex] + " % WP";
        }
        
        dataClone->SetTitle(wpValue.c_str());
        dataClone->Draw("LP");
        signalClone->Draw("LP SAME");
        bkgClone->Draw("LP SAME");
        
        // Add a small legend for the first subplot only
        if (wpIndex == 0) {
            auto legendSmall = new TLegend(0.55, 0.70, 0.95, 0.92);
            legendSmall->SetBorderSize(0);
            legendSmall->SetFillStyle(0);
            legendSmall->SetTextSize(0.05);
            legendSmall->AddEntry(dataClone, "Data", "lp");
            legendSmall->AddEntry(signalClone, "Z(#rightarrowbb)+ #gamma", "lp");
            legendSmall->AddEntry(bkgClone, "#gamma + jets", "lp");
            legendSmall->Draw();
        }
    }
    
    canvasAll->SaveAs("../Plots/DB_Deco/CutFlow/EfficiencyComparison_AllWPs_Grid.pdf");
    //canvasAll->SaveAs("../Plots/DB_Deco/CutFlow/EfficiencyComparison_AllWPs_Grid.png");
    
    delete canvasAll;
    // Clean up
    for (auto hist : dataEfficiencyHistograms) delete hist;
    for (auto hist : signalEfficiencyHistograms) delete hist;
    for (auto hist : bkgEfficiencyHistograms) delete hist;
    //delete canvasComp;
    //delete legendComp;
    //delete canvasSignal;
    //delete legendSignal;
}

void CompareEfficienciesData() {
    bool atlasStyleLoaded = true;
    try {
        gROOT->LoadMacro("../../Util/AtlasStyle.C");
        SetAtlasStyle();
    } catch (...) {
        std::cerr << "Warning: AtlasStyle.C not found. Using default ROOT style." << std::endl;
        atlasStyleLoaded = false;
        // Set a reasonable default style
        gStyle->SetOptStat(0);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gStyle->SetLegendBorderSize(0);
    }
    
    // Define working points to compare
    //std::vector<std::string> workingPoints = {"25", "30", "37", "46", "58", "74", "94", "125", "155"};
    std::vector<std::string> workingPoints = {"25", "46", "74"};
    
    // Get username for EOS path
    std::string username = std::getenv("USER") ? std::getenv("USER") : "tamezza";
    
    // Define EOS path
    //std::string eosPath = "root://eosuser.cern.ch//eos/user/t/" + username + "/Documents/EasyBjets/ZbbJets/";
    std::string eosPath = "../../";
    
    
    // Define colors for different working points
    std::vector<int> colors = {kRed, kBlue, kGreen+2, kMagenta, kCyan+2, kOrange+7, kViolet+3, kPink+2, kSpring+4};
    
    // Vectors to store efficiency histograms for data, signal, and background
    std::vector<TH1F*> efficiencyHistograms;
    std::vector<TH1F*> signalEfficiencyHistograms;
    std::vector<TH1F*> bkgEfficiencyHistograms;
    
    // Create efficiency histograms for all samples
    for (size_t i = 0; i < workingPoints.size(); ++i) {
        const auto& wp = workingPoints[i];

        std::string wpLabel;
        if (wp == "125") {
            wpLabel = "1p25";
        } else if (wp == "155") {
            wpLabel = "1p55";
        } else {
            wpLabel = "0p" + wp;
        }
        
        // Data histogram
        TH1F* hist = new TH1F(("h_eff_" + wpLabel).c_str(), 
                              ("Efficiency (WP: " + wpLabel + ");Cut;Efficiency").c_str(), 
                              8, 0, 8);
        hist->SetLineColor(colors[i]);
        hist->SetLineWidth(2);
        hist->SetMarkerStyle(20 + i);
        hist->SetMarkerColor(colors[i]);
        efficiencyHistograms.push_back(hist);
        
        // Signal histogram
        TH1F* histSignal = new TH1F(("h_eff_signal_" + wpLabel).c_str(), 
                                    ("Signal Efficiency (WP: " + wpLabel + ");Cut;Efficiency").c_str(), 
                                    8, 0, 8);
        histSignal->SetLineColor(colors[i]);
        histSignal->SetLineWidth(2);
        histSignal->SetMarkerStyle(20 + i);
        histSignal->SetMarkerColor(colors[i]);
        signalEfficiencyHistograms.push_back(histSignal);
        
        // Background histogram
        TH1F* histBkg = new TH1F(("h_eff_bkg_" + wpLabel).c_str(), 
                                 ("Background Efficiency (WP: " + wpLabel + ");Cut;Efficiency").c_str(), 
                                 8, 0, 8);
        histBkg->SetLineColor(colors[i]);
        histBkg->SetLineWidth(2);
        histBkg->SetMarkerStyle(20 + i);
        histBkg->SetMarkerColor(colors[i]);
        bkgEfficiencyHistograms.push_back(histBkg);
    }
    
    // Define cut names
    std::vector<std::string> cutNames = {
        "Initial Events",          // Cut 0
        "D_{xbb} WP",              // Cut 1
        "p_{T}^{J} > 450",         // Cut 2
        "|#eta^{J}| < 1.2",        // Cut 3
        "Small-R p_{T} & #eta",    // Cut 4
        "JVT > 0.59",              // Cut 5
        "|#alpha - #pi| < 0.3",    // Cut 6
        "#Delta R (J1, DB) > 1.5" // Cut 7
    };
    
    // Process each working point for all samples
    for (size_t i = 0; i < workingPoints.size(); ++i) {
        const auto& wp = workingPoints[i];

        std::string fileWP, labelWP;
        if (wp == "125") {
            fileWP = "1p25";
            labelWP = "1p25";
        } else if (wp == "155") {
            fileWP = "1p55";
            labelWP = "1p55";
        } else {
            fileWP = "0p" + wp;
            labelWP = "0p" + wp;
        }
        
        // Process DATA file
        std::string filePath = eosPath + "NtupleSlim/DB_Deco_" + fileWP + "wp_slim/data_deco_" + fileWP + "wp_slim.root";
        TFile* file = TFile::Open(filePath.c_str(), "READ");
        
        if (!file || file->IsZombie()) {
            std::cerr << "Error: Could not open file " << filePath << std::endl;
            if (file) delete file;
            continue;
        }
        
        TH1F* cutflowHist = static_cast<TH1F*>(file->Get("h_cutflow"));
        if (!cutflowHist) {
            std::cerr << "Error: Could not find h_cutflow in file " << filePath << std::endl;
            file->Close();
            delete file;
            continue;
        }
        
        double initialCount = cutflowHist->GetBinContent(1);
        for (int bin = 1; bin <= 12 && bin <= cutflowHist->GetNbinsX(); ++bin) {
            double count = cutflowHist->GetBinContent(bin);
            double efficiency = (initialCount > 0) ? count / initialCount : 0;
            efficiencyHistograms[i]->SetBinContent(bin, efficiency);

            if (bin == 2) {
                efficiencyHistograms[i]->GetXaxis()->SetBinLabel(bin, ("D_{xbb} " + labelWP + " WP").c_str());
            } else {
                efficiencyHistograms[i]->GetXaxis()->SetBinLabel(bin, cutNames[bin-1].c_str());
            }
        }
        
        file->Close();
        delete file;
        
        // Process SIGNAL file
        std::string signalPath = eosPath + "NtupleSlim/DB_Deco_" + fileWP + "wp_slim/Zbby_deco_" + fileWP + "wp_slim.root";
        TFile* signalFile = TFile::Open(signalPath.c_str(), "READ");
        
        if (!signalFile || signalFile->IsZombie()) {
            std::cerr << "Error: Could not open signal file " << signalPath << std::endl;
            if (signalFile) delete signalFile;
            continue;
        }
        
        TH1F* signalCutflow = static_cast<TH1F*>(signalFile->Get("h_cutflow"));
        if (!signalCutflow) {
            std::cerr << "Error: Could not find h_cutflow in signal file " << signalPath << std::endl;
            signalFile->Close();
            delete signalFile;
            continue;
        }
        
        double signalInitialCount = signalCutflow->GetBinContent(1);
        for (int bin = 1; bin <= 12 && bin <= signalCutflow->GetNbinsX(); ++bin) {
            double count = signalCutflow->GetBinContent(bin);
            double efficiency = (signalInitialCount > 0) ? count / signalInitialCount : 0;
            signalEfficiencyHistograms[i]->SetBinContent(bin, efficiency);

            if (bin == 2) {
                signalEfficiencyHistograms[i]->GetXaxis()->SetBinLabel(bin, ("D_{xbb} " + labelWP + " WP").c_str());
            } else {
                signalEfficiencyHistograms[i]->GetXaxis()->SetBinLabel(bin, cutNames[bin-1].c_str());
            }
        }
        
        signalFile->Close();
        delete signalFile;
        
        // Process BACKGROUND file
        std::string bkgPath = eosPath + "NtupleSlim/DB_Deco_" + fileWP + "wp_slim/gamma_jets_deco_" + fileWP + "wp_slim.root";
        TFile* bkgFile = TFile::Open(bkgPath.c_str(), "READ");
        
        if (!bkgFile || bkgFile->IsZombie()) {
            std::cerr << "Error: Could not open background file " << bkgPath << std::endl;
            if (bkgFile) delete bkgFile;
            continue;
        }
        
        TH1F* bkgCutflow = static_cast<TH1F*>(bkgFile->Get("h_cutflow"));
        if (!bkgCutflow) {
            std::cerr << "Error: Could not find h_cutflow in background file " << bkgPath << std::endl;
            bkgFile->Close();
            delete bkgFile;
            continue;
        }
        
        double bkgInitialCount = bkgCutflow->GetBinContent(1);
        for (int bin = 1; bin <= 12 && bin <= bkgCutflow->GetNbinsX(); ++bin) {
            double count = bkgCutflow->GetBinContent(bin);
            double efficiency = (bkgInitialCount > 0) ? count / bkgInitialCount : 0;
            bkgEfficiencyHistograms[i]->SetBinContent(bin, efficiency);

            if (bin == 2) {
                bkgEfficiencyHistograms[i]->GetXaxis()->SetBinLabel(bin, ("D_{xbb} " + labelWP + " WP").c_str());
            } else {
                bkgEfficiencyHistograms[i]->GetXaxis()->SetBinLabel(bin, cutNames[bin-1].c_str());
            }
        }
        
        bkgFile->Close();
        delete bkgFile;
    }
    
    // ==================== DATA EFFICIENCY CANVAS (LOG SCALE) ====================
    auto canvasEff = new TCanvas("canvasEff", "Data Efficiency Comparison", 800, 600);
    canvasEff->SetLogy();
    
    auto legend = new TLegend(0.60, 0.55, 0.9, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.030);
    
    bool firstDrawn = false;
    double minEff = 1.0, maxEff = 0.0;
    
    for (size_t i = 0; i < efficiencyHistograms.size(); ++i) {
        if (efficiencyHistograms[i]->GetEntries() > 0) {
            for (int bin = 1; bin <= efficiencyHistograms[i]->GetNbinsX(); ++bin) {
                double eff = efficiencyHistograms[i]->GetBinContent(bin);
                if (eff > 0) {
                    if (eff < minEff) minEff = eff;
                    if (eff > maxEff) maxEff = eff;
                }
            }
        }
    }
    
    double yMin = (minEff > 0) ? minEff * 0.5 : 1e-7;
    double yMax = (maxEff > 0) ? maxEff * 3.0 : 3.0;
    if (yMin < 1e-7) yMin = 1e-7;
    if (yMax > 3.0) yMax = 3.0;
    
    for (size_t i = 0; i < efficiencyHistograms.size(); ++i) {
        if (efficiencyHistograms[i]->GetEntries() > 0) {
            efficiencyHistograms[i]->GetYaxis()->SetRangeUser(yMin, yMax);
            efficiencyHistograms[i]->GetYaxis()->SetTitle("Data Selection Efficiency");
            efficiencyHistograms[i]->GetYaxis()->SetTitleOffset(1.3);
            efficiencyHistograms[i]->GetYaxis()->SetTitleSize(0.045);
            efficiencyHistograms[i]->GetYaxis()->SetLabelSize(0.04);
            
            efficiencyHistograms[i]->GetXaxis()->SetTitle("Cuts");
            efficiencyHistograms[i]->GetXaxis()->SetTitleSize(0.045);
            efficiencyHistograms[i]->GetXaxis()->SetLabelSize(0.04);
            efficiencyHistograms[i]->GetXaxis()->SetTitleOffset(1.1);
            
            efficiencyHistograms[i]->SetLineWidth(3);
            efficiencyHistograms[i]->SetMarkerSize(1.5);
            
            efficiencyHistograms[i]->Draw(firstDrawn ? "LP SAME" : "LP");
            
            std::string legendLabel;
            if (workingPoints[i] == "125") {
                legendLabel = "1.25 % QCD Eff. WP";
            } else if (workingPoints[i] == "155") {
                legendLabel = "1.55 % QCD Eff. WP";
            } else {
                double wpValue = std::stod(workingPoints[i]) / 100.0;
                std::stringstream ss;
                ss << std::fixed << std::setprecision(2) << wpValue << " % QCD Eff. WP";
                legendLabel = ss.str();
            }
            
            legend->AddEntry(efficiencyHistograms[i], legendLabel.c_str(), "lp");
            firstDrawn = true;
        }
    }
    
    legend->Draw();
    DrawAtlasLabel();
    
    canvasEff->SaveAs("../Plots/DB_Deco/CutFlow/DataEfficiencyComparison_log.pdf");
    //canvasEff->SaveAs("../Plots/DB_Deco/CutFlow/DataEfficiencyComparison_log.png");
    
    // ==================== SIGNAL EFFICIENCY CANVAS (LOG SCALE) ====================
    auto canvasSignal = new TCanvas("canvasSignal", "Signal Efficiency Comparison", 800, 600);
    canvasSignal->SetLogy();
    
    auto legendSignal = new TLegend(0.60, 0.55, 0.9, 0.9);
    legendSignal->SetBorderSize(0);
    legendSignal->SetFillStyle(0);
    legendSignal->SetTextSize(0.030);
    
    double minEffSignal = 1.0, maxEffSignal = 0.0;
    for (size_t i = 0; i < signalEfficiencyHistograms.size(); ++i) {
        if (signalEfficiencyHistograms[i]->GetEntries() > 0) {
            for (int bin = 1; bin <= signalEfficiencyHistograms[i]->GetNbinsX(); ++bin) {
                double eff = signalEfficiencyHistograms[i]->GetBinContent(bin);
                if (eff > 0) {
                    if (eff < minEffSignal) minEffSignal = eff;
                    if (eff > maxEffSignal) maxEffSignal = eff;
                }
            }
        }
    }
    
    double yMinSignal = (minEffSignal > 0) ? minEffSignal * 0.5 : 1e-7;
    double yMaxSignal = (maxEffSignal > 0) ? maxEffSignal * 3.0 : 3.0;
    if (yMinSignal < 1e-7) yMinSignal = 1e-7;
    if (yMaxSignal > 3.0) yMaxSignal = 3.0;
    
    bool firstSignalDrawn = false;
    for (size_t i = 0; i < signalEfficiencyHistograms.size(); ++i) {
        if (signalEfficiencyHistograms[i]->GetEntries() > 0) {
            signalEfficiencyHistograms[i]->GetYaxis()->SetRangeUser(yMinSignal, yMaxSignal);
            signalEfficiencyHistograms[i]->GetYaxis()->SetTitle("Signal Selection Efficiency");
            signalEfficiencyHistograms[i]->GetYaxis()->SetTitleOffset(1.3);
            signalEfficiencyHistograms[i]->GetYaxis()->SetTitleSize(0.045);
            signalEfficiencyHistograms[i]->GetYaxis()->SetLabelSize(0.04);
            
            signalEfficiencyHistograms[i]->GetXaxis()->SetTitle("Cuts");
            signalEfficiencyHistograms[i]->GetXaxis()->SetTitleSize(0.045);
            signalEfficiencyHistograms[i]->GetXaxis()->SetLabelSize(0.04);
            signalEfficiencyHistograms[i]->GetXaxis()->SetTitleOffset(1.1);
            
            signalEfficiencyHistograms[i]->SetLineWidth(3);
            signalEfficiencyHistograms[i]->SetMarkerSize(1.5);
            
            signalEfficiencyHistograms[i]->Draw(firstSignalDrawn ? "LP SAME" : "LP");
            
            std::string legendLabel;
            if (workingPoints[i] == "125") {
                legendLabel = "1.25 % QCD Eff. WP";
            } else if (workingPoints[i] == "155") {
                legendLabel = "1.55 % QCD Eff. WP";
            } else {
                double wpValue = std::stod(workingPoints[i]) / 100.0;
                std::stringstream ss;
                ss << std::fixed << std::setprecision(2) << wpValue << " % QCD Eff. WP";
                legendLabel = ss.str();
            }
            
            legendSignal->AddEntry(signalEfficiencyHistograms[i], legendLabel.c_str(), "lp");
            firstSignalDrawn = true;
        }
    }
    
    legendSignal->Draw();
    DrawAtlasLabel();
    
    canvasSignal->SaveAs("../Plots/DB_Deco/CutFlow/SignalEfficiencyComparison_log.pdf");
    //canvasSignal->SaveAs("../Plots/DB_Deco/CutFlow/SignalEfficiencyComparison_log.png");
    
    // ==================== BACKGROUND EFFICIENCY CANVAS (LOG SCALE) ====================
    auto canvasBkg = new TCanvas("canvasBkg", "Background Efficiency Comparison", 800, 600);
    canvasBkg->SetLogy();
    
    auto legendBkg = new TLegend(0.60, 0.55, 0.9, 0.9);
    legendBkg->SetBorderSize(0);
    legendBkg->SetFillStyle(0);
    legendBkg->SetTextSize(0.030);
    
    double minEffBkg = 1.0, maxEffBkg = 0.0;
    for (size_t i = 0; i < bkgEfficiencyHistograms.size(); ++i) {
        if (bkgEfficiencyHistograms[i]->GetEntries() > 0) {
            for (int bin = 1; bin <= bkgEfficiencyHistograms[i]->GetNbinsX(); ++bin) {
                double eff = bkgEfficiencyHistograms[i]->GetBinContent(bin);
                if (eff > 0) {
                    if (eff < minEffBkg) minEffBkg = eff;
                    if (eff > maxEffBkg) maxEffBkg = eff;
                }
            }
        }
    }
    
    double yMinBkg = (minEffBkg > 0) ? minEffBkg * 0.5 : 1e-7;
    double yMaxBkg = (maxEffBkg > 0) ? maxEffBkg * 3.0 : 3.0;
    if (yMinBkg < 1e-7) yMinBkg = 1e-7;
    if (yMaxBkg > 3.0) yMaxBkg = 3.0;
    
    bool firstBkgDrawn = false;
    for (size_t i = 0; i < bkgEfficiencyHistograms.size(); ++i) {
        if (bkgEfficiencyHistograms[i]->GetEntries() > 0) {
            bkgEfficiencyHistograms[i]->GetYaxis()->SetRangeUser(yMinBkg, yMaxBkg);
            bkgEfficiencyHistograms[i]->GetYaxis()->SetTitle("Background Selection Efficiency");
            bkgEfficiencyHistograms[i]->GetYaxis()->SetTitleOffset(1.3);
            bkgEfficiencyHistograms[i]->GetYaxis()->SetTitleSize(0.045);
            bkgEfficiencyHistograms[i]->GetYaxis()->SetLabelSize(0.04);
            
            bkgEfficiencyHistograms[i]->GetXaxis()->SetTitle("Cuts");
            bkgEfficiencyHistograms[i]->GetXaxis()->SetTitleSize(0.045);
            bkgEfficiencyHistograms[i]->GetXaxis()->SetLabelSize(0.04);
            bkgEfficiencyHistograms[i]->GetXaxis()->SetTitleOffset(1.1);
            
            bkgEfficiencyHistograms[i]->SetLineWidth(3);
            bkgEfficiencyHistograms[i]->SetMarkerSize(1.5);
            
            bkgEfficiencyHistograms[i]->Draw(firstBkgDrawn ? "LP SAME" : "LP");
            
            std::string legendLabel;
            if (workingPoints[i] == "125") {
                legendLabel = "1.25 % QCD Eff. WP";
            } else if (workingPoints[i] == "155") {
                legendLabel = "1.55 % QCD Eff. WP";
            } else {
                double wpValue = std::stod(workingPoints[i]) / 100.0;
                std::stringstream ss;
                ss << std::fixed << std::setprecision(2) << wpValue << " % QCD Eff. WP";
                legendLabel = ss.str();
            }
            
            legendBkg->AddEntry(bkgEfficiencyHistograms[i], legendLabel.c_str(), "lp");
            firstBkgDrawn = true;
        }
    }
    
    legendBkg->Draw();
    DrawAtlasLabel();
    
    canvasBkg->SaveAs("../Plots/DB_Deco/CutFlow/BackgroundEfficiencyComparison_log.pdf");
    //canvasBkg->SaveAs("../Plots/DB_Deco/CutFlow/BackgroundEfficiencyComparison_log.png");
    
    // ==================== EXTENDED EFFICIENCY TABLE ====================
    std::cout << "\n+----------+------------------+------------------+------------------+\n";
    std::cout << "| WP       | Data Final Eff.  | Signal Final Eff.| Bkg Final Eff.   |\n";
    std::cout << "+----------+------------------+------------------+------------------+\n";
    
    for (size_t i = 0; i < efficiencyHistograms.size(); ++i) {
        if (efficiencyHistograms[i]->GetEntries() > 0) {
            double dataFinalEff = efficiencyHistograms[i]->GetBinContent(8);
            double signalFinalEff = signalEfficiencyHistograms[i]->GetBinContent(8);
            double bkgFinalEff = bkgEfficiencyHistograms[i]->GetBinContent(8);
            
            std::string wpDisplay;
            if (workingPoints[i] == "125") {
                wpDisplay = "1p25";
            } else if (workingPoints[i] == "155") {
                wpDisplay = "1p55";
            } else {
                wpDisplay = "0p" + workingPoints[i];
            }
            
            std::cout << "| " << std::left << std::setw(8) << wpDisplay << " | " 
                     << std::right << std::setw(16) << std::scientific << std::setprecision(4) << dataFinalEff << " | " 
                     << std::setw(16) << signalFinalEff << " | "
                     << std::setw(16) << bkgFinalEff << " |\n";
        }
    }
    
    std::cout << "+----------+------------------+------------------+------------------+\n";
    
    // ==================== DATA LINEAR SCALE CANVAS ====================
    auto canvasEffLinear = new TCanvas("canvasEffLinear", "Data Efficiency Comparison (Linear Scale)", 800, 600);
    canvasEffLinear->SetGridy();
    canvasEffLinear->SetGridx();
    
    auto legendLinear = new TLegend(0.60, 0.15, 0.9, 0.45);
    legendLinear->SetBorderSize(0);
    legendLinear->SetFillStyle(0);
    legendLinear->SetTextSize(0.030);
    
    std::vector<TH1F*> linearHistograms;
    bool firstLinearDrawn = false;
    
    for (size_t i = 0; i < efficiencyHistograms.size(); ++i) {
        if (efficiencyHistograms[i]->GetEntries() > 0) {
            TH1F* hLinear = static_cast<TH1F*>(efficiencyHistograms[i]->Clone(("linear_" + std::string(efficiencyHistograms[i]->GetName())).c_str()));
            
            hLinear->GetYaxis()->SetRangeUser(0, 1.1);
            hLinear->GetYaxis()->SetTitle("Data Selection Efficiency");
            hLinear->GetYaxis()->SetTitleOffset(1.3);
            hLinear->GetYaxis()->SetTitleSize(0.045);
            hLinear->GetYaxis()->SetLabelSize(0.04);
            
            hLinear->GetXaxis()->SetTitle("Cuts");
            hLinear->GetXaxis()->SetTitleSize(0.045);
            hLinear->GetXaxis()->SetLabelSize(0.04);
            hLinear->GetXaxis()->SetTitleOffset(1.1);
            
            hLinear->SetLineWidth(3);
            hLinear->SetMarkerSize(1.5);
            
            hLinear->Draw(firstLinearDrawn ? "LP SAME" : "LP");
            
            std::string legendLabel;
            if (workingPoints[i] == "125") {
                legendLabel = "1.25 % QCD Eff. WP";
            } else if (workingPoints[i] == "155") {
                legendLabel = "1.55 % QCD Eff. WP";
            } else {
                double wpValue = std::stod(workingPoints[i]) / 100.0;
                std::stringstream ss;
                ss << std::fixed << std::setprecision(2) << wpValue << " % QCD Eff. WP";
                legendLabel = ss.str();
            }
            legendLinear->AddEntry(hLinear, legendLabel.c_str(), "lp");
            
            linearHistograms.push_back(hLinear);
            firstLinearDrawn = true;
        }
    }
    
    legendLinear->Draw();
    DrawAtlasLabel();
    
    TLatex infoLabelLinear;
    infoLabelLinear.SetNDC();
    infoLabelLinear.SetTextFont(42);
    infoLabelLinear.SetTextSize(0.035);
    infoLabelLinear.DrawLatex(0.2, 0.78, "Efficiency relative to initial events");
    
    canvasEffLinear->SaveAs("../Plots/DB_Deco/CutFlow/DataEfficiencyComparison_Linear.pdf");
    //canvasEffLinear->SaveAs("../Plots/DB_Deco/CutFlow/DataEfficiencyComparison_Linear.png");
    
    // Clean up
    for (auto hist : linearHistograms) {
        delete hist;
    }
    delete canvasEffLinear;
    delete legendLinear;
    
    // Clean up efficiency histograms
    for (auto hist : efficiencyHistograms) {
        delete hist;
    }
    delete canvasEff;
    delete legend;
}


void CompareSignalToBackground() {
    gROOT->LoadMacro("../../Util/AtlasStyle.C");
    SetAtlasStyle();
    
    // Define working points to compare
    //std::vector<std::string> workingPoints = {"25", "30", "37", "46", "58", "74", "94", "125", "155"};
    std::vector<std::string> workingPoints = {"25", "46", "74"};
    
    // Get username for EOS path
    std::string username = std::getenv("USER") ? std::getenv("USER") : "tamezza";
    
    // Define EOS path
    //std::string basePath = "root://eosuser.cern.ch//eos/user/t/" + username + "/Documents/EasyBjets/ZbbJets/";
    std::string basePath = "../../";
    
    // Vectors to store S/B ratios for each cut
    std::vector<double> signalToBkg;
    std::vector<std::string> wpLabels;
    
    // Create canvas
    auto canvasSB = new TCanvas("canvasSB", "Signal to Background Ratio", 800, 600);
    canvasSB->SetGridy();
    
    // Create graph 
    auto graphSB = new TGraph();
    graphSB->SetTitle("Signal to Background Ratio;Working Point;S/B Ratio (Final Selection)");
    graphSB->SetMarkerStyle(20);
    graphSB->SetMarkerSize(1.5);
    graphSB->SetMarkerColor(kBlue);
    graphSB->SetLineColor(kBlue);
    graphSB->SetLineWidth(2);
    
    // Process each working point
    int pointIndex = 0;
    for (const auto& wp : workingPoints) {
        
        std::string fileWP, labelWP;
        if (wp == "125") {
            fileWP = "1p25";
            labelWP = "1p25";
        } else if (wp == "155") {
            fileWP = "1p55";
            labelWP = "1p55";
        } else {
            fileWP = "0p" + wp;
            labelWP = "0p" + wp;
        }
        
        // Construct file paths for this working point
        std::string dataPath = basePath + "NtupleSlim/DB_Deco_" + fileWP + "wp_slim/data_deco_" + fileWP + "wp_slim.root";
        std::string zbbPath = basePath + "NtupleSlim/DB_Deco_" + fileWP + "wp_slim/Zbby_deco_" + fileWP + "wp_slim.root";
        std::string mjPath = basePath + "NtupleSlim/DB_Deco_" + fileWP + "wp_slim/gamma_jets_deco_" + fileWP + "wp_slim.root";
        
        // Open ROOT files
        TFile* dataFile = TFile::Open(dataPath.c_str(), "READ");
        TFile* zbbFile = TFile::Open(zbbPath.c_str(), "READ");
        TFile* mjFile = TFile::Open(mjPath.c_str(), "READ");
        
        if (!dataFile || dataFile->IsZombie() || !zbbFile || zbbFile->IsZombie() || !mjFile || mjFile->IsZombie()) {
            std::cerr << "Error: Could not open files for working point " << labelWP << std::endl;
            if (dataFile) delete dataFile;
            if (zbbFile) delete zbbFile;
            if (mjFile) delete mjFile;
            continue;
        }
        
        // Get the cutflow histograms
        TH1F* dataCF = static_cast<TH1F*>(dataFile->Get("h_cutflow"));
        TH1F* zbbCF = static_cast<TH1F*>(zbbFile->Get("h_cutflow"));
        TH1F* mjCF = static_cast<TH1F*>(mjFile->Get("h_cutflow"));
        
        if (!dataCF || !zbbCF || !mjCF) {
            std::cerr << "Error: Could not find cutflow histograms for working point " << labelWP << std::endl;
            dataFile->Close(); zbbFile->Close(); mjFile->Close();
            delete dataFile; delete zbbFile; delete mjFile;
            continue;
        }
        
        // Calculate scale factor for gamma_jets_ (normalization bin is 4)
        double normBin = 4;
        double dataVal = dataCF->GetBinContent(normBin);
        double mjVal = mjCF->GetBinContent(normBin);
        double scaleFactor = (mjVal > 0) ? dataVal / mjVal : 1.0;
        mjCF->Scale(scaleFactor);
        
        // Get final bin values (bin 8 - after all cuts)
        int finalBin = 8;
        double signalVal = zbbCF->GetBinContent(finalBin);
        double bkgVal = mjCF->GetBinContent(finalBin);
        
        // Calculate S/B ratio
        double sbRatio = (bkgVal > 0) ? signalVal / bkgVal : 0;
        
        signalToBkg.push_back(sbRatio);
        wpLabels.push_back(labelWP);
        
        // Add point to graph - handle different WP value formats
        double wpValue;
        if (wp == "125") {
            wpValue = 1.25;
        } else if (wp == "155") {
            wpValue = 1.55;
        } else {
            wpValue = std::stod(wp) / 100.0; // Convert to decimal
        }
        graphSB->SetPoint(pointIndex++, wpValue, sbRatio);
        
        // Clean up
        dataFile->Close(); zbbFile->Close(); mjFile->Close();
        delete dataFile; delete zbbFile; delete mjFile;
    }
    
    // Draw the graph
    graphSB->Draw("APL");
    
    // Customize the graph appearance
    graphSB->GetXaxis()->SetTitle("Working Point Value");
    graphSB->GetYaxis()->SetTitle("Signal to Background Ratio");
    
    // Draw ATLAS label
    DrawAtlasLabel();
    
    // Save the plot
    canvasSB->SaveAs("../Plots/DB_Deco/CutFlow/SignalToBackgroundRatio.pdf");
    //canvasSB->SaveAs("../Plots/DB_Deco/CutFlow/SignalToBackgroundRatio.png");
    
    // Print S/B table
    std::cout << "\n+---------------+----------------+\n";
    std::cout << "| Working Point | S/B Ratio      |\n";
    std::cout << "+---------------+----------------+\n";
    
    for (size_t i = 0; i < wpLabels.size(); ++i) {
        std::cout << "| " << std::setw(13) << wpLabels[i] << " | " 
                  << std::right << std::setw(14) << std::fixed << std::setprecision(4) << signalToBkg[i] << " |\n";
    }
    
    std::cout << "+---------------+----------------+\n";
    
    // Clean up
    delete canvasSB;
    delete graphSB;
}

void CompareDataToMC() {
    gROOT->LoadMacro("../../Util/AtlasStyle.C");
    SetAtlasStyle();
    
    // Define working points to compare
    //std::vector<std::string> workingPoints = {"25", "30", "37", "46", "58", "74", "94", "125", "155"};
    std::vector<std::string> workingPoints = {"25", "46", "74"};
    
    // Get username for EOS path
    std::string username = std::getenv("USER") ? std::getenv("USER") : "tamezza";
    
    // Define EOS path
    //std::string eosPath = "root://eosuser.cern.ch//eos/user/t/" + username + "/Documents/EasyBjets/ZbbJets/";
    std::string eosPath = "../../";
    
    // Create a histogram to store the data/MC ratios for final cut
    TH1F* h_dataMCRatio = new TH1F("h_dataMCRatio", "Data/MC Ratio at Final Selection;Working Point;Data/MC", 
                                   workingPoints.size(), 0, workingPoints.size());
    
    // FIXED: Set bin labels to working point values correctly
    for (int i = 0; i < workingPoints.size(); ++i) {
        std::string binLabel;
        if (workingPoints[i] == "125") {
            binLabel = "1p25";
        } else if (workingPoints[i] == "155") {
            binLabel = "1p55";
        } else {
            binLabel = "0p" + workingPoints[i];
        }
        h_dataMCRatio->GetXaxis()->SetBinLabel(i+1, binLabel.c_str());
    }
    
    // Process each working point
    for (int i = 0; i < workingPoints.size(); ++i) {
        const auto& wp = workingPoints[i];

        std::string fileWP, labelWP;
        if (wp == "125") {
            fileWP = "1p25";
            labelWP = "1p25";
        } else if (wp == "155") {
            fileWP = "1p55";
            labelWP = "1p55";
        } else {
            fileWP = "0p" + wp;
            labelWP = "0p" + wp;
        }
        
        // Try to open the files for this working point
        std::string dataPath = eosPath + "NtupleSlim/DB_Deco_" + fileWP + "wp_slim/data_deco_" + fileWP + "wp_slim.root";
        std::string zbbPath = eosPath + "NtupleSlim/DB_Deco_" + fileWP + "wp_slim/Zbby_deco_" + fileWP + "wp_slim.root";
        std::string mjPath = eosPath + "NtupleSlim/DB_Deco_" + fileWP + "wp_slim/gamma_jets_deco_" + fileWP + "wp_slim.root";
        
        TFile* dataFile = TFile::Open(dataPath.c_str(), "READ");
        TFile* zbbFile = TFile::Open(zbbPath.c_str(), "READ");
        TFile* mjFile = TFile::Open(mjPath.c_str(), "READ");
        
        if (!dataFile || dataFile->IsZombie() || !zbbFile || zbbFile->IsZombie() || !mjFile || mjFile->IsZombie()) {
            std::cerr << "Error: Could not open files for working point " << labelWP << std::endl;
            if (dataFile) delete dataFile;
            if (zbbFile) delete zbbFile;
            if (mjFile) delete mjFile;
            continue;
        }
        
        // Get the cutflow histograms
        TH1F* dataCF = static_cast<TH1F*>(dataFile->Get("h_cutflow"));
        TH1F* zbbCF = static_cast<TH1F*>(zbbFile->Get("h_cutflow"));
        TH1F* mjCF = static_cast<TH1F*>(mjFile->Get("h_cutflow"));
        
        if (!dataCF || !zbbCF || !mjCF) {
            std::cerr << "Error: Could not find cutflow histograms for working point " << labelWP << std::endl;
            dataFile->Close(); zbbFile->Close(); mjFile->Close();
            delete dataFile; delete zbbFile; delete mjFile;
            continue;
        }
        
        // Calculate scale factor for gamma_jets_ (normalization bin is 4)
        double normBin = 4;
        double dataVal = dataCF->GetBinContent(normBin);
        double mjVal = mjCF->GetBinContent(normBin);
        double scaleFactor = (mjVal > 0) ? dataVal / mjVal : 1.0;
        mjCF->Scale(scaleFactor);
        
        // Get final bin values (bin 8 - after all cuts)
        int finalBin = 8;
        double dataFinal = dataCF->GetBinContent(finalBin);
        double zbbFinal = zbbCF->GetBinContent(finalBin);
        double mjFinal = mjCF->GetBinContent(finalBin);
        
        // Calculate total MC and data/MC ratio
        double totalMC = zbbFinal + mjFinal;
        double dataMCRatio = (totalMC > 0) ? dataFinal / totalMC : 0;
        
        // Fill the histogram
        h_dataMCRatio->SetBinContent(i+1, dataMCRatio);
        
        // Set error
        double dataError = sqrt(dataFinal); // Assume Poisson error for data
        double mcError = sqrt(zbbFinal*zbbFinal + mjFinal*mjFinal); // Sum in quadrature
        double ratioError = (totalMC > 0) ? 
                           sqrt(pow(dataError/totalMC, 2) + pow(dataFinal*mcError/(totalMC*totalMC), 2)) : 0;
        h_dataMCRatio->SetBinError(i+1, ratioError);
        
        // Clean up
        dataFile->Close(); zbbFile->Close(); mjFile->Close();
        delete dataFile; delete zbbFile; delete mjFile;
    }
    
    // Style the histogram
    h_dataMCRatio->SetMarkerStyle(21);
    h_dataMCRatio->SetMarkerColor(kBlue);
    h_dataMCRatio->SetLineColor(kBlue);
    h_dataMCRatio->SetLineWidth(2);
    h_dataMCRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
    
    // Create canvas
    auto canvasDMC = new TCanvas("canvasDMC", "Data/MC Ratio", 1000, 800);
    canvasDMC->SetGridy();
    
    // Draw histogram
    h_dataMCRatio->Draw("EP");
    
    // Add a horizontal line at 1.0
    TLine* line = new TLine(0, 1.0, workingPoints.size(), 1.0);
    line->SetLineColor(kRed);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw();
    
    // Draw ATLAS label
    DrawAtlasLabel();
    
    // Save the plot
    canvasDMC->SaveAs("../Plots/DB_Deco/CutFlow/DataMCRatio.pdf");
    //canvasDMC->SaveAs("../Plots/DB_Deco/CutFlow/DataMCRatio.png");
    
    // FIXED: Print Data/MC ratio table with correct labels
    std::cout << "\n+--------------+----------------+\n";
    std::cout << "| Working Point | Data/MC Ratio |\n";
    std::cout << "+--------------+----------------+\n";
    
    for (int i = 0; i < workingPoints.size(); ++i) {
        std::string wpDisplay;
        if (workingPoints[i] == "125") {
            wpDisplay = "1p25";
        } else if (workingPoints[i] == "155") {
            wpDisplay = "1p55";
        } else {
            wpDisplay = "0p" + workingPoints[i];
        }
        
        std::cout << "| " << std::left << std::setw(12) << wpDisplay << " | " 
                 << std::right << std::setw(14) << std::fixed << std::setprecision(4) 
                 << h_dataMCRatio->GetBinContent(i+1) << " ± " 
                 << std::fixed << std::setprecision(4) << h_dataMCRatio->GetBinError(i+1) << " |\n";
    }
    
    std::cout << "+--------------+----------------+\n";
    
    // Clean up
    delete canvasDMC;
    delete h_dataMCRatio;
    delete line;
}


/* --------------------------------- CutFlow --------------------------------------------------------- */

void IndividualStackedCutflows() {
    // Try to load AtlasStyle.C, but handle the case if it's not found
    bool atlasStyleLoaded = true;
    try {
        gROOT->LoadMacro("../../Util/AtlasStyle.C");
        SetAtlasStyle();
    } catch (...) {
        std::cerr << "Warning: AtlasStyle.C not found. Using default ROOT style." << std::endl;
        atlasStyleLoaded = false;
        // Set a reasonable default style
        gStyle->SetOptStat(0);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gStyle->SetLegendBorderSize(0);
    }
    
    // Define working points to compare
    //std::vector<std::string> workingPoints = {"25", "30", "37", "46", "58", "74", "94", "125", "155"};
    std::vector<std::string> workingPoints = {"25", "46", "74"};
    
    // Get username for EOS path
    std::string username = std::getenv("USER") ? std::getenv("USER") : "tamezza";
    
    // Define EOS path
    //std::string eosPath = "root://eosuser.cern.ch//eos/user/t/" + username + "/Documents/EasyBjets/ZbbJets/";
    std::string eosPath = "../../";
    
    // Define cut names
    std::vector<std::string> cutNames = {
        "Initial Events",          // Cut 0
        "D_{xbb} WP",              // Cut 1
        "p_{T}^{J} > 450",         // Cut 2
        "|#eta^{J}| < 1.2",        // Cut 3
        "Small-R p_{T} & #eta",    // Cut 4
        "JVT > 0.59",              // Cut 5
        "|#alpha - #pi| < 0.3",    // Cut 6
        "#Delta R (J1, DB) > 1.5" // Cut 7
    };
    
    // Process each working point
    for (const auto& wp : workingPoints) {
        std::string fileWP, labelWP;
        if (wp == "125") {
            fileWP = "1p25";
            labelWP = "1p25";
        } else if (wp == "155") {
            fileWP = "1p55";
            labelWP = "1p55";
        } else {
            fileWP = "0p" + wp;
            labelWP = "0p" + wp;
        }
        
        std::cout << "Creating stacked cutflow for working point " << labelWP << std::endl;
        
        // Try to open the files for this working point
        std::string dataPath = eosPath + "NtupleSlim/DB_Deco_" + fileWP + "wp_slim/data_deco_" + fileWP + "wp_slim.root";
        std::string zbbPath = eosPath + "NtupleSlim/DB_Deco_" + fileWP + "wp_slim/Zbby_deco_" + fileWP + "wp_slim.root";
        std::string mjPath = eosPath + "NtupleSlim/DB_Deco_" + fileWP + "wp_slim/gamma_jets_deco_" + fileWP + "wp_slim.root";
        
        TFile* dataFile = TFile::Open(dataPath.c_str(), "READ");
        TFile* zbbFile = TFile::Open(zbbPath.c_str(), "READ");
        TFile* mjFile = TFile::Open(mjPath.c_str(), "READ");
        
        if (!dataFile || dataFile->IsZombie() || !zbbFile || zbbFile->IsZombie() || !mjFile || mjFile->IsZombie()) {
            std::cerr << "Error: Could not open files for working point " << labelWP << std::endl;
            if (dataFile) delete dataFile;
            if (zbbFile) delete zbbFile;
            if (mjFile) delete mjFile;
            continue;
        }
        
        // Get the cutflow histograms
        TH1F* dataCF = static_cast<TH1F*>(dataFile->Get("h_cutflow"));
        TH1F* zbbCF = static_cast<TH1F*>(zbbFile->Get("h_cutflow"));
        TH1F* mjCF = static_cast<TH1F*>(mjFile->Get("h_cutflow"));
        
        if (!dataCF || !zbbCF || !mjCF) {
            std::cerr << "Error: Could not find cutflow histograms for working point " << labelWP << std::endl;
            dataFile->Close(); zbbFile->Close(); mjFile->Close();
            delete dataFile; delete zbbFile; delete mjFile;
            continue;
        }
        
        // Calculate scale factor for gamma_jets_ (normalization bin is 4)
        double normBin = 4;
        double dataVal = dataCF->GetBinContent(normBin);
        double mjVal = mjCF->GetBinContent(normBin);
        double scaleFactor = (mjVal > 0) ? dataVal / mjVal : 1.0;
        //mjCF->Scale(scaleFactor);
        mjCF->Scale(1.);
        
        // Create copies of histograms to work with
        TH1F* hData = static_cast<TH1F*>(dataCF->Clone("h_data_cutflow"));
        TH1F* hZbb = static_cast<TH1F*>(zbbCF->Clone("h_zbb_cutflow"));
        TH1F* hMJ = static_cast<TH1F*>(mjCF->Clone("h_mj_cutflow"));
        
        // Set histogram styles
        hData->SetMarkerStyle(20);
        hData->SetMarkerSize(1.2);
        hData->SetMarkerColor(kBlack);
        hData->SetLineColor(kBlack);
        hData->SetLineWidth(2);
        
        hZbb->SetFillColor(kRed);
        hZbb->SetLineColor(kRed);
        hZbb->SetLineWidth(1);
        
        hMJ->SetFillColor(kBlue);
        hMJ->SetLineColor(kBlue);
        hMJ->SetLineWidth(1);

        // Corrected logic to display working point as "0.25 % QCD Eff. WP" for bin labels
        std::string WP_Label_display_for_bins;
        if (wp == "125") {
            WP_Label_display_for_bins = "1.25 % QCD Eff. WP";
        } else if (wp == "155") {
            WP_Label_display_for_bins = "1.55 % QCD Eff. WP";
        } else {
            double wpValue = std::stod(wp) / 100.0;
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << wpValue << " % QCD Eff. WP";
            WP_Label_display_for_bins = ss.str();
        }
        
        // FIXED: Set bin labels with correct working point format
        for (int bin = 1; bin <= cutNames.size() && bin <= hData->GetNbinsX(); ++bin) {
            if (bin == 2) { // Special case for b-tagging WP
                /*hData->GetXaxis()->SetBinLabel(bin, ("D_{xbb} " + labelWP + " WP").c_str());
                hZbb->GetXaxis()->SetBinLabel(bin, ("D_{xbb} " + labelWP + " WP").c_str());
                hMJ->GetXaxis()->SetBinLabel(bin, ("D_{xbb} " + labelWP + " WP").c_str());*/
                hData->GetXaxis()->SetBinLabel(bin, ("D_{xbb} " + WP_Label_display_for_bins).c_str());
                hZbb->GetXaxis()->SetBinLabel(bin, ("D_{xbb} " + WP_Label_display_for_bins).c_str());
                hMJ->GetXaxis()->SetBinLabel(bin, ("D_{xbb} " + WP_Label_display_for_bins).c_str());
            } else {
                hData->GetXaxis()->SetBinLabel(bin, cutNames[bin-1].c_str());
                hZbb->GetXaxis()->SetBinLabel(bin, cutNames[bin-1].c_str());
                hMJ->GetXaxis()->SetBinLabel(bin, cutNames[bin-1].c_str());
            }
        }


        // Corrected logic to display working point as "0.25 % QCD Eff. WP"
        /*std::string WP_Label;
        if (wp == "125") {
            WP_Label = "1.25 % QCD Eff. WP";
        } else if (wp == "155") {
            WP_Label = "1.55 % QCD Eff. WP";
        } else {
            double wpValue = std::stod(wp) / 100.0; // Use 'wp' (current string)
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << wpValue << " % QCD Eff. WP";
            WP_Label = ss.str();
        }*/

        // FIXED: Create canvas with correct naming
        TCanvas* canvas = new TCanvas(("canvasStackedCF_" + labelWP).c_str(), 
                                      ("Stacked Cutflow (WP: " + labelWP + ")").c_str(), 800, 600);
        canvas->cd();
        
        // Create upper pad for stacked plot
        TPad* padTop = new TPad("padTop", "padTop", 0, 0.3, 1, 1);
        padTop->SetBottomMargin(0.02);
        padTop->SetLogy();
        padTop->Draw();
        padTop->cd();
        
        // Create stacked histogram
        THStack* stack = new THStack(("stackCF_" + labelWP).c_str(), 
                                     ("Cutflow (WP: " + labelWP + ");Cut;Events").c_str());
        
        // Add histograms to stack (signal first, then background)
        stack->Add(hZbb);
        stack->Add(hMJ);
        
        // Find max Y value for plot range
        double maxY = hData->GetMaximum();
        
        // Draw stack and data
        stack->Draw("HIST");
        stack->SetMaximum(maxY * 1000); // Give some room at the top
        stack->SetMinimum(0.1); // Small value for log scale
        
        // Fix axis labels (they get overwritten in the stack)
        stack->GetXaxis()->SetLabelSize(0);
        
        // Draw data points
        hData->Draw("EP SAME");
        
        // Create legend
        TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->AddEntry(hData, "Data", "lep");
        legend->AddEntry(hZbb, "Z(#rightarrowbb)+ #gamma", "f");
        legend->AddEntry(hMJ, "#gamma + jets", "f");
        legend->Draw();
        
        // Draw ATLAS label
        DrawAtlasLabel();
        
        // Add working point label
        TLatex wpLabel;
        wpLabel.SetNDC();
        wpLabel.SetTextFont(42);
        wpLabel.SetTextSize(0.04);
        //wpLabel.DrawLatex(0.2, 0.78, ("Working Point: " + labelWP).c_str());
        
        // Corrected logic to display working point as "0.25 % QCD Eff. WP"
        std::string WP_Label;
        if (wp == "125") {
            WP_Label = "1.25 % QCD Eff. WP";
        } else if (wp == "155") {
            WP_Label = "1.55 % QCD Eff. WP";
        } else {
            double wpValue = std::stod(wp) / 100.0; // Use 'wp' (current string)
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << wpValue << " % QCD Eff. WP";
            WP_Label = ss.str();
        }

        wpLabel.DrawLatex(0.2, 0.78, ("Working Point: " + WP_Label).c_str());
        
        // Create lower pad for ratio
        canvas->cd();
        TPad* padBottom = new TPad("padBottom", "padBottom", 0, 0.05, 1, 0.3);
        padBottom->SetTopMargin(0.05);
        padBottom->SetBottomMargin(0.3);
        padBottom->Draw();
        padBottom->cd();
        
        // Create data/MC ratio histogram
        TH1F* hRatio = static_cast<TH1F*>(hData->Clone("h_ratio"));
        
        // Create sum of all MC for denominator
        TH1F* hSumMC = static_cast<TH1F*>(hZbb->Clone("h_sumMC"));
        hSumMC->Add(hMJ);
        
        // Calculate ratio
        hRatio->Divide(hSumMC);
        
        // Configure ratio plot
        hRatio->SetTitle("");
        hRatio->GetYaxis()->SetTitle("Data/MC");
        hRatio->GetYaxis()->SetTitleSize(0.10);
        hRatio->GetYaxis()->SetTitleOffset(0.5);
        hRatio->GetYaxis()->SetLabelSize(0.08);
        hRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
        hRatio->GetYaxis()->SetNdivisions(505);
        
        hRatio->GetXaxis()->SetTitle("Cuts");
        hRatio->GetXaxis()->SetTitleSize(0.10);
        hRatio->GetXaxis()->SetTitleOffset(1.2);
        hRatio->GetXaxis()->SetLabelSize(0.10);
        
        // Draw ratio
        hRatio->Draw("EP");
        
        // Add a reference line at 1
        TLine* line = new TLine(0, 1, cutNames.size(), 1);
        line->SetLineStyle(2);
        line->SetLineColor(kBlack);
        line->SetLineWidth(2);
        line->Draw("SAME");
        
        // FIXED: Save canvas with correct file names
        canvas->SaveAs(("../Plots/DB_Deco/CutFlow/StackedCutflow_WP" + labelWP + ".pdf").c_str());
        //canvas->SaveAs(("../Plots/DB_Deco/CutFlow/StackedCutflow_WP" + labelWP + ".png").c_str());
        
        // FIXED: Print cutflow table with correct working point
        std::cout << "\nCutflow table for working point " << labelWP << ":" << std::endl;
        std::cout << "+------+----------------------+-----------------+--------------------+---------------------+----------------+\n";
        std::cout << "| Cut  | Description          | Data            | Z->bb + #gamma     | #gamma + jets       | Data/MC        |\n";
        std::cout << "+------+----------------------+-----------------+--------------------+---------------------+----------------+\n";
        
        // Open a file to save the cutflow table
        std::ofstream cutflowFile;
        cutflowFile.open(("../Plots/DB_Deco/CutFlow/CutflowTable_WP" + labelWP + ".txt").c_str());
        
        cutflowFile << "Cutflow table for working point " << labelWP << ":\n";
        cutflowFile << "+------+----------------------+-----------------+--------------------+---------------------+----------------+\n";
        cutflowFile << "| Cut  | Description          | Data            | Z->bb + #gamma     | #gamma + jets       | Data/MC        |\n";
        cutflowFile << "+------+----------------------+-----------------+--------------------+---------------------+----------------+\n";
        
        for (int i = 1; i <= cutNames.size() && i <= hData->GetNbinsX(); ++i) {
            double dataContent = hData->GetBinContent(i);
            double zbbContent = hZbb->GetBinContent(i);
            double mjContent = hMJ->GetBinContent(i);
            double totalMC = zbbContent + mjContent;
            double dataToMC = (totalMC > 0) ? dataContent / totalMC : 0.0;
            
            // FIXED: Use correct cut label format
            std::string cutLabel = (i == 2) ? "D_{xbb} " + labelWP + " WP" : cutNames[i-1];
            
            // Format for console output
            std::ostringstream line_stream;
            line_stream << "| " << std::setw(4) << i-1 << " | " 
                       << std::left << std::setw(20) << cutLabel << " | " 
                       << std::right << std::setw(9) << std::fixed << std::setprecision(1) << dataContent << " | " 
                       << std::setw(12) << std::fixed << std::setprecision(1) << zbbContent << " | " 
                       << std::setw(13) << std::fixed << std::setprecision(1) << mjContent << " | " 
                       << std::setw(8) << std::fixed << std::setprecision(3) << dataToMC << " |";
                       
            std::cout << line_stream.str() << std::endl;
            cutflowFile << line_stream.str() << "\n";
        }
        
        std::cout << "+------+----------------------+-----------+--------------+---------------+----------+\n";
        cutflowFile << "+------+----------------------+-----------+--------------+---------------+----------+\n";
        
        // Calculate some useful metrics
        int finalBin = cutNames.size();
        double finalDataEvents = hData->GetBinContent(finalBin);
        double finalZbbEvents = hZbb->GetBinContent(finalBin);
        double finalMjEvents = hMJ->GetBinContent(finalBin);
        double finalTotalMC = finalZbbEvents + finalMjEvents;
        double finalDataToMC = (finalTotalMC > 0) ? finalDataEvents / finalTotalMC : 0.0;
        double finalStoB = (finalMjEvents > 0) ? finalZbbEvents / finalMjEvents : 0.0;
        double efficiencyData = (hData->GetBinContent(1) > 0) ? finalDataEvents / hData->GetBinContent(1) : 0.0;
        
        // Format summary output
        std::ostringstream summary;
        summary << "\nSummary for working point " << labelWP << ":\n";
        summary << "  - Final data events: " << finalDataEvents << "\n";
        summary << "  - Final Z->bb events: " << finalZbbEvents << "\n";
        summary << "  - Final gamma_jets_ events: " << finalMjEvents << "\n";
        summary << "  - Final data/MC ratio: " << finalDataToMC << "\n";
        summary << "  - Final S/B ratio: " << finalStoB << "\n";
        summary << "  - Overall efficiency: " << std::fixed << std::setprecision(6) << efficiencyData << "\n";
        
        std::cout << summary.str();
        cutflowFile << summary.str();
        
        // Save relative efficiencies for each Cuts
        cutflowFile << "\nRelative efficiencies at each step (normalized to initial events):\n";
        cutflowFile << "+------+----------------------+-----------------+-----------------+-----------------+\n";
        cutflowFile << "| Cut  | Description          | Data            | Z->bb + #gamma  | #gamma + jets   |\n";
        cutflowFile << "+------+----------------------+-----------------+-----------------+-----------------+\n";
        
        double initialDataCount = hData->GetBinContent(1);
        double initialZbbCount = hZbb->GetBinContent(1);
        double initialMjCount = hMJ->GetBinContent(1);
        
        for (int i = 1; i <= cutNames.size() && i <= hData->GetNbinsX(); ++i) {
            double dataEff = (initialDataCount > 0) ? hData->GetBinContent(i) / initialDataCount : 0.0;
            double zbbEff = (initialZbbCount > 0) ? hZbb->GetBinContent(i) / initialZbbCount : 0.0;
            double mjEff = (initialMjCount > 0) ? hMJ->GetBinContent(i) / initialMjCount : 0.0;
            
            std::string cutLabel = (i == 2) ? "D_{xbb} " + labelWP + " WP" : cutNames[i-1];
            
            cutflowFile << "| " << std::setw(4) << i-1 << " | " 
                        << std::left << std::setw(20) << cutLabel << " | " 
                        << std::right << std::setw(9) << std::fixed << std::setprecision(6) << dataEff << " | " 
                        << std::setw(9) << std::fixed << std::setprecision(6) << zbbEff << " | " 
                        << std::setw(9) << std::fixed << std::setprecision(6) << mjEff << " |\n";
        }
        
        cutflowFile << "+------+----------------------+-----------+-----------+-----------+\n";
        
        // Close the file
        cutflowFile.close();
        std::cout << "Cutflow table saved to: " << ("../Plots/DB_Deco/CutFlow/CutflowTable_WP" + labelWP + ".txt") << std::endl;
        
        // Create efficiency histograms
        TH1F* hDataEff = static_cast<TH1F*>(hData->Clone("h_data_eff"));
        TH1F* hZbbEff = static_cast<TH1F*>(hZbb->Clone("h_zbb_eff"));
        TH1F* hMJEff = static_cast<TH1F*>(hMJ->Clone("h_mj_eff"));
        
        // Calculate efficiency relative to initial bin
        double dataInitial = hData->GetBinContent(1);
        double zbbInitial = hZbb->GetBinContent(1);
        double mjInitial = hMJ->GetBinContent(1);
        
        if (dataInitial > 0 && zbbInitial > 0 && mjInitial > 0) {
            for (int bin = 1; bin <= hData->GetNbinsX(); ++bin) {
                hDataEff->SetBinContent(bin, hData->GetBinContent(bin) / dataInitial);
                hZbbEff->SetBinContent(bin, hZbb->GetBinContent(bin) / zbbInitial);
                hMJEff->SetBinContent(bin, hMJ->GetBinContent(bin) / mjInitial);
                
                // Set reasonable errors
                hDataEff->SetBinError(bin, hData->GetBinError(bin) / dataInitial);
                hZbbEff->SetBinError(bin, hZbb->GetBinError(bin) / zbbInitial);
                hMJEff->SetBinError(bin, hMJ->GetBinError(bin) / mjInitial);
            }
            
            // Create efficiency plot
            TCanvas* canvasEff = new TCanvas(("canvasEfficiency_" + labelWP).c_str(), 
                                          ("Efficiency (WP: " + labelWP + ")").c_str(), 800, 600);
            canvasEff->SetGridy();
            canvasEff->SetLogy(); // Log scale for efficiency
            
            // Create a legend
            TLegend* legendEff = new TLegend(0.7, 0.7, 0.9, 0.9);
            legendEff->SetBorderSize(0);
            legendEff->SetFillStyle(0);
            
            // Set axis titles
            hDataEff->SetTitle(("Data Selection Efficiency (WP: " + labelWP + ")").c_str());
            hDataEff->GetYaxis()->SetTitle("Efficiency (log scale)");
            hDataEff->GetYaxis()->SetRangeUser(1e-6, 2); // Reasonable range for log scale
            
            // Draw histograms
            hDataEff->Draw("EP");
            hZbbEff->Draw("HIST SAME");
            hMJEff->Draw("HIST SAME");
            
            // Add to legend
            legendEff->AddEntry(hDataEff, "Data", "lep");
            legendEff->AddEntry(hZbbEff, "Z(#rightarrowbb)+ #gamma", "l");
            legendEff->AddEntry(hMJEff, "#gamma + jets", "l");
            legendEff->Draw();
            
            // Draw ATLAS label
            DrawAtlasLabel();
            
            // Add working point label
            TLatex wpLabelEff;
            wpLabelEff.SetNDC();
            wpLabelEff.SetTextFont(42);
            wpLabelEff.SetTextSize(0.04);
            wpLabelEff.DrawLatex(0.2, 0.78, ("Working Point: " + labelWP).c_str());
            
            // Save the plot
            canvasEff->SaveAs(("../Plots/DB_Deco/CutFlow/Efficiency_WP" + labelWP + ".pdf").c_str());
            //canvasEff->SaveAs(("../Plots/DB_Deco/CutFlow/Efficiency_WP" + labelWP + ".png").c_str());
            
            // Clean up
            delete canvasEff;
            delete legendEff;
        }
        
        // Clean up all histograms and objects
        delete hDataEff;
        delete hZbbEff;
        delete hMJEff;
        delete hData;
        delete hZbb;
        delete hMJ;
        delete hRatio;
        delete hSumMC;
        delete canvas;
        delete legend;
        delete line;
        
        dataFile->Close(); 
        zbbFile->Close(); 
        mjFile->Close();
        delete dataFile; 
        delete zbbFile; 
        delete mjFile;
    }
}


// Add to the main function
void CutFlow_Deco() {
    // Try to load AtlasStyle.C, but handle the case if it's not found
    bool atlasStyleLoaded = true;
    try {
        gROOT->LoadMacro("../../Util/AtlasStyle.C");
        SetAtlasStyle();
    } catch (...) {
        std::cerr << "Warning: AtlasStyle.C not found. Using default ROOT style." << std::endl;
        atlasStyleLoaded = false;
        // Set a reasonable default style
        gStyle->SetOptStat(0);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gStyle->SetLegendBorderSize(0);
    }
    
    // Define the output path for plots
    std::string outputPath = "../Plots/DB_Deco/CutFlow/";
    
    // Create output directory if it doesn't exist
    gSystem->Exec(("mkdir -p " + outputPath).c_str());
    
    CompareEfficiencies();
    // Run the comparison analyses
    std::cout << "Comparing Efficiencies..." << std::endl;
    CompareEfficienciesData();
    
    std::cout << "Comparing Signal to Background Ratios..." << std::endl;
    CompareSignalToBackground();
    
    std::cout << "Comparing Data to MC Ratios..." << std::endl;
    CompareDataToMC();
    
    std::cout << "Creating Individual Stacked Cutflows..." << std::endl;
    IndividualStackedCutflows();

    //std::cout << "Creating Working Point Comparison Summary..." << std::endl;
    //CreateCutflowSummary();
    
    std::cout << "All comparisons complete." << std::endl;
}
