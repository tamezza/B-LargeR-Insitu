#include "../../Util/AtlasStyle.C"

void DrawAtlasLabel() {
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(72);
    latex.SetTextColor(kBlack);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.2, 0.88, "ATLAS");

    latex.SetTextFont(42);
    latex.DrawLatex(0.3, 0.88, "Internal");

    latex.DrawLatex(0.2, 0.84, "#sqrt{s} = 13 TeV, 140 fb^{-1}");
}


void CompareEfficiencies() {
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
    //std::vector<int> colors = {kRed, kBlue, kGreen+2, kMagenta, kCyan+2, kOrange+7, kViolet+3, kPink+2, kSpring+4};
    std::vector<int> colors = {kRed, kMagenta, kOrange+7};
    
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
        "D_{xbb} Cuts",              // Cut 1
        "p_{T}^{J} > 450",         // Cut 2
        "|#eta^{J}| < 1.2",        // Cut 3
        "Small-R p_{T} & #eta",    // Cut 4
        "JVT > 0.59",              // Cut 5
        "|#alpha - #pi| < 0.3",    // Cut 6
        "#Delta R (J1, MJB) > 1.5" // Cut 7
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
        std::string filePath = eosPath + "NtupleSlim/MJB_Deco_" + fileWP + "wp_slim/data_deco_" + fileWP + "wp_slim.root";
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
        for (int bin = 1; bin <= 8 && bin <= cutflowHist->GetNbinsX(); ++bin) {
            double count = cutflowHist->GetBinContent(bin);
            double efficiency = (initialCount > 0) ? count / initialCount : 0;
            efficiencyHistograms[i]->SetBinContent(bin, efficiency);

            if (bin == 2) {
                efficiencyHistograms[i]->GetXaxis()->SetBinLabel(bin, ("D_{xbb} ")); // ("D_{xbb} " + labelWP + " WP").c_str()
            } else {
                efficiencyHistograms[i]->GetXaxis()->SetBinLabel(bin, cutNames[bin-1].c_str());
            }
        }
        
        file->Close();
        delete file;
        
        // Process SIGNAL file
        std::string signalPath = eosPath + "NtupleSlim/MJB_Deco_" + fileWP + "wp_slim/ZbbJets_deco_" + fileWP + "wp_slim.root";
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
        for (int bin = 1; bin <= 8 && bin <= signalCutflow->GetNbinsX(); ++bin) {
            double count = signalCutflow->GetBinContent(bin);
            double efficiency = (signalInitialCount > 0) ? count / signalInitialCount : 0;
            signalEfficiencyHistograms[i]->SetBinContent(bin, efficiency);

            if (bin == 2) {
                signalEfficiencyHistograms[i]->GetXaxis()->SetBinLabel(bin, ("D_{xbb} ")); // ("D_{xbb} " + labelWP + " WP").c_str()
            } else {
                signalEfficiencyHistograms[i]->GetXaxis()->SetBinLabel(bin, cutNames[bin-1].c_str());
            }
        }
        
        signalFile->Close();
        delete signalFile;
        
        // Process BACKGROUND file
        std::string bkgPath = eosPath + "NtupleSlim/MJB_Deco_" + fileWP + "wp_slim/Multijets_deco_" + fileWP + "wp_slim.root";
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
        for (int bin = 1; bin <= 8 && bin <= bkgCutflow->GetNbinsX(); ++bin) {
            double count = bkgCutflow->GetBinContent(bin);
            double efficiency = (bkgInitialCount > 0) ? count / bkgInitialCount : 0;
            bkgEfficiencyHistograms[i]->SetBinContent(bin, efficiency);

            if (bin == 2) {
                bkgEfficiencyHistograms[i]->GetXaxis()->SetBinLabel(bin, ("D_{xbb} ")); // ("D_{xbb} " + labelWP + " WP").c_str()
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
    
    canvasEff->SaveAs("../Plots/MJB_Deco/Efficiency/DataEfficiencyComparison_log.pdf");
    //canvasEff->SaveAs("../Plots/MJB_Deco/Efficiency/DataEfficiencyComparison_log.png");
    
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
    
    canvasSignal->SaveAs("../Plots/MJB_Deco/Efficiency/SignalEfficiencyComparison_log.pdf");
    //canvasSignal->SaveAs("../Plots/MJB_Deco/Efficiency/SignalEfficiencyComparison_log.png");
    
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
    
    canvasBkg->SaveAs("../Plots/MJB_Deco/Efficiency/BackgroundEfficiencyComparison_log.pdf");
    //canvasBkg->SaveAs("../Plots/MJB_Deco/Efficiency/BackgroundEfficiencyComparison_log.png");
    
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
    auto canvasEffLinear = new TCanvas("canvasEffLinear", "Data Efficiency Comparison (Linear Scale)", 1200, 900);
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
    infoLabelLinear.DrawLatex(0.2, 0.77, "Efficiency relative to initial events");
    
    canvasEffLinear->SaveAs("../Plots/MJB_Deco/Efficiency/DataEfficiencyComparison_Linear.pdf");
    //canvasEffLinear->SaveAs("../Plots/MJB_Deco/Efficiency/DataEfficiencyComparison_Linear.png");
    
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

// Add to the main function
void Efficiency_Deco() {
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
    std::string outputPath = "../Plots/MJB_Deco/Efficiency/";
    
    // Create output directory if it doesn't exist
    gSystem->Exec(("mkdir -p " + outputPath).c_str());
    
    // Run the comparison analyses
    std::cout << "Comparing Efficiencies..." << std::endl;
    CompareEfficiencies();
    
    std::cout << "All comparisons complete." << std::endl;
}
