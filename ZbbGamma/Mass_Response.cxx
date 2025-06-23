
/*
 * b-Jets Plotting Code
 * 
 * This script analyzes large-R (bb)-jets samples and data
 * It processes jet mass distributions, 
 * And direct balance ratio, creating comparison plots with proper uncertainty bands.
 */
#include "../../Util/AtlasStyle.C"

/**
 * Draws the standard ATLAS label on plots
 */
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

/**
 * Normalizes all histograms to unit area
 * @param histograms Vector of histograms to normalize
 */
void NormalizeHistograms(std::vector<TH1F*>& histograms) {
    for (auto& hist : histograms) {
        if (hist->Integral() != 0) {
            hist->Scale(1.0 / hist->Integral());
        }
    }
}

/**
 * Structure to manage working point configurations
 */
/*
struct WorkingPoint {
    std::string name;                       // Working point name (e.g., "46")
    std::vector<std::string> inputFiles;    // List of input ROOT files
    std::string outputDir;                  // Output directory for plots
*/
    /**
     * Constructor
     * @param wp Working point identifier
     */
    /*
    WorkingPoint(const std::string& wp) : name(wp) {
        // Define input files for the working point
        inputFiles = {
            "../../Data/DB/data_0p" + wp + "wp_slim.root",          
            "../../Data/DB/Multijets_0p" + wp + "wp_slim.root",      
            "../../Data/DB/ZbbJets_0p" + wp + "wp_slim.root"         
            
        };
        outputDir = "../Plots/DB/" + wp + "/MassRdb/";
    }
};
*/

struct WorkingPoint {
    std::string name;                       // Working point name (e.g., "46")
    std::vector<std::string> inputFiles;    // List of input ROOT files
    std::string outputDir;                  // Output directory for plots

    /**
     * Constructor
     * @param wp Working point identifier
     */
    WorkingPoint(const std::string& wp) : name(wp) {
        // Create proper ROOT URL for remote access
        // Format: root://eosuser.cern.ch//eos/user/t/tamezza/...
        // Get current username for path construction 
        std::string username = std::getenv("USER") ? std::getenv("USER") : "tamezza";
        
        // APPROACH 1: When running directly on lxplus in the same session /eos/user/t/tamezza/
        //std::string eosPath = "root://eosuser.cern.ch//eos/user/t/" + username + "/Documents/EasyBjets/ZbbJets/";
        std::string eosPath = "../../";
        
        // Define input files for the working point with proper ROOT URLs
        /*inputFiles = {
            eosPath + "NtupleSlim/DB_0p" + wp + "wp_slim/data_0p" + wp + "wp_slim.root",          
            eosPath + "NtupleSlim/DB_0p" + wp + "wp_slim/Multijets_0p" + wp + "wp_slim.root",      
            eosPath + "NtupleSlim/DB_0p" + wp + "wp_slim/ZbbJets_0p" + wp + "wp_slim.root"         
        };
        
        outputDir = "../Plots/DB/" + wp + "/MassRdb/";*/
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

        inputFiles = {
            eosPath + "NtupleSlim/DB_" + fileWP + "wp_slim/data_" + fileWP + "wp_slim.root",          
            eosPath + "NtupleSlim/DB_" + fileWP + "wp_slim/gamma_jets_" + fileWP + "wp_slim.root",      
            eosPath + "NtupleSlim/DB_" + fileWP + "wp_slim/Zbby_" + fileWP + "wp_slim.root"         
        };
        //outputDir = "../Plots/DB/" + fileWP + "/MassResponse/";
        outputDir = "../Plots/DB/MassResponse/" + fileWP + "/";
    }
}; //kinit tamezza@CERN.CH



/**
 * Creates directory if it does not exist
 * @param path Directory path to create
 */
void CreateDirectoryIfNotExists(const std::string& path) {
    std::string command = "mkdir -p " + path;
    system(command.c_str());
}

/**
 * Main class for histogram management and analysis
 */
class HistogramManager {
private:
    std::vector<TH1F*> h_ljet_m;
    std::vector<TH1F*> h_ljet_m_50_350;
    std::vector<TH1F*> h_ljet_m_50_150;
    std::vector<TH1F*> h_ljet_m_50_150_2;
    std::vector<TH1F*> h_Asym;
    std::vector<TH1F*> h_Asym_2;
    std::vector<TH1F*> h_R_DB_1;
    std::vector<TH1F*> h_R_DB_2;
    std::vector<TH1F*> h_R_DB_3;
    std::vector<TH1F*> h_R_DB;

    // Background composition histograms
    std::unique_ptr<TH1D> h_background_composition;
    std::vector<std::unique_ptr<TH1F>> h_R_DB_bins;
    std::vector<std::unique_ptr<TH1F>> h_mass_bins;

    std::string workingPoint;   // Current working point being analyzed
    std::string outputDir;      // Output directory for plots

    double  p0;  // Scale factor intercept parameter
    double  p1;  // Scale factor slope parameter
    //std::vector<double> fitParameters;  // Store fit parameters [intercept, slope]
    //std::vector<double> fitErrors;      // Store fit parameter errors



    /**
     * Initialize histograms for a specific variable
     * @param histograms Vector to store the created histograms
     * @param variable Variable name for histogram naming
     * @param title Title of the histogram
     * @param nbins Number of bins
     * @param min Minimum x-axis value
     * @param max Maximum x-axis value
     */
    void InitializeHistogramGroup(std::vector<TH1F*>& histograms, const std::string& variable, 
                                const std::string& title, int nbins, float min, float max) {
        // Create histograms for data, MultiJets, and Z->bb jets
        std::vector<std::string> names = { "data", "gamma_jets", "Zbby" };
        for (const auto& name : names) {
            histograms.push_back(new TH1F(("h_" + name + "_" + variable).c_str(), 
                                        (title + ";Entries").c_str(), nbins, min, max));
        }
    }

    /**
     * Sets the color scheme for the histograms
     * @param histograms Vector of histograms to color
     */
    void SetHistogramColors(std::vector<TH1F*>& histograms) {
        std::vector<int> colors = { kBlack, kBlue, kRed };  // Data, MultiJets, Z->bb
        for (size_t i = 0; i < histograms.size(); ++i) {
            histograms[i]->SetLineColor(colors[i]);
            histograms[i]->SetFillColor(colors[i]);
        }
    }

public:
    /**
     * Constructor
     * @param wp Working point identifier
     */
    HistogramManager(const std::string& wp, const std::string& dir) : workingPoint(wp), outputDir(dir), p0(1.0), p1(0.0)
        //fitParameters(2, 0.0),  // Initialize with 2 zeros
        //fitParameters{1.0, 0.0}, // Initialize with 1.0 (intercept) and 0.0 (slope)
        //fitErrors(2, 0.0)       // Initialize with 2 zeros
    {
        InitializeHistograms();
    }

    /**
     * Destructor - clean up all dynamically allocated histograms
     */
    ~HistogramManager() {
        // Cleanup all histograms
        for (auto hist : h_ljet_m) delete hist;
        for (auto hist : h_ljet_m_50_350) delete hist;
        for (auto hist : h_ljet_m_50_150) delete hist;
        for (auto hist : h_ljet_m_50_150_2) delete hist;
        for (auto hist : h_Asym) delete hist;
        for (auto hist : h_Asym_2) delete hist;
        for (auto hist : h_R_DB_1) delete hist;
        for (auto hist : h_R_DB_2) delete hist;
        for (auto hist : h_R_DB_3) delete hist;
        for (auto hist : h_R_DB) delete hist;
        // ... cleanup other histogram vectors
    }

    /**
     * Initialize all histograms used in the analysis
     */
    void InitializeHistograms() {
        InitializeHistogramGroup(h_ljet_m, "ljet_m", "Large-R Jet Mass;M^{J} [GeV];Events", 50, 50, 200);
        InitializeHistogramGroup(h_ljet_m_50_350, "ljet_m_50_350", "Large-R Jet Mass;M^{J} [GeV];Events", 50, 50, 350);
        InitializeHistogramGroup(h_ljet_m_50_150, "ljet_m_50_150", "Large-R Jet Mass;M^{J} [GeV];Events", 20, 50, 150);
        InitializeHistogramGroup(h_ljet_m_50_150_2, "ljet_m_50_150_2", "Large-R Jet Mass;M^{J} [GeV];Events", 50, 50, 150);
        InitializeHistogramGroup(h_Asym, "Asym", "Asymmetry;Asymmetry;Events", 50, -1, 1);
        InitializeHistogramGroup(h_Asym_2, "Asym_2", "Asymmetry;Asymmetry;Events", 20, -0.35, 0.35);
        InitializeHistogramGroup(h_R_DB_1, "R_DB_1", "The Direct Balance Ratio;R_{MJB};Events", 50, 0, 2);
        InitializeHistogramGroup(h_R_DB_2, "R_DB_2", "The Direct Balance Ratio;R_{MJB};Events", 55, 0, 2.5);
        InitializeHistogramGroup(h_R_DB_3, "R_DB_3", "The Direct Balance Ratio;R_{MJB};Events", 55, 0, 2);
        InitializeHistogramGroup(h_R_DB, "R_DB", "The Direct Balance Ratio;R_{MJB};Events", 60, 0, 100);

        // Set colors scheme for all histogram groups
        SetHistogramColors(h_ljet_m);
        SetHistogramColors(h_ljet_m_50_350);
        SetHistogramColors(h_ljet_m_50_150);
        SetHistogramColors(h_ljet_m_50_150_2);
        SetHistogramColors(h_Asym);
        SetHistogramColors(h_Asym_2);
        SetHistogramColors(h_R_DB_1);
        SetHistogramColors(h_R_DB_2);
        SetHistogramColors(h_R_DB_3);
        SetHistogramColors(h_R_DB);

        // Initialize background composition histogram for flavor categories
        h_background_composition = std::make_unique<TH1D>("h_background_composition", 
            "Background Composition;Category;Events", 6, 0, 6);

        // Initialize binned histograms
        for (int i = 0; i < 6; ++i) {
            h_R_DB_bins.push_back(std::make_unique<TH1F>(
                ("h_R_DB_bin_" + std::to_string(i)).c_str(),
                "R_DB;R_DB;Events", 60, 0, 3.5));
            h_mass_bins.push_back(std::make_unique<TH1F>(
                ("h_mass_bin_" + std::to_string(i)).c_str(),
                "Mass;Mass [GeV];Events", 50, 50, 500));
        }
    }

    /**
     * Setup branch addresses for reading from the TTree
     * @param tree The TTree to read from
     * @param ljet_m Output vector for jet mass
     * @param weight Output for event weight
     * @param Asym Output for asymmetry value
     * @param R_DB Output for direct balance ratio
     * @param ljet_flavour_label Output for jet flavor labels
     * @param ljet_truth_label Output for jet truth labels
     */
    void SetBranchAddresses(TTree* tree, std::vector<float>*& ljet_m, double& weight, 
                           double& Asym, double& R_DB, 
                           std::vector<int>*& ljet_flavour_label, 
                           std::vector<int>*& ljet_truth_label) {
        // Set addresses for basic variables
        tree->SetBranchAddress("ljet_m", &ljet_m);
        tree->SetBranchAddress("total_weight", &weight);
        tree->SetBranchAddress("Asym", &Asym);
        tree->SetBranchAddress("R_DB", &R_DB);

        // Check if flavor labeling branches exist (may not be in data)
        if (!tree->GetBranch("ljet_flavour_label") || !tree->GetBranch("ljet_truth_label")) 
            return;
        
        // Set addresses for flavor information
        tree->SetBranchAddress("ljet_flavour_label", &ljet_flavour_label);
        tree->SetBranchAddress("ljet_truth_label", &ljet_truth_label);
    }

    /**
     * Fill histograms from a TTree
     * @param tree The TTree to read from
     * @param histograms Vector of histograms to fill
     * @param isData Flag indicating if the tree contains data (not MC)
     */
    void FillHistograms(TTree* tree, std::vector<TH1F*>& histograms, bool isData) {
        // Variables to read from the tree
        std::vector<float>* ljet_m = nullptr;
        double weight = 0.0;
        double R_DB = 0.0;
        double Asym = 0.0;
        std::vector<int>* ljet_flavour_label = nullptr;
        std::vector<int>* ljet_truth_label = nullptr;

        // Setup the branch addresses
        SetBranchAddresses(tree, ljet_m, weight, Asym, R_DB, 
                          ljet_flavour_label, ljet_truth_label);
        
        // Loop over all entries in the tree
        Long64_t nentries = tree->GetEntries();
        for (Long64_t i = 0; i < nentries; ++i) {
            tree->GetEntry(i);
            double eventWeight = isData ? 1.0 : weight; // Data has weight=1

            // Fill all histograms with appropriate variables
            histograms[0]->Fill((ljet_m->at(0))/1000.0, eventWeight);
            histograms[1]->Fill((ljet_m->at(0))/1000.0, eventWeight);
            histograms[2]->Fill((ljet_m->at(0))/1000.0, eventWeight);
            histograms[3]->Fill((ljet_m->at(0))/1000.0, eventWeight);
            histograms[4]->Fill(Asym, eventWeight);
            histograms[5]->Fill(Asym, eventWeight);
            histograms[6]->Fill(R_DB, eventWeight);
            histograms[7]->Fill(R_DB, eventWeight);
            histograms[8]->Fill(R_DB, eventWeight);
            histograms[9]->Fill(R_DB, eventWeight);

            // For MC samples, fill flavor composition info
            if (!isData && ljet_flavour_label && ljet_truth_label) {
                // Count number of b and c quarks in the jet
                int number_of_b = 0, number_of_c = 0;
                for (const auto& label : *ljet_flavour_label) {
                    if (label == 55) number_of_b += 2;
                    else if (label == 54) { number_of_b += 1; number_of_c += 1; }
                    else if (label == 51) number_of_b += 1;
                    else if (number_of_b == 0) {
                        if (label == 44) number_of_c += 2;
                        else if (label == 41) number_of_c += 1;
                        else if (label == 11) number_of_c = 0;
                    }
                }

                // Determine jet flavor category
                int jet_bin = (number_of_b >= 2) ? 0 :              // >=2 b-quarks
                    (number_of_b == 1 && number_of_c >= 1) ? 1 :    // 1 b-quark, >=1 c-quark
                    (number_of_b == 1 && number_of_c == 0) ? 2 :    // 1 b-quark, 0 c-quarks
                    (number_of_b == 0 && number_of_c >= 2) ? 3 :    // 0 b-quarks, >=2 c-quarks
                    (number_of_b == 0 && number_of_c == 1) ? 4 : 5; // 0 b-quarks, 1 c-quark or light jets

                // Fill jet flavor category histograms
                h_background_composition->Fill(jet_bin);
                if (jet_bin < h_R_DB_bins.size()) {
                    h_R_DB_bins[jet_bin]->Fill(R_DB, weight);
                    h_mass_bins[jet_bin]->Fill((ljet_m->at(0))/1000.0, weight);
                }
            }
        }

        if (isData) {
            std::cout << "Data histogram entries: " << histograms[0]->GetEntries() << std::endl;
        }
    }

    /**
     * Process all input files and fill histograms
     * @param files List of input ROOT files
     */
    void ProcessFiles(const std::vector<std::string>& files) {
        //std::vector<std::unique_ptr<TFile>> rootFiles;
         //std::vector<TTree*> trees;

         std::vector<TFile*> rootFiles;
         std::vector<TTree*> trees;
 
         // Open all files and get their trees
         for (const auto& file : files) {
             std::cout << "Opening file: " << file << std::endl;
             // Use TFile::Open directly for remote file access
             TFile* rootFile = TFile::Open(file.c_str(), "READ");

             //auto rootFile = std::make_unique<TFile>(file.c_str(), "READ");

             //if (!rootFile || rootFile->IsZombie()) {
             //   std::cerr << "Error: Could not open file " << file << std::endl;
             //   continue;
            //}
            
            if (rootFile && !rootFile->IsZombie()) {
                TTree* tree = (TTree*)rootFile->Get("nominal");
                if (tree) {
                    trees.push_back(tree);
                    rootFiles.push_back(std::move(rootFile));
                }
            }
        }

        // Process each file/tree
        for (size_t i = 0; i < trees.size(); ++i) {
            bool isData = (i == 0); // First file is assumed to be data
            // Group histograms for filling
            std::vector<TH1F*> histograms = {
                h_ljet_m[i], h_ljet_m_50_350[i], h_ljet_m_50_150[i], 
                h_ljet_m_50_150_2[i], h_Asym[i], h_Asym_2[i], 
                h_R_DB_1[i], h_R_DB_2[i], h_R_DB_3[i], h_R_DB[i]
            };
            FillHistograms(trees[i], histograms, isData);
        }

        // Calculate and apply scale factor from sidebands
        CalculateAndApplyScaleFactor();
    }

    /**
     * Calculate data/MC scale factor from sideband regions and apply it to MC
     */
    void CalculateAndApplyScaleFactor() {
        // Define sideband regions (avoiding signal region)
        int binLow1 = h_ljet_m_50_150[1]->FindBin(50);
        int binHigh1 = h_ljet_m_50_150[1]->FindBin(65);
        int binLow2 = h_ljet_m_50_150[1]->FindBin(110);
        int binHigh2 = h_ljet_m_50_150[1]->FindBin(150);
        
        // Create arrays to store x and y values for the fit
        std::vector<double> x_values;
        std::vector<double> y_values;
        std::vector<double> y_errors;
        
        // Extract bin contents for first sideband region
        for (int bin = binLow1; bin <= binHigh1; ++bin) {
            double bin_center = h_ljet_m_50_150[1]->GetBinCenter(bin);
            double data_content = h_ljet_m_50_150[0]->GetBinContent(bin);
            double mc_content = h_ljet_m_50_150[1]->GetBinContent(bin);
            
            // Avoid division by zero
            if (mc_content > 0) {
                x_values.push_back(bin_center);
                y_values.push_back(data_content / mc_content);
                
                // Calculate error on the ratio
                double data_error = h_ljet_m_50_150[0]->GetBinError(bin);
                double mc_error = h_ljet_m_50_150[1]->GetBinError(bin);
                double ratio_error = sqrt(pow(data_error/mc_content, 2) + 
                                         pow(data_content*mc_error/(mc_content*mc_content), 2));
                y_errors.push_back(ratio_error);
            }
        }
        
        // Extract bin contents for second sideband region
        for (int bin = binLow2; bin <= binHigh2; ++bin) {
            double bin_center = h_ljet_m_50_150[1]->GetBinCenter(bin);
            double data_content = h_ljet_m_50_150[0]->GetBinContent(bin);
            double mc_content = h_ljet_m_50_150[1]->GetBinContent(bin);
            
            // Avoid division by zero
            if (mc_content > 0) {
                x_values.push_back(bin_center);
                y_values.push_back(data_content / mc_content);
                
                // Calculate error on the ratio
                double data_error = h_ljet_m_50_150[0]->GetBinError(bin);
                double mc_error = h_ljet_m_50_150[1]->GetBinError(bin);
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
        //double p0 = fit_func->GetParameter(0);  // Intercept
        //double p1 = fit_func->GetParameter(1);  // Slope
        p0 = fit_func->GetParameter(0);  // Instead of double p0 = ...
        p1 = fit_func->GetParameter(1);  // Instead of double p1 = ...
        double p0_err = fit_func->GetParError(0);
        double p1_err = fit_func->GetParError(1);
        
        printf("Linear fit results: SF(m) = %.4f (± %.4f) + %.8f (± %.8f) * m\n", 
               p0, p0_err, p1, p1_err);
        
        // Save the fit result to a canvas
        TCanvas* fit_canvas = new TCanvas("fit_canvas", "Scale Factor Fit", 800, 600);
        ratio_graph->SetTitle("Data/MC Ratio in Sidebands;Large-R Jet Mass [MeV];Data/MC");
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
        
        DrawAtlasLabel();
        fit_canvas->SaveAs((this->outputDir + "scale_factor_fit.pdf").c_str());
        
        // Apply the mass-dependent scale factor to multijet histograms
        std::vector<TH1F*> multijetHistograms = {
            h_ljet_m[1], h_ljet_m_50_350[1], h_ljet_m_50_150[1],
            h_ljet_m_50_150_2[1]
        };
        
        // First scale mass histograms with mass-dependent scale factor
        for (TH1F* hist : multijetHistograms) {
            if (hist) {
                for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
                    double bin_center = hist->GetBinCenter(bin);
                    double content = hist->GetBinContent(bin);
                    double error = hist->GetBinError(bin);
                    
                    // Calculate the bin-by-bin scale factor
                    double scale_factor = p0 + p1 * bin_center;
                    //double scale_factor = 1;
                    
                    // Apply the scale factor
                    hist->SetBinContent(bin, content * scale_factor);
                    hist->SetBinError(bin, error * scale_factor);
                }
            }
        }
        
        // Calculate average scale factor for histograms that do not have mass on x-axis
        double avgScaleFactor = 0.0;
        int nBins = 0;
        
        // Calculate average SF in the mass range of interest (50-150 GeV)
        for (int bin = binLow1; bin <= binHigh2; ++bin) {
            double bin_center = h_ljet_m_50_150[1]->GetBinCenter(bin);
            double scale_factor = p0 + p1 * bin_center;
            //double scale_factor = 1;
            avgScaleFactor += scale_factor;
            nBins++;
        }
        
        avgScaleFactor /= nBins;
        printf("Average Scale Factor: %.4f\n", avgScaleFactor);
        
        // Apply the average scale factor to non-mass histograms
        std::vector<TH1F*> otherHistograms = {
            h_Asym[1], h_Asym_2[1],
            h_R_DB_1[1], h_R_DB_2[1], h_R_DB_3[1], h_R_DB[1]
        };
        
        for (TH1F* hist : otherHistograms) {
            if (hist) {
                hist->Scale(avgScaleFactor);
            }
        }
        
        // Clean up
        delete fit_func;
        delete ratio_graph;
        delete fit_canvas;
    }

    /**
     * Save all histograms to output files and create plots
     * @param outputDir Directory to save the output files
     * @param createRatioPlots Flag to create ratio panels in plots
     */
    void SaveHistograms(bool createRatioPlots = true) {
        CreateDirectoryIfNotExists(outputDir);

        // List of canvas names for saving
        std::vector<std::string> canvasNames = {
            "Large_R_Jet_m", "Large_R_Jet_m_50_350", "Large_R_Jet_m_50_150",
            "Large_R_Jet_m_50_150_2", "Asymmetry", "Asymmetry_2",
            "R_DB_1", "R_DB_2", "R_DB_3", "R_DB"
        };

        // Group histograms for saving
        std::vector<std::vector<TH1F*>> histogramGroups = {
            h_ljet_m, h_ljet_m_50_350, h_ljet_m_50_150, h_ljet_m_50_150_2,
            h_Asym, h_Asym_2, h_R_DB_1, h_R_DB_2, h_R_DB_3, h_R_DB
        };

        // Save each histogram group
        for (size_t i = 0; i < canvasNames.size(); ++i) {
            SaveHistogramGroup(histogramGroups[i], canvasNames[i], outputDir, createRatioPlots);
        }

        // Save background composition
        TFile outputFile((outputDir + "background_composition.root").c_str(), "RECREATE");
        h_background_composition->Write();
        for (const auto& hist : h_R_DB_bins) {
            hist->Write();
        }
        for (const auto& hist : h_mass_bins) {
            hist->Write();
        }
        outputFile.Close();

        // Create and save additional plots
        SaveBackgroundComposition(outputDir);
        SaveSplitBackgroundHistograms(outputDir);
    }

private:
    /**
     * Save a group of histograms to a canvas with optional ratio plot
     * @param histograms Vector of histograms to save
     * @param name Name for the canvas and output file
     * @param outputDir Directory to save the output file
     * @param createRatioPlot Flag to create a ratio panel in the plot
     */
    void SaveHistogramGroup(std::vector<TH1F*>& histograms, const std::string& name, 
                          const std::string& outputDir, bool createRatioPlot) {
        TCanvas* canvas = new TCanvas(("canvas_" + name).c_str(), name.c_str(), 800, 800);
        
        if (createRatioPlot) {
            // Create upper plot pad (70% of canvas)
            TPad* pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
            pad1->SetBottomMargin(0.015);
            if (name.find("_log") != std::string::npos) {
                pad1->SetLogy(true);
            }
            pad1->Draw();
            pad1->cd();

            // THStack* hs = new THStack("hs", (name + ";Events").c_str());
            // Create stack and draw histograms
            THStack* hs = new THStack("hs", "");  // Empty title, we will set it manually
            for (size_t i = 1; i < histograms.size(); ++i) {
                hs->Add(histograms[i]);
            }
            hs->SetMaximum(hs->GetMaximum() * 1.5);

            // Set proper axis titles
            hs->Draw("HIST");
            hs->GetXaxis()->SetTitle(histograms[0]->GetXaxis()->GetTitle());
            hs->GetYaxis()->SetTitle(histograms[0]->GetYaxis()->GetTitle());
            hs->GetXaxis()->SetTitleSize(0.05);
            hs->GetYaxis()->SetTitleSize(0.05);
            hs->GetXaxis()->SetTitleOffset(1.2);
            hs->GetYaxis()->SetTitleOffset(1.2);
            hs->GetXaxis()->SetLabelSize(0);  // Hide x-axis labels for upper pad
            hs->SetMaximum(hs->GetMaximum() * 1.5);
            
            histograms[0]->Draw("EP same");

            //hs->GetXaxis()->SetLabelSize(0);
            //histograms[0]->GetXaxis()->SetLabelSize(0);

            // Calculate and draw MC statistical uncertainty band
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

            // Draw legend
            TLegend* legend = new TLegend(0.75, 0.6, 0.9, 0.93);
            legend->SetBorderSize(0);
            legend->SetTextSize(0.03);
            legend->AddEntry(histograms[0], "Data", "lep");
            legend->AddEntry(histograms[2], "Z (#rightarrow bb) + #gamma", "f");
            legend->AddEntry(histograms[1], "#gamma + jets", "f");
            legend->AddEntry(err_band, "MC Stat. Unc.", "f");
            legend->Draw();

            DrawAtlasLabel();

            // Ratio plot pad
            canvas->cd();
            TPad* pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
            pad2->SetTopMargin(0.05);
            pad2->SetBottomMargin(0.35);
            pad2->Draw();
            pad2->cd();

            // Create and draw ratio plot
            TH1F* h_ratio = (TH1F*)histograms[0]->Clone();
            h_ratio->SetLineColor(kBlack);
            h_ratio->SetMinimum(0.5);
            h_ratio->SetMaximum(1.5);
            h_ratio->Divide(h_mc_total);
            //h_ratio->Draw("EP");
            h_ratio->SetMarkerStyle(21);
            // Use the original histograms axis titles
            h_ratio->GetXaxis()->SetTitle(histograms[0]->GetXaxis()->GetTitle());
            //h_ratio->GetXaxis()->SetTitle(name.c_str());
            h_ratio->GetXaxis()->SetTitleSize(0.12);
            h_ratio->GetXaxis()->SetLabelSize(0.12);
            h_ratio->GetXaxis()->SetTitleOffset(1.25);
            h_ratio->GetYaxis()->SetTitle("Data/MC");
            h_ratio->GetYaxis()->SetTitleSize(0.15);
            h_ratio->GetYaxis()->SetLabelSize(0.12);
            h_ratio->GetYaxis()->SetTitleOffset(0.3);
            h_ratio->GetYaxis()->SetNdivisions(505);
            h_ratio->Draw("El same");

            // Add the scale factor fit to the ratio plot if this is a mass histogram 
            /*
            if (name.find("Large_R_Jet_m") != std::string::npos) {
                // Create a function for the scale factor in the ratio pad
                //TF1* sf_line = new TF1("sf_line", "[0] + [1] * x", h_ratio->GetXaxis()->GetXmin(), h_ratio->GetXaxis()->GetXmax());
                double xmin = h_ratio->GetXaxis()->GetXmin();
                double xmax = h_ratio->GetXaxis()->GetXmax();
                TF1* sf_line = new TF1("sf_line", "[0] + [1]*x", xmin, xmax);

                sf_line->SetParameters(p0, p1);
                sf_line->SetLineColor(kBlue);
                //sf_line->SetLineWidth(2);
                sf_line->SetLineStyle(2);
                sf_line->Draw("SAME");
    
                TLatex latex;
                latex.SetTextSize(0.12);
                latex.SetTextColor(kBlue);
                latex.DrawLatexNDC(0.175, 0.85, Form("Fit: %.2f + %.3fx", p0, p1));
                
            }
            */
    

            // Draw ratio uncertainty band
            TGraphAsymmErrors* ratio_err_band = new TGraphAsymmErrors(h_mc_total);
            for (int i = 0; i < ratio_err_band->GetN(); ++i) {
                ratio_err_band->SetPointY(i, 1);
                //double rel_error = h_mc_total->GetBinError(i + 1) / h_mc_total->GetBinContent(i + 1);
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
            //line->Draw();
            line->Draw("same");

            pad2->Update();
            canvas->Update();
        }

        canvas->SaveAs((outputDir + name + ".pdf").c_str());
        delete canvas;
    }

    /**
     * Save background composition histogram
     * @param outputDir Directory to save the output file
     */
    void SaveBackgroundComposition(const std::string& outputDir) {
        TCanvas* canvas = new TCanvas("background_composition", "Background Composition", 800, 600);
        h_background_composition->Draw();

        // Add labels for each bin
        h_background_composition->GetXaxis()->SetBinLabel(1, ">=2 b");
        h_background_composition->GetXaxis()->SetBinLabel(2, "1b, >=1c");
        h_background_composition->GetXaxis()->SetBinLabel(3, "1b, 0c");
        h_background_composition->GetXaxis()->SetBinLabel(4, "0b, >=2c");
        h_background_composition->GetXaxis()->SetBinLabel(5, "0b, 1c");
        h_background_composition->GetXaxis()->SetBinLabel(6, "Light");

        DrawAtlasLabel();
        canvas->SaveAs((outputDir + "background_composition.pdf").c_str());
        delete canvas;
    }

    /**
     * Save histograms split by background type
     * @param outputDir Directory to save the output files
     */
    void SaveSplitBackgroundHistograms(const std::string& outputDir) {
        // Save R_DB split histograms
        TCanvas* canvas_rdb = new TCanvas("R_DB_split", "R_DB Split", 800, 600);
        DrawSplitHistograms(h_R_DB_bins, canvas_rdb);
        canvas_rdb->SaveAs((outputDir + "R_DB_split.pdf").c_str());
        delete canvas_rdb;

        // Save mass split histograms
        TCanvas* canvas_mass = new TCanvas("mass_split", "Mass Split", 800, 600);
        DrawSplitHistograms(h_mass_bins, canvas_mass);
        canvas_mass->SaveAs((outputDir + "mass_split.pdf").c_str());
        delete canvas_mass;
    }

    /**
     * Draw a set of histograms on one canvas with different colors
     * @param histograms Vector of histograms to draw
     * @param canvas Canvas to draw on
     */
    void DrawSplitHistograms(const std::vector<std::unique_ptr<TH1F>>& histograms, TCanvas* canvas) {
        // Create legend
        TLegend* legend = new TLegend(0.75, 0.6, 0.9, 0.93);
        legend->SetBorderSize(0);
        legend->SetTextSize(0.03);

        // Find maximum y-value across all histograms
        double max = 0;
        for (const auto& hist : histograms) {
            if (hist->GetMaximum() > max) {
                max = hist->GetMaximum();
            }
        }
        max *= 1.75;

        // Names for the legend
        std::vector<std::string> binNames = {
            ">=2 b-quarks", 
            "1 b-quark, >= c-quark", 
            "1 b-quark, 0 c-quarks",
            "0 b-quarks, >= c-quarks",
            "0 b-quarks, 1 c-quark",
            "Light jets"
        };

        // Draw all histograms with different colors
        for (size_t i = 0; i < histograms.size(); ++i) {
            histograms[i]->SetMaximum(max);
            histograms[i]->SetLineWidth(2);
            histograms[i]->SetLineColor(i + 1);
            histograms[i]->Draw(i == 0 ? "HIST" : "HIST SAME");
            legend->AddEntry(histograms[i].get(), ("Bin " + std::to_string(i)).c_str(), "l");
        }

        legend->Draw();
        DrawAtlasLabel();
    }
};

/**
 * Main function to run the analysis
 */
void Mass_Response() {
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();

    // Define working points to process
    //std::vector<std::string> workingPoints = {"50", "55", "60", "65"};
    //std::vector<std::string> workingPoints = {"25", "30", "37", "46", "58", "74", "94", "125", "155"};
    std::vector<std::string> workingPoints = {"25", "46", "74"};

    // Process each working point
    for (const auto& wp : workingPoints) {
        std::cout << "Processing working point " << wp << std::endl;
        
        // Configure working point
        WorkingPoint workingPoint(wp);
        CreateDirectoryIfNotExists(workingPoint.outputDir);

        // Create histogram manager and process data
        HistogramManager histManager(wp, workingPoint.outputDir); // Pass outputDir to construc
        histManager.ProcessFiles(workingPoint.inputFiles);
        //histManager.SaveHistograms(workingPoint.outputDir);
        histManager.SaveHistograms(); // No parameters needed now
        

        std::cout << "Completed processing working point " << wp << std::endl;
    }
}

int main() {
    // Disable ROOT global graphics   
    gROOT->SetBatch(kTRUE);

    // Run the analysis
    Mass_Response();

    // Cleanup ROOT
    gROOT->GetListOfCanvases()->Delete();
    gDirectory->DeleteAll();

    // Force application termination
    gSystem->Exit(0);

    return 0;
}