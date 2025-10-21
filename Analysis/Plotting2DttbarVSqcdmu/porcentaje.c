#include <TChain.h>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>

// Funciones auxiliares
std::vector<std::pair<float, float>> getCuts() {
    std::vector<std::pair<float, float>> cuts;
    for (int i = 1; i <= 12; ++i) {
        cuts.emplace_back(0.05f * i, 5.f * i);
    }
    return cuts;
}

std::vector<std::string> getCutLabels(const std::vector<std::pair<float, float>>& cuts) {
    std::vector<std::string> labels;
    for (const auto& cut : cuts) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(2) << cut.first << "/" << static_cast<int>(cut.second);
        labels.push_back(oss.str());
    }
    return labels;
}

std::vector<float> porcentajeEnVentana(TChain* chain, const std::string& xvar, const std::string& yvar, const std::vector<std::pair<float, float>>& cuts) {
    std::vector<float> porcentajes(cuts.size(), 0.0f);

    float xval = 0.0f;
    float yval = 0.0f;

    chain->SetBranchAddress(xvar.c_str(), &xval);
    chain->SetBranchAddress(yvar.c_str(), &yval);

    Long64_t nEntries = chain->GetEntries();
    std::vector<int> dentro(cuts.size(), 0);

    for (Long64_t i = 0; i < nEntries; ++i) {
        chain->GetEntry(i);
        for (size_t j = 0; j < cuts.size(); ++j) {
            if (xval < cuts[j].first && yval < cuts[j].second) {
                dentro[j]++;
            }
        }
    }

    for (size_t j = 0; j < cuts.size(); ++j) {
        porcentajes[j] = (nEntries > 0) ? 100.f * dentro[j] / nEntries : 0.f;
    }

    return porcentajes;
}

// Función principal ROOT-interpretable
void porcentaje() {
    TChain ttbar("Events");
    TChain qcd("Events");

    ttbar.Add("/eos/home-c/cmunozdi/AnalysisSamples/2dQCDcuts/TTtoLNu2Q/*.root");
    qcd.Add("/eos/home-c/cmunozdi/AnalysisSamples/2dQCDcuts/QCDMu/**/*.root");

    std::vector<int> thresholds = {15, 20, 25, 30};
    auto cuts = getCuts();
    auto cutLabels = getCutLabels(cuts);

    std::cout << std::setw(18) << "Cut dR/ptRel";
    for (int thr : thresholds) {
        std::cout << std::setw(25) << ("QCDMu (thr " + std::to_string(thr) + " GeV)");
        std::cout << std::setw(25) << ("TTbar (thr " + std::to_string(thr) + " GeV)");
    }
    std::cout << "\n";

    for (size_t i = 0; i < cuts.size(); ++i) {
        std::cout << std::setw(18) << cutLabels[i];

        for (int thr : thresholds) {
            std::string xvar = "lepton_" + std::to_string(thr) + ".dR_to_jet";
            std::string yvar = "lepton_" + std::to_string(thr) + ".pt_rel_to_jet";

            auto qcd_pct = porcentajeEnVentana(&qcd, xvar, yvar, cuts);
            auto ttbar_pct = porcentajeEnVentana(&ttbar, xvar, yvar, cuts);

            std::cout << std::setw(24) << std::fixed << std::setprecision(2) << qcd_pct[i] << "%";
            std::cout << std::setw(24) << std::fixed << std::setprecision(2) << ttbar_pct[i] << "%";
        }
        std::cout << "\n";
    }
}
