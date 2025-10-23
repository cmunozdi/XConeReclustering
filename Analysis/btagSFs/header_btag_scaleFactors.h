#ifndef __nanoPF_helpers_h__
#define __nanoPF_helpers_h__

#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"
#include <vector>
#include "TVector3.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include <TObject.h>


#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include <cassert>
#include <filesystem>
#include <iostream>
#include <string>
#include "correction.h"

using namespace ROOT::VecOps;
using namespace std;
namespace fs = std::filesystem;

using rvec_f = const RVec<float>;
using rvec_i = const RVec<int>;
using rvec_b = const RVec<bool>;

// Global variable for correction files
std::unique_ptr<correction::CorrectionSet> sf_btagset;
std::unique_ptr<correction::CorrectionSet> eff_btagset;
std::unique_ptr<correction::CorrectionSet> sf_muoset;
std::unique_ptr<correction::CorrectionSet> sf_purew;

// Initialize `cset` at the begining of the header
inline void initializeBTagCorrectionSet() {
    fs::path fname_sf_btag_2023 = "./btagSFs/btagging.json.gz"; //"/eos/user/c/cmunozdi/SWAN_projects/BtagScaleFactors/btagging.json.gz"; 
    fs::path fname_eff_btag_2023 = "./btagSFs/btag_efficiencies_combined.json"; //"/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightIDNoLepVeto/output_efficiencies_rdf/btag_efficiencies_combined.json"; 
    // fs::path fname_btag_2023 = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2023_Summer23/btagging.json.gz";
    if (!fs::exists(fname_sf_btag_2023)) {
        throw std::runtime_error("Correction file does not exist: " + fname_sf_btag_2023.string());
    }
    sf_btagset = correction::CorrectionSet::from_file(fname_sf_btag_2023.string());
    std::cout << "Correction file loaded susscessfully: " << fname_sf_btag_2023.string() << std::endl;

    if (!fs::exists(fname_eff_btag_2023)) {
        throw std::runtime_error("Correction file does not exist: " + fname_eff_btag_2023.string());
    }
    eff_btagset = correction::CorrectionSet::from_file(fname_eff_btag_2023.string());
    std::cout << "Correction file loaded susscessfully: " << fname_eff_btag_2023.string() << std::endl;
}

inline void initializeMUOCorrectionSet() {
    fs::path fname_sf_MUO_2023 = "./muoSFs/muon_HighPt.json.gz";
    if (!fs::exists(fname_sf_MUO_2023)) {
        throw std::runtime_error("Correction file does not exist: " + fname_sf_MUO_2023.string());
    }
    sf_muoset = correction::CorrectionSet::from_file(fname_sf_MUO_2023.string());
    std::cout << "Correction file loaded susscessfully: " << fname_sf_MUO_2023.string() << std::endl;
}

inline void initializePUReweightingCorrectionSet() {
    fs::path fname_puReweighting_2023 = "./puRew/puWeights.json.gz";
    if (!fs::exists(fname_puReweighting_2023)) {
        throw std::runtime_error("Correction file does not exist: " + fname_puReweighting_2023.string());
    }
    sf_purew = correction::CorrectionSet::from_file(fname_puReweighting_2023.string());
    std::cout << "Correction file loaded susscessfully: " << fname_puReweighting_2023.string() << std::endl;
}

// Funtion to compute btagWeight of the event given the jet collection
inline float compute_btagWeight(const rvec_f &jet_pts, const rvec_f &jet_etas, const ROOT::VecOps::RVec<UChar_t>& jet_flavs, const rvec_f &jet_btags) {
    float weight = 1.;
    float jet_i_eta = 0.;
    for (size_t i = 0; i < jet_pts.size(); i++) {

        jet_i_eta = jet_etas[i];
        if (jet_pts[i] < 30. || abs(jet_etas[i]) > 2.5 /*|| jet_pts[i] > 1000*/) continue;

        if(abs(jet_etas[i])==2.5) jet_i_eta=2.499; // to avoid the out of range error in the btag SF evaluation

        float sf = 1.;
        try {
            if (static_cast<int>(jet_flavs[i]) == 0) {
                sf = sf_btagset->at("deepJet_light")->evaluate(
                    std::vector<correction::Variable::Type>{"central", "T", static_cast<int>(jet_flavs[i]), abs(jet_i_eta), jet_pts[i]}
                );
            } else {
                sf = sf_btagset->at("deepJet_comb")->evaluate(
                    std::vector<correction::Variable::Type>{"central", "T", static_cast<int>(jet_flavs[i]), abs(jet_i_eta), jet_pts[i]}
                );
            }
        } catch (const std::exception &e) {
            
            std::cout << "Jet " << i << ": pt = " << jet_pts[i] << ", eta = " << jet_etas[i] 
                  << ", flav = " << static_cast<int>(jet_flavs[i]) << ", btag = " << jet_btags[i] << std::endl;
            std::cerr << "Error evaluationg the SF for jet " << i << ": " << e.what() << ". Value of sf=" << sf << std::endl;
            continue;
        }

        float eff = 1.;
        try {
            eff = eff_btagset->at("btag_efficiency")->evaluate(
                std::vector<correction::Variable::Type>{jet_i_eta, jet_pts[i], static_cast<int>(jet_flavs[i])}
            );
        } catch (const std::exception &e) {
            
            std::cout << "Jet " << i << ": pt = " << jet_pts[i] << ", eta = " << jet_etas[i] 
                  << ", flav = " << static_cast<int>(jet_flavs[i]) << ", btag = " << jet_btags[i] << std::endl;
            std::cerr << "Error evaluating efficiency for jet " << i << ": " << e.what() << ". Value of eff=" << eff << std::endl;
            continue;
        }

        float jet_w = 1;
        if (jet_btags[i] > 0.6553) {
            jet_w = sf;
        } else {
            jet_w = (1 - sf * eff) / (1 - eff);
        }

        weight *= jet_w;
    }
    return weight;
}

// Function to compute MUO weight of the event with just one muon
inline float compute_MUOWeight(const float &muo_pt, const float &muo_eta){
    float weight =1.;
    float muo_p = muo_pt*std::cosh(muo_eta);

    float w_global = 1., w_highptid = 1., w_iso = 1., w_hlt = 1.;

    try {
        w_global = sf_muoset->at("NUM_GlobalMuons_DEN_TrackerMuonProbes")->evaluate({muo_eta, muo_p, "nominal"});
    } catch (const std::exception &e) {
        std::cerr << "Error in w_global: " << e.what() << std::endl;
    }

    try {
        w_highptid = sf_muoset->at("NUM_HighPtID_DEN_GlobalMuonProbes")->evaluate({muo_eta, muo_pt, "nominal"});
    } catch (const std::exception &e) {
        std::cerr << "Error in w_highptid: " << e.what() << std::endl;
    }

    try {
        w_iso = sf_muoset->at("NUM_probe_LooseRelTkIso_DEN_HighPtProbes")->evaluate({muo_eta, muo_pt, "nominal"});
    } catch (const std::exception &e) {
        std::cerr << "Error in w_iso: " << e.what() << std::endl;
    }

    try {
        w_hlt = sf_muoset->at("NUM_HLT_DEN_HighPtLooseRelIsoProbes")->evaluate({muo_eta, muo_pt, "nominal"});
    } catch (const std::exception &e) {
        std::cerr << "Error in w_hlt: " << e.what() << std::endl;
    }

    weight = w_global * w_highptid * w_iso * w_hlt;
    return weight;

}

//Function to compute PU reweighting correction
inline float compute_PUReweight(const float &PU_NumTrueInteractions){
    float weight = 1.;

    try{
        weight = sf_purew->at("Collisions2023_366403_369802_eraBC_GoldenJson")->evaluate({PU_NumTrueInteractions, "nominal"});
    } catch (const std::exception &e) {
        std::cerr << "Error in w_pu: " << e.what() << std::endl;
    }

    return weight;
}

struct Lepton{
    std::vector<float> pt;
    std::vector<float> eta;
    std::vector<float> phi;
    std::vector<int> pdgId;
    std::vector<float> dR_to_jet;
    std::vector<float> pt_rel_to_jet;
    int n_lep=0;
    std::vector<float> closest_jet_pt; // To store the pt of the closest jet to the lepton, if needed
    std::vector<float> closest_jet_eta; // To store the eta of the closest jet to the lepton, if needed
    std::vector<float> closest_jet_phi; // To store the phi of the closest jet to the lepton, if needed
    std::vector<float> closest_jet_mass; // To store the mass of the closest jet to the lepton, if needed
    std::vector<int> closest_jet_idx; // To store the index of the closest jet to the lepton, if needed
    float closest_jet_threshold = -1; // To store the threshold for the closest jet, if needed

    ClassDef(Lepton, 1) // ROOT necesita esta macro para generar diccionarios
};

struct JetReclus {
    int n_jets = 0;
    float R = 0.0;
    float beta = 0.0;
    float ptmin = 0.0;
    float etamax = 0.0;
    std::vector<bool> findLepton;
    std::vector<float> pt;
    std::vector<float> eta;
    std::vector<float> phi;
    std::vector<float> mass;

//     std::vector<int> topjet_list;
//     std::vector<int> topjet_list2;

    // Constructor explícito para garantizar que se inicialicen a cero
    JetReclus() : n_jets(0), R(0.0), beta(0.0), ptmin(0.0), etamax(0.0), findLepton(false) {}

    ClassDef(JetReclus, 1) // ROOT necesita esta macro para generar diccionarios
};

struct TopJetReclus {
    int n_subjets = 0;
    float ptmin = 0.0;
    float etamax = 0.0;
    float pt;
    float eta;
    float phi;
    float mass;
    std::vector<bool> subjets_in_topjet;
    float pt_W;
    float eta_W;
    float phi_W;
    float mass_W;
    std::vector<bool> subjets_in_W;

    TopJetReclus(): n_subjets(0), ptmin(0.0), etamax(0.0), pt(0.0), eta(0.0), phi(0.0), mass(0.0), subjets_in_topjet(false), pt_W(0.0), eta_W(0.0), phi_W(0.0), mass_W(0.0), subjets_in_W(false) {}

    ClassDef(TopJetReclus, 1) // ROOT necesita esta macro para generar diccionarios

};

struct XConeReclusteredJets {
    JetReclus fatjets;
    JetReclus subjets;
    TopJetReclus topjets;

    std::vector<int> fastjet_CandsList;
    std::vector<int> subjet_CandsList;

    // Constructor por defecto que inicializa las estructuras internas
    XConeReclusteredJets() : fatjets(), subjets(), topjets() {}

    ClassDef(XConeReclusteredJets, 1) // ROOT necesita esta macro para generar diccionarios
};

#endif