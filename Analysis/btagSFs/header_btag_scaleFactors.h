#ifndef __nanoPF_helpers_h__
#define __nanoPF_helpers_h__

#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"
#include <vector>
#include "TVector3.h"
#include "TMath.h"
#include "../deltaR.h"
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
std::unique_ptr<correction::CorrectionSet> cset;
std::unique_ptr<correction::CorrectionSet> sf_btagset;
std::unique_ptr<correction::CorrectionSet> eff_btagset;
std::unique_ptr<correction::CorrectionSet> sf_muoset;
std::unique_ptr<correction::CorrectionSet> sf_purew;
std::unique_ptr<correction::CorrectionSet> jet_veto_map_set;

// Initialize `cset` at the begining of the header
inline void initializeCorrectionSet() {
    fs::path fname_ak4_2023 = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/jet_jerc.json.gz";//"./2023jec/jet_jerc.json.gz";
    // fs::path fname_btag_2023 = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2023_Summer23/btagging.json.gz";
    if (!fs::exists(fname_ak4_2023)) {
        throw std::runtime_error("El archivo de corrección no existe: " + fname_ak4_2023.string());
    }
    cset = correction::CorrectionSet::from_file(fname_ak4_2023.string());
    std::cout << "Archivo de corrección cargado exitosamente: " << fname_ak4_2023.string() << std::endl;
    // if (!fs::exists(fname_btag_2023)) {
    //     throw std::runtime_error("El archivo de corrección de b-tagging no existe: " + fname_btag_2023.string());
    // }
    // cset_btag = correction::CorrectionSet::from_file(fname_btag_2023.string());
    // std::cout << "Archivo de corrección de b-tagging cargado exitosamente: " << fname_btag_2023.string() << std::endl;
    // std::cout << "Corrections set initialized." << std::endl;
}

inline void initializeBTagCorrectionSet() {
    fs::path fname_sf_btag_2023 = "/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/btagging.json.gz";//"./btagSFs/btagging.json.gz"; //"/eos/user/c/cmunozdi/SWAN_projects/BtagScaleFactors/btagging.json.gz"; 
    fs::path fname_eff_btag_2023 = "/eos/project/r/rtu-topanalysis/cmunozdi/AnalysisSamples_JetTightIDNoLepVeto_Full/output_efficiencies_rdf_topSemi_2022postEE/btag_efficiencies_combined.json";//"./btagSFs/btag_efficiencies_combined.json"; //"/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightIDNoLepVeto/output_efficiencies_rdf/btag_efficiencies_combined.json"; 
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
    fs::path fname_sf_MUO_2023 = "/cvmfs/cms-griddata.cern.ch/cat/metadata/MUO/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/muon_HighPt.json.gz";//"./muoSFs/muon_HighPt.json.gz";
    if (!fs::exists(fname_sf_MUO_2023)) {
        throw std::runtime_error("Correction file does not exist: " + fname_sf_MUO_2023.string());
    }
    sf_muoset = correction::CorrectionSet::from_file(fname_sf_MUO_2023.string());
    std::cout << "Correction file loaded susscessfully: " << fname_sf_MUO_2023.string() << std::endl;
}

inline void initializePUReweightingCorrectionSet() {
    fs::path fname_puReweighting_2023 = "/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/puWeights.json.gz";//"./puRew/puWeights.json.gz";
    if (!fs::exists(fname_puReweighting_2023)) {
        throw std::runtime_error("Correction file does not exist: " + fname_puReweighting_2023.string());
    }
    sf_purew = correction::CorrectionSet::from_file(fname_puReweighting_2023.string());
    std::cout << "Correction file loaded susscessfully: " << fname_puReweighting_2023.string() << std::endl;
}

//Initialize jet veto map json file
inline void initializeJetVetoMap(){
    fs::path fname_jetVetoMap_2023 = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/jetvetomaps.json.gz";//"/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-23CSep23-Summer23-NanoAODv12/latest/jetvetomaps.json.gz";
    if (!fs::exists(fname_jetVetoMap_2023)) {
        throw std::runtime_error("Correction file does not exist: " + fname_jetVetoMap_2023.string());
    }
    jet_veto_map_set = correction::CorrectionSet::from_file(fname_jetVetoMap_2023.string());
    std::cout << "Correction file loaded susscessfully: " << fname_jetVetoMap_2023.string() << std::endl;
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
        // std::cerr << "Error in w_global: " << e.what() << std::endl;
    }

    try {
        w_highptid = sf_muoset->at("NUM_HighPtID_DEN_GlobalMuonProbes")->evaluate({muo_eta, muo_pt, "nominal"});
    } catch (const std::exception &e) {
        // std::cerr << "Error in w_highptid: " << e.what() << std::endl;
    }

    try {
        w_iso = sf_muoset->at("NUM_probe_LooseRelTkIso_DEN_HighPtProbes")->evaluate({muo_eta, muo_pt, "nominal"});
    } catch (const std::exception &e) {
        // std::cerr << "Error in w_iso: " << e.what() << std::endl;
    }

    try {
        w_hlt = sf_muoset->at("NUM_HLT_DEN_HighPtLooseRelIsoProbes")->evaluate({muo_eta, muo_pt, "nominal"});
    } catch (const std::exception &e) {
        // std::cerr << "Error in w_hlt: " << e.what() << std::endl;
    }

    weight = w_global * w_highptid * w_iso * w_hlt;
    return weight;

}

//Function to compute PU reweighting correction
inline float compute_PUReweight(const float &PU_NumTrueInteractions){
    float weight = 1.;

    try{
        weight = sf_purew->at("Collisions2022_359022_362760_eraEFG_GoldenJson")->evaluate({PU_NumTrueInteractions, "nominal"}); //at("Collisions2023_366403_369802_eraBC_GoldenJson")->evaluate({PU_NumTrueInteractions, "nominal"});
    } catch (const std::exception &e) {
        // std::cerr << "Error in w_pu: " << e.what() << std::endl;
    }

    return weight;
}

//Funtion to Veto or not events based on jet veto map
inline bool is_event_not_vetoed_by_jetVetoMap(const rvec_f &jet_etas, const rvec_f &jet_phis){

    bool is_vetoed = false;
    //Implement the jet veto map evaluation here for each jet in the event until one vetoes the event or all jets are checked
    for(size_t i=0; i<jet_etas.size(); i++){
        float veto_map_value = 0.;
        try{
            veto_map_value = jet_veto_map_set->at("Summer22EE_23Sep2023_RunEFG_V1")->evaluate({"jetvetomap",jet_etas[i], jet_phis[i]}); // Summer23Prompt23_RunC_V1
        } catch (const std::exception &e) {
            std::cerr << "Error in jet veto map evaluation for jet " << i << ": " << e.what() << std::endl;
            continue;
        }
        if(veto_map_value != 0.) return false; //Event is vetoed
    }

    return true; // Event is not vetoed
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

inline double singleLevelEnergyCorr(const map<string, correction::Variable::Type>& example,
                                    const unique_ptr<correction::CorrectionSet>& cset,
                                    string jec, string lvl, string algo) {
    string key = jec + '_' + lvl + '_' + algo;
    // std::cout << "JSON access to key: " << key << std::endl;

    // correction::Correction::Ref sf = cset->at(key);
    correction::CompoundCorrection::Ref sf = cset->compound().at(key);

    // std::cout << "Inputs:";
    vector<correction::Variable::Type> inputs;
    for (const correction::Variable& input : sf->inputs()) {
        // std::cout << ' ' << input.name();
        inputs.push_back(example.at(input.name()));
        // std::cout << " = " << get<double>(example.at(input.name()));
    }
    // std::cout << std::endl;

    double result = sf->evaluate(inputs);
    // std::cout << "JSON result: " << result << std::endl;

    return result;
}

inline float compute_ptrel_new(const TLorentzVector& muon,
                                const TLorentzVector& jet_corr,
                                float& dR, float jet_threshold) {
    dR = muon.DeltaR(jet_corr);
    if (jet_corr.Pt() < jet_threshold) {
        return -1.0; // Indicating no valid ptrel
    }
    return muon.Vect().Cross(jet_corr.Vect().Unit()).Mag();
}

inline float compute_ptrel_new(const TLorentzVector& muon,
                                TLorentzVector& jet_corr,
                                float& dR,
                                float jet_threshold,
                                float jet_rawFactor,
                                double jet_Area,
                                float rho_PU,
                                float run=-1.) {
    // Raw Jet
    TLorentzVector jet_raw = jet_corr * (1.0 - jet_rawFactor);
    // std::cout << "Jet_raw: Pt = " << jet_raw.Pt() << ", jet_corr: Pt = " << jet_corr.Pt() << ", jet_rawFactor = " << jet_rawFactor << std::endl;

    // Muon
    TLorentzVector muon_raw = muon;

    // Raw jet without the muon
    TLorentzVector jet_raw_no_mu = jet_raw - muon_raw;

    //Variable map with the raw jet
    map<string, correction::Variable::Type> jet_raw_map {
        {"JetPt", jet_raw_no_mu.Pt()},
        {"JetEta", jet_raw_no_mu.Eta()},
        {"JetPhi", jet_raw_no_mu.Phi()},
        {"JetA", jet_Area},
        {"Rho", rho_PU},
        {"run", run}
    };

    // Then the conditions to apply the correction
    string jec = "Summer23Prompt23_V2_DATA",
           lvl = "L1L2L3Res",
           algo = "AK4PFPuppi";

    // Getting the correction factor
    double jet_corrFactor=0.;
    if(jet_Area==0 && rho_PU==0){
        jet_corrFactor = 1.;
        // std::cout << "Jet area and rho_PU are zero, setting jet_corrFactor to 1.0" << std::endl;
    } else{
        jet_corrFactor = singleLevelEnergyCorr(jet_raw_map, cset, jec, lvl, algo);
        // //Imprimir el jet_raw_map
        // std::cout << "Jet raw map for correction: " << std::endl;
        // for (const auto& pair : jet_raw_map) {
        //     std::cout << "\t" << pair.first << ": " << get<double>(pair.second) << std::endl;
        // }    
    }

    // Corrected jet without the muon
    TLorentzVector jet_corr_no_mu = jet_raw_no_mu * jet_corrFactor;
    


    // If the corrected jet without the muon has a pt lower than 15 GeV, return -1 (indicating no valid ptrel)
    if (jet_corr_no_mu.Pt() < jet_threshold) {
        return -1.0;
    }


    // Calculate dR between the muon and the corrected jet without the muon
    dR = muon.DeltaR(jet_corr_no_mu);


    std::cout << "Jet_raw_no_mu: Pt = " << jet_raw_no_mu.Pt() << ", Jet_corr_no_mu: Pt = " << jet_corr_no_mu.Pt() << ", jet_corrFactor = " << jet_corrFactor << ", muonpT: " << muon.Pt() << ", dR: " << dR << ", ptrel: " << muon.Vect().Cross(jet_corr_no_mu.Vect().Unit()).Mag() << std::endl << std::endl;

    //Save the corrected jet without the muon
    jet_corr.SetPtEtaPhiM(jet_corr_no_mu.Pt(),
                          jet_corr_no_mu.Eta(),
                          jet_corr_no_mu.Phi(),
                          jet_corr_no_mu.M());

    // Calculate ptrel
    return muon.Vect().Cross(jet_corr_no_mu.Vect().Unit()).Mag();
}

inline double get_mass(int pdgId) {
    TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(pdgId);
    if (particle) {
        return particle->Mass(); // Devuelve la masa en GeV/c^2
    } else {
        return -1.0; // Indica que la partícula no fue encontrada
    }
}

inline Lepton triggerLepton(const rvec_f &pt_le, const rvec_f &eta_le, const rvec_f &phi_le, const rvec_i &pdgId_le, const rvec_i &jetIdx_le, const rvec_b &tightId_le, /*const rvec_f &ak4_lesubDeltaeta, const rvec_f &ak4_lesubDeltaphi, const rvec_f &ak4_lesubfactor,*/ const rvec_f &PFCands_pt, const rvec_f &PFCands_eta, const rvec_f &PFCands_phi, const rvec_i &PFCands_pdgId, const rvec_f &ak4_pt, const rvec_f &ak4_eta, const rvec_f &ak4_phi, const rvec_f &ak4_mass, const rvec_b &ak4_passJetIdTightLepVeto, const rvec_f &ak4_rawFactor, const rvec_f &ak4_area, const float &rho_PU, const float run, float ak4_pt_threshold = 25.0, const bool isMuon = true, bool AnalysisCuts = false, bool GenLevel = false, const float thr_dR = 0.3, const float thr_ptrel = 30.0) {
    // std::cout << "Is GenLevel: " << GenLevel << std::endl;
    // std::cout << "Trigger Lepton selection: isMuon = " << isMuon << ", AnalysisCuts = " << AnalysisCuts << std::endl;
    // std::cout << "Jet threshold: " << ak4_pt_threshold << std::endl;
    ROOT::VecOps::RVec<bool> tightId_le_local;// = tightId_le; // Crear una copia modificable

    // std::cout << "Redefining ..." << std::endl;

    if (tightId_le.size() == 0) {
        // std::cout << "Warning: tightId_le is empty, assuming all leptons are tight." << std::endl;
        tightId_le_local = rvec_b(pt_le.size(), true); // Asumimos que todos son tight
    }
    else{
        tightId_le_local = tightId_le; // Usar la copia local
    }
    Lepton lepton;
    double pt_lim;
    double eta_max = 99.; //3.
    double ptrel_max = 20.; //30.0
    if (AnalysisCuts) {
        eta_max = 2.4;
        ptrel_max = 30.;
    }

    //Create the correction object
    // fs::path fname_ak4_2023 = "./2023jec/jet_jerc.json.gz";
//TENGO QUE LLEVAR ESTO AL PYTHON FUERA DEL C++#################################
    // fs::path fname_ak4_2023 = "./jet_jerc.json.gz";
    // assert(fs::exists(fname_ak4_2023));
    // unique_ptr<correction::CorrectionSet> cset = correction::CorrectionSet::from_file(fname_ak4_2023.string());

    for (size_t i = 0; i < pt_le.size(); i++) {
        

        // std::cout << "Es un tight lepton: " << i << ", valor --> " << tightId_le_local[i]  << std::endl;
        if((abs(pdgId_le[i]) != 11 && abs(pdgId_le[i]) != 13) || !tightId_le_local[i]) continue;

        
        double min_dR_to_jet = 999.;
        double min_pt_rel_to_jet = -1.;

        if ((isMuon || AnalysisCuts) && abs(pdgId_le[i]) == 13){
            pt_lim = 45.0;
            if (AnalysisCuts) pt_lim = 55.0;   
        } else if ((!isMuon || AnalysisCuts) && abs(pdgId_le[i]) == 11){
            pt_lim = 120.0;
        }
        
        TLorentzVector tem_closest_jet;
        int tmp_closest_jet_idx;

        // Muones y electrones con pt > pt_lim GeV: aplicar criterio de P_T^rel
        if ((pt_le[i] > pt_lim) && (abs(eta_le[i]) < eta_max)) {
            // if(jetIdx_le[i] < 0) continue; // Si el jetIdx es negativo, no pertenece a ningún jet
            
            bool isInsideJet = false;
            TLorentzVector lepton_p4;
            lepton_p4.SetPtEtaPhiM(pt_le[i], eta_le[i], phi_le[i], get_mass(pdgId_le[i]));
            // TVector3 lepton_p3(pt_le[i] * cos(phi_le[i]), pt_le[i] * sin(phi_le[i]), pt_le[i] * sinh(eta_le[i]));

            // std::cout << endl << "Lepton " << i << ": ptx = " << lepton_p3.Px() << ", pty = " << lepton_p3.Py() << ", ptz = " << lepton_p3.Pz() << ", pt = " << pt_le[i] << ", eta = " << eta_le[i] << ", phi = " << phi_le[i] << ", pdgId = " << pdgId_le[i] << endl;
            for (size_t j = 0; j < ak4_eta.size(); j++) {

                //Continue if the jet does not pass the tight ID or is below the threshold
                if (!ak4_passJetIdTightLepVeto[j] ) continue; // Si el jet no pasa el ID tight, lo ignoramos
                //if(jetIdx_le[i] != static_cast<int>(j)) continue; // Only consider the jet that the lepton belongs to
                TLorentzVector jet_corr_p4;
                jet_corr_p4.SetPtEtaPhiM(ak4_pt[j], ak4_eta[j], ak4_phi[j], ak4_mass[j]); // Assuming mass is 0 for jets
                // TVector3 jet_p3(ak4_pt[j] * cos(ak4_phi[j]), ak4_pt[j] * sin(ak4_phi[j]), ak4_pt[j] * sinh(ak4_eta[j]));
                // TVector3 adjusted_jet_p3 = jet_p3;
                // if(deltaR(lepton_p3, jet_p3)<0.4) adjusted_jet_p3 -= lepton_p3; //Only substract the lepton when dR is <0.4

                // std::cout << endl << "\tJet " << j << ": ptx = " << jet_p3.Px() << ", pty = " << jet_p3.Py() << ", ptz = " << jet_p3.Pz();

                float tmp_ptrel;

                float dR;

                // std::cout << endl << "########Jet number: " << j << ", muon_idxJet = " << jetIdx_le[i] << "##########" << endl;
                // std::cout << "Corrected Jet (without muon): Pt = " << (- ak4_lesubfactor[j] + 1) << ", eta = " << ak4_lesubDeltaeta+ak4_eta[j] << ", phi = " << ak4_lesubDeltaphi+ak4_phi[j] << endl; 
                bool muon_belongs_to_jet = (jetIdx_le[i] == static_cast<int>(j));
                //std::cout << "\t muon_belongs_to_jet = " << muon_belongs_to_jet << std::endl;

                // dR = lepton_p4.DeltaR(jet_corr_p4);
                tmp_ptrel = compute_ptrel_new(lepton_p4, jet_corr_p4, dR, ak4_pt_threshold);

                if(muon_belongs_to_jet && dR < 0.4) {

                    tmp_ptrel = compute_ptrel_new(lepton_p4, jet_corr_p4, dR, ak4_pt_threshold, ak4_rawFactor[j], ak4_area[j], rho_PU, run);
                }

                // std::cout << " min_dR_to_jet = " << min_dR_to_jet << ", deltaR = " << deltaR(eta_le[i], phi_le[i], ak4_eta[j], ak4_phi[j])  << ", deltaRpart = " << dR << endl;
                // std::cout << "\tdeltaRadj = " << deltaR(lepton_p3, adjusted_jet_p3)  << endl;

                // If this jet closer than all the previous jets, update the minimum dR and pt_rel
                if ((min_dR_to_jet > dR)&&(tmp_ptrel != -1.0)) {

                    min_pt_rel_to_jet = tmp_ptrel;
                    min_dR_to_jet = dR;
                    tem_closest_jet = jet_corr_p4; // Save the closest jet
                    tmp_closest_jet_idx = j; // Save the index of the closest jet
                }
            }
            // Comment this three lines if you do not want to apply 2D cuts (deltaR and ptrel)
            // if(min_dR_to_jet < 0.4 && min_pt_rel_to_jet < ptrel_max) {
            //     isInsideJet = true;
            // }   
            if ((min_dR_to_jet>thr_dR || min_pt_rel_to_jet>thr_ptrel) && (min_dR_to_jet!= 999. && min_pt_rel_to_jet != -1.)) {
                lepton.pt.push_back(pt_le[i]);
                lepton.eta.push_back(eta_le[i]);
                lepton.phi.push_back(phi_le[i]);
                lepton.pdgId.push_back(pdgId_le[i]);
                lepton.dR_to_jet.push_back(min_dR_to_jet);
                lepton.pt_rel_to_jet.push_back(min_pt_rel_to_jet);
                lepton.closest_jet_pt.push_back(tem_closest_jet.Pt());
                lepton.closest_jet_eta.push_back(tem_closest_jet.Eta());
                lepton.closest_jet_phi.push_back(tem_closest_jet.Phi());
                lepton.closest_jet_mass.push_back(tem_closest_jet.M());
                lepton.closest_jet_idx.push_back(tmp_closest_jet_idx); // Store the index of the closest jet
                lepton.closest_jet_threshold = ak4_pt_threshold; // Store the threshold for the closest jet
                std::cout << "#########dR: " << min_dR_to_jet << ", ptrel: " << min_pt_rel_to_jet << ",########## ptMuon: " << pt_le[i] << ", corrJetPt: " << tem_closest_jet.Pt() << ", beforeCorrJetPt: " << ak4_pt[tmp_closest_jet_idx] << std::endl;
                if(min_dR_to_jet<0.05 && min_pt_rel_to_jet>50){
                    // std::cout << endl << "Lepton " << i << ": ptx = " << lepton_p3.Px() << ", pty = " << lepton_p3.Py() << ", ptz = " << lepton_p3.Pz() << ", eta = " << eta_le[i] << ", phi = " << phi_le[i] << endl;
                    // std::cout << endl << "\tJet " << ": ptx = " << tem_jet_p3.Px() << ", pty = " << tem_jet_p3.Py() << ", ptz = " << tem_jet_p3.Pz() << ", eta = " << tmp_jet_eta << ", phi = " << tmp_jet_phi;
                    // std::cout << ", min_dR_to_jet = " << min_dR_to_jet << ", deltaR = " << std::sqrt(reco::deltaR2(eta_le[i], phi_le[i], tmp_jet_eta, tmp_jet_phi)) << endl;
                    // std::cout << "\t\t tmp_ptrel = " << min_pt_rel_to_jet << endl;
                }
            }
        }

        // Electrones con 45 < pt < 120 GeV: aplicar criterio de aislamiento
        if ((!isMuon || AnalysisCuts) && abs(pdgId_le[i]) == 11 && pt_le[i] < 120) {
            if(AnalysisCuts && pt_le[i] < 60.0) continue;
            else if (!AnalysisCuts && pt_le[i] < 45.0) continue;
            // Calcular aislamiento
            float I_ch = 0.0, I_gamma = 0.0, I_neutral = 0.0;
            for (size_t j = 0; j < PFCands_pt.size(); j++) {
                float dR = deltaR(eta_le[i], phi_le[i], PFCands_eta[j], PFCands_phi[j]);
                if (dR > 0.3) continue; // Fuera del cono de aislamiento

                if (abs(PFCands_pdgId[j]) == 211) { // Hadrón cargado
                    I_ch += PFCands_pt[j];
                } else if (PFCands_pdgId[j] == 22) { // Fotón
                    I_gamma += PFCands_pt[j];
                } else if (abs(PFCands_pdgId[j]) == 130) { // Hadrón neutro
                    I_neutral += PFCands_pt[j];
                }
            }

            // Aislamiento combinado
            float I_combined = I_ch + I_gamma + I_neutral;
            float E_T = pt_le[i];
            float isolation_threshold = 0.0;

            if (abs(eta_le[i]) < 1.479) { // Barrel
                isolation_threshold = 0.029 + 0.51 / E_T;
            } else { // Endcaps
                isolation_threshold = 0.0445 + 0.963 / E_T;
            }

            // Aplicar criterio de aislamiento
            if ((I_combined / E_T) < isolation_threshold) {
                lepton.pt.push_back(pt_le[i]);
                lepton.eta.push_back(eta_le[i]);
                lepton.phi.push_back(phi_le[i]);
                lepton.pdgId.push_back(pdgId_le[i]);
            }
        }
    }

    lepton.n_lep = lepton.pt.size();
    return lepton;
}

//Function that prints the first lepton's variables for debugging
inline void printLeptonDebug(const Lepton& lepton) {
    if (lepton.n_lep > 0) {
        std::cout << "First lepton details:" << std::endl;
        std::cout << "Pt: " << lepton.pt[0] << std::endl;
        std::cout << "Eta: " << lepton.eta[0] << std::endl;
        std::cout << "Phi: " << lepton.phi[0] << std::endl;
        std::cout << "PdgId: " << lepton.pdgId[0] << std::endl;
        std::cout << "dR to jet: " << lepton.dR_to_jet[0] << std::endl;
        std::cout << "Pt rel to jet: " << lepton.pt_rel_to_jet[0] << std::endl;
    } else {
        std::cout << "No leptons found." << std::endl;
    }
}


inline Lepton CombineLeptons(const Lepton& muons, const Lepton& electrons) {
    Lepton combined;

    // Combinar los vectores de propiedades
    combined.pt.insert(combined.pt.end(), muons.pt.begin(), muons.pt.end());
    combined.pt.insert(combined.pt.end(), electrons.pt.begin(), electrons.pt.end());

    combined.eta.insert(combined.eta.end(), muons.eta.begin(), muons.eta.end());
    combined.eta.insert(combined.eta.end(), electrons.eta.begin(), electrons.eta.end());

    combined.phi.insert(combined.phi.end(), muons.phi.begin(), muons.phi.end());
    combined.phi.insert(combined.phi.end(), electrons.phi.begin(), electrons.phi.end());

    combined.pdgId.insert(combined.pdgId.end(), muons.pdgId.begin(), muons.pdgId.end());
    combined.pdgId.insert(combined.pdgId.end(), electrons.pdgId.begin(), electrons.pdgId.end());

    combined.dR_to_jet.insert(combined.dR_to_jet.end(), muons.dR_to_jet.begin(), muons.dR_to_jet.end());
    combined.dR_to_jet.insert(combined.dR_to_jet.end(), electrons.dR_to_jet.begin(), electrons.dR_to_jet.end());

    combined.pt_rel_to_jet.insert(combined.pt_rel_to_jet.end(), muons.pt_rel_to_jet.begin(), muons.pt_rel_to_jet.end());
    combined.pt_rel_to_jet.insert(combined.pt_rel_to_jet.end(), electrons.pt_rel_to_jet.begin(), electrons.pt_rel_to_jet.end());

    combined.closest_jet_pt.insert(combined.closest_jet_pt.end(), muons.closest_jet_pt.begin(), muons.closest_jet_pt.end());
    combined.closest_jet_pt.insert(combined.closest_jet_pt.end(), electrons.closest_jet_pt.begin(), electrons.closest_jet_pt.end());

    combined.closest_jet_eta.insert(combined.closest_jet_eta.end(), muons.closest_jet_eta.begin(), muons.closest_jet_eta.end());
    combined.closest_jet_eta.insert(combined.closest_jet_eta.end(), electrons.closest_jet_eta.begin(), electrons.closest_jet_eta.end());

    combined.closest_jet_phi.insert(combined.closest_jet_phi.end(), muons.closest_jet_phi.begin(), muons.closest_jet_phi.end());
    combined.closest_jet_phi.insert(combined.closest_jet_phi.end(), electrons.closest_jet_phi.begin(), electrons.closest_jet_phi.end());

    combined.closest_jet_mass.insert(combined.closest_jet_mass.end(), muons.closest_jet_mass.begin(), muons.closest_jet_mass.end());
    combined.closest_jet_mass.insert(combined.closest_jet_mass.end(), electrons.closest_jet_mass.begin(), electrons.closest_jet_mass.end());

    combined.closest_jet_idx.insert(combined.closest_jet_idx.end(), muons.closest_jet_idx.begin(), muons.closest_jet_idx.end());
    combined.closest_jet_idx.insert(combined.closest_jet_idx.end(), electrons.closest_jet_idx.begin(), electrons.closest_jet_idx.end());


    // Actualizar el número total de leptones
    combined.n_lep = muons.n_lep + electrons.n_lep;

    // Actualizar el umbral del jet más cercano
    combined.closest_jet_threshold = std::max(muons.closest_jet_threshold, electrons.closest_jet_threshold);

    //Print the first lepton's variables for debugging
    // printLeptonDebug(combined);

    return combined;
}





#endif