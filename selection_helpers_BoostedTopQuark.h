#ifndef __nanoPF_helpers_h__
#define __nanoPF_helpers_h__

#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"
#include <vector>
#include "TVector3.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "Rivet/Tools/ParticleIdUtils.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include <fastjet/PseudoJet.hh>
// #include "/afs/cern.ch/user/c/cmunozdi/Analysis/GeneratePFNanos/CMSSW_13_0_13/src/fastjet-install/include/fastjet/contrib/XConePlugin.hh" //RUN LOCAL
// #include "/afs/cern.ch/user/c/cmunozdi/Analysis/GeneratePFNanos/CMSSW_13_0_13/src/fastjet-install/include/fastjet/contrib/Nsubjettiness.hh" //RUN LOCAL
#include "fastjet/contrib/XConePlugin.hh" //RUN CRAB
#include "fastjet/contrib/Nsubjettiness.hh" //RUN CRAB
#include "deltaR.h"


#include <TObject.h>


#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include <cassert>
#include <filesystem>
#include <iostream>
#include <string>
#include "correction.h"

using namespace ROOT::VecOps;
using namespace fastjet;
using namespace contrib;
using namespace std;
namespace fs = std::filesystem;

using rvec_f = const RVec<float>;
using rvec_i = const RVec<int>;
using rvec_b = const RVec<bool>;

// Variable global para el archivo de corrección de jets
extern std::unique_ptr<correction::CorrectionSet> cset;
// extern std::unique_ptr<correction::CorrectionSet> cset_btag;
extern std::unique_ptr<correction::CorrectionSet> sf_btagset;
extern std::unique_ptr<correction::CorrectionSet> eff_btagset;
extern std::unique_ptr<correction::CorrectionSet> sf_muoset;
extern std::unique_ptr<correction::CorrectionSet> sf_purew;

// NEW: global tags for which correction to pick inside the JSON. Change values as needed in selection_helpers_BoostedTopQuark.cpp
extern std::string g_jec_tag;  //= "Summer23Prompt23_V2_MC";
extern std::string g_lvl_tag;  //= "L1L2L3Res";
extern std::string g_algo_tag; // = "AK4PFPuppi";

extern bool g_isMC;

// Inicializar `cset` al comienzo del header
inline void initializeCorrectionSet(bool isMC = true, const std::string &jec_override = "") {
    // choose which key to use inside the single JSON (can be overridden)
    std::cout << "Initializing JEC CorrectionSet with isMC = " << isMC << " and jec_override = '" << jec_override << "'" << std::endl;

    g_isMC = isMC;
    if (!jec_override.empty()) {
        g_jec_tag = jec_override;
    } else {
        g_jec_tag = g_isMC ? "Summer23BPixPrompt23_V3_MC" : "Summer23BPixPrompt23_V3_DATA"; //"Summer22EE_22Sep2023_V3_MC" : "Summer22EE_22Sep2023_RunE_V3_DATA"; //"Summer23Prompt23_V2_MC" : "Summer23Prompt23_V2_DATA";
    }
    // lvl/algo kept as defaults but can be changed from Python if desired:
    std::cout << "Using JEC tag: " << g_jec_tag << std::endl;
    g_lvl_tag  = "L1L2L3Res";
    g_algo_tag = "AK4PFPuppi";

    fs::path fname_ak4_2023 = "./jet_jerc.json.gz"; //"/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-23CSep23-Summer23-NanoAODv12/latest/jet_jerc.json.gz"; For local, use the whole path. For CRAB, use relative path and upload the file with the job.
    if (!fs::exists(fname_ak4_2023)) {
        throw std::runtime_error("El archivo de corrección no existe: " + fname_ak4_2023.string());
    }
    cset = correction::CorrectionSet::from_file(fname_ak4_2023.string());
    std::cout << "Archivo de corrección cargado exitosamente: " << fname_ak4_2023.string() << std::endl;
}

extern bool debug;
// Declarar `cset` como una variable global externa
// extern std::unique_ptr<correction::CorrectionSet> cset;

/**
 * @brief Returns the particle's mass from its pdgId
 *
 * @param pdgId
 * @return double particle mass in GeV/c^2. Returns -1.0 if it does not find the particle
 */
inline double get_mass(int pdgId) {
    TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(pdgId);
    if (particle) {
        return particle->Mass(); // Devuelve la masa en GeV/c^2
    } else {
        return -1.0; // Indica que la partícula no fue encontrada
    }
}

struct Lepton{
    std::vector<float> pt;
    std::vector<float> eta;
    std::vector<float> phi;
    std::vector<int> pdgId;
    std::vector<int> leptonIdx; // Index of the lepton in the original collection
    std::vector<float> dR_to_jet;
    std::vector<float> pt_rel_to_jet;
    int n_lep=0;
    int n_lep_after2Dcuts=0; // New variable to count leptons after 2D cuts
    std::vector<float> closest_jet_pt; // To store the pt of the closest jet to the lepton, if needed
    std::vector<float> closest_jet_eta; // To store the eta of the closest jet to the lepton, if needed
    std::vector<float> closest_jet_phi; // To store the phi of the closest jet to the lepton, if needed
    std::vector<float> closest_jet_mass; // To store the mass of the closest jet to the lepton, if needed
    std::vector<int> closest_jet_idx; // To store the index of the closest jet to the lepton, if needed
    float closest_jet_threshold = -1; // To store the threshold for the closest jet, if needed

    ClassDef(Lepton, 1) // ROOT necesita esta macro para generar diccionarios
};

inline double singleLevelEnergyCorr(const map<string, correction::Variable::Type>& example) {
    string key = g_jec_tag + '_' + g_lvl_tag + '_' + g_algo_tag;
    // string key = jec + '_' + lvl + '_' + algo;
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
        {"Rho", rho_PU}
    };

    if(!g_isMC){
        // std::cout << "Data event: adding run number to JEC inputs." << std::endl;
        if (run < 0) {
            std::cerr << "compute_ptrel_new: run not provided for DATA but required by correction. Using run=0 as fallback." << std::endl;
            jet_raw_map["run"] = 0.;
        } else {
            // std::cout << "compute_ptrel_new: using run=" << run << " for DATA correction." << std::endl;
            jet_raw_map["run"] = run;
        }
    }

    // Then the conditions to apply the correction
    // string jec = "Summer23Prompt23_V2_DATA",
    //        lvl = "L1L2L3Res",
    //        algo = "AK4PFPuppi";

    // Getting the correction factor
    double jet_corrFactor=0.;
    if(jet_Area==0 && rho_PU==0){
        jet_corrFactor = 1.;
        // std::cout << "Jet area and rho_PU are zero, setting jet_corrFactor to 1.0" << std::endl;
    } else{
        jet_corrFactor = singleLevelEnergyCorr(jet_raw_map);
        // //Imprimir el jet_raw_map
        // std::cout << "Jet raw map for correction: " << std::endl;
        // for (const auto& pair : jet_raw_map) {
        //     std::cout << "\t" << pair.first << ": " << get<double>(pair.second) << std::endl;
        // }    
    }

    // Corrected jet without the muon
    TLorentzVector jet_corr_no_mu = jet_raw_no_mu * jet_corrFactor;
    // std::cout << "Jet_raw_no_mu: Pt = " << jet_raw_no_mu.Pt() << ", Jet_corr_no_mu: Pt = " << jet_corr_no_mu.Pt() << ", jet_corrFactor = " << jet_corrFactor << ", jet_area = " << jet_Area << ", rho_PU = " << rho_PU << std::endl << std::endl;


    // If the corrected jet without the muon has a pt lower than 15 GeV, return -1 (indicating no valid ptrel)
    if (jet_corr_no_mu.Pt() < jet_threshold) {
        return -1.0;
    }
    // Calculate dR between the muon and the corrected jet without the muon
    dR = muon.DeltaR(jet_corr_no_mu);

    //Save the corrected jet without the muon
    jet_corr.SetPtEtaPhiM(jet_corr_no_mu.Pt(),
                          jet_corr_no_mu.Eta(),
                          jet_corr_no_mu.Phi(),
                          jet_corr_no_mu.M());

    // Calculate ptrel
    return muon.Vect().Cross(jet_corr_no_mu.Vect().Unit()).Mag();
}

inline ROOT::VecOps::RVec<Short_t> calculateGenCandsJetIdx(const ROOT::VecOps::RVec<int>& GenJetCands_genCandsIdx,
                                                           const ROOT::VecOps::RVec<int>& GenJetCands_jetIdx,
                                                           int nGenCands, int nGenJet) {
    // Inicializa el vector de índices con -1 (no pertenece a ningún jet)
    ROOT::VecOps::RVec<Short_t> GenCands_jetIdx(nGenCands, -1);

    // Itera sobre las GenJetCands para asignar los índices de los jets
    for (size_t i = 0; i < GenJetCands_genCandsIdx.size(); ++i) {
        int genCandIdx = GenJetCands_genCandsIdx[i];
        if (genCandIdx >= 0 && genCandIdx < nGenCands) {
            int jetIdx = GenJetCands_jetIdx[i];
            if (jetIdx >= 0 && jetIdx < nGenJet) {
                GenCands_jetIdx[genCandIdx] = static_cast<Short_t>(jetIdx);
            }
        }
    }

    return GenCands_jetIdx;
}

inline Lepton pickLepton(const rvec_f &pt_mu, const rvec_f &eta_mu, const rvec_f &phi_mu, const rvec_i &pdgId_mu, const rvec_f &pt_e, const rvec_f &eta_e, const rvec_f &phi_e, const rvec_i &pdgId_e){
    //Count how many leptons with pt>55 and eta<2.4 are in the event. If more than one, or none, return empty Lepton struct. Otherwise, return the lepton.
    Lepton lepton;

    for (size_t i = 0; i < pt_mu.size(); i++) {
        if (pt_mu[i] > 55.0 && abs(eta_mu[i]) < 2.4 && abs(pdgId_mu[i]) == 13) {
            lepton.pt.push_back(pt_mu[i]);
            lepton.eta.push_back(eta_mu[i]);
            lepton.phi.push_back(phi_mu[i]);
            lepton.pdgId.push_back(pdgId_mu[i]);
            lepton.leptonIdx.push_back(i);
        }
    }

    for (size_t i = 0; i < pt_e.size(); i++) {
        if (pt_e[i] > 55.0 && abs(eta_e[i]) < 2.4 && abs(pdgId_e[i]) == 11) {
            lepton.pt.push_back(pt_e[i]);
            lepton.eta.push_back(eta_e[i]);
            lepton.phi.push_back(phi_e[i]);
            lepton.pdgId.push_back(pdgId_e[i]);
            lepton.leptonIdx.push_back(i);
        }
    }

    lepton.n_lep = lepton.pt.size();
    return lepton;
}

//Function to the lepton struct the information related to the closest jet
inline Lepton addLeptonJetInfo(const Lepton &lepton_in, const rvec_i &jetIdx_mu, const rvec_i &jetIdx_e, const rvec_f &pt_pfcands, const rvec_f &eta_pfcands, const rvec_f &phi_pfcands, const rvec_i &pdgId_pfcands, const rvec_f &ak4_pt, const rvec_f &ak4_eta, const rvec_f &ak4_phi, const rvec_f &ak4_mass, const rvec_i &ak4_jetId, const rvec_f &ak4_rawFactor, const rvec_f &ak4_area, const float &rho_PU, const float run, float ak4_pt_threshold = 15.0) {
    //First thing is to check that lepton_in only contains 1 single lepton. If not, return empty struct.
    Lepton lepton_out;  // Create empty output struct
    if (lepton_in.n_lep != 1) {
        return lepton_out;  // Return empty if input doesn't have exactly 1 lepton
    }

    //Secondly, for electrons with pt<120 GeV, we do not add the jet info because the isolation is computed differently.
    if (abs(lepton_in.pdgId[0]) == 11 && lepton_in.pt[0] <= 120.0) {

        float I_ch = 0.0, I_gamma = 0.0, I_neutral = 0.0;
        for (size_t j = 0; j < pt_pfcands.size(); j++) {
            float dR = deltaR(lepton_in.eta[0], lepton_in.phi[0], eta_pfcands[j], phi_pfcands[j]);
            if (dR > 0.3) continue; // Outside the isolation cone

            if (abs(pdgId_pfcands[j]) == 211) { // Charged hadron
                I_ch += pt_pfcands[j];
            } else if (pdgId_pfcands[j] == 22) { // Photon
                I_gamma += pt_pfcands[j];
            } else if (abs(pdgId_pfcands[j]) == 130) { // Neutral hadron
                I_neutral += pt_pfcands[j];
            }
        }

        // Combined isolation
        float I_combined = I_ch + I_gamma + I_neutral;
        float E_T = lepton_in.pt[0];
        float isolation_threshold = 0.0;

        if (abs(lepton_in.eta[0]) < 1.479) { // Barrel
            isolation_threshold = 0.029 + 0.51 / E_T;
        } else { // Endcaps
            isolation_threshold = 0.0445 + 0.963 / E_T;
        }

        // Apply isolation criterion
        if ((I_combined / E_T) < isolation_threshold) {
            lepton_out.pt.push_back(lepton_in.pt[0]);
            lepton_out.eta.push_back(lepton_in.eta[0]);
            lepton_out.phi.push_back(lepton_in.phi[0]);
            lepton_out.pdgId.push_back(lepton_in.pdgId[0]);
            lepton_out.leptonIdx.push_back(lepton_in.leptonIdx[0]);
            lepton_out.n_lep_after2Dcuts = 1; // Mark that the lepton passed the 2D cuts
            lepton_out.n_lep = 1;
        }
    }

    // For muons and high-pt electrons, we add the jet info
    else if ( (abs(lepton_in.pdgId[0]) == 13) || (abs(lepton_in.pdgId[0]) == 11 && lepton_in.pt[0] > 120.0) ) {
        //Get the index of the jet according to the BTVnano
        int leptonIdx = lepton_in.leptonIdx[0];
        int jetIdx = -1;
        if (abs(lepton_in.pdgId[0]) == 13) { // Muon
            jetIdx = jetIdx_mu[leptonIdx];
        } else if (abs(lepton_in.pdgId[0]) == 11) { // Electron
            jetIdx = jetIdx_e[leptonIdx];
        }

        // Lorentz vector of the lepton
        TLorentzVector lepton_p4;
        lepton_p4.SetPtEtaPhiM(lepton_in.pt[0], lepton_in.eta[0], lepton_in.phi[0], get_mass(lepton_in.pdgId[0]));

        // Variables to store the closest jet info
        double min_dR_to_jet = 999.;
        double min_pt_rel_to_jet = -1.;
        TLorentzVector tem_closest_jet;
        int tmp_closest_jet_idx;

        //Run over all AK4 jets to find the closest one
        for (size_t j = 0; j < ak4_pt.size(); j++) {

            //First, if jet does not pass the tight lepton Vecto, skip it
            if (ak4_jetId[j] < 6) continue;

            // Second, build lorentz vector of the jet
            TLorentzVector jet_p4;
            jet_p4.SetPtEtaPhiM(ak4_pt[j], ak4_eta[j], ak4_phi[j], ak4_mass[j]);

            float tmp_ptrel = -1.;
            float dR = -1.;

            // The lepton belongs to the jet if its index matches the jetIdx from BTVnano
            bool leptonBelongsToJet = (static_cast<int>(j) == jetIdx);

            //Compute a first stimated ptrel and dR
            tmp_ptrel = compute_ptrel_new(lepton_p4, jet_p4, dR, ak4_pt_threshold);

            if (leptonBelongsToJet && dR < 0.4) {
                // If the lepton belongs to the jet, recompute ptrel with the new method, substracting the lepton from the jet
                tmp_ptrel = compute_ptrel_new(lepton_p4, jet_p4, dR, ak4_pt_threshold, ak4_rawFactor[j], ak4_area[j], rho_PU, run);
            }

            // Store the closest jet info
            if ((min_dR_to_jet > dR) && (tmp_ptrel != -1.0)){
                min_dR_to_jet = dR;
                min_pt_rel_to_jet = tmp_ptrel;
                tem_closest_jet = jet_p4; // Save the closest jet p4
                tmp_closest_jet_idx = j; // Save the index of the closest jet
            }
        }

        //Save the lepton info
        if((min_dR_to_jet != 999.)&&(min_pt_rel_to_jet != -1.0)){
            lepton_out.pt.push_back(lepton_in.pt[0]);
            lepton_out.eta.push_back(lepton_in.eta[0]);
            lepton_out.phi.push_back(lepton_in.phi[0]);
            lepton_out.pdgId.push_back(lepton_in.pdgId[0]);
            lepton_out.leptonIdx.push_back(lepton_in.leptonIdx[0]);
            lepton_out.dR_to_jet.push_back(min_dR_to_jet);
            lepton_out.pt_rel_to_jet.push_back(min_pt_rel_to_jet);
            lepton_out.closest_jet_pt.push_back(tem_closest_jet.Pt());;
            lepton_out.closest_jet_eta.push_back(tem_closest_jet.Eta());;
            lepton_out.closest_jet_phi.push_back(tem_closest_jet.Phi());;
            lepton_out.closest_jet_mass.push_back(tem_closest_jet.M());;
            lepton_out.closest_jet_idx.push_back(tmp_closest_jet_idx);
            lepton_out.closest_jet_threshold = ak4_pt_threshold;
            lepton_out.n_lep = 1;
            lepton_out.n_lep_after2Dcuts = 0; // Mark that the lepton no passed the 2D cuts
            if(min_dR_to_jet > 0.3 || min_pt_rel_to_jet > 30.0){
                lepton_out.n_lep_after2Dcuts = 1; // Mark that the lepton passed the 2D cuts
            }
        }
    }

    return lepton_out;
}

//Same function pickLepton but for gen level particles, i.e. there are no GenMuons or GenElectrons, only GenCands
inline Lepton pickGenLepton(const rvec_f &pt_gencands, const rvec_f &eta_gencands, const rvec_f &phi_gencands, const rvec_i &pdgId_gencands){
    //Count how many leptons with pt>55 and eta<2.4 are in the event. If more than one, or none, return empty Lepton struct. Otherwise, return the lepton.
    Lepton lepton;

    for (size_t i = 0; i < pt_gencands.size(); i++) {
        if (pt_gencands[i] > 55.0 && abs(eta_gencands[i]) < 2.4 && (abs(pdgId_gencands[i]) == 11 || abs(pdgId_gencands[i]) == 13)) {
            lepton.pt.push_back(pt_gencands[i]);
            lepton.eta.push_back(eta_gencands[i]);
            lepton.phi.push_back(phi_gencands[i]);
            lepton.pdgId.push_back(pdgId_gencands[i]);
            lepton.leptonIdx.push_back(i);
        }
    }

    lepton.n_lep = lepton.pt.size();
    return lepton;
}

// Same function as addLeptonJetInfo but for gen level leptons. As we do not have to check trigger selections on gen level, this function should unify both: pickLepton and addLeptonJetInfo, i.e. it should pick events with exactly one gen lepton, and then pass the 2D isolation if muon or high-pt electron, while for low-pt electrons it should compute the isolation as in addLeptonJetInfo.
inline Lepton GenLepton(const rvec_f &pt_gencands, const rvec_f &eta_gencands, const rvec_f &phi_gencands, const rvec_i &pdgId_gencands, const rvec_i &jetIdx_gencands, const rvec_f &genak4_pt, const rvec_f &genak4_eta, const rvec_f &genak4_phi, const rvec_f &genak4_mass) {
    //First thing is to build the lepton struct with events with exactly one lepton
    Lepton pick_lepton = pickGenLepton(pt_gencands, eta_gencands, phi_gencands, pdgId_gencands);

    //Once we have the single lepton struct, we can reuse the same function addLeptonJetInfo, passing trivial values for Jet_rawFactor, Jet_area, rho and run, as they are not needed at gen level. For gen level we do not have separate muon/electron jet indices, so we pass the same for both: in the function, we use the index we saved in the lepton struct, so it does not matter.
    Lepton lepton_out = addLeptonJetInfo(pick_lepton, jetIdx_gencands, jetIdx_gencands, pt_gencands, eta_gencands, phi_gencands, pdgId_gencands, genak4_pt, genak4_eta, genak4_phi, genak4_mass, rvec_i(genak4_pt.size(), 7), rvec_f(genak4_pt.size(), 0.0), rvec_f(genak4_pt.size(), 0.0), 0.0, 0.0);

    return lepton_out;
}

inline ROOT::RVec<bool> pass_JetIdTightLepVeto(const ROOT::RVec<float>& Jet_eta,
                                   const ROOT::RVec<float>& Jet_neHEF, const ROOT::RVec<float>& Jet_neEmEF,
                                   const ROOT::RVec<int>& Jet_chMultiplicity, const ROOT::RVec<int>& Jet_neMultiplicity,
                                   const ROOT::RVec<float>& Jet_chHEF, const ROOT::RVec<float>& Jet_muEF,
                                   const ROOT::RVec<float>& Jet_chEmEF, bool LepVeto = false){
    ROOT::RVec<bool> Jet_passJetIdTight = rvec_b(Jet_eta.size(), false);
    ROOT::RVec<bool> Jet_passJetIdTightLepVeto = rvec_b(Jet_eta.size(), false);
    for(size_t jet_idx = 0; jet_idx < Jet_eta.size(); jet_idx++) {
        if(abs(Jet_eta[jet_idx]) <= 2.6){
            Jet_passJetIdTight[jet_idx] = (Jet_neHEF[jet_idx] < 0.99) && (Jet_neEmEF[jet_idx] < 0.9) && (Jet_chMultiplicity[jet_idx]+Jet_neMultiplicity[jet_idx] > 1) && (Jet_chHEF[jet_idx] > 0.01) && (Jet_chMultiplicity[jet_idx] > 0);
        }else if((abs(Jet_eta[jet_idx]) > 2.6) && (abs(Jet_eta[jet_idx]) <= 2.7)){
            Jet_passJetIdTight[jet_idx] = (Jet_neHEF[jet_idx] < 0.90) && (Jet_neEmEF[jet_idx] < 0.99);
        }else if((abs(Jet_eta[jet_idx]) > 2.7) && (abs(Jet_eta[jet_idx]) <= 3.0)){
            Jet_passJetIdTight[jet_idx] = (Jet_neHEF[jet_idx] < 0.90);
        }else if((abs(Jet_eta[jet_idx]) > 3.0) ){
            Jet_passJetIdTight[jet_idx] = (Jet_neMultiplicity[jet_idx] >=2 ) && (Jet_neEmEF[jet_idx] < 0.4);
        }

        // Now apply the lepton veto
        if(LepVeto){
            if(abs(Jet_eta[jet_idx]) <= 2.7){
                Jet_passJetIdTightLepVeto[jet_idx] = Jet_passJetIdTight[jet_idx] && (Jet_muEF[jet_idx] < 0.8) && (Jet_chEmEF[jet_idx] < 0.8);
            }else Jet_passJetIdTightLepVeto[jet_idx] = Jet_passJetIdTight[jet_idx];
        }
    }
    if(LepVeto) return Jet_passJetIdTightLepVeto;
    else return Jet_passJetIdTight;
    // return Jet_passJetIdTightLepVeto;
}

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
    std::vector<float> area;

//     std::vector<int> topjet_list;
//     std::vector<int> topjet_list2;

    // Constructor explícito para garantizar que se inicialicen a cero
    JetReclus() : n_jets(0), R(0.0), beta(0.0), ptmin(0.0), etamax(0.0), findLepton(false) {}

    ClassDef(JetReclus, 2) // ROOT necesita esta macro para generar diccionarios (version 2: añadido area)
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



/*BUILDING XCONE RECLUSTERING FUNCTIONS*/
inline void initPlugin(std::unique_ptr<NjettinessPlugin> & ptr, int N, float R0, float beta, bool usePseudoXCone){
    if(usePseudoXCone){
        ptr.reset(new PseudoXConePlugin(N, R0, beta));
    } else {
        ptr.reset(new XConePlugin(N, R0, beta));
    }
}

inline XConeReclusteredJets buildXConeJets(const Lepton &lep, const rvec_f &pt, const rvec_f &eta, const rvec_f &phi, const rvec_f &mass, const rvec_i &pdgId, int NJets = 2, float RJets = 1.2, float BetaJets = 2.0, float ptminJets = 200., float etamaxJets = 3., int NSubJets = 3, float RSubJets = 0.4, float BetaSubJets = 2., float ptminSubJets = 22.5, float etamaxSubJets = 3., bool doLeptonSpecific = true, bool usePseudoXCone = false){

    //Return empty XConeReclusteredJets if lepton struct has more than one lepton
    if (doLeptonSpecific && lep.n_lep != 1){
        if(debug) std::cout << "More than one lepton found (or none) in buildXConeJets" << std::endl;
        XConeReclusteredJets XConeJetsCollection;
        return XConeJetsCollection;
    }
    // else std::cout << "Number of leptons found in buildXConeJets: " << lep.n_lep << std::endl;
    // const reco::Candidate *lepton(nullptr);
    std::vector<fastjet::PseudoJet>::iterator lepton_iter;// = _psj.end();
    float lepton_min_pt = 55;
    float DRLeptonJet_ = 999.;

    //Convert particles to PseudoJets
    std::vector<fastjet::PseudoJet> _psj;
    lepton_iter = _psj.end();
    int i_gl=-1;
    int lepton_index = -1;

    fastjet::PseudoJet lepton;
    for (size_t i = 0; i < pt.size(); i++){
        i_gl++;
        ROOT::Math::PtEtaPhiMVector p4(pt[i], eta[i], phi[i], mass[i]);
        if (std::isnan(p4.px())||
            std::isnan(p4.py())||
            std::isnan(p4.pz())||
            std::isinf(p4.px())||
            std::isinf(p4.py())||
            std::isinf(p4.pz())||
            (p4.pt() == 0)) continue;

        fastjet::PseudoJet tmp_particle = fastjet::PseudoJet(p4.px(), p4.py(), p4.pz(), p4.energy());
        tmp_particle.set_user_index(i);
        _psj.push_back(tmp_particle);
        // if(debug) std::cout << "Direccion de memoria de _psj[0]: " << &_psj[0] << std::endl;

        if (doLeptonSpecific){
            uint pdgid = abs(pdgId[i]);
            auto candPt = p4.pt();
            // if((pdgid==11 || pdgid==13) && candPt > 55.0){
            //     std::cout << "\tAdded PFCand with pt: " << p4.pt() << ", eta: " << p4.eta() << ", phi: " << p4.phi() << ", mass: " << p4.mass() << ", pdgId: " << pdgId[i] << std::endl;
            // }
            if ((pdgid == 11 || pdgid == 13) && (candPt > lepton_min_pt)){
                // lepton_iter = std::prev(_psj.end());
                lepton_index = i_gl;
                lepton_min_pt = candPt;
                lepton = _psj.back();
            }
            // if(puppiw[i] >= 0 && i<50){
            //     // std::cout << "\tAdded PFCand with pt: " << p4.pt() << ", eta: " << p4.eta() << ", phi: " << p4.phi() << ", mass: " << p4.mass() << ", pdgId: " << pdgId[i] << ", puppiW: " << puppiw[i] << std::endl;
            //     // std::cout << p4.pt() << ";" << p4.eta() << ";" << p4.phi() << ";" << p4.mass() << ";" << pdgId[i] << ";" << puppiw[i] << ";" << puppiwnolep[i] << ";" << std::endl;
            // }
        }//doLeptonSpecific
    }//for loop over pfcands of the event

    bool foundLepton = false;

    if (doLeptonSpecific && lepton != 0.){
        if(abs(lepton_min_pt - lep.pt[0])/lep.pt[0] > 0.00000001){
            // std::cout << "PT del lepton PFCand es: " << lepton_min_pt << ", eta: " << _psj[lepton_index].eta() << ", phi: " << _psj[lepton_index].phi() <<  ", pdgId: " << pdgId[lepton_index] << ", puppiW: " << puppiw[lepton_index] << std::endl;
            // std::cout << "El pT del lepton Muon/Electron es: " << lep.pt[0] << ", eta: " << lep.eta[0] << ", phi: " << lep.phi[0] << ", pdgId: " << lep.pdgId[0] << std::endl;
        }
        if(debug) std::cout << "PT del lepton es: " << lepton_min_pt << std::endl;
        // lepton = *lepton_iter;
        if(debug) std::cout << "PT del lepton es1: " << lepton.pt() << std::endl;
        foundLepton = true;
        if(debug) std::cout << "PT del lepton es: " << lepton_min_pt << std::endl;
    }

    // if(!foundLepton){
    //     //std::cout << "foundLepton: " << foundLepton << std::endl;
    //     //std::cout << "doLeptonSpecific: " << doLeptonSpecific << std::endl;
    // }

    if (doLeptonSpecific && !foundLepton){
        if(debug) std::cout << "foundLepton: " << foundLepton << std::endl;
        if(debug) std::cout << "No lepton found in buildXConeJets" << std::endl;
        XConeReclusteredJets XConeJetsCollection;
        return XConeJetsCollection;
    }

    //Make sure to have enough particles with non-zero momentum in the event
    unsigned int minParticleCount = NJets * NSubJets;
    if (doLeptonSpecific && foundLepton) minParticleCount -= 1; //since one jet will use NSubJets-1 (cause of the neutrino)



    //Declare jet collections
    XConeReclusteredJets XConeJetsCollection;
    if(_psj.size() < minParticleCount){
        return XConeJetsCollection;
    }

    //Run clustering of fatjets
    std::unique_ptr<NjettinessPlugin> plugin_xcone;
    initPlugin(plugin_xcone, NJets, RJets, BetaJets, usePseudoXCone);
    JetDefinition jet_def_xcone(plugin_xcone.get());
    // Add GhostedAreaSpec for jet area calculation
    fastjet::GhostedAreaSpec ghost_spec_xcone(5.0);  // eta_max = 5.0
    fastjet::ClusterSequenceArea clust_seq_xcone(_psj, jet_def_xcone, ghost_spec_xcone);
    Selector selector_fatjets = SelectorPtMin(10); //pt larger than 10GeV (this is the loosest criteria, for the jet where the lepton is clustered; the main fatjet would have to be larger than 400GeV)
    std::vector<fastjet::PseudoJet> fatjets = sorted_by_pt(selector_fatjets(clust_seq_xcone.inclusive_jets()));

    //Print some info about the fatjets like number, pt, eta, phi, mass, area, etc.
    // std::cout << "Number of fatjets clustered: " << fatjets.size() << std::endl;
    // for (unsigned int i = 0; i < fatjets.size(); i++){
    //     const auto &fj = fatjets[i];
    //     std::cout << "\tFatjet " << i << ": pt = " << fj.pt() << ", eta = " << fj.eta() << ", phi = " << fj.phi_std() << ", mass = " << fj.m() << ", area = " << fj.area() << std::endl;
    // }



    //First fet and write list of particles in fatjets:
    //If particle i is clustered in jet j, the i-th entry of the list == j
    vector<int> list_fat;
    list_fat.clear();
    list_fat = clust_seq_xcone.particle_jet_indices(fatjets);
    //This is the mapping of subjet to hard jet
    std::vector<std::vector<int>> indices;
    indices.resize(fatjets.size());

    for(const auto &fjet: fatjets){
        XConeJetsCollection.fatjets.pt.push_back(fjet.pt());
        XConeJetsCollection.fatjets.eta.push_back(fjet.eta());
        XConeJetsCollection.fatjets.phi.push_back(fjet.phi_std());
        XConeJetsCollection.fatjets.mass.push_back(fjet.m());
        XConeJetsCollection.fatjets.area.push_back(fjet.area());
    }
    XConeJetsCollection.fatjets.n_jets = XConeJetsCollection.fatjets.pt.size();
    XConeJetsCollection.fatjets.R = RJets;
    XConeJetsCollection.fatjets.beta = BetaJets;
    XConeJetsCollection.fatjets.ptmin = 10;
    XConeJetsCollection.fatjets.etamax = 999;
    XConeJetsCollection.fastjet_CandsList = list_fat;

    //Check if subjets should be clustered
    bool doSubjets = true;
    if(NSubJets == 0) doSubjets = false;

    //For lepton specific clustering, we find the fatjet closest to the lepton,
    //and then use N-1 on its constituents
    uint lepton_jet_ind = 999;
    if (doLeptonSpecific && foundLepton){
        float dr_min = 999;
        for (uint iFj = 0; iFj < fatjets.size(); iFj++){
            const auto &fj = fatjets[iFj];
            float this_dr = deltaR(fj, lepton);
            if (this_dr < dr_min && this_dr < DRLeptonJet_){
                dr_min = this_dr;
                lepton_jet_ind = iFj;
            }
        }

        //Fill findLepton vector for fatjets
        for(uint iFj=0; iFj<fatjets.size(); iFj++){
            if(iFj == lepton_jet_ind) XConeJetsCollection.fatjets.findLepton.push_back(true);
            else XConeJetsCollection.fatjets.findLepton.push_back(false);
        }
    }

    // XConeJetsCollection.fatjets.findLepton[lepton_jet_ind] = true;

    //Return if the fatjet where the lepton is not clusted has a pt smaller than 400GeV
    // if (doLeptonSpecific && lepton_jet_ind != 999){
    //     for(uint i=0; i<fatjets.size(); i++){
    //         if(i == lepton_jet_ind) continue;
    //         if(fatjets[i].pt() < ptminJets){
    //             XConeReclusteredJets XConeJetsCollection;
    //             return XConeJetsCollection;
    //         }
    //     }
    // }
    //For baseline selection, at least one of the fatjets should have a pt larger than 200 GeV
    bool atLeastOneFastjet = false;
    for(uint i=0; i<fatjets.size(); i++){
        if(fatjets[i].pt() > ptminJets) atLeastOneFastjet = true;
    }
    if(!atLeastOneFastjet){
        XConeReclusteredJets XConeJetsCollection;
        return XConeJetsCollection;
    }

    vector<int> subjet_list, subjet_list2, topjet_list;
    subjet_list.clear();
    subjet_list2.clear();
    topjet_list.clear();
    //Loop over all fatjets and cluster subjets
    for(unsigned int i=0; i<fatjets.size(); i++){
        if(debug) std::cout << "\tNumero de fatjet: " << i << std::endl;

        //Get particles in fatjet
        std::vector<fastjet::PseudoJet> particle_in_fatjet;
        for (unsigned int ipart=0; ipart < _psj.size(); ipart++){
            if ( list_fat[ipart] < 0 ) continue; //If particle is not clustered in any jet
            else if (abs(list_fat[ipart]) == i){
                particle_in_fatjet.push_back(_psj.at(ipart));
            }
        }


        auto thisNSubJets = NSubJets;
        if (doLeptonSpecific && (i == lepton_jet_ind)){
            --thisNSubJets;
        }


        //Check if there are more particles than required subjets
        bool enoughParticles = (particle_in_fatjet.size() >= thisNSubJets);
        if(debug) std::cout << "\t\tNumero de particulas en el fatjet: " << particle_in_fatjet.size() << std::endl;
        if(debug) std::cout << "\t\tNumero de subjets: " << thisNSubJets << std::endl;
        if (!enoughParticles){
            XConeReclusteredJets XConeJetsCollection;
            // // XConeJetsCollection.clear();
            return XConeJetsCollection;
        }



        //Run second clustering step (subjets) for each fat jet
        std::vector<fastjet::PseudoJet> subjets;
        std::unique_ptr<fastjet::ClusterSequenceArea> clust_seq_sub;
        vector<int> list_sub;
        if(enoughParticles && doSubjets){
            std::unique_ptr<NjettinessPlugin> plugin_xcone_sub;
            initPlugin(plugin_xcone_sub, thisNSubJets, RSubJets, BetaSubJets, usePseudoXCone);
            JetDefinition jet_def_sub(plugin_xcone_sub.get());
            // Add GhostedAreaSpec for jet area calculation
            fastjet::GhostedAreaSpec ghost_spec_sub(5.0);  // eta_max = 5.0
            clust_seq_sub.reset(new fastjet::ClusterSequenceArea(particle_in_fatjet, jet_def_sub, ghost_spec_sub));
            Selector selector_subjets;
            if(i != lepton_jet_ind) selector_subjets = SelectorPtMin(ptminSubJets) && SelectorAbsRapMax(etamaxSubJets);
            else selector_subjets = SelectorPtMin(0.);
            subjets = sorted_by_pt(selector_subjets(clust_seq_sub->inclusive_jets()));
            list_sub.clear();
            list_sub = clust_seq_sub->particle_jet_indices(subjets);

            //Print some info about the subjets like number, pt, eta, phi, mass, area, etc.
            // std::cout << "\tNumber of subjets clustered in fatjet " << i << ": " << subjets.size() << std::endl;
            // for (unsigned int j = 0; j < subjets.size(); j++){
            //     const auto &sj = subjets[j];
            //     std::cout << "\t\tSubjet " << j << ": pt = " << sj.pt() << ", eta = " << sj.eta() << ", phi = " << sj.phi_std() << ", mass = " << sj.m() << ", area = " << sj.area() << std::endl;
            // }
        }


//         for(int s = 0; s < list_sub.size(); s++){
//             if(list_sub[s]==-1) continue
//             list_sub[s] = list_sub[s] + 3*i;
//         }
//         list_global_sub.insert(list_global_sub.end(),std::begin(list_sub),std::end(list_sub));


        for(const auto &subjet: subjets){
            XConeJetsCollection.subjets.pt.push_back(subjet.pt());
            XConeJetsCollection.subjets.eta.push_back(subjet.eta());
            XConeJetsCollection.subjets.phi.push_back(subjet.phi_std());
            XConeJetsCollection.subjets.mass.push_back(subjet.m());
            XConeJetsCollection.subjets.area.push_back(subjet.area());
        }
        XConeJetsCollection.subjets.n_jets = XConeJetsCollection.subjets.pt.size();
        XConeJetsCollection.subjets.R = RSubJets;
        XConeJetsCollection.subjets.beta = BetaSubJets;
        XConeJetsCollection.subjets.ptmin = ptminSubJets;
        XConeJetsCollection.subjets.etamax = etamaxSubJets;

        //For lepton specific clustering, we find the subjet closest to the lepton, I'm not sure if this is 100% correct
        uint lepton_subjet_ind = 999;
        if (doLeptonSpecific && foundLepton){
            float dr_min = 999;
            for (uint iSj = 0; iSj < subjets.size(); iSj++){
                if(lepton_jet_ind!=i) continue; //We are in the hadronic fatjet
                const auto &sj = subjets[iSj];
                float this_dr = deltaR(sj, lepton);
                if (this_dr < dr_min && this_dr < DRLeptonJet_){
                    dr_min = this_dr;
                    lepton_subjet_ind = iSj;
                }
            }

            //Fill findLepton vector for subjets
            for(uint iSj=0; iSj<subjets.size(); iSj++){
                if(iSj == lepton_subjet_ind) XConeJetsCollection.subjets.findLepton.push_back(true);
                else XConeJetsCollection.subjets.findLepton.push_back(false);
            }
        }

        //Return if the subjets belonging to the fatjet where the lepton is not clusted have a pt smaller than 30GeV or abs(eta) larger than 2.5
        // if (doLeptonSpecific && lepton_jet_ind!=i){
        //     for(uint j=0; j<subjets.size(); j++){
        //         if(subjets[j].pt() < ptminSubJets || abs(subjets[j].eta()) > etamaxSubJets){
        //             XConeReclusteredJets XConeJetsCollection;
        //             return XConeJetsCollection;
        //         }
        //     }

        // }




        if(i==0){
            subjet_list = list_sub;
        }else{
            subjet_list2 = list_sub;
        }


        //Check we got the number of subjets we asked for
        if(doSubjets && subjets.size() != thisNSubJets){
            if(debug) std::cout << "Only found " << subjets.size() << " subjets but requested " << thisNSubJets << ". "
                      << "Fatjet had " << particle_in_fatjet.size() << " constituents." << std::endl;
            // XConeReclusteredJets XConeJetsCollection;
            // // XConeJetsCollection.clear();
            // return XConeJetsCollection;
            if(subjets.size()==0){
                XConeReclusteredJets XConeJetsCollection;
                return XConeJetsCollection;
            }
        }


        //Get particles from subjets
        if(enoughParticles && doSubjets){
            std::vector<fastjet::PseudoJet> particles_from_subjets;
            particles_from_subjets.clear();
            for (uint s=0; s<subjets.size();s++){
                for(const auto &particle: subjets[s].constituents()){
                    particles_from_subjets.push_back(particle);
                }
                indices[i].push_back(XConeJetsCollection.subjets.pt.size());
            }
        }

        if(debug) std::cout << "LLEGO AQUI" << std::endl;

        
        //Topjets
        //We will save here the jet with the hadronic W boson
        //which we will use to get the mass of the top quark.
        //The topjets objetcs are just the sum of the corresponding subjets
        if(lepton_jet_ind!=i && lepton_jet_ind!=999){//hadronic decay of W boson scenario
            ROOT::Math::PtEtaPhiMVector topjet_p4(0,0,0,0);
            for (const auto &subjet: subjets){
                topjet_p4 += ROOT::Math::PtEtaPhiMVector(subjet.pt(), subjet.eta(), subjet.phi_std(), subjet.m());
                XConeJetsCollection.topjets.subjets_in_topjet.push_back(true);
            }
            XConeJetsCollection.topjets.pt = (topjet_p4.Pt());
            XConeJetsCollection.topjets.eta = (topjet_p4.Eta());
            XConeJetsCollection.topjets.phi = (topjet_p4.Phi());
            XConeJetsCollection.topjets.mass = (topjet_p4.M());
            XConeJetsCollection.topjets.n_subjets = subjets.size();
            XConeJetsCollection.topjets.ptmin = ptminJets;
            XConeJetsCollection.topjets.etamax = etamaxJets;

            //Fill here the W branches of the topjets.
            //Find the two subjets that are closer to the W mass
            //If there are less than two subjets, we will define that as the W
            float mass_W = 80.385;
            float mass_diff = 999.;
            int subjet1 = 0;
            int subjet2 = -1;
            ROOT::Math::PtEtaPhiMVector W_p4;
            for(int iSj=0; iSj<subjets.size(); iSj++){
                if(subjets.size()<2) break;
                for(int jSj=iSj+1; jSj<subjets.size(); jSj++){
                    ROOT::Math::PtEtaPhiMVector tmp_W_p4 = ROOT::Math::PtEtaPhiMVector(subjets[iSj].pt(), subjets[iSj].eta(), subjets[iSj].phi_std(), subjets[iSj].m()) + ROOT::Math::PtEtaPhiMVector(subjets[jSj].pt(), subjets[jSj].eta(), subjets[jSj].phi_std(), subjets[jSj].m());
                    if(abs(tmp_W_p4.M()-mass_W)<mass_diff){
                        mass_diff = abs(tmp_W_p4.M()-mass_W);
                        W_p4 = tmp_W_p4;
                        subjet1 = iSj;
                        subjet2 = jSj;
                    }
                }
            }
            if(debug) std::cout << "Subjet1: " << subjet1 << " Subjet2: " << subjet2 << std::endl;
            if(subjet2==-1){//There is only one subjet for the hadronic decay of W, and that jet will be consider as the W (You can apply cuts on the number of subjets in the hadronic decay)
                XConeJetsCollection.topjets.pt_W = (subjets[subjet1].pt());
                XConeJetsCollection.topjets.eta_W = (subjets[subjet1].eta());
                XConeJetsCollection.topjets.phi_W = (subjets[subjet1].phi_std());
                XConeJetsCollection.topjets.mass_W = (subjets[subjet1].m());
                XConeJetsCollection.topjets.subjets_in_W.push_back(true);
            }else{
                XConeJetsCollection.topjets.pt_W = W_p4.Pt();
                XConeJetsCollection.topjets.eta_W = W_p4.Eta();
                XConeJetsCollection.topjets.phi_W = W_p4.Phi();
                XConeJetsCollection.topjets.mass_W = W_p4.M();
                for(int iSj=0; iSj<subjets.size(); iSj++){
                    if(iSj == subjet1 || iSj == subjet2) XConeJetsCollection.topjets.subjets_in_W.push_back(true);
                    else XConeJetsCollection.topjets.subjets_in_W.push_back(false);
                }
            }

        }else{//leptonic decay of W boson scenario
            for(int iSj=0; iSj<subjets.size(); iSj++){
                XConeJetsCollection.topjets.subjets_in_topjet.push_back(false);
                XConeJetsCollection.topjets.subjets_in_W.push_back(false);
            }
        }
//         std::unique_ptr<NjettinessPlugin> plugin_xcone_top;
//         initPlugin(plugin_xcone_top, 1, 1.2, BetaJets, usePseudoXCone);
//         JetDefinition jet_def_top(plugin_xcone_top.get());
//         ClusterSequence clust_seq_top(particles_from_subjets, jet_def_top);
// //         Selector selector_topjets = SelectorPtMin(ptminJets) && SelectorAbsRapMax(etamaxJets);
// //         std::vector<fastjet::PseudoJet> topjets = sorted_by_pt(selector_topjets(clust_seq_top.inclusive_jets()));
//         std::vector<fastjet::PseudoJet> topjets = sorted_by_pt(clust_seq_top.inclusive_jets(0));



//         if(i==0){
//             topjet_list = clust_seq_top.particle_jet_indices(topjets);
//         }else{
//             topjet_list2 = clust_seq_top.particle_jet_indices(topjets);
//         }

//         for (const auto &topjet: topjets){
//             XConeJetsCollection.topjets.pt.push_back(topjet.pt());
//             XConeJetsCollection.topjets.eta.push_back(topjet.eta());
//             XConeJetsCollection.topjets.phi.push_back(topjet.phi_std());
//             XConeJetsCollection.topjets.mass.push_back(topjet.m());
//         }

//        /* //Adding the subjets properties from subjets
//         float sum_pt = 0.0, sum_mass = 0.0;
//         float weighted_eta = 0.0, weighted_phi = 0.0;

//         for (const auto& subjet : subjets) {
//             float pt = subjet.pt();
//             float eta = subjet.eta();
//             float phi = subjet.phi();
//             float mass = subjet.m();

//             sum_pt += pt;
//             sum_mass += mass;
//             weighted_eta += pt * eta;
//             weighted_phi += pt * phi;
//         }

//         float eta_topjet = (sum_pt > 0) ? (weighted_eta / sum_pt) : 0.0;
//         float phi_topjet = (sum_pt > 0) ? (weighted_phi / sum_pt) : 0.0;

//         XConeJetsCollection.topjets.pt.push_back(sum_pt);
//         XConeJetsCollection.topjets.eta.push_back(eta_topjet);
//         XConeJetsCollection.topjets.phi.push_back(phi_topjet);
//         XConeJetsCollection.topjets.mass.push_back(sum_mass);*/



//         XConeJetsCollection.topjets.n_jets = XConeJetsCollection.topjets.pt.size();
//         XConeJetsCollection.topjets.R = 1.2;
//         XConeJetsCollection.topjets.beta = BetaJets;
//         XConeJetsCollection.topjets.ptmin = ptminJets;
//         XConeJetsCollection.topjets.etamax = etamaxJets;

//         //Check we got the number of topjets we asked for
//         if(topjets.size() != 1){
//             //std::cout << "Found " << topjets.size() << " jets but only requested one top-jet.\n"
//                         << "Added at least one jet." << std::endl;
//         }
    }//loop over fatjets

    //Return if the invariant mass of the fatjet where the lepton is clustered is greater than the invariant mass of the fatjet where the lepton is not clustered
    // if (doLeptonSpecific && foundLepton){
    //     if(lepton_jet_ind == 0 && fatjets[lepton_jet_ind].m() > fatjets[1].m()){
    //         XConeReclusteredJets XConeJetsCollection;
    //         return XConeJetsCollection;
    //     }else if(lepton_jet_ind == 1 && fatjets[lepton_jet_ind].m() > fatjets[0].m()){
    //         XConeReclusteredJets XConeJetsCollection;
    //         return XConeJetsCollection;
    //     }
    // }


//     vector<int> list_global_sub;
//     list_global_sub.clear();
    int s1=0, s2=0;
    for(int i = 0; i<XConeJetsCollection.fastjet_CandsList.size(); i++){
        if(XConeJetsCollection.fastjet_CandsList[i]==-1){
            if(debug) std::cout << "\t\t\tNumero PFCand: " << i << std::endl;
            XConeJetsCollection.subjet_CandsList.push_back(-1);
        }else if(XConeJetsCollection.fastjet_CandsList[i]==0){
            if(debug) std::cout << "\t\t\tNumero PFCand(2): " << i << std::endl;
            XConeJetsCollection.subjet_CandsList.push_back(subjet_list[s1]);
            s1++;
        }else if(XConeJetsCollection.fastjet_CandsList[i]==1){
            if(subjet_list2[s2]==-1) {
                if(debug) std::cout << "\t\t\tNumero PFCand(3): " << i << std::endl;
                XConeJetsCollection.subjet_CandsList.push_back(-1);
            }else{
                if(debug) std::cout << "\t\t\tNumero PFCand(4): " << i << std::endl;
                XConeJetsCollection.subjet_CandsList.push_back(subjet_list2[s2] + 1 + *std::max_element(subjet_list.begin(), subjet_list.end()));
            }
            s2++;
        }

    }
    // s1=0, s2=0;
    // for(int i = 0; i<XConeJetsCollection.subjet_CandsList.size(); i++){
    //     if(XConeJetsCollection.subjet_CandsList[i]==-1){
    //         XConeJetsCollection.topjet_list.push_back(-1);
    //     }else if(XConeJetsCollection.subjet_CandsList[i]>-1 && XConeJetsCollection.subjet_CandsList[i]<3){
    //         XConeJetsCollection.topjet_list.push_back(topjet_list[s1]);
    //         s1++;
    //     }else if(XConeJetsCollection.subjet_CandsList[i]>2 && XConeJetsCollection.subjet_CandsList[i]<6){
    //         if(topjet_list2[s2]==-1){
    //             XConeJetsCollection.topjet_list.push_back(-1);
    //         }else{
    //             XConeJetsCollection.topjet_list.push_back(topjet_list2[s2]+1);
    //         }
    //         s2++;
    //     }
    // }


//     vector<int> list_global_


//     XConeJetsCollection.subjet_CandsList = list_global_sub;

    return XConeJetsCollection;


}







#endif