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

// Inicializar `cset` al comienzo del header
inline void initializeCorrectionSet() {
    fs::path fname_ak4_2023 = "./jet_jerc.json.gz";
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

// Inicializar `cset` al comienzo del header
inline void initializeBTagCorrectionSet() {
    fs::path fname_sf_btag_2023 = "./btagging.json.gz"; //"/eos/user/c/cmunozdi/SWAN_projects/BtagScaleFactors/btagging.json.gz"; 
    fs::path fname_eff_btag_2023 = "./btag_efficiencies_combined.json"; //"/eos/project/r/rtu-topanalysis/AnalysisSamples_JetTightIDNoLepVeto/output_efficiencies_rdf/btag_efficiencies_combined.json"; 
    // fs::path fname_btag_2023 = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2023_Summer23/btagging.json.gz";
    if (!fs::exists(fname_sf_btag_2023)) {
        throw std::runtime_error("El archivo de corrección no existe: " + fname_sf_btag_2023.string());
    }
    sf_btagset = correction::CorrectionSet::from_file(fname_sf_btag_2023.string());
    std::cout << "Archivo de corrección cargado exitosamente: " << fname_sf_btag_2023.string() << std::endl;
    if (!fs::exists(fname_eff_btag_2023)) {
        throw std::runtime_error("El archivo de corrección no existe: " + fname_eff_btag_2023.string());
    }
    eff_btagset = correction::CorrectionSet::from_file(fname_eff_btag_2023.string());
    std::cout << "Archivo de corrección cargado exitosamente: " << fname_eff_btag_2023.string() << std::endl;
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

// Funtion to compute btagWeight of the event given the jet collection
inline float compute_btagWeight(const rvec_f &jet_pts, const rvec_f &jet_etas, const ROOT::VecOps::RVec<UChar_t>& jet_flavs, const rvec_f &jet_btags) {
    float weight = 1.;
    for (size_t i = 0; i < jet_pts.size(); i++) {

        if (jet_pts[i] < 30. || abs(jet_etas[i]) > 2.5 /*|| jet_pts[i] > 1000*/) continue;

        float sf = 1.;
        try {
            if (static_cast<int>(jet_flavs[i]) == 0) {
                sf = sf_btagset->at("deepJet_light")->evaluate(
                    std::vector<correction::Variable::Type>{"central", "T", static_cast<int>(jet_flavs[i]), abs(jet_etas[i]), jet_pts[i]}
                );
            } else {
                sf = sf_btagset->at("deepJet_comb")->evaluate(
                    std::vector<correction::Variable::Type>{"central", "T", static_cast<int>(jet_flavs[i]), abs(jet_etas[i]), jet_pts[i]}
                );
            }
        } catch (const std::exception &e) {
            
            std::cout << "Jet " << i << ": pt = " << jet_pts[i] << ", eta = " << jet_etas[i] 
                  << ", flav = " << static_cast<int>(jet_flavs[i]) << ", btag = " << jet_btags[i] << std::endl;
            std::cerr << "Error al evaluar el SF para el jet " << i << ": " << e.what() << ". Value of sf=" << sf << std::endl;
            continue;
        }

        float eff = 1.;
        try {
            eff = eff_btagset->at("btag_efficiency")->evaluate(
                std::vector<correction::Variable::Type>{jet_etas[i], jet_pts[i], static_cast<int>(jet_flavs[i])}
            );
        } catch (const std::exception &e) {
            
            std::cout << "Jet " << i << ": pt = " << jet_pts[i] << ", eta = " << jet_etas[i] 
                  << ", flav = " << static_cast<int>(jet_flavs[i]) << ", btag = " << jet_btags[i] << std::endl;
            std::cerr << "Error al evaluar la eficiencia para el jet " << i << ": " << e.what() << ". Value of eff=" << eff << std::endl;
            continue;
        }

        float jet_w = 1;
        if (jet_btags[i] > 0.4648) {
            jet_w = sf;
        } else {
            jet_w = (1 - sf * eff) / (1 - eff);
        }

        weight *= jet_w;
    }
    return weight;
}


// /**
//     @short select gen candidates (for now only charged)
// */
// rvec_b selectGenCandidates(const rvec_i &pdgId) {

//     std::vector<bool> isValidCand(pdgId.size(),false);

//     //select charged particles only
//     for(size_t i=0; i<pdgId.size(); i++) {

//         if( !Rivet::PID::isCharged(pdgId[i]) ) continue;

//         //... do other selections here

//         isValidCand[i]=0;
//     }

//     return rvec_b(isValidCand.begin(), isValidCand.end());
// }


/**
    @short runs fastjet on a collection of particles and returns jets with pt, eta, phi
*/

struct JetAntikTReclus {
    vector<float> pt;
    vector<float> eta;
    vector<float> phi;
    vector<float> mass;
    int n_jets;
    float R;
    float ptmin;
};

inline JetAntikTReclus buildJets(const ROOT::VecOps::RVec<float> &pt,
                  const ROOT::VecOps::RVec<float> &eta,
                  const ROOT::VecOps::RVec<float> &phi,
                  const ROOT::VecOps::RVec<int> &pdgId,
                  float R = .8, float ptmin = 400.) {

    using namespace fastjet;

    std::vector<PseudoJet> particles;
    for (size_t i = 0; i < pt.size(); i++) {
        float p_mass = get_mass(pdgId[i]);
        if (p_mass < 0) continue;
        ROOT::Math::PtEtaPhiMVector p4(pt[i], eta[i], phi[i], p_mass);
        particles.emplace_back(p4.px(), p4.py(), p4.pz(), p4.energy());
    }

    JetDefinition jet_def(antikt_algorithm, R);
    ClusterSequence cs(particles, jet_def);
    Selector selector_jets = SelectorPtMin(ptmin) && SelectorAbsRapMax(2.5);
    auto jets = sorted_by_pt(selector_jets(cs.inclusive_jets()));

    JetAntikTReclus jet_data;
    for (const auto &jet : jets) {
        jet_data.pt.push_back(jet.pt());
        jet_data.eta.push_back(jet.eta());
        jet_data.phi.push_back(jet.phi_std());
        jet_data.mass.push_back(jet.m());
    }
    jet_data.n_jets = jet_data.pt.size();
    jet_data.R = R;
    jet_data.ptmin = ptmin;
//     //std::cout << "Number of particles: " << jet_data.pt.size() << std::endl;

    return jet_data;
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
                                const unique_ptr<correction::CorrectionSet>& cset) {
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

    // Then the conditions to apply the correction
    string jec = "Summer23Prompt23_V2_MC",
           lvl = "L1L2L3Res",
           algo = "AK4PFPuppi";

    // Getting the correction factor
    double jet_corrFactor=0.;
    if(jet_Area==0 && rho_PU==0) jet_corrFactor = 1.;
    else jet_corrFactor = singleLevelEnergyCorr(jet_raw_map, cset, jec, lvl, algo);

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

inline float compute_ptrel(const TLorentzVector& muon,
                            const TLorentzVector& jet_corr,
                            float jet_rawFactor, float& dR, bool muonBelongsToJet = false) {
    // Jet with no correction
    TLorentzVector jet_raw = jet_corr * (1.0 - jet_rawFactor);


    // std::cout << "Muon pt: " << muon.Pt() << ", eta: " << muon.Eta() << ", phi: " << muon.Phi() << std::endl;
    // std::cout << "Jet corr pt: " << jet_corr.Pt() << ", eta: " << jet_corr.Eta() << ", phi: " << jet_corr.Phi() << ", rawFactor: " << jet_rawFactor << std::endl;
    // std::cout << "Jet raw pt: " << jet_raw.Pt() << ", eta: " << jet_raw.Eta() << ", phi: " << jet_raw.Phi() << std::endl;
    // std::cout << "muon.DeltaR(jet_corr): " << muon.DeltaR(jet_corr) << std::endl;
    // std::cout << "deltaR(muon, jet_corr): " << deltaR(muon.Eta(), muon.Phi(), jet_corr.Eta(), jet_corr.Phi())  << std::endl;

    // If deltaR(lep, jet)>0.4, do not sustract the muon from the jet
    if (muon.DeltaR(jet_corr) > 0.4 || !muonBelongsToJet) {
        dR = muon.DeltaR(jet_corr);
        // std::cout << "\t DeltaR is greater than 0.4 (" << dR << "), returning ptrel: " << muon.Vect().Cross(jet_corr.Vect().Unit()).Mag() << std::endl;
        return muon.Vect().Cross(jet_corr.Vect().Unit()).Mag();
    }

    // Verify that the jet and the lepton are not the same vector
    // if (muon.DeltaR(jet_raw) < 1e-4) {
    //     //std::cout << "\t Jet and muon are the same vector, returning -1" << std::endl;
    //     return -1;  // will look for the next closest jet
    // }

    //Uncorrected the muon (it was embedded in the raw jet)
    float l1corrFactor;
    if(abs(1.0-jet_rawFactor)>1e-5){
        l1corrFactor = 1.0 / (1.0 - jet_rawFactor);
    }else{
        l1corrFactor = 1.0; // If the jet_rawFactor is very close to 1, we assume no correction is needed
    }
    TLorentzVector muon_raw = muon;// * (1.0 / l1corrFactor);

    //Jet without the muon, in RAW scale
    TLorentzVector jet_raw_no_mu = jet_raw - muon_raw;

    // Jet corrected without the muon
    float scale = jet_corr.Pt() / jet_raw.Pt();  // same factor applied to the whole jet
    TLorentzVector jet_corr_no_mu = jet_raw_no_mu * scale;

    // Get the dR between the muon and the corrected jet without the muon
    dR = muon.DeltaR(jet_corr_no_mu);

    // std::cout << "\t l1corrFactor: " << l1corrFactor << std::endl;
    // std::cout << "\t muon_raw_pt: " << muon_raw.Pt() << ", eta: " << muon_raw.Eta() << ", phi: " << muon_raw.Phi() << std::endl;
    // std::cout << "\t jet_raw_no_mu_pt: " << jet_raw_no_mu.Pt() << ", eta: " << jet_raw_no_mu.Eta() << ", phi: " << jet_raw_no_mu.Phi() << std::endl;
    // std::cout << "\t scale: " << scale << std::endl;
    // std::cout << "\t jet_corr_no_mu_pt: " << jet_corr_no_mu.Pt() << ", eta: " << jet_corr_no_mu.Eta() << ", phi: " << jet_corr_no_mu.Phi() << std::endl;
    // std::cout << "\t dR(muon, jet_corr_no_mu): " << dR << std::endl;

    // If the corrected jet without the muon has a pt lower than 15 GeV, return -1 (indicating no valid ptrel)
    if (jet_corr_no_mu.Pt() < 15) {
        return -1.0;
    }

    // Final jet = corrected jet without the muon + corrected muon
    // with this, we can easily get the corrected jet without the muon by sustracting the muon
    TLorentzVector jet_final = jet_corr_no_mu + muon;

   // Get the dR between the muon and the corrected jet without the muon
    dR = muon.DeltaR(jet_corr_no_mu);


    // std::cout << "\t jet_final_pt: " << jet_final.Pt() << ", eta: " << jet_final.Eta() << ", phi: " << jet_final.Phi() << std::endl;
    // std::cout << "\t dR(muon, jet_final): " << dR << std::endl;

    // pTrel
    float ptrel = muon.Vect().Cross((/*jet_final - muon*/jet_corr_no_mu).Vect().Unit()).Mag();

    // std::cout << "\t ptrel: " << ptrel << std::endl;

    return ptrel;
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

inline Lepton triggerLepton(const rvec_f &pt_le, const rvec_f &eta_le, const rvec_f &phi_le, const rvec_i &pdgId_le, const rvec_i &jetIdx_le, const rvec_b &tightId_le, /*const rvec_f &ak4_lesubDeltaeta, const rvec_f &ak4_lesubDeltaphi, const rvec_f &ak4_lesubfactor,*/ const rvec_f &PFCands_pt, const rvec_f &PFCands_eta, const rvec_f &PFCands_phi, const rvec_i &PFCands_pdgId, const rvec_f &ak4_pt, const rvec_f &ak4_eta, const rvec_f &ak4_phi, const rvec_f &ak4_mass, const rvec_b &ak4_passJetIdTightLepVeto, const rvec_f &ak4_rawFactor, const rvec_f &ak4_area, const float &rho_PU, float ak4_pt_threshold = 25.0, const bool isMuon = true, bool AnalysisCuts = false, bool GenLevel = false, const float thr_dR = 0.3, const float thr_ptrel = 30.0) {
    // std::cout << "Is GenLevel: " << GenLevel << std::endl;
    // std::cout << "Trigger Lepton selection: isMuon = " << isMuon << ", AnalysisCuts = " << AnalysisCuts << std::endl;
    // std::cout << "Jet threshold: " << ak4_pt_threshold << std::endl;
    ROOT::VecOps::RVec<bool> tightId_le_local;// = tightId_le; // Crear una copia modificable

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

                    tmp_ptrel = compute_ptrel_new(lepton_p4, jet_corr_p4, dR, ak4_pt_threshold, ak4_rawFactor[j], ak4_area[j], rho_PU, cset);
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
                // std::cout << "#########dR: " << min_dR_to_jet << ", ptrel: " << min_pt_rel_to_jet << ", ptrelNANO: " << le_jetptrel[i] << "##########" << endl << std::endl;
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

    return combined;
}

//Funtion to calculate the missing transverse momentum (if we return -sqrt((px)^2+(py)^2) we get the missing energy, just in the opposite direction of the rest of PFCands)
inline double Get_pTmiss (const rvec_f &pt, const rvec_f &phi){
    double px = 0.0;
    double py = 0.0;
    for(size_t i=0; i<pt.size(); i++){
        px += pt[i]*cos(phi[i]);
        py += pt[i]*sin(phi[i]);
    }
    return sqrt(px*px + py*py);

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



/*BUILDING XCONE RECLUSTERING FUNCTIONS*/
inline void initPlugin(std::unique_ptr<NjettinessPlugin> & ptr, int N, float R0, float beta, bool usePseudoXCone){
    if(usePseudoXCone){
        ptr.reset(new PseudoXConePlugin(N, R0, beta));
    } else {
        ptr.reset(new XConePlugin(N, R0, beta));
    }
}

inline XConeReclusteredJets buildXConeJets(const rvec_f &pt, const rvec_f &eta, const rvec_f &phi, const rvec_f &mass, const rvec_i &pdgId, int NJets = 2, float RJets = 1.2, float BetaJets = 2.0, float ptminJets = 200., float etamaxJets = 3., int NSubJets = 3, float RSubJets = 0.4, float BetaSubJets = 2., float ptminSubJets = 22.5, float etamaxSubJets = 3., bool doLeptonSpecific = true, bool usePseudoXCone = false){

    // const reco::Candidate *lepton(nullptr);
    std::vector<fastjet::PseudoJet>::iterator lepton_iter;// = _psj.end();
    float lepton_min_pt = 00;
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
            if ((pdgid == 11 || pdgid == 13) && (candPt > lepton_min_pt)){
                // lepton_iter = std::prev(_psj.end());
                lepton_index = i_gl;
                lepton_min_pt = candPt;
                lepton = _psj.back();
            }
        }//doLeptonSpecific
    }//for loop over pfcands of the event

    bool foundLepton = false;

    if (doLeptonSpecific && lepton != 0.){
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
    ClusterSequence clust_seq_xcone(_psj, jet_def_xcone);
    Selector selector_fatjets = SelectorPtMin(10); //pt larger than 10GeV (this is the loosest criteria, for the jet where the lepton is clustered; the main fatjet would have to be larger than 400GeV)
    std::vector<fastjet::PseudoJet> fatjets = sorted_by_pt(selector_fatjets(clust_seq_xcone.inclusive_jets()));



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
        std::unique_ptr<ClusterSequence> clust_seq_sub;
        vector<int> list_sub;
        if(enoughParticles && doSubjets){
            std::unique_ptr<NjettinessPlugin> plugin_xcone_sub;
            initPlugin(plugin_xcone_sub, thisNSubJets, RSubJets, BetaSubJets, usePseudoXCone);
            JetDefinition jet_def_sub(plugin_xcone_sub.get());
            clust_seq_sub.reset(new ClusterSequence(particle_in_fatjet, jet_def_sub));
            Selector selector_subjets;
            if(i != lepton_jet_ind) selector_subjets = SelectorPtMin(ptminSubJets) && SelectorAbsRapMax(etamaxSubJets);
            else selector_subjets = SelectorPtMin(0.);
            subjets = sorted_by_pt(selector_subjets(clust_seq_sub->inclusive_jets()));
            list_sub.clear();
            list_sub = clust_seq_sub->particle_jet_indices(subjets);
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
        if (doLeptonSpecific && lepton_jet_ind!=i){
            for(uint j=0; j<subjets.size(); j++){
                if(subjets[j].pt() < ptminSubJets || abs(subjets[j].eta()) > etamaxSubJets){
                    XConeReclusteredJets XConeJetsCollection;
                    return XConeJetsCollection;
                }
            }

        }




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
    if (doLeptonSpecific && foundLepton){
        if(lepton_jet_ind == 0 && fatjets[lepton_jet_ind].m() > fatjets[1].m()){
            XConeReclusteredJets XConeJetsCollection;
            return XConeJetsCollection;
        }else if(lepton_jet_ind == 1 && fatjets[lepton_jet_ind].m() > fatjets[0].m()){
            XConeReclusteredJets XConeJetsCollection;
            return XConeJetsCollection;
        }
    }


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