#ifndef __nanoPF_helpers_h__
#define __nanoPF_helpers_h__

#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"
#include <vector>

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

using namespace ROOT::VecOps;
using namespace fastjet;
using namespace contrib;
using namespace std;

using rvec_f = const RVec<float>;
using rvec_i = const RVec<int>;
using rvec_b = const RVec<bool>;

bool debug = false;

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


/**
    @short select gen candidates (for now only charged)
*/
rvec_b selectGenCandidates(const rvec_i &pdgId) {

    std::vector<bool> isValidCand(pdgId.size(),false);

    //select charged particles only
    for(size_t i=0; i<pdgId.size(); i++) {

        if( !Rivet::PID::isCharged(pdgId[i]) ) continue;

        //... do other selections here

        isValidCand[i]=0;
    }

    return rvec_b(isValidCand.begin(), isValidCand.end());
}


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

JetAntikTReclus buildJets(const ROOT::VecOps::RVec<float> &pt,
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
//     std::cout << "Number of particles: " << jet_data.pt.size() << std::endl;

    return jet_data;
}

struct Lepton{
    std::vector<float> pt;
    std::vector<float> eta;
    std::vector<float> phi;
    std::vector<int> pdgId;
    int n_lep=0;

};

Lepton triggerLepton(const rvec_f &pt_le, const rvec_f &eta_le, const rvec_f &phi_le, const rvec_i &pdgId_le, const rvec_f &ak4_eta, const rvec_f &ak4_phi, bool isMuon = true){
    Lepton lepton;
    if(isMuon){
        //Selection for muon: pt > 60, |eta| < 2.4. Also, has not to be inside a AK4 jet and has a pt_rel<40 being pt_rel the relative pt of the muon in the orthogonal direction to the jet
        if(pt_le.size()>0){
            for(size_t i=0; i<pt_le.size(); i++){
                if(pt_le[i]>60 && abs(eta_le[i])<2.4){
                    bool isInsideJet = false;
                    for(size_t j=0; j<ak4_eta.size(); j++){
                        if(deltaR(eta_le[i],ak4_eta[j],phi_le[i],ak4_phi[j])<0.4){
                            if(pt_le[i]*sin(deltaPhi(phi_le[i],ak4_phi[j]))<40) isInsideJet = true;
                        }
                    }
                    if(!isInsideJet){
                        lepton.pt.push_back(pt_le[i]);
                        lepton.eta.push_back(eta_le[i]);
                        lepton.phi.push_back(phi_le[i]);
                        lepton.pdgId.push_back(pdgId_le[i]);
                    }
                }
            }
        }
    }else{
        //Selection for electron: pt > 120, |eta| < 2.4. Also, has not to be inside a AK4 jet and has a pt_rel<40 being pt_rel the relative pt of the electron in the orthogonal direction to the jet
        if(pt_le.size()>0){
            for(size_t i=0; i<pt_le.size(); i++){
                if(pt_le[i]>120 && abs(eta_le[i])<2.4){
                    bool isInsideJet = false;
                    for(size_t j=0; j<ak4_eta.size(); j++){
                        if(deltaR(eta_le[i],ak4_eta[j],phi_le[i],ak4_phi[j])<0.4){
                            if(pt_le[i]*sin(deltaPhi(phi_le[i],ak4_phi[j]))<40) isInsideJet = true;
                        }
                    }
                    if(!isInsideJet){
                        lepton.pt.push_back(pt_le[i]);
                        lepton.eta.push_back(eta_le[i]);
                        lepton.phi.push_back(phi_le[i]);
                        lepton.pdgId.push_back(pdgId_le[i]);
                    }
                }
            }
        }
    }
    lepton.n_lep = lepton.pt.size();
    return lepton;
}

Lepton CombineLeptons(const Lepton& muons, const Lepton& electrons) {
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

    // Actualizar el número total de leptones
    combined.n_lep = muons.n_lep + electrons.n_lep;

    return combined;
}

//Funtion to calculate the missing transverse momentum (if we return -sqrt((px)^2+(py)^2) we get the missing energy, just in the opposite direction of the rest of PFCands)
double Get_pTmiss (const rvec_f &pt, const rvec_f &phi){
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
void initPlugin(std::unique_ptr<NjettinessPlugin> & ptr, int N, float R0, float beta, bool usePseudoXCone){
    if(usePseudoXCone){
        ptr.reset(new PseudoXConePlugin(N, R0, beta));
    } else {
        ptr.reset(new XConePlugin(N, R0, beta));
    }
}

XConeReclusteredJets buildXConeJets(const rvec_f &pt, const rvec_f &eta, const rvec_f &phi, const rvec_f &mass, const rvec_i &pdgId, int NJets = 2, float RJets = 1.2, float BetaJets = 2.0, float ptminJets = 400., float etamaxJets = 2.4, int NSubJets = 3, float RSubJets = 0.4, float BetaSubJets = 2., float ptminSubJets = 30., float etamaxSubJets = 2.5, bool doLeptonSpecific = true, bool usePseudoXCone = false){

    // const reco::Candidate *lepton(nullptr);
    std::vector<fastjet::PseudoJet>::iterator lepton_iter;// = _psj.end();
    float lepton_max_pt = 00;
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
            if ((pdgid == 11 || pdgid == 13) && (candPt > lepton_max_pt)){
                // lepton_iter = std::prev(_psj.end());
                lepton_index = i_gl;
                lepton_max_pt = candPt;
                lepton = _psj.back();
            }
        }//doLeptonSpecific
    }//for loop over pfcands of the event

    bool foundLepton = false;

    if (doLeptonSpecific && lepton != 0.){
        if(debug) std::cout << "PT del lepton es: " << lepton_max_pt << std::endl;
        // lepton = *lepton_iter;
        if(debug) std::cout << "PT del lepton es1: " << lepton.pt() << std::endl;
        foundLepton = true;
        if(debug) std::cout << "PT del lepton es: " << lepton_max_pt << std::endl;
    }

    // if(!foundLepton){
    //     std::cout << "foundLepton: " << foundLepton << std::endl;
    //     std::cout << "doLeptonSpecific: " << doLeptonSpecific << std::endl;
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
    if (doLeptonSpecific && lepton_jet_ind != 999){
        for(uint i=0; i<fatjets.size(); i++){
            if(i == lepton_jet_ind) continue;
            if(fatjets[i].pt() < ptminJets){
                XConeReclusteredJets XConeJetsCollection;
                return XConeJetsCollection;
            }
        }
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
            float mass_W = 80.3602;
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
//             std::cout << "Found " << topjets.size() << " jets but only requested one top-jet.\n"
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




/**
    @short checks for opposite charge, compatible flavor quark pairs for W boson candidates
*/
rvec_b dijetCands(const rvec_i &pdgId, const rvec_f &pt, const rvec_f &eta, const rvec_f &phi, float m0=80.4, float deltam=31.0)
{
    std::vector<bool> isValidCand(pdgId.size(), false);

    // Loop over the list of quarks
    for (size_t i = 0; i < pdgId.size(); i++) {

        // Select light quarks (u, d, s, c)
        if (!(abs(pdgId[i]) == 1 || abs(pdgId[i]) == 2 || abs(pdgId[i]) == 3 || abs(pdgId[i]) == 4)) continue;

        // Try to make a pair with charge +-1
        for (size_t j = i + 1; j < pdgId.size(); j++) {

            // Ensure the second quark is also light
            if (!(abs(pdgId[j]) == 1 || abs(pdgId[j]) == 2 || abs(pdgId[j]) == 3 || abs(pdgId[j]) == 4)) continue;

            // Require opposite pdgId sign (q+anti-q)
            if (pdgId[i] * pdgId[j] > 0) continue;

            //assign quark mass
            float mass1(get_mass(pdgId[i]));
//             if(abs(pdgId[i])==1) mass1=0.00470;
//             if(abs(pdgId[i])==2) mass1=0.00216;
//             if(abs(pdgId[i])==3) mass1=0.09350;
//             if(abs(pdgId[i])==4) mass1=1.27300;
            if(mass1<0) continue;

            float mass2(get_mass(pdgId[j]));
//             if(abs(pdgId[j])==1) mass2=0.00470;
//             if(abs(pdgId[j])==2) mass2=0.00216;
//             if(abs(pdgId[j])==3) mass2=0.09350;
//             if(abs(pdgId[j])==4) mass2=1.27300;
            if(mass1<0) continue;

            // Compute mass of the quark pair system
            ROOT::Math::PtEtaPhiMVector pi(pt[i], eta[i], phi[i], mass1);
            ROOT::Math::PtEtaPhiMVector pj(pt[j], eta[j], phi[j], mass2);
            float mqq = (pi + pj).M();

            // Check compatibility with the W mass window
            if (fabs(mqq - m0) > deltam) continue;

            // Mark the pair as valid
            isValidCand[i] = true;
            isValidCand[j] = true;
        }
    }

    return rvec_b(isValidCand.begin(), isValidCand.end());
}



/**
    @short checks if a set of objects are isolated with respect to a reference in the eta-phi plane
*/
rvec_b crossClean(const rvec_f &eta,const rvec_f &phi, const rvec_f &eta_ref, const rvec_f &phi_ref,float cone=0.4)
{
    std::vector<bool> isIso;

    //loop over the first list of objects
    for(size_t i=0; i<eta.size(); i++) {

        float minDR(9999.);
        for(size_t j=0; i<eta_ref.size(); i++) {
            minDR = min(minDR,ROOT::VecOps::DeltaR(eta[i],eta_ref[j],phi[i],phi_ref[j]));
        }
        isIso.push_back( minDR>cone );
    }

    return rvec_b(isIso.begin(), isIso.end());
}

/**
   @returns a kinematics feature of a two body system
*/
float kinematics2lq(const int &pdgId1, const float &pt1, const float &eta1, const float &phi1,
                   const int &pdgId2, const float &pt2, const float &eta2, const float &phi2,
                   std::string kin="mass")
{
    float m1(get_mass(pdgId1)), m2(get_mass(pdgId2));
//     if(abs(pdgId1)==1) m1 = 0.00470;
//     else if(abs(pdgId1)==2) m1 = 0.00216;
//     else if(abs(pdgId1)==3) m1 = 0.09350;
//     else if(abs(pdgId1)==4) m1 = 1.27300;

//     if(abs(pdgId2)==1) m2 = 0.00470;
//     else if(abs(pdgId2)==2) m2 = 0.00216;
//     else if(abs(pdgId2)==3) m2 = 0.09350;
//     else if(abs(pdgId2)==4) m2 = 1.27300;

    ROOT::Math::PtEtaPhiMVector p1(pt1, eta1, phi1, m1);
    ROOT::Math::PtEtaPhiMVector p2(pt2, eta2, phi2, m2);
    if(kin=="pt") return (p1+p2).Pt();
    if(kin=="eta") return (p1+p2).Eta();
    if(kin=="phi") return (p1+p2).Phi();
    return (p1+p2).M();
}

#include <vector>
#include "ROOT/RVec.hxx"

// struct Lepton {
//     float pt;
//     float phi;
//     float eta;
//     int pdgId;
// };

// std::vector<Lepton> EventSelectionLepton(
//     const ROOT::VecOps::RVec<float>& lepton_pt,
//     const ROOT::VecOps::RVec<float>& lepton_eta,
//     const ROOT::VecOps::RVec<float>& lepton_phi,
//     const ROOT::VecOps::RVec<int>& lepton_pdgId,
//     const ROOT::VecOps::RVec<int>& muon_mask,
//     const ROOT::VecOps::RVec<float>& electron_pt,
//     const ROOT::VecOps::RVec<float>& electron_phi,
//     const ROOT::VecOps::RVec<float>& electron_eta,
//     const ROOT::VecOps::RVec<int>& electron_pdgId,
//     const ROOT::VecOps::RVec<int>& electron_mask) {

//     std::vector<Lepton> leptonInfo;

//     // Añadir información de muones seleccionados
//     for (size_t i = 0; i < muon_pt.size(); ++i) {
//         if (muon_mask[i]) {
//             leptonInfo.push_back({muon_pt[i], muon_phi[i], muon_eta[i], muon_pdgId[i]});
//         }
//     }

//     // Añadir información de electrones seleccionados
//     for (size_t i = 0; i < electron_pt.size(); ++i) {
//         if (electron_mask[i]) {
//             leptonInfo.push_back({electron_pt[i], electron_phi[i], electron_eta[i], electron_pdgId[i]});
//         }
//     }

//     return leptonInfo;
// }

#endif