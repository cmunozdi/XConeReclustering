#!/usr/bin/env python3

import ROOT
from tqdm import tqdm
import os
import numpy as np
import pandas as pd

def get_root_files(directory):
    root_files = []
    for dirpath, _, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith('.root'):
                root_files.append(os.path.join(dirpath, filename))
    return root_files

def runSelectionOn(infile='/eos/user/c/cmunozdi/PFNano_Run3/mc_summer23/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v2_BTV_Run3_2023_Comm_MINIAODv4/250317_151402/0000/MC_defaultAK4_2023_143.root',
                   outfile='XCone_reclusteringTTtoLNu2Q_2023test1.root'):
    
    """
    runs a simple selection on a NanoGEN file and saves a snapshot of the event
    the selection makes use of RDataFrames and additional functions coded in selection_helpers.h
    """

    # ROOT.gSystem.Load("$LCIO/lib/libRivet.so")
    # ROOT.gSystem.Load("$LCIO/lib/libHepMC.so")
    # ROOT.gSystem.Load("$LCIO/lib/libfastjet.so")
    
#     ROOT.gInterpreter.Declare('#include "Particle.h"')

    ROOT.gInterpreter.Declare('#include "selection_helpers_BoostedTopQuark.h"')

    ROOT.ROOT.EnableImplicitMT()
    rdf = ROOT.RDataFrame('Events', infile) 
    # rdf = rdf.Range(300)
    #lepton selection: select in the kinematics region of interest 
    rdf = rdf.Define('muon', 'triggerLepton(Muon_pt, Muon_eta, Muon_phi, Muon_pdgId, Jet_eta, Jet_phi, true)') \
             .Define('electron', 'triggerLepton(Electron_pt, Electron_eta, Electron_phi, Electron_pdgId, Jet_eta, Jet_phi, false)') \
             .Define('lepton', 'CombineLeptons(muon, electron)') \
             .Define('n_leptons', 'lepton.n_lep') \
             .Filter('n_leptons==1')
    # print("lepton selection")

    #b-jets selection:
    #Tight working points:
         #for 2022_Summer22 is 0.7183
         #for 2022_Summer22EE is 0.73
         #for 2023_Summer23 is 0.6553
         #for 2023_Summer23BPrix is 0.7994

    #I think is better to check the btagDeepFlavB value after de reclustering (different values for different campaigns) && Jet_btagDeepFlavB>0.6553

    rdf = rdf.Define('jet', 'Jet_pt>30 && abs(Jet_eta) < 2.4') \
             .Define('n_jets', 'Sum(jet)') \
             .Filter('n_jets>0')
#     for c in ['partonFlavour', 'pt', 'eta', 'phi']:
#         rdf = rdf.Define(f'PFCandsbJet_{c}',f'Jet_{c}[pfcandsbjet]')
    # print("jet selection")

    #Supresion of multijet backgrounds from the production of light-flavour quarks and gluons
    rdf = rdf.Define('pt_miss', 'Get_pTmiss(PFCands_pt, PFCands_phi)') \
             .Filter('pt_miss>50')

    # Build jets from PFCands#, 2, 1.2, 2.0, 400., 2.4, 3, 0.4, 2., 25., 2.4, true, false)') \
    rdf = rdf.Define('jet_XConefromPFCands','buildXConeJets(PFCands_pt, PFCands_eta, PFCands_phi, PFCands_mass, PFCands_pdgId)') \
             .Filter('jet_XConefromPFCands.fatjets.n_jets > 0') \
# #              .Define('n_jets_XConefromPFCands', 'jet_XConefromPFCands.topjets.n_jets') \
# #              .Define('jet_AntikTfromPFCands_pt', 'jet_AntikTfromPFCands.pt') \
# #              .Define('jet_AntikTfromPFCands_eta', 'jet_AntikTfromPFCands.eta') \
# #              .Define('jet_AntikTfromPFCands_phi', 'jet_AntikTfromPFCands.phi')

    # rdf = rdf.Define('jet_XConefromGenCands','buildXConeJets(GenCands_pt, GenCands_eta, GenCands_phi, GenCands_mass, GenCands_pdgId)') \
            #  .Filter('jet_XConefromGenCands.fatjets.n_jets > 0') \

#     # Define properties of jetºs for each index
# #     for i in range(1, 3):
# #         for c in ['pt', 'eta', 'phi']:
# #             rdf = rdf.Define(
# #                 f'jet_{i}_AntikTfromPFCands_{c}',
# #                 f'n_jets_AntikTfromPFCands > {i-1} ? jet_AntikTfromPFCands[{i-1}].{c} : -100'
# #             )

#     print("Jet antikt pfcands")
    
#    # rdf = rdf.Define( 'jet_XConefromPFCands', 'buildXConeJets(PFCands_pt, PFCands_eta, PFCands_phi,PFCands_pdgId, 2, 1.2, 2., false)') \
#     #         .Define('n_jets_XConefromPFCands', 'jet_XConefromPFCands.size()')



    
    #save the selection
    columns=['event', 'lepton', 'n_jets', 'Jet_btagDeepFlavB', 'pt_miss'] #,'n_fatjets'] 'muon', 'electron',#'n_jets_AntikTfromPFCands', , 'jet_1_AntikTfromPFCands_pt', 'jet_2_AntikTfromPFCands_pt', 'jet_1_AntikTfromPFCands_eta', 'jet_2_AntikTfromPFCands_eta', 'jet_1_AntikTfromPFCands_phi', 'jet_2_AntikTfromPFCands_phi']
    columns+=['jet_XConefromPFCands'] #, 'jet_XConefromGenCands'
    columns+=['PFCands_pt', 'PFCands_phi', 'PFCands_eta', 'PFCands_pdgId']
    # columns+=['GenCands_pt', 'GenCands_phi', 'GenCands_eta', 'GenCands_pdgId']
    # columns+=['Muon_pt','Muon_phi','Muon_eta', 'Muon_pdgId', 'Electron_pt','Electron_phi','Electron_eta', 'Electron_pdgId', 'lepton_pt','lepton_phi','lepton_eta', 'lepton_pdgId']
#     for i in range(1,3):
#     columns+=['jet_AntikTfromPFCands_pt','jet_AntikTfromPFCands_eta','jet_AntikTfromPFCands_phi']
    # columns+=['PFCands_pt', 'PFCands_phi', 'PFCands_eta', 'PFCands_pdgId']
    # print("creating required columns")
    # print("Available columns in DataFrame:", rdf.GetColumnNames())
    
#     rdf.Foreach(lambda ev: print(f"Event: {ev}"), ["event"])
    

    # Verificar que columns es una lista de cadenas de texto
    if not isinstance(columns, list):
        raise TypeError("columns debe ser una lista de cadenas de texto")
    for col in columns:
        if not isinstance(col, str):
            raise TypeError(f"Cada elemento en columns debe ser una cadena de texto. Encontrado: {type(col)}")

    # print("creating required columns")
    # print("Available columns in DataFrame:", rdf.GetColumnNames())
    # for col in tqdm(columns, desc="Checking columns"):
    #     pass 


    rdf.Snapshot('Events', outfile, columns)
    # print("snapshot")
    ROOT.ROOT.DisableImplicitMT()

# directory = '/eos/user/c/cmunozdi/PFNano_Run3/mc_summer23/TTtoLNu2Q'
directory = '/eos/user/c/cmunozdi/PFNano_Run3/data_2023/Muon1'

input_files = get_root_files(directory)

runSelectionOn(input_files, '/eos/user/c/cmunozdi/AnalysisMCSamples/2023/data/XCone_reclustering_data_Muon1_2023.root')
# runSelectionOn(infile="../btvnano-prod/MC_defaultAK4_2023_TTtoLNu2Q_MT-166p5.root", outfile='XCone_reclusteringTTtoLNu2Q_MT-166p5_2023test1.root') #    infile= input_files infile='nano_mc_Run3_TTtoLNu2Q_MT_171p5_NANO_Skim_haddnano.root', outfile='XCone_reclusteringTTtoLNu2Q.root'

