import argparse
import os
import ROOT

# Input argument parsing
parser = argparse.ArgumentParser(description="Process PFNANO with XCone algo for MC or Data.")
parser.add_argument("--input", required=True, help="Input file or directory")
parser.add_argument("--output", required=True, help="Output file")
parser.add_argument("--isMC", action="store_true", help="Are we processing MC? (buy default is Data)")

args = parser.parse_args()

# Verify if the input is a file or a directory
if os.path.isdir(args.input):
    input_files = [os.path.join(args.input, f) for f in os.listdir(args.input) if f.endswith(".root")]
elif os.path.isfile(args.input):
    input_files = [args.input]
else:
    raise ValueError(f"Input is not a valid file or direcroty: {args.input}")

output_file = args.output
isMC = args.isMC


# ROOT.gSystem.Load("$LCIO/lib/libRivet.so")
# ROOT.gSystem.Load("$LCIO/lib/libHepMC.so")
# ROOT.gSystem.Load("$LCIO/lib/libfastjet.so")

# ROOT.gInterpreter.Declare('#include "Particle.h"')
# ROOT.gInterpreter.Declare('#include "/afs/cern.ch/user/c/cmunozdi/Analysis/GeneratePFNanos/CMSSW_13_0_13/src/XConeReclustering/selection_helpers_BoostedTopQuark.h"') ##RUN LOCAL
ROOT.gInterpreter.Declare('#include "selection_helpers_BoostedTopQuark.h"') ##RUN CRAB
verbosity = ROOT.Experimental.RLogScopedVerbosity(ROOT.Detail.RDF.RDFLogChannel(), ROOT.Experimental.ELogLevel.kInfo)

ROOT.ROOT.EnableImplicitMT()
# Creates the RDataFrame
rdf = ROOT.RDataFrame("Events", input_files)

# Lepton selection: select in the kinematic region of interest
rdf = rdf.Define('muon', 'triggerLepton(Muon_pt, Muon_eta, Muon_phi, Muon_pdgId, Jet_eta, Jet_phi, true)') \
         .Define('electron', 'triggerLepton(Electron_pt, Electron_eta, Electron_phi, Electron_pdgId, Jet_eta, Jet_phi, false)') \
         .Define('lepton', 'CombineLeptons(muon, electron)') \
         .Define('n_leptons', 'lepton.n_lep') \
         .Filter('n_leptons ==1')


# b-jets selection:
    #Tight working points:
         #for 2022_Summer22 is 0.7183
         #for 2022_Summer22EE is 0.73
         #for 2023_Summer23 is 0.6553
         #for 2023_Summer23BPrix is 0.7994
    #I think is better to check the btagDeepFlavB value after de reclustering (different values for different campaigns) && Jet_btagDeepFlavB>0.6553
rdf = rdf.Define('jet', 'Jet_pt>30 && abs(Jet_eta) <2.4') \
         .Define('n_jets', 'Sum(jet)') \
         .Filter('n_jets > 0')


# Supression of multijet backgrounds from the production of light-flavor quarks and gluons
rdf = rdf.Define('pt_miss', 'Get_pTmiss(PFCands_pt, PFCands_phi)') \
         .Filter('pt_miss>50')


# XCone jet clustering for PFCands
rdf = rdf.Define('jet_XConefromPFCands','buildXConeJets(PFCands_pt, PFCands_eta, PFCands_phi, PFCands_mass, PFCands_pdgId)') \
         .Filter('jet_XConefromPFCands.fatjets.n_jets > 0') \

if isMC:
    # Only execute this line if we are processing MC
    rdf = rdf.Define('jet_XConefromGenCands', 'buildXConeJets(GenCands_pt, GenCands_eta, GenCands_phi, GenCands_mass, GenCands_pdgId)') \
            #  .Filter('jet_XConefromGenCands.fatjets.n_jets > 0')

# Save the columns needed (if running this on crab, maybe better to save all the columns, since selecting the interesting events will notably reduce the size of the file)
# columns = ['event', 'lepton', 'n_jets', 'Jet_btagDeepFlavB', 'pt_miss']
# columns += ['jet_XConefromPFCands']
# columns += ['PFCands_pt', 'PFCands_phi', 'PFCands_eta', 'PFCands_pdgId']

# if isMC:
#     # Only add these columns if we are processing MC
#     columns += ['jet_XConefromGenCands']
#     columns += ['GenCands_pt', 'GenCands_phi', 'GenCands_eta', 'GenCands_pdgId']


# # Verify that columns is a list of strings
# if not isinstance(columns, list):
#     raise TypeError("columns must be a list of strings")
# for col in columns:
#     if not isinstance(col, str):
#         raise TypeError(f"Each element in columns must be a string. Found: {type(col)}")

# Create a snapshot with the selected columns
rdf.Snapshot('Events', output_file)#, columns)

ROOT.ROOT.DisableImplicitMT()

print(f"Processing completed. Output file: {output_file}")