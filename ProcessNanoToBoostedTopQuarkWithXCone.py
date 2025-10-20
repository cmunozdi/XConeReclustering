import argparse
import os
import ROOT
import subprocess


# Verify if FASTJET_CONTRIB_BASE is set
fjc_base = os.environ.get("FASTJET_CONTRIB_BASE")
if not fjc_base:
    raise RuntimeError("FASTJET_CONTRIB_BASE is not set. Ensure the environment is configured correctly.")

# Agregar la ruta de los headers al intérprete de ROOT
ROOT.gSystem.AddIncludePath(f'-I{fjc_base}/include')




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

output_file_events = args.output.replace(".root", "_events.root")
isMC = args.isMC

# Crear un nombre de archivo para rdf_runs basado en args.output
output_file_runs = args.output.replace(".root", "_runs.root")
output_file_lumi = args.output.replace(".root", "_lumi.root")

# ROOT.gSystem.Load("$LCIO/lib/libRivet.so")
# ROOT.gSystem.Load("$LCIO/lib/libHepMC.so")
# ROOT.gSystem.Load("$LCIO/lib/libfastjet.so")

ROOT.gInterpreter.Declare('#include "selection_helpers_BoostedTopQuark.h"') ##RUN CRAB
# Llamar a la función para inicializar `cset`
ROOT.gInterpreter.ProcessLine("initializeCorrectionSet();")
# ROOT.gInterpreter.ProcessLine("initializeBTagCorrectionSet();")
# verbosity = ROOT.Experimental.RLogScopedVerbosity(ROOT.Detail.RDF.RDFLogChannel(), ROOT.Experimental.ELogLevel.kInfo)

ROOT.ROOT.EnableImplicitMT()
# Creates the RDataFrame
rdf_runs = ROOT.RDataFrame("Runs", input_files)
rdf_runs.Snapshot('Runs', output_file_runs)
rdf_lumi = ROOT.RDataFrame("LuminosityBlocks", input_files)
rdf_lumi.Snapshot('LuminosityBlocks', output_file_lumi)
rdf = ROOT.RDataFrame("Events", input_files)


# Define Jet_passTightID
rdf = rdf.Define("Jet_passJetIdTight", "pass_JetIdTightLepVeto(Jet_eta, Jet_neHEF, Jet_neEmEF, Jet_chMultiplicity, Jet_neMultiplicity, \
                 Jet_chHEF, Jet_muEF, Jet_chEmEF, false)") \
         .Define("Jet_passJetIdTightLepVeto", "pass_JetIdTightLepVeto(Jet_eta, Jet_neHEF, Jet_neEmEF, Jet_chMultiplicity, Jet_neMultiplicity, \
                 Jet_chHEF, Jet_muEF, Jet_chEmEF, true)") \
        #  .Redefine("Jet_passJetIdTightLepVeto", "std::vector<bool>(Jet_passJetIdTightLepVeto.begin(), Jet_passJetIdTightLepVeto.end())") \

# Lepton selection: select in the kinematic region of interest 
        #  .Redefine('Electron_tightId', 'std::vector<bool>(Electron_tightId.begin(), Electron_tightId.end())') \

rdf = rdf.Define('muon', 'triggerLepton(Muon_pt, Muon_eta, Muon_phi, Muon_pdgId, Muon_jetIdx, Muon_tightId, PFCands_pt, PFCands_eta, PFCands_phi, \
                 PFCands_pdgId, Jet_pt, Jet_eta, Jet_phi, Jet_mass, Jet_passJetIdTight, Jet_rawFactor, Jet_area, Rho_fixedGridRhoFastjetAll, 25, true, false)') \
         .Define('Electron_tightId', 'Electron_cutBased >= 4') \
         .Define('electron', 'triggerLepton(Electron_pt, Electron_eta, Electron_phi, Electron_pdgId, Electron_jetIdx, Electron_tightId, \
                 PFCands_pt, PFCands_eta, PFCands_phi, PFCands_pdgId, Jet_pt, Jet_eta, Jet_phi, Jet_mass, Jet_passJetIdTight, Jet_rawFactor, Jet_area, Rho_fixedGridRhoFastjetAll, 25, false, false)') \
         .Define('lepton', 'CombineLeptons(muon, electron)') \
         .Define('n_leptons', 'lepton.n_lep') \
         .Define('Num_lep_above_15', "Sum(Muon_pt > 15 && abs(Muon_eta) < 2.4 && Muon_tightId) + Sum(Electron_pt > 15 && abs(Electron_eta) < 2.4 && Electron_tightId)") \
        #  .Filter('Num_lep_above_15 == 1')
        #  .Filter('n_leptons >=1')

# b-jets selection:
    #Tight working points:
         #for 2022_Summer22 is 0.7183
         #for 2022_Summer22EE is 0.73
         #for 2023_Summer23 is 0.6553
         #for 2023_Summer23BPrix is 0.7994
    #I think is better to check the btagDeepFlavB value after de reclustering (different values for different campaigns) && Jet_btagDeepFlavB>0.6553
    # Lowering down the cuts (original && (Jet_pt>22.5) && (abs(Jet_eta) <3) Jet_passJetIdTightLepVeto
rdf = rdf.Define('pseudo_bjet', 'Jet_pt>22.5 && abs(Jet_eta) <3 && Jet_passJetIdTightLepVeto') \
         .Define('n_pseudo_bjets', 'Sum(pseudo_bjet)') \
        #  .Filter('n_jets > 0')


# Suppression of multijet backgrounds from the production of light-flavor quarks and gluons Get_pTmiss(PFCands_pt, PFCands_phi)
# rdf = rdf.Redefine('pt_miss', 'PuppiMET_pt') \
        #  .Filter('pt_miss>50')


# XCone jet clustering for PFCands
rdf = rdf.Define('jet_XConefromPFCands','buildXConeJets(PFCands_pt, PFCands_eta, PFCands_phi, PFCands_mass, PFCands_pdgId)') \
         .Define('n_fatjets', 'jet_XConefromPFCands.fatjets.n_jets') \
        #  .Filter('jet_XConefromPFCands.fatjets.n_jets > 0') \

# # XCone for PFCands: pass detector selection flag: #Reducing cuts as well   && (n_jets > 0) && (pt_miss > 37.5) && (n_fatjets > 0) || lepton_20.n_lep >= 1 || lepton_25.n_lep >= 1 || lepton_30.n_lep >= 1  && (n_pseudo_bjets > 0) && (PuppiMET_pt > 37.5)
rdf = rdf.Define('pass_detector_selection', '''
                return  ((lepton.n_lep >= 1) && (Num_lep_above_15 >= 1) && (n_pseudo_bjets > 0) && (PuppiMET_pt > 37.5) && (n_fatjets > 0) && (jet_XConefromPFCands.topjets.n_subjets ==3) && (jet_XConefromPFCands.topjets.pt > 200));
                ''')

if isMC:

    # Only for ttbar samples, add the top pt reweighting
    rdf = rdf.Define("GenPartTop_pt", "GenPart_pt[(abs(GenPart_pdgId) == 6) && (GenPart_statusFlags & (1<<13)) != 0]") \
             .Define("TopPtWeight_dataPowheg", "GenPartTop_pt.size() == 2 ? sqrt( exp(0.0615 - 0.0005*GenPartTop_pt[0]) * exp(0.0615 - 0.0005*GenPartTop_pt[1]) ) : 1.0") \
             .Define("TopPtWeight_NNLOpNLOEW", "GenPartTop_pt.size() == 2 ? sqrt((0.103*exp(-0.0118*GenPartTop_pt[0])-0.000134*GenPartTop_pt[0]+0.973)*(0.103*exp(-0.0118*GenPartTop_pt[1])-0.000134*GenPartTop_pt[1]+0.973)) : 1.0") \
    #Define b-tagging weights
    # rdf = rdf.Define("btagWeight", "compute_btagWeight(Jet_pt, Jet_eta, Jet_hadronFlavour, Jet_btagUParTAK4B)") \
    # Only execute this line if we are processing MC  Get_pTmiss(GenCands_pt, GenCands_phi) && (n_gen_fatjets > 0
    rdf = rdf.Define('GenCands_jetIdx', 'calculateGenCandsJetIdx(GenJetCands_genCandsIdx, GenJetCands_jetIdx, nGenCands, nGenJet)') \
             .Define('gen_muon', 'triggerLepton(GenCands_pt, GenCands_eta, GenCands_phi, GenCands_pdgId, GenCands_jetIdx, ROOT::VecOps::RVec<bool>(), GenCands_pt, GenCands_eta, GenCands_phi, GenCands_pdgId, GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass, ROOT::VecOps::RVec<bool>(GenJet_eta.size(), true), ROOT::VecOps::RVec<float>(GenJet_pt.size(), 0.0), ROOT::VecOps::RVec<float>(GenJet_pt.size(), 0.0), 0, 25, true, false, true)') \
             .Define('gen_electron', 'triggerLepton(GenCands_pt, GenCands_eta, GenCands_phi, GenCands_pdgId, GenCands_jetIdx, ROOT::VecOps::RVec<bool>(), GenCands_pt, GenCands_eta, GenCands_phi, GenCands_pdgId, GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass, ROOT::VecOps::RVec<bool>(GenJet_eta.size(), true), ROOT::VecOps::RVec<float>(GenJet_pt.size(), 0.0), ROOT::VecOps::RVec<float>(GenJet_pt.size(), 0.0), 0, 25, false, false, true)') \
             .Define('gen_lepton', 'CombineLeptons(gen_muon, gen_electron)') \
             .Define('n_gen_leptons', 'gen_lepton.n_lep') \
             .Define('num_gen_lep_above_15', "Sum(GenCands_pt > 15 && abs(GenCands_eta) < 2.4 && (abs(GenCands_pdgId)==13 || abs(GenCands_pdgId)==11))") \
             .Define('gen_pseudo_bjet', 'GenJet_pt>22.5 && abs(GenJet_eta) <3') \
             .Define('n_gen_pseudo_bjets', 'Sum(gen_pseudo_bjet)') \
            #  .Redefine('gen_pt_miss', 'GenMET_pt') \
    rdf = rdf.Define('jet_XConefromGenCands', 'buildXConeJets(GenCands_pt, GenCands_eta, GenCands_phi, GenCands_mass, GenCands_pdgId)') \
             .Define('n_gen_fatjets', 'jet_XConefromGenCands.fatjets.n_jets') \
             .Define('pass_particle_selection', '''
                    return ((n_gen_leptons >= 1) && (num_gen_lep_above_15 >= 1) && (n_gen_pseudo_bjets > 0) && (GenMET_pt > 37.5) && (n_gen_fatjets > 0) && (jet_XConefromGenCands.topjets.n_subjets ==3) && (jet_XConefromGenCands.topjets.pt > 200));
                    ''') \
             .Redefine("GenCands_jetIdx", "std::vector<Short_t>(GenCands_jetIdx.begin(), GenCands_jetIdx.end())") \
            #  .Filter('jet_XConefromGenCands.fatjets.n_jets > 0')
else:
    # If not MC, we don't need to define gen_* variables
    rdf = rdf.Define('pass_particle_selection', 'false')

# Filter the events with at least one flag set to true
rdf = rdf.Filter('pass_detector_selection || pass_particle_selection') 

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


# opts = ROOT.RDF.RSnapshotOptions()
# opts.fMode = "RECREATE"
# opts.fOverwriteIfExists = True

# Create a snapshot with the selected columns
rdf.Snapshot('Events', output_file_events)#, "", opts)#, columns)
# output_file.WriteObject(rdf, "Events")
# output_file.WriteObject(rdf_runs, "Runs")
# Copiar el árbol "Runs" desde el primer input file
# f_in = ROOT.TFile.Open(input_files[0])
# t_runs = f_in.Get("Runs")
# f_out = ROOT.TFile.Open(output_file, "UPDATE")
# if t_runs:
#     f_out.cd()
#     t_runs.CloneTree(-1, "fast").Write("Runs")
# f_out.Close()
# f_in.Close()

# rdf_runs.Snapshot('Runs', output_file, "", opts)

ROOT.ROOT.DisableImplicitMT()

print(f"Processing completed. Output file: {args.output}")

# Merge the output files if there are multiple input files
output_files = [output_file_events, output_file_runs, output_file_lumi]

combined_output_file = args.output

hadd_command = ["hadd", "-f", combined_output_file] + output_files

try: 
    subprocess.run(hadd_command, check=True)
    print(f"Successfully merged output files into: {combined_output_file}")
    subprocess.run(["rm"] + output_files, check=True)
except subprocess.CalledProcessError as e:
    print(f"Error merging output files: {e}")
    

