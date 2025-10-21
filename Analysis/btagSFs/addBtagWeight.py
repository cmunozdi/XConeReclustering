import argparse
import os
import ROOT
import subprocess

run_small_test = False  # Set to True for local testing, False for production

# # Verify if FASTJET_CONTRIB_BASE is set
# fjc_base = os.environ.get("FASTJET_CONTRIB_BASE")
# if not fjc_base:
#     raise RuntimeError("FASTJET_CONTRIB_BASE is not set. Ensure the environment is configured correctly.")

# # Agregar la ruta de los headers al intérprete de ROOT
# ROOT.gSystem.AddIncludePath(f'-I{fjc_base}/include')




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

ROOT.gInterpreter.Declare('#include "../selection_helpers_BoostedTopQuark.h"') ##RUN CRAB
# Llamar a la función para inicializar `cset`
ROOT.gInterpreter.ProcessLine("initializeBTagCorrectionSet();")
# verbosity = ROOT.Experimental.RLogScopedVerbosity(ROOT.Detail.RDF.RDFLogChannel(), ROOT.Experimental.ELogLevel.kInfo)

if not run_small_test: ROOT.ROOT.EnableImplicitMT()
# Creates the RDataFrame
rdf_runs = ROOT.RDataFrame("Runs", input_files)
rdf_runs.Snapshot('Runs', output_file_runs)
rdf_lumi = ROOT.RDataFrame("LuminosityBlocks", input_files)
rdf_lumi.Snapshot('LuminosityBlocks', output_file_lumi)
rdf = ROOT.RDataFrame("Events", input_files)
if run_small_test: rdf = rdf.Range(10000)

rdf = rdf.Redefine("btagWeight", "compute_btagWeight(Jet_pt, Jet_eta, Jet_hadronFlavour, Jet_btagDeepFlavB)") \


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

if not run_small_test: ROOT.ROOT.DisableImplicitMT()

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
    

