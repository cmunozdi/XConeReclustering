// selection_helpers_BoostedTopQuark.cc
#include "selection_helpers_BoostedTopQuark.h"

// Implementación de funciones o clases declaradas en el archivo de cabecera
bool debug = false;
std::unique_ptr<correction::CorrectionSet> cset;
std::unique_ptr<correction::CorrectionSet> sf_btagset;
std::unique_ptr<correction::CorrectionSet> eff_btagset;
std::unique_ptr<correction::CorrectionSet> sf_muoset;
std::unique_ptr<correction::CorrectionSet> sf_purew;

std::string g_jec_tag  = "Summer23Prompt23_V2_MC";
std::string g_lvl_tag  = "L1L2L3Res";
std::string g_algo_tag = "AK4PFPuppi";
// std::string jec_key; // = g_jec_tag + '_' + g_lvl_tag + '_' + g_algo_tag;
bool g_isMC = true;

void exampleFunction() {
    if (debug) {
        std::cout << "Debugging enabled in selection_helpers_BoostedTopQuark.cc" << std::endl;
    }
}