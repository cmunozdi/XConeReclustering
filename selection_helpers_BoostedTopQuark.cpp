// selection_helpers_BoostedTopQuark.cc
#include "selection_helpers_BoostedTopQuark.h"

// Implementación de funciones o clases declaradas en el archivo de cabecera
bool debug = false;
std::unique_ptr<correction::CorrectionSet> cset;
std::unique_ptr<correction::CorrectionSet> sf_btagset;
std::unique_ptr<correction::CorrectionSet> eff_btagset;

void exampleFunction() {
    if (debug) {
        std::cout << "Debugging enabled in selection_helpers_BoostedTopQuark.cc" << std::endl;
    }
}