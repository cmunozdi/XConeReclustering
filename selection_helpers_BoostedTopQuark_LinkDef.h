#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

// Declara las estructuras y funciones que deseas incluir en el diccionario
#pragma link C++ struct Lepton+;
#pragma link C++ struct XConeReclusteredJets+;
#pragma link C++ struct JetReclus+;
#pragma link C++ struct TopJetReclus+;
#pragma link C++ function pickLepton+;
#pragma link C++ function addLeptonJetInfo+;
#pragma link C++ function singleLevelEnergyCorr+;
#pragma link C++ function compute_ptrel_new+;
#pragma link C++ function calculateGenCandsJetIdx+;
#pragma link C++ function pass_JetIdTightLepVeto+;
#pragma link C++ function buildXConeJets+;
#pragma link C++ function initPlugin+;


#endif