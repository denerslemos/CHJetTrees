#include "CallLibraries.h"  // call libraries from ROOT and C++

// Reconstructed Jets
std::unique_ptr<TTreeReaderArray<int>> JetRecoType;
std::unique_ptr<TTreeReaderArray<float>> JetRecoE;
std::unique_ptr<TTreeReaderArray<float>> JetRecoPx;
std::unique_ptr<TTreeReaderArray<float>> JetRecoPy;
std::unique_ptr<TTreeReaderArray<float>> JetRecoPz;
std::unique_ptr<TTreeReaderArray<float>> JetRecoM;
std::unique_ptr<TTreeReaderArray<int>> JetRecoPDG;
std::unique_ptr<TTreeReaderArray<unsigned int>> JetRecoCBegin;
std::unique_ptr<TTreeReaderArray<unsigned int>> JetRecoCEnd;
std::unique_ptr<TTreeReaderArray<int>> JetRecoCIdx;

// Reconstructed charged particles (Tracks)
std::unique_ptr<TTreeReaderArray<float>> TrkRecoE;
std::unique_ptr<TTreeReaderArray<float>> TrkRecoPx;
std::unique_ptr<TTreeReaderArray<float>> TrkRecoPy;
std::unique_ptr<TTreeReaderArray<float>> TrkRecoPz;
std::unique_ptr<TTreeReaderArray<float>> TrkRecoM;
std::unique_ptr<TTreeReaderArray<int>> TrkRecoPDG;

// Reco to Gen track association
std::unique_ptr<TTreeReaderArray<unsigned int>> TrkPartAssocRec;
std::unique_ptr<TTreeReaderArray<unsigned int>> TrkPartAssocSim;
std::unique_ptr<TTreeReaderArray<float>> TrkPartAssocWeight;

// Generated Jets
std::unique_ptr<TTreeReaderArray<int>> JetGenType;
std::unique_ptr<TTreeReaderArray<float>> JetGenE;
std::unique_ptr<TTreeReaderArray<float>> JetGenPx;
std::unique_ptr<TTreeReaderArray<float>> JetGenPy;
std::unique_ptr<TTreeReaderArray<float>> JetGenPz;
std::unique_ptr<TTreeReaderArray<float>> JetGenM;
std::unique_ptr<TTreeReaderArray<int>> JetGenPDG;
std::unique_ptr<TTreeReaderArray<unsigned int>> JetGenCBegin;
std::unique_ptr<TTreeReaderArray<unsigned int>> JetGenCEnd;
std::unique_ptr<TTreeReaderArray<int>> JetGenCIdx;

// Generated particles -> only stable particles from generator level -> Not available ?
std::unique_ptr<TTreeReaderArray<float>> TrkGenPx;
std::unique_ptr<TTreeReaderArray<float>> TrkGenPy;
std::unique_ptr<TTreeReaderArray<float>> TrkGenPz;
std::unique_ptr<TTreeReaderArray<float>> TrkGenE;
std::unique_ptr<TTreeReaderArray<float>> TrkGenM;
std::unique_ptr<TTreeReaderArray<int>> TrkGenPDG;

// MCParticles particles -> full generator information + sec decay from GEANT
std::unique_ptr<TTreeReaderArray<double>> TrkMCGenPx;
std::unique_ptr<TTreeReaderArray<double>> TrkMCGenPy;
std::unique_ptr<TTreeReaderArray<double>> TrkMCGenPz;
std::unique_ptr<TTreeReaderArray<double>> TrkMCGenM;
std::unique_ptr<TTreeReaderArray<int>> TrkMCGenPDG;
std::unique_ptr<TTreeReaderArray<int>> TrkMCGenStatus;



// Reads jet tree variables
/*
Arguments:
chain: TChain of input files
tree_reader: the object for TTreeReader
*/
void JetTreeReader(TChain* chain, std::unique_ptr<TTreeReader>& tree_reader) {

    tree_reader = std::make_unique<TTreeReader>(chain);

	// Reconstructed Jets
    JetRecoType = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "ReconstructedChargedJets.type");
    JetRecoE = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "ReconstructedChargedJets.energy");
    JetRecoPx = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "ReconstructedChargedJets.momentum.x");
    JetRecoPy = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "ReconstructedChargedJets.momentum.y");
    JetRecoPz = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "ReconstructedChargedJets.momentum.z");
    JetRecoM = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "ReconstructedChargedJets.mass");
    JetRecoPDG = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "ReconstructedChargedJets.PDG");
    JetRecoCBegin = std::make_unique<TTreeReaderArray<unsigned int>>(*tree_reader, "ReconstructedChargedJets.particles_begin");
    JetRecoCEnd = std::make_unique<TTreeReaderArray<unsigned int>>(*tree_reader, "ReconstructedChargedJets.particles_end");
    JetRecoCIdx = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "_ReconstructedChargedJets_particles.index");

	// Reconstructed charged particles (Tracks)
    TrkRecoE = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "ReconstructedChargedParticles.energy");
    TrkRecoPx = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "ReconstructedChargedParticles.momentum.x");
    TrkRecoPy = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "ReconstructedChargedParticles.momentum.y");
    TrkRecoPz = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "ReconstructedChargedParticles.momentum.z");
    TrkRecoM = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "ReconstructedChargedParticles.mass");
    TrkRecoPDG = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "ReconstructedChargedParticles.PDG");
	TrkPartAssocRec = std::make_unique<TTreeReaderArray<unsigned int>>(*tree_reader, "ReconstructedChargedParticleAssociations.recID"); // Reco <-> MCParticle
	TrkPartAssocSim = std::make_unique<TTreeReaderArray<unsigned int>>(*tree_reader, "ReconstructedChargedParticleAssociations.simID");
	TrkPartAssocWeight = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "ReconstructedChargedParticleAssociations.weight");

	// Generated Jets
	JetGenType = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "GeneratedChargedJets.type");
    JetGenE = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "GeneratedChargedJets.energy");
    JetGenPx = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "GeneratedChargedJets.momentum.x");
    JetGenPy = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "GeneratedChargedJets.momentum.y");
    JetGenPz = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "GeneratedChargedJets.momentum.z");
    JetGenM = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "GeneratedChargedJets.mass");
    JetGenPDG = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "GeneratedChargedJets.PDG");
    JetGenCBegin = std::make_unique<TTreeReaderArray<unsigned int>>(*tree_reader, "GeneratedChargedJets.particles_begin");
    JetGenCEnd = std::make_unique<TTreeReaderArray<unsigned int>>(*tree_reader, "GeneratedChargedJets.particles_end");
    JetGenCIdx = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "_GeneratedChargedJets_particles.index");

	// MC Gen -> add GEANT Stuff
	TrkMCGenStatus = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "MCParticles.generatorStatus");
	TrkMCGenPx = std::make_unique<TTreeReaderArray<double>>(*tree_reader, "MCParticles.momentum.x");
	TrkMCGenPy = std::make_unique<TTreeReaderArray<double>>(*tree_reader, "MCParticles.momentum.y");
	TrkMCGenPz = std::make_unique<TTreeReaderArray<double>>(*tree_reader, "MCParticles.momentum.z");
	TrkMCGenM = std::make_unique<TTreeReaderArray<double>>(*tree_reader, "MCParticles.mass");
	TrkMCGenPDG = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "MCParticles.PDG");
	TrkMCGenCharge = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "MCParticles.charge");

	// Full gen level
	/*
	TrkGenE = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "GeneratedParticles.energy");
	TrkGenPx = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "GeneratedParticles.momentum.x");
	TrkGenPy = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "GeneratedParticles.momentum.y");
	TrkGenPz = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "GeneratedParticles.momentum.z");
	TrkGenM = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "GeneratedParticles.mass");
	TrkGenPDG = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "GeneratedParticles.PDG");
	TrkGenCharge = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "GeneratedParticles.charge");
	*/

}
