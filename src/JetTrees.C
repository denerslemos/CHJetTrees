#include "JetTreeReader.h"  // call libraries from ROOT and C++

void JetTrees(TString InputFileList, TString OutputFile){

	typedef ROOT::Math::PxPyPzEVector LorentzVector;

	int NEVENTS = 0;
	// Reco Jets (Variable-length vectors for multiple jets per event)
	std::vector<float> RecoJet_pt;
	std::vector<float> RecoJet_eta;
	std::vector<float> RecoJet_phi;
	std::vector<float> RecoJet_E;
	std::vector<float> RecoJet_M;
	std::vector<bool> RecoJet_hasElectron;
	std::vector<float> RecoJet_maxPtPart_pt; 
	// Gen Jets (Variable-length vectors for multiple jets per event)
	std::vector<float> GenJet_pt;
	std::vector<float> GenJet_eta;
	std::vector<float> GenJet_phi;
	std::vector<float> GenJet_E;
	std::vector<float> GenJet_M;
	std::vector<bool> GenJet_hasElectron;
	std::vector<float> GenJet_maxPtPart_pt;

	// Read the list of input file(s)
	fstream FileList;
	FileList.open(Form("%s",InputFileList.Data()), ios::in);
	if(!FileList.is_open()){cout << "List of input files not founded!" << endl; return;}{cout << "List of input files founded! --> " << InputFileList.Data() << endl;}

	// Make a chain and a vector of file names
	std::vector<TString> FileListVector;
	string FileChain;
	while(getline(FileList, FileChain)){FileListVector.push_back(FileChain.c_str());}
	FileList.close();	
	TChain *mychain = new TChain("events");
	for (std::vector<TString>::iterator listIterator = FileListVector.begin(); listIterator != FileListVector.end(); listIterator++){
		TFile *testfile = TFile::Open(*listIterator,"READ");
		if(testfile && !testfile->IsZombie() && !testfile->TestBit(TFile::kRecovered)){ // safety against corrupted files
			cout << "Adding file " << *listIterator << " to the chains" << endl; // adding files to the chains for each step
			mychain->Add(*listIterator);
		}else{cout << "File: " << *listIterator << " failed!" << endl;}
	}

	// Reading trees
	std::unique_ptr<TTreeReader> tree_reader;
	JetTreeReader(mychain, tree_reader); 
	JetTree->Branch("NEVENTS", &NEVENTS, "NEVENTS/I"); // To keep track of the event number	
	
	// Make Output
	TFile *OutFile = TFile::Open(Form("%s",OutputFile.Data()),"RECREATE");	
	TTree *JetTree = new TTree("JetTree", "JetTree");	
	// Event-level branches
	OutTree->Branch("NEVENTS", &NEVENTS, "NEVENTS/I"); 
	// Reco Jet Branches
	OutTree->Branch("RecoJet_pt", &RecoJet_pt);
	OutTree->Branch("RecoJet_eta", &RecoJet_eta);
	OutTree->Branch("RecoJet_phi", &RecoJet_phi);
	OutTree->Branch("RecoJet_E", &RecoJet_E);
	OutTree->Branch("RecoJet_M", &RecoJet_M);
	OutTree->Branch("RecoJet_hasElectron", &RecoJet_hasElectron);
    OutTree->Branch("RecoJet_maxPtPart_pt", &RecoJet_maxPtPart_pt); 
	// Gen Jet Branches
	OutTree->Branch("GenJet_pt", &GenJet_pt);
	OutTree->Branch("GenJet_eta", &GenJet_eta);
	OutTree->Branch("GenJet_phi", &GenJet_phi);
	OutTree->Branch("GenJet_E", &GenJet_E);
	OutTree->Branch("GenJet_M", &GenJet_M);
	OutTree->Branch("GenJet_hasElectron", &GenJet_hasElectron);
    OutTree->Branch("GenJet_maxPtPart_pt", &GenJet_maxPtPart_pt); 
	// --- END TTree setup ---
	
	// Loop over events	
	while(tree_reader->Next()) {	
	
	    if(NEVENTS%50000 == 0) cout << "Events Processed: " << NEVENTS << endl;

		// Clear all vectors for the new event
		RecoJet_pt.clear();
		RecoJet_eta.clear();
		RecoJet_phi.clear();
		RecoJet_E.clear();
		RecoJet_M.clear();
		RecoJet_hasElectron.clear();
		RecoJet_maxPtPart_pt.clear();
		
		GenJet_pt.clear();
		GenJet_eta.clear();
		GenJet_phi.clear();
		GenJet_E.clear();
		GenJet_M.clear();
		GenJet_hasElectron.clear();
		GenJet_maxPtPart_pt.clear();

	    // Analyze Reconstructed Jets
		for(unsigned int ijet = 0; ijet < JetRecoType->GetSize(); ijet++) {
			// Make 4-vector
			LorentzVector JetReco((*JetRecoPx)[ijet], (*JetRecoPy)[ijet], (*JetRecoPz)[ijet], (*JetRecoE)[ijet]);
			RecoJet_pt.push_back(RecoJetPt);
			RecoJet_eta.push_back(JetReco.Eta());
			RecoJet_phi.push_back(JetReco.Phi());
			RecoJet_E.push_back(JetReco.E());
			RecoJet_M.push_back((*JetRecoM)[ijet]);

			bool hasElectron = false; // Check if Jet Contains an Electron - Use Particle Matching to Find True PID
			float maxPtReco = -1.0; // check max track pT in the constituents
			for(unsigned int icjet = (*JetRecoCBegin)[ijet]; icjet < (*JetRecoCEnd)[ijet]; icjet++) {// Loop over jet constituents (particles within the jet)
		   		int chargePartIndex = (*JetRecoCIdx)[icjet]; // ReconstructedChargedParticle Index for m'th Jet Component
				// Calculate constituent pT
				float px = (*TrkRecoPx)[chargePartIndex];
				float py = (*TrkRecoPy)[chargePartIndex];
				float pt = std::sqrt(px*px + py*py);
				// Update Max Pt Particle
				if (pt > maxPtReco) { maxPtReco = pt; }
				// Find electron in a jet
			    int elecIndex = -1;
			    float elecIndexWeight = -1.0;
		    	for(unsigned int itrkass = 0; itrkass < TrkPartAssocRec->GetSize(); itrkass++){ // Loop Over All ReconstructedChargedParticleAssociations
					if((*TrkPartAssocRec)[itrkass] == chargePartIndex){ // Select Entry Matching the ReconstructedChargedParticle Index
					    if((*TrkPartAssocWeight)[itrkass] > elecIndexWeight){ // Find Particle with Greatest Weight = Contributed Most Hits to Track
							elecIndex = (*TrkPartAssocSim)[itrkass]; // Get Index of MCParticle Associated with ReconstructedChargedParticle
							elecIndexWeight = (*TrkPartAssocWeight)[itrkass];
			      		}
			  		}
		      	}
				if((*TrkMCGenPDG)[elecIndex] == 11) hasElectron = true;
	  		}

			RecoJet_hasElectron.push_back(hasElectron);
			RecoJet_maxPtPart_pt.push_back(maxPtReco);
    
		}
		
		// Analyze Gen Jets
		for(unsigned int igjet = 0; igjet < JetGenType->GetSize(); igjet++) {
			LorentzVector JetGen((*JetGenPx)[igjet], (*JetGenPy)[igjet], (*JetGenPz)[igjet], (*JetGenE)[igjet]);
			GenJet_pt.push_back(JetGen.Pt());
			GenJet_eta.push_back(JetGen.Eta());
			GenJet_phi.push_back(JetGen.Phi());
			GenJet_E.push_back(JetGen.E());
			GenJet_M.push_back((*JetGenM)[igjet]);
			// -> Check for electrons
			// Loop over jet constituents (particles within the jet)
			bool hasGenElectron = false; 
			float maxPtGen = -1.0;
			for(unsigned int icgjet = (*JetGenCBegin)[igjet]; icgjet < (*JetGenCEnd)[igjet]; icgjet++) { 
				int genPartIndex = (*JetGenCIdx)[icgjet];
				// Calculate constituent pT
				float px = (*TrkGenPx)[genPartIndex];
				float py = (*TrkGenPy)[genPartIndex];
				float pt = std::sqrt(px*px + py*py);
				// Update Max Pt Particle
				if (pt > maxPtGen) { maxPtGen = pt; }
				int gTrkPDG = (*TrkGenPDG)[genPartIndex];	    		
			    if(gTrkPDG == 11) hasGenElectron = true;					
			}
			
			GenJet_hasElectron.push_back(hasGenElectron);
			GenJet_maxPtPart_pt.push_back(maxPtGen);			

		}

		NEVENTS++;

		OutTree->Fill();
		
	}
	
	cout << "Total number of events: " << NEVENTS << endl;
	OutFile->Write(); 	
	OutFile->Close();
	

}

