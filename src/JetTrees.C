#include "JetTreeReader.h"  // call libraries from ROOT and C++

void JetTrees(TString InputFileList, TString OutputFile){

	typedef ROOT::Math::PxPyPzEVector LorentzVector;

	int NEVENTS = 0;
	int EVETMULTRECO = 0; 
	int EVETMULTGEN = 0;
	// Reco Jets (Variable-length vectors for multiple jets per event)
	std::vector<float> RecoJet_pt;
	std::vector<float> RecoJet_eta;
	std::vector<float> RecoJet_phi;
	std::vector<float> RecoJet_E;
	std::vector<float> RecoJet_M;
	std::vector<bool> RecoJet_hasElectron;
	std::vector<float> RecoJet_maxPtPart_pt; 
	// Reco jet constituents
	std::vector<std::vector<float>> RecoJet_constituent_pt;
	std::vector<std::vector<float>> RecoJet_constituent_eta;
	std::vector<std::vector<float>> RecoJet_constituent_phi;
	std::vector<std::vector<int>> RecoJet_constituent_nhits;

	// Gen Jets (Variable-length vectors for multiple jets per event)
	std::vector<float> GenJet_pt;
	std::vector<float> GenJet_eta;
	std::vector<float> GenJet_phi;
	std::vector<float> GenJet_E;
	std::vector<float> GenJet_M;
	std::vector<bool> GenJet_hasElectron;
	std::vector<bool> GenJet_hasNeutral;
	std::vector<float> GenJet_maxPtPart_pt;
	// Gen jet constituents
	std::vector<std::vector<float>> GenJet_constituent_pt;
	std::vector<std::vector<float>> GenJet_constituent_eta;
	std::vector<std::vector<float>> GenJet_constituent_phi;

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
	
	// Make Output
	TFile *OutFile = TFile::Open(Form("%s",OutputFile.Data()),"RECREATE");	
	TTree *JetTree = new TTree("JetTree", "JetTree");	
	// Event-level branches
	JetTree->Branch("NEVENTS", &NEVENTS, "NEVENTS/I"); 
	JetTree->Branch("EVETMULTRECO", &EVETMULTRECO, "EVETMULTRECO/I");
	JetTree->Branch("EVETMULTGEN", &EVETMULTGEN, "EVETMULTGEN/I");
	// Reco Jet Branches
	JetTree->Branch("RecoJet_pt", &RecoJet_pt);
	JetTree->Branch("RecoJet_eta", &RecoJet_eta);
	JetTree->Branch("RecoJet_phi", &RecoJet_phi);
	JetTree->Branch("RecoJet_E", &RecoJet_E);
	JetTree->Branch("RecoJet_M", &RecoJet_M);
	JetTree->Branch("RecoJet_hasElectron", &RecoJet_hasElectron);
    JetTree->Branch("RecoJet_maxPtPart_pt", &RecoJet_maxPtPart_pt); 
	// Reco Jet constituent branches
	JetTree->Branch("RecoJet_constituent_pt", &RecoJet_constituent_pt);
	JetTree->Branch("RecoJet_constituent_eta", &RecoJet_constituent_eta);
	JetTree->Branch("RecoJet_constituent_phi", &RecoJet_constituent_phi);
	JetTree->Branch("RecoJet_constituent_nhits", &RecoJet_constituent_nhits);

	// Gen Jet Branches
	JetTree->Branch("GenJet_pt", &GenJet_pt);
	JetTree->Branch("GenJet_eta", &GenJet_eta);
	JetTree->Branch("GenJet_phi", &GenJet_phi);
	JetTree->Branch("GenJet_E", &GenJet_E);
	JetTree->Branch("GenJet_M", &GenJet_M);
	JetTree->Branch("GenJet_hasElectron", &GenJet_hasElectron);
	JetTree->Branch("GenJet_hasNeutral", &GenJet_hasNeutral);
    JetTree->Branch("GenJet_maxPtPart_pt", &GenJet_maxPtPart_pt); 
	// Gen Jet constituent branches
	JetTree->Branch("GenJet_constituent_pt", &GenJet_constituent_pt);
	JetTree->Branch("GenJet_constituent_eta", &GenJet_constituent_eta);
	JetTree->Branch("GenJet_constituent_phi", &GenJet_constituent_phi);

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
        RecoJet_constituent_pt.clear(); 
        RecoJet_constituent_eta.clear();
        RecoJet_constituent_phi.clear(); 
        RecoJet_constituent_nhits.clear();
		
		GenJet_pt.clear();
		GenJet_eta.clear();
		GenJet_phi.clear();
		GenJet_E.clear();
		GenJet_M.clear();
		GenJet_hasElectron.clear();
		GenJet_hasNeutral.clear();
		GenJet_maxPtPart_pt.clear();
        GenJet_constituent_pt.clear(); 
        GenJet_constituent_eta.clear();
        GenJet_constituent_phi.clear(); 

		EVETMULTRECO = TrkRecoPx->GetSize();
		EVETMULTGEN = TrkGenPx->GetSize();

	    // Analyze Reconstructed Jets
		for(unsigned int ijet = 0; ijet < JetRecoType->GetSize(); ijet++) {
			// Make 4-vector
			LorentzVector JetReco((*JetRecoPx)[ijet], (*JetRecoPy)[ijet], (*JetRecoPz)[ijet], (*JetRecoE)[ijet]);
			RecoJet_pt.push_back(JetReco.Pt());
			RecoJet_eta.push_back(JetReco.Eta());
			RecoJet_phi.push_back(JetReco.Phi());
			RecoJet_E.push_back(JetReco.E());
			RecoJet_M.push_back((*JetRecoM)[ijet]);

			std::vector<float> const_pt;
			std::vector<float> const_eta;
			std::vector<float> const_phi;
			std::vector<int> const_nhits;
            const_pt.clear();
            const_eta.clear();
            const_phi.clear();
            const_nhits.clear();

			bool hasElectron = false; // Check if Jet Contains an Electron - Use Particle Matching to Find True PID
			float maxPtReco = -1.0; // check max track pT in the constituents
			for(unsigned int icjet = (*JetRecoCBegin)[ijet]; icjet < (*JetRecoCEnd)[ijet]; icjet++) {// Loop over jet constituents (particles within the jet)
		   		int chargePartIndex = (*JetRecoCIdx)[icjet]; // ReconstructedChargedParticle Index for m'th Jet Component
				// Calculate constituent pT
				TVector3 TrkVec((*TrkRecoPx)[chargePartIndex], (*TrkRecoPy)[chargePartIndex], (*TrkRecoPz)[chargePartIndex]);
				const_pt.push_back(TrkVec.Pt());
                const_eta.push_back(TrkVec.Eta());
                const_phi.push_back(TrkVec.Phi());
                const_nhits.push_back((*TrkRecoNhits)[chargePartIndex]);
				// Update Max Pt Particle
				if (TrkVec.Pt() > maxPtReco) { maxPtReco = TrkVec.Pt(); }
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
            RecoJet_constituent_pt.push_back(const_pt);
            RecoJet_constituent_eta.push_back(const_eta);
            RecoJet_constituent_phi.push_back(const_phi);   
            RecoJet_constituent_nhits.push_back(const_nhits);   

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
            std::vector<float> gconst_pt;
            std::vector<float> gconst_eta;
            std::vector<float> gconst_phi;
            gconst_pt.clear();
            gconst_eta.clear();
            gconst_phi.clear();   
            
			bool hasGenElectron = false; 
			bool hasGenNeutral = false; 
			float maxPtGen = -1.0;
			for(unsigned int icgjet = (*JetGenCBegin)[igjet]; icgjet < (*JetGenCEnd)[igjet]; icgjet++) { 
				int genPartIndex = (*JetGenCIdx)[icgjet];
				// Calculate constituent pT
				TVector3 gTrkVec((*TrkGenPx)[genPartIndex], (*TrkGenPy)[genPartIndex], (*TrkGenPz)[genPartIndex]);
                gconst_pt.push_back(gTrkVec.Pt());
                gconst_eta.push_back(gTrkVec.Eta());
                gconst_phi.push_back(gTrkVec.Phi());
				float charge = (*TrkGenCharge)[genPartIndex];
				// Update Max Pt Particle
				if (gTrkVec.Pt() > maxPtGen) { maxPtGen = gTrkVec.Pt(); }
				int gTrkPDG = (*TrkGenPDG)[genPartIndex];	    		
			    if(gTrkPDG == 11) hasGenElectron = true;
			    if(charge == 0)	hasGenNeutral = true;		
			}
			
			GenJet_hasElectron.push_back(hasGenElectron);
			GenJet_hasNeutral.push_back(hasGenNeutral);
			GenJet_maxPtPart_pt.push_back(maxPtGen);			
            GenJet_constituent_pt.push_back(gconst_pt);
            GenJet_constituent_eta.push_back(gconst_eta);
            GenJet_constituent_phi.push_back(gconst_phi);         
            
		}

		NEVENTS++;

		JetTree->Fill();
		
	}
	
	cout << "Total number of events: " << NEVENTS << endl;
	OutFile->Write(); 	
	OutFile->Close();
	

}

