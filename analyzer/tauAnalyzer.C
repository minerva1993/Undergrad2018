#define tauAnalyzer_cxx
#include "tauAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
using namespace std;

void tauAnalyzer::Loop(const std::string outFileName)
{
//   In a ROOT session, you can do:
//      root> .L tauAnalyzer.C
//      root> tauAnalyzer t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

  if (fChain == 0) return;

  TFile* f = new TFile(outFileName.c_str(), "recreate");
  TTree* tree = new TTree("tree", "ntuple");

  std::vector<int> b_matched_tauTag;

	int nCuts = 7;
	int jets_n, bJets_n, tauJets_n;

  tree->Branch("matched_tauTag", "std::vector<int>", &b_matched_tauTag);
 
	TH1F* h_tauTag_matched[7];
	TH1F* h_jets_n[7];
	TH1F* h_bJets_n[7];
	TH1F* h_tauJets_n[7];

	for(int i=0; i < 7; i++){
		h_tauTag_matched[i] = new TH1F(Form("nTauTag_matched_S%i", i),"",2,0,2);
		h_jets_n[i] = new TH1F(Form("nJets_S%i", i),"",10,0,10);
		h_bJets_n[i] = new TH1F(Form("nbJets_S%i", i),"",5,0,5);
		h_tauJets_n[i] = new TH1F(Form("ntauJets_S%i", i),"",5,0,5);

		h_tauTag_matched[i]->Sumw2();
		h_jets_n[i]->Sumw2();
		h_bJets_n[i]->Sumw2();
		h_tauJets_n[i]->Sumw2();
	}

	//Event Loop
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;

  for(Long64_t jentry=0; jentry<nentries; jentry++) {
    fChain->GetEntry(jentry);
    std::cout << jentry << " / " << nentries << '\r';

    vector<int> matchedIdx;
    TLorentzVector jet;

		jets_n = 0;
		bJets_n = 0;
		tauJets_n = 0;

		for(int i=0; i<nJet; i++){
			if(Jet_pt[i] > 30 && abs(Jet_eta[i]) < 2.4){
				jet.SetPtEtaPhiM(Jet_pt[i],Jet_eta[i],Jet_phi[i], Jet_m[i]);
				jets_n++;
 
				TLorentzVector genParticle;
        for(int j=0; j<nGenParticle; j++){
          if(std::abs(GenParticle_pdgId[j]) != 15) continue;
          genParticle.SetPtEtaPhiM(GenParticle_pt[j], GenParticle_eta[j], GenParticle_phi[j], GenParticle_m[j]);
          if(jet.DeltaR(genParticle) < 0.5) matchedIdx.push_back(j);
        }

				if(Jet_bTag[i] == 1) bJets_n++;
				if(Jet_tauTag[i] == 1) tauJets_n++;
			}
    }

		//Selection Loop
		bool eventSelection[nCuts];
		for(int bcut = 0; bcut < nCuts; bcut++) eventSelection[bcut] = false;

		eventSelection[0] = true;
		eventSelection[1] = nMuon + nElectron == 1;
		eventSelection[2] = nMuon + nElectron == 2;
		eventSelection[3] = ( nElectron + nMuon >= 3 );
		eventSelection[4] = ( (nMuon + nElectron == 1) && jets_n >= 2 );
		eventSelection[5] = ( (nMuon + nElectron == 2) && jets_n >= 2 && bJets_n == 1 );
		eventSelection[6] = ( (nMuon + nElectron == 1) && jets_n >= 2 && tauJets_n >= 1 );

		for(int cut = 0; cut < nCuts; cut++) {
			if( eventSelection[cut] ) {
				h_jets_n[cut]->Fill(jets_n);
				h_bJets_n[cut]->Fill(bJets_n);
				h_tauJets_n[cut]->Fill(tauJets_n);

				for(vector<int>::iterator iter=matchedIdx.begin(); iter!=matchedIdx.end(); ++iter){
					b_matched_tauTag.push_back(Jet_tauTag[*iter]);
					h_tauTag_matched[cut]->Fill(Jet_tauTag[*iter]);
				}
				b_matched_tauTag.clear();
			}
		}
		tree->Fill();
	}
 
	tree->Write();

	for(int i = 0; i < nCuts; i++){
		h_jets_n[i]->Write();
		h_bJets_n[i]->Write();
		h_tauJets_n[i]->Write();
		h_tauTag_matched[i]->Write();
	}
	f->Close();
}
