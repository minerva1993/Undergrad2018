#define tauAnalyzer_cxx
#include "tauAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
using namespace std;

void tauAnalyzer::Loop()
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

	std::vector<int> b_matched_tauTag;

   TFile* f = new TFile("out.root", "recreate");
   TTree* tree = new TTree("tree", "ntuple");
	tree->Branch("matched_tauTag", "std::vector<int>", &b_matched_tauTag);
	//TH1F* h_tauTag_matched = new TH1F("tauTag","",2,0,2);

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries; jentry++) {
		fChain->GetEntry(jentry);

		vector<int> matchedIdx;

		TLorentzVector jet;
		for(int i=0; i<nJet; i++){
			if(Jet_pt[i] > 30 && abs(Jet_eta[i]) < 2.4){
				jet.SetPtEtaPhiM(Jet_pt[i],Jet_eta[i],Jet_phi[i], Jet_m[i]);
			    TLorentzVector genParticle;
			    for(int j=0; j<nGenParticle; j++){
				  if(std::abs(GenParticle_pdgId[j]) != 15) continue;
				  genParticle.SetPtEtaPhiM(GenParticle_pt[j], GenParticle_eta[j], GenParticle_phi[j], GenParticle_m[j]);
				  if(jet.DeltaR(genParticle) < 0.5) matchedIdx.push_back(j);
			}
		  }
		}
	
		for(vector<int>::iterator iter=matchedIdx.begin(); iter!=matchedIdx.end(); ++iter){
			b_matched_tauTag.push_back(Jet_tauTag[*iter]);
		}

		tree->Fill();
        b_matched_tauTag.clear();

	}
	

//		if(muons_n*electrons_n!=0) continue;
/*		
   	for(Int_t j=0; j<jets_n; j++){
         if(jets_pt[j]>30 && abs(jets_eta[j])<2.5){
				nJets++;
				if(jets_bTag[j]==1) nBJets++;
				if(nJets==1) jet1.SetPtEtaPhiM(jets_pt[j], jets_eta[j], jets_phi[j], jets_m[j]);
				if(nJets==2){
					jet2.SetPtEtaPhiM(jets_pt[j], jets_eta[j], jets_phi[j], jets_m[j]);
					dijet_m=(jet1+jet2).M();	
				}
				j_pt.push_back(jets_pt[j]);
				j_eta.push_back(jets_eta[j]);
 	 		}
		}
*/
	tree->Write();
   f->Close();
}
