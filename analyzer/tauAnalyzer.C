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
  TH1F* h_lep_DR[7];
  TH1F* h_leptau_DR[7];

  for(int i=0; i < 7; i++){
    h_tauTag_matched[i] = new TH1F(Form("h_nTauTag_matched_S%i", i),";#Tau Tagged",2,0,2);
    h_jets_n[i] = new TH1F(Form("h_nJets_S%i", i),";Jet Multiplicity",10,0,10);
    h_bJets_n[i] = new TH1F(Form("h_nbJets_S%i", i),";bJet Multiplicity",5,0,5);
    h_tauJets_n[i] = new TH1F(Form("h_ntauJets_S%i", i),";#Tau Jet Multiplicity",5,0,5);
    h_lep_DR[i] = new TH1F(Form("h_lep_DR_S%i", i), ";#Delta R_{ll}", 40, 0, 4);
    h_leptau_DR[i] = new TH1F(Form("h_leptau_DR_S%i", i), ";#Delta R_{#Tau l}", 40, 0, 4);

    h_tauTag_matched[i]->Sumw2();
    h_jets_n[i]->Sumw2();
    h_bJets_n[i]->Sumw2();
    h_tauJets_n[i]->Sumw2();
    h_lep_DR[i]->Sumw2();
    h_leptau_DR[i]->Sumw2();
  }

  //Event Loop
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;

  for(Long64_t jentry=0; jentry<nentries; jentry++) {
    fChain->GetEntry(jentry);
    std::cout << jentry << " / " << nentries << '\r';

    std::vector<int> GoodMuIdx, GoodElecIdx;
    int nGoodMuon = 0;
    int nGoodElectron= 0;
    int tau_idx;

    for(int i = 0; i < nMuon; i++)
      if(Muon_pt[i] > 30 && std::abs(Muon_eta[i]) < 2.4) GoodMuIdx.push_back(i);
    for(int i = 0; i < nElectron; i++)
      if(Electron_pt[i] > 35 && std::abs(Electron_eta[i]) < 2.4) GoodElecIdx.push_back(i);

    nGoodMuon = GoodMuIdx.size();
    nGoodElectron = GoodElecIdx.size();

    std::vector<int> matchedIdx, tauIdx;
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
        if(Jet_tauTag[i] == 1){
          tauJets_n++;
          tauIdx.push_back(i);
        }
      }
    }

    //Selection Loop
    bool eventSelection[nCuts];
    for(int bcut = 0; bcut < nCuts; bcut++) eventSelection[bcut] = false;

    eventSelection[0] = true;
    eventSelection[1] = nGoodMuon + nGoodElectron == 1;
    eventSelection[2] = nGoodMuon + nGoodElectron == 2;
    eventSelection[3] = ( nGoodElectron + nGoodMuon >= 3 );
    eventSelection[4] = ( (nGoodMuon + nGoodElectron == 1) && jets_n >= 2 );
    eventSelection[5] = ( (nGoodMuon + nGoodElectron == 2) && jets_n >= 2 && bJets_n == 1 );
    eventSelection[6] = ( (nGoodMuon + nGoodElectron == 1) && jets_n >= 2 && tauJets_n == 1 );

    TLorentzVector lepton[nGoodMuon + nGoodElectron];

    for(int cut = 0; cut < nCuts; cut++) {
      if( eventSelection[cut] ) {
        h_jets_n[cut]->Fill(jets_n);
        h_bJets_n[cut]->Fill(bJets_n);
        h_tauJets_n[cut]->Fill(tauJets_n);

        for(int i = 0; i < nGoodMuon; i++){
          int idx = GoodMuIdx.at(i);
          lepton[i].SetPtEtaPhiM(Muon_pt[idx], Muon_eta[idx], Muon_phi[idx], Muon_m[idx]);
        }
        for(int i = nGoodMuon; i < nGoodMuon + nGoodElectron; i++){
          int idx = GoodElecIdx.at(i-nGoodMuon);
          lepton[i].SetPtEtaPhiM(Electron_pt[idx], Electron_eta[idx], Electron_phi[idx], Electron_m[idx]);
        }

				for(vector<int>::iterator iter=matchedIdx.begin(); iter!=matchedIdx.end(); ++iter){
          b_matched_tauTag.push_back(Jet_tauTag[*iter]);
          h_tauTag_matched[cut]->Fill(Jet_tauTag[*iter]);
        }
        b_matched_tauTag.clear();

        TLorentzVector tau;
        if( cut == 2 || cut == 5 ) h_lep_DR[cut]->Fill(lepton[0].DeltaR(lepton[1]));
        if( cut == 6 ){
          for(vector<int>::iterator iter=tauIdx.begin(); iter!=tauIdx.end(); ++iter){
            tau.SetPtEtaPhiM(Jet_pt[*iter], Jet_phi[*iter], Jet_eta[*iter], Jet_m[*iter]);
          }
          h_leptau_DR[cut]->Fill(lepton[0].DeltaR(tau));
        }
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
    h_lep_DR[i]->Write();
    h_leptau_DR[i]->Write();
  }
  f->Close();
}
