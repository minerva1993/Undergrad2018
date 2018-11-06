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

  int nCuts = 13;
  int nCh = 8;
  int jets_n, bJets_n, cJets_n, tauJets_n;

  tree->Branch("matched_tauTag", "std::vector<int>", &b_matched_tauTag);
 
  TH1F* h_tauTag_matched[8][13];
  TH1F* h_jets_n[8][13];
  TH1F* h_bJets_n[8][13];
  TH1F* h_cJets_n[8][13];
  TH1F* h_tauJets_n[8][18];
  TH1F* h_lep_DR[8][13];
  TH1F* h_leptau_DR[8][13];

  for(int ch=0; ch < nCh; ch++){
    for(int i=0; i < nCuts; i++){
      h_tauTag_matched[ch][i] = new TH1F(Form("h_nTauTag_matched_Ch%i_S%i",ch,i),";#Tau Tagged",2,0,2);
      h_jets_n[ch][i] = new TH1F(Form("h_nJets_Ch%i_S%i",ch,i),";Jet Multiplicity",10,0,10);
      h_bJets_n[ch][i] = new TH1F(Form("h_nbJets_Ch%i_S%i",ch,i),";bJet Multiplicity",5,0,5);
      h_cJets_n[ch][i] = new TH1F(Form("h_ncJets_Ch%i_S%i",ch,i),";bJet Multiplicity",5,0,5);
      h_tauJets_n[ch][i] = new TH1F(Form("h_ntauJets_Ch%i_S%i",ch,i),";#Tau Jet Multiplicity",5,0,5);
      h_lep_DR[ch][i] = new TH1F(Form("h_lep_DR_Ch%i_S%i",ch,i), ";#Delta R_{ll}", 40, 0, 4);
      h_leptau_DR[ch][i] = new TH1F(Form("h_leptau_DR_Ch%i_S%i",ch,i), ";#Delta R_{#Tau l}", 40, 0, 4);

      h_tauTag_matched[ch][i]->Sumw2();
      h_jets_n[ch][i]->Sumw2();
      h_bJets_n[ch][i]->Sumw2();
      h_cJets_n[ch][i]->Sumw2();
      h_tauJets_n[ch][i]->Sumw2();
      h_lep_DR[ch][i]->Sumw2();
      h_leptau_DR[ch][i]->Sumw2();
    }
  }

  TH1F* EventInfo;
  EventInfo = new TH1F("EventInfo","EventInfo",2,0,2);
  EventInfo->GetXaxis()->SetBinLabel(1,"Number of Events");
  EventInfo->GetXaxis()->SetBinLabel(2,"Sum of Weights");


  //Event Loop
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;

  for(Long64_t jentry=0; jentry<nentries; jentry++) {
    fChain->GetEntry(jentry);
    //std::cout << jentry << " / " << nentries << '\r';

    EventInfo->Fill(0.5, 1.0);
    std::vector<int> GoodMuIdx, GoodElecIdx;
    int nGoodMuon = 0;
    int nGoodElectron= 0;
    int tau_idx;

    for(int i = 0; i < nMuon; i++){
      if(GoodMuIdx.size() == 0 && Muon_pt[i] > 30 && std::abs(Muon_eta[i]) < 2.4 && Muon_relIso[i] < 0.15)
        GoodMuIdx.push_back(i);
      else if(GoodMuIdx.size() > 0 && Muon_pt[i] > 25 && std::abs(Muon_eta[i]) < 2.4 && Muon_relIso[i] < 0.15)
        GoodMuIdx.push_back(i);
    }
    for(int i = 0; i < nElectron; i++){
      if(GoodElecIdx.size() == 0 && Electron_pt[i] > 35 && std::abs(Electron_eta[i]) < 2.4 && Electron_relIso[i] < 0.15)
        GoodElecIdx.push_back(i);
      else if(GoodElecIdx.size() > 0 && Electron_pt[i] > 30 && std::abs(Electron_eta[i]) < 2.4 && Electron_relIso[i] < 0.15)
        GoodElecIdx.push_back(i);
    }

    nGoodMuon = GoodMuIdx.size();
    nGoodElectron = GoodElecIdx.size();

    std::vector<int> matchedIdx, tauIdx;
    TLorentzVector jet;

    jets_n = 0;
    bJets_n = 0;
    cJets_n = 0;
    tauJets_n = 0;

    //lepton selection
    bool pass_lep[8];
    pass_lep[0] = (nGoodMuon==1) && (nGoodElectron==0);
    pass_lep[1] = (nGoodMuon==0) && (nGoodElectron==1);
    pass_lep[2] = (pass_lep[0] || pass_lep[1]); 
    pass_lep[3] = (nGoodMuon==2) && (nGoodElectron==0);
    pass_lep[4] = (nGoodMuon==0) && (nGoodElectron==2);
    pass_lep[5] = (nGoodMuon==1) && (nGoodElectron==1);
    pass_lep[6] = (pass_lep[3] || pass_lep[4] || pass_lep[5]);
    pass_lep[7] = (nGoodMuon+nGoodElectron==3);

    //Jet multiplicity
    for(int i=0; i<nJet; i++){
      if(Jet_pt[i] > 30 && abs(Jet_eta[i]) < 2.4){
        jet.SetPtEtaPhiM(Jet_pt[i],Jet_eta[i],Jet_phi[i], Jet_m[i]);
        jets_n++;

        TLorentzVector genParticle;
        for(int j=0; j<nGenParticle; j++){
          if(std::abs(GenParticle_pdgId[j]) != 15 and Jet_tauTag[i] != 1) continue;
          genParticle.SetPtEtaPhiM(GenParticle_pt[j], GenParticle_eta[j], GenParticle_phi[j], GenParticle_m[j]);
          if(jet.DeltaR(genParticle) < 0.5) matchedIdx.push_back(j);
        }

        if(Jet_bTag[i] == 1) bJets_n++;
        if(Jet_cTag[i] == 2) cJets_n++;
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
    eventSelection[1] = jets_n >= 2;
    eventSelection[2] = jets_n >= 2 && bJets_n == 1;
    eventSelection[3] = jets_n >= 2 && tauJets_n == 1;
    eventSelection[4] = jets_n >= 2 && tauJets_n == 2;
    eventSelection[5] = jets_n >= 2 && bJets_n == 1 && tauJets_n == 1;
    eventSelection[6] = jets_n >= 2 && bJets_n == 1 && tauJets_n == 2;
    eventSelection[7] = jets_n >= 3;
    eventSelection[8] = jets_n >= 3 && bJets_n == 1;
    eventSelection[9] = jets_n >= 3 && tauJets_n == 1;
    eventSelection[10] = jets_n >= 3 && tauJets_n == 2;
    eventSelection[11] = jets_n >= 3 && bJets_n == 1 && tauJets_n == 1;
    eventSelection[12] = jets_n >= 3 && bJets_n == 1 && tauJets_n == 2;

    TLorentzVector lepton[nGoodMuon + nGoodElectron];

    for(int MODE = 0; MODE< 8; MODE++){
      for(int cut = 0; cut < nCuts; cut++) {
        if( eventSelection[cut] && pass_lep[MODE] ) {
          h_jets_n[MODE][cut]->Fill(jets_n);
          h_bJets_n[MODE][cut]->Fill(bJets_n);
          h_cJets_n[MODE][cut]->Fill(cJets_n);
          h_tauJets_n[MODE][cut]->Fill(tauJets_n);

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
            h_tauTag_matched[MODE][cut]->Fill(Jet_tauTag[*iter]);
          }
          b_matched_tauTag.clear();

          TLorentzVector tau;
          if( MODE >= 3 ) h_lep_DR[MODE][cut]->Fill(lepton[0].DeltaR(lepton[1]));
/*
          if( cut == 6 ){
            for(vector<int>::iterator iter=tauIdx.begin(); iter!=tauIdx.end(); ++iter){
              tau.SetPtEtaPhiM(Jet_pt[*iter], Jet_phi[*iter], Jet_eta[*iter], Jet_m[*iter]);
            }
            h_leptau_DR[MODE][cut]->Fill(lepton[0].DeltaR(tau));
          }
*/
        }
      }//selection loop
    }
    tree->Fill();
  }
  tree->Write();

  for(int ch=0; ch < nCh; ch++){
    for(int i = 0; i < nCuts; i++){
      h_jets_n[ch][i]->Write();
      h_bJets_n[ch][i]->Write();
      h_cJets_n[ch][i]->Write();
      h_tauJets_n[ch][i]->Write();
      h_tauTag_matched[ch][i]->Write();
      h_lep_DR[ch][i]->Write();
      h_leptau_DR[ch][i]->Write();
    }
  }
  EventInfo->Write();
  f->Close();
}
