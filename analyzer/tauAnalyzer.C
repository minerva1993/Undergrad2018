#define tauAnalyzer_cxx
#include "tauAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
using namespace std;

void tauAnalyzer::Loop(const std::string outFileName)
{

  if (fChain == 0) return;

  TFile* f = new TFile(outFileName.c_str(), "recreate");
  TTree* tree = new TTree("tree", "ntuple");

  std::vector<int> b_matched_tauTag;

  int nCuts = 83;
  int nCh = 11;
  int jets_n, bJets_n, cJets_n, tauJets_n;
  float jet_pt_sum;

  tree->Branch("matched_tauTag", "std::vector<int>", &b_matched_tauTag);
 
  TH1F* h_tauTag_matched[11][83];
  TH1F* h_jets_n[11][83];
  TH1F* h_bJets_n[11][83];
  TH1F* h_cJets_n[11][83];
  TH1F* h_tauJets_n[11][83];
  TH1F* h_lepDR[11][83];
  TH1F* h_leptau_DR[11][83];
  TH1F* h_tautauDR[11][83];
  TH1F* h_MET[11][83];
  TH1F* h_jetptsum[11][83];
  TH1F* h_smTop[11][83];
  TH1F* h_sigTop[11][83];

  for(int ch=0; ch < nCh; ch++){
    for(int i=0; i < nCuts; i++){
      h_tauTag_matched[ch][i] = new TH1F(Form("h_nTauTagmatched_Ch%i_S%i",ch,i),";#tau Tagged",2,0,2);
      h_jets_n[ch][i] = new TH1F(Form("h_nJets_Ch%i_S%i",ch,i),";Jet Multiplicity",10,0,10);
      h_bJets_n[ch][i] = new TH1F(Form("h_nbJets_Ch%i_S%i",ch,i),";b jet Multiplicity",5,0,5);
      h_cJets_n[ch][i] = new TH1F(Form("h_ncJets_Ch%i_S%i",ch,i),";c jet Multiplicity",5,0,5);
      h_tauJets_n[ch][i] = new TH1F(Form("h_ntauJets_Ch%i_S%i",ch,i),";#tau Jet Multiplicity",5,0,5);
      h_lepDR[ch][i] = new TH1F(Form("h_lepDR_Ch%i_S%i",ch,i), ";#Delta R_{ll}", 40, 0, 4);
      h_leptau_DR[ch][i] = new TH1F(Form("h_leptauDR_Ch%i_S%i",ch,i), ";#Delta R_{#tau l}", 40, 0, 4);
      h_tautauDR[ch][i] = new TH1F(Form("h_tautauDR_Ch%i_S%i",ch,i), ";#Delta R_{#tau #tau}", 40, 0, 4);
      h_MET[ch][i] = new TH1F(Form("h_MET_Ch%i_S%i",ch,i), ";MET (GeV)", 40, 0, 300);
      h_jetptsum[ch][i] = new TH1F(Form("h_jetPtsum_Ch%i_S%i",ch,i), ";H_{T} (GeV)",40,0,1000);
      h_smTop[ch][i] = new TH1F(Form("h_smTop_Ch%i_S%i",ch,i), ";SM Top Mass (GeV)",40,0,300);
      h_sigTop[ch][i] = new TH1F(Form("h_sigTop_Ch%i_S%i",ch,i), ";Signal Top Mass (GeV)",40,0,300);

      h_tauTag_matched[ch][i]->Sumw2();
      h_jets_n[ch][i]->Sumw2();
      h_bJets_n[ch][i]->Sumw2();
      h_cJets_n[ch][i]->Sumw2();
      h_tauJets_n[ch][i]->Sumw2();
      h_lepDR[ch][i]->Sumw2();
      h_leptau_DR[ch][i]->Sumw2();
      h_MET[ch][i]->Sumw2();
      h_jetptsum[ch][i]->Sumw2();
      h_smTop[ch][i]->Sumw2();
      h_sigTop[ch][i]->Sumw2();
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
      if(GoodElecIdx.size() == 0 && Electron_pt[i] > 35 && std::abs(Electron_eta[i]) < 2.4 && Electron_relIso[i] < 0.25)
        GoodElecIdx.push_back(i);
      else if(GoodElecIdx.size() > 0 && Electron_pt[i] > 30 && std::abs(Electron_eta[i]) < 2.4 && Electron_relIso[i] < 0.25)
        GoodElecIdx.push_back(i);
    }

    nGoodMuon = GoodMuIdx.size();
    nGoodElectron = GoodElecIdx.size();

    std::vector<int> matchedIdx, tauIdx, bIdx, cIdx, jetIdx;
    TLorentzVector jet;

    jets_n = 0;
    bJets_n = 0;
    cJets_n = 0;
    tauJets_n = 0;
    jet_pt_sum = 0;

    //lepton selection
    bool pass_lep[nCh];
    pass_lep[0] = (nGoodMuon==1) && (nGoodElectron==0);
    pass_lep[1] = (nGoodMuon==0) && (nGoodElectron==1);
    pass_lep[2] = (pass_lep[0] || pass_lep[1]); 
    pass_lep[3] = (nGoodMuon==2) && (nGoodElectron==0);
    pass_lep[4] = (nGoodMuon==0) && (nGoodElectron==2);
    pass_lep[5] = (nGoodMuon==1) && (nGoodElectron==1);
    pass_lep[6] = (pass_lep[3] || pass_lep[4] || pass_lep[5]);
    pass_lep[7] = (nGoodMuon+nGoodElectron == 3);
    pass_lep[8] = (nGoodMuon+nGoodElectron >= 1);
    pass_lep[9] = (nGoodMuon+nGoodElectron >= 2);
    pass_lep[10] = (nGoodMuon+nGoodElectron >= 2 && nGoodMuon > 0);

    //Jet multiplicity
    for(int i=0; i<nJet; i++){
      if(Jet_pt[i] > 30 && abs(Jet_eta[i]) < 2.4){
        jet.SetPtEtaPhiM(Jet_pt[i],Jet_eta[i],Jet_phi[i], Jet_m[i]);
        jets_n++;
        jet_pt_sum += Jet_pt[i];
        jetIdx.push_back(i);

        TLorentzVector genParticle;
        for(int j=0; j<nGenParticle; j++){
          if(std::abs(GenParticle_pdgId[j]) != 15 and Jet_tauTag[i] != 1) continue;
          genParticle.SetPtEtaPhiM(GenParticle_pt[j], GenParticle_eta[j], GenParticle_phi[j], GenParticle_m[j]);
          if(jet.DeltaR(genParticle) < 0.5) matchedIdx.push_back(j);
        }

        if(Jet_bTag[i] == 1){
          bJets_n++;
          bIdx.push_back(i);
        }
        if(Jet_cTag[i] == 2){
          cJets_n++;
          cIdx.push_back(i);
        }
        if(Jet_tauTag[i] == 1){
          tauJets_n++;
          tauIdx.push_back(i);
        }
      }
    }

    //Selection Loop
    bool eventSelection[nCuts];
    for(int bcut = 0; bcut < nCuts; bcut++) eventSelection[bcut] = false;

    eventSelection[0] =  true;
    eventSelection[1] =  jets_n >= 2;
    eventSelection[2] =  jets_n >= 2 && bJets_n == 1;
    eventSelection[3] =  jets_n >= 2 && bJets_n >= 1;
    eventSelection[4] =  jets_n >= 2 && cJets_n == 1;
    eventSelection[5] =  jets_n >= 2 && cJets_n >= 1;
    eventSelection[6] =  jets_n >= 2 && tauJets_n == 1;
    eventSelection[7] =  jets_n >= 2 && tauJets_n == 2;
    eventSelection[8] =  jets_n >= 2 && tauJets_n >= 1;
    eventSelection[9] =  jets_n >= 2 && jet_pt_sum > 200;
    eventSelection[10] = jets_n >= 2 && bJets_n == 1 && tauJets_n == 1;
    eventSelection[11] = jets_n >= 2 && bJets_n == 1 && tauJets_n == 2;
    eventSelection[12] = jets_n >= 2 && bJets_n == 1 && tauJets_n >= 1;
    eventSelection[13] = jets_n >= 2 && bJets_n >= 1 && tauJets_n == 1;
    eventSelection[14] = jets_n >= 2 && bJets_n >= 1 && tauJets_n == 2;
    eventSelection[15] = jets_n >= 2 && bJets_n >= 1 && tauJets_n >= 1;
    eventSelection[16] = jets_n >= 2 && bJets_n == 1 && tauJets_n == 1 && jet_pt_sum > 200;
    eventSelection[17] = jets_n >= 2 && bJets_n == 1 && tauJets_n == 2 && jet_pt_sum > 200;
    eventSelection[18] = jets_n >= 2 && bJets_n == 1 && tauJets_n >= 1 && jet_pt_sum > 200;
    eventSelection[19] = jets_n >= 2 && bJets_n >= 1 && tauJets_n == 1 && jet_pt_sum > 200;
    eventSelection[20] = jets_n >= 2 && bJets_n >= 1 && tauJets_n == 2 && jet_pt_sum > 200;
    eventSelection[21] = jets_n >= 2 && bJets_n >= 1 && tauJets_n >= 1 && jet_pt_sum > 200;
    eventSelection[22] = jets_n >= 2 && bJets_n == 1 && cJets_n == 1;
    eventSelection[23] = jets_n >= 2 && bJets_n >= 1 && cJets_n == 1;
    eventSelection[24] = jets_n >= 2 && bJets_n == 1 && cJets_n >= 1;
    eventSelection[25] = jets_n >= 2 && bJets_n >= 1 && cJets_n >= 1;
    eventSelection[26] = jets_n >= 2 && bJets_n == 1 && cJets_n == 1 && jet_pt_sum > 200;
    eventSelection[27] = jets_n >= 2 && bJets_n >= 1 && cJets_n == 1 && jet_pt_sum > 200;
    eventSelection[28] = jets_n >= 2 && bJets_n == 1 && cJets_n >= 1 && jet_pt_sum > 200;
    eventSelection[29] = jets_n >= 2 && bJets_n >= 1 && cJets_n >= 1 && jet_pt_sum > 200;
    eventSelection[30] = jets_n >= 2 && cJets_n == 1 && tauJets_n == 1;
    eventSelection[31] = jets_n >= 2 && cJets_n == 1 && tauJets_n == 2;
    eventSelection[32] = jets_n >= 2 && cJets_n >= 1 && tauJets_n == 1;
    eventSelection[33] = jets_n >= 2 && cJets_n >= 1 && tauJets_n == 2;
    eventSelection[34] = jets_n >= 2 && cJets_n >= 1 && tauJets_n >= 1;
    eventSelection[35] = jets_n >= 2 && cJets_n >= 1 && tauJets_n >= 2;
    eventSelection[36] = jets_n >= 2 && bJets_n == 1 && MET_pt > 80;
    eventSelection[37] = jets_n >= 2 && cJets_n == 1 && MET_pt > 80;
    eventSelection[38] = jets_n >= 2 && bJets_n == 1 && cJets_n == 1 && MET_pt > 80;
    eventSelection[39] = jets_n >= 2 && bJets_n == 1 && MET_pt > 80 && jet_pt_sum > 200;
    eventSelection[40] = jets_n >= 2 && cJets_n == 1 && MET_pt > 80 && jet_pt_sum > 200;
    eventSelection[41] = jets_n >= 2 && bJets_n == 1 && cJets_n == 1 && MET_pt > 80 && jet_pt_sum > 200;
    eventSelection[42] = jets_n >= 3;
    eventSelection[43] = jets_n >= 3 && bJets_n == 1;
    eventSelection[44] = jets_n >= 3 && bJets_n >= 1;
    eventSelection[45] = jets_n >= 3 && cJets_n == 1;
    eventSelection[46] = jets_n >= 3 && cJets_n >= 1;
    eventSelection[47] = jets_n >= 3 && tauJets_n == 1;
    eventSelection[48] = jets_n >= 3 && tauJets_n == 2;
    eventSelection[49] = jets_n >= 3 && tauJets_n >= 1;
    eventSelection[50] = jets_n >= 3 && jet_pt_sum > 200;
    eventSelection[51] = jets_n >= 3 && bJets_n == 1 && tauJets_n == 1;
    eventSelection[52] = jets_n >= 3 && bJets_n == 1 && tauJets_n == 2;
    eventSelection[53] = jets_n >= 3 && bJets_n == 1 && tauJets_n >= 1;
    eventSelection[54] = jets_n >= 3 && bJets_n >= 1 && tauJets_n == 1;
    eventSelection[55] = jets_n >= 3 && bJets_n >= 1 && tauJets_n == 2;
    eventSelection[56] = jets_n >= 3 && bJets_n >= 1 && tauJets_n >= 1;
    eventSelection[57] = jets_n >= 3 && bJets_n == 1 && tauJets_n == 1 && jet_pt_sum > 200;
    eventSelection[58] = jets_n >= 3 && bJets_n == 1 && tauJets_n == 2 && jet_pt_sum > 200;
    eventSelection[59] = jets_n >= 3 && bJets_n == 1 && tauJets_n >= 1 && jet_pt_sum > 200;
    eventSelection[60] = jets_n >= 3 && bJets_n >= 1 && tauJets_n == 1 && jet_pt_sum > 200;
    eventSelection[61] = jets_n >= 3 && bJets_n >= 1 && tauJets_n == 2 && jet_pt_sum > 200;
    eventSelection[62] = jets_n >= 3 && bJets_n >= 1 && tauJets_n >= 1 && jet_pt_sum > 200;
    eventSelection[63] = jets_n >= 3 && bJets_n == 1 && cJets_n == 1;
    eventSelection[64] = jets_n >= 3 && bJets_n >= 1 && cJets_n == 1;
    eventSelection[65] = jets_n >= 3 && bJets_n == 1 && cJets_n >= 1;
    eventSelection[66] = jets_n >= 3 && bJets_n >= 1 && cJets_n >= 1;
    eventSelection[67] = jets_n >= 3 && bJets_n == 1 && cJets_n == 1 && jet_pt_sum > 200;
    eventSelection[68] = jets_n >= 3 && bJets_n >= 1 && cJets_n == 1 && jet_pt_sum > 200;
    eventSelection[69] = jets_n >= 3 && bJets_n == 1 && cJets_n >= 1 && jet_pt_sum > 200;
    eventSelection[70] = jets_n >= 3 && bJets_n >= 1 && cJets_n >= 1 && jet_pt_sum > 200;
    eventSelection[71] = jets_n >= 3 && cJets_n == 1 && tauJets_n == 1;
    eventSelection[72] = jets_n >= 3 && cJets_n == 1 && tauJets_n == 2;
    eventSelection[73] = jets_n >= 3 && cJets_n >= 1 && tauJets_n == 1;
    eventSelection[74] = jets_n >= 3 && cJets_n >= 1 && tauJets_n == 2;
    eventSelection[75] = jets_n >= 3 && cJets_n >= 1 && tauJets_n >= 1;
    eventSelection[76] = jets_n >= 3 && cJets_n >= 1 && tauJets_n >= 2;
    eventSelection[77] = jets_n >= 3 && bJets_n == 1 && MET_pt > 80;
    eventSelection[78] = jets_n >= 3 && cJets_n == 1 && MET_pt > 80;
    eventSelection[79] = jets_n >= 3 && bJets_n == 1 && cJets_n == 1 && MET_pt > 80;
    eventSelection[80] = jets_n >= 3 && bJets_n == 1 && MET_pt > 80 && jet_pt_sum > 200;
    eventSelection[81] = jets_n >= 3 && cJets_n == 1 && MET_pt > 80 && jet_pt_sum > 200;
    eventSelection[82] = jets_n >= 3 && bJets_n == 1 && cJets_n == 1 && MET_pt > 80 && jet_pt_sum > 200;



    TLorentzVector lepton[nGoodMuon + nGoodElectron];
    for(int i = 0; i < nGoodMuon; i++){
      int idx = GoodMuIdx.at(i);
      lepton[i].SetPtEtaPhiM(Muon_pt[idx], Muon_eta[idx], Muon_phi[idx], Muon_m[idx]);
    }
    for(int i = nGoodMuon; i < nGoodMuon + nGoodElectron; i++){
      int idx = GoodElecIdx.at(i-nGoodMuon);
      lepton[i].SetPtEtaPhiM(Electron_pt[idx], Electron_eta[idx], Electron_phi[idx], Electron_m[idx]);
    }

    TLorentzVector met;
    double MET_x = TMath::Abs(MET_pt)*TMath::Cos(MET_phi);
    double MET_y = TMath::Abs(MET_pt)*TMath::Sin(MET_phi);
    met.SetPxPyPzE( MET_x, MET_y, 0, MET_pt );

    TLorentzVector tauJet[2];
    for( unsigned int i = 0; i < tauIdx.size(); i++){
      int j = tauIdx.at(i);
      if( tauJet[0].Pt() == 0 ) tauJet[0].SetPtEtaPhiM(Jet_pt[j],Jet_eta[j],Jet_phi[j], Jet_m[j]);
      else if( tauJet[1].Pt() == 0) tauJet[1].SetPtEtaPhiM(Jet_pt[j],Jet_eta[j],Jet_phi[j], Jet_m[j]);
      else continue;
    }
    TLorentzVector bJet[2];
    for( unsigned int i = 0; i < bIdx.size(); i++){
      int j = bIdx.at(i);
      if( bJet[0].Pt() == 0 ) bJet[0].SetPtEtaPhiM(Jet_pt[j],Jet_eta[j],Jet_phi[j], Jet_m[j]);
      else if( bJet[1].Pt() == 0) bJet[1].SetPtEtaPhiM(Jet_pt[j],Jet_eta[j],Jet_phi[j], Jet_m[j]);
      else continue;
    }

    for(int MODE = 0; MODE < nCh; MODE++){
      for(int cut = 0; cut < nCuts; cut++) {
        if( eventSelection[cut] && pass_lep[MODE] ) {
          h_jets_n[MODE][cut]->Fill(jets_n);
          h_bJets_n[MODE][cut]->Fill(bJets_n);
          h_cJets_n[MODE][cut]->Fill(cJets_n);
          h_tauJets_n[MODE][cut]->Fill(tauJets_n);
          h_MET[MODE][cut]->Fill(MET_pt);
          h_jetptsum[MODE][cut]->Fill(jet_pt_sum);

          for(vector<int>::iterator iter=matchedIdx.begin(); iter!=matchedIdx.end(); ++iter){
            b_matched_tauTag.push_back(Jet_tauTag[*iter]);
            h_tauTag_matched[MODE][cut]->Fill(Jet_tauTag[*iter]);
          }
          b_matched_tauTag.clear();

          if( nGoodMuon + nGoodElectron > 1  ) h_lepDR[MODE][cut]->Fill(lepton[0].DeltaR(lepton[1]));
          if( tauIdx.size() > 0 && nGoodMuon > 0 ){
            float leptaudR = lepton[0].DeltaR(tauJet[0]);
            if( leptaudR < 2.2 ) continue;
            h_leptau_DR[MODE][cut]->Fill(leptaudR);
          }
          if( tauIdx.size() > 1 ) h_tautauDR[MODE][cut]->Fill(tauJet[1].DeltaR(tauJet[0]));

          if(tauIdx.size() > 1 && cIdx.size() > 0){
            TLorentzVector tempC;
            int k = cIdx.at(0);
            tempC.SetPtEtaPhiM(Jet_pt[k],Jet_eta[k],Jet_phi[k], Jet_m[k]);
            float topmass = (tempC+tauJet[0]+tauJet[1]).M();
            //if(topmass > 200 || topmass < 150) continue;
            h_sigTop[MODE][cut]->Fill(topmass);
          }
        }
      }//selection loop
    }
    tree->Fill();
  }
  tree->Write();

  for(int ch=0; ch < nCh; ch++){
    for(int i = 0; i < nCuts-1; i++){
      if( ch == 0 || ch == 1 || ch == 3 || ch == 4 || ch == 5 || ch == 7 || ch == 8) continue; 
      h_jets_n[ch][i]->Write();
      //h_bJets_n[ch][i]->Write();
      //h_cJets_n[ch][i]->Write();
      //h_tauJets_n[ch][i]->Write();
      //h_tauTag_matched[ch][i]->Write();
      //h_lepDR[ch][i]->Write();
      h_leptau_DR[ch][i]->Write();
      //h_tautauDR[ch][i]->Write();
      //h_MET[ch][i]->Write();
      //h_jetptsum[ch][i]->Write();
      //h_smTop[ch][i]->Write();
      h_sigTop[ch][i]->Write();
    }
  }
  EventInfo->Write();
  f->Close();
}
