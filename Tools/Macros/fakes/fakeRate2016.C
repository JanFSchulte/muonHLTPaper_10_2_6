#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TTree.h"
#include "TBranch.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include "DataFormats/Math/interface/deltaR.h"
#include "MuonHLTNtuples/Analyzers/src/MuonTree.h"
#include "TLorentzVector.h"

double muonmass = 0.10565837;
bool debug = false;

enum Sig { 
  Prompt = 0,
  DiMuon,
  LowPt,
  DisplacedOld,
  DisplacedNew,
};

bool selectTagMuon  (MuonCand, TH1F* );
bool selectProbeMuon(MuonCand, MuonCand, TH1F* );
bool selectMuon     (MuonCand);
bool selectGenMuon  (GenParticleCand);
bool matchMuon      (MuonCand, std::vector<HLTObjCand>, std::string);
bool firedL1        (          std::vector<HLTObjCand>, std::string);
//bool matchMuonWithL3(MuonCand, std::vector<HLTMuonCand>);
bool  matchMuonWithL3 (MuonCand, std::vector<HltTrackCand>);
bool  matchTrackWithGen(HltTrackCand, std::vector<GenParticleCand>);
bool  matchMuonWithGen (HLTMuonCand, std::vector<GenParticleCand>);

std::string getProbeFilter(int,bool);
float getLeadingPtCut(int);
float getTrailingPtCut(int);

void printProgBar(int);

double pt_bins[20]  = { 5, 7, 9, 12, 16,  20 ,  24 ,  27 ,   30,   35,   40,   45,   50,  60, 70 ,  90, 150, 250,500,1000 };
double dz_bins[11]  = {-15, -8, -6, -4, -2, 0, 2, 4, 6, 8, 15};
double eta_bins[16] = {-2.4, -2.1, -1.6, -1.2, -1.04, -0.9, -0.3, -0.2,  0.2, 0.3, 0.9, 1.04, 1.2, 1.6, 2.1, 2.4};
double iso_bins[12] = { 0  , 0.02, 0.04, 0.06, 0.08,  0.1, 0.12, 0.16, 0.2, 0.3, 0.6, 1   };
double offlineIsoCut = 0.15;

double weights[120] = {
0,
0,
0,
0,
0,
0,
3962.58,
1922.78,
770.845,
499.675,
210.538,
117.538,
63.4768,
37.1985,
23.0937,
14.7195,
10.3597,
7.20071,
5.35578,
4.11288,
3.27684,
2.64972,
2.23488,
1.9275,
1.69965,
1.51698,
1.36991,
1.26135,
1.16694,
1.08709,
1.01763,
0.938948,
0.877911,
0.811432,
0.747691,
0.685075,
0.627142,
0.570469,
0.516718,
0.4713,
0.427329,
0.387521,
0.350465,
0.319454,
0.290395,
0.266889,
0.248097,
0.233828,
0.219835,
0.20588,
0.197742,
0.191238,
0.183904,
0.175513,
0.175626,
0.177029,
0.174181,
0.172184,
0.175612,
0.182704,
0.182884,
0.178116,
0.188324,
0.198227,
0.192363,
0.197881,
0.212715,
0.212715,
0.217636,
0.214742,
0.239151,
0.228956,
0.216051,
0.266567,
0.296731,
0.281813,
0.265021,
0.270535,
0.315188,
0.292378,
0.260119,
0.330573,
0.355631,
0.335217,
0.3565,
0.418806,
0.380833,
0.296492,
0.607797,
0.40312,
0.437765,
0.433877,
0.672169,
0.320663,
0.655824,
0.669958,
0.893277,
0.409419,
0.527845,
0.691569,
0.571697,
1.33992,
0.638055,
0.466057,
0.687136,
1.19104,
0.198506,
1.19104,
2.97759,
0.297759,
1.78655,
0.595518,
1.78655,
0,
0,
0,
0,
0,
0,
0};
/// TAG-DEFINITION: 
std::string isofilterTag  = "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07::HLT";

/// for PROMPT-MUONS   (close-by and far-away) 
//std::string L1filter      = "hltL1fL1sMu22or25L1Filtered0::TEST"; 
//std::string L2filter      = "hltL2fL1sMu22or25L1f0L2Filtered10Q::TEST";
//std::string L3filter      = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q::TEST"; 
std::string L1filter      = "hltL1fL1sMu22or25L1Filtered0::MYHLT";
std::string L2filter      = "hltL2fL1sMu22or25L1f0L2Filtered10Q::MYHLT";
std::string L3filter      = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q::MYHLT"; 


// ******************************************
//       T&P definitions                    *
//                                          *
std::string thepassfilter  = L3filter;
//std::string theprobefilter = L1filter; 
float offlinePtCut         = 24.;
//                                          *
//                                          *
// ******************************************

void fakeRate2016(TString inputfilename="/eos/uscms/store/user/bmahakud/ProductionHLTAN_LPC_IterL3HighStat/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ProductionHLTAN_LPC_IterL3HighStat/181130_193653/0000/muonNtupleIterL3.root", std::string effmeasured="MC2018",bool isMC=false){

  int flavor=Sig::Prompt;

  bool doingL1 = thepassfilter.find("L1fL1") != std::string::npos; 

  TFile* outfile = TFile::Open(Form("%s_IterL3preFilter.root", effmeasured.c_str()),"RECREATE");
  std::cout << "output file: " << outfile -> GetName() << std::endl;

  TH1F* dimuon_mass             = new TH1F("h_dimuon_mass"          ,"dimuon_mass"      , 1500,  0,  150 );
  TH1F* tagiso                  = new TH1F("h_tagiso"               ,"tagiso"           ,  100,  0,  1   );  
  TH1F* ProbePt                 = new TH1F("h_ProbePt"              ,"ProbeMuonPt"      ,10000,  0,  1000   );  

  TEfficiency* fakeL2Pt           = new TEfficiency("fakeL2Pt"          ,"fakeL2Pt"           ,   19,  pt_bins );
  TEfficiency* fakeL2Eta          = new TEfficiency("fakeL2Eta"         ,"fakeL2Eta"          ,   15, eta_bins );
  TEfficiency* fakeL2Phi          = new TEfficiency("fakeL2Phi"         ,"fakeL2Phi"          ,   20, -3.2, 3.2);
  TEfficiency* fakeL2nVtx         = new TEfficiency("fakeL2nVtx"        ,"fakeL2nVtx"         ,   50, 0, 100);

  TEfficiency* fakeL2EtaPt          = new TEfficiency("fakeL2EtaPt"         ,"fakeL2EtaPt"          ,   15, eta_bins );
  TEfficiency* fakeL2PhiPt          = new TEfficiency("fakeL2PhiPt"         ,"fakeL2PhiPt"          ,   20, -3.2, 3.2);
  TEfficiency* fakeL2nVtxPt         = new TEfficiency("fakeL2nVtxPt"        ,"fakeL2nVtxPt"         ,   50, 0, 100);

  TEfficiency* fakeOIPt           = new TEfficiency("fakeOIPt"          ,"fakeOIPt"           ,   19,  pt_bins );
  TEfficiency* fakeOIEta          = new TEfficiency("fakeOIEta"         ,"fakeOIEta"          ,   15, eta_bins );
  TEfficiency* fakeOIPhi          = new TEfficiency("fakeOIPhi"         ,"fakeOIPhi"          ,   20, -3.2, 3.2);
  TEfficiency* fakeOInVtx         = new TEfficiency("fakeOInVtx"        ,"fakeOInVtx"         ,   50, 0, 100);

  TEfficiency* fakeOIEtaPt          = new TEfficiency("fakeOIEtaPt"         ,"fakeOIEtaPt"          ,   15, eta_bins );
  TEfficiency* fakeOIPhiPt          = new TEfficiency("fakeOIPhiPt"         ,"fakeOIPhiPt"          ,   20, -3.2, 3.2);
  TEfficiency* fakeOInVtxPt         = new TEfficiency("fakeOInVtxPt"        ,"fakeOInVtxPt"         ,   50, 0, 100);

  TEfficiency* fakeIOPt           = new TEfficiency("fakeIOPt"          ,"fakeIOPt"           ,   19,  pt_bins );
  TEfficiency* fakeIOEta          = new TEfficiency("fakeIOEta"         ,"fakeIOEta"          ,   15, eta_bins );
  TEfficiency* fakeIOPhi          = new TEfficiency("fakeIOPhi"         ,"fakeIOPhi"          ,   20, -3.2, 3.2);
  TEfficiency* fakeIOnVtx         = new TEfficiency("fakeIOnVtx"        ,"fakeIOnVtx"         ,   50, 0, 100);

  TEfficiency* fakeIOEtaPt          = new TEfficiency("fakeIOEtaPt"         ,"fakeIOEtaPt"          ,   15, eta_bins );
  TEfficiency* fakeIOPhiPt          = new TEfficiency("fakeIOPhiPt"         ,"fakeIOPhiPt"          ,   20, -3.2, 3.2);
  TEfficiency* fakeIOnVtxPt         = new TEfficiency("fakeIOnVtxPt"        ,"fakeIOnVtxPt"         ,   50, 0, 100);





  TEfficiency* fakePt           = new TEfficiency("fakePt"          ,"fakePt"           ,   19,  pt_bins );
  TEfficiency* fakeEta          = new TEfficiency("fakeEta"         ,"fakeEta"          ,   15, eta_bins );
  TEfficiency* fakePhi          = new TEfficiency("fakePhi"         ,"fakePhi"          ,   20, -3.2, 3.2);
  TEfficiency* fakenVtx         = new TEfficiency("fakenVtx"        ,"fakenVtx"         ,   50, 0, 100);

  TEfficiency* fakeEtaPt          = new TEfficiency("fakeEtaPt"         ,"fakeEtaPt"          ,   15, eta_bins );
  TEfficiency* fakePhiPt          = new TEfficiency("fakePhiPt"         ,"fakePhiPt"          ,   20, -3.2, 3.2);
  TEfficiency* fakenVtxPt         = new TEfficiency("fakenVtxPt"        ,"fakenVtxPt"         ,   50, 0, 100);

  TEfficiency* fakeChi2         = new TEfficiency("fakeChi2"        , "fakezhi2"        ,  60,    0.,  7.);
  TEfficiency* fakeDxy          = new TEfficiency("fakeDxy"         , "fakeDxy"         ,  100,  -0.3, 0.3);
  TEfficiency* fakeDz           = new TEfficiency("fakeDz"          , "fakeDz"          ,   10,   dz_bins);
  TEfficiency* fakePixHit       = new TEfficiency("fakePixHit"      , "fakePixHit"      ,   20,  -0.5,19.5);
  TEfficiency* fakeLayHit       = new TEfficiency("fakeLayHit"      , "fakeLayHit"      ,   16,   2.5,18.5);
  TEfficiency* fakePixLay       = new TEfficiency("fakePixLay"      , "fakePixLay"      ,    7,  -0.5, 6.5);

  TEfficiency* fakeChi2Pt         = new TEfficiency("fakeChi2Pt"        , "fakezhi2Pt"        ,  60,    0.,  7.);
  TEfficiency* fakeDxyPt          = new TEfficiency("fakeDxyPt"         , "fakeDxyPt"         ,  100,  -0.3, 0.3);
  TEfficiency* fakeDzPt           = new TEfficiency("fakeDzPt"          , "fakeDzPt"          ,   10,   dz_bins);
  TEfficiency* fakePixHitPt       = new TEfficiency("fakePixHitPt"      , "fakePixHitPt"      ,   20,  -0.5,19.5);
  TEfficiency* fakeLayHitPt       = new TEfficiency("fakeLayHitPt"      , "fakeLayHitPt"      ,   16,   2.5,18.5);
  TEfficiency* fakePixLayPt       = new TEfficiency("fakePixLayPt"      , "fakePixLayPt"      ,    7,  -0.5, 6.5);



  double offlineiso04 = 100;
  
  TChain *tree = new TChain("muonNtuples/muonTree");
  tree->Add(inputfilename); 

 
  
  if (!tree) {
    std::cout << " *** tree not found *** " << std::endl;
    return;
  }
    
  MuonEvent* ev      = new MuonEvent(); 
  //TBranch*  evBranch = tree->GetBranch("event"); 
  //evBranch -> SetAddress(&ev);
  TBranch*  evBranch; 
  tree-> SetBranchAddress("event",&ev,&evBranch);


  int nentries = tree->GetEntries();
//  int nentries = 10000000;// tree->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;

  bool flagfile = false;

  for (Int_t eventNo=0; eventNo < nentries; eventNo++)     {
    Int_t IgetEvent   = tree   -> GetEvent(eventNo);
    printProgBar((int)(eventNo*100./nentries));
 

    double weight = 1;
    //if (isMC) weight = weights[ev -> nVtx];
    unsigned int nmuons = ev->hltmuons.size(); 
    for (int imu = 0; imu < nmuons; imu++){ 
        bool isFake = true;
	// select the probe muon  & match to the probe:  
	if (matchMuonWithGen(ev->hltmuons.at(imu),ev->genParticles)) isFake = false;
	fakePt   -> FillWeighted( isFake, weight, ev -> hltmuons.at(imu).pt );
	//if (ev -> muons.at(jmu).pt < offlinePtCut) continue;
	
	 fakeEta  -> FillWeighted( isFake, weight, ev -> hltmuons.at(imu).eta);
	 fakePhi  -> FillWeighted( isFake, weight, ev -> hltmuons.at(imu).phi);
	 fakenVtx -> FillWeighted( isFake, weight, ev -> nVtx                );
		
/*	 fakeChi2   -> FillWeighted( pass, weight, ev -> muons.at(imu).innerchi2 );
	 fakeDxy    -> FillWeighted( pass, weight, ev -> muons.at(imu).innerdxy);
	 fakeDz     -> FillWeighted( pass, weight, ev -> muons.at(imu).innerdz);
	 fakePixHit -> FillWeighted( pass, weight, ev -> muons.at(imu).innerpixelHits);
	 fakeLayHit -> FillWeighted( pass, weight, ev -> muons.at(imu).innerlayerHits);
	 fakePixLay -> FillWeighted( pass, weight, ev -> muons.at(imu).innerpixelLayers);
*/
	if (ev -> hltmuons.at(imu).pt < offlinePtCut) continue;

	 fakeEtaPt  -> FillWeighted( isFake, weight, ev -> hltmuons.at(imu).eta);
	 fakePhiPt  -> FillWeighted( isFake, weight, ev -> hltmuons.at(imu).phi);
	 fakenVtxPt -> FillWeighted( isFake, weight, ev -> nVtx                );
/*
	 fakeChi2Pt -> FillWeighted( pass, weight, ev -> muons.at(imu).innerchi2 );
	 fakeDxyPt  -> FillWeighted( pass, weight, ev -> muons.at(imu).innerdxy);
	 fakeDzPt   -> FillWeighted( pass, weight, ev -> muons.at(imu).innerdz);
	 fakePixHitPt-> FillWeighted( pass, weight, ev -> muons.at(imu).innerpixelHits);
	 fakeLayHitPt-> FillWeighted( pass, weight, ev -> muons.at(imu).innerlayerHits);
	 fakePixLayPt-> FillWeighted( pass, weight, ev -> muons.at(imu).innerpixelLayers);
*/

    }


    nmuons = ev->hltTrackOI.size(); 
    for (int imu = 0; imu < nmuons; imu++){ 
        bool isFake = true;
	// select the probe muon  & match to the probe:  
	if (matchTrackWithGen(ev->hltTrackOI.at(imu),ev->genParticles)) isFake = false;
	 fakeOIPt   -> FillWeighted( isFake, weight, ev -> hltTrackOI.at(imu).pt );
	 fakeOIEta  -> FillWeighted( isFake, weight, ev -> hltTrackOI.at(imu).eta);
	 fakeOIPhi  -> FillWeighted( isFake, weight, ev -> hltTrackOI.at(imu).phi);
	 fakeOInVtx -> FillWeighted( isFake, weight, ev -> nVtx                );
		
	if (ev -> hltTrackOI.at(imu).pt < offlinePtCut) continue;

	 fakeOIEtaPt  -> FillWeighted( isFake, weight, ev -> hltTrackOI.at(imu).eta);
	 fakeOIPhiPt  -> FillWeighted( isFake, weight, ev -> hltTrackOI.at(imu).phi);
	 fakeOInVtxPt -> FillWeighted( isFake, weight, ev -> nVtx                );


    }

    nmuons = ev->hltTrackIOL2.size(); 
    for (int imu = 0; imu < nmuons; imu++){ 
        bool isFake = true;
	// select the probe muon  & match to the probe:  
	if (matchTrackWithGen(ev->hltTrackIOL2.at(imu),ev->genParticles)) isFake = false;
	 fakeIOPt   -> FillWeighted( isFake, weight, ev -> hltTrackIOL2.at(imu).pt );
	 fakeIOEta  -> FillWeighted( isFake, weight, ev -> hltTrackIOL2.at(imu).eta);
	 fakeIOPhi  -> FillWeighted( isFake, weight, ev -> hltTrackIOL2.at(imu).phi);
	 fakeIOnVtx -> FillWeighted( isFake, weight, ev -> nVtx                );
		
	if (ev -> hltTrackIOL2.at(imu).pt < offlinePtCut) continue;

	 fakeIOEtaPt  -> FillWeighted( isFake, weight, ev -> hltTrackIOL2.at(imu).eta);
	 fakeIOPhiPt  -> FillWeighted( isFake, weight, ev -> hltTrackIOL2.at(imu).phi);
	 fakeIOnVtxPt -> FillWeighted( isFake, weight, ev -> nVtx                );


    } 
    nmuons = ev->L2muons.size(); 
    for (int imu = 0; imu < nmuons; imu++){ 
        bool isFake = true;
	// select the probe muon  & match to the probe:  
	if (matchMuonWithGen(ev->L2muons.at(imu),ev->genParticles)) isFake = false;
	 fakeL2Pt   -> FillWeighted( isFake, weight, ev -> L2muons.at(imu).pt );
	 fakeL2Eta  -> FillWeighted( isFake, weight, ev -> L2muons.at(imu).eta);
	 fakeL2Phi  -> FillWeighted( isFake, weight, ev -> L2muons.at(imu).phi);
	 fakeL2nVtx -> FillWeighted( isFake, weight, ev -> nVtx                );
		
	if (ev -> L2muons.at(imu).pt < offlinePtCut) continue;

	 fakeL2EtaPt  -> FillWeighted( isFake, weight, ev -> L2muons.at(imu).eta);
	 fakeL2PhiPt  -> FillWeighted( isFake, weight, ev -> L2muons.at(imu).phi);
	 fakeL2nVtxPt -> FillWeighted( isFake, weight, ev -> nVtx                );


    }
  }
  fakePt->Write();
  fakeEta->Write();
  fakePhi->Write();
  fakenVtx->Write();
  fakeEtaPt->Write();
  fakePhiPt->Write();
  fakenVtxPt->Write();


  fakeOIPt->Write();
  fakeOIEta->Write();
  fakeOIPhi->Write();
  fakeOInVtx->Write();
  fakeOIEtaPt->Write();
  fakeOIPhiPt->Write();
  fakeOInVtxPt->Write();

  fakeIOPt->Write();
  fakeIOEta->Write();
  fakeIOPhi->Write();
  fakeIOnVtx->Write();
  fakeIOEtaPt->Write();
  fakeIOPhiPt->Write();
  fakeIOnVtxPt->Write();




  fakeL2Pt->Write();
  fakeL2Eta->Write();
  fakeL2Phi->Write();
  fakeL2nVtx->Write();
  fakeL2EtaPt->Write();
  fakeL2PhiPt->Write();
  fakeL2nVtxPt->Write();


  fakeChi2->Write();
  fakeDxy->Write();
  fakeDz->Write();
  fakePixHit->Write();
  fakeLayHit->Write();
  fakePixLay->Write();

  fakeChi2Pt->Write();
  fakeDxyPt->Write();
  fakeDzPt->Write();
  fakePixHitPt->Write();
  fakeLayHitPt->Write();
  fakePixLayPt->Write();


  outfile          -> Close();  
 
  return;
}
bool firedL1( std::vector<HLTObjCand> toc, std::string L1FilterName){ 
  int ntoc = toc.size();
  for ( std::vector<HLTObjCand>::const_iterator it = toc.begin(); it != toc.end(); ++it ) { 
    if ( it->filterTag.compare(L1FilterName) == 0) return true;
  }
  return false;
}
bool matchMuon(MuonCand mu, std::vector<HLTObjCand> toc, std::string tagFilterName){

  bool match = false;
  int ntoc = toc.size();

  float minDR = 0.1; 
  if (tagFilterName.find("L1fL1") != std::string::npos) minDR = 0.3;
  float theDR = 100;
  for ( std::vector<HLTObjCand>::const_iterator it = toc.begin(); it != toc.end(); ++it ) {
    if ( it->filterTag.compare(tagFilterName) == 0) {
      theDR = deltaR(it -> eta, it -> phi, mu.eta, mu.phi);
      if (theDR < minDR){
        minDR = theDR;
        match = true;
      }
    }
  }
  
  return match;
}

bool selectTagMuon(MuonCand mu, TH1F* tagh){
  
  if (!( mu.pt         > 26.)) return false; 
  if (!( fabs(mu.eta)  < 2.4 )) return false; 
  if (!( mu.isTight    == 1  )) return false; 
  
  //add isolation cut
  float offlineiso04 = mu.chargedDep_dR04 + std::max(0., mu.photonDep_dR04 + mu.neutralDep_dR04 - 0.5*mu.puPt_dR04);
  offlineiso04       = offlineiso04 / mu.pt;
  tagh -> Fill(offlineiso04);
  if (offlineiso04   > offlineIsoCut) return false; 

  return true;
}

float getLeadingPtCut(int signature){ 
  float ptcut = 0.;
  if (signature == Sig::Prompt) ptcut = 26.;
  if (signature == Sig::DiMuon) ptcut = 18.;
  if (signature == Sig::LowPt ) ptcut = 0.;
  return ptcut;
}

float getTrailingPtCut(int signature){ 
  float ptcut = 0.;
  if (signature == Sig::Prompt) ptcut = 24.;
  if (signature == Sig::DiMuon) ptcut = 8. ;
  if (signature == Sig::LowPt ) ptcut = 0. ;
  return ptcut;
}


bool selectMuon(MuonCand mu){  
  if (!( mu.pt         > offlinePtCut  )) return false; 
  if (!( fabs(mu.eta)  < 2.4 )) return false;
  if (!( mu.isLoose    == 1  )) return false; 
  return true;
}

bool selectGenMuon(GenParticleCand mu){
  if (!( fabs(mu.pdgId) == 13)) return false;
  if (!( mu.pt         > offlinePtCut  )) return false; 
  if (!( fabs(mu.eta)  < 2.4 )) return false;
  return true;
}

//select the probe muon
bool selectProbeMuon(MuonCand mu, MuonCand tagMu, TH1F* dimuon_mass){
  
  if (mu.pt == tagMu.pt  && 
      mu.eta == tagMu.eta &&
      mu.phi == tagMu.phi ) 
    return false;
  if ( deltaR(tagMu.eta,tagMu.phi, mu.eta, mu.phi) <= 0.3 ) return false;
  
  if (!( mu.pt          > 0  )) return false; 
  if (!( fabs(mu.eta)  < 2.4 )) return false; 
  if (!( mu.isTight    == 1  )) return false; 
  if (mu.charge * tagMu.charge > 0) return false;
  //add isolation cut
  float offlineiso04 = mu.chargedDep_dR04 + std::max(0., mu.photonDep_dR04 + mu.neutralDep_dR04 - 0.5*mu.puPt_dR04);
  offlineiso04       = offlineiso04 / mu.pt;
  if (offlineiso04   > offlineIsoCut) return false; 
  
  TLorentzVector mu1, mu2;
  mu1.SetPtEtaPhiM (mu.pt   , mu.eta   , mu.phi   , muonmass);
  mu2.SetPtEtaPhiM (tagMu.pt, tagMu.eta, tagMu.phi, muonmass);
  double mumumass = (mu1 + mu2).M();
  dimuon_mass -> Fill(mumumass); 
  if (! (mumumass > 81. && mumumass < 101. )) return false;
  
  return true;
}

bool matchMuonWithL3(MuonCand mu, std::vector<HltTrackCand> L3cands){

  bool match = false;
  float minDR = 0.3;
  float theDR = 100;
  for ( std::vector<HltTrackCand>::const_iterator it = L3cands.begin(); it != L3cands.end(); ++it ) {
    theDR = deltaR(it -> eta, it -> phi, mu.eta, mu.phi);
    if (theDR < minDR){
      minDR = theDR;
      match = true;
    }
  }
  return match;
}
bool matchMuonWithGen(HLTMuonCand mu, std::vector<GenParticleCand> genParticles){

  bool match = false;
  float minDR = 0.3;
  float theDR = 100;
  for ( std::vector<GenParticleCand>::const_iterator it = genParticles.begin(); it != genParticles.end(); ++it ) {
    //if (!(fabs(pdgId->pdgId) == 13)) continue;
    theDR = deltaR(it -> eta, it -> phi, mu.eta, mu.phi);
    if (theDR < minDR){
      minDR = theDR;
      match = true;
    }
  }
  return match;
}
bool matchTrackWithGen(HltTrackCand mu, std::vector<GenParticleCand> genParticles){

  bool match = false;
  float minDR = 0.3;
  float theDR = 100;
  for ( std::vector<GenParticleCand>::const_iterator it = genParticles.begin(); it != genParticles.end(); ++it ) {
    //if (!(fabs(pdgId->pdgId) == 13)) continue;
    theDR = deltaR(it -> eta, it -> phi, mu.eta, mu.phi);
    if (theDR < minDR){
      minDR = theDR;
      match = true;
    }
  }
  return match;
}


std::string getProbeFilter(int signature,bool isMC){
  if (signature == Sig::Prompt) { 
    if (isMC) return "hltL1fL1sMu22or25L1Filtered0::MYHLT"; //Prompt
    else return "hltL1fL1sMu22or25L1Filtered0::TEST"; //Prompt
  }
  if (signature == Sig::DiMuon) { 
    return "hltL1fL1sDoubleMu155L1Filtered0::TEST"; //Dimuon
  }
  if (signature == Sig::LowPt ) {
    return "hltL1fL1sL1sDoubleMu4SQOSdRMax1p2L1Filtered0::TEST";  //JPsi
  }
  if (signature == Sig::DisplacedOld ) { 
    return "hltL1fDimuonL1Filtered0::TEST"; //Displaced OLD
  }
  if (signature == Sig::DisplacedNew ) {
    return "hltDimuon3L1Filtered0::TEST"; //Displaced NEW
  }
  return "none";
}
void printProgBar( int percent ){
  std::string bar;  
  for(int i = 0; i < 50; i++){
    if( i < (percent/2)){
      bar.replace(i,1,"=");
    }else if( i == (percent/2)){
      bar.replace(i,1,">");
    }else{
      bar.replace(i,1," ");
    }
  }

  std::cout<< "\r" "[" << bar << "] ";
  std::cout.width( 3 );
  std::cout<< percent << "%     " << std::flush;
}



