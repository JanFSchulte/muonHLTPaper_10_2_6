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

bool selectTagMuon  (MuonCand);
bool selectProbeMuon(MuonCand, MuonCand );
bool selectMuon     (MuonCand);
bool selectGenMuon  (GenParticleCand);
bool matchMuon      (MuonCand, std::vector<HLTObjCand>, std::string);
bool firedL1        (          std::vector<HLTObjCand>, std::string);
//bool matchMuonWithL3(MuonCand, std::vector<HLTMuonCand>);
bool  matchMuonWithL3 (MuonCand, std::vector<HltTrackCand>);

std::string getProbeFilter(int);
float getLeadingPtCut(int);
float getTrailingPtCut(int);

void printProgBar(int);

double pt_bins[20]  = { 5, 7, 9, 12, 16,  20 ,  24 ,  27 ,   30,   35,   40,   45,   50,  60, 70 ,  90, 150, 250,500,1000 };
double dz_bins[11]  = {-15, -8, -6, -4, -2, 0, 2, 4, 6, 8, 15};
double eta_bins[16] = {-2.4, -2.1, -1.6, -1.2, -1.04, -0.9, -0.3, -0.2,  0.2, 0.3, 0.9, 1.04, 1.2, 1.6, 2.1, 2.4};
double iso_bins[12] = { 0  , 0.02, 0.04, 0.06, 0.08,  0.1, 0.12, 0.16, 0.2, 0.3, 0.6, 1   };
double offlineIsoCut = 0.15;


/// TAG-DEFINITION: 
std::string isofilterTag  = "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07::HLT";

/// for PROMPT-MUONS   (close-by and far-away) 
//std::string L1filter      = "hltL1fL1sMu22or25L1Filtered0::TEST"; 
//std::string L2filter      = "hltL2fL1sMu22or25L1f0L2Filtered10Q::TEST";
//std::string L3filter      = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q::TEST"; 
std::string L1filter      = "hltL1fL1sMu22or25L1Filtered0::TEST";
std::string L2filter      = "hltL2fL1sMu22or25L1f0L2Filtered10Q::TEST";
std::string L3filter      = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q::TEST"; 


// ******************************************
//       T&P definitions                    *
//                                          *
std::string thepassfilter  = L3filter;
//std::string theprobefilter = L1filter; 
float offlinePtCut         = 26.;
//                                          *
//                                          *
// ******************************************

void getVtxDist(TString inputfilename="/eos/uscms/store/user/bmahakud/ProductionHLTAN_LPC_IterL3HighStat/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ProductionHLTAN_LPC_IterL3HighStat/181130_193653/0000/muonNtupleIterL3.root", std::string effmeasured="MC2018"){

  int flavor=Sig::Prompt;

  bool doingL1 = thepassfilter.find("L1fL1") != std::string::npos; 

  TFile* outfile = TFile::Open(Form("%s_IterL3preFilter.root", effmeasured.c_str()),"RECREATE");
  std::cout << "output file: " << outfile -> GetName() << std::endl;

  TH1F* nVtx=new TH1F("nVtx","nVtx",120,0,120 );



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
  offlinePtCut = getLeadingPtCut(flavor);
  float ptcut1 = getLeadingPtCut(flavor);
  float ptcut2 = getTrailingPtCut(flavor);
	
  for (Int_t eventNo=0; eventNo < nentries; eventNo++)     {
    Int_t IgetEvent   = tree   -> GetEvent(eventNo);
    printProgBar((int)(eventNo*100./nentries));
   	
	nVtx->Fill(ev->nVtx); 
  
  } 
 //Writing the histograms in a file.
  outfile           -> cd();


  nVtx->Write();


 
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

bool selectTagMuon(MuonCand mu){
  
  if (!( mu.pt         > offlinePtCut)) return false; 
  if (!( fabs(mu.eta)  < 2.4 )) return false; 
  if (!( mu.isTight    == 1  )) return false; 
  
  //add isolation cut
  float offlineiso04 = mu.chargedDep_dR04 + std::max(0., mu.photonDep_dR04 + mu.neutralDep_dR04 - 0.5*mu.puPt_dR04);
  offlineiso04       = offlineiso04 / mu.pt;
  if (offlineiso04   > offlineIsoCut) return false; 

  return true;
}

float getLeadingPtCut(int signature){ 
  float ptcut = 0.;
  if (signature == Sig::Prompt) ptcut = 29.;
  if (signature == Sig::DiMuon) ptcut = 18.;
  if (signature == Sig::LowPt ) ptcut = 0.;
  return ptcut;
}

float getTrailingPtCut(int signature){ 
  float ptcut = 0.;
  if (signature == Sig::Prompt) ptcut = 27.;
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
bool selectProbeMuon(MuonCand mu, MuonCand tagMu){
  
  if (mu.pt == tagMu.pt  && 
      mu.eta == tagMu.eta &&
      mu.phi == tagMu.phi ) 
    return false;
  
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
  if (! (mumumass > 81. && mumumass < 101. )) return false;
  
  return true;
}

bool matchMuonWithL3(MuonCand mu, std::vector<HltTrackCand> L3cands){

  bool match = false;
  float minDR = 0.1;
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




std::string getProbeFilter(int signature){
  if (signature == Sig::Prompt) { 
    return "hltL1fL1sMu22or25L1Filtered0::TEST"; //Prompt
    //return "hltL1fL1sMu22or25L1Filtered0::MYHLT"; //Prompt
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
