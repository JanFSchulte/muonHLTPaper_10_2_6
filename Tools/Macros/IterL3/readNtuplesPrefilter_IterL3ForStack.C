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

void readNtuplesPrefilter_IterL3ForStack(TString inputfilename="/eos/uscms/store/user/bmahakud/ProductionHLTAN_LPC_IterL3HighStat/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ProductionHLTAN_LPC_IterL3HighStat/181130_193653/0000/muonNtupleIterL3.root", std::string effmeasured="MC2018",bool isMC=false){

  int flavor=Sig::Prompt;

  bool doingL1 = thepassfilter.find("L1fL1") != std::string::npos; 

  TFile* outfile = TFile::Open(Form("%s_IterL3preFilter.root", effmeasured.c_str()),"RECREATE");
  std::cout << "output file: " << outfile -> GetName() << std::endl;

  TH1F* dimuon_mass             = new TH1F("h_dimuon_mass"          ,"dimuon_mass"      , 1500,  0,  150 );
  TH1F* tagiso                  = new TH1F("h_tagiso"               ,"tagiso"           ,  100,  0,  1   );  
  TH1F* ProbePt                 = new TH1F("h_ProbePt"              ,"ProbeMuonPt"      ,10000,  0,  1000   );  

  TH1F* OIeffPt=new TH1F("OIeffPt","OIeffPt",19,  pt_bins);
  TH1F* OIeffEta=new TH1F("OIeffEta","OIeffEta",15, eta_bins); 
  TH1F* OIeffPhi=new TH1F("OIeffPhi","OIeffPhi",20, -3.2, 3.2);
  TH1F* OIeffnVtx=new TH1F("OIeffnVtx","OIeffnVtx",50,0,100);

  TH1F* effPtTot=new TH1F("effPtTot","effPtTot",19,  pt_bins);
  TH1F* effEtaTot=new TH1F("effEtaTot","effEtaTot",15, eta_bins);
  TH1F* effPhiTot=new TH1F("effPhiTot","effPhiTot",20, -3.2, 3.2);
  TH1F* effnVtxTot=new TH1F("effnVtxTot","effnVtxTot",50,0,100);



  TH1F* OIplusIOL2effPt=new TH1F("OIplusIOL2effPt","OIplusIOL2effPt",19,  pt_bins);
  TH1F* OIplusIOL2effEta=new TH1F("OIplusIOL2effEta","OIplusIOL2effEta",15, eta_bins);
  TH1F* OIplusIOL2effPhi=new TH1F("OIplusIOL2effPhi","OIplusIOL2effPhi",20, -3.2, 3.2);
  TH1F* OIplusIOL2effnVtx=new TH1F("OIplusIOL2effnVtx","OIplusIOL2effnVtx",50,0,100);

  TH1F* OIplusIOL2plusIOL1effPt=new TH1F("OIplusIOL2plusIOL1effPt","OIplusIOL2plusIOL1effPt",19,  pt_bins);
  TH1F* OIplusIOL2plusIOL1effEta=new TH1F("OIplusIOL2plusIOL1effEta","OIplusIOL2plusIOL1effEta",15, eta_bins);
  TH1F* OIplusIOL2plusIOL1effPhi=new TH1F("OIplusIOL2plusIOL1effPhi","OIplusIOL2plusIOL1effPhi",20, -3.2, 3.2);
  TH1F* OIplusIOL2plusIOL1effnVtx=new TH1F("OIplusIOL2plusIOL1effnVtx","OIplusIOL2plusIOL1effnVtx",50,0,100); 

  OIeffPt->Sumw2(); 
  OIeffEta->Sumw2(); 
  OIeffPhi->Sumw2(); 
  OIeffnVtx->Sumw2(); 

  effPtTot->Sumw2(); 
  effEtaTot->Sumw2(); 
  effPhiTot->Sumw2(); 
  effnVtxTot->Sumw2(); 

  OIplusIOL2effPt->Sumw2();
  OIplusIOL2effEta->Sumw2();
  OIplusIOL2effPhi->Sumw2();
  OIplusIOL2effnVtx->Sumw2();


  OIplusIOL2plusIOL1effPt->Sumw2();
  OIplusIOL2plusIOL1effEta->Sumw2();
  OIplusIOL2plusIOL1effPhi->Sumw2();
  OIplusIOL2plusIOL1effnVtx->Sumw2();


  TEfficiency* OIPt    = new TEfficiency("OIPt"   ,"OIPt"    ,   19,  pt_bins );
  TEfficiency* OIIOPt    = new TEfficiency("OIIOPt"   ,"OIIOPt"    ,   19,  pt_bins );
  TEfficiency* OIIOIOPt    = new TEfficiency("OIIOIOPt"   ,"OIIOIOPt"    ,   19,  pt_bins );

  TEfficiency* OIEta          = new TEfficiency("OIEta"         ,"OIEta"          ,   15, eta_bins );
  TEfficiency* OIIOEta          = new TEfficiency("OIIOEta"         ,"OIIOEta"          ,   15, eta_bins );
  TEfficiency* OIIOIOEta          = new TEfficiency("OIIOIOEta"         ,"OIIOIOEta"          ,   15, eta_bins );

  TEfficiency* OIPhi          = new TEfficiency("OIPhi"         ,"OIPhi"          ,   20, -3.2, 3.2);
  TEfficiency* OIIOPhi          = new TEfficiency("OIIOPhi"         ,"OIIOPhi"          ,   20, -3.2, 3.2);
  TEfficiency* OIIOIOPhi          = new TEfficiency("OIIOIOPhi"         ,"OIIOIOPhi"          ,   20, -3.2, 3.2);
 
  TEfficiency* OInVtx          = new TEfficiency("OInVtx"         ,"OInVtx"          ,   50, 0, 100);
  TEfficiency* OIIOnVtx          = new TEfficiency("OIIOnVtx"         ,"OIIOnVtx"          ,   50, 0, 100);
  TEfficiency* OIIOIOnVtx          = new TEfficiency("OIIOIOnVtx"         ,"OIIOIOnVtx"          ,   50, 0, 100);

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
  std::string theprobefilter = getProbeFilter(flavor,isMC);
  offlinePtCut = getLeadingPtCut(flavor);
  float ptcut1 = getLeadingPtCut(flavor);
  float ptcut2 = getTrailingPtCut(flavor);
	
  for (Int_t eventNo=0; eventNo < nentries; eventNo++)     {
    Int_t IgetEvent   = tree   -> GetEvent(eventNo);
    printProgBar((int)(eventNo*100./nentries));
 

    double weight = 1;
    if (isMC) weight = weights[ev -> nVtx];
   
    unsigned int nmuons = ev->muons.size(); 
    if (nmuons < 2) continue; 
    for (int imu = 0; imu < nmuons; imu++){ 
      // select the tag muon        
      if (debug) cout <<"select Tag muon" << endl;
      if (! selectTagMuon(ev -> muons.at(imu), tagiso)) continue; 
      if (! matchMuon(ev -> muons.at(imu), ev -> hltTag.objects, isofilterTag)) continue;
      for (int jmu = 0; jmu < nmuons; jmu++){
	bool passOI   = false;
	bool passIOL2 = false;
        bool passIOL1 = false;
	// select the probe muon  & match to the probe:  
	if (!selectProbeMuon(ev -> muons.at(jmu), ev -> muons.at(imu), dimuon_mass)) continue;
	if (!(matchMuon(ev -> muons.at(jmu), ev -> hlt.objects, theprobefilter))) continue;
		// select the pass muon
		//
	ProbePt->Fill(ev -> muons.at(jmu).pt);
	if (matchMuonWithL3(ev->muons.at(jmu),ev->hltTrackOI)) passOI = true;
        if (matchMuonWithL3(ev->muons.at(jmu),ev->hltTrackIOL2)) passIOL2 = true;	
        if (matchMuonWithL3(ev->muons.at(jmu),ev->hltTrackIOL1)) passIOL1 = true;
	if (isMC){
		OIPt   -> FillWeighted( passOI, weight, ev -> muons.at(jmu).pt );
		OIIOPt   -> FillWeighted( (passOI || passIOL2), weight, ev -> muons.at(jmu).pt );
		OIIOIOPt   -> FillWeighted( (passOI || passIOL2 || passIOL1), weight, ev -> muons.at(jmu).pt );
	}
	else{
		OIPt   -> Fill( passOI, ev -> muons.at(jmu).pt );
		OIIOPt   -> Fill( (passOI || passIOL2), ev -> muons.at(jmu).pt );
		OIIOIOPt   -> Fill( (passOI || passIOL2 || passIOL1), ev -> muons.at(jmu).pt );
	}
         effPtTot->Fill(ev -> muons.at(jmu).pt);	

         if(passOI)OIeffPt->Fill(ev -> muons.at(jmu).pt);
	 if(passOI || passIOL2)OIplusIOL2effPt->Fill(ev -> muons.at(jmu).pt);
	 if(passOI || passIOL2 || passIOL1)OIplusIOL2plusIOL1effPt->Fill(ev -> muons.at(jmu).pt);

	 if (ev -> muons.at(jmu).pt < offlinePtCut) continue;

	if (isMC){
		OIEta  -> FillWeighted( passOI, weight, ev -> muons.at(jmu).eta);
		OIPhi  -> FillWeighted( passOI, weight, ev -> muons.at(jmu).phi);
		OInVtx -> FillWeighted( passOI, weight, ev -> nVtx             );

		OIIOEta  -> FillWeighted( (passOI || passIOL2), weight, ev -> muons.at(jmu).eta);
		OIIOPhi  -> FillWeighted( (passOI || passIOL2), weight, ev -> muons.at(jmu).phi);
		OIIOnVtx -> FillWeighted( (passOI || passIOL2), weight, ev -> nVtx             );

		OIIOIOEta  -> FillWeighted( (passOI || passIOL2 || passIOL1), weight, ev -> muons.at(jmu).eta);
		OIIOIOPhi  -> FillWeighted( (passOI || passIOL2 || passIOL1), weight, ev -> muons.at(jmu).phi);
		OIIOIOnVtx -> FillWeighted( (passOI || passIOL2 || passIOL1), weight, ev -> nVtx             );
	}
	else{
		OIEta  -> Fill( passOI, ev -> muons.at(jmu).eta);
		OIPhi  -> Fill( passOI, ev -> muons.at(jmu).phi);
		OInVtx -> Fill( passOI, ev -> nVtx             );

		OIIOEta  -> Fill( (passOI || passIOL2), ev -> muons.at(jmu).eta);
		OIIOPhi  -> Fill( (passOI || passIOL2), ev -> muons.at(jmu).phi);
		OIIOnVtx -> Fill( (passOI || passIOL2), ev -> nVtx             );

		OIIOIOEta  -> Fill( (passOI || passIOL2 || passIOL1), ev -> muons.at(jmu).eta);
		OIIOIOPhi  -> Fill( (passOI || passIOL2 || passIOL1), ev -> muons.at(jmu).phi);
		OIIOIOnVtx -> Fill( (passOI || passIOL2 || passIOL1), ev -> nVtx             );
	}


         effEtaTot->Fill(ev -> muons.at(jmu).eta);
         effPhiTot->Fill(ev -> muons.at(jmu).phi);
         effnVtxTot->Fill(ev ->nVtx);


         if(passOI)OIeffEta->Fill(ev -> muons.at(jmu).eta);
         if(passOI)OIeffPhi->Fill(ev -> muons.at(jmu).phi);
         if(passOI)OIeffnVtx->Fill(ev -> nVtx);

         if(passOI || passIOL2)OIplusIOL2effEta->Fill(ev -> muons.at(jmu).eta);
         if(passOI || passIOL2)OIplusIOL2effPhi->Fill(ev -> muons.at(jmu).phi);
         if(passOI || passIOL2)OIplusIOL2effnVtx->Fill(ev -> nVtx);  

         if(passOI || passIOL2 || passIOL1)OIplusIOL2plusIOL1effEta->Fill(ev -> muons.at(jmu).eta);
         if(passOI || passIOL2 || passIOL1)OIplusIOL2plusIOL1effPhi->Fill(ev -> muons.at(jmu).phi);
         if(passOI || passIOL2 || passIOL1)OIplusIOL2plusIOL1effnVtx->Fill(ev -> nVtx); 





      } // nmuons
    }
  
  } 
  OIeffPt->Divide(effPtTot);
  OIeffEta->Divide(effEtaTot);
  OIeffPhi->Divide(effPhiTot);
  OIeffnVtx->Divide(effnVtxTot);

  OIplusIOL2effPt->Divide(effPtTot);
  OIplusIOL2effEta->Divide(effEtaTot);
  OIplusIOL2effPhi->Divide(effPhiTot);
  OIplusIOL2effnVtx->Divide(effnVtxTot);
  
  OIplusIOL2plusIOL1effPt->Divide(effPtTot);
  OIplusIOL2plusIOL1effEta->Divide(effEtaTot);
  OIplusIOL2plusIOL1effPhi->Divide(effPhiTot);
  OIplusIOL2plusIOL1effnVtx->Divide(effnVtxTot);



 
  //Writing the histograms in a file.
  outfile           -> cd();
  OIeffPt->Write(); 
  OIeffEta->Write();
  OIeffPhi->Write();
  OIeffnVtx->Write();
  OIplusIOL2effPt->Write();
  OIplusIOL2effEta->Write();
  OIplusIOL2effPhi->Write();
  OIplusIOL2effnVtx->Write();
  OIplusIOL2plusIOL1effPt->Write();
  OIplusIOL2plusIOL1effEta->Write();
  OIplusIOL2plusIOL1effPhi->Write();
  OIplusIOL2plusIOL1effnVtx->Write();
  


  OIPt->Write();
  OIIOPt->Write();
  OIIOIOPt->Write();

  OIEta->Write();
  OIIOEta->Write();
  OIIOIOEta->Write();

  OIPhi->Write();
  OIIOPhi->Write();
  OIIOIOPhi->Write();

  OInVtx->Write();
  OIIOnVtx->Write();
  OIIOIOnVtx->Write();

  ProbePt->Write();

 
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



