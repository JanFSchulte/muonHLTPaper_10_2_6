#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TChain.h"
#include "TH2F.h"
#include "TTree.h"
#include "TBranch.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include "DataFormats/Math/interface/deltaR.h"
#include "MuonHLTNtuples/Analyzers/src/MuonTree.h"
#include "MuonHLTNtuples/Analyzers/src/MuonTreeLinkDef.h"
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
int getSign(int);
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

void pTResolution(TString inputfilename="/eos/uscms/store/user/bmahakud/ProductionHLTAN_LPC_IterL3HighStat/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ProductionHLTAN_LPC_IterL3HighStat/181130_193653/0000/muonNtupleIterL3.root", std::string effmeasured="MC2018"){

  int flavor=Sig::Prompt;

  bool doingL1 = thepassfilter.find("L1fL1") != std::string::npos; 

  TFile* outfile = TFile::Open(Form("%s_IterL3preFilter.root", effmeasured.c_str()),"RECREATE");
  std::cout << "output file: " << outfile -> GetName() << std::endl;


  TH1F* pTChargeAllL1 =new TH1F("pTChargeAllL1","pTChargeAllL1",19,  pt_bins );
  TH1F* pTChargeMisL1 =new TH1F("pTChargeMisL1","pTChargeMisL1",19,  pt_bins );
  TH1F* pTChargeAllL2 =new TH1F("pTChargeAllL2","pTChargeAllL2",19,  pt_bins );
  TH1F* pTChargeMisL2 =new TH1F("pTChargeMisL2","pTChargeMisL2",19,  pt_bins );
  TH1F* pTChargeAllL3 =new TH1F("pTChargeAllL3","pTChargeAllL3",19,  pt_bins );
  TH1F* pTChargeMisL3 =new TH1F("pTChargeMisL3","pTChargeMisL3",19,  pt_bins );

  TH1F* etaChargeAllL1 =new TH1F("etaChargeAllL1","etaChargeAllL1",15,  eta_bins );
  TH1F* etaChargeMisL1 =new TH1F("etaChargeMisL1","etaChargeMisL1",15,  eta_bins );
  TH1F* etaChargeAllL2 =new TH1F("etaChargeAllL2","etaChargeAllL2",15,  eta_bins );
  TH1F* etaChargeMisL2 =new TH1F("etaChargeMisL2","etaChargeMisL2",15,  eta_bins );
  TH1F* etaChargeAllL3 =new TH1F("etaChargeAllL3","etaChargeAllL3",15,  eta_bins );
  TH1F* etaChargeMisL3 =new TH1F("etaChargeMisL3","etaChargeMisL3",15,  eta_bins );

  TH1F* pTChargeAllGenL1 =new TH1F("pTChargeAllGenL1","pTChargeAllGenL1",19,  pt_bins );
  TH1F* pTChargeMisGenL1 =new TH1F("pTChargeMisGenL1","pTChargeMisGenL1",19,  pt_bins );
  TH1F* pTChargeAllGenL2 =new TH1F("pTChargeAllGenL2","pTChargeAllGenL2",19,  pt_bins );
  TH1F* pTChargeMisGenL2 =new TH1F("pTChargeMisGenL2","pTChargeMisGenL2",19,  pt_bins );
  TH1F* pTChargeAllGenL3 =new TH1F("pTChargeAllGenL3","pTChargeAllGenL3",19,  pt_bins );
  TH1F* pTChargeMisGenL3 =new TH1F("pTChargeMisGenL3","pTChargeMisGenL3",19,  pt_bins );

  TH1F* etaChargeAllGenL1 =new TH1F("etaChargeAllGenL1","etaChargeAllGenL1",15,  eta_bins );
  TH1F* etaChargeMisGenL1 =new TH1F("etaChargeMisGenL1","etaChargeMisGenL1",15,  eta_bins );
  TH1F* etaChargeAllGenL2 =new TH1F("etaChargeAllGenL2","etaChargeAllGenL2",15,  eta_bins );
  TH1F* etaChargeMisGenL2 =new TH1F("etaChargeMisGenL2","etaChargeMisGenL2",15,  eta_bins );
  TH1F* etaChargeAllGenL3 =new TH1F("etaChargeAllGenL3","etaChargeAllGenL3",15,  eta_bins );
  TH1F* etaChargeMisGenL3 =new TH1F("etaChargeMisGenL3","etaChargeMisGenL3",15,  eta_bins );







  TH1F* ptResL1=new TH1F("ptResL1","ptResL1",1000,-5,5 );
  TH2F* ptResVsPtL1 = new TH2F("ptResVsPtL1","ptResVsPtL1",1000,-1,1,100,0,1000);
  TH2F* ptResVsEtaL1 = new TH2F("ptResVsEtaL1","ptResVsEtaL1",1000,-5,5,100,-2.4,2.4);

  TH1F* ptResL2=new TH1F("ptResL2","ptResL2",1000,-5,5 );
  TH2F* ptResVsPtL2 = new TH2F("ptResVsPtL2","ptResVsPtL2",1000,-5,5,100,0,1000);
  TH2F* ptResVsEtaL2 = new TH2F("ptResVsEtaL2","ptResVsEtaL2",1000,-5,5,100,-2.4,2.4);

  TH1F* ptResL3=new TH1F("ptResL3","ptResL3",1000,-5,5 );
  TH2F* ptResVsPtL3 = new TH2F("ptResVsPtL3","ptResVsPtL3",1000,-2.5,2.5,100,0,1000);
  TH2F* ptResVsEtaL3 = new TH2F("ptResVsEtaL3","ptResVsEtaL3",1000,-2.5,2.5,100,-2.4,2.4);

  TH1F* ptResL1Inv=new TH1F("ptResL1Inv","ptResL1Inv",1000,-5,5 );
  TH2F* ptResVsPtL1Inv = new TH2F("ptResVsPtL1Inv","ptResVsPtL1Inv",1000,-2.5,2.5,100,0,1000);
  TH2F* ptResVsEtaL1Inv = new TH2F("ptResVsEtaL1Inv","ptResVsEtaL1Inv",1000,-2.5,2.5,100,-2.4,2.4);

  TH1F* ptResL2Inv=new TH1F("ptResL2Inv","ptResL2Inv",1000,-5,5 );
  TH2F* ptResVsPtL2Inv = new TH2F("ptResVsPtL2Inv","ptResVsPtL2Inv",1000,-2.5,2.5,100,0,1000);
  TH2F* ptResVsEtaL2Inv = new TH2F("ptResVsEtaL2Inv","ptResVsEtaL2Inv",1000,-2.5,2.5,100,-2.4,2.4);

  TH1F* ptResL3Inv=new TH1F("ptResL3Inv","ptResL3Inv",1000,-5,5 );
  TH2F* ptResVsPtL3Inv = new TH2F("ptResVsPtL3Inv","ptResVsPtL3Inv",1000,-2.5,2.5,100,0,1000);
  TH2F* ptResVsEtaL3Inv = new TH2F("ptResVsEtaL3Inv","ptResVsEtaL3Inv",1000,-2.5,2.5,100,-2.4,2.4);

  TH1F* ptResL1QPt=new TH1F("ptResL1QPt","ptResL1QPt",1000,-5,5 );
  TH2F* ptResVsPtL1QPt = new TH2F("ptResVsPtL1QPt","ptResVsPtL1QPt",1000,-2.5,2.5,100,0,1000);
  TH2F* ptResVsEtaL1QPt = new TH2F("ptResVsEtaL1QPt","ptResVsEtaL1QPt",1000,-2.5,2.5,100,-2.4,2.4);

  TH1F* ptResL2QPt=new TH1F("ptResL2QPt","ptResL2QPt",1000,-5,5 );
  TH2F* ptResVsPtL2QPt = new TH2F("ptResVsPtL2QPt","ptResVsPtL2QPt",1000,-2.5,2.5,100,0,1000);
  TH2F* ptResVsEtaL2QPt = new TH2F("ptResVsEtaL2QPt","ptResVsEtaL2QPt",1000,-2.5,2.5,100,-2.4,2.4);

  TH1F* ptResL3QPt=new TH1F("ptResL3QPt","ptResL3QPt",1000,-5,5 );
  TH2F* ptResVsPtL3QPt = new TH2F("ptResVsPtL3QPt","ptResVsPtL3QPt",1000,-2.5,2.5,100,0,1000);
  TH2F* ptResVsEtaL3QPt = new TH2F("ptResVsEtaL3QPt","ptResVsEtaL3QPt",1000,-2.5,2.5,100,-2.4,2.4);


  TH1F* ptResL1Gen=new TH1F("ptResL1Gen","ptResL1Gen",1000,-5,5 );
  TH2F* ptResVsPtL1Gen = new TH2F("ptResVsPtL1Gen","ptResVsPtL1Gen",1000,-1,1,100,0,1000);
  TH2F* ptResVsEtaL1Gen = new TH2F("ptResVsEtaL1Gen","ptResVsEtaL1Gen",1000,-5,5,100,-2.4,2.4);

  TH1F* ptResL2Gen = new TH1F("ptResL2Gen","ptResL2Gen",1000,-5,5 );
  TH2F* ptResVsPtL2Gen = new TH2F("ptResVsPtL2Gen","ptResVsPtL2Gen",1000,-5,5,100,0,1000);
  TH2F* ptResVsEtaL2Gen = new TH2F("ptResVsEtaL2Gen","ptResVsEtaL2Gen",1000,-5,5,100,-2.4,2.4);

  TH1F* ptResL3Gen = new TH1F("ptResL3Gen","ptResL3Gen",1000,-5,5 );
  TH2F* ptResVsPtL3Gen = new TH2F("ptResVsPtL3Gen","ptResVsPtL3Gen",1000,-2.5,2.5,100,0,1000);
  TH2F* ptResVsEtaL3Gen = new TH2F("ptResVsEtaL3Gen","ptResVsEtaL3Gen",1000,-2.5,2.5,100,-2.4,2.4);

  TH1F* ptResL1InvGen = new TH1F("ptResL1InvGen","ptResL1InvGen",1000,-5,5 );
  TH2F* ptResVsPtL1InvGen = new TH2F("ptResVsPtL1InvGen","ptResVsPtL1InvGen",1000,-2.5,2.5,100,0,1000);
  TH2F* ptResVsEtaL1InvGen = new TH2F("ptResVsEtaL1InvGen","ptResVsEtaL1InvGen",1000,-2.5,2.5,100,-2.4,2.4);

  TH1F* ptResL2InvGen = new TH1F("ptResL2InvGen","ptResL2InvGen",1000,-5,5 );
  TH2F* ptResVsPtL2InvGen = new TH2F("ptResVsPtL2InvGen","ptResVsPtL2InvGen",1000,-2.5,2.5,100,0,1000);
  TH2F* ptResVsEtaL2InvGen = new TH2F("ptResVsEtaL2InvGen","ptResVsEtaL2InvGen",1000,-2.5,2.5,100,-2.4,2.4);

  TH1F* ptResL3InvGen = new TH1F("ptResL3InvGen","ptResL3InvGen",1000,-5,5 );
  TH2F* ptResVsPtL3InvGen = new TH2F("ptResVsPtL3InvGen","ptResVsPtL3InvGen",1000,-2.5,2.5,100,0,1000);
  TH2F* ptResVsEtaL3InvGen = new TH2F("ptResVsEtaL3InvGen","ptResVsEtaL3InvGen",1000,-2.5,2.5,100,-2.4,2.4);

  TH1F* ptResL1QPtGen=new TH1F("ptResL1QPtGen","ptResL1QPtGen",1000,-5,5 );
  TH2F* ptResVsPtL1QPtGen = new TH2F("ptResVsPtL1QPtGen","ptResVsPtL1QPtGen",1000,-2.5,2.5,100,0,1000);
  TH2F* ptResVsEtaL1QPtGen = new TH2F("ptResVsEtaL1QPtGen","ptResVsEtaL1QPtGen",1000,-2.5,2.5,100,-2.4,2.4);

  TH1F* ptResL2QPtGen=new TH1F("ptResL2QPtGen","ptResL2QPtGen",1000,-5,5 );
  TH2F* ptResVsPtL2QPtGen = new TH2F("ptResVsPtL2QPtGen","ptResVsPtL2QPtGen",1000,-2.5,2.5,100,0,1000);
  TH2F* ptResVsEtaL2QPtGen = new TH2F("ptResVsEtaL2QPtGen","ptResVsEtaL2QPtGen",1000,-2.5,2.5,100,-2.4,2.4);

  TH1F* ptResL3QPtGen=new TH1F("ptResL3QPtGen","ptResL3QPtGen",1000,-5,5 );
  TH2F* ptResVsPtL3QPtGen = new TH2F("ptResVsPtL3QPtGen","ptResVsPtL3QPtGen",1000,-2.5,2.5,100,0,1000);
  TH2F* ptResVsEtaL3QPtGen = new TH2F("ptResVsEtaL3QPtGen","ptResVsEtaL3QPtGen",1000,-2.5,2.5,100,-2.4,2.4);





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
    
    unsigned int nmuons = ev->muons.size(); 
    if (nmuons < 2) continue; 
    for (int imu = 0; imu < nmuons; imu++){ 
      // select a good offline muon        
      if (debug) cout <<"select Tag muon" << endl;
      if (! selectTagMuon(ev -> muons.at(imu))) continue; 

      for (int jmu = 0; jmu < nmuons; jmu++){

	if (!selectProbeMuon(ev -> muons.at(jmu), ev -> muons.at(imu))) continue;

	// for L3
	bool match = false;
        float minDR = 0.1;
        float theDR = 100;
	HLTMuonCand closestL3;
        for ( std::vector<HLTMuonCand>::const_iterator it = ev -> hltmuons.begin(); it != ev -> hltmuons.end(); ++it ) {
    		theDR = deltaR(it -> eta, it -> phi, ev -> muons.at(jmu).eta,ev -> muons.at(jmu).phi);
   		 if (theDR < minDR){
      			minDR = theDR;
      			match = true;	
			closestL3 = *it;
    		}
  	}
	if (match){

		double ptResL3InvVal = (1./ ev -> muons.at(jmu).pt - 1. / closestL3.pt ) / (1./ev -> muons.at(jmu).pt);
		double ptResL3QPtVal = (ev -> muons.at(jmu).charge/ ev -> muons.at(jmu).pt - closestL3.charge / closestL3.pt ) / (ev -> muons.at(jmu).charge /ev -> muons.at(jmu).pt);
		double ptResL3Val = (ev -> muons.at(jmu).pt - closestL3.pt ) / ev -> muons.at(jmu).pt;
		ptResL3Inv->Fill(ptResL3InvVal);
		ptResL3QPt->Fill(ptResL3QPtVal);
		ptResL3->Fill(ptResL3Val);
		ptResVsPtL3Inv->Fill(ptResL3InvVal,ev -> muons.at(jmu).pt);
		ptResVsPtL3QPt->Fill(ptResL3QPtVal,ev -> muons.at(jmu).pt);
		ptResVsPtL3->Fill(ptResL3Val,ev -> muons.at(jmu).pt);
		ptResVsEtaL3Inv->Fill(ptResL3InvVal,ev -> muons.at(jmu).eta);
		ptResVsEtaL3QPt->Fill(ptResL3QPtVal,ev -> muons.at(jmu).eta);
		ptResVsEtaL3->Fill(ptResL3Val,ev -> muons.at(jmu).eta);


		pTChargeAllL3->Fill(ev -> muons.at(jmu).pt);
		if (! (ev -> muons.at(jmu).charge == closestL3.charge)) pTChargeMisL3->Fill(ev -> muons.at(jmu).pt);
		etaChargeAllL3->Fill(ev -> muons.at(jmu).eta);
		if (! (ev -> muons.at(jmu).charge == closestL3.charge)) etaChargeMisL3->Fill(ev -> muons.at(jmu).eta);
	}
	
	// for L2
	match = false;
        minDR = 0.1;
        theDR = 100;
	HLTMuonCand closestL2;
        for ( std::vector<HLTMuonCand>::const_iterator it = ev -> L2muons.begin(); it != ev -> L2muons.end(); ++it ) {
    		theDR = deltaR(it -> eta, it -> phi, ev -> muons.at(jmu).eta,ev -> muons.at(jmu).phi);
   		 if (theDR < minDR){
      			minDR = theDR;
      			match = true;	
			closestL2 = *it;
    		}
  	}
	if (match){
		double ptResL2InvVal = (1./ ev -> muons.at(jmu).pt - 1. / closestL2.pt ) / (1./ev -> muons.at(jmu).pt);
		double ptResL2QPtVal = (ev -> muons.at(jmu).charge / ev -> muons.at(jmu).pt - closestL2.charge / closestL2.pt ) / (ev -> muons.at(jmu).charge / ev -> muons.at(jmu).pt);
		double ptResL2Val = (ev -> muons.at(jmu).pt - closestL2.pt ) / ev -> muons.at(jmu).pt;
                ptResL2Inv->Fill(ptResL2InvVal);
                ptResL2QPt->Fill(ptResL2QPtVal);
                ptResL2->Fill(ptResL2Val);
                ptResVsPtL2Inv->Fill(ptResL2InvVal,ev -> muons.at(jmu).pt);
                ptResVsPtL2QPt->Fill(ptResL2QPtVal,ev -> muons.at(jmu).pt);
                ptResVsPtL2->Fill(ptResL2Val,ev -> muons.at(jmu).pt);
                ptResVsEtaL2Inv->Fill(ptResL2InvVal,ev -> muons.at(jmu).eta);
                ptResVsEtaL2QPt->Fill(ptResL2QPtVal,ev -> muons.at(jmu).eta);
                ptResVsEtaL2->Fill(ptResL2Val,ev -> muons.at(jmu).eta);

		pTChargeAllL2->Fill(ev -> muons.at(jmu).pt);
		if (! (ev -> muons.at(jmu).charge == closestL2.charge)) pTChargeMisL2->Fill(ev -> muons.at(jmu).pt);
		etaChargeAllL2->Fill(ev -> muons.at(jmu).eta);
		if (! (ev -> muons.at(jmu).charge == closestL2.charge)) etaChargeMisL2->Fill(ev -> muons.at(jmu).eta);

	}
	

	// for L2
	match = false;
        minDR = 0.3;
        theDR = 100;
	L1MuonCand closestL1;
        for ( std::vector<L1MuonCand>::const_iterator it = ev -> L1muons.begin(); it != ev -> L1muons.end(); ++it ) {
		if (it->quality < 12) continue;
    		theDR = deltaR(it -> eta, it -> phi, ev -> muons.at(jmu).eta,ev -> muons.at(jmu).phi);
   		 if (theDR < minDR){
      			minDR = theDR;
      			match = true;	
			closestL1 = *it;
    		}
  	}
	if (match){
		double ptResL1InvVal = (1./ ev -> muons.at(jmu).pt - 1. / closestL1.pt ) / (1./ev -> muons.at(jmu).pt);
		double ptResL1QPtVal = (ev -> muons.at(jmu).charge / ev -> muons.at(jmu).pt - closestL1.charge / closestL1.pt ) / (ev -> muons.at(jmu).charge / ev -> muons.at(jmu).pt);
		double ptResL1Val = (ev -> muons.at(jmu).pt - closestL1.pt ) / ev -> muons.at(jmu).pt;
                ptResL1Inv->Fill(ptResL1InvVal);
                ptResL1QPt->Fill(ptResL1QPtVal);
                ptResL1->Fill(ptResL1Val);
                ptResVsPtL1Inv->Fill(ptResL1InvVal,ev -> muons.at(jmu).pt);
                ptResVsPtL1QPt->Fill(ptResL1QPtVal,ev -> muons.at(jmu).pt);
                ptResVsPtL1->Fill(ptResL1Val,ev -> muons.at(jmu).pt);
                ptResVsEtaL1Inv->Fill(ptResL1InvVal,ev -> muons.at(jmu).eta);
                ptResVsEtaL1QPt->Fill(ptResL1QPtVal,ev -> muons.at(jmu).eta);
                ptResVsEtaL1->Fill(ptResL1Val,ev -> muons.at(jmu).eta);

		pTChargeAllL1->Fill(ev -> muons.at(jmu).pt);
		if (! (ev -> muons.at(jmu).charge == closestL1.charge)) pTChargeMisL1->Fill(ev -> muons.at(jmu).pt);
		etaChargeAllL1->Fill(ev -> muons.at(jmu).eta);
		if (! (ev -> muons.at(jmu).charge == closestL1.charge)) etaChargeMisL1->Fill(ev -> muons.at(jmu).eta);

	}
	

	
     }


    }
       int nGenMuons = ev->genParticles.size();
      for (int jmu = 0; jmu < nGenMuons; jmu++){

	//if (!selectProbeMuon(ev -> muons.at(jmu), ev -> muons.at(imu))) continue;
	if (ev -> genParticles.at(jmu).pt < 26) continue;
	// for L3
	bool match = false;
        float minDR = 0.1;
        float theDR = 100;
	HLTMuonCand closestL3;
        for ( std::vector<HLTMuonCand>::const_iterator it = ev -> hltmuons.begin(); it != ev -> hltmuons.end(); ++it ) {
    		theDR = deltaR(it -> eta, it -> phi, ev -> genParticles.at(jmu).eta,ev -> genParticles.at(jmu).phi);
   		 if (theDR < minDR){
      			minDR = theDR;
      			match = true;	
			closestL3 = *it;
    		}
  	}
	if (match){
		double ptResL3InvVal = (1./ ev -> genParticles.at(jmu).pt - 1. / closestL3.pt ) / (1./ev -> genParticles.at(jmu).pt);
		double ptResL3QPtVal = (getSign(ev -> genParticles.at(jmu).pdgId) / ev -> genParticles.at(jmu).pt - closestL3.charge / closestL3.pt ) / (getSign(ev-> genParticles.at(jmu).pdgId) / ev -> genParticles.at(jmu).pt);
		double ptResL3Val = (ev -> genParticles.at(jmu).pt - closestL3.pt ) / ev -> genParticles.at(jmu).pt;
		ptResL3InvGen->Fill(ptResL3InvVal);
		ptResL3QPtGen->Fill(ptResL3QPtVal);
		ptResL3Gen->Fill(ptResL3Val);
		ptResVsPtL3InvGen->Fill(ptResL3InvVal,ev -> genParticles.at(jmu).pt);
		ptResVsPtL3QPtGen->Fill(ptResL3QPtVal,ev -> genParticles.at(jmu).pt);
		ptResVsPtL3Gen->Fill(ptResL3Val,ev -> genParticles.at(jmu).pt);
		ptResVsEtaL3InvGen->Fill(ptResL3InvVal,ev -> genParticles.at(jmu).eta);
		ptResVsEtaL3QPtGen->Fill(ptResL3QPtVal,ev -> genParticles.at(jmu).eta);
		ptResVsEtaL3Gen->Fill(ptResL3Val,ev -> genParticles.at(jmu).eta);

		pTChargeAllGenL3->Fill(ev -> genParticles.at(jmu).pt);
		if (! (getSign(ev -> genParticles.at(jmu).pdgId) == closestL3.charge)) pTChargeMisGenL3->Fill(ev -> genParticles.at(jmu).pt);
		etaChargeAllGenL3->Fill(ev -> genParticles.at(jmu).eta);
		if (! (getSign(ev -> genParticles.at(jmu).pdgId) == closestL3.charge)) etaChargeMisGenL3->Fill(ev -> genParticles.at(jmu).eta);

	}
	
	// for L2
	match = false;
        minDR = 0.1;
        theDR = 100;
	HLTMuonCand closestL2;
        for ( std::vector<HLTMuonCand>::const_iterator it = ev -> L2muons.begin(); it != ev -> L2muons.end(); ++it ) {
    		theDR = deltaR(it -> eta, it -> phi, ev -> genParticles.at(jmu).eta,ev -> genParticles.at(jmu).phi);
   		 if (theDR < minDR){
      			minDR = theDR;
      			match = true;	
			closestL2 = *it;
    		}
  	}
	if (match){
		double ptResL2InvVal = (1./ ev -> genParticles.at(jmu).pt - 1. / closestL2.pt ) / (1./ev -> genParticles.at(jmu).pt);
		double ptResL2QPtVal = (getSign(ev -> genParticles.at(jmu).pdgId) / ev -> genParticles.at(jmu).pt - closestL2.charge / closestL2.pt ) / (getSign(ev -> genParticles.at(jmu).pdgId) /ev -> genParticles.at(jmu).pt);
		double ptResL2Val = (ev -> genParticles.at(jmu).pt - closestL2.pt ) / ev -> genParticles.at(jmu).pt;
                ptResL2InvGen->Fill(ptResL2InvVal);
                ptResL2QPtGen->Fill(ptResL2QPtVal);
                ptResL2Gen->Fill(ptResL2Val);
                ptResVsPtL2InvGen->Fill(ptResL2InvVal,ev -> genParticles.at(jmu).pt);
                ptResVsPtL2QPtGen->Fill(ptResL2QPtVal,ev -> genParticles.at(jmu).pt);
                ptResVsPtL2Gen->Fill(ptResL2Val,ev -> genParticles.at(jmu).pt);
                ptResVsEtaL2InvGen->Fill(ptResL2InvVal,ev -> genParticles.at(jmu).eta);
                ptResVsEtaL2QPtGen->Fill(ptResL2QPtVal,ev -> genParticles.at(jmu).eta);
                ptResVsEtaL2Gen->Fill(ptResL2Val,ev -> genParticles.at(jmu).eta);


		pTChargeAllGenL2->Fill(ev -> genParticles.at(jmu).pt);
		if (! (getSign(ev -> genParticles.at(jmu).pdgId) == closestL2.charge)) pTChargeMisGenL2->Fill(ev -> genParticles.at(jmu).pt);
		etaChargeAllGenL2->Fill(ev -> genParticles.at(jmu).eta);
		if (! (getSign(ev -> genParticles.at(jmu).pdgId) == closestL2.charge)) etaChargeMisGenL2->Fill(ev -> genParticles.at(jmu).eta);


	}
	

	// for L2
	match = false;
        minDR = 0.3;
        theDR = 100;
	L1MuonCand closestL1;
        for ( std::vector<L1MuonCand>::const_iterator it = ev -> L1muons.begin(); it != ev -> L1muons.end(); ++it ) {
		if (it->quality < 12) continue;
    		theDR = deltaR(it -> eta, it -> phi, ev -> genParticles.at(jmu).eta,ev -> genParticles.at(jmu).phi);
   		 if (theDR < minDR){
      			minDR = theDR;
      			match = true;	
			closestL1 = *it;
    		}
  	}
	if (match){
		double ptResL1InvVal = (1./ ev -> genParticles.at(jmu).pt - 1. / closestL1.pt ) / (1./ev -> genParticles.at(jmu).pt);
		double ptResL1QPtVal = ( getSign(ev -> genParticles.at(jmu).pdgId) / ev -> genParticles.at(jmu).pt - closestL2.charge / closestL1.pt ) / (getSign( ev -> genParticles.at(jmu).pdgId) / ev -> genParticles.at(jmu).pt);
		double ptResL1Val = (ev -> genParticles.at(jmu).pt - closestL1.pt ) / ev -> genParticles.at(jmu).pt;
                ptResL1InvGen->Fill(ptResL1InvVal);
                ptResL1QPtGen->Fill(ptResL1QPtVal);
                ptResL1Gen->Fill(ptResL1Val);
                ptResVsPtL1InvGen->Fill(ptResL1InvVal,ev -> genParticles.at(jmu).pt);
                ptResVsPtL1QPtGen->Fill(ptResL1QPtVal,ev -> genParticles.at(jmu).pt);
                ptResVsPtL1Gen->Fill(ptResL1Val,ev -> genParticles.at(jmu).pt);
                ptResVsEtaL1InvGen->Fill(ptResL1InvVal,ev -> genParticles.at(jmu).eta);
                ptResVsEtaL1QPtGen->Fill(ptResL1QPtVal,ev -> genParticles.at(jmu).eta);
                ptResVsEtaL1Gen->Fill(ptResL1Val,ev -> genParticles.at(jmu).eta);

		pTChargeAllGenL1->Fill(ev -> genParticles.at(jmu).pt);
		if (! (getSign(ev -> genParticles.at(jmu).pdgId) == closestL1.charge)) pTChargeMisGenL1->Fill(ev -> genParticles.at(jmu).pt);
		etaChargeAllGenL1->Fill(ev -> genParticles.at(jmu).eta);
		if (! (getSign(ev -> genParticles.at(jmu).pdgId) == closestL1.charge)) etaChargeMisGenL1->Fill(ev -> genParticles.at(jmu).eta);


	}
	

	
     }
 
  } 
 //Writing the histograms in a file.
  outfile           -> cd();


  ptResL1->Write();
  ptResVsPtL1->Write();
  ptResVsEtaL1->Write();

  ptResL2->Write();
  ptResVsPtL2->Write();
  ptResVsEtaL2->Write();

  ptResL3->Write();
  ptResVsPtL3->Write();
  ptResVsEtaL3->Write();

  ptResL1Inv->Write();
  ptResVsPtL1Inv->Write();
  ptResVsEtaL1Inv->Write();

  ptResL2Inv->Write();
  ptResVsPtL2Inv->Write();
  ptResVsEtaL2Inv->Write();

  ptResL3Inv->Write();
  ptResVsPtL3Inv->Write();
  ptResVsEtaL3Inv->Write();

  ptResL1QPt->Write();
  ptResVsPtL1QPt->Write();
  ptResVsEtaL1QPt->Write();

  ptResL2QPt->Write();
  ptResVsPtL2QPt->Write();
  ptResVsEtaL2QPt->Write();

  ptResL3QPt->Write();
  ptResVsPtL3QPt->Write();
  ptResVsEtaL3QPt->Write();


  ptResL1Gen->Write();
  ptResVsPtL1Gen->Write();
  ptResVsEtaL1Gen->Write();

  ptResL2Gen->Write();
  ptResVsPtL2Gen->Write();
  ptResVsEtaL2Gen->Write();

  ptResL3Gen->Write();
  ptResVsPtL3Gen->Write();
  ptResVsEtaL3Gen->Write();

  ptResL1InvGen->Write();
  ptResVsPtL1InvGen->Write();
  ptResVsEtaL1InvGen->Write();

  ptResL2InvGen->Write();
  ptResVsPtL2InvGen->Write();
  ptResVsEtaL2InvGen->Write();

  ptResL3InvGen->Write();
  ptResVsPtL3InvGen->Write();
  ptResVsEtaL3InvGen->Write();

  ptResL1QPtGen->Write();
  ptResVsPtL1QPtGen->Write();
  ptResVsEtaL1QPtGen->Write();

  ptResL2QPtGen->Write();
  ptResVsPtL2QPtGen->Write();
  ptResVsEtaL2QPtGen->Write();

  ptResL3QPtGen->Write();
  ptResVsPtL3QPtGen->Write();
  ptResVsEtaL3QPtGen->Write();



  pTChargeAllL1->Write();
  pTChargeMisL1->Write();
  pTChargeAllL2->Write();
  pTChargeMisL2->Write();
  pTChargeAllL3->Write();
  pTChargeMisL3->Write();

  etaChargeAllL1->Write();
  etaChargeMisL1->Write();
  etaChargeAllL2->Write();
  etaChargeMisL2->Write();
  etaChargeAllL3->Write();
  etaChargeMisL3->Write();

  pTChargeAllGenL1->Write();
  pTChargeMisGenL1->Write();
  pTChargeAllGenL2->Write();
  pTChargeMisGenL2->Write();
  pTChargeAllGenL3->Write();
  pTChargeMisGenL3->Write();

  etaChargeAllGenL1->Write();
  etaChargeMisGenL1->Write();
  etaChargeAllGenL2->Write();
  etaChargeMisGenL2->Write();
  etaChargeAllGenL3->Write();
  etaChargeMisGenL3->Write();



  outfile->Close();

 
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

int getSign(int pdgId){

if (pdgId < 0) return 1;
else return -1;

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
