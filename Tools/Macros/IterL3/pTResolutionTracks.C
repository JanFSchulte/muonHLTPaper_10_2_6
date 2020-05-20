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

void pTResolutionTracks(TString inputfilename="/eos/uscms/store/user/bmahakud/ProductionHLTAN_LPC_IterL3HighStat/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ProductionHLTAN_LPC_IterL3HighStat/181130_193653/0000/muonNtupleIterL3.root", std::string effmeasured="MC2018"){

  int flavor=Sig::Prompt;

  bool doingL1 = thepassfilter.find("L1fL1") != std::string::npos; 

  TFile* outfile = TFile::Open(Form("%s_IterL3preFilter.root", effmeasured.c_str()),"RECREATE");
  std::cout << "output file: " << outfile -> GetName() << std::endl;

  TH1F* ptResL1=new TH1F("ptResL1","ptResL1",1000,-5,5 );
  TH2F* ptResVsPtL1 = new TH2F("ptResVsPtL1","ptResVsPtL1",1000,-5,5,100,0,1000);
  TH2F* ptResVsEtaL1 = new TH2F("ptResVsEtaL1","ptResVsEtaL1",1000,-5,5,100,-2.4,2.4);

  TH1F* ptResL2=new TH1F("ptResL2","ptResL2",1000,-5,5 );
  TH2F* ptResVsPtL2 = new TH2F("ptResVsPtL2","ptResVsPtL2",1000,-5,5,100,0,1000);
  TH2F* ptResVsEtaL2 = new TH2F("ptResVsEtaL2","ptResVsEtaL2",1000,-5,5,100,-2.4,2.4);

  TH1F* ptResL3OI=new TH1F("ptResL3OI","ptResL3OI",1000,-5,5 );
  TH2F* ptResVsPtL3OI = new TH2F("ptResVsPtL3OI","ptResVsPtL3OI",1000,-5,5,100,0,1000);
  TH2F* ptResVsEtaL3OI = new TH2F("ptResVsEtaL3OI","ptResVsEtaL3OI",1000,-5,5,100,-2.4,2.4);

  TH1F* ptResL3IOL2=new TH1F("ptResL3IOL2","ptResL3IOL2",1000,-5,5 );
  TH2F* ptResVsPtL3IOL2 = new TH2F("ptResVsPtL3IOL2","ptResVsPtL3IOL2",1000,-5,5,100,0,1000);
  TH2F* ptResVsEtaL3IOL2 = new TH2F("ptResVsEtaL3IOL2","ptResVsEtaL3IOL2",1000,-5,5,100,-2.4,2.4);

  TH1F* ptResL3IOL1=new TH1F("ptResL3IOL1","ptResL3IOL1",1000,-5,5 );
  TH2F* ptResVsPtL3IOL1 = new TH2F("ptResVsPtL3IOL1","ptResVsPtL3IOL1",1000,-5,5,100,0,1000);
  TH2F* ptResVsEtaL3IOL1 = new TH2F("ptResVsEtaL3IOL1","ptResVsEtaL3IOL1",1000,-5,5,100,-2.4,2.4);

  TH1F* ptResL1Inv=new TH1F("ptResL1Inv","ptResL1Inv",1000,-5,5 );
  TH2F* ptResVsPtL1Inv = new TH2F("ptResVsPtL1Inv","ptResVsPtL1Inv",1000,-2.5,2.5,100,0,1000);
  TH2F* ptResVsEtaL1Inv = new TH2F("ptResVsEtaL1Inv","ptResVsEtaL1Inv",1000,-2.5,2.5,100,-2.4,2.4);

  TH1F* ptResL2Inv=new TH1F("ptResL2Inv","ptResL2Inv",1000,-5,5 );
  TH2F* ptResVsPtL2Inv = new TH2F("ptResVsPtL2Inv","ptResVsPtL2Inv",1000,-2.5,2.5,100,0,1000);
  TH2F* ptResVsEtaL2Inv = new TH2F("ptResVsEtaL2Inv","ptResVsEtaL2Inv",1000,-2.5,2.5,100,-2.4,2.4);

  TH1F* ptResL3OIInv=new TH1F("ptResL3OIInv","ptResL3OIInv",1000,-5,5 );
  TH2F* ptResVsPtL3OIInv = new TH2F("ptResVsPtL3OIInv","ptResVsPtL3OIInv",1000,-2.5,2.5,100,0,1000);
  TH2F* ptResVsEtaL3OIInv = new TH2F("ptResVsEtaL3OIInv","ptResVsEtaL3OIInv",1000,-2.5,2.5,100,-2.4,2.4);

  TH1F* ptResL3IOL2Inv=new TH1F("ptResL3IOL2Inv","ptResL3IOL2Inv",1000,-5,5 );
  TH2F* ptResVsPtL3IOL2Inv = new TH2F("ptResVsPtL3IOL2Inv","ptResVsPtL3IOL2Inv",1000,-2.5,2.5,100,0,1000);
  TH2F* ptResVsEtaL3IOL2Inv = new TH2F("ptResVsEtaL3IOL2Inv","ptResVsEtaL3IOL2Inv",1000,-2.5,2.5,100,-2.4,2.4);

  TH1F* ptResL3IOL1Inv=new TH1F("ptResL3IOL1Inv","ptResL3IOL1Inv",1000,-5,5 );
  TH2F* ptResVsPtL3IOL1Inv = new TH2F("ptResVsPtL3IOL1Inv","ptResVsPtL3IOL1Inv",1000,-2.5,2.5,100,0,1000);
  TH2F* ptResVsEtaL3IOL1Inv = new TH2F("ptResVsEtaL3IOL1Inv","ptResVsEtaL3IOL1Inv",1000,-2.5,2.5,100,-2.4,2.4);






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
	HltTrackCand closestL3;
        for ( std::vector<HltTrackCand>::const_iterator it = ev -> hltTrackOI.begin(); it != ev -> hltTrackOI.end(); ++it ) {
    		theDR = deltaR(it -> eta, it -> phi, ev -> muons.at(jmu).eta,ev -> muons.at(jmu).phi);
   		 if (theDR < minDR){
      			minDR = theDR;
      			match = true;	
			closestL3 = *it;
    		}
  	}
	if (match){
		double ptResL3InvVal = (1./ ev -> muons.at(jmu).pt - 1. / closestL3.pt ) / (1./ev -> muons.at(jmu).pt);
		double ptResL3Val = (ev -> muons.at(jmu).pt - closestL3.pt ) / ev -> muons.at(jmu).pt;
		ptResL3OIInv->Fill(ptResL3InvVal);
		ptResL3OI->Fill(ptResL3Val);
		ptResVsPtL3OIInv->Fill(ptResL3InvVal,ev -> muons.at(jmu).pt);
		ptResVsPtL3OI->Fill(ptResL3Val,ev -> muons.at(jmu).pt);
		ptResVsEtaL3OIInv->Fill(ptResL3InvVal,ev -> muons.at(jmu).eta);
		ptResVsEtaL3OI->Fill(ptResL3Val,ev -> muons.at(jmu).eta);
	}
	// for L3 IOL2
	match = false;
        minDR = 0.1;
        theDR = 100;
	HltTrackCand closestL3IOL2;
        for ( std::vector<HltTrackCand>::const_iterator it = ev -> hltTrackIOL2.begin(); it != ev -> hltTrackIOL2.end(); ++it ) {
    		theDR = deltaR(it -> eta, it -> phi, ev -> muons.at(jmu).eta,ev -> muons.at(jmu).phi);
   		 if (theDR < minDR){
      			minDR = theDR;
      			match = true;	
			closestL3IOL2 = *it;
    		}
  	}
	if (match){
		double ptResL3InvVal = (1./ ev -> muons.at(jmu).pt - 1. / closestL3IOL2.pt ) / (1./ev -> muons.at(jmu).pt);
		double ptResL3Val = (ev -> muons.at(jmu).pt - closestL3IOL2.pt ) / ev -> muons.at(jmu).pt;
		ptResL3IOL2Inv->Fill(ptResL3InvVal);
		ptResL3IOL2->Fill(ptResL3Val);
		ptResVsPtL3IOL2Inv->Fill(ptResL3InvVal,ev -> muons.at(jmu).pt);
		ptResVsPtL3IOL2->Fill(ptResL3Val,ev -> muons.at(jmu).pt);
		ptResVsEtaL3IOL2Inv->Fill(ptResL3InvVal,ev -> muons.at(jmu).eta);
		ptResVsEtaL3IOL2->Fill(ptResL3Val,ev -> muons.at(jmu).eta);
	}
	// for L3 IOL1
	match = false;
        minDR = 0.1;
        theDR = 100;
	HltTrackCand closestL3IOL1;
        for ( std::vector<HltTrackCand>::const_iterator it = ev -> hltTrackIOL1.begin(); it != ev -> hltTrackIOL1.end(); ++it ) {
    		theDR = deltaR(it -> eta, it -> phi, ev -> muons.at(jmu).eta,ev -> muons.at(jmu).phi);
   		 if (theDR < minDR){
      			minDR = theDR;
      			match = true;	
			closestL3IOL1 = *it;
    		}
  	}
	if (match){
		double ptResL3InvVal = (1./ ev -> muons.at(jmu).pt - 1. / closestL3IOL1.pt ) / (1./ev -> muons.at(jmu).pt);
		double ptResL3Val = (ev -> muons.at(jmu).pt - closestL3IOL1.pt ) / ev -> muons.at(jmu).pt;
		ptResL3IOL1Inv->Fill(ptResL3InvVal);
		ptResL3IOL1->Fill(ptResL3Val);
		ptResVsPtL3IOL1Inv->Fill(ptResL3InvVal,ev -> muons.at(jmu).pt);
		ptResVsPtL3IOL1->Fill(ptResL3Val,ev -> muons.at(jmu).pt);
		ptResVsEtaL3IOL1Inv->Fill(ptResL3InvVal,ev -> muons.at(jmu).eta);
		ptResVsEtaL3IOL1->Fill(ptResL3Val,ev -> muons.at(jmu).eta);
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
		double ptResL2Val = (ev -> muons.at(jmu).pt - closestL2.pt ) / ev -> muons.at(jmu).pt;
                ptResL2Inv->Fill(ptResL2InvVal);
                ptResL2->Fill(ptResL2Val);
                ptResVsPtL2Inv->Fill(ptResL2InvVal,ev -> muons.at(jmu).pt);
                ptResVsPtL2->Fill(ptResL2Val,ev -> muons.at(jmu).pt);
                ptResVsEtaL2Inv->Fill(ptResL2InvVal,ev -> muons.at(jmu).eta);
                ptResVsEtaL2->Fill(ptResL2Val,ev -> muons.at(jmu).eta);

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
		double ptResL1Val = (ev -> muons.at(jmu).pt - closestL1.pt ) / ev -> muons.at(jmu).pt;
                ptResL1Inv->Fill(ptResL1InvVal);
                ptResL1->Fill(ptResL1Val);
                ptResVsPtL1Inv->Fill(ptResL1InvVal,ev -> muons.at(jmu).pt);
                ptResVsPtL1->Fill(ptResL1Val,ev -> muons.at(jmu).pt);
                ptResVsEtaL1Inv->Fill(ptResL1InvVal,ev -> muons.at(jmu).eta);
                ptResVsEtaL1->Fill(ptResL1Val,ev -> muons.at(jmu).eta);

	}
	

	
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

  ptResL3OI->Write();
  ptResVsPtL3OI->Write();
  ptResVsEtaL3OI->Write();

  ptResL3IOL2->Write();
  ptResVsPtL3IOL2->Write();
  ptResVsEtaL3IOL2->Write();

  ptResL3IOL1->Write();
  ptResVsPtL3IOL1->Write();
  ptResVsEtaL3IOL1->Write();

  ptResL1Inv->Write();
  ptResVsPtL1Inv->Write();
  ptResVsEtaL1Inv->Write();

  ptResL2Inv->Write();
  ptResVsPtL2Inv->Write();
  ptResVsEtaL2Inv->Write();

  ptResL3OIInv->Write();
  ptResVsPtL3OIInv->Write();
  ptResVsEtaL3OIInv->Write();

  ptResL3IOL2Inv->Write();
  ptResVsPtL3IOL2Inv->Write();
  ptResVsEtaL3IOL2Inv->Write();

  ptResL3IOL1Inv->Write();
  ptResVsPtL3IOL1Inv->Write();
  ptResVsEtaL3IOL1Inv->Write();


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
