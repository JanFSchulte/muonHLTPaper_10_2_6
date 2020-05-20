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
bool matchMuonWithL3(MuonCand, std::vector<HLTMuonCand>);
std::string getProbeFilter(int, bool);
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
4226.67,
1823.04,
787.477,
486.122,
211.872,
118.316,
63.4982,
37.9888,
23.0261,
14.9799,
10.2863,
7.27979,
5.34811,
4.10273,
3.26837,
2.64915,
2.23506,
1.93052,
1.69373,
1.51683,
1.36889,
1.26224,
1.16777,
1.08684,
1.01703,
0.939873,
0.876546,
0.810688,
0.748066,
0.685921,
0.625705,
0.571451,
0.516831,
0.471159,
0.42725,
0.387711,
0.351214,
0.32003,
0.290533,
0.267097,
0.248765,
0.234077,
0.2196,
0.205689,
0.197416,
0.191746,
0.183693,
0.175746,
0.175846,
0.176275,
0.173514,
0.172328,
0.175078,
0.181157,
0.182818,
0.177181,
0.1876,
0.196788,
0.192088,
0.197586,
0.213389,
0.212861,
0.219671,
0.216015,
0.237371,
0.227794,
0.216879,
0.261989,
0.28839,
0.283026,
0.264886,
0.267633,
0.315637,
0.295262,
0.266044,
0.337673,
0.35408,
0.311411,
0.369371,
0.394094,
0.382813,
0.290877,
0.582273,
0.413041,
0.433895,
0.446835,
0.618919,
0.29157,
0.669854,
0.714608,
0.672572,
0.423471,
0.635207,
0.635207,
0.554362,
1.34514,
0.635207,
0.557742,
0.680579,
0.846942,
0.190562,
1.5245,
2.11736,
0.423471,
0.879517,
0.635207,
1.27041,
2.85843,
1.90562,
2.28674,
0,
0.95281,
0,
7.62248,
};

/// TAG-DEFINITION: 
std::string isofilterTag  = "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07::HLT";

/// for PROMPT-MUONS   (close-by and far-away) 
//std::string L1filter      = "hltL1fL1sMu22or25L1Filtered0::TEST"; 
//std::string L2filter      = "hltL2fL1sMu22or25L1f0L2Filtered10Q::TEST";
//std::string L3filter      = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q::TEST"; 
std::string L1filter      = "hltL1fL1sMu22or25L1Filtered0::MYHLT"; 
std::string L2filter      = "hltL2fL1sMu22or25L1f0L2Filtered10Q::MYHLT";
std::string L3filter      = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q::MYHLT"; 


/* 
/// for DIMUONS
std::string L1filter      = "hltL1fL1sDoubleMu155L1Filtered0::TEST"; 
std::string L2filter      = "hltL2pfL1sDoubleMu155L1f0L2PreFiltered0::TEST";
std::string L3filter      = "hltL3fL1DoubleMu155fFiltered17::TEST"; 

/// for DIMUONS (JPSI)
std::string L1filter      = "hltL1fL1sL1sDoubleMu4SQOSdRMax1p2L1Filtered0::TEST"; 
std::string L2filter      = "hltL2fL1sL1sDoubleMu4SQOSdRMax1p2L1f0L2PreFiltered0::TEST";
std::string L3filter      = "hltDimuon25JpsiL3fL3Filtered::TEST"; 

/// for DISPLACED: 
std::string L1filter      = "hltL1fDimuonL1Filtered0::TEST"; 
std::string L2filter      = "hltL2fDimuonL1f0L2NoVtxFiltered16::TEST";
std::string L3filter      = "hltL3fDimuonL1f0L2NVf16L3NoFiltersNoVtxFiltered43::TEST"; 

std::string L1filter      = "hltDimuon3L1Filtered0::TEST"; 
std::string L2filter      = "hltDimuon3L2PreFiltered0::TEST";
std::string L3filter      = "hltDoubleMu3L3FilteredNoVtx::TEST"; 
*/
// ******************************************
//       T&P definitions                    *
//                                          *
std::string thepassfilter  = L3filter;
//std::string theprobefilter = L1filter; 
float offlinePtCut         = 24.;
//                                          *
//                                          *
// ******************************************

void readNtuplesPrefilter_IterL3(TString inputfilename="/eos/uscms/store/user/bmahakud/ProductionHLTAN_LPC_IterL3HighStat/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ProductionHLTAN_LPC_IterL3HighStat/181130_193653/0000/muonNtupleIterL3.root", std::string effmeasured="MC2018", bool isMC = false){

  int flavor=Sig::Prompt;

  bool doingL1 = thepassfilter.find("L1fL1") != std::string::npos; 

  TFile* outfile = TFile::Open(Form("%s_IterL3preFilter.root", effmeasured.c_str()),"RECREATE");
  std::cout << "output file: " << outfile -> GetName() << std::endl;

  //Create histograms  
  TH1F* dimuon_mass             = new TH1F("h_dimuon_mass"          ,"dimuon_mass"      , 1500,  0,  150 );
  TH1F* tagiso                  = new TH1F("h_tagiso"               ,"tagiso"           ,  100,  0,  1   );
  TH1F* tagMuonPt               = new TH1F("h_tagMuonPt"            ,"tagMuonPt"        ,  150,  0,  150 );
  TH1F* nvtx_event              = new TH1F("h_nvtx_event"           ,"nvtx_event"       ,   60,  0,   60 );
 
  TEfficiency* muonPt_barrel    = new TEfficiency("muonPt_barrel"   ,"muonPt_barrel"    ,   19,  pt_bins );
  TEfficiency* muonPt_endcap    = new TEfficiency("muonPt_endcap"   ,"muonPt_endcap"    ,   19,  pt_bins );
  TEfficiency* muonPt           = new TEfficiency("muonPt"          ,"muonPt"           ,   19,  pt_bins ); 
  TEfficiency* muonPtTurnOn     = new TEfficiency("muonPtTurnOn"    ,"muonPtTurnOn"     ,   19,  pt_bins ); 
  TEfficiency* muonEta          = new TEfficiency("muonEta"         ,"muonEta"          ,   15, eta_bins );
  TEfficiency* muonPhi          = new TEfficiency("muonPhi"         ,"muonPhi"          ,   20, -3.2, 3.2);
  TEfficiency* muonEff          = new TEfficiency("muonEff"         ,"muonEff"          ,    1,   0., 1.0);
  TEfficiency* muonDeltaR       = new TEfficiency("muonDeltaR"      ,"muonDeltaR"       ,   30,   0., 3.0);
  TEfficiency* muonDeltaPhi     = new TEfficiency("muonDeltaPhi"    ,"muonDeltaPhi"     ,   30,   0., 3.2); 

  // GLOBAL QUANTITIES//
  TEfficiency* muonchi2         = new TEfficiency("muonchi2"        , "muonchi2"        ,  60,    0.,  7.);
  TEfficiency* muondxy          = new TEfficiency("muondxy"         , "muondxy"         ,  100,  -0.3, 0.3);
  TEfficiency* muondz           = new TEfficiency("muondz"          , "muondz"          ,   10,   dz_bins);
  TEfficiency* muonPixHit       = new TEfficiency("muonPixHit"      , "muonPixHit"      ,   20,  -0.5,19.5);
  TEfficiency* muonLayHit       = new TEfficiency("muonLayHit"      , "muonLayHit"      ,   16,   2.5,18.5);
  TEfficiency* muonPixLay       = new TEfficiency("muonPixLay"      , "muonPixLay"      ,    7,  -0.5, 6.5);
  TEfficiency* muoninnerPt      = new TEfficiency("muoninnerPt"     , "muoninnerPt"     ,   19,   pt_bins );
  TEfficiency* muoninnerEta     = new TEfficiency("muoninnerEta"    , "muoninnerEta"    ,   15,   eta_bins);
  TEfficiency* muoninnerPhi     = new TEfficiency("muoninnerPhi"    , "muoninnerPhi"    ,   20,  -3.2, 3.2);
  
  TEfficiency* failingMuonPt    = new TEfficiency("failingMuonPt"   ,"failingMuonPt"    ,   19,  pt_bins ); 
  TEfficiency* failingMuonEta   = new TEfficiency("failingMuonEta"  ,"failingMuonEta"   ,   15, eta_bins );
  TEfficiency* failingMuonPhi   = new TEfficiency("failingMuonPhi"  ,"failingMuonPhi"   ,   20, -3.2, 3.2);
  TEfficiency* failingMuonEff   = new TEfficiency("failingMuonEff"  ,"failingMuonEff"   ,   1 ,   0., 1.0);

  TH1F* ProbePt                 = new TH1F("h_ProbePt"              ,"ProbeMuonPt"      ,  10000, 0, 1000 );
  TH1F* PassingProbePt          = new TH1F("h_PassingProbePt"       ,"PassingMuonPt"    ,  19,  pt_bins );
  TH1F* PassingProbeEta         = new TH1F("h_PassingProbeEta"      ,"PassingMuonEta"   ,  15, eta_bins );
  TH1F* PassingProbePhi         = new TH1F("h_PassingProbePhi"      ,"PassingMuonPhi"   ,  20, -3.2, 3.2);
  TH1F* PassingProbeMll         = new TH1F("h_PassingProbeMll"      ,"PassingMuonMll"   ,  20,  86., 96.); 

  TH1F* FailingProbePt          = new TH1F("h_FailingProbePt"       ,"FailingMuonPt"    ,  19,  pt_bins );
  TH1F* FailingProbeEta         = new TH1F("h_FailingProbeEta"      ,"FailingMuonEta"   ,  15, eta_bins );
  TH1F* FailingProbePhi         = new TH1F("h_FailingProbePhi"      ,"FailingMuonPhi"   ,  20, -3.2, 3.2);
  TH1F* FailingProbeMll         = new TH1F("h_FailingProbeMll"      ,"FailingMuonMll"   ,  20,  86., 96.);

  // di-muon efficiencies 
  TEfficiency* diMuonPt         = new TEfficiency("diMuonPt"      ,"diMuonPt"       ,   19,  pt_bins  , 19,  pt_bins ); 
  TEfficiency* diMuonEta        = new TEfficiency("diMuonEta"     ,"diMuonEta"      ,   15, eta_bins  , 15, eta_bins );
  TEfficiency* diMuonPhi        = new TEfficiency("diMuonPhi"     ,"diMuonPhi"      ,   20, -3.2, 3.2 , 20, -3.2, 3.2);
  TEfficiency* diMuonEff        = new TEfficiency("diMuonEff"     ,"diMuonEff"      ,    1,   0., 1.0);
  TEfficiency* diMuonDeltaR     = new TEfficiency("diMuonDeltaR"  ,"diMuonDeltaR"   ,   30,   0., 3.0);
  TEfficiency* diMuonLeadPt     = new TEfficiency("diMuonLeadPt"  ,"diMuonLeadPt"   ,   19,  pt_bins ); 
  TEfficiency* diMuonLeadEta    = new TEfficiency("diMuonLeadEta" ,"diMuonLeadEta"  ,   15,  eta_bins );
  TEfficiency* diMuonLeadPhi    = new TEfficiency("diMuonLeadPhi" ,"diMuonLeadPhi"  ,   20, -3.2, 3.2);
  TEfficiency* diMuonTrailPt    = new TEfficiency("diMuonTrailPt" ,"diMuonTrailPt"  ,   19,  pt_bins ); 
  TEfficiency* diMuonTrailEta   = new TEfficiency("diMuonTrailEta","diMuonTrailEta" ,   15, eta_bins );
  TEfficiency* diMuonTrailPhi   = new TEfficiency("diMuonTrailPhi","diMuonTrailPhi" ,   20, -3.2, 3.2);

  TEfficiency* nvtx             = new TEfficiency("nvtx"             ,"nvtx"             ,   60,    0,  60);
  TEfficiency* nvtx_barrel      = new TEfficiency("nvtx_barrel"      ,"nvtx_barrel"      ,   60,    0,  60);
  TEfficiency* nvtx_endcap      = new TEfficiency("nvtx_endcap"      ,"nvtx_endcap"      ,   60,    0,  60);
   
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
  std::cout << "Number of entries = " << nentries << std::endl;

  bool flagfile = false;
  std::string theprobefilter = getProbeFilter(flavor, isMC);
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
    nvtx_event-> Fill( ev -> nVtx   ); 
    for (int imu = 0; imu < nmuons; imu++){ 
      // select the tag muon        
      if (debug) cout <<"select Tag muon" << endl;
      if (! selectTagMuon(ev -> muons.at(imu), tagiso)) continue; 
      
      if (! matchMuon(ev -> muons.at(imu), ev -> hltTag.objects, isofilterTag)) continue;
      tagMuonPt -> Fill ( ev -> muons.at(imu).pt) ; 
      
      for (int jmu = 0; jmu < nmuons; jmu++){
	bool pass   = false;
	
	// select the probe muon  & match to the probe:  
	if (!selectProbeMuon(ev -> muons.at(jmu), ev -> muons.at(imu), dimuon_mass)) continue;
	if (!doingL1 && !(matchMuon(ev -> muons.at(jmu), ev -> hlt.objects, theprobefilter))) continue;	
	// select the pass muon
	if (matchMuonWithL3(ev->muons.at(jmu),ev->hltmuons)) pass = true;
	
	muonPtTurnOn -> Fill( pass, ev -> muons.at(jmu).pt );
        ProbePt      -> Fill( ev -> muons.at(jmu).pt ); 
	if (ev -> muons.at(jmu).pt < offlinePtCut) continue;
	
	TLorentzVector mu1, mu2;
	mu1.SetPtEtaPhiM (ev->muons.at(imu).pt,ev->muons.at(imu).eta,ev->muons.at(imu).phi, muonmass); 
	mu2.SetPtEtaPhiM (ev->muons.at(jmu).pt,ev->muons.at(jmu).eta,ev->muons.at(jmu).phi, muonmass);
	double mumumass = (mu1 + mu2).M();
	double DeltaR   = mu1.DeltaR(mu2);
	double DeltaPhi = mu1.DeltaPhi(mu2);

	if (pass) { 
	  PassingProbePt  -> Fill( ev -> muons.at(jmu).pt  );
	  PassingProbeEta -> Fill( ev -> muons.at(jmu).eta );
	  PassingProbePhi -> Fill( ev -> muons.at(jmu).phi );
	  PassingProbeMll -> Fill( mumumass                );
	}	      
	else {       
	  FailingProbePt  -> Fill( ev -> muons.at(jmu).pt  );
	  FailingProbeEta -> Fill( ev -> muons.at(jmu).eta );
	  FailingProbePhi -> Fill( ev -> muons.at(jmu).phi );
	  FailingProbeMll -> Fill( mumumass                );
	}
	if (isMC){	
		muonPt       -> FillWeighted( pass, weight, ev -> muons.at(jmu).pt );
		muonEta      -> FillWeighted( pass, weight, ev -> muons.at(jmu).eta);
		muonPhi      -> FillWeighted( pass, weight, ev -> muons.at(jmu).phi);
		muonEff      -> FillWeighted( pass, weight, 0.5                    );
		muonDeltaR   -> FillWeighted( pass, weight, DeltaR);
		muonDeltaPhi -> FillWeighted( pass, weight, DeltaPhi);
		nvtx	     -> FillWeighted( pass, weight, ev -> nVtx);

		muoninnerPt  -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerpt);
		muoninnerEta -> FillWeighted( pass, weight, ev -> muons.at(jmu).innereta); 
		muoninnerPhi -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerphi); 
		
		muonchi2   -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerchi2 );
		muondxy    -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerdxy);
		muondz     -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerdz);
		muonPixHit -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerpixelHits);
		muonLayHit -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerlayerHits);
		muonPixLay -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerpixelLayers);

		failingMuonPt  -> FillWeighted( !pass, weight, ev->muons.at(jmu).pt );
		failingMuonEta -> FillWeighted( !pass, weight, ev->muons.at(jmu).eta);
		failingMuonPhi -> FillWeighted( !pass, weight, ev->muons.at(jmu).phi);
		failingMuonEff -> FillWeighted( !pass, weight, 0.5                    );
	}
	else{
		muonPt       -> Fill( pass, ev -> muons.at(jmu).pt );
		muonEta      -> Fill( pass, ev -> muons.at(jmu).eta);
		muonPhi      -> Fill( pass, ev -> muons.at(jmu).phi);
		muonEff      -> Fill( pass, 0.5                    );
		muonDeltaR   -> Fill( pass, DeltaR);
		muonDeltaPhi -> Fill( pass, DeltaPhi);
		nvtx	     -> Fill( pass, ev -> nVtx);

		muoninnerPt  -> Fill( pass, ev -> muons.at(jmu).innerpt);
		muoninnerEta -> Fill( pass, ev -> muons.at(jmu).innereta); 
		muoninnerPhi -> Fill( pass, ev -> muons.at(jmu).innerphi); 
		
		muonchi2   -> Fill( pass, ev -> muons.at(jmu).innerchi2 );
		muondxy    -> Fill( pass, ev -> muons.at(jmu).innerdxy);
		muondz     -> Fill( pass, ev -> muons.at(jmu).innerdz);
		muonPixHit -> Fill( pass, ev -> muons.at(jmu).innerpixelHits);
		muonLayHit -> Fill( pass, ev -> muons.at(jmu).innerlayerHits);
		muonPixLay -> Fill( pass, ev -> muons.at(jmu).innerpixelLayers);

		failingMuonPt  -> Fill( !pass, ev->muons.at(jmu).pt );
		failingMuonEta -> Fill( !pass, ev->muons.at(jmu).eta);
		failingMuonPhi -> Fill( !pass, ev->muons.at(jmu).phi);
		failingMuonEff -> Fill( !pass, 0.5                    );
	}
      } // nmuons
    }
    
    /// NOW FOR DIMUONS:
    // both should pass the L1 filter! 
    if (firedL1(ev -> hlt.objects, theprobefilter)) { 
      for (int imu = 0; imu < nmuons; imu++){
	if (!selectMuon(ev -> muons.at(imu))) continue; 
	for (int jmu = imu+1; jmu < nmuons; jmu++){
	  if (!selectMuon(ev -> muons.at(jmu))) continue; 
	  bool pass = false;
	  if (matchMuonWithL3(ev->muons.at(jmu),ev->hltmuons) && 
	      matchMuonWithL3(ev->muons.at(imu),ev->hltmuons)) pass = true;
	  
	  double DeltaR = deltaR(ev->muons.at(imu).eta,ev->muons.at(imu).phi,ev->muons.at(jmu).eta,ev->muons.at(jmu).phi);
	  diMuonPt       -> Fill( pass, ev->muons.at(imu).pt , ev->muons.at(jmu).pt  ); 
	  diMuonLeadPt   -> Fill( pass, ev->muons.at(imu).pt ); 
	  diMuonTrailPt  -> Fill( pass, ev->muons.at(jmu).pt ); 
	  
	  if (ev->muons.at(imu).pt < ptcut1) continue;
	  if (ev->muons.at(jmu).pt < ptcut2) continue;
	  diMuonEta      -> Fill( pass, ev->muons.at(imu).eta, ev->muons.at(jmu).eta );
	  diMuonPhi      -> Fill( pass, ev->muons.at(imu).phi, ev->muons.at(jmu).phi );
	  diMuonDeltaR   -> Fill( pass, DeltaR);
	  diMuonEff      -> Fill( pass, 0.5);
	  
	  diMuonLeadEta  -> Fill( pass, ev->muons.at(imu).eta ); 
	  diMuonTrailEta -> Fill( pass, ev->muons.at(jmu).eta ); 
	  diMuonLeadPhi  -> Fill( pass, ev->muons.at(imu).phi ); 
	  diMuonTrailPhi -> Fill( pass, ev->muons.at(jmu).phi ); 
	}
      }
    }
  }  
  
  //Writing the histograms in a file.
  outfile           -> cd();
  tagMuonPt         -> Write();
  
  muonPt            -> Write();
  muonPtTurnOn      -> Write();
  muonEta           -> Write();
  muonPhi           -> Write();
  muonEff           -> Write();
  muonDeltaR        -> Write();
  muonDeltaPhi      -> Write();


  muonchi2     -> Write();
  muondxy      -> Write();
  muondz       -> Write();
  muonPixHit   -> Write();
  muonLayHit   -> Write();
  muonPixLay   -> Write();
  muoninnerPt  -> Write();
  muoninnerEta -> Write();
  muoninnerPhi -> Write();

  failingMuonPt   -> Write();
  failingMuonEta  -> Write();
  failingMuonPhi  -> Write();
  failingMuonEff  -> Write();

  ProbePt -> Write();

  PassingProbePt  -> Write();
  PassingProbeEta -> Write();
  PassingProbePhi -> Write();
  PassingProbeMll -> Write();
  
  FailingProbePt  -> Write();
  FailingProbeEta -> Write();
  FailingProbePhi -> Write();
  FailingProbeMll -> Write();
  
  /// Per-event (di-muon) efficiency.
  diMuonPt         -> Write();
  diMuonEta        -> Write();
  diMuonPhi        -> Write();
  diMuonEff        -> Write();
  diMuonDeltaR     -> Write();
  diMuonLeadPt     -> Write(); 
  diMuonLeadEta    -> Write();
  diMuonLeadPhi    -> Write();
  diMuonTrailPt    -> Write(); 
  diMuonTrailEta   -> Write();
  diMuonTrailPhi   -> Write();
  nvtx_event       -> Write();
  nvtx             -> Write();
  
  dimuon_mass      -> Write();
  tagiso           -> Write();
  
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
  if (signature == Sig::Prompt) ptcut = 26.;
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

bool matchMuonWithL3(MuonCand mu, std::vector<HLTMuonCand> L3cands){

  bool match = false;
  float minDR = 0.3;
  float theDR = 100;
  for ( std::vector<HLTMuonCand>::const_iterator it = L3cands.begin(); it != L3cands.end(); ++it ) { 
    theDR = deltaR(it -> eta, it -> phi, mu.eta, mu.phi); 
    if (theDR < minDR){ 
      minDR = theDR;
      match = true;
    }
  }
  return match;
}
std::string getProbeFilter(int signature, bool isMC){
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
