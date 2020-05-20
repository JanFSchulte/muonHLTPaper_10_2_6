

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

bool selectTagMuon  (MuonCand, TH1F* );
bool selectProbeMuon(MuonCand, MuonCand, TH1F* );
bool selectMuon     (MuonCand);
bool matchMuon      (MuonCand, std::vector<HLTObjCand>, std::string);
//HLTMuonCand  matchL3        (MuonCand, std::vector<HLTMuonCand>);
//L1MuonCand   matchL1        (MuonCand, std::vector<L1MuonCand>);
bool  matchMuonWithL3 (MuonCand, std::vector<HLTMuonCand>);
void printProgBar(int);

//std::string L1filter      =  "hltL1fL1sMu22Or25L1Filtered0::TEST"; 
//std::string L2filter      = "hltL2fL1sMu22or25L1f0L2Filtered10Q::TEST";
//std::string L3filter      =  "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q::TEST"; 
//std::string isofilterTag  = "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p09::HLT";
/// TAG-DEFINITION: 
std::string isofilterTag  = "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07::HLT";

/// for PROMPT-MUONS   (close-by and far-away) 
std::string L1filter      = "hltL1fL1sMu22or25L1Filtered0::TEST"; 
std::string L2filter      = "hltL2fL1sMu22or25L1f0L2Filtered10Q::TEST";
std::string L3filter      = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q::TEST"; 


double pt_bins[20]  = { 5, 7, 9, 12, 16,  20 ,  24 ,  27 ,   30,   35,   40,   45,   50,  60, 70 ,  90, 150, 250,500,1000 };
double dz_bins[11]  = { -15., -8., -6., -4., -2.,  0.,  2.,  4.,   6.,   8.,  15.};
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
0.63520,
0.63520,
0.55436,
1.3451,
0.63520,
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



// ******************************************
//                                          *
//                                          *
std::string hltname        = "HLT_IsoMu27_v10"; 
std::string thepassfilter  = L3filter;
std::string the2016filter  = "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07::HLT";
std::string theprobefilter = L1filter; 
float offlinePtCut         = 24.; 
//                                          *
//                                          *
// ******************************************

void readNtuplesPre_CascadeAndTkMu(TString inputfilename="files/files/",/*/eos/uscms/store/user/bmahakud/ProductionCasTest_v1/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ProductionCasTest_v1/181205_133338/0000/muonNtupleCorrCasTk.root",*/ std::string effmeasured="CascadeORTkMu_", bool isMC=false){

  ///eos/uscms/store/user/bmahakud/ProductionCasTest_v1/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ProductionCasTest_v1/181205_133338/0000/NtupleTMP.root
  ///eos/uscms/store/user/bmahakud/TestCascade_LPC_v3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/TestCascade_LPC_v3/181203_142416/0000/muonNtupleCasNew.root

  bool doingL1 = thepassfilter.find("L1fL1") != std::string::npos; 

  TFile* outfile = TFile::Open(Form("%s_efficiency_prefilter.root", effmeasured.c_str()),"RECREATE");
  std::cout << "output file: " << outfile -> GetName() << std::endl;

  //Create histograms  
  TH1F* dimuon_mass             = new TH1F("h_dimuon_mass"          ,"dimuon_mass"      , 1500,  0,  150 );
  TH2F* dimuon_mass_vs_pt_probePass       = new TH2F("h_dimuon_mass_vs_pt_probePass"    ,"dimuon mass vs pt" ,1500,  0,  150, 200, 0,  2000 );
  TH2F* dimuon_mass_vs_pt_probeFail       = new TH2F("h_dimuon_mass_vs_pt_probeFail"    ,"dimuon mass vs pt" ,1500,  0,  150, 200, 0, 2000 );


  TH1F* tagiso                  = new TH1F("h_tagiso"               ,"tagiso"           ,  100,  0,  1   );
  TH1F* tagMuonPt               = new TH1F("h_tagMuonPt"            ,"tagMuonPt"        ,  150,  0,  150 );
  TH1F* nvtx_event              = new TH1F("h_nvtx_event"           ,"nvtx_event"       ,   50,  0,   100 );



  TEfficiency* muonPt_barrel_cascade    = new TEfficiency("muonPt_barrel_cascade"   ,"muonPt_barrel_cascade"    ,   19,  pt_bins );
  TEfficiency* muonPt_endcap_cascade    = new TEfficiency("muonPt_endcap_cascade"   ,"muonPt_endcap_cascade"    ,   19,  pt_bins );
  TEfficiency* muonPt_cascade           = new TEfficiency("muonPt_cascade"          ,"muonPt_cascade"           ,   19,  pt_bins ); 
  TEfficiency* muonPtTurnOn_cascade     = new TEfficiency("muonPtTurnOn_cascade"    ,"muonPtTurnOn_cascade"     ,   19,  pt_bins ); 
  TEfficiency* muonEta_cascade          = new TEfficiency("muonEta_cascade"         ,"muonEta_cascade"          ,   15, eta_bins );
  TEfficiency* muonPhi_cascade          = new TEfficiency("muonPhi_cascade"         ,"muonPhi_cascade"          ,   20, -3.2, 3.2);
  TEfficiency* muonEff_cascade          = new TEfficiency("muonEff_cascade"         ,"muonEff_cascade"          ,    1,   0., 1.0);

  // GLOBAL QUANTITIES//
  TEfficiency* muonchi2_cascade         = new TEfficiency("muonchi2_cascade"        , "muonchi2_cascade"        ,  20,    0.,  7.);
  TEfficiency* muondxy_cascade          = new TEfficiency("muondxy_cascade"         , "muondxy_cascade"         ,  100,  -0.3, 0.3);
  TEfficiency* muondz_cascade           = new TEfficiency("muondz_cascade"          , "muondz_cascade"          ,   10,   dz_bins);
  TEfficiency* muonPixHit_cascade       = new TEfficiency("muonPixHit_cascade"      , "muonPixHit_cascade"      ,   20,  -0.5,19.5);
  TEfficiency* muonLayHit_cascade       = new TEfficiency("muonLayHit_cascade"      , "muonLayHit_cascade"      ,   16,   2.5,18.5);
  TEfficiency* muonPixLay_cascade       = new TEfficiency("muonPixLay_cascade"      , "muonPixLay_cascade"      ,    7,  -0.5, 6.5);
  TEfficiency* muoninnerPt_cascade      = new TEfficiency("muoninnerPt_cascade"     , "muoninnerPt_cascade"     ,   19,   pt_bins );
  TEfficiency* muoninnerEta_cascade     = new TEfficiency("muoninnerEta_cascade"    , "muoninnerEta_cascade"    ,   15,   eta_bins);
  TEfficiency* muoninnerPhi_cascade     = new TEfficiency("muoninnerPhi_cascade"    , "muoninnerPhi_cascade"    ,   20,  -3.2, 3.2);
  TEfficiency* muonValHits_cascade      = new TEfficiency("muonValHits_cascade"     , "muonValHits_cascade"     ,   20, -0.5, 19.5);

  TEfficiency* muonPt_barrel_tkmu    = new TEfficiency("muonPt_barrel_tkmu"   ,"muonPt_barrel_tkmu"    ,   19,  pt_bins );
  TEfficiency* muonPt_endcap_tkmu    = new TEfficiency("muonPt_endcap_tkmu"   ,"muonPt_endcap_tkmu"    ,   19,  pt_bins );
  TEfficiency* muonPt_tkmu           = new TEfficiency("muonPt_tkmu"          ,"muonPt_tkmu"           ,   19,  pt_bins ); 
  TEfficiency* muonPtTurnOn_tkmu     = new TEfficiency("muonPtTurnOn_tkmu"    ,"muonPtTurnOn_tkmu"     ,   19,  pt_bins ); 
  TEfficiency* muonEta_tkmu          = new TEfficiency("muonEta_tkmu"         ,"muonEta_tkmu"          ,   15, eta_bins );
  TEfficiency* muonPhi_tkmu          = new TEfficiency("muonPhi_tkmu"         ,"muonPhi_tkmu"          ,   20, -3.2, 3.2);
  TEfficiency* muonEff_tkmu          = new TEfficiency("muonEff_tkmu"         ,"muonEff_tkmu"          ,    1,   0., 1.0);

  // GLOBAL QUANTITIES//
  TEfficiency* muonchi2_tkmu         = new TEfficiency("muonchi2_tkmu"        , "muonchi2_tkmu"        ,  20,    0.,  7.);
  TEfficiency* muondxy_tkmu          = new TEfficiency("muondxy_tkmu"         , "muondxy_tkmu"         ,  100,  -0.3, 0.3);
  TEfficiency* muondz_tkmu           = new TEfficiency("muondz_tkmu"          , "muondz_tkmu"          ,   10,   dz_bins);
  TEfficiency* muonPixHit_tkmu       = new TEfficiency("muonPixHit_tkmu"      , "muonPixHit_tkmu"      ,   20,  -0.5,19.5);
  TEfficiency* muonLayHit_tkmu       = new TEfficiency("muonLayHit_tkmu"      , "muonLayHit_tkmu"      ,   16,   2.5,18.5);
  TEfficiency* muonPixLay_tkmu       = new TEfficiency("muonPixLay_tkmu"      , "muonPixLay_tkmu"      ,    7,  -0.5, 6.5);
  TEfficiency* muoninnerPt_tkmu      = new TEfficiency("muoninnerPt_tkmu"     , "muoninnerPt_tkmu"     ,   19,   pt_bins );
  TEfficiency* muoninnerEta_tkmu     = new TEfficiency("muoninnerEta_tkmu"    , "muoninnerEta_tkmu"    ,   15,   eta_bins);
  TEfficiency* muoninnerPhi_tkmu     = new TEfficiency("muoninnerPhi_tkmu"    , "muoninnerPhi_tkmu"    ,   20,  -3.2, 3.2);
  TEfficiency* muonValHits_tkmu      = new TEfficiency("muonValHits_tkmu"     , "muonValHits_tkmu"     ,   20, -0.5, 19.5);



 
  TEfficiency* muonPt_barrel    = new TEfficiency("muonPt_barrel"   ,"muonPt_barrel"    ,   19,  pt_bins );
  TEfficiency* muonPt_endcap    = new TEfficiency("muonPt_endcap"   ,"muonPt_endcap"    ,   19,  pt_bins );
  TEfficiency* muonPt           = new TEfficiency("muonPt"          ,"muonPt"           ,   19,  pt_bins ); 
  TEfficiency* muonPtTurnOn     = new TEfficiency("muonPtTurnOn"    ,"muonPtTurnOn"     ,   19,  pt_bins ); 
  TEfficiency* muonEta          = new TEfficiency("muonEta"         ,"muonEta"          ,   15, eta_bins );
  TEfficiency* muonPhi          = new TEfficiency("muonPhi"         ,"muonPhi"          ,   20, -3.2, 3.2);
  TEfficiency* muonEff          = new TEfficiency("muonEff"         ,"muonEff"          ,    1,   0., 1.0);

  // GLOBAL QUANTITIES//
  TEfficiency* muonchi2         = new TEfficiency("muonchi2"        , "muonchi2"        ,  20,    0.,  7.);
  TEfficiency* muondxy          = new TEfficiency("muondxy"         , "muondxy"         ,  100,  -0.3, 0.3);
  TEfficiency* muondz           = new TEfficiency("muondz"          , "muondz"          ,   10,   dz_bins);
  TEfficiency* muonPixHit       = new TEfficiency("muonPixHit"      , "muonPixHit"      ,   20,  -0.5,19.5);
  TEfficiency* muonLayHit       = new TEfficiency("muonLayHit"      , "muonLayHit"      ,   16,   2.5,18.5);
  TEfficiency* muonPixLay       = new TEfficiency("muonPixLay"      , "muonPixLay"      ,    7,  -0.5, 6.5);
  TEfficiency* muoninnerPt      = new TEfficiency("muoninnerPt"     , "muoninnerPt"     ,   19,   pt_bins );
  TEfficiency* muoninnerEta     = new TEfficiency("muoninnerEta"    , "muoninnerEta"    ,   15,   eta_bins);
  TEfficiency* muoninnerPhi     = new TEfficiency("muoninnerPhi"    , "muoninnerPhi"    ,   20,  -3.2, 3.2);
  TEfficiency* muonValHits      = new TEfficiency("muonValHits"     , "muonValHits"     ,   20, -0.5, 19.5);

  // BARREL//
  TEfficiency* muoninnerPt_barrel    = new TEfficiency("muoninnerPt_barrel"   , "muoninnerPt_barrel"   ,   19,   pt_bins );
  TEfficiency* muoninnerEta_barrel   = new TEfficiency("muoninnerEta_barrel"  , "muoninnerEta_barrel"  ,   15,   eta_bins);
  TEfficiency* muoninnerPhi_barrel   = new TEfficiency("muoninnerPhi_barrel"  , "muoninnerPhi_barrel"  ,   20,  -3.2, 3.2);
  TEfficiency* muonchi2_barrel       = new TEfficiency("muonchi2_barrel"      , "muonchi2_barrel"      ,  20,    0.,   7.);
  TEfficiency* muondxy_barrel        = new TEfficiency("muondxy_barrel"       , "muondxy_barrel"       ,  100, -0.3,  0.3);
  TEfficiency* muondz_barrel         = new TEfficiency("muondz_barrel"        , "muondz_barrel"        ,   10,  dz_bins);
  TEfficiency* muonPixHit_barrel     = new TEfficiency("muonPixHit_barrel"    , "muonPixHit_barrel"    ,   20,  -0.5, 19.5);
  TEfficiency* muonLayHit_barrel     = new TEfficiency("muonLayHit_barrel"    , "muonLayHit_barrel"    ,   16,   2.5, 18.5);
  TEfficiency* muonPixLay_barrel     = new TEfficiency("muonPixLay_barrel"    , "muonPixLay_barrel"    ,    7,  -0.5,  6.5);
  TEfficiency* muonValHits_barrel    = new TEfficiency("muonValHits_barrel"   , "muonValHits_barrel"   ,   20, -0.5, 19.5);

  // INTERMEDIATE REGION//
  TEfficiency* muonPt_int       = new TEfficiency("muonPt_int"       , "muonPt_int"        ,   19,  pt_bins );
  TEfficiency* nvtx_int         = new TEfficiency("nvtx_int"         , "nvtx_int"          ,   60,    0,   60);
  TEfficiency* muoninnerPt_int  = new TEfficiency("muoninnerPt_int"  , "muoninnerPt_int"   ,   19,   pt_bins );
  TEfficiency* muoninnerEta_int = new TEfficiency("muoninnerEta_int" , "muoninnerEta_int"  ,   15,   eta_bins);
  TEfficiency* muoninnerPhi_int = new TEfficiency("muoninnerPhi_int" , "muoninnerPhi_int"  ,   20,  -3.2, 3.2);
  TEfficiency* muonchi2_int     = new TEfficiency("muonchi2_int"     , "muonchi2_int"      ,  20,   0.,  15.);
  TEfficiency* muondxy_int      = new TEfficiency("muondxy_int"      , "muondxy_int"       ,  100, -0.3,  0.3);
  TEfficiency* muondz_int       = new TEfficiency("muondz_int"       , "muondz_int"        ,   10, dz_bins);
  TEfficiency* muonPixHit_int   = new TEfficiency("muonPixHit_int"   , "muonPixHit_int"    ,   20, -0.5, 19.5);
  TEfficiency* muonLayHit_int   = new TEfficiency("muonLayHit_int"   , "muonLayHit_int"    ,   16,  2.5, 18.5);
  TEfficiency* muonPixLay_int   = new TEfficiency("muonPixLay_int"   , "muonPixLay_int"    ,    7, -0.5,  6.5);
  TEfficiency* muonValHits_int  = new TEfficiency("muonValHits_int"  , "muonValHits_int"   ,   20, -0.5, 19.5);

  // ENDCAP//
  TEfficiency* muoninnerPt_endcap  = new TEfficiency("muoninnerPt_endcap"  , "muoninnerPt_endcap"   ,   19,   pt_bins );
  TEfficiency* muoninnerEta_endcap = new TEfficiency("muoninnerEta_endcap" , "muoninnerEta_endcap"  ,   15,   eta_bins);
  TEfficiency* muoninnerPhi_endcap = new TEfficiency("muoninnerPhiendcap"  , "muoninnerPhi_endcap"  ,   20,  -3.2, 3.2);
  TEfficiency* muonchi2_endcap     = new TEfficiency("muonchi2_endcap"     , "muonchi2_endcap"      ,  20,   0.,   7.);
  TEfficiency* muondxy_endcap      = new TEfficiency("muondxy_endcap"      , "muondxy_endcap"       ,  100, -0.3,  0.3);
  TEfficiency* muondz_endcap       = new TEfficiency("muondz_endcap"       , "muondz_endcap"        ,   10, dz_bins);
  TEfficiency* muonPixHit_endcap   = new TEfficiency("muonPixHit_endcap"   , "muonPixHit_endcap"    ,   20, -0.5, 19.5);
  TEfficiency* muonLayHit_endcap   = new TEfficiency("muonLayHit_endcap"   , "muonLayHit_endcap"    ,   16,  2.5, 18.5);
  TEfficiency* muonPixLay_endcap   = new TEfficiency("muonPixLay_endcap"   , "muonPixLay_endcap"    ,    7, -0.5,  6.5);
  TEfficiency* muonValHits_endcap  = new TEfficiency("muonValHits_endcap"  , "muonValHits_endcap"   ,   20, -0.5, 19.5);

  TEfficiency* muonOver16Pt     = new TEfficiency("muonOver16Pt"    ,"muonOver16Pt"     ,   19,  pt_bins ); 
  TEfficiency* muonOver16Eta    = new TEfficiency("muonOver16Eta"   ,"muonOver16Eta"    ,   15, eta_bins );
  TEfficiency* muonOver16Phi    = new TEfficiency("muonOver16Phi"   ,"muonOver16Phi"    ,   20, -3.2, 3.2);
  TEfficiency* muonOver16Eff    = new TEfficiency("muonOver16Eff"   ,"muonOver16Eff"    ,   1 ,   0., 1.0);
  
  TEfficiency* failingMuonPt    = new TEfficiency("failingMuonPt"   ,"failingMuonPt"    ,   19,  pt_bins ); 
  TEfficiency* failingMuonEta   = new TEfficiency("failingMuonEta"  ,"failingMuonEta"   ,   15, eta_bins );
  TEfficiency* failingMuonPhi   = new TEfficiency("failingMuonPhi"  ,"failingMuonPhi"   ,   20, -3.2, 3.2);
  TEfficiency* failingMuonEff   = new TEfficiency("failingMuonEff"  ,"failingMuonEff"   ,   1 ,   0., 1.0);

  TH1F* ProbePt                 = new TH1F("h_ProbePt"              ,"ProbeMuonPt"      ,10000, 0, 1000 );
  TH1F* PassingProbePt          = new TH1F("h_PassingProbePt"       ,"PassingMuonPt"    ,  19,  pt_bins );
  TH1F* PassingProbeEta         = new TH1F("h_PassingProbeEta"      ,"PassingMuonEta"   ,  15, eta_bins );
  TH1F* PassingProbePhi         = new TH1F("h_PassingProbePhi"      ,"PassingMuonPhi"   ,  20, -3.2, 3.2);
  TH1F* PassingProbeMll         = new TH1F("h_PassingProbeMll"      ,"PassingMuonMll"   ,  20,  86., 96.);

  TH1F* FailingProbePt          = new TH1F("h_FailingProbePt"       ,"FailingMuonPt"    ,  19,  pt_bins );
  TH1F* FailingProbeEta         = new TH1F("h_FailingProbeEta"      ,"FailingMuonEta"   ,  15, eta_bins );
  TH1F* FailingProbePhi         = new TH1F("h_FailingProbePhi"      ,"FailingMuonPhi"   ,  20, -3.2, 3.2);
  TH1F* FailingProbeMll         = new TH1F("h_FailingProbeMll"      ,"FailingMuonMll"   ,  20,  86., 96.);

  // for dimuon eff: 
  TEfficiency* diMuonPt         = new TEfficiency("diMuonPt"        ,"diMuonPt"         ,   19,  pt_bins ); 
  TEfficiency* diMuonEta        = new TEfficiency("diMuonEta"       ,"diMuonEta"        ,   15, eta_bins );
  TEfficiency* diMuonPhi        = new TEfficiency("diMuonPhi"       ,"diMuonPhi"        ,   20, -3.2, 3.2);

  // for noFilter matching: 
  TEfficiency* hltmuonPt        = new TEfficiency("hltmuonPt"       ,"hltmuonPt"        ,   19,  pt_bins ); 
  TEfficiency* hltmuonEta       = new TEfficiency("hltmuonEta"      ,"hltmuonEta"       ,   15, eta_bins );
  TEfficiency* hltmuonPhi       = new TEfficiency("hltmuonPhi"      ,"hltmuonPhi"       ,   20, -3.2, 3.2);
  TEfficiency* hltmuonEff       = new TEfficiency("hltmuonEff"      ,"hltmuonEff"       ,   1 ,   0.,  1.);

  TEfficiency* hltmuonFromL2Pt  = new TEfficiency("hltmuonFromL2Pt"  ,"hltmuonFromL2Pt"  ,   19,  pt_bins ); 
  TEfficiency* hltmuonFromL2Eta = new TEfficiency("hltmuonFromL2Eta" ,"hltmuonFromL2Eta" ,   15, eta_bins );
  TEfficiency* hltmuonFromL2Phi = new TEfficiency("hltmuonFromL2Phi" ,"hltmuonFromL2Phi" ,   20, -3.2, 3.2);
  TEfficiency* hltmuonFromL2Eff = new TEfficiency("hltmuonFromL2Eff" ,"hltmuonFromL2Eff" ,   1 ,   0.,  1.);

  TEfficiency* nvtx             = new TEfficiency("nvtx"             ,"nvtx"             ,   50,    0,  100);
  TEfficiency* nvtx_cascade     = new TEfficiency("nvtx_cascade"     ,"nvtx_cascade"     ,   50,    0,  100);
  TEfficiency* nvtx_tkmu        = new TEfficiency("nvtx_tkmu"        ,"nvtx_tkmu"        ,   50,    0,  100);
  TEfficiency* nvtx_barrel      = new TEfficiency("nvtx_barrel"      ,"nvtx_barrel"      ,   60,    0,  60);
  TEfficiency* nvtx_endcap      = new TEfficiency("nvtx_endcap"      ,"nvtx_endcap"      ,   60,    0,  60);
   
  double offlineiso04 = 100;
  TChain *tree = new TChain("muonNtuples/muonTree");
  tree->Add(inputfilename); 


//  TFile* inputfile = TFile::Open(inputfilename, "READ");
//  std::cout << "input file: " << inputfile -> GetName() << std::endl;

//  TTree *tree = (TTree*) inputfile -> Get("muonNtuples/muonTree"); 
  
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
  //int nentries = 1000000;
  std::cout << "Number of entries = " << nentries << std::endl;

  bool flagfile = false;
  for (Int_t eventNo=0; eventNo < nentries; eventNo++) 
    {
      Int_t IgetEvent   = tree   -> GetEvent(eventNo); 
      printProgBar((int)(eventNo*100./nentries));
      double weight = 1;
      if (isMC) weight = weights[ev -> nVtx];
 
      unsigned int nmuons = ev->muons.size();
      if (nmuons < 2) continue;
      unsigned int nhltmuons = ev->hltmuons.size(); 
    
      //    if (!ev-> hltTag.find(hltname)) continue;
      //nvtx_event-> Fill( ev -> nVtx   ); 

      for (int imu = 0; imu < nmuons; imu++){ 
	// select the tag muon        
	if (! selectTagMuon(ev -> muons.at(imu), tagiso)) continue;
	if (! matchMuon(ev -> muons.at(imu), ev -> hltTag.objects, isofilterTag)) continue;
	tagMuonPt -> Fill ( ev -> muons.at(imu).pt) ; 

	for (int jmu = 0; jmu < nmuons; jmu++){
	  bool pass   = false;
	  bool matchWith2016 = false;

	  // select the probe muon
	  if (!selectProbeMuon(ev -> muons.at(jmu), ev -> muons.at(imu), dimuon_mass)) continue;
	   //cout<< "probe muon"<< endl;
	  if (!doingL1 && !(matchMuon(ev -> muons.at(jmu), ev -> hlt.objects, theprobefilter))) continue;

	  //PRE FILTER
          bool passhlt=false;
          bool passTkMu=false;
          if(matchMuonWithL3(ev->muons.at(jmu),ev->hltmuons))passhlt=true;
          if(matchMuonWithL3(ev->muons.at(jmu),ev->tkmuons))passTkMu=true;

	  //if ( matchMuonWithL3(ev->muons.at(jmu),ev->hltmuons) || matchMuonWithL3(ev->muons.at(jmu),ev->tkmuons)) pass=true;
          if(passTkMu  || passhlt)pass=true;




	  //	  if ( matchMuonWithL3(ev->muons.at(jmu),ev->hltmuons)) pass=true;

	  muonPtTurnOn -> Fill( pass, ev -> muons.at(jmu).pt ); 
	  ProbePt->Fill(ev -> muons.at(jmu).pt); 
	  // now require pT cut: 
	  if (ev -> muons.at(jmu).pt < offlinePtCut) continue; 

	  TLorentzVector mu1, mu2;
	  mu1.SetPtEtaPhiM (ev->muons.at(imu).pt,ev->muons.at(imu).eta,ev->muons.at(imu).phi, muonmass); 
	  mu2.SetPtEtaPhiM (ev->muons.at(jmu).pt,ev->muons.at(jmu).eta,ev->muons.at(jmu).phi, muonmass);
	  double mumumass = (mu1 + mu2).M();

         if (passTkMu) dimuon_mass_vs_pt_probePass->Fill(mumumass, ev->muons.at(jmu).pt);
	 else dimuon_mass_vs_pt_probeFail->Fill(mumumass,ev->muons.at(jmu).pt);

	 if(isMC){
	  muonPt       -> FillWeighted( pass, weight, ev -> muons.at(jmu).pt );
	  muonEta      -> FillWeighted( pass, weight, ev -> muons.at(jmu).eta);
	  muonPhi      -> FillWeighted( pass, weight, ev -> muons.at(jmu).phi);
	  muonEff      -> FillWeighted( pass, weight, 0.5                    );

	  muoninnerPt  -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerpt);
	  muoninnerEta -> FillWeighted( pass, weight, ev -> muons.at(jmu).innereta); 
	  muoninnerPhi -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerphi); 

	  muonchi2   -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerchi2 );
	  muondxy    -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerdxy);
	  muondz     -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerdz);
	  muonPixHit -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerpixelHits);
	  muonLayHit -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerlayerHits);
	  muonPixLay -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerpixelLayers);
	  muonValHits-> FillWeighted( pass, weight, ev -> muons.at(jmu).innervalidHits);

	  //if (passhlt){

		  muonPt_cascade       -> FillWeighted( passhlt, weight, ev -> muons.at(jmu).pt );
		  muonEta_cascade      -> FillWeighted( passhlt, weight, ev -> muons.at(jmu).eta);
		  muonPhi_cascade      -> FillWeighted( passhlt, weight, ev -> muons.at(jmu).phi);
		  muonEff_cascade      -> FillWeighted( passhlt, weight, 0.5                    );

		  muoninnerPt_cascade  -> FillWeighted( passhlt, weight, ev -> muons.at(jmu).innerpt);
		  muoninnerEta_cascade -> FillWeighted( passhlt, weight, ev -> muons.at(jmu).innereta); 
		  muoninnerPhi_cascade -> FillWeighted( passhlt, weight, ev -> muons.at(jmu).innerphi); 

		  muonchi2_cascade   -> FillWeighted( passhlt, weight, ev -> muons.at(jmu).innerchi2 );
		  muondxy_cascade    -> FillWeighted( passhlt, weight, ev -> muons.at(jmu).innerdxy);
		  muondz_cascade     -> FillWeighted( passhlt, weight, ev -> muons.at(jmu).innerdz);
		  muonPixHit_cascade -> FillWeighted( passhlt, weight, ev -> muons.at(jmu).innerpixelHits);
		  muonLayHit_cascade -> FillWeighted( passhlt, weight, ev -> muons.at(jmu).innerlayerHits);
		  muonPixLay_cascade -> FillWeighted( passhlt, weight, ev -> muons.at(jmu).innerpixelLayers);
		  muonValHits_cascade-> FillWeighted( passhlt, weight, ev -> muons.at(jmu).innervalidHits);



	  //}
	  //if (passTkMu){


		  muonPt_tkmu       -> FillWeighted( passTkMu, weight, ev -> muons.at(jmu).pt );
		  muonEta_tkmu      -> FillWeighted( passTkMu, weight, ev -> muons.at(jmu).eta);
		  muonPhi_tkmu      -> FillWeighted( passTkMu, weight, ev -> muons.at(jmu).phi);
		  muonEff_tkmu      -> FillWeighted( passTkMu, weight, 0.5                    );

		  muoninnerPt_tkmu  -> FillWeighted( passTkMu, weight, ev -> muons.at(jmu).innerpt);
		  muoninnerEta_tkmu -> FillWeighted( passTkMu, weight, ev -> muons.at(jmu).innereta); 
		  muoninnerPhi_tkmu -> FillWeighted( passTkMu, weight, ev -> muons.at(jmu).innerphi); 

		  muonchi2_tkmu   -> FillWeighted( passTkMu, weight, ev -> muons.at(jmu).innerchi2 );
		  muondxy_tkmu    -> FillWeighted( passTkMu, weight, ev -> muons.at(jmu).innerdxy);
		  muondz_tkmu     -> FillWeighted( passTkMu, weight, ev -> muons.at(jmu).innerdz);
		  muonPixHit_tkmu -> FillWeighted( passTkMu, weight, ev -> muons.at(jmu).innerpixelHits);
		  muonLayHit_tkmu -> FillWeighted( passTkMu, weight, ev -> muons.at(jmu).innerlayerHits);
		  muonPixLay_tkmu -> FillWeighted( passTkMu, weight, ev -> muons.at(jmu).innerpixelLayers);
		  muonValHits_tkmu-> FillWeighted( passTkMu, weight, ev -> muons.at(jmu).innervalidHits);


	  nvtx           -> FillWeighted( pass, weight, ev -> nVtx             );
	  nvtx_cascade   -> FillWeighted( passhlt, weight, ev -> nVtx             );
	  nvtx_tkmu      -> FillWeighted( passTkMu, weight, ev -> nVtx             );

	  // BARREL//
	  if (fabs(ev -> muons.at(jmu).eta) <= 0.9){
	    muonPt_barrel       -> FillWeighted( pass, weight, ev -> muons.at(jmu).pt );
	    muoninnerPt_barrel  -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerpt); 
	    muoninnerEta_barrel -> FillWeighted( pass, weight, ev -> muons.at(jmu).innereta); 
	    muoninnerPhi_barrel -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerphi); 
	    nvtx_barrel         -> FillWeighted( pass, weight, ev -> nVtx );
	    muonchi2_barrel     -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerchi2 );
	    muondxy_barrel      -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerdxy);
	    muondz_barrel       -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerdz); 
	    muonPixHit_barrel   -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerpixelHits); 
	    muonLayHit_barrel   -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerlayerHits); 
	    muonPixLay_barrel   -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerpixelLayers);  
	    muonValHits_barrel  -> FillWeighted( pass, weight, ev -> muons.at(jmu).innervalidHits);

	  }

	  // INTERMEDIATE REGION//
	  if (fabs(ev -> muons.at(jmu).eta)>0.9 && fabs(ev -> muons.at(jmu).eta)<1.6){
	    muonPt_int       -> FillWeighted( pass, weight, ev -> muons.at(jmu).pt );
	    muoninnerPt_int  -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerpt); 
	    muoninnerEta_int -> FillWeighted( pass, weight, ev -> muons.at(jmu).innereta); 
	    muoninnerPhi_int -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerphi); 
	    nvtx_int         -> FillWeighted( pass, weight, ev -> nVtx );
	    muonchi2_int     -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerchi2 );
	    muondxy_int      -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerdxy);
	    muondz_int       -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerdz); 
	    muonPixHit_int   -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerpixelHits); 
	    muonLayHit_int   -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerlayerHits); 
	    muonPixLay_int   -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerpixelLayers);  
	    muonValHits_int  -> FillWeighted( pass, weight, ev -> muons.at(jmu).innervalidHits);
	  }
	  // ENDCAP//
	  if ( fabs(ev -> muons.at(jmu).eta)>=1.6){
	    muonPt_endcap       -> FillWeighted( pass, weight, ev -> muons.at(jmu).pt );
	    muoninnerPt_endcap  -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerpt); 
	    muoninnerEta_endcap -> FillWeighted( pass, weight, ev -> muons.at(jmu).innereta); 
	    muoninnerPhi_endcap -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerphi); 
	    nvtx_endcap         -> FillWeighted( pass, weight, ev -> nVtx );
	    muonchi2_endcap     -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerchi2 );
	    muondxy_endcap      -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerdxy);
	    muondz_endcap       -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerdz); 
	    muonPixHit_endcap   -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerpixelHits); 
	    muonLayHit_endcap   -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerlayerHits); 
	    muonPixLay_endcap   -> FillWeighted( pass, weight, ev -> muons.at(jmu).innerpixelLayers); 
	    muonValHits_endcap  -> FillWeighted( pass, weight, ev -> muons.at(jmu).innervalidHits);
	   }
	 }
	else{
	  muonPt       -> Fill( pass, ev -> muons.at(jmu).pt );
	  muonEta      -> Fill( pass, ev -> muons.at(jmu).eta);
	  muonPhi      -> Fill( pass, ev -> muons.at(jmu).phi);
	  muonEff      -> Fill( pass, 0.5                    );

	  muoninnerPt  -> Fill( pass, ev -> muons.at(jmu).innerpt);
	  muoninnerEta -> Fill( pass, ev -> muons.at(jmu).innereta); 
	  muoninnerPhi -> Fill( pass, ev -> muons.at(jmu).innerphi); 

	  muonchi2   -> Fill( pass, ev -> muons.at(jmu).innerchi2 );
	  muondxy    -> Fill( pass, ev -> muons.at(jmu).innerdxy);
	  muondz     -> Fill( pass, ev -> muons.at(jmu).innerdz);
	  muonPixHit -> Fill( pass, ev -> muons.at(jmu).innerpixelHits);
	  muonLayHit -> Fill( pass, ev -> muons.at(jmu).innerlayerHits);
	  muonPixLay -> Fill( pass, ev -> muons.at(jmu).innerpixelLayers);
	  muonValHits-> Fill( pass, ev -> muons.at(jmu).innervalidHits);

	  //if (passhlt){

		  muonPt_cascade       -> Fill( passhlt, ev -> muons.at(jmu).pt );
		  muonEta_cascade      -> Fill( passhlt, ev -> muons.at(jmu).eta);
		  muonPhi_cascade      -> Fill( passhlt, ev -> muons.at(jmu).phi);
		  muonEff_cascade      -> Fill( passhlt, 0.5                    );

		  muoninnerPt_cascade  -> Fill( passhlt, ev -> muons.at(jmu).innerpt);
		  muoninnerEta_cascade -> Fill( passhlt, ev -> muons.at(jmu).innereta); 
		  muoninnerPhi_cascade -> Fill( passhlt, ev -> muons.at(jmu).innerphi); 

		  muonchi2_cascade   -> Fill( passhlt, ev -> muons.at(jmu).innerchi2 );
		  muondxy_cascade    -> Fill( passhlt, ev -> muons.at(jmu).innerdxy);
		  muondz_cascade     -> Fill( passhlt, ev -> muons.at(jmu).innerdz);
		  muonPixHit_cascade -> Fill( passhlt, ev -> muons.at(jmu).innerpixelHits);
		  muonLayHit_cascade -> Fill( passhlt, ev -> muons.at(jmu).innerlayerHits);
		  muonPixLay_cascade -> Fill( passhlt, ev -> muons.at(jmu).innerpixelLayers);
		  muonValHits_cascade-> Fill( passhlt, ev -> muons.at(jmu).innervalidHits);



	  //}
	  //if (passTkMu){


		  muonPt_tkmu       -> Fill( passTkMu, ev -> muons.at(jmu).pt );
		  muonEta_tkmu      -> Fill( passTkMu, ev -> muons.at(jmu).eta);
		  muonPhi_tkmu      -> Fill( passTkMu, ev -> muons.at(jmu).phi);
		  muonEff_tkmu      -> Fill( passTkMu, 0.5                    );

		  muoninnerPt_tkmu  -> Fill( passTkMu, ev -> muons.at(jmu).innerpt);
		  muoninnerEta_tkmu -> Fill( passTkMu, ev -> muons.at(jmu).innereta); 
		  muoninnerPhi_tkmu -> Fill( passTkMu, ev -> muons.at(jmu).innerphi); 

		  muonchi2_tkmu   -> Fill( passTkMu, ev -> muons.at(jmu).innerchi2 );
		  muondxy_tkmu    -> Fill( passTkMu, ev -> muons.at(jmu).innerdxy);
		  muondz_tkmu     -> Fill( passTkMu, ev -> muons.at(jmu).innerdz);
		  muonPixHit_tkmu -> Fill( passTkMu, ev -> muons.at(jmu).innerpixelHits);
		  muonLayHit_tkmu -> Fill( passTkMu, ev -> muons.at(jmu).innerlayerHits);
		  muonPixLay_tkmu -> Fill( passTkMu, ev -> muons.at(jmu).innerpixelLayers);
		  muonValHits_tkmu-> Fill( passTkMu, ev -> muons.at(jmu).innervalidHits);


	  nvtx           -> Fill( pass, ev -> nVtx             );
	  nvtx_cascade   -> Fill( passhlt, ev -> nVtx             );
	  nvtx_tkmu      -> Fill( passTkMu, ev -> nVtx             );

	  // BARREL//
	  if (fabs(ev -> muons.at(jmu).eta) <= 0.9){
	    muonPt_barrel       -> Fill( pass, ev -> muons.at(jmu).pt );
	    muoninnerPt_barrel  -> Fill( pass, ev -> muons.at(jmu).innerpt); 
	    muoninnerEta_barrel -> Fill( pass, ev -> muons.at(jmu).innereta); 
	    muoninnerPhi_barrel -> Fill( pass, ev -> muons.at(jmu).innerphi); 
	    nvtx_barrel         -> Fill( pass, ev -> nVtx );
	    muonchi2_barrel     -> Fill( pass, ev -> muons.at(jmu).innerchi2 );
	    muondxy_barrel      -> Fill( pass, ev -> muons.at(jmu).innerdxy);
	    muondz_barrel       -> Fill( pass, ev -> muons.at(jmu).innerdz); 
	    muonPixHit_barrel   -> Fill( pass, ev -> muons.at(jmu).innerpixelHits); 
	    muonLayHit_barrel   -> Fill( pass, ev -> muons.at(jmu).innerlayerHits); 
	    muonPixLay_barrel   -> Fill( pass, ev -> muons.at(jmu).innerpixelLayers);  
	    muonValHits_barrel  -> Fill( pass, ev -> muons.at(jmu).innervalidHits);

	  }

	  // INTERMEDIATE REGION//
	  if (fabs(ev -> muons.at(jmu).eta)>0.9 && fabs(ev -> muons.at(jmu).eta)<1.6){
	    muonPt_int       -> Fill( pass, ev -> muons.at(jmu).pt );
	    muoninnerPt_int  -> Fill( pass, ev -> muons.at(jmu).innerpt); 
	    muoninnerEta_int -> Fill( pass, ev -> muons.at(jmu).innereta); 
	    muoninnerPhi_int -> Fill( pass, ev -> muons.at(jmu).innerphi); 
	    nvtx_int         -> Fill( pass, ev -> nVtx );
	    muonchi2_int     -> Fill( pass, ev -> muons.at(jmu).innerchi2 );
	    muondxy_int      -> Fill( pass, ev -> muons.at(jmu).innerdxy);
	    muondz_int       -> Fill( pass, ev -> muons.at(jmu).innerdz); 
	    muonPixHit_int   -> Fill( pass, ev -> muons.at(jmu).innerpixelHits); 
	    muonLayHit_int   -> Fill( pass, ev -> muons.at(jmu).innerlayerHits); 
	    muonPixLay_int   -> Fill( pass, ev -> muons.at(jmu).innerpixelLayers);  
	    muonValHits_int  -> Fill( pass, ev -> muons.at(jmu).innervalidHits);
	  }
	  // ENDCAP//
	  if ( fabs(ev -> muons.at(jmu).eta)>=1.6){
	    muonPt_endcap       -> Fill( pass, ev -> muons.at(jmu).pt );
	    muoninnerPt_endcap  -> Fill( pass, ev -> muons.at(jmu).innerpt); 
	    muoninnerEta_endcap -> Fill( pass, ev -> muons.at(jmu).innereta); 
	    muoninnerPhi_endcap -> Fill( pass, ev -> muons.at(jmu).innerphi); 
	    nvtx_endcap         -> Fill( pass, ev -> nVtx );
	    muonchi2_endcap     -> Fill( pass, ev -> muons.at(jmu).innerchi2 );
	    muondxy_endcap      -> Fill( pass, ev -> muons.at(jmu).innerdxy);
	    muondz_endcap       -> Fill( pass, ev -> muons.at(jmu).innerdz); 
	    muonPixHit_endcap   -> Fill( pass, ev -> muons.at(jmu).innerpixelHits); 
	    muonLayHit_endcap   -> Fill( pass, ev -> muons.at(jmu).innerlayerHits); 
	    muonPixLay_endcap   -> Fill( pass, ev -> muons.at(jmu).innerpixelLayers); 
	    muonValHits_endcap  -> Fill( pass, ev -> muons.at(jmu).innervalidHits);
	   }

	}
	} // nmuons
      }
    }  

 
  //Writing the histograms in a file.
  outfile           -> cd();
  tagMuonPt         -> Write();
 
  nvtx->Write();
  nvtx_cascade->Write();
  nvtx_tkmu->Write();

  dimuon_mass_vs_pt_probePass -> Write();
  dimuon_mass_vs_pt_probeFail -> Write();
 
  muonPt            -> Write();
  muonPtTurnOn      -> Write();
  muonEta           -> Write();
  muonPhi           -> Write();
  muonEff           -> Write();

  muonchi2     -> Write();
  muondxy      -> Write();
  muondz       -> Write();
  muonPixHit   -> Write();
  muonLayHit   -> Write();
  muonPixLay   -> Write();
  muoninnerPt  -> Write();
  muoninnerEta -> Write();
  muoninnerPhi -> Write();
  muonValHits  -> Write();

  muonPt_cascade            -> Write();
  muonPtTurnOn_cascade      -> Write();
  muonEta_cascade           -> Write();
  muonPhi_cascade           -> Write();
  muonEff_cascade           -> Write();

  muonchi2_cascade     -> Write();
  muondxy_cascade      -> Write();
  muondz_cascade       -> Write();
  muonPixHit_cascade   -> Write();
  muonLayHit_cascade   -> Write();
  muonPixLay_cascade   -> Write();
  muoninnerPt_cascade  -> Write();
  muoninnerEta_cascade -> Write();
  muoninnerPhi_cascade -> Write();
  muonValHits_cascade  -> Write();

  muonPt_tkmu            -> Write();
  muonPtTurnOn_tkmu      -> Write();
  muonEta_tkmu           -> Write();
  muonPhi_tkmu           -> Write();
  muonEff_tkmu           -> Write();

  muonchi2_tkmu     -> Write();
  muondxy_tkmu      -> Write();
  muondz_tkmu       -> Write();
  muonPixHit_tkmu   -> Write();
  muonLayHit_tkmu   -> Write();
  muonPixLay_tkmu   -> Write();
  muoninnerPt_tkmu  -> Write();
  muoninnerEta_tkmu -> Write();
  muoninnerPhi_tkmu -> Write();
  muonValHits_tkmu  -> Write();



  muonPt_barrel       -> Write();
  muoninnerPt_barrel  -> Write();
  muoninnerEta_barrel -> Write();
  muoninnerPhi_barrel -> Write();
  nvtx_barrel         -> Write();
  muonchi2_barrel     -> Write();
  muondxy_barrel      -> Write();
  muondz_barrel       -> Write();
  muonPixHit_barrel   -> Write();
  muonLayHit_barrel   -> Write();
  muonPixLay_barrel   -> Write();
  muonValHits_barrel  -> Write();

  muonPt_int       -> Write();
  muoninnerPt_int  -> Write();
  muoninnerEta_int -> Write();
  muoninnerPhi_int -> Write();
  nvtx_int         -> Write();
  muonchi2_int     -> Write();
  muondxy_int      -> Write();
  muondz_int       -> Write();
  muonPixHit_int   -> Write();
  muonLayHit_int   -> Write();
  muonPixLay_int   -> Write();
  muonValHits_int  -> Write();

  muonPt_endcap       -> Write();
  muoninnerPt_endcap  -> Write();
  muoninnerEta_endcap -> Write();
  muoninnerPhi_endcap -> Write();
  nvtx_endcap         -> Write();
  muonchi2_endcap     -> Write();
  muondxy_endcap      -> Write();
  muondz_endcap       -> Write();
  muonPixHit_endcap   -> Write();
  muonLayHit_endcap   -> Write();
  muonPixLay_endcap   -> Write();
  muonValHits_endcap  -> Write();

  failingMuonPt   -> Write();
  failingMuonEta  -> Write();
  failingMuonPhi  -> Write();
  failingMuonEff  -> Write();

  ProbePt->Write();

  PassingProbePt  -> Write();
  PassingProbeEta -> Write();
  PassingProbePhi -> Write();
  PassingProbeMll -> Write();
  
  FailingProbePt  -> Write();
  FailingProbeEta -> Write();
  FailingProbePhi -> Write();
  FailingProbeMll -> Write();

  hltmuonPt       -> Write();
  hltmuonEta      -> Write();
  hltmuonPhi      -> Write();
  hltmuonEff      -> Write();

  hltmuonFromL2Pt  -> Write();
  hltmuonFromL2Eta -> Write();
  hltmuonFromL2Phi -> Write();
  hltmuonFromL2Eff -> Write();

  diMuonPt         -> Write();
  diMuonEta        -> Write();
  diMuonPhi        -> Write();

  nvtx_event       -> Write();
//  nvtx             -> Write();
  
  dimuon_mass      -> Write();
  tagiso           -> Write();
  
  outfile          -> Close();  
  
  return;
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
  
  if (!( mu.pt         > 26  )) return false; 
  if (!( fabs(mu.eta)  < 2.4 )) return false; 
  if (!( mu.isTight    == 1  )) return false; 
  
  //add isolation cut
  float offlineiso04 = mu.chargedDep_dR04 + std::max(0., mu.photonDep_dR04 + mu.neutralDep_dR04 - 0.5*mu.puPt_dR04);
  offlineiso04       = offlineiso04 / mu.pt;
  tagh -> Fill(offlineiso04);
  if (offlineiso04   > offlineIsoCut) return false; 

  return true;
}


bool selectMuon(MuonCand mu){  
  if (!( mu.pt         > offlinePtCut  )) return false; 
  if (!( fabs(mu.eta)  < 2.4 )) return false;
  if (!( mu.isLoose    == 1  )) return false; 
  return true;
}


//select the probe muon
bool selectProbeMuon(MuonCand mu, MuonCand tagMu, TH1F* dimuon_mass){
  
  if (mu.pt == tagMu.pt  && 
      mu.pt == tagMu.eta &&
      mu.pt == tagMu.phi ) 
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
//  if (! (mumumass > 81. && mumumass < 101. )) return false;
  
  return true;
}

//*********************************************************************************************************************
//HLTMuonCand matchL3(MuonCand mu, std::vector<HLTMuonCand> L3cands){

//  bool match = false;
//  int nL3 = L3cands.size();

//  float minDR = 0.1;
//  float theDR = 100;
//  HLTMuonCand theL3;
//  theL3.pt        = -1000;
//  theL3.eta       = -1000;
//  theL3.phi       = -1000;
//  theL3.trkpt     = -1000;
//  theL3.ecalDep   = -1000;
//  theL3.hcalDep   = -1000;
//  theL3.trkDep    = -1000;
//  theL3.ecalDep05 = -1000;
//  theL3.hcalDep05 = -1000;
//  theL3.ecalDep1  = -1000;
//  theL3.hcalDep1  = -1000;
  
//  for ( std::vector<HLTMuonCand>::const_iterator it = L3cands.begin(); it != L3cands.end(); ++it ) {
//    theDR = deltaR(it -> eta, it -> phi, mu.eta, mu.phi);
//    if (theDR < minDR){
//      minDR = theDR;
//      match = true;
//      theL3 = *it;
//    }
//  }
//  return theL3;
//}


bool matchMuonWithL3(MuonCand mu, std::vector<HLTMuonCand> L3cands){

  bool match = false;
  float minDR = 0.1;
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


//L1MuonCand matchL1(MuonCand mu, std::vector<L1MuonCand> L1cands){

//  bool match = false;
//  int nL1 = L1cands.size();

//  float minDR = 0.3;
//  float theDR = 100;
//  L1MuonCand theL1;
//  theL1.pt        = -1000;
//  theL1.eta       = -1000;
//  theL1.phi       = -1000;
  
//  for ( std::vector<L1MuonCand>::const_iterator it = L1cands.begin(); it != L1cands.end(); ++it ) {
//    theDR = deltaR(it -> eta, it -> phi, mu.eta, mu.phi);
//    if (theDR < minDR){
//      minDR = theDR;
//      match = true;
//      theL1 = *it;
//    }
//  }
//  return theL1;
//}

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
