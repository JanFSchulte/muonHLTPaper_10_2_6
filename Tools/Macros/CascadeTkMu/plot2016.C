#include "SetStyle.C"
void plot2016(){

TStyle *gStyle = new TStyle("gStyle","Style for P-TDR");
SetStyle st;
st.SetPars(gStyle);



gStyle->SetTitleXSize(0.05); 
gStyle->SetTitleYOffset(1.25); 
gStyle->SetTitleXOffset(1.0); 
gStyle->SetTitleYSize(0.05);
gStyle->SetLabelColor(1, "XYZ");
gStyle->SetLabelFont(42, "XYZ");
gStyle->SetLabelOffset(0.005, "XYZ");
gStyle->SetLabelSize(0.035, "XYZ");

 
	TFile *f=new TFile("CTk2018v5_efficiency_prefilter.root","READ");
/*
  TH1F *OIeffPt = (TH1F*)f->Get("muonPt");
  TH1F *OIeffEta = (TH1F*)f->Get("muonEta");
  TH1F *OIeffPhi = (TH1F*)f->Get("muonPhi");
  TH1F *OIeffnVtx = (TH1F*)f->Get("nVtx");
  TH1F *OIplusIOL2effPt = (TH1F*)f->Get("muonPt_cascade");
  TH1F *OIplusIOL2effEta = (TH1F*)f->Get("muonEta_cascade");
  TH1F *OIplusIOL2effPhi = (TH1F*)f->Get("muonPhi_cascade");
  TH1F *OIplusIOL2effnVtx = (TH1F*)f->Get("nVtx_cascade");
  TH1F *OIplusIOL2plusIOL1effPt = (TH1F*)f->Get("muonPt_tkmu");
  TH1F *OIplusIOL2plusIOL1effEta = (TH1F*)f->Get("muonEta_tkmu");
  TH1F *OIplusIOL2plusIOL1effPhi = (TH1F*)f->Get("muonPhi_tkmu");
  TH1F *OIplusIOL2plusIOL1effnVtx = (TH1F*)f->Get("nVtx_tkmu");

*/
  TEfficiency *OIeffPt = (TEfficiency*)f->Get("muonPt");
  TEfficiency *OIeffEta = (TEfficiency*)f->Get("muonEta");
  TEfficiency *OIeffPhi = (TEfficiency*)f->Get("muonPhi");
  TEfficiency *OIeffnVtx = (TEfficiency*)f->Get("nvtx");
  TEfficiency *OIplusIOL2effPt = (TEfficiency*)f->Get("muonPt_cascade");
  TEfficiency *OIplusIOL2effEta = (TEfficiency*)f->Get("muonEta_cascade");
  TEfficiency *OIplusIOL2effPhi = (TEfficiency*)f->Get("muonPhi_cascade");
  TEfficiency *OIplusIOL2effnVtx = (TEfficiency*)f->Get("nvtx_cascade");
  TEfficiency *OIplusIOL2plusIOL1effPt = (TEfficiency*)f->Get("muonPt_tkmu");
  TEfficiency *OIplusIOL2plusIOL1effEta = (TEfficiency*)f->Get("muonEta_tkmu");
  TEfficiency *OIplusIOL2plusIOL1effPhi = (TEfficiency*)f->Get("muonPhi_tkmu");
  TEfficiency *OIplusIOL2plusIOL1effnVtx = (TEfficiency*)f->Get("nvtx_tkmu");
std::cout << "banana" << std::endl;
/*  OIplusIOL2plusIOL1effEta->GetYaxis()->SetRangeUser(0.80,1.05);
  OIplusIOL2plusIOL1effEta->GetYaxis()->SetTitle(" HLT reconstruction efficiency"); 
  OIplusIOL2plusIOL1effEta->GetYaxis()->SetTitleSize(0.05); 
  OIplusIOL2plusIOL1effEta->GetYaxis()->SetTitleOffset(1.25); 
  OIplusIOL2plusIOL1effEta->GetXaxis()->SetTitleSize(0.05); 
  OIplusIOL2plusIOL1effEta->GetXaxis()->SetTitle("#eta (#mu) ");

  OIplusIOL2plusIOL1effPt->GetYaxis()->SetRangeUser(0.80,1.05);
  OIplusIOL2plusIOL1effPt->GetXaxis()->SetRangeUser(30,500);
  OIplusIOL2plusIOL1effPt->GetYaxis()->SetTitle(" HLT reconstruction efficiency");
  OIplusIOL2plusIOL1effPt->GetYaxis()->SetTitleSize(0.05); 
  OIplusIOL2plusIOL1effPt->GetYaxis()->SetTitleOffset(1.25); 
  OIplusIOL2plusIOL1effPt->GetXaxis()->SetTitleSize(0.05); 
  OIplusIOL2plusIOL1effPt->GetXaxis()->SetTitle("p_{T} (#mu) [GeV] ");

  OIplusIOL2plusIOL1effPhi->GetYaxis()->SetRangeUser(0.80,1.05);
  OIplusIOL2plusIOL1effPhi->GetYaxis()->SetTitle(" HLT reconstruction efficiency");
  OIplusIOL2plusIOL1effPhi->GetYaxis()->SetTitleSize(0.05); 
  OIplusIOL2plusIOL1effPhi->GetYaxis()->SetTitleOffset(1.25); 
  OIplusIOL2plusIOL1effPhi->GetXaxis()->SetTitleSize(0.05); 
  OIplusIOL2plusIOL1effPhi->GetXaxis()->SetTitle("#phi (#mu) ");

  OIplusIOL2plusIOL1effnVtx->GetYaxis()->SetRangeUser(0.80,1.05);
  OIplusIOL2plusIOL1effnVtx->GetXaxis()->SetRangeUser(0,60);
  OIplusIOL2plusIOL1effnVtx->GetYaxis()->SetTitle(" HLT reconstruction efficiency");
  OIplusIOL2plusIOL1effnVtx->GetYaxis()->SetTitleSize(0.05); 
  OIplusIOL2plusIOL1effnVtx->GetYaxis()->SetTitleOffset(1.25); 
  OIplusIOL2plusIOL1effnVtx->GetXaxis()->SetTitleSize(0.05); 
  OIplusIOL2plusIOL1effnVtx->GetXaxis()->SetTitle("Number of reconstructed primary vertices ");

*/
std::cout << "kiwi" << std::endl; 

TLatex * latex = new TLatex();
latex->SetTextFont(42);
latex->SetTextAlign(31);
latex->SetTextSize(0.04);
latex->SetNDC(true);
TLatex * latexCMS = new TLatex();
latexCMS->SetTextFont(61);
latexCMS->SetTextSize(0.055);
latexCMS->SetNDC(true);
TLatex *latexCMSExtra = new TLatex();
latexCMSExtra->SetTextFont(52);
latexCMSExtra->SetTextSize(0.03);
latexCMSExtra->SetNDC(true) ;



  char Legname1[100];


  TLegend *leg_1D[24];

  for(int k0=0;k0<24;k0++){
  sprintf(Legname1,"leg_1D%i",k0);

  leg_1D[k0]=new TLegend(0.2086957,0.204878,0.7341137,0.3520557,NULL,"brNDC");
  leg_1D[k0]->SetTextFont(62);
  leg_1D[k0]->SetLineColor(1);
  leg_1D[k0]->SetLineStyle(1);
  leg_1D[k0]->SetLineWidth(3);
  leg_1D[k0]->SetFillColor(0);
  leg_1D[k0]->SetFillStyle(1001);
  leg_1D[k0]->SetShadowColor(0);
  leg_1D[k0]->SetDrawOption(0);
  leg_1D[k0]->SetBorderSize(0);
  leg_1D[k0]->SetTextSize(0.04);

  }


 TCanvas *c1=new TCanvas("c1","c1");
 c1->cd();
 TPad *plotPad = new TPad("plotPad","plotPad",0,0,1,1);
// plotPad->UseCurrentStyle();
 plotPad->Draw();
 plotPad->cd();
 plotPad->DrawFrame(30,0.8,500,1.05,";p_{T} (#mu) [GeV] ; HLT reconstruction efficiency");
 //plotPad->SetLogx();
 leg_1D[0]->AddEntry(OIeffPt,"cascade || tracker muon","lp");
 leg_1D[0]->AddEntry(OIplusIOL2effPt,"cascade","lp");
 leg_1D[0]->AddEntry(OIplusIOL2plusIOL1effPt,"tracker muon","lp");
 OIplusIOL2plusIOL1effPt->SetLineColor(kBlue);
 OIplusIOL2effPt->SetLineColor(kRed);
 OIeffPt->SetLineColor(kBlack);
 OIplusIOL2plusIOL1effPt->SetFillColor(kRed-10);
 OIplusIOL2effPt->SetFillColor(kRed);
 OIeffPt->SetFillColor(kRed+2);
 OIplusIOL2plusIOL1effPt->SetMarkerColor(kBlue);
 OIplusIOL2effPt->SetMarkerColor(kRed);
 OIeffPt->SetMarkerColor(kBlack);
 OIplusIOL2plusIOL1effPt->SetMarkerStyle(20);
 OIplusIOL2effPt->SetMarkerStyle(21);
 OIeffPt->SetMarkerStyle(22);
 

 OIplusIOL2plusIOL1effPt->Draw("pe same");
 OIplusIOL2effPt->Draw("pe same ");
 OIeffPt->Draw("pe same ");
 leg_1D[0]->Draw();

latex->DrawLatex(0.95, 0.96, "X.X fb^{-1} (2018) (13 TeV)");
latexCMS->DrawLatex(0.15,0.955,"CMS");
latexCMSExtra->DrawLatex(0.275,0.955,"Preliminary");



 gPad->RedrawAxis();
 c1->SaveAs("CascadeTkMuEffvsPt_2018.pdf");
 c1->SaveAs("CascadeTkMuEffvsPt_2018.png");






 TCanvas *c2=new TCanvas("c2","c2"); 
 c2->cd();
 plotPad->DrawFrame(-2.4,0.8,2.4,1.05,";#eta (#mu) ; HLT reconstruction efficiency");
 leg_1D[1]->AddEntry(OIeffEta,"cascade || tracker muon","lp");
 leg_1D[1]->AddEntry(OIplusIOL2effEta,"cascade","lp");
 leg_1D[1]->AddEntry(OIplusIOL2plusIOL1effEta,"tracker muon","lp");
 OIplusIOL2plusIOL1effEta->SetLineColor(kBlue);
 OIplusIOL2effEta->SetLineColor(kRed);
 OIeffEta->SetLineColor(kBlack);
 OIplusIOL2plusIOL1effEta->SetFillColor(kRed-10);
 OIplusIOL2effEta->SetFillColor(kRed);
 OIeffEta->SetFillColor(kRed+2);
 OIplusIOL2plusIOL1effEta->SetMarkerColor(kBlue);
 OIplusIOL2effEta->SetMarkerColor(kRed);
 OIeffEta->SetMarkerColor(kBlack);
 OIplusIOL2plusIOL1effEta->SetMarkerStyle(20);
 OIplusIOL2effEta->SetMarkerStyle(21);
 OIeffEta->SetMarkerStyle(22);
 

 OIplusIOL2plusIOL1effEta->Draw("pe same");
 OIplusIOL2effEta->Draw("pe same");
 OIeffEta->Draw("pe same");
 leg_1D[1]->Draw();
 gPad->RedrawAxis();
latex->DrawLatex(0.95, 0.96, "X.X fb^{-1} (2018) (13 TeV)");
latexCMS->DrawLatex(0.15,0.955,"CMS");
latexCMSExtra->DrawLatex(0.275,0.955,"Preliminary");

c2->SaveAs("CascadeTkMuEffvsEta_2018.pdf");
 c2->SaveAs("CascadeTkMuEffvsEta_2018.png");



 TCanvas *c3=new TCanvas("c3","c3");
 c3->cd();
 plotPad->DrawFrame(-3.14,0.8,3.14,1.05,";#phi (#mu) ; HLT Reconstruction Effciency");
 leg_1D[2]->AddEntry(OIeffPhi,"cascade || tracker muon","lp");
 leg_1D[2]->AddEntry(OIplusIOL2effPhi,"cascade","lp");
 leg_1D[2]->AddEntry(OIplusIOL2plusIOL1effPhi,"tracker muon","lp");
 OIplusIOL2plusIOL1effPhi->SetLineColor(kBlue);
 OIplusIOL2effPhi->SetLineColor(kRed);
 OIeffPhi->SetLineColor(kBlack);
 OIplusIOL2plusIOL1effPhi->SetFillColor(kRed-10);
 OIplusIOL2effPhi->SetFillColor(kRed);
 OIeffPhi->SetFillColor(kRed+2);
 OIplusIOL2plusIOL1effPhi->SetMarkerColor(kBlue);
 OIplusIOL2effPhi->SetMarkerColor(kRed);
 OIeffPhi->SetMarkerColor(kBlack);
 OIplusIOL2plusIOL1effPhi->SetMarkerStyle(20);
 OIplusIOL2effPhi->SetMarkerStyle(21);
 OIeffPhi->SetMarkerStyle(22);
 
 OIplusIOL2plusIOL1effPhi->Draw("pe same");
 OIplusIOL2effPhi->Draw("pe same");
 OIeffPhi->Draw("pe same");
 leg_1D[2]->Draw();
 gPad->RedrawAxis();
latex->DrawLatex(0.95, 0.96, "X.X fb^{-1} (2018) (13 TeV)");
latexCMS->DrawLatex(0.15,0.955,"CMS");
latexCMSExtra->DrawLatex(0.275,0.955,"Preliminary");

c3->SaveAs("CascadeTkMuEffvsPhi_2018.pdf");
 c3->SaveAs("CascadeTkMuEffvsPhi_2018.png");



 TCanvas *c4=new TCanvas("c4","c4");
 c4->cd();
 plotPad->DrawFrame(0,0.8,60,1.05,"; Number of reconstructed primary vertices ; HLT reconstruction efficiency");

 leg_1D[3]->AddEntry(OIeffnVtx,"cascade || tracker muon","lp");
 leg_1D[3]->AddEntry(OIplusIOL2effnVtx,"cascade","lp");
 leg_1D[3]->AddEntry(OIplusIOL2plusIOL1effnVtx,"tracker muon","lp");
 OIplusIOL2plusIOL1effnVtx->SetLineColor(kBlue);
 OIplusIOL2effnVtx->SetLineColor(kRed);
 OIeffnVtx->SetLineColor(kBlack);
 OIplusIOL2plusIOL1effnVtx->SetMarkerColor(kBlue);
 OIplusIOL2effnVtx->SetMarkerColor(kRed);
 OIeffnVtx->SetMarkerColor(kBlack);
 OIplusIOL2plusIOL1effnVtx->SetMarkerStyle(20);
 OIplusIOL2effnVtx->SetMarkerStyle(21);
 OIeffnVtx->SetMarkerStyle(22);
 
 //OIplusIOL2plusIOL1effnVtx->Rebin(2);
 //OIplusIOL2effnVtx->Rebin(2);
 //OIeffnVtx->Rebin(2);
 OIplusIOL2plusIOL1effnVtx->Draw("pe same");
 OIplusIOL2effnVtx->Draw("pe same");
 OIeffnVtx->Draw("pe same");
 leg_1D[3]->Draw();
latex->DrawLatex(0.95, 0.96, "X.X fb^{-1} (2018) (13 TeV)");
latexCMS->DrawLatex(0.15,0.955,"CMS");
latexCMSExtra->DrawLatex(0.275,0.955,"Preliminary");

gPad->RedrawAxis();
 c4->SaveAs("CascadeTkMuEffvsnVtx_2018.pdf");
 c4->SaveAs("CascadeTkMuEffvsnVtx_2018.png");




















}
