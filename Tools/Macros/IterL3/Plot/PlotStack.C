#include "SetStyle.C"
void PlotStack(){

TStyle *gStyle = new TStyle("gStyle","Style for P-TDR");
SetStyle st;
st.SetPars(gStyle);



	TFile *f=new TFile("../PromptRecoD_IterL3preFilter_save.root","READ");


  TH1F *OIeffPt = (TH1F*)f->Get("OIeffPt");
  TH1F *OIeffEta = (TH1F*)f->Get("OIeffEta");
  TH1F *OIeffPhi = (TH1F*)f->Get("OIeffPhi");
  TH1F *OIeffnVtx = (TH1F*)f->Get("OIeffnVtx");
  TH1F *OIplusIOL2effPt = (TH1F*)f->Get("OIplusIOL2effPt");
  TH1F *OIplusIOL2effEta = (TH1F*)f->Get("OIplusIOL2effEta");
  TH1F *OIplusIOL2effPhi = (TH1F*)f->Get("OIplusIOL2effPhi");
  TH1F *OIplusIOL2effnVtx = (TH1F*)f->Get("OIplusIOL2effnVtx");
  TH1F *OIplusIOL2plusIOL1effPt = (TH1F*)f->Get("OIplusIOL2plusIOL1effPt");
  TH1F *OIplusIOL2plusIOL1effEta = (TH1F*)f->Get("OIplusIOL2plusIOL1effEta");
  TH1F *OIplusIOL2plusIOL1effPhi = (TH1F*)f->Get("OIplusIOL2plusIOL1effPhi");
  TH1F *OIplusIOL2plusIOL1effnVtx = (TH1F*)f->Get("OIplusIOL2plusIOL1effnVtx");

  OIplusIOL2plusIOL1effEta->GetYaxis()->SetRangeUser(0.80,1.05);
  OIplusIOL2plusIOL1effEta->GetYaxis()->SetTitle(" HLT Reconstruction Efficiency"); 
  OIplusIOL2plusIOL1effEta->GetXaxis()->SetTitle("#eta (#mu) ");

  OIplusIOL2plusIOL1effPt->GetYaxis()->SetRangeUser(0.80,1.05);
  OIplusIOL2plusIOL1effPt->GetXaxis()->SetRangeUser(30,500);
  OIplusIOL2plusIOL1effPt->GetYaxis()->SetTitle(" HLT Reconstruction Efficiency");
  OIplusIOL2plusIOL1effPt->GetXaxis()->SetTitle("p_{T} (#mu) [GeV] ");

  OIplusIOL2plusIOL1effPhi->GetYaxis()->SetRangeUser(0.80,1.05);
  OIplusIOL2plusIOL1effPhi->GetYaxis()->SetTitle(" HLT Reconstruction Efficiency");
  OIplusIOL2plusIOL1effPhi->GetXaxis()->SetTitle("#phi (#mu) ");

  OIplusIOL2plusIOL1effnVtx->GetYaxis()->SetRangeUser(0.80,1.05);
  OIplusIOL2plusIOL1effnVtx->GetXaxis()->SetRangeUser(0,60);
  OIplusIOL2plusIOL1effnVtx->GetYaxis()->SetTitle(" HLT Reconstruction Efficiency");
  OIplusIOL2plusIOL1effnVtx->GetXaxis()->SetTitle("Number of reconstructed primary vertices ");







  TPaveText *pCMS = new TPaveText(0.1238796,0.95,0.8645485,0.99,"brNDC");
  pCMS->SetBorderSize(0);
  pCMS->SetFillStyle(0);
  pCMS->SetTextAlign(11);
  pCMS->SetTextFont(42);
  pCMS->SetTextSize(0.04);
  pCMS->AddText("#bf{CMS} #it{Preliminary}                           #sqrt{s}= 13 TeV (2018)");

  char Legname1[100];


  TLegend *leg_1D[24];

  for(int k0=0;k0<24;k0++){
  sprintf(Legname1,"leg_1D%i",k0);

  leg_1D[k0]=new TLegend(0.6086957,0.804878,0.7341137,0.9320557,NULL,"brNDC");
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
 leg_1D[0]->AddEntry(OIeffPt,"OI","f");
 leg_1D[0]->AddEntry(OIplusIOL2effPt,"OI+IOL2","f");
 leg_1D[0]->AddEntry(OIplusIOL2plusIOL1effPt,"OI+IOL2+IOL1","f");
 OIplusIOL2plusIOL1effPt->SetFillColor(kRed-10);
 OIplusIOL2effPt->SetFillColor(kRed);
 OIeffPt->SetFillColor(kRed+2);
 OIplusIOL2plusIOL1effPt->Draw("hist ");
 OIplusIOL2effPt->Draw("hist same ");
 OIeffPt->Draw("hist same ");
 pCMS->Draw();
 leg_1D[0]->Draw();
 gPad->RedrawAxis();
 c1->SaveAs("StackEffvsPt.pdf");
 c1->SaveAs("StackEffvsPt.png");






 TCanvas *c2=new TCanvas("c2","c2"); 
 c2->cd();
 leg_1D[1]->AddEntry(OIeffEta,"OI","f");
 leg_1D[1]->AddEntry(OIplusIOL2effEta,"OI+IOL2","f");
 leg_1D[1]->AddEntry(OIplusIOL2plusIOL1effEta,"OI+IOL2+IOL1","f");
 OIplusIOL2plusIOL1effEta->SetFillColor(kRed-10);
 OIplusIOL2effEta->SetFillColor(kRed);
 OIeffEta->SetFillColor(kRed+2);
 OIplusIOL2plusIOL1effEta->Draw("hist");
 OIplusIOL2effEta->Draw("hist same");
 OIeffEta->Draw("hist same");
 pCMS->Draw();
 leg_1D[1]->Draw();
 gPad->RedrawAxis();
 c2->SaveAs("StackEffvsEta.pdf");
 c2->SaveAs("StackEffvsEta.png");



 TCanvas *c3=new TCanvas("c3","c3");
 c3->cd();
 leg_1D[2]->AddEntry(OIeffPhi,"OI","f");
 leg_1D[2]->AddEntry(OIplusIOL2effPhi,"OI+IOL2","f");
 leg_1D[2]->AddEntry(OIplusIOL2plusIOL1effPhi,"OI+IOL2+IOL1","f");
 OIplusIOL2plusIOL1effPhi->SetFillColor(kRed-10);
 OIplusIOL2effPhi->SetFillColor(kRed);
 OIeffPhi->SetFillColor(kRed+2);
 OIplusIOL2plusIOL1effPhi->Draw("hist");
 OIplusIOL2effPhi->Draw("hist same");
 OIeffPhi->Draw("hist same");
 pCMS->Draw();
 leg_1D[2]->Draw();
 gPad->RedrawAxis();
 c3->SaveAs("StackEffvsPhi.pdf");
 c3->SaveAs("StackEffvsPhi.png");



 TCanvas *c4=new TCanvas("c4","c4");
 c4->cd();
 leg_1D[3]->AddEntry(OIeffnVtx,"OI","f");
 leg_1D[3]->AddEntry(OIplusIOL2effnVtx,"OI+IOL2","f");
 leg_1D[3]->AddEntry(OIplusIOL2plusIOL1effnVtx,"OI+IOL2+IOL1","f");
 OIplusIOL2plusIOL1effnVtx->SetFillColor(kRed-10);
 OIplusIOL2effnVtx->SetFillColor(kRed);
 OIeffnVtx->SetFillColor(kRed+2);
 //OIplusIOL2plusIOL1effnVtx->Rebin(2);
 //OIplusIOL2effnVtx->Rebin(2);
 //OIeffnVtx->Rebin(2);
 OIplusIOL2plusIOL1effnVtx->Draw("hist");
 OIplusIOL2effnVtx->Draw("hist same");
 OIeffnVtx->Draw("hist same");
 pCMS->Draw();
 leg_1D[3]->Draw();
 gPad->RedrawAxis();
 c4->SaveAs("StackEffvsnVtx.pdf");
 c4->SaveAs("StackEffvsnVtx.png");




















}
