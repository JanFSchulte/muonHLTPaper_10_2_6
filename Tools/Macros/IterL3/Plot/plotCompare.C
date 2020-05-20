#include "SetStyle.C"
void plotCompare(){

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

 
	TFile *f=new TFile("../../CascadeTkMu/CTk2018v7_efficiency_prefilter.root","READ");
  TEfficiency *effPt2016 = (TEfficiency*)f->Get("muonPt");
  TEfficiency *effEta2016 = (TEfficiency*)f->Get("muonEta");
  TEfficiency *effPhi2016 = (TEfficiency*)f->Get("muonPhi");
  TEfficiency *effnVtx2016 = (TEfficiency*)f->Get("nvtx");

   TFile *f2=new TFile("2018v2_IterL3preFilter.root","READ");


  TEfficiency *effPt2018 = (TEfficiency*)f2->Get("muonPt");
  TEfficiency *effEta2018 = (TEfficiency*)f2->Get("muonEta");
  TEfficiency *effPhi2018 = (TEfficiency*)f2->Get("muonPhi");
  TEfficiency *effnVtx2018 = (TEfficiency*)f2->Get("nvtx");

   TFile *f3=new TFile("../PromptRecoDv2_IterL3preFilter.root","READ");

//  TH1F *effPt2018NoID = (TH1F*)f3->Get("OIplusIOL2plusIOL1effPt");
//  TH1F *effEta2018NoID = (TH1F*)f3->Get("OIplusIOL2plusIOL1effEta");
//  TH1F *effPhi2018NoID = (TH1F*)f3->Get("OIplusIOL2plusIOL1effPhi");
//  TH1F *effnVtx2018NoID = (TH1F*)f3->Get("OIplusIOL2plusIOL1effnVtx");

  TEfficiency *effPt2018NoID = (TEfficiency*)f3->Get("OIIOIOPt");
  TEfficiency *effEta2018NoID = (TEfficiency*)f3->Get("OIIOIOEta");
  TEfficiency *effPhi2018NoID = (TEfficiency*)f3->Get("OIIOIOPhi");
  TEfficiency *effnVtx2018NoID = (TEfficiency*)f3->Get("OIIOIOnVtx");





TLatex * latex = new TLatex();
latex->SetTextFont(42);
latex->SetTextAlign(31);
latex->SetTextSize(0.03);
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

  leg_1D[k0]=new TLegend(0.2586957,0.804878,0.7341137,0.9320557,NULL,"brNDC");
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
 plotPad->DrawFrame(30,0.9,500,1.05,";p_{T} (#mu) [GeV] ; HLT reconstruction efficiency");
 leg_1D[0]->AddEntry(effPt2016,"cascade || tracker muon","lp");
 leg_1D[0]->AddEntry(effPt2018NoID,"iterative","lp");
 leg_1D[0]->AddEntry(effPt2018,"iterative + muon ID","lp");
 effPt2016->SetLineColor(kBlue);
 effPt2016->SetMarkerColor(kBlue);
 effPt2018->SetMarkerColor(kRed);
 effPt2018->SetLineColor(kRed);
 effPt2016->SetMarkerStyle(21);
 effPt2018->SetMarkerStyle(22);
 effPt2018NoID->SetMarkerStyle(23);
 effPt2016->Draw("pe same ");
 effPt2018->Draw("pe same ");
 effPt2018NoID->Draw("pe same ");
 leg_1D[0]->Draw();


latex->DrawLatex(0.95, 0.96, "X.X fb^{-1} (2018), 2.5 fb^{-1} (2018) (13 TeV)");
latexCMS->DrawLatex(0.15,0.955,"CMS");
latexCMSExtra->DrawLatex(0.275,0.955,"Preliminary");



 gPad->RedrawAxis();
 c1->SaveAs("CompareEffvsPt_2018.pdf");
 c1->SaveAs("CompareEffvsPt_2018.png");






 TCanvas *c2=new TCanvas("c2","c2"); 
 c2->cd();
 plotPad->DrawFrame(-2.4,0.9,2.4,1.05,";#eta (#mu) ; HLT reconstruction efficiency");
 leg_1D[1]->AddEntry(effEta2016,"cascade || tracker muon","l");
 leg_1D[1]->AddEntry(effEta2018NoID,"iterative","l");
 leg_1D[1]->AddEntry(effEta2018,"iterative + muon ID","l");
 effEta2016->SetLineColor(kBlue);
 effEta2018->SetLineColor(kRed);
 effEta2016->SetMarkerColor(kBlue);
 effEta2018->SetMarkerColor(kRed);
 effEta2016->SetMarkerStyle(21);
 effEta2018->SetMarkerStyle(22);
 effEta2018NoID->SetMarkerStyle(22);

 effEta2016->Draw("pe same");
 effEta2018->Draw("pe same");
 effEta2018NoID->Draw("pe same");
 leg_1D[1]->Draw();
 gPad->RedrawAxis();
latex->DrawLatex(0.95, 0.96, "X.X fb^{-1} (2018), 2.5 fb^{-1} (2018) (13 TeV)");
latexCMS->DrawLatex(0.15,0.955,"CMS");
latexCMSExtra->DrawLatex(0.275,0.955,"Preliminary");

 c2->SaveAs("CompareEffvsEta_2018.pdf");
 c2->SaveAs("CompareEffvsEta_2018.png");



 TCanvas *c3=new TCanvas("c3","c3");
 c3->cd();
 plotPad->DrawFrame(-3.14,0.9,3.14,1.05,";#phi (#mu) ; HLT reconstruction effciency");
 leg_1D[2]->AddEntry(effPhi2016,"cascade || tracker muon","l");
 leg_1D[2]->AddEntry(effPhi2018NoID,"iterative","l");
 leg_1D[2]->AddEntry(effPhi2018,"iterative + muon ID","l");
 effPhi2016->SetLineColor(kBlue);
 effPhi2018->SetLineColor(kRed);
 effPhi2016->SetMarkerColor(kBlue);
 effPhi2018->SetMarkerColor(kRed);
 effPhi2016->SetMarkerStyle(21);
 effPhi2018->SetMarkerStyle(22);
 effPhi2018NoID->SetMarkerStyle(23);


 effPhi2016->Draw("pe same");
 effPhi2018->Draw("pe same");
 effPhi2018NoID->Draw("pe same");
 leg_1D[2]->Draw();
 gPad->RedrawAxis();
latex->DrawLatex(0.95, 0.96, "X.X fb^{-1} (2018), 2.5 fb^{-1} (2018) (13 TeV)");
latexCMS->DrawLatex(0.15,0.955,"CMS");
latexCMSExtra->DrawLatex(0.275,0.955,"Preliminary");

 c3->SaveAs("CompareEffvsPhi_2018.pdf");
 c3->SaveAs("CompareEffvsPhi_2018.png");



 TCanvas *c4=new TCanvas("c4","c4");
 c4->cd();
 plotPad->DrawFrame(0,0.9,60,1.05,"; Number of reconstructed primary vertices ; HLT reconstruction efficiency");
 leg_1D[3]->AddEntry(effnVtx2016,"cascade || tracker muon","l");
 leg_1D[3]->AddEntry(effnVtx2018NoID,"iterative","l");
 leg_1D[3]->AddEntry(effnVtx2018,"iterative + muon ID","l");
 effnVtx2016->SetLineColor(kBlue);
 effnVtx2018->SetLineColor(kRed);
 effnVtx2016->SetMarkerColor(kBlue);
 effnVtx2018->SetMarkerColor(kRed);
 effnVtx2016->SetMarkerStyle(21);
 effnVtx2018->SetMarkerStyle(22);
 effnVtx2018NoID->SetMarkerStyle(23);


 effnVtx2016->Draw("pe same");
 effnVtx2018->Draw("pe same");
 effnVtx2018NoID->Draw("pe same");
 leg_1D[3]->Draw();
latex->DrawLatex(0.95, 0.96, "X.X fb^{-1} (2018), 2.5 fb^{-1} (2018) (13 TeV)");
latexCMS->DrawLatex(0.15,0.955,"CMS");
latexCMSExtra->DrawLatex(0.275,0.955,"Preliminary");

gPad->RedrawAxis();
 c4->SaveAs("CompareEffvsnVtx_2018.pdf");
 c4->SaveAs("CompareEffvsnVtx_2018.png");




















}
