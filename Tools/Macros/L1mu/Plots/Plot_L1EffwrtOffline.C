

#include "SetStyle.C"

void Plot_L1EffwrtOffline(){

TStyle *gStyle = new TStyle("gStyle","Style for P-TDR");
SetStyle st;
st.SetPars(gStyle);


TFile *f1=new TFile("../PromptRecoA1_PostFilter_L1EffwrtOffline_L1Qua12_pt22.root","READ");
TFile *f2=new TFile("../PromptRecoA1_PostFilter_L1EffwrtOffline_L1Qua8_pt15.root","READ");

TEfficiency *Pt_v1=(TEfficiency*)f1->Get("muonPt");
TEfficiency *Pt_v2=(TEfficiency*)f2->Get("muonPt");

TEfficiency *Eta_v1=(TEfficiency*)f1->Get("muonEta");
TEfficiency *Eta_v2=(TEfficiency*)f2->Get("muonEta");

TEfficiency *Phi_v1=(TEfficiency*)f1->Get("muonPhi");
TEfficiency *Phi_v2=(TEfficiency*)f2->Get("muonPhi");

TEfficiency *nVtx_v1=(TEfficiency*)f1->Get("muonEffnVtx");//muonEffnVtx
TEfficiency *nVtx_v2=(TEfficiency*)f2->Get("muonEffnVtx");



Pt_v1->SetLineColor(4);
Pt_v1->SetLineWidth(1);
Pt_v1->SetMarkerStyle(20);
Pt_v1->SetMarkerColor(4);
Pt_v1->SetMarkerSize(1.5);
Pt_v2->SetLineColor(2);
Pt_v2->SetLineWidth(1);
Pt_v2->SetMarkerStyle(22);
Pt_v2->SetMarkerColor(2);
Pt_v2->SetMarkerSize(1.5);





Eta_v1->SetLineColor(4);
Eta_v1->SetLineWidth(1);
Eta_v1->SetMarkerStyle(20);
Eta_v1->SetMarkerColor(4);
Eta_v1->SetMarkerSize(1.5);
Eta_v2->SetLineColor(2);
Eta_v2->SetLineWidth(1);
Eta_v2->SetMarkerStyle(22);
Eta_v2->SetMarkerColor(2);
Eta_v2->SetMarkerSize(1.5);


Phi_v1->SetLineColor(4);
Phi_v1->SetLineWidth(1);
Phi_v1->SetMarkerStyle(20);
Phi_v1->SetMarkerColor(4);
Phi_v1->SetMarkerSize(1.5);
Phi_v2->SetLineColor(2);
Phi_v2->SetLineWidth(1);
Phi_v2->SetMarkerStyle(22);
Phi_v2->SetMarkerColor(2);
Phi_v2->SetMarkerSize(1.5);


nVtx_v1->SetLineColor(4);
nVtx_v1->SetLineWidth(1);
nVtx_v1->SetMarkerStyle(20);
nVtx_v1->SetMarkerColor(4);
nVtx_v1->SetMarkerSize(1.5);
nVtx_v2->SetLineColor(2);
nVtx_v2->SetLineWidth(1);
nVtx_v2->SetMarkerStyle(22);
nVtx_v2->SetMarkerColor(2);
nVtx_v2->SetMarkerSize(1.5);









char Legname1[100];


TLegend *leg_1D[24]; 

for(int k0=0;k0<24;k0++){
sprintf(Legname1,"leg_1D%i",k0);

leg_1D[k0]=new TLegend(0.4,0.4,0.60,0.6);
leg_1D[k0]->SetTextFont(62);
leg_1D[k0]->SetLineColor(1);
leg_1D[k0]->SetLineStyle(1);
leg_1D[k0]->SetLineWidth(3);
leg_1D[k0]->SetFillColor(0);
leg_1D[k0]->SetFillStyle(1001);
leg_1D[k0]->SetShadowColor(0);
leg_1D[k0]->SetDrawOption(0);
leg_1D[k0]->SetBorderSize(0);
leg_1D[k0]->SetTextSize(0.05);

}


TPaveText *pCMS = new TPaveText(0.1638796,0.8101045,0.8645485,0.9721254,"brNDC");
  pCMS->SetBorderSize(0);
  pCMS->SetFillStyle(0);
  pCMS->SetTextAlign(11);
  pCMS->SetTextFont(42);
  pCMS->SetTextSize(0.04);
  pCMS->AddText("                                        #sqrt{s}= 13 TeV (2018)");
  pCMS->AddText("               ");
  pCMS->AddText("               ");
  pCMS->AddText("               ");
  pCMS->AddText("               ");
  pCMS->AddText("#bf{CMS} ");
  pCMS->AddText("               ");
  pCMS->AddText("               ");
  pCMS->AddText("#it{Preliminary}");

Pt_v1->SetTitle(" ; p_{T} (#mu) [GeV]  ; L1 Efficiency");
Eta_v1->SetTitle(" ; #eta (#mu)   ; L1 Efficiency");
Phi_v1->SetTitle(" ; #phi (#mu) ; L1 Efficiency");
nVtx_v1->SetTitle(" ; n_{reco. vertices}  ; L1 Efficiency");


leg_1D[0]->AddEntry(Pt_v1,"Q >=12, p_{T}^{L1} >22 GeV","l");
leg_1D[0]->AddEntry(Pt_v2,"Q >=8, p_{T}^{L1} >15 GeV","l");

leg_1D[1]->AddEntry(Eta_v1,"Q >=12, p_{T}^{L1} >22 GeV","l");
leg_1D[1]->AddEntry(Eta_v2,"Q >=8, p_{T}^{L1} >15 GeV","l");

leg_1D[2]->AddEntry(Phi_v1,"Q >=12, p_{T}^{L1} >22 GeV","l");
leg_1D[2]->AddEntry(Phi_v2,"Q >=8, p_{T}^{L1} >15 GeV","l");

leg_1D[3]->AddEntry(nVtx_v1,"Q >=12, p_{T}^{L1} >22 GeV","l");
leg_1D[3]->AddEntry(nVtx_v2,"Q >=8, p_{T}^{L1} >15 GeV","l");




TCanvas *c1=new TCanvas("c1","c1");
c1->cd();

Pt_v1->Draw();
Pt_v2->Draw("same");
gPad->Update();
auto graph1A = Pt_v1->TEfficiency::GetPaintedGraph();
auto graph2A = Pt_v2->TEfficiency::GetPaintedGraph();
graph1A->SetMinimum(0.5);
graph1A->SetMaximum(1.08);
graph1A->GetXaxis()->SetRangeUser(25,150.);
graph2A->SetMinimum(0.5);
graph2A->SetMaximum(1.08);
graph2A->GetXaxis()->SetRangeUser(25,150.);
gPad->Update();
leg_1D[0]->Draw();
pCMS->Draw();
c1->SaveAs("Pt_L1EffwrtOffline.pdf");
c1->SaveAs("Pt_L1EffwrtOffline.png");





TCanvas *c2=new TCanvas("c2","c2");
c2->cd();

Eta_v1->Draw();
Eta_v2->Draw("same");
gPad->Update();
auto graph1B = Eta_v1->TEfficiency::GetPaintedGraph();
auto graph2B = Eta_v2->TEfficiency::GetPaintedGraph();
graph1B->SetMinimum(0.5);
graph1B->SetMaximum(1.08);
graph1B->GetXaxis()->SetRangeUser(-2.4,2.4);
graph2B->SetMinimum(0.5);
graph2B->SetMaximum(1.08);
graph2B->GetXaxis()->SetRangeUser(-2.4,2.4);
gPad->Update();
leg_1D[1]->Draw();
pCMS->Draw();
c2->SaveAs("Eta_L1EffwrtOffline.pdf");
c2->SaveAs("Eta_L1EffwrtOffline.png");

TCanvas *c3=new TCanvas("c3","c3");
c3->cd();

Phi_v1->Draw();
Phi_v2->Draw("same");
gPad->Update();
auto graph1C = Phi_v1->TEfficiency::GetPaintedGraph();
auto graph2C = Phi_v2->TEfficiency::GetPaintedGraph();
graph1C->SetMinimum(0.5);
graph1C->SetMaximum(1.08);
graph1C->GetXaxis()->SetRangeUser(-3.14,3.14);
graph2C->SetMinimum(0.5);
graph2C->SetMaximum(1.08);
graph2C->GetXaxis()->SetRangeUser(-3.14,3.14);
gPad->Update();
leg_1D[2]->Draw();
pCMS->Draw();
c3->SaveAs("Phi_L1EffwrtOffline.pdf");
c3->SaveAs("Phi_L1EffwrtOffline.png");






TCanvas *c4=new TCanvas("c4","c4");
c4->cd();

nVtx_v1->Draw();
nVtx_v2->Draw("same");
gPad->Update();
auto graph1D = Phi_v1->TEfficiency::GetPaintedGraph();
auto graph2D = Phi_v2->TEfficiency::GetPaintedGraph();
graph1D->SetMinimum(0.5);
graph1D->SetMaximum(1.08);
graph1D->GetXaxis()->SetRangeUser(10,40);
graph2D->SetMinimum(0.5);
graph2D->SetMaximum(1.08);
graph2D->GetXaxis()->SetRangeUser(10,40);
gPad->Update();
leg_1D[3]->Draw();
pCMS->Draw();
c4->SaveAs("nVtx_L1EffwrtOffline.pdf");
c4->SaveAs("nVtx_L1EffwrtOffline.png");



























}
