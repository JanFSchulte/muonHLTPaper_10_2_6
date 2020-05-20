

#include "SetStyle.C"

void Plot_L1EffwrtOffline(){

TStyle *gStyle = new TStyle("gStyle","Style for P-TDR");
SetStyle st;
st.SetPars(gStyle);

TFile *f1=new TFile("../L1EffFullStats_PostFilter_L1EffwrtOffline.root","READ");
TFile *f2=new TFile("../L1EffFullStatsQ8_PostFilter_L1EffwrtOffline.root","READ");


//TFile *f1=new TFile("../PromptRecoA1_PostFilter_L1EffwrtOffline_L1Qua12_pt22.root","READ");
//TFile *f2=new TFile("../PromptRecoA1_PostFilter_L1EffwrtOffline_L1Qua8_pt15.root","READ");

TEfficiency *Pt_v1=(TEfficiency*)f1->Get("muonPt");
TEfficiency *Pt_v2=(TEfficiency*)f2->Get("muonPt");

TEfficiency *Eta_v1=(TEfficiency*)f1->Get("muonEta");
TEfficiency *Eta_v2=(TEfficiency*)f2->Get("muonEta");

TEfficiency *Phi_v1=(TEfficiency*)f1->Get("muonPhi");
TEfficiency *Phi_v2=(TEfficiency*)f2->Get("muonPhi");

TEfficiency *nVtx_v1=(TEfficiency*)f1->Get("muonEffnVtx");//muonEffnVtx
TEfficiency *nVtx_v2=(TEfficiency*)f2->Get("muonEffnVtx");

gStyle->SetTitleXSize(0.05); 
gStyle->SetTitleYOffset(1.25); 
gStyle->SetTitleXOffset(1.0); 
gStyle->SetTitleYSize(0.05);
gStyle->SetLabelColor(1, "XYZ");
gStyle->SetLabelFont(42, "XYZ");
gStyle->SetLabelOffset(0.005, "XYZ");
gStyle->SetLabelSize(0.035, "XYZ");



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

leg_1D[k0]=new TLegend(0.4,0.7,0.60,0.9);
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




Pt_v1->SetTitle(" ; p_{T} (#mu) [GeV]  ; L1 Efficiency");
Eta_v1->SetTitle(" ; #eta (#mu)   ; L1 Efficiency");
Phi_v1->SetTitle(" ; #phi (#mu) ; L1 Efficiency");
nVtx_v1->SetTitle(" ; Number of reconstructed primary vertices; L1 Efficiency");

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
graph1A->SetMinimum(0.8);
graph1A->SetMaximum(1.1);
graph1A->GetXaxis()->SetRangeUser(25,150.);
graph1A->GetYaxis()->SetTitleSize(0.05); 
graph1A->GetYaxis()->SetTitleOffset(1.25); 
//graph1A->GetXaxis()->SetTitleOffset(1.35); 
graph1A->GetXaxis()->SetTitleSize(0.05); 

graph2A->SetMinimum(0.8);
graph2A->SetMaximum(1.1);
graph2A->GetXaxis()->SetRangeUser(25,150.);
gPad->Update();
leg_1D[0]->Draw();

 gPad->RedrawAxis();
latex->DrawLatex(0.95, 0.96, "5 fb^{-1} (2018) (13 TeV)");
latexCMS->DrawLatex(0.15,0.955,"CMS");
latexCMSExtra->DrawLatex(0.275,0.955,"Preliminary");


c1->SaveAs("Pt_L1EffwrtOffline.pdf");
c1->SaveAs("Pt_L1EffwrtOffline.png");





TCanvas *c2=new TCanvas("c2","c2");
c2->cd();

Eta_v1->Draw();
Eta_v2->Draw("same");
gPad->Update();
auto graph1B = Eta_v1->TEfficiency::GetPaintedGraph();
auto graph2B = Eta_v2->TEfficiency::GetPaintedGraph();
graph1B->SetMinimum(0.8);
graph1B->SetMaximum(1.1);
graph1B->GetXaxis()->SetRangeUser(-2.4,2.4);
graph1B->GetYaxis()->SetTitleSize(0.05); 
graph1B->GetYaxis()->SetTitleOffset(1.25); 
//graph1B->GetXaxis()->SetTitleOffset(1.35); 
graph1B->GetXaxis()->SetTitleSize(0.05); 

graph2B->SetMinimum(0.8);
graph2B->SetMaximum(1.1);
graph2B->GetXaxis()->SetRangeUser(-2.4,2.4);
gPad->Update();
leg_1D[1]->Draw();
 gPad->RedrawAxis();
latex->DrawLatex(0.95, 0.96, "5 fb^{-1} (2018) (13 TeV)");
latexCMS->DrawLatex(0.15,0.955,"CMS");
latexCMSExtra->DrawLatex(0.275,0.955,"Preliminary");
c2->SaveAs("Eta_L1EffwrtOffline.pdf");
c2->SaveAs("Eta_L1EffwrtOffline.png");

TCanvas *c3=new TCanvas("c3","c3");
c3->cd();

Phi_v1->Draw();
Phi_v2->Draw("same");
gPad->Update();
auto graph1C = Phi_v1->TEfficiency::GetPaintedGraph();
auto graph2C = Phi_v2->TEfficiency::GetPaintedGraph();
graph1C->SetMinimum(0.8);
graph1C->SetMaximum(1.1);
graph1C->GetXaxis()->SetRangeUser(-3.14,3.14);
graph1C->GetYaxis()->SetTitleSize(0.05); 
graph1C->GetYaxis()->SetTitleOffset(1.25); 
//graph1C->GetXaxis()->SetTitleOffset(1.35); 
graph1C->GetXaxis()->SetTitleSize(0.05); 

graph2C->SetMinimum(0.8);
graph2C->SetMaximum(1.1);
graph2C->GetXaxis()->SetRangeUser(-3.14,3.14);
gPad->Update();
leg_1D[2]->Draw();
 gPad->RedrawAxis();
latex->DrawLatex(0.95, 0.96, "5 fb^{-1} (2018) (13 TeV)");
latexCMS->DrawLatex(0.15,0.955,"CMS");
latexCMSExtra->DrawLatex(0.275,0.955,"Preliminary");
c3->SaveAs("Phi_L1EffwrtOffline.pdf");
c3->SaveAs("Phi_L1EffwrtOffline.png");






TCanvas *c4=new TCanvas("c4","c4");
c4->cd();

nVtx_v1->Draw();
nVtx_v2->Draw("same");
gPad->Update();
auto graph1D = nVtx_v1->TEfficiency::GetPaintedGraph();
auto graph2D = nVtx_v2->TEfficiency::GetPaintedGraph();
graph1D->SetMinimum(0.8);
graph1D->SetMaximum(1.1);
graph1D->GetXaxis()->SetRangeUser(0,60);
graph1D->GetYaxis()->SetTitleSize(0.05); 
graph1D->GetYaxis()->SetTitleOffset(1.25); 
//graph1D->GetXaxis()->SetTitleOffset(1.35); 
graph1D->GetXaxis()->SetTitleSize(0.05); 

graph2D->SetMinimum(0.8);
graph2D->SetMaximum(1.1);
graph2D->GetXaxis()->SetRangeUser(0,60);
gPad->Update();
leg_1D[3]->Draw();
 gPad->RedrawAxis();
latex->DrawLatex(0.95, 0.96, "5 fb^{-1} (2018) (13 TeV)");
latexCMS->DrawLatex(0.15,0.955,"CMS");
latexCMSExtra->DrawLatex(0.275,0.955,"Preliminary");
c4->SaveAs("nVtx_L1EffwrtOffline.pdf");
c4->SaveAs("nVtx_L1EffwrtOffline.png");



























}
