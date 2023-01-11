#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TColor.h"
#include "TAxis.h"
#include "TLine.h"
#include "TStyle.h"

void testHBT()
{
    gStyle->SetOptStat(0);

    const TString filepath = "/home/jedkol/lustre/hades/user/kjedrzej/iHKM/ic/output/urqmd/femtopp19a.root";
    const double hXmin = 0.0, hXmax = 0.2;
    
    TH1D *hDen,*hNumQ,*hNumQSC;
    TLine *lLine;
    TFile *inpFile;

    TCanvas *canv = new TCanvas("canv","",800,450);

    inpFile = TFile::Open(filepath);
    hDen = (TH1D*) inpFile->Get("den1d");
    hNumQ = (TH1D*) inpFile->Get("num1d");
    hNumQSC = (TH1D*) inpFile->Get("num1dqsc");

    hNumQ->Divide(hDen);
    hNumQ->SetMarkerColor(TColor::GetColor("#2d7f9d"));
    hNumQ->SetLineColor(TColor::GetColor("#2d7f9d"));
    hNumQ->GetXaxis()->SetRangeUser(hXmin,hXmax);
    hNumQ->SetTitle("1D p - p correlation funciton");

    hNumQSC->Divide(hDen);
    hNumQSC->SetMarkerColor(TColor::GetColor("#dc7684"));
    hNumQSC->SetLineColor(TColor::GetColor("#dc7684"));

    lLine = new TLine(hXmin,1.,hXmax,1);
    lLine->SetLineColor(kBlack);
    lLine->SetLineStyle(kDashed);

    hNumQ->Draw();
    hNumQSC->Draw("same");
    lLine->Draw("same");

    canv->SaveAs("./output/canv1Dpp.png");
}