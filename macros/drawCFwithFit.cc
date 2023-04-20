#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TColor.h"
#include "TAxis.h"
#include "TLine.h"
#include "TStyle.h"
#include "TLegend.h"

TH1D* createFitHisto()
{
    const float fitYval[] = {0.186139,1.03439,1.70836,1.89326,1.74052,1.49529,1.28736,1.14161,1.05321,1.00036,0.975295,0.963088,0.960411,0.962584,0.965645,0.96831,0.972307,0.974299,0.978558,0.981151,0.982222,0.984651,0.987045,0.986598,0.987324,0.989202,0.990442,0.990345,0.991446,0.99055,0.992298,0.992094,0.992617,0.993748,0.993622,0.994411,0.993687,0.993381,0.993762,0.994058};
    const int arrSize = 40;
    TH1D* hOut = new TH1D("hout","",arrSize,0.,0.4);
    for (int i = 1; i <= arrSize; i++)
        hOut->SetBinContent(i,fitYval[i-1]);
    
    return hOut;
}

void prepareHistogram(TH1D *hist,int colour, int marker, int rebin = 0)
{
    hist->SetMarkerColor(colour);
    hist->SetMarkerStyle(marker);
    hist->SetMarkerSize(1.2);
    hist->SetLineColor(colour);
    hist->SetLineWidth(3);

    if(rebin)
    {
        hist->Rebin(rebin);
        hist->Scale(1./rebin);
    }
}

void drawCFwithFit()
{
    gStyle->SetOptStat(0);

    const TString filePath = "/home/jedkol/Downloads/HADES/iHKM/macros/output/femtoProton.root";
    const double hXmin = 0.0, hXmax = 0.2;
    const int tcBlue = TColor::GetColor("#2d7f9d"), tcRed = TColor::GetColor("#dc7684");
    
    TH1D *hRatQSC,*hFit;
    TLine *lLine;
    TFile *inpFile;

    TLegend *leg = new TLegend(0.6,0.7,0.9,0.9,"","NDC");

    TCanvas *canv = new TCanvas("canv","",800,800);

    inpFile = TFile::Open(filePath,"read");
    hRatQSC = (TH1D*) inpFile->Get("rat1dqsc");
    hRatQSC->SetTitle("1D p-p Model CF with CorrFit");
    prepareHistogram(hRatQSC,tcBlue,20);
    hRatQSC->Draw("p");

    hFit = createFitHisto();
    prepareHistogram(hFit,tcRed,20);
    hFit->Draw("c same");

    leg->AddEntry(hRatQSC,"iHKM","p");
    leg->AddEntry(hFit,"fit","l");
    leg->Draw("same");

    canv->SaveAs("./output/canv1DppFit.png");
}