#include <iostream>
#include <fstream>
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
    TString fitFileTxt = "output.txt";
    ifstream fitFile;
    double val;
    std::vector<double> fitYval;

    fitFile.open(fitFileTxt);
    while (fitFile >> val)
        fitYval.push_back(val);
    fitFile.close();
    
    TH1D* hOut = new TH1D("hout","",fitYval.size(),0.005,0.405);
    for (int i = 1; i <= fitYval.size(); i++)
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