#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TColor.h"
#include "TAxis.h"
#include "TLine.h"
#include "TStyle.h"
#include "TMath.h"

void SetErrors(TH1D *hout, TH1D *hNum, TH1D *hDen)
{
    const int iterMax = hout->GetNbinsX();
    double vErr = 0, vNum = 0, vDen = 0, eNum = 0, eDen = 0;
    for (int i = 1; i <= iterMax; i++)
    {
        vErr = 0;
        vNum = hNum->GetBinContent(i);
        eNum = hNum->GetBinError(i);
        vDen = hDen->GetBinContent(i);
        eDen = hDen->GetBinError(i);

        // propagacja błędów dla funcji num/den z uwzględnieniem korelacji pomiędzy składowymi
        if (fabs(vDen) > std::numeric_limits<double>::epsilon())
            vErr = TMath::Sqrt((eNum*eNum)/(vDen*vDen) + ((vNum*vNum)*(eDen*eDen))/(vDen*vDen*vDen*vDen) - (2*vNum*eNum*eDen)/(vDen*vDen*vDen));
        hout->SetBinError(i,vErr);
    }
}

void PrepareHistogram(TH1D *hist,int colour, int marker, int rebin = 0)
{
    hist->SetMarkerColor(colour);
    hist->SetMarkerStyle(marker);
    hist->SetLineColor(colour);
    hist->SetFillColor(colour);

    if(rebin)
    {
        hist->Rebin(rebin);
        hist->Scale(1./rebin);
    }
}

void testHBT()
{
    gStyle->SetOptStat(0);

    const TString filepath = "/home/jedkol/Downloads/HADES/data/iHKM/femtopp19a.root";
    const TString outFilePath = "/home/jedkol/Downloads/HADES/iHKM/macros/output/femtoProton.root";
    const double hXmin = 0.0, hXmax = 0.2;
    const int tcBlue = TColor::GetColor("#2d7f9d"), tcRed = TColor::GetColor("#dc7684");
    
    TH1D *hDen,*hNumQ,*hNumQSC,*hRatQ,*hRatQSC;
    TLine *lLine;
    TFile *inpFile, *outFile;

    TCanvas *canv = new TCanvas("canv","",800,450);

    inpFile = TFile::Open(filepath);
    hDen = (TH1D*) inpFile->Get("den1d");
    hNumQ = (TH1D*) inpFile->Get("num1d");
    hRatQ = new TH1D(*hNumQ);
    hRatQ->SetName("rat1d");
    hNumQSC = (TH1D*) inpFile->Get("num1dqsc");
    hRatQSC = new TH1D(*hNumQSC);
    hRatQSC->SetName("rat1dqsc");

    hRatQ->Divide(hDen);
    SetErrors(hRatQ,hNumQ,hDen);
    PrepareHistogram(hRatQ,tcBlue,20,2);
    hRatQ->GetXaxis()->SetRangeUser(hXmin,hXmax);
    hRatQ->GetYaxis()->SetRangeUser(0.,2.);
    hRatQ->SetTitle("1D p - p correlation funciton"); // 1D p - p correlation funciton, 1D #pi^{+} - #pi^{+} correlation funciton

    hRatQSC->Divide(hDen);
    SetErrors(hRatQSC,hNumQSC,hDen);
    PrepareHistogram(hRatQSC,tcRed,20,2);

    lLine = new TLine(hXmin,1.,hXmax,1);
    lLine->SetLineColor(kBlack);
    lLine->SetLineStyle(kDashed);

    hRatQ->Draw("p");
    hRatQSC->Draw("p same");
    lLine->Draw("same");

    outFile = TFile::Open(outFilePath,"recreate");
    hRatQ->Write();
    hRatQSC->Write();
    canv->SaveAs("./output/canv1Dpp.png");
}