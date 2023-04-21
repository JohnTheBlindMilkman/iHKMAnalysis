#include <iostream>
#include "TStyle.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH3.h"
#include "TLine.h"
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

void prepareGraph(TH1D* hist, int col)
{
    hist->SetMarkerColor(col);
    hist->SetMarkerStyle(20);
    hist->SetLineColor(col);

    //hist->GetXaxis()->SetLabelSize();
    //hist->GetXaxis()->SetLabelOffset();
    //hist->GetXaxis()->SetTitleSize();
    //hist->GetXaxis()->SetTitleOffset();

    hist->GetYaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetLabelOffset(0.005);
    hist->GetYaxis()->SetLabelFont(62);
    hist->GetYaxis()->SetTitleSize(0.04);
    hist->GetYaxis()->SetTitleOffset(1.8);
    hist->GetYaxis()->SetTitleFont(62);
    hist->GetYaxis()->SetNdivisions(510);
}

void drawProton3DCF()
{
    gStyle->SetOptStat(0);
    const TString filePath = "/home/jedkol/Downloads/STAR/data/MuDst39_230127.root";
    const TString sProj[] = {"x","y","z"};
    const TString sProjName[] = {"out","side","long"};
    const int cent = 1;
    const double drawMin = 0.0, drawMax = 0.27;
    const int rebin = 4;
    const int wbin = 10;

    double norm = 1.;
    int binc = 0, binm = 0;
    TFile *impFile,*outFile;
    TCanvas *canv;
    TH3D *hNum3D, *hDen3D;
    TH1D *hRatPP[3],*hNumPP[3],*hDenPP[3];
    TLine *lline;

    canv = new TCanvas("canv","",1800,600);
    canv->Divide(3,1);

    lline = new TLine(drawMin,1.,drawMax,1.);
    lline->SetLineStyle(kDashed);

    impFile = TFile::Open(filePath);

    hNum3D = (TH3D*) impFile->Get(Form("NumProtonProton%dBPLCMS3DCF",cent));
    hDen3D = (TH3D*) impFile->Get(Form("DenProtonProton%dBPLCMS3DCF",cent));

    for(int i : {0,1,2}) // loop over projections: out, side, long
    {
        binc = hNum3D->GetXaxis()->FindFixBin(0.0);
        binm = binc + wbin;

        hNum3D->GetXaxis()->SetRange(binc, (i == 0) ? hNum3D->GetNbinsX() : binm); // if this works I will declare myself the king of C++
        hNum3D->GetYaxis()->SetRange(binc, (i == 1) ? hNum3D->GetNbinsY() : binm);
        hNum3D->GetZaxis()->SetRange(binc, (i == 2) ? hNum3D->GetNbinsZ() : binm);
        hDen3D->GetXaxis()->SetRange(binc, (i == 0) ? hNum3D->GetNbinsX() : binm);
        hDen3D->GetYaxis()->SetRange(binc, (i == 1) ? hNum3D->GetNbinsY() : binm);
        hDen3D->GetZaxis()->SetRange(binc, (i == 2) ? hNum3D->GetNbinsZ() : binm);

        hNumPP[i] = (TH1D*) hNum3D->Project3D(sProj[i]);
        hDenPP[i] = (TH1D*) hDen3D->Project3D(sProj[i]);
        hRatPP[i] = new TH1D(*hNumPP[i]);
        hRatPP[i]->SetName(Form("RatProtonProton%dBPLCMS3DCF_%s",cent,sProjName[i].Data()));
        //hRatPP[i]->Sumw2();
        //hRatPP[i]->Reset("ICE");
        hRatPP[i]->Divide(hDenPP[i]);
        hRatPP[i]->Rebin(rebin);
        norm = getNorm(hRatPP[i],0.15,0.3);
        hRatPP[i]->GetXaxis()->SetRangeUser(drawMin,drawMax);
        //hRatPP[i]->GetYaxis()->SetRangeUser(0.9,1.2);
        hRatPP[i]->SetTitle(Form("Projection of 3D p-p correlation for %d-%d%% ceontrality;q_{%s};CF(q_{%s})",cent*10-10,cent*10,sProjName[i].Data(),sProjName[i].Data()));
        prepareGraph(hRatPP[i],TColor::GetColor("#965f77"));
        if (norm > 0)
            hRatPP[i]->Scale(1./norm);
        canv->cd(i+1)->SetMargin(0.15,0.01,0.15,0.1);
        hRatPP[i]->Draw();
        lline->Draw("same");
    }
    
    canv->SaveAs("./output/CF3Dcent0to10.png");
    //outFile = TFile::Open("./output/ProtonProton1BPLCMS3DCF.root","recreate");
    //canv->Write();
    //for( int i : {0,1,2})
        //hRatPP[i]->Write();

    //outFile->Close();
    //impFile->Close();
}