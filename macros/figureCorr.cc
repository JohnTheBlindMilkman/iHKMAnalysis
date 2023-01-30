// ==========================================================================
// TherminaCorr - an extension of THERMINATOR 2 therm2_hbtfit functionality
// Author: Jedrzej Kolas
// Code at: https://github.com/JohnTheBlindMilkman/TherminaCorr
// ==========================================================================

#include <iostream>
#include <TFile.h>
#include <TCanvas.h>
#include <TString.h>
#include <TFile.h>
#include <TH1.h>
#include <TColor.h>
#include <TStyle.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TLegend.h>

using namespace std;

void figureCorr()
{
    gStyle->SetOptStat(0);

    TString filePath = "/home/jedkol/lustre/hades/user/kjedrzej/iHKM/ic/output/urqmd/";
    const int kTset[] = {0};
    const int len = sizeof(kTset)/sizeof(kTset[1]);
    const TString projName[] = {"out", "side", "long"};
    const Color_t cGreen = TColor::GetColor(26,201,33);
    const Color_t cRed = TColor::GetColor(252,90,58);
    TDirectory *dir = new TDirectory();

    TCanvas *c[len];
    TFile *inFile;
    TH1D *hCF[3], *hFIT[3], *hFITcopy[3];
    TLegend *leg = new TLegend(0.75,0.75,0.95,0.95,"","NDC");
    int iter = 0;

    if(! dir->GetDirectory(filePath + "outGraphs/"))
        gSystem->Exec("mkdir " + filePath + "outGraphs/");

    for(int i : kTset)
    {
        c[iter] = new TCanvas(Form("c%d",iter),"",1500,500);
        c[iter]->Divide(3,1);
        inFile = TFile::Open(Form("%shbtfitpipi%da.root",filePath.Data(),i));

        for(auto j : {0,1,2})
        {
            c[iter]->cd(j+1)->SetMargin(0.15,0.05,0.15,0.01);
            hCF[j] = (TH1D*) inFile->Get(Form("CFproj_%s_3",projName[j].Data()));
            hCF[j]->GetYaxis()->SetTitleOffset(1.25);
            hCF[j]->GetYaxis()->SetNdivisions(507);
            hCF[j]->GetXaxis()->SetTitleOffset(1.1);
            hCF[j]->GetXaxis()->SetNdivisions(507);
            hCF[j]->SetLineColor(cGreen);
            hCF[j]->SetLineWidth(5);
            hCF[j]->SetFillColor(cGreen);
            hCF[j]->SetFillStyle(3390);
            hCF[j]->SetMarkerColor(cGreen);
            hCF[j]->SetMarkerStyle(8);                

            hFIT[j] = (TH1D*) inFile->Get(Form("FITproj_%s_3",projName[j].Data()));
            hFIT[j]->SetLineColor(cRed);
            hFIT[j]->SetLineStyle(1);
            hFITcopy[j] = (TH1D*) hFIT[j]->Clone(Form("FITproj_%s_3_clone",projName[j].Data()));

            hFIT[j]->SetFillColor(cRed);
            hFIT[j]->SetFillStyle(3344);

            hCF[j]->Draw("p e4");
            hFIT[j]->Draw("e4 same");
            hFITcopy[j]->Draw("c hist same");
            if(j == 0)
            {
                if(iter == 0)
                {
                    leg->AddEntry(hCF[j],"Data","pf");
                    leg->AddEntry(hFIT[j],"Fit","lf");
                }
                leg->Draw("same");
            }
        }
        c[iter]->SaveAs(Form("%soutGraphs/canvaspipi%da.png",filePath.Data(),i));

        iter++;
    }
}