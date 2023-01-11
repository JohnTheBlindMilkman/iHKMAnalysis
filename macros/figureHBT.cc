// ==========================================================================
// TherminaCorr - an extension of THERMINATOR 2 therm2_hbtfit functionality
// Author: Jedrzej Kolas
// Code at: https://github.com/JohnTheBlindMilkman/TherminaCorr
// ==========================================================================

#include <iostream>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TString.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TAxis.h>
#include <TPaveText.h>
#include <TText.h>
#include <TLegend.h>
#include <TColor.h>
#include <TVirtualPad.h>

using namespace std;

const int nPads = 6;
const EColor color[] = {kGreen,kRed,kOrange};

void preparegraph(TGraphErrors *aGE)
{
    aGE->SetTitle("");
    aGE->SetLineWidth(2);
    aGE->SetMarkerStyle(20);
    //aGE->SetMarkerSize(2);
    aGE->GetXaxis()->SetLimits(101.,599.999);
}

void preparegraph(TGraphErrors *aGE, int col, bool isLambda = false)
{
    double min,max;

    if(isLambda)
    {
        min = -0.2;
        max = 1.7;
    }
    else
    {
        min = 1.7;
        max = 8.3;
    }
    aGE->SetTitle("");
    aGE->SetMarkerStyle(20);
    if(col == 0)
    {
        aGE->SetMarkerColor(color[col]+2);
        aGE->SetLineColor(color[col]+2);
    }
    else
    {
        aGE->SetMarkerColor(color[col]);
        aGE->SetLineColor(color[col]);
    }
    aGE->SetMinimum(min);
    aGE->SetMaximum(max);
    aGE->GetXaxis()->SetLimits(101.,599.999);
}

void figureHBT()
{  
    gStyle->SetPalette(kLake);

    const TString folderPath = "/home/jedkol/lustre/hades/user/kjedrzej/iHKM/ic/output/urqmd";
    const TString tName[] = {"#pi^{+} - #pi^{+}","#pi^{-} - #pi^{-}","#pi^{0} - #pi^{0}"};
    const TString pName[] = {"pipi","piMpiM","pi0pi0"};
    const TString fName[] = {"pipi","pimpim","pi0pi0"};
    const TString epsType[] = {"E0","E1","E2","E3","E4","E5","E6"};
    const TString modelName[] = {"gRinv_0_3","gRout_0_3","gRside_0_3","gRlong_0_3","gLambdaInv_0_3","gLambdaOSL_0_3"};
    const TString expName[] = {"geExpInv","geExpO","geExpS","geExpL","geExpLamInv","geLamosl"};
    const TString rName[] = {"R_{inv}","R_{out}","R_{side}","R_{long}","#lambda_{inv}","#lambda_{osl}"};
    const double epsVal[] = {0.0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6};
    const int NoExp = 0, NoMod = 1;

    TGraphErrors *gModel[NoMod][nPads];
    TGraphAsymmErrors *gaeExp[NoExp][nPads];
    TFile *fExpData,*fModelData;
    TLegend *legExp,*legMod;
    TCanvas *c1;
    TVirtualPad *tvp;
    int maxExpGraph;
    TString drawCom;

    for(int i = 0; i < NoExp; i++)  //input experimental data
    {
        fExpData = TFile::Open(Form("/home/jedkol/Downloads/modelSR/macrosForHBT/%sR3DExp.root",fName[i].Data()),"read");
        if(fExpData->IsOpen())
        {
            cout << "File " << Form("%sR3DExp.root",fName[i].Data()) << " is open" << endl;

            for(int j = 0; j < nPads; j++)
                gaeExp[i][j] = (TGraphAsymmErrors*) fExpData->Get(expName[j].Data());
        }

        fExpData->Close();
    }

    for(int i = 0; i < NoMod; i++)  //input model data
    {
        fModelData = TFile::Open(Form("%s/parameterFit%s.root",folderPath.Data(),pName[0].Data()),"read");
        if(fModelData->IsOpen())
        {
            cout << "File " << Form("parameterFit%s.root",pName[0].Data()) << " is open" << endl;

            for(int j = 0; j < nPads; j++)
                gModel[i][j] = (TGraphErrors*) fModelData->Get(modelName[j].Data());
        }

        fModelData->Close();
    }

    c1 = new TCanvas("c1","",900,900);
    c1->SetMargin(0.,0.,0.2,0.2);
    c1->Divide(2,3,0,0);

    legMod = new TLegend(0.6,0.5,0.8,1.);
    legMod->SetBorderSize(0);
    legMod->SetFillStyle(0);

    for(int i = 0; i < nPads; i++)
    {
        tvp = c1->cd(i+1);
        tvp->SetGrid(1,1);
        gModel[0][i]->GetXaxis()->SetNdivisions(508);
        gModel[0][i]->GetYaxis()->SetNdivisions(505);

        gModel[0][i]->GetYaxis()->SetLabelFont(42);
        gModel[0][i]->GetYaxis()->SetLabelOffset();

        gModel[0][i]->GetYaxis()->SetTitle(Form("%s %s",rName[i].Data(),"[fm]"));
        gModel[0][i]->GetYaxis()->SetTitleOffset(0.8);
        gModel[0][i]->GetYaxis()->CenterTitle(true);
            
        if(i > 3)
        {
            gModel[0][i]->GetXaxis()->SetTitle("k_{T} [MeV]");
            gModel[0][i]->GetXaxis()->SetTitleSize(0.08);
            gModel[0][i]->GetXaxis()->SetTitleOffset(1.1);
            gModel[0][i]->GetXaxis()->SetLabelSize(0.08);
            gModel[0][i]->GetXaxis()->SetLabelFont(42);
            gModel[0][i]->GetXaxis()->SetLabelOffset();

            gModel[0][i]->GetYaxis()->SetTitleSize(0.09);
            gModel[0][i]->GetYaxis()->SetTitleOffset(0.85);
            gModel[0][i]->GetYaxis()->SetLabelSize(0.08);
            gModel[0][i]->GetYaxis()->SetLabelOffset(0.01);
        }
        else
        {
            gModel[0][i]->GetYaxis()->SetTitleSize(0.1);
            gModel[0][i]->GetYaxis()->SetTitleOffset(0.85);
            gModel[0][i]->GetYaxis()->SetLabelSize(0.1);
            gModel[0][i]->GetYaxis()->SetLabelOffset(0.01);
        }

        if(i == 0 || i == 2 || i == 4)
            tvp->SetLeftMargin(0.2);
        else
            tvp->SetRightMargin(0.2);

        if(i < 4)
            preparegraph(gModel[0][i],0);
        else 
            preparegraph(gModel[0][i],0,true);

        drawCom = "ap";
        if(i == 1 || i == 3 || i == 5)
            drawCom += " y+";

        gModel[0][i]->Draw(drawCom.Data());

        if(i == 0)
            legMod->AddEntry(gModel[0][i],tName[0].Data(),"p");
    }

    c1->Draw();
    c1->SaveAs(folderPath + "/outGraphs/HBT.png");
    c1->SaveAs(folderPath + "/outGraphs/HBT.eps");
}