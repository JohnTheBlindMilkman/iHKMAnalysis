#include <iostream>

#include "TBenchmark.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TH1.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TVirtualPad.h"

void prepareGraph(TH1F *h, TString title, int color)
{
    h->SetTitle(title.Data());
    h->SetMarkerStyle(20);
    h->SetMarkerColor(color);
    h->SetLineColor(color);
}

void testiHKM()
{
    gStyle->SetOptStat(0);

    const TString filePath = "/home/jedkol/lustre/hades/user/kjedrzej/iHKM/ic/output/urqmd/";
    const TString treeName = "treeini"; //treeini & treefin -> difference?
    const int noFiles = 100;

    TCanvas *canv;
    TVirtualPad *vp;
    float pt,y,widthPt,widthY;

    TH1F *hMomentum = new TH1F("hMomentum","",50,0,5);
    TH1F *hRapidity = new TH1F("hRapidity","",100,-5,5);    // w STARze nie robią rozkładu rapidity, więc chyba nie ma sensu tego robić
    TChain *chMain = new TChain(treeName.Data());

    for (int i = 0; i < noFiles; i++)
        chMain->Add(Form("%s%d.root",filePath.Data(),i));

    TTreeReader ttReader(chMain);
    TTreeReaderValue<int> npart(ttReader,"npart");
    TTreeReaderArray<float> px(ttReader,"px");
    TTreeReaderArray<float> py(ttReader,"py");
    TTreeReaderArray<float> pz(ttReader,"pz");
    TTreeReaderArray<float> E(ttReader,"E");
    TTreeReaderArray<int> id(ttReader,"id");
    TTreeReaderArray<int> mid(ttReader,"mid");

    gBenchmark->Start("testiHKM");

    widthPt = hMomentum->GetBinWidth(10);
    widthY = hRapidity->GetBinWidth(10);

    while (ttReader.Next())
    {
        for (int i = 0; i < *npart; i++)
        {
            if (id[i] != 2212) // if not proton
                continue;
            if (mid[i] != 0) // if not primordial
                continue;

            y = 0.5 * TMath::Log((E[i] + pz[i])/(E[i] - pz[i]));
            hRapidity->Fill(y);
            if (TMath::Abs(y) < widthY) // midrapidity; has to be equal to rapidity bin width
            {
                pt = sqrtf(px[i]*px[i] + py[i]*py[i]);
                hMomentum->Fill(pt,1./(pt));
            }  
        }
    }

    hMomentum->Scale(1./(hMomentum->GetEntries()*2*TMath::Pi()*widthPt*widthY));
    hRapidity->Scale(1./(hRapidity->GetEntries()*widthY));

    canv = new TCanvas("c","c",1600,800);
    canv->SetMargin(0.3,0.05,0.2,0.05);
    canv->Divide(2,1);

    vp = canv->cd(1);
    vp->SetLogy();
    vp->SetMargin(0.2,0.05,0.1,0.1);
    prepareGraph(hMomentum,"p_{t} distribution for protons;p_{t} [GeV/c]; #frac{1}{2 #pi p_{t}} #frac{d^{2}N}{dp_{t} dy}|_{y=0} [(GeV/c)^{-2}]",2);
    hMomentum->Draw("c hist");

    vp = canv->cd(2);
    vp->SetMargin(0.2,0.05,0.1,0.1);
    prepareGraph(hRapidity,"y distribution for protons;y [GeV/c];#frac{dN}{dy}",4);
    hRapidity->Draw("c hist");

    gBenchmark->Show("testiHKM");

    canv->SaveAs("./output/testiHKM.eps");
    canv->SaveAs("./output/testiHKM.png");
}