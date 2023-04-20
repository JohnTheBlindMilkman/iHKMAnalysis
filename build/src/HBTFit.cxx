#include <TPad.h>

#include "HBTFit.h"
#include "THGlobal.h"

HBTFit::HBTFit()
{
}

HBTFit::~HBTFit()
{
}

TH3D* HBTFit::getfitprojc(TH3D *expden, TF3 *funqk)
{
    TH3D *projhist = new TH3D(*expden);

    for (int q1int = 1; q1int <= expden->GetNbinsX(); q1int++)
        for (int q2int = 1; q2int <= expden->GetNbinsY(); q2int++) 
            for (int q3int = 1; q3int <= expden->GetNbinsZ(); q3int++) 
                projhist->SetBinContent(q1int, q2int, q3int,expden->GetBinContent(q1int, q2int, q3int) * funqk->Eval(expden->GetXaxis()->GetBinCenter(q1int),expden->GetYaxis()->GetBinCenter(q2int),expden->GetZaxis()->GetBinCenter(q3int)));

    return projhist;
}

TH1D* HBTFit::getfitprojc(TH1D *expden, TF1 *funqk)
{
    TH1D *projhist = new TH1D(*expden);

    for(int iter = 1; iter <= expden->GetNbinsX(); iter++)
        projhist->SetBinContent(iter,expden->GetBinContent(iter) * funqk->Eval(expden->GetBinCenter(iter)));

    return projhist;
}

TH1D* HBTFit::getFitHisto(TH1D *expden, TF1 *funqk)
{
    TH1D *projhist = new TH1D(*expden);

    for(int iter = 1; iter <= expden->GetNbinsX(); iter++)
        projhist->SetBinContent(iter,funqk->Eval(expden->GetBinCenter(iter)));

    return projhist;
}

void HBTFit::preparepad(TH1D *hExp, TH1D* hFit)
{
    gPad->SetFillStyle(4000);
    gPad->SetFillColor(0);
    gPad->SetRightMargin(0.02);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.12);
    gPad->SetGridy();

    hExp->Draw("HISTPE1");
    hFit->Draw("SAMEHISTL");
}

void HBTFit::preparehist(TH1D *hist, int projType, int wType, TString type)
{
    if(projType < 0 || projType > 3 || wType < 1 || wType > 3)
    {
        PRINT_MESSAGE("<HBTFit::preparehist>: Unknown projection type or width");
        exit(THGlobal::Error::generalunsuportedValue);
    }

    if (!type.CompareTo("FIT")) 
    {
        hist->SetLineColor(8);
        hist->SetLineWidth(2);
        hist->SetLineStyle(2);
    }
    else if(!type.CompareTo("CF"))
    {
        hist->SetMarkerSize(1.0);
        hist->SetMarkerColor(2);
        hist->SetMarkerStyle(kOpenCircle);
        hist->SetTitle(Form(";q_{%s} [GeV/c];C(q_{%s})",THGlobal::sProjNames[projType].Data(),THGlobal::sProjNames[projType].Data()));
        hist->SetMinimum(hist->GetMinimum()*0.9);
        hist->SetMaximum(hist->GetMaximum()*1.1);
        hist->GetXaxis()->SetLabelSize(0.055);
        hist->GetYaxis()->SetLabelSize(0.055);
        hist->GetXaxis()->SetTitleSize(0.055);
        hist->GetYaxis()->SetTitleSize(0.055);
        hist->GetYaxis()->SetTitleOffset(0.8);
    }
    else
    {
        PRINT_MESSAGE("<HBTFit::preparehist>: Unknown histogram type");
        exit(THGlobal::Error::generalunsuportedValue);
    }

    hist->SetName(Form("%sproj_%s_%d",type.Data(),THGlobal::sProjNames[projType].Data(),wType));
}

void HBTFit::setErrors(TH1D *hout, TH1D *hNum, TH1D *hDen)
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

        // if F = x/y dF = sqrt(dx^2/y^2 + x^2dy^2/y^4 - 2xdxdy/y^3)
        if (fabs(vDen) > std::numeric_limits<double>::epsilon()) // make sure we don't divide by 0
            vErr = TMath::Sqrt((eNum*eNum)/(vDen*vDen) + ((vNum*vNum)*(eDen*eDen))/(vDen*vDen*vDen*vDen) - (2*vNum*eNum*eDen)/(vDen*vDen*vDen));
        hout->SetBinError(i,vErr);
    }
}

void HBTFit::setErrors(TH3D *hout, TH3D *hNum, TH3D *hDen)
{
    const int iterX = hout->GetNbinsX(), iterY = hout->GetNbinsY(), iterZ = hout->GetNbinsZ();
    double vErr = 0, vNum = 0, vDen = 0, eNum = 0, eDen = 0;
    for (int xx = 1; xx <= iterX; xx++)
        for (int yy = 1; yy <= iterY; yy++)
            for (int zz = 1; zz <= iterZ; zz++)
            {
                vErr = 0;
                vNum = hNum->GetBinContent(xx,yy,zz);
                eNum = hNum->GetBinError(xx,yy,zz);
                vDen = hDen->GetBinContent(xx,yy,zz);
                eDen = hDen->GetBinError(xx,yy,zz);

                // if F = x/y dF = sqrt(dx^2/y^2 + x^2dy^2/y^4 - 2xdxdy/y^3)
                if (fabs(vDen) > std::numeric_limits<double>::epsilon()) // make sure we don't divide by 0
                    vErr = TMath::Sqrt((eNum*eNum)/(vDen*vDen) + ((vNum*vNum)*(eDen*eDen))/(vDen*vDen*vDen*vDen) - (2*vNum*eNum*eDen)/(vDen*vDen*vDen));
                hout->SetBinError(xx,yy,zz,vErr);
            }
}

TH1D* HBTFit::getproj(TH3D *numq, TH3D *denq, int nproj, int wbin, double norm)
{
    Double_t intnexp = norm;
    Double_t intdexp = 1.0;
    TString sProj = nullptr;
    TH1D *hproj,*denbuf,*numbuf;

    Int_t binc = numq->GetXaxis()->FindFixBin(0.0);
    Int_t binm = binc + wbin;

    numq->GetXaxis()->SetRange(binc, binm);
    numq->GetYaxis()->SetRange(binc, binm);
    numq->GetZaxis()->SetRange(binc, binm);
    denq->GetXaxis()->SetRange(binc, binm);
    denq->GetYaxis()->SetRange(binc, binm);
    denq->GetZaxis()->SetRange(binc, binm);

    switch (nproj) 
    {
        case 0:
            sProj = "x";
            denq->GetXaxis()->SetRange(1,numq->GetNbinsX());
            numq->GetXaxis()->SetRange(1,numq->GetNbinsX());
            break;
        case 1:
            sProj = "y";
            denq->GetYaxis()->SetRange(1,numq->GetNbinsY());
            numq->GetYaxis()->SetRange(1,numq->GetNbinsY());
            break;
        case 2:
            sProj = "z";
            denq->GetZaxis()->SetRange(1,numq->GetNbinsZ());
            numq->GetZaxis()->SetRange(1,numq->GetNbinsZ());
            break;
        default:
            PRINT_MESSAGE("<HBTFit::getproj>: Unknown projection type");
            exit(THGlobal::Error::generalunsuportedValue);
    }
    
    denbuf = new TH1D(*((TH1D *) denq->Project3D(sProj)));
    numbuf = new TH1D(*((TH1D *) numq->Project3D(sProj)));

    hproj = new TH1D(*numbuf);
    hproj->Sumw2();
    hproj->Reset("ICE");
    hproj->Divide(numbuf, denbuf, 1.0, 1.0, "");
    //this->setErrors(hproj,numbuf,denbuf);

    /*for (int iter=1; iter<hproj->GetNbinsX(); iter++)
        if (numbuf->GetBinContent(iter)) 
        {
            Double_t dn = numbuf->GetBinError(iter);
            Double_t an = numbuf->GetBinContent(iter);
            Double_t dd = denbuf->GetBinError(iter);
            Double_t ad = denbuf->GetBinContent(iter);
            hproj->SetBinError(iter, TMath::Sqrt((dn*dn*ad*ad + dd*dd*an*an + dd*dd*dn*dn)/(ad*ad*ad*ad))); // change this to match testHBT (set global erros and don't do them on projection) - JJ
        }*/

    hproj->Scale(intdexp/intnexp);

    return hproj;
}

Double_t HBTFit::fungek(Double_t *x, Double_t *par)
{
    Double_t qosq = x[0]*x[0];
    Double_t qssq = x[1]*x[1];
    Double_t qlsq = x[2]*x[2];
    Double_t lam =  TMath::Abs(par[1]);

    Double_t gpart = exp((-par[2]*par[2]*qosq-par[3]*par[3]*qssq-par[4]*par[4]*qlsq)/0.038938); // make it a constant? - JJ

    return (par[0] * (1 + lam*gpart));
}

Double_t HBTFit::fungek1D(Double_t *x, Double_t *par)
{
    Double_t qinv = x[0]*x[0];
    Double_t lam =  TMath::Abs(par[0]);
    Double_t gpart = exp((-par[1]*par[1]*qinv)/0.038938);

    return (par[2] * (1 + lam*gpart));
}