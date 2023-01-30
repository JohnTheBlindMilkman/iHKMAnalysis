#ifndef _TH2_HBTFIT_H_
    #define _TH2_HBTFIT_H_

    #include <TH3D.h>
    #include <TF3.h>

    class HBTFit 
    {
        public:
            HBTFit();
            ~HBTFit();

            TH3D *getfitprojc(TH3D *expden, TF3 *funqk);
            TH1D *getfitprojc(TH1D *expden, TF1 *funqk);
            TH1D *getFitHisto(TH1D *expden, TF1 *funqk);
            void preparepad(TH1D *hExp, TH1D* hFit);
            void preparehist(TH1D *hist, int projType, int wType, TString type);
            void setErrors(TH1D *hout, TH1D *hNum, TH1D *hDen);
            void setErrors(TH3D *hout, TH3D *hNum, TH3D *hDen);
            TH1D *getproj(TH3D *numq, TH3D *denq, int nproj, int wbin, double norm);
            Double_t fungek(Double_t *x, Double_t *par);
            Double_t fungek1D(Double_t *x, Double_t *par);
    };

#endif