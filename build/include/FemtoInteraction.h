#ifndef _FEMTO_INTERACTION_
    #define _FEMTO_INTERACTION_

    #include <iostream>
    #include <complex.h>
    #include <TROOT.h>
    //#include <RDataFrame.hxx>
    #include <TMath.h>
    #include "THGlobal.h"

    class FemtoInteraction
    {
        private:
            std::complex<long double> GetG(long double eta, long double rho, long double hfun);
            std::tuple<std::complex<long double>,std::complex<long double>> GetFC(long double kstar, long double eta, long double hfun);
            std::tuple<long double, long double> BFunPFun(long double eta, long double rho);
            std::complex<long double> GetFFsingle(int sign);
            std::tuple<std::complex<long double>,std::complex<long double>> GetFFdouble();
            void InitGamov();
            double Gamov(double arg);
            long double ChiIm(long double eta);
            long double GetH(long double eta);

            std::complex<long double> d0s,f0s,d0t,f0t; // to-Do: check wheather we need to use long double - JJ
            long double pionac,coulqscpart;
            double mKStarSigned,mKStarOut,mKStarSide,mKStarLong,mROS,mRSS,mRLS,mRSt,mKO,mKS,mKL,mDE,mRO,mRS,mRL,mDT,mKR,mKC,mKP,mRSide,mROut,mRLong,mDTime,mRTrans,mRSidePairCMS,mRLongPairCMS,mDTimePairLCMS,mROutPairCMS,mDTimePairCMS,mRStar;
            int twoSpin,partpid,partpid2,pairtype;

        public:
            FemtoInteraction();
            ~FemtoInteraction();
            double PairKinematics(THGlobal::ParticleCoor part1, THGlobal::ParticleCoor part2, THGlobal::ReferenceFrame refFrame, double ktmin, double ktMax, double btmin, double btMax);
            double GetCoulomb();
            double GetQuantum(float tDiff);
            double GetQuantumCoulomb(float tDiff);
            double GetQuantumCoulombStrong(float tDiff);
            double GetNewQuantumCoulombStrong(float tDiff);

            //template <typename RDF>
            //auto PairKinematics(RDF df, TString sRefFrame);
    };

    inline long double FemtoInteraction::ChiIm(long double eta) { return Gamov(1.0 / (eta * pionac)) / (2.0 * eta);}
    
    /*template <typename RDF>
    auto FemtoInteraction::PairKinematics(RDF df, TString sRefFrame)
    {
        auto sqrtSum = [](float x, float y){return sqrt(x*x + y*y)};
        auto sqrtDiff = [](float t, float z){return sqrt(t*t - z*z)};
        auto dfOut = df.Define("tPt",sqrtSum,{"px","py"}).Define("tMt",sqrtDiff,{"E","pz"}); //this seems too complicated for now, gonna make a simple macro to see how this works - JJ
    }*/

#endif