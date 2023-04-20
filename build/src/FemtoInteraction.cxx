#include "FemtoInteraction.h"

FemtoInteraction::FemtoInteraction() : d0s(0.,0.), f0s(0.,0.), d0t(0.,0.), f0t(0.,0.), twoSpin(0), pairtype(2), pionac(0), coulqscpart(0), mKStarSigned(0), mKStarOut(0), mKStarSide(0), mKStarLong(0), mRO(0), mROS(0), mRS(0), mRSS(0), mRLS(0), mRSt(0), mKO(0), mKS(0), mKL(0), mDE(0), mRL(0), mDT(0), mKR(0), mKC(0), mKP(0), mRSide(0), mROut(0), mRLong(0), mDTime(0), mRTrans(0), mRSidePairCMS(0), mROutPairCMS(0), mRLongPairCMS(0), mDTimePairLCMS(0), mDTimePairCMS(0), mRStar(0)
{
    InitGamov();
}

FemtoInteraction::~FemtoInteraction() 
{
}

double FemtoInteraction::PairKinematics(THGlobal::ParticleCoor part1, THGlobal::ParticleCoor part2, THGlobal::ReferenceFrame refFrame, double ktMin, double ktMax, double btMin, double btMax)
{
    // Calculate pair variables
    double tPx = part1.px+part2.px;
    double tPy = part1.py+part2.py;
    double tPz = part1.pz+part2.pz;
    double tE  = part1.e +part2.e;
    double tPt = tPx*tPx + tPy*tPy;
    double tMt = tE*tE - tPz*tPz;//mCVK;
    double tM  = sqrt(tMt - tPt);
    tMt = sqrt(tMt);
    tPt = sqrt(tPt);
    double mKT = tPt/2.0;
    double pairphi = TMath::ATan2(tPy, tPx);

    if ((mKT < ktMin) || (mKT > ktMax))
        return 0.;

    double mBetat = tPt/tMt;
    if ((mBetat < btMin) || (mBetat > btMax)) 
    {
        mKT = -1.0;
        return 0.;
    }

    THGlobal::ParticleCoor particle1lcms,particle2lcms,particle1prf,particle2prf;

    // Boost to LCMS
    double tBeta = tPz/tE;
    double tGamma = tE/tMt;	    
    mKStarLong = tGamma * (part1.pz - tBeta * part1.e);
    double tE1L = tGamma * (part1.e  - tBeta * part1.pz);
    // Transform to LCMS
    particle1lcms.z = tGamma * (part1.z - tBeta * part1.t);
    particle1lcms.t = tGamma * (part1.t - tBeta * part1.z);
    
    particle2lcms.z = tGamma * (part2.z - tBeta * part2.t);
    particle2lcms.t = tGamma * (part2.t - tBeta * part2.z);
    
    particle1prf.pz = particle1lcms.pz = tGamma * (part1.pz - tBeta * part1.e);
    particle1lcms.e  = tGamma * (part1.e  - tBeta * part1.pz);
    
    particle2prf.pz = particle2lcms.pz = tGamma * (part2.pz - tBeta * part2.e);
    particle2lcms.e  = tGamma * (part2.e  - tBeta * part2.pz);
    
    // Rotate in transverse plane
    mKStarOut  = ( part1.px*tPx + part1.py*tPy)/tPt;
    mKStarSide = (-part1.px*tPy + part1.py*tPx)/tPt;
        
    particle1lcms.px = mKStarOut;
    particle1lcms.py = mKStarSide;
    
    particle2lcms.px = (part2.px*tPx + part2.py*tPy)/tPt;
    particle2lcms.py = (part2.py*tPx - part2.px*tPy)/tPt;;
    
    mKO = particle1lcms.px - particle2lcms.px;
    mKS = particle1lcms.py - particle2lcms.py;
    mKL = particle1lcms.pz - particle2lcms.pz;
    mDE = particle1lcms.e  - particle2lcms.e;

    // save the rotated coordinates in LCMS variables
    particle1lcms.x = ( part1.x*tPx + part1.y*tPy)/tPt;
    particle1lcms.y = (-part1.x*tPy + part1.y*tPx)/tPt;

    particle2lcms.x = ( part2.x*tPx + part2.y*tPy)/tPt;
    particle2lcms.y = (-part2.x*tPy + part2.y*tPx)/tPt;

    // Boost to pair cms
    mKStarOut = tMt/tM * (mKStarOut - tPt/tMt * tE1L);
    
    Double_t tBetat = tPt/tMt;
    Double_t tGammat = 1.0/TMath::Sqrt(1.0-tBetat*tBetat);
    
    particle1prf.x = tGammat*(particle1lcms.x - tBetat*particle1lcms.t);
    particle1prf.t = tGammat*(particle1lcms.t - tBetat*particle1lcms.x);
    
    particle2prf.x = tGammat*(particle2lcms.x - tBetat*particle2lcms.t);
    particle2prf.t = tGammat*(particle2lcms.t - tBetat*particle2lcms.x);
    
    mRO = (particle1lcms.x - particle2lcms.x) / THGlobal::GeVtoFm;
    mRS = (particle1lcms.y - particle2lcms.y) / THGlobal::GeVtoFm;
    mRL = (particle1lcms.z - particle2lcms.z) / THGlobal::GeVtoFm;
    mDT = (particle1lcms.t - particle2lcms.t) / THGlobal::GeVtoFm;

    mKO = particle1lcms.px - particle2lcms.px;

    mKStarSigned = mKStarOut>0.? 1. : -1.;
    mKStarSigned *= sqrt(mKStarSide*mKStarSide + mKStarOut*mKStarOut + mKStarLong*mKStarLong);

    mKR = fabs(mKStarSigned);
    if ( mKR < 1e-10 || fabs(mKR) < 1e-10 ) 
        mKC = 0.0;
    else 
        mKC=mKStarLong/mKR;

    mKP=atan2(mKStarSide,mKStarOut);

    double tDX = part1.x-part2.x;
    double tDY = part1.y-part2.y;
    mRLong = part1.z-part2.z;
    mDTime = part1.t-part2.t;

    mRTrans = tDX>0.? ::sqrt(tDX*tDX+tDY*tDY) : -1.*::sqrt(tDX*tDX+tDY*tDY);
    mROut = (tDX*tPx + tDY*tPy)/tPt;
    mRSide = (-tDX*tPy + tDY*tPx)/tPt;

    mRSidePairCMS = mRSide;
    mRSS = mRSidePairCMS / THGlobal::GeVtoFm;

    mRLongPairCMS = tGamma*(mRLong - tBeta* mDTime);
    mDTimePairLCMS = tGamma*(mDTime - tBeta* mRLong);

    mRLS = mRLongPairCMS / THGlobal::GeVtoFm;
    tBeta = tPt/tMt;
    tGamma = tMt/tM;

    mROutPairCMS = tGamma*(mROut - tBeta* mDTimePairLCMS);
    mROS = mROutPairCMS / THGlobal::GeVtoFm;
    mDTimePairCMS = tGamma*(mDTimePairLCMS - tBeta* mROut);

    mRStar = sqrt(mROutPairCMS*mROutPairCMS + mRSidePairCMS*mRSidePairCMS + mRLongPairCMS*mRLongPairCMS);
    mRSt = mRStar / THGlobal::GeVtoFm;

    return mKStarSigned;

    /*if (writesourcehists) 
    {
        tROutLCMS ->Fill(mROut);
        tRSideLCMS->Fill(mRSide);
        tRLongLCMS->Fill(mRLong);
        tRTimeLCMS->Fill(mDTimePairLCMS);
        tROutSideLongLCMS->Fill(mROut, mRSide, mRLong);
        
        tROutPRF ->Fill(mROutPairCMS);
        tRSidePRF->Fill(mRSidePairCMS);
        tRLongPRF->Fill(mRLongPairCMS);
        tRTimePRF->Fill(mDTimePairCMS);
        tROutSideLongPRF->Fill(mROutPairCMS, mRSidePairCMS, mRLongPairCMS);
        
        tROutLCMSQinv->Fill(mROut, fabs(mKStarSigned)*2.0);
        tRSideLCMSQinv->Fill(mRSide, fabs(mKStarSigned)*2.0);
        tRLongLCMSQinv->Fill(mRLong, fabs(mKStarSigned)*2.0);
        tRTimeLCMSQinv->Fill(mDTimePairLCMS, fabs(mKStarSigned)*2.0);
        
        tROutPRFQinv->Fill(mROutPairCMS, fabs(mKStarSigned)*2.0);
        tRSidePRFQinv->Fill(mRSidePairCMS, fabs(mKStarSigned)*2.0);
        tRLongPRFQinv->Fill(mRLongPairCMS, fabs(mKStarSigned)*2.0);
        tRTimePRFQinv->Fill(mDTimePairCMS, fabs(mKStarSigned)*2.0);
        
        tRInvPRF->Fill(mRStar, 1.0/(mRStar*mRStar));
    }*/
}
// Calculate quantum weight according to B-E or F-D statistics
double FemtoInteraction::GetQuantum(float tDiff)
{
    static int pcount = 0;
    if (tDiff > 500.0)
        return 1.0;

    if (twoSpin == 0) 
        return 1.0 + TMath::Cos(-mKO*mRO - mKS*mRS - mKL*mRL + mDE*mDT);
    else if (twoSpin == 1) 
    {
        if (pcount ==3) 
            return 1.0 + TMath::Cos(-mKO*mRO - mKS*mRS - mKL*mRL + mDE*mDT);
        else 
            return 1.0 - TMath::Cos(-mKO*mRO - mKS*mRS - mKL*mRL + mDE*mDT);
        
        pcount++;
        if (pcount == 4) 
            pcount=0;
    }
}

// Calculate the wave-function modulus sqaured
// for non-identical particles
// that is the Coulomb interaction
double FemtoInteraction::GetCoulomb()
{
    double kstrst = mKStarOut*mROS + mKStarSide*mRSS + mKStarLong*mRLS;

    // Classical limit - if distance is larger than Coulomb radius, 
    // the interaction does not matter
    if (fabs(mRSt) > fabs(pionac)) 
        return (1.0);
    if (fabs(mRSt) == 0.0) 
        return (Gamov(fabs(mKStarSigned)));

    // Classical limit - in the case of large k*r* product we go to 
    // classical coulomb interaction
    if (fabs(mKStarSigned) * mRSt * (1.0 + kstrst/(mRSt*fabs(mKStarSigned)))> 15.0)
        return (1.0 - 1.0/(mRSt*pionac*mKStarSigned*mKStarSigned));
    
    // Calculate the F function
    std::complex<long double> ffplus = GetFFsingle(1);

    if (!finite(ffplus.real())) 
    {
        std::cout << "FFPlus Re not a number !" << " " << std::endl;
        std::cout << mRSt << " " << mKStarSigned << " " << mROS << " " << mRSS << " " << mRLS << std::endl;
    }

    if (!finite(ffplus.imag()))
        std::cout << "FFPlus Im not a number !" << " " << std::endl;

    return (Gamov(fabs(mKStarSigned)) * norm(ffplus) * norm(ffplus));
}

// Calculate the wave-function modulus sqaured
// for identical bosons (symmetrized)
// with Coulomb interaction included
double FemtoInteraction::GetQuantumCoulomb(float tDiff) //this has to be turned into a functor - JJ
{
    if (mRSt < 0.0000000001) // why this arbitrary number?? - JJ
        return 1.0;

    if (tDiff > 500.0)
        return 1.0;

    double kstrst = mKStarOut*mROS + mKStarSide*mRSS + mKStarLong*mRLS;
    int ccase = 0;
    static int pcount = 0;
    int wavesign = 1;

    if (twoSpin == 1) 
    {
        if (pcount == 3)
            wavesign = 1;
        else
            wavesign = -1;

        pcount++;
        if (pcount == 4) 
            pcount = 0;
    }

    // Classical limit - if distance is larger than Coulomb radius, 
    // the interaction does not matter
    if (fabs(mRSt) > fabs(pionac)) 
        return (1.0 + wavesign*cos(2*kstrst));

    // Classical limit - in the case of large k* we go to 
    // classical coulomb interaction
    long double testp = fabs(mKStarSigned) * mRSt * (1.0 + kstrst/(mRSt*fabs(mKStarSigned)));
    long double testm = fabs(mKStarSigned) * mRSt * (1.0 - kstrst/(mRSt*fabs(mKStarSigned)));

    if ((testp> 15.0) && (testm> 15.0))
    {
        double fasymplus  = (1.0 - 1.0/((mRSt+kstrst)*pionac*mKStarSigned*mKStarSigned));
        double fasymminus = (1.0 - 1.0/((mRSt-kstrst)*pionac*mKStarSigned*mKStarSigned));
        return 0.5 * ((fasymplus + fasymminus)*cos(2*kstrst) + (2.0*sqrt(fasymplus*fasymminus)));
    }
    
    std::complex<long double> ffplus, ffminus;
    // Check for the classical limit in both functions separately
    if (((testp< 15.0) && (testm< 15.0))) // ||
    {
        // Calculate the F function
        std::tie(ffplus, ffminus) = GetFFdouble();
        ccase = 1;
    }
    else if (testp < 15.0)
    {
        double asym;
        ffplus = GetFFsingle(1);
        asym = sqrt(1.0 - 1.0/(mRSt*(1.0 - kstrst/(mRSt*fabs(mKStarSigned))*pionac*mKStarSigned*mKStarSigned)))/Gamov(fabs(mKStarSigned));
        
        ffminus = {(asym < 1.0 ? 1.0 + (asym -1.0) *2.0 : 1.0 + (asym -1.0) /2.0),sqrt(asym*asym - ffminus.real()*ffminus.real())};

        ccase = 2;
    }
    else 
    {
        double asym;
        ffminus = GetFFsingle(-1);
        asym = sqrt(1.0 - 1.0/(mRSt*(1.0 + kstrst/(mRSt*fabs(mKStarSigned))*pionac*mKStarSigned*mKStarSigned)))/Gamov(fabs(mKStarSigned));

        ffplus = {(asym < 1.0 ? 1.0 + (asym -1.0) *2.0 : 1.0 + (asym -1.0) /2.0),sqrt(asym*asym - ffplus.real()*ffplus.real())};

        ccase = 3;
    }

    std::complex<long double> expikr = {cos(kstrst),sin(kstrst)};
    std::complex<long double> sterm = expikr * expikr * ffplus * conj(ffminus);
    std::complex<long double> tterm = conj(expikr * expikr) * conj(ffplus) * ffminus;


    if (!finite(ffminus.real())) // for the love of god standarise the error messages - JJ
        std::cout << "FFMinus Re not a number !" << " " << ccase<< std::endl;
    
    if (!finite(ffminus.imag()))
        std::cout << "FFMinus Im not a number !" << " " << ccase<< std::endl;

    if ((ffplus.real() > 2.0) || (ffplus.real() < -2.0))
        std::cout << "FFplus Re wild !" << ffplus.real() << std::endl;

    if ((ffplus.imag() > 2.0) || (ffplus.imag() < -2.0))
        std::cout << "FFplus Im wild !" << ffplus.imag() << std::endl;

    if ((ffminus.real() > 2.0) || (ffminus.real() < -2.0))
        std::cout << "FFminus Re wild !" << ffminus.real() << " " << ccase << std::endl;
    
    if ((ffminus.imag() > 2.0) || (ffminus.imag() < -2.0))
        std::cout << "FFminus Im wild !" << ffminus.imag() << " " << ccase << std::endl;

    coulqscpart = 0.5 * Gamov(fabs(mKStarSigned)) * (norm(ffplus)*norm(ffplus) + norm(ffminus)*norm(ffminus));

    return (0.5 * Gamov(fabs(mKStarSigned)) * (norm(ffplus)*norm(ffplus) + wavesign*sterm.real() + wavesign*tterm.real() + norm(ffminus)*norm(ffminus)));
}

// Calculate the wave-function modulus sqaured
// for identical bosons (symmetrized) or fermions (antisymmetrized)
// with Coulomb interaction and Strong interaction included
double FemtoInteraction::GetQuantumCoulombStrong(float tDiff)
{
    if (mRSt < 0.0000000001)
        return 1.0;

    if (mRSt < 1.0/0.197327) 
        return GetQuantumCoulomb(tDiff);
    
    double tKstRst = mKStarOut*mROS + mKStarSide*mRSS + mKStarLong*mRLS;
    long double kstar = fabs(mKStarSigned);
    long double rho = mRSt * kstar;
    
    int ccase = 0;
    int wavesign = 1;

    // Classical limit - if distance is larger than Coulomb radius, 
    // the interaction does not matter
    //  if (fabs(mRSt) > fabs(pionac)) return (1.0 + wavesign*cos(2*tKstRst));

    // Classical limit - in the case of large k* we go to 
    // classical coulomb interaction
    long double testp = rho * (1.0 + tKstRst/(rho));
    long double testm = rho * (1.0 - tKstRst/(rho));

    std::complex<long double> ffplus, ffminus;
    long double ffmRe, ffpRe;
    if ((testp> 15.0) && (testm> 15.0))
    {
        double asym;
        asym = (1.0 - 1.0/(mRSt*(1.0 - tKstRst/rho)*pionac*kstar*kstar))/Gamov(kstar);
        //      cout << "as1 " << asym << endl;
        asym = sqrt(asym);
        if (asym < 1.0) 
            ffmRe = 1.0 + (asym -1.0) *2.0;
        else
            ffmRe = 1.0 + (asym -1.0) /2.0;

        ffminus = {ffmRe,sqrt(asym*asym - ffminus.real()*ffminus.real())};

        asym = (1.0 - 1.0/(mRSt*(1.0 + tKstRst/rho)*pionac*kstar*kstar))/Gamov(kstar);
        //      cout << "as2 " << asym << endl;
        asym = sqrt(asym);
        if (asym < 1.0) 
            ffpRe = 1.0 + (asym -1.0) *2.0;
        else
            ffpRe = 1.0 + (asym -1.0) /2.0;

        ffplus = {ffpRe,sqrt(asym*asym - ffplus.real()*ffplus.real())};
      
    }
    // Check for the classical limit in both functions separately
    else if (((testp< 15.0) && (testm< 15.0))) // ||
    {
        // Calculate the F function
        std::tie(ffplus,ffminus) = GetFFdouble();
        ccase = 1;
    }
    else if (testp< 15.0)
    {
        double asym;
        ffplus = GetFFsingle(1);
        ffminus = GetFFsingle(-1);
        if ((fabs(ffminus.real()) > 2.0) || fabs(ffminus.imag()) > 2.0) 
        {
            asym = (1.0 - 1.0/(mRSt*(1.0 - tKstRst/(rho)*pionac*kstar*kstar)))/Gamov(kstar);

            asym = sqrt(asym);
            if (asym < 1.0) 
                ffmRe = 1.0 + (asym -1.0) *2.0;
            else
                ffmRe = 1.0 + (asym -1.0) /2.0;

            ffminus = {ffmRe,sqrt(asym*asym - ffminus.real()*ffminus.real())};
        }
        ccase = 2;
    }
    else 
    {
        double asym;
        ffminus = GetFFsingle(-1);
        ffplus = GetFFsingle(1);
        if ((fabs(ffplus.real()) > 2.0) || fabs(ffplus.imag()) > 2.0) 
        {
            asym = (1.0 - 1.0/(mRSt*(1.0 + tKstRst/(rho)*pionac*kstar*kstar)))/Gamov(kstar);
            asym = sqrt(asym);
            if (asym < 1.0) 
                ffpRe = 1.0 + (asym -1.0) *2.0;
            else
                ffpRe = 1.0 + (asym -1.0) /2.0;

            ffplus = {ffpRe,sqrt(asym*asym - ffplus.real()*ffplus.real())};
        }
        ccase = 3;
    }

    long double eta  = 1.0/(kstar*pionac);
    long double hfun = GetH(eta);
    std::complex<long double> gtilde  = GetG(eta, rho, hfun);
    std::complex<long double> gtilor  = gtilde / (long double) mRSt;

    std::complex<long double> fcouls, fcoult;
    std::tie(fcouls,fcoult) = GetFC(kstar, eta, hfun);
    
    std::complex<long double> fgs = gtilor * fcouls;
    long double fgmods = norm(fgs);

    std::complex<long double> fgt = gtilor * fcoult;
    long double fgmodt = norm(fgt);

    std::complex<long double> expikr = {cos(tKstRst),-sin(tKstRst)};
    std::complex<long double> expikrc = conj(expikr);
    std::complex<long double> ffplusc = conj(ffplus);
    std::complex<long double> ffminusc = conj(ffminus);

    std::complex<long double> expikr2 = expikr * expikr;
    std::complex<long double> expikrc2 = conj(expikr2);
    std::complex<long double> sterm = expikr2 * ffplus * ffminusc;
    std::complex<long double> tterm = expikrc2 * ffminus * ffplusc;

    std::complex<long double> epfpc = expikr * ffplus;
    std::complex<long double> emfmc = expikrc * ffminus;

    long double fcgefhs = (fgs.real()*emfmc.real() + fgs.imag()*emfmc.imag());
    long double fcgefgs = (fgs.real()*epfpc.real() + fgs.imag()*epfpc.imag());

    long double fcgefht = (fgt.real()*emfmc.real() + fgt.imag()*emfmc.imag());
    long double fcgefgt = (fgt.real()*epfpc.real() + fgt.imag()*epfpc.imag());

    long double smult = 1+wavesign;

    if (!finite(ffminus.real()))
        std::cout << "FFMinus Re not a number ! " << testp << " " << testm << " " << ccase<< std::endl;
  
    if (!finite(ffminus.imag()))
        std::cout << "FFMinus Im not a number !" << testp << " " << testm << " " << ccase<< std::endl;

    if ((ffplus.real() > 2.0) || (ffplus.real() < -2.0))
        std::cout << "FFplus Re wild !" << ffplus.real() << " case " << ccase << " " << testp << " " << testm << std::endl;

    if ((ffplus.imag() > 2.0) || (ffplus.imag() < -2.0))
        std::cout << "FFplus Im wild !" << ffplus.imag() << " case " << ccase << " " << testp << " " << testm << std::endl;

    if ((ffminus.real() > 2.0) || (ffminus.real() < -2.0))
        std::cout << "FFminus Re wild !" << ffminus.real() << " case " << ccase << std::endl;
    
    if ((ffminus.imag() > 2.0) || (ffminus.imag() < -2.0))
        std::cout << "FFminus Im wild !" << ffminus.imag() << " case " << ccase << std::endl;

    if (twoSpin == 1) 
    {
        wavesign = 1;
        smult = 2;
        long double singlet = (0.5 * Gamov(kstar) * (2.0 * fgmods * smult + norm(ffplus) + norm(ffminus) + wavesign*sterm.real() + wavesign*tterm.real() + smult * 2 * (fcgefhs + fcgefgs)));
        wavesign = -1;
        smult = 0;
        long double triplet = (0.5 * Gamov(kstar) * (2.0 * fgmodt * smult + norm(ffplus) + norm(ffminus) + wavesign*sterm.real() + wavesign*tterm.real() + smult * 2 * (fcgefht + fcgefgt)));
        return (0.25 * singlet + 0.75 * triplet);
    }
    else
        return (0.5 * Gamov(kstar) * (2.0 * fgmods * smult + norm(ffplus) + norm(ffminus) + wavesign*sterm.real() + wavesign*tterm.real() + smult * 2 * (fcgefhs + fcgefgs)));
}

// Calculate the wave-function modulus sqaured
// for identical bosons (symmetrized) or fermions (antisymmetrized)
// with Coulomb interaction and Strong interaction included
double FemtoInteraction::GetNewQuantumCoulombStrong(float tDiff)
{
    if (mRSt < 0.0000000001)
        return 1.0;

    if (mRSt < 1.0/THGlobal::GeVtoFm) // once fucking again, what are those arbitrary numbers? - JJ
        return GetQuantumCoulomb(tDiff);

    if (tDiff > 500.0) // does this make me a hypocryte? - JJ
        return 1.0;
    
    double tKstRst = mKStarOut*mROS + mKStarSide*mRSS + mKStarLong*mRLS;
    long double kstar = fabs(mKStarSigned);
    long double rho = mRSt * kstar;
    
    int ccase = 0;
    int wavesign = 1;

    // Classical limit - if distance is larger than Coulomb radius, 
    // the interaction does not matter
    //  if (fabs(mRSt) > fabs(pionac)) return (1.0 + wavesign*cos(2*tKstRst));

    // Classical limit - in the case of large k* we go to 
    // classical coulomb interaction
    long double testp = rho * (1.0 + tKstRst/(rho));
    long double testm = rho * (1.0 - tKstRst/(rho));

    std::complex<long double> ffplus, ffminus;
    if ((testp> 15.0) && (testm> 15.0))
    {
        double asym;
        asym = sqrt((1.0 - 1.0 / (mRSt * (1.0 - tKstRst / rho) * pionac * kstar * kstar)) / Gamov(kstar));
        ffminus = {(asym < 1.0 ? 1.0 + (asym -1.0) *2.0 : 1.0 + (asym -1.0) /2.0),sqrt(asym*asym - ffminus.real()*ffminus.real())};

        asym = sqrt((1.0 - 1.0 / (mRSt * (1.0 + tKstRst / rho) * pionac * kstar * kstar)) / Gamov(kstar));
        ffplus = {(asym < 1.0 ? 1.0 + (asym -1.0) *2.0 : 1.0 + (asym -1.0) /2.0),sqrt(asym*asym - ffplus.real()*ffplus.real())};
    }
    // Check for the classical limit in both functions separately
    else if (((testp< 15.0) && (testm< 15.0))) // ||
    {
        // Calculate the F function
        std::tie(ffplus,ffminus) = GetFFdouble();
        ccase = 1;
    }
    else if (testp < 15.0)
    {
        double asym;
        ffplus = GetFFsingle(1);
        ffminus = GetFFsingle(-1);
        if ((fabs(ffminus.real()) > 2.0) || fabs(ffminus.imag()) > 2.0) 
        {
            asym = sqrt((1.0 - 1.0 / (mRSt * (1.0 - tKstRst / (rho) * pionac * kstar * kstar))) / Gamov(kstar));
            ffminus = {(asym < 1.0 ? 1.0 + (asym -1.0) *2.0 : 1.0 + (asym -1.0) /2.0),sqrt(asym*asym - ffminus.real()*ffminus.real())};
        }
        ccase = 2;
    }
    else 
    {
        double asym;
        ffplus = GetFFsingle(1);
        ffminus = GetFFsingle(-1);
        if ((fabs(ffplus.real()) > 2.0) || fabs(ffplus.imag()) > 2.0) 
        {
            asym = sqrt((1.0 - 1.0 / (mRSt * (1.0 + tKstRst / (rho) * pionac * kstar * kstar))) / Gamov(kstar));
            ffplus = {(asym < 1.0 ? 1.0 + (asym -1.0) *2.0 : 1.0 + (asym -1.0) /2.0),sqrt(asym*asym - ffplus.real()*ffplus.real())};
        }
        ccase = 3;
    }

    long double eta  = 1.0/(kstar*pionac);
    long double hfun = GetH(eta);
    std::complex<long double> gtilde  = GetG(eta, rho, hfun) / (long double) mRSt;
    std::complex<long double> fcouls, fcoult;
    std::tie(fcouls, fcoult) = GetFC(kstar, eta, hfun);
    
    std::complex<long double> fgs = gtilde * fcouls;
    std::complex<long double> fgt = gtilde * fcoult;
    long double fgmods = abs(fgs)*abs(fgs);
    long double fgmodt = abs(fgt)*abs(fgt);

    std::complex<long double> expikr = {cos(tKstRst),-sin(tKstRst)};
    std::complex<long double> sterm = expikr * expikr * ffplus * conj(ffminus);
    std::complex<long double> tterm = conj(expikr * expikr) * conj(ffplus) * ffminus; // isn't this just conj(sterm) ? - JJ

    std::complex<long double> epfpc = expikr * ffplus;
    std::complex<long double> emfmc = conj(expikr) * ffminus;

    long double fcgefhs = (fgs.real()*emfmc.real() + fgs.imag()*emfmc.imag());
    long double fcgefgs = (fgs.real()*epfpc.real() + fgs.imag()*epfpc.imag());

    long double fcgefht = (fgt.real()*emfmc.real() + fgt.imag()*emfmc.imag());
    long double fcgefgt = (fgt.real()*epfpc.real() + fgt.imag()*epfpc.imag());

    long double smult = 1 + wavesign;

    if (!finite(ffminus.real()))
        std::cout << "FFMinus Re not a number ! " << testp << " " << testm << " " << ccase<< std::endl;
  
    if (!finite(ffminus.imag()))
        std::cout << "FFMinus Im not a number !" << testp << " " << testm << " " << ccase<< std::endl;

    if ((ffplus.real() > 2.0) || (ffplus.real() < -2.0))
        std::cout << "FFplus Re wild !" << ffplus.real() << " case " << ccase << " " << testp << " " << testm << std::endl;

    if ((ffplus.imag() > 2.0) || (ffplus.imag() < -2.0))
        std::cout << "FFplus Im wild !" << ffplus.imag() << " case " << ccase << " " << testp << " " << testm << std::endl;

    if ((ffminus.real() > 2.0) || (ffminus.real() < -2.0))
        std::cout << "FFminus Re wild !" << ffminus.real() << " case " << ccase << std::endl;
    
    if ((ffminus.imag() > 2.0) || (ffminus.imag() < -2.0))
        std::cout << "FFminus Im wild !" << ffminus.imag() << " case " << ccase << std::endl;

    if (twoSpin == 1) 
    {
        wavesign = 1;
        smult = 2;
        long double singlet = (0.5 * Gamov(kstar) * (2.0 * fgmods * smult + norm(ffplus) + norm(ffminus) + wavesign*sterm.real() + wavesign*tterm.real() + smult * 2 * (fcgefhs + fcgefgs)));
        wavesign = -1;
        smult = 0;
        long double triplet = (0.5 * Gamov(kstar) * (2.0 * fgmodt * smult + norm(ffplus) + norm(ffminus) + wavesign*sterm.real() + wavesign*tterm.real() + smult * 2 * (fcgefht + fcgefgt)));
        return (0.25 * singlet + 0.75 * triplet);
    }
    else
        return (0.5 * Gamov(kstar) * (2.0 * fgmods * smult + norm(ffplus) + norm(ffminus) + wavesign*sterm.real() + wavesign*tterm.real() + smult * 2 * (fcgefhs + fcgefgs)));
}

// Returns the std::tuple<bsum,psum>
std::tuple<long double, long double> FemtoInteraction::BFunPFun(long double eta, long double rho)
{
    long double bsum = 1. + eta * rho;
    long double bnpu, bn = eta * rho, bnmu = 1.;
    long double psum = 1.;
    long double pnpu, pn = 0., pnmu = 1.;

    if (rho > TMath::Pi()*4) 
        return std::make_tuple(sin(rho)/rho,cos(rho));

    for (int iter = 1; iter < 100000; iter++) 
    {
        bnpu = 2 * eta*rho *bn - rho*rho*bnmu;
        bnpu /= (1.0*iter+1.0)*(1.0*iter+2.0);
        bsum += bnpu;

        pnpu = 2 * eta*rho *pn - rho*rho*pnmu - (2.0*iter+1.0)*2.0*eta*rho*bn;
        pnpu /= (1.0*iter)*(1.0*iter+1.0);
        psum += pnpu;

        bnmu = bn;
        bn   = bnpu;

        pnmu = pn;
        pn   = pnpu;
        if ((fabs(pnmu) + fabs (bnmu) + fabs(bnpu) + fabs(pnpu)) < 1.0e-20) // whyyyyyyy this random value? - JJ
            break;
    }

    return std::make_tuple(bsum,psum);
}

// Returns std::tuple<fcs,fct>
std::tuple<std::complex<long double>,std::complex<long double>> FemtoInteraction::GetFC(long double kstar, long double eta, long double hfun)
{
    std::complex<long double> fcs,fct,cia;

    cia = {hfun,ChiIm(eta)};
    cia *= (long double) 2. / pionac;
    fcs = ((long double) 1. / f0s) + (d0s * (0.5 * kstar * kstar)) - cia;
    fct = ((long double) 1. / f0t) + (d0t * (0.5 * kstar * kstar)) - cia;

    return std::make_tuple(((long double) 1. / fcs),((long double) 1. / fct));
}

// Calculate the Gtilde function
std::complex<long double> FemtoInteraction::GetG(long double eta, long double rho, long double hfun)
{
    long double euler = 0.577215665;
    std::complex<long double> gtemp;
    long double bres, pres, bmult;

    std::tie(bres, pres) = BFunPFun(eta, rho);
    bmult = 2.0*eta*rho*bres;

    gtemp = {pres + bmult * (log(fabs(2.0 * eta * rho)) + 2.0 * euler - 1.0 + hfun), bmult * ChiIm(eta)};
    return gtemp;
}

// Calculates the h function for the strong interaction
long double FemtoInteraction::GetH(long double eta) // it looks weird, but ok I'll allow it - JJ
{
    long double euler = 0.577215665;
    long double etasum = log(1.0/eta) - euler;
    long double series = 0.0;
    long double x2inv = (eta*eta);
    long double element;
    long double save;

    for (int iter = 1; iter < 1000000; iter++) 
    {
        element = ((1.0*iter)*(1.0*iter) + x2inv)*(1.0*iter);
        element = 1.0/element;
        if (iter == 1) 
            save = 1.0e-10 * element;

        series += element;
        if (element < save) 
            break;
    }
    series *= x2inv;
    etasum += series;

    return etasum;
}

void FemtoInteraction::InitGamov() 
{
    twoSpin = 0;
    switch(pairtype) // all the numbers should be, in my opinion, constants - JJ
    {
        case 0:
            pionac = 387.5 / THGlobal::GeVtoFm;
            partpid = 211;
            partpid2 = 0;
            break;
        case 1:
            pionac = 109.55 / THGlobal::GeVtoFm;
            partpid = 321;
            partpid2 = 0;
            break;
        case 2:
            pionac = 57.63975274 / THGlobal::GeVtoFm;
            partpid = 2212;
            partpid2 = 0;
            twoSpin = 1;
            d0s = 2.77 / THGlobal::GeVtoFm;
            f0s = 7.77 / THGlobal::GeVtoFm;
            d0t = 1.7  / THGlobal::GeVtoFm;
            f0t = -5.4 / THGlobal::GeVtoFm;
            break;
        case 3:
            pionac = 248.519 / THGlobal::GeVtoFm;
            partpid = 221;
            partpid2 = 321;
            break;
        case 4:
            pionac = -248.519 / THGlobal::GeVtoFm;
            partpid = -211;
            partpid2 = 321;
            break;
        case 5:
            pionac = 222.564 / THGlobal::GeVtoFm;
            partpid = 211;
            partpid2 = 2212;
            break;
        case 6:
            pionac = -222.564 / THGlobal::GeVtoFm;
            partpid = -211;
            partpid2 = 2212;
            break;
        case 7:
            pionac = 83.59432006 / THGlobal::GeVtoFm;
            partpid = 321;
            partpid2 = 2212;
            break;
        default:
            std::cout << "Unknown pair type " << pairtype << std::endl; // this should be some standarised error message + it should not occur - JJ
            exit(0);
    }
/*
    long double oneoveracsq = 1.0/(pionac * pionac);
    long double twopioverac = 2.0*TMath::Pi()/pionac;
    double tpaoverk;
    double gamov[2000];
    
    for (int iter=0; iter<2000; iter++) 
    {
        tpaoverk = twopioverac/(iter*0.0002 + 0.0001);
        gamov[iter] = tpaoverk * 1.0 / (exp(tpaoverk) - 1);
    }*/ // leaving this commented for now, the Gamov mehtod does exactly this and this part is never used in the origina code - JJ
}

double FemtoInteraction::Gamov(double arg)
{
    long double eta = (2.0*TMath::Pi()/pionac) / arg;
    return (eta) * 1.0/(exp(eta) - 1.0);
}

// Calculates the confluent hypergeometric function F
// from single orientation of cos(theta*)
// For non-symmetrized wave-function (non-identical particles)
std::complex<long double> FemtoInteraction::GetFFsingle(int sign)
{
    long double eta, ksip, kstar, kstrst, coskr, tcomp;
    std::complex<long double> alfa, zetp, sump, fcomp, scompp;
    int nsteps;

    kstar = fabs(mKStarSigned);
    kstrst = mKStarOut * mROS + mKStarSide * mRSS + mKStarLong * mRLS;
    coskr = sign * kstrst / (kstar * mRSt);
    ksip = kstar * mRSt * (1 + coskr);

    if (ksip > 15.0)
        nsteps = 170;
    else if (ksip > 10.0)
        nsteps = 45;
    else if (ksip > 5.0)
        nsteps = 35;
    else
        nsteps = 15;

    eta = 1.0 / (kstar * pionac);
    alfa = {0.0,-eta};

    zetp = {0.0,ksip};
        
    sump = {1.0,0.0};
    tcomp = 1.0;
        
    for (int istep = 0; istep < nsteps; istep++) 
    {
        if (istep > 0) 
            sump += (fcomp * scompp) / (tcomp * tcomp);
        
        fcomp *= (alfa + (long double) istep);
        scompp *= zetp;
        tcomp *= (istep+1);

        if (norm(sump) < 1.0e-14) // why is this exactly 1.0e-14? - JJ
            break;
    }
    
    return sump;
}

// Calculates the confluent hypergeometric function
// For two orientations of cos(theta*) 
// For symmetrized wave-function (identical particles)
// Returns std::tuple<ffp,ffm>
std::tuple<std::complex<long double>,std::complex<long double>> FemtoInteraction::GetFFdouble()
{
    long double eta, ksip, ksim, kstar, kstrst, coskr, tcomp;
    std::complex<long double> alfa, zetp, zetm, sump, summ, fcomp, scompp, scompm;
    int nsteps;

    kstar = fabs(mKStarSigned);
    kstrst = mKStarOut*mROS + mKStarSide*mRSS + mKStarLong*mRLS;
    coskr = kstrst/(fabs(mKStarSigned) * mRSt);
    ksip = kstar*mRSt*(1 + coskr);
    ksim = kstar*mRSt*(1 - coskr);

    if ((ksip < 5.0) && (ksim < 5.0))
        nsteps = 25;
    else if ((ksip < 10.0) && (ksim < 10.0))
        nsteps = 45;
    else if ((ksip < 15.0) && (ksim < 15.0))
        nsteps = 110;
    else
        nsteps = 150;
    
    eta = 1.0/(kstar * pionac);
    alfa = {0.0,-eta};

    zetp = {0.0,ksip};
    zetm = {0.0,ksim};
        
    sump = {1.0,0.0};
    summ = {1.0,0.0};
    tcomp = 1.0;
        
    for (int istep=0; istep<nsteps; istep++) 
    { 
        if (istep > 0) 
        {
            sump += (fcomp * scompp) / (tcomp * tcomp);
            summ += (fcomp * scompm) / (tcomp * tcomp);
        }
        
        fcomp *= (alfa + (long double) istep);
        scompp *= zetp;
        scompm *= zetm;
        tcomp *= (istep+1);
    }
    
    return std::make_tuple(sump,summ);
}