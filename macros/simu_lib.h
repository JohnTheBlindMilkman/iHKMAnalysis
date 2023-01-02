#include <iostream>
#include <random>

#include "TBenchmark.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TVirtualPad.h"
#include <math.h>

#define PIPID 211
#define KPID 321
#define PPID 2212
#define ETAABS 0.35
#define PTMIN 0.12
#define PTMAX 0.9
#define TMAX 500
#define ABSRAP 1.0// temporary value that i used once

#define COULOMBSTEPS 170

#define KAONMASS 0.493677
#define PIONMASS 0.13957
#define PROTONMASS 0.93827

#define FOURPI 1.25663706143591725e+01

using namespace std;

//complex *************************************************************
struct _dcomplex 
{
  long double re;
  long double im;
};

typedef struct _dcomplex dcomplex;

dcomplex mult(dcomplex arga, dcomplex argb)
{
  dcomplex res;
  
  res.re = arga.re * argb.re - arga.im * argb.im;
  res.im = arga.re * argb.im + argb.re * arga.im;

  return res;
}
dcomplex conj(dcomplex arg)
{
  dcomplex res;
  
  res.re = arg.re;
  res.im = -arg.im;

  return res;
}

dcomplex mult(dcomplex arga, long double argb)
{
  dcomplex res;
  
  res.re = arga.re * argb;
  res.im = arga.im * argb;

  return res;
}

long double modl2(dcomplex arg)
{
  return arg.re*arg.re + arg.im*arg.im;
}

long double modl(dcomplex arg)
{
  return hypot(arg.re, arg.im);
}

dcomplex invr(dcomplex arg)
{
  dcomplex res;
  
  res.re = arg.re/modl2(arg);
  res.im = - arg.im/modl2(arg);

  return res;
}

//end complex*************************************************************



void		InitializeGamow();
double		Gamow(double arg);
dcomplex	mult(dcomplex arga, dcomplex argb);
dcomplex	conj(dcomplex arg);
dcomplex	mult(dcomplex arga, long double argb);
long double	modl2(dcomplex arg);
long double	modl(dcomplex arg);
dcomplex	invr(dcomplex arg);
void		GetFFsingle(dcomplex *ffp, int sign = 1);
void		GetFFdouble(dcomplex *ffp, dcomplex *ffm);
double		GetQuantumCoulomb();
double		GetCoulomb();

double mRL, mRS, mRO, mDT, mKO, mKS, mKL, mDE, mROS, mRSS, mRLS, mRSt;
double mKStarLong, mKStarOut, mKStarSide, mRSidePairCMS, mRLongPairCMS, mDTimePairLCMS, mROutPairCMS, mDTimePairCMS, mRLong, mROut, mRSide, mDTime, mRTrans, mBetat;
double mRStar, mKStarSigned;
  
dcomplex d0s;
dcomplex f0s;
dcomplex d0t;
dcomplex f0t;

int part1pid = PIPID;
int part2pid = PIPID;

int pairtype = 0;
double pionac = 387.5 / 0.197327;
int twospin = 0;// 0 for pion-pion, kaon-kaon, 1 for proton-proton

int pcount = 0;  
double gamov[2000];
long double oneoveracsq;
long double twopioverac;
long double coulqscpart;  


void InitializeGamow()
{
  twospin = 0;
  if (pairtype == 0) {
    pionac = 387.5 / 0.197327;
    part1pid = PIPID;
    part2pid = PIPID;
  }
  else if (pairtype == 1) {
    pionac = 109.55 / 0.197327;
    part1pid = KPID;
    part2pid = KPID;
  }
  else if (pairtype == 2) {
    pionac = 57.63975274 / 0.197327;
    part1pid = PPID;
    part2pid = PPID;
    twospin = 1;
    d0s.re = 2.77 / 0.197327;
    f0s.re = 7.77 / 0.197327;
    d0t.re = 1.7  / 0.197327;
    f0t.re = -5.4 / 0.197327;
  }
  else if (pairtype == 3) {
    pionac = 248.519 / 0.197327;
    part1pid = PIPID;
    part2pid = KPID;
  }
  else if (pairtype == 4) {
    pionac = -248.519 / 0.197327;
    part1pid = -PIPID;
    part2pid = KPID;
  }
  else if (pairtype == 5) {
    pionac = 222.564 / 0.197327;
    part1pid = PIPID;
    part2pid = PPID;
  }
  else if (pairtype == 6) {
    pionac = -222.564 / 0.197327;
    part1pid = -PIPID;
    part2pid = PPID;
  }
  else if (pairtype == 7) {
    pionac = 83.59432006 / 0.197327;
    part1pid = KPID;
    part2pid = PPID;
  }
  else {
    
    cout << "Unknown pair type " << pairtype << endl;
    exit(0);
  }

  oneoveracsq = 1.0/(pionac * pionac);
  twopioverac = 2.0*TMath::Pi()/pionac;
  double tpaoverk;
  
  for (int iter=0; iter<2000; iter++) {
    tpaoverk = twopioverac/(iter*0.0002 + 0.0001);
    gamov[iter] = tpaoverk * 1.0 / (exp(tpaoverk) - 1);
  }
}

double Gamow(double arg)
{
  long double eta = twopioverac / arg;
  return (eta) * 1.0/(exp(eta) - 1.0);
}

// Calculates the confluent hypergeometric function F
// from single orientation of cos(theta*)
// For non-symmetrized wave-function (non-identical particles)
void GetFFsingle(dcomplex *ffp, int sign)
{
  double comprep[COULOMBSTEPS];
  double compimp[COULOMBSTEPS];
  double eta, ksip;
  dcomplex alfa, zetp;
  
  int nsteps;

  double kstar = fabs(mKStarSigned);
  double kstrst = mKStarOut*mROS + mKStarSide*mRSS + mKStarLong*mRLS;
  double coskr = sign * kstrst/(fabs(mKStarSigned) * mRSt);

  if (kstar*mRSt*(1+coskr) > 15.0)
    nsteps = 170;
  else if (kstar*mRSt*(1+coskr) > 10.0)
    nsteps = 45;
  else if (kstar*mRSt*(1+coskr) > 5.0)
    nsteps = 35;
  else
    nsteps = 15;

  eta = 1.0/(kstar * pionac);
  alfa.re = 0.0;
  alfa.im = -eta;

  dcomplex fcomp, scompp;
  double tcomp;
  dcomplex sump;
  dcomplex fcmult;

  double rad = mRSt;

  ksip = kstar*rad*(1+coskr);

  zetp.re = 0.0;
  zetp.im = ksip;
      
  fcomp.re = 1.0;
  fcomp.im = 0.0;
  scompp.re = 1.0; 
  scompp.im = 0.0;
  tcomp = 1.0;
      
  for (int istep=0; istep<nsteps; istep++) {
    sump = mult(fcomp, scompp);

    sump = mult(sump, 1.0/(tcomp*tcomp));
	
    if (istep == 0) {
      comprep[istep] = sump.re;
      compimp[istep] = sump.im;
    }
    else {
      comprep[istep] = comprep[istep-1] + sump.re;
      compimp[istep] = compimp[istep-1] + sump.im;
    }
    
    fcmult.re = alfa.re + istep;
    fcmult.im = alfa.im;
	
    fcomp = mult(fcomp, fcmult);
    scompp = mult(scompp, zetp);
    tcomp *= (istep+1);

    if ((sump.re*sump.re+sump.im*sump.im) < 1.0e-14) {
      nsteps = istep;
      break;
    }
  }
  
  ffp->re = comprep[nsteps-1];
  ffp->im = compimp[nsteps-1]; 
}

// Calculates the confluent hypergeometric function
// For two orientations of cos(theta*) 
// For symmetrized wave-function (identical particles)
void GetFFdouble(dcomplex *ffp, dcomplex *ffm)
{
  long double comprep[COULOMBSTEPS];
  long double compimp[COULOMBSTEPS];
  long double comprem[COULOMBSTEPS];
  long double compimm[COULOMBSTEPS];
  long double eta, ksip, ksim;
  dcomplex alfa, zetp, zetm;
  
  int nsteps;

  long double kstar = fabs(mKStarSigned);
  long double kstrst = mKStarOut*mROS + mKStarSide*mRSS + mKStarLong*mRLS;
  long double coskr = kstrst/(fabs(mKStarSigned) * mRSt);

  if ((kstar*mRSt*(1+coskr) < 5.0) &&
      (kstar*mRSt*(1-coskr) < 5.0))
    nsteps = 25;
  else if ((kstar*mRSt*(1+coskr) < 10.0) &&
	   (kstar*mRSt*(1-coskr) < 10.0))
    nsteps = 45;
  else if ((kstar*mRSt*(1+coskr) < 15.0) &&
	   (kstar*mRSt*(1-coskr) < 15.0))
    nsteps = 110;
  else
    nsteps = 150;
  
  eta = 1.0/(kstar * pionac);
  alfa.re = 0.0;
  alfa.im = -eta;

  dcomplex fcomp, scompp, scompm;
  long double tcomp;
  dcomplex sump, summ;
  dcomplex fcmult;

  long double rad = mRSt;

  ksip = kstar*rad*(1+coskr);
  ksim = kstar*rad*(1-coskr);

  zetp.re = 0.0;
  zetp.im = ksip;
      
  zetm.re = 0.0;
  zetm.im = ksim;

  fcomp.re = 1.0;
  fcomp.im = 0.0;
  scompp.re = 1.0; 
  scompp.im = 0.0;
  scompm.re = 1.0; 
  scompm.im = 0.0;
  tcomp = 1.0;
      
  for (int istep=0; istep<nsteps; istep++) {
    sump = mult(fcomp, scompp);
    summ = mult(fcomp, scompm);

    sump = mult(sump, 1.0/(tcomp*tcomp));
    summ = mult(summ, 1.0/(tcomp*tcomp));
	
	
    if (istep == 0) {
      comprep[istep] = sump.re;
      compimp[istep] = sump.im;
      
      comprem[istep] = summ.re;
      compimm[istep] = summ.im;
    }
    else {
      comprep[istep] = comprep[istep-1] + sump.re;
      compimp[istep] = compimp[istep-1] + sump.im;
      
      comprem[istep] = comprem[istep-1] + summ.re;
      compimm[istep] = compimm[istep-1] + summ.im;
    }
    
    fcmult.re = alfa.re + istep;
    fcmult.im = alfa.im;
	
    fcomp = mult(fcomp, fcmult);
    scompp = mult(scompp, zetp);
    scompm = mult(scompm, zetm);
    tcomp *= (istep+1);
  }
  
  ffp->re = comprep[nsteps-1];
  ffp->im = compimp[nsteps-1];

  ffm->re = comprem[nsteps-1];
  ffm->im = compimm[nsteps-1];
}



double GetCoulomb()
{
  double kstrst = mKStarOut*mROS + mKStarSide*mRSS + mKStarLong*mRLS;

  // Classical limit - if distance is larger than Coulomb radius, 
  // the interaction does not matter
  if (fabs(mRSt) > fabs(pionac)) return (1.0);
  if (fabs(mRSt) == 0.0){
      if (Gamow(fabs(mKStarSigned)) > 0.001)
	  return (Gamow(fabs(mKStarSigned)));
  }

  // Classical limit - in the case of large k*r* product we go to 
  // classical coulomb interaction
  if (fabs(mKStarSigned) * mRSt * (1.0 + kstrst/(mRSt*fabs(mKStarSigned)))> 15.0)
    return (1.0 - 1.0/(mRSt*pionac*mKStarSigned*mKStarSigned));
  
  // Calculate the F function
  dcomplex ffplus;
  GetFFsingle(&ffplus);

  if(modl2(ffplus)<3)
  	return (Gamow(fabs(mKStarSigned)) * modl2(ffplus));
  else
  	return 0;
  	
}

double GetQuantumCoulomb()
{
  if (mRSt < 0.0000000001)
    return 1.0;

  double kstrst = mKStarOut*mROS + mKStarSide*mRSS + mKStarLong*mRLS;
  int ccase = 0;
  static int pcount = 0;
  int wavesign = 1;

  if (twospin == 1) {
    if (pcount == 3)
      wavesign = 1;
    else
      wavesign = -1;
    pcount++;
    if (pcount == 4) pcount = 0;
  }

  if (fabs(mRSt) > fabs(pionac)){
      return (1.0 + wavesign*cos(2*kstrst));
  }
  
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

  dcomplex ffplus, ffminus;
  // Check for the classical limit in both functions separately
  if (((testp< 15.0) && (testm< 15.0))) // ||
  {
     // Calculate the F function
     GetFFdouble(&ffplus, &ffminus);
     ccase = 1;
  }
  else if (testp< 15.0)
  {
     double asym;
     GetFFsingle(&ffplus);
     asym = (1.0 - 1.0/(mRSt*(1.0 - kstrst/(mRSt*fabs(mKStarSigned))*pionac*mKStarSigned*mKStarSigned)))/Gamow(fabs(mKStarSigned));
     asym = sqrt(asym);
     if (asym < 1.0) 
	ffminus.re = 1.0 + (asym -1.0) *2.0;
     else
	ffminus.re = 1.0 + (asym -1.0) /2.0;
     ffminus.im = sqrt(asym*asym - ffminus.re*ffminus.re);
     ccase = 2;
  }
  else 
  {
      double asym;
      GetFFsingle(&ffminus, -1);
      asym = (1.0 - 1.0/(mRSt*(1.0 + kstrst/(mRSt*fabs(mKStarSigned))*pionac*mKStarSigned*mKStarSigned)))/Gamow(fabs(mKStarSigned));
      asym = sqrt(asym);
      if (asym < 1.0) 
	ffplus.re = 1.0 + (asym -1.0) *2.0;
      else
	ffplus.re = 1.0 + (asym -1.0) /2.0;
      ffplus.im = sqrt(asym*asym - ffplus.re*ffplus.re);
      ccase = 3;

   }

  dcomplex expikr;
  expikr.re = cos(kstrst);
  expikr.im = sin(kstrst);
  dcomplex expikrc = conj(expikr);
  dcomplex ffplusc = conj(ffplus);
  dcomplex ffminusc = conj(ffminus);

  dcomplex expikr2 = mult(expikr, expikr);
  dcomplex expikrc2 = conj(expikr2);
  dcomplex sterm = mult(expikr2,  mult(ffplus, ffminusc));
  dcomplex tterm = mult(expikrc2, mult(ffminus, ffplusc));

  if((modl2(ffminus)<3) && (modl2(ffplus)<3))
     return (0.5 * Gamow(fabs(mKStarSigned)) * 
	  (modl2(ffplus) + wavesign*sterm.re + wavesign*tterm.re + modl2(ffminus)));
  else
     return 0;
}


void PairKinematics(int ii, int jj, TTreeReaderArray<float> &px, TTreeReaderArray<float> &py, TTreeReaderArray<float> &pz, TTreeReaderArray<float> &E, TTreeReaderArray<float> &t, TTreeReaderArray<float> &x, TTreeReaderArray<float> &y, TTreeReaderArray<float> &z){

  // Calculate pair variables
  double tPx = px[ii]+px[jj];
  double tPy = py[ii]+py[jj];
  double tPz = pz[ii]+pz[jj];
  double tE  = E[ii] +E[jj];
  
  double tPt = tPx*tPx + tPy*tPy;
  double tMt = tE*tE - tPz*tPz;//mCVK;
  double tM  = sqrt(tMt - tPt);
  
  tMt = sqrt(tMt);
  tPt = sqrt(tPt);
  
  double mKT = tPt/2.0;
  if ((mKT<0.15) || (mKT>0.25))
    return;

  double mBetat = tPt/tMt;  
  if ((mBetat<0.0) || (mBetat>1.0)) {
    mKT = -1.0;
    return;
  }

  double tbm = tPz/tE;
  double tgm = tE/tMt;
      
  double mkl = tgm * (pz[ii] - tbm*E[ii]);
  double met = tgm * (E[ii] - tbm*pz[ii]);

  double mko = ( px[ii]*tPx + py[ii]*tPy)/tPt;
  double mks = (-px[ii]*tPy + py[ii]*tPx)/tPt;

  mko = (tMt/tM) * (mko - (tPt/tMt) * met);

  double mkv = sqrt(mko*mko + mks*mks + mkl*mkl);
  mKStarSigned=mkv*0.5;
  
  
  // Boost to LCMS
  double tBeta = tPz/tE;
  double tGamma = tE/tMt;	    
  mKStarLong = tGamma * (pz[ii] - tBeta * E[ii]);
  double tE1L = tGamma * (E[ii]  - tBeta * pz[ii]);
  
  // Transform to LCMS
  double particle1lcms_z = tGamma * (z[ii] - tBeta * t[ii]);
  double particle1lcms_t = tGamma * (t[ii] - tBeta * z[ii]);
  
  double particle2lcms_z = tGamma * (z[jj] - tBeta * t[jj]);
  double particle2lcms_t = tGamma * (t[jj] - tBeta * z[jj]);
  
  double particle1prf_pz = tGamma * (pz[ii] - tBeta * E[ii]);
  double particle1lcms_pz = tGamma * (pz[ii] - tBeta * E[ii]);
  double particle1lcms_e  = tGamma * (E[ii]  - tBeta * pz[ii]);
  
  double particle2prf_pz = tGamma * (pz[jj] - tBeta * E[jj]);
  double particle2lcms_pz = tGamma * (pz[jj] - tBeta * E[jj]);
  double particle2lcms_e  = tGamma * (E[jj]  - tBeta * pz[jj]);
  
  // Rotate in transverse plane
  mKStarOut  = ( px[ii]*tPx + py[ii]*tPy)/tPt;
  mKStarSide = (-px[ii]*tPy + py[ii]*tPx)/tPt;
      
  double particle1lcms_px = mKStarOut;
  double particle1lcms_py = mKStarSide;
  
  double particle2lcms_px = (px[jj]*tPx + py[jj]*tPy)/tPt;
  double particle2lcms_py = (py[jj]*tPx - px[jj]*tPy)/tPt;;
  
  mKO = particle1lcms_px - particle2lcms_px;
  mKS = particle1lcms_py - particle2lcms_py;
  mKL = particle1lcms_pz - particle2lcms_pz;
  mDE = particle1lcms_e  - particle2lcms_e;
  
  // save the rotated coordinates in LCMS variables
  double particle1lcms_x = ( x[ii]*tPx + y[ii]*tPy)/tPt;
  double particle1lcms_y = (-x[ii]*tPy + y[ii]*tPx)/tPt;

  double particle2lcms_x = ( x[jj]*tPx + y[jj]*tPy)/tPt;
  double particle2lcms_y = (-x[jj]*tPy + y[jj]*tPx)/tPt;
  
  // Boost to pair cms
  mKStarOut = tMt/tM * (mKStarOut - tPt/tMt * tE1L);
  
  Double_t tBetat = tPt/tMt;
  Double_t tGammat = 1.0/TMath::Sqrt(1.0-tBetat*tBetat);
  
  double particle1prf_x = tGammat*(particle1lcms_x - tBetat*particle1lcms_t);
  double particle1prf_t = tGammat*(particle1lcms_t - tBetat*particle1lcms_x);
  
  double particle2prf_x = tGammat*(particle2lcms_x - tBetat*particle2lcms_t);
  double particle2prf_t = tGammat*(particle2lcms_t - tBetat*particle2lcms_x);
  
  mRO = (particle1lcms_x - particle2lcms_x)/0.197327;
  mRS = (particle1lcms_y - particle2lcms_y)/0.197327;
  mRL = (particle1lcms_z - particle2lcms_z)/0.197327;
  mDT = (particle1lcms_t - particle2lcms_t)/0.197327;
 
  mKO = particle1lcms_px - particle2lcms_px;
  
  double tDX = x[ii]-x[jj];
  double tDY = y[ii]-y[jj];
  mRLong = z[ii]-z[jj];
  mDTime = t[ii]-t[jj];

  mRTrans = tDX>0.? ::sqrt(tDX*tDX+tDY*tDY) : -1.*::sqrt(tDX*tDX+tDY*tDY);
  mROut = (tDX*tPx + tDY*tPy)/tPt;
  mRSide = (-tDX*tPy + tDY*tPx)/tPt;

  mRSidePairCMS = mRSide;
  mRSS = mRSidePairCMS/0.197327;

  mRLongPairCMS = tGamma*(mRLong - tBeta* mDTime);
  mDTimePairLCMS = tGamma*(mDTime - tBeta* mRLong);

  mRLS = mRLongPairCMS/0.197327;
  tBeta = tPt/tMt;
  tGamma = tMt/tM;

  mROutPairCMS = tGamma*(mROut - tBeta* mDTimePairLCMS);
  mROS = mROutPairCMS/0.197327;
  mDTimePairCMS = tGamma*(mDTimePairLCMS - tBeta* mROut);

  mRStar = ::sqrt(mROutPairCMS*mROutPairCMS + mRSidePairCMS*mRSidePairCMS +
		  mRLongPairCMS*mRLongPairCMS);
  mRSt = mRStar/0.197327;  
  
}



void PairKinematics2(int ii, int jj, TTreeReaderArray<float> &px, TTreeReaderArray<float> &py, TTreeReaderArray<float> &pz, TTreeReaderArray<float> &E, TTreeReaderArray<float> &pxM, TTreeReaderArray<float> &pyM, TTreeReaderArray<float> &pzM, TTreeReaderArray<float> &EM){
  // Calculate pair variables
  double tPx = px[ii]+pxM[jj];
  double tPy = py[ii]+pyM[jj];
  double tPz = pz[ii]+pzM[jj];

  double tE  = E[ii] +EM[jj];
  double tPt = tPx*tPx + tPy*tPy;
 
  double tMt = tE*tE - tPz*tPz;//mCVK;
  double tM  = sqrt(tMt - tPt);
  tMt = sqrt(tMt);
  tPt = sqrt(tPt);
  
  double tbm = tPz/tE;
  double tgm = tE/tMt;
      
  double mkl = tgm * (pz[ii] - tbm*E[ii]);
  double met = tgm * (E[ii] - tbm*pz[ii]);

  double mko = ( px[ii]*tPx + py[ii]*tPy)/tPt;
  double mks = (-px[ii]*tPy + py[ii]*tPx)/tPt;

  mko = (tMt/tM) * (mko - (tPt/tMt) * met);

  double mkv = sqrt(mko*mko + mks*mks + mkl*mkl);
  mKStarSigned=mkv*0.5;
}



TH1D *getCF(TH1D *hnum, TH1D *hden, TString name, Float_t normmin, Float_t normmax)
{

    Int_t binmin = hnum->FindBin(normmin);
    Int_t binmax = hnum->FindBin(normmax);
    Double_t intN = hnum->Integral(binmin, binmax);
    Double_t intD = hden->Integral(binmin, binmax);
    Double_t normalization = intN / intD;
    TH1D *hdensc = (TH1D *)hden->Clone("hdensc");
    hdensc->Scale(normalization);
    TH1D *h = (TH1D *)hnum->Clone(name);
    h->Divide(hdensc);
    return h;
};


