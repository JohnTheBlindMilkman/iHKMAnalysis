#include "MixingManager.h"

MixingManager::MixingManager()
{
}

MixingManager::MixingManager(int argc, char *argv[])
{
    std::tie(mIniFile,mPpid) = IOManager::Parse(argc,argv);
    //mIOManager = std::unique_ptr<IOManager>(new IOManager(mIniFile,"femto_mixer.tmp"));
    mIOManager = std::make_unique<IOManager>(mIniFile,"femto_mixer.tmp");
}

MixingManager::~MixingManager()
{
}

void MixingManager::Init()
{
    // get necessary parameters from ini file
    // get reader
    /*mCoul = mIOManager->GetIniParamter(THGlobal::EuseCoulomb);
    mPrim = mIOManager->GetIniParamter(THGlobal::EisPrimordial);
    mPairType = mIOManager->GetIniParamter(THGlobal::EpairType);
    mEvtsToMix = mIOManager->GetIniParamter(THGlobal::EeventToMix);*/
    mCoul = true;
    mPrim = false;
    mPairType = 2;
    mEvtsToMix = 50;
    mKtMin = 0.2;
    mKtMax = 2.0;
}

void MixingManager::Start()
{
    // do the mixing process
    // fill histograms
    const int noFiles = 10;
    int evNum = 0;
    TString inpDir = "/home/jedkol/lustre/hades/user/kjedrzej/iHKM/14p5GeV/0to10cent/";
    TH1D *num1dqsc = new TH1D("num1dqsc", "num1dqsc;2k* [GeV/c];C(k*)", 200, 0.0, 0.4);
    TH1D *den1d = new TH1D("den1d", "den1d;2k* [GeV/c];C(k*)", 200, 0.0, 0.4);
    TH1D *wght = new TH1D("wght","wght",50,0,5);

    TChain *evtCh = new TChain("treefin");
    for (int i = 0; i < noFiles; i++)
        evtCh->Add(Form("%s%d.root",inpDir.Data(),i));

    TTreeReader reader(evtCh);
    std::vector<THGlobal::ParticleCoor> part;
    std::vector<std::array<double,2>> weights;
    THGlobal::ParticleCoor particle;
    TTreeReaderValue<int> npart(reader,"npart");
    TTreeReaderArray<float> px(reader,"px");
    TTreeReaderArray<float> py(reader,"py");
    TTreeReaderArray<float> pz(reader,"pz");
    TTreeReaderArray<float> E(reader,"E");
    TTreeReaderArray<int> id(reader,"id");
    TTreeReaderArray<int> mid(reader,"mid");
    TTreeReaderArray<float> t(reader,"t");
    TTreeReaderArray<float> x(reader,"x");
    TTreeReaderArray<float> y(reader,"y");
    TTreeReaderArray<float> z(reader,"z");
    while (reader.Next())
    {
        for (int iter = 0; iter < *npart; iter++) 
        {
            particle = {px[iter],py[iter],pz[iter],E[iter],x[iter],y[iter],z[iter],t[iter],id[iter],mid[iter]};
            if(TrackCut(particle))
                part.push_back(particle);
        }
        weights = Mix(part);
        AppendToHistogram(weights,num1dqsc,den1d,wght);

        part.clear();
        part.resize(0);
        weights.clear();
        weights.resize(0);
        evNum++;

        std::cout << "Event:\t" << evNum << "/" << 600 * noFiles << "\r";
        std::cout.flush();
    }

    TFile *fOut = TFile::Open("test.root","recreate");
    num1dqsc->Write();
    den1d->Write();
    wght->Write();
    fOut->Close();

    delete num1dqsc;
    delete den1d;
    delete evtCh;
    delete fOut;
}

void MixingManager::Finish()
{
    // call deconstructors ?
    // inform IOManager to finish
}

bool MixingManager::TrackCut(const THGlobal::ParticleCoor &part)
{
    if (part.id != 2212) // if not proton
        return false;

    float rap = (0.5 * log((part.e + part.pz) / (part.e - part.pz))); 
    if (abs(rap) >= 0.7) // if rapidity not good
        return false;

    float pt = hypot(part.px,part.py);
    if (pt <= 0.4 || pt >= 1.8) // if transverse momentum not good
        return false;

    float peta = -log(tan(atan2(pt,part.pz)/2.));
    if (peta >= 0.35) // if pseudorapidity (?) not good
        return false;

    return true;
}

bool MixingManager::PairCut(const THGlobal::ParticleCoor &part1, const THGlobal::ParticleCoor &part2)
{
    float kT = sqrt(pow(part1.px + part2.px,2) + pow(part1.py + part2.py,2)) / 2.;
    if (kT < mKtMin || kT > mKtMax)
        return false;

    return true;
}

std::vector<std::array<double,2>> MixingManager::Mix(const std::vector<THGlobal::ParticleCoor> &partList) // To-Do: pass weights by reference to omit copying data? - JJ
{
    std::vector<std::array<double,2>> weights;
    std::array<double,2> tmp = {0.,0.};

    for (std::size_t it = partList.size(); it > 0; it--)
    {
        for (std::size_t jt = 0; jt < it; jt++)
        {
            if (!PairCut(partList[it],partList[jt]))
                continue;

            tmp[0] = PairKinematics(partList[it],partList[jt],THGlobal::LCMS,mKtMin,mKtMax,0.,1.);
            if (tmp[0] > 0.)
            {
                tmp[1] = GetQuantumCoulombStrong(fabs(partList[it].t - partList[jt].t));
                weights.push_back(tmp);
            }
        }
    }

    return weights;
}

void MixingManager::AppendToHistogram(const std::vector<std::array<double,2>> &vecOfArr, TH1D *hNum, TH1D *hDen, TH1D *hWeight)
{
    if (vecOfArr.size() > mEvtsToMix)
        for (auto iter : vecOfArr)
        {
            hNum->Fill(fabs(iter[0])*2,iter[1]);
            hDen->Fill(fabs(iter[0])*2,1.0);
            hWeight->Fill(iter[1]);
        }
    else    
        return;
}
/*
// this would be more universal, but needs general refinement, i.e. filtr can't use the same obj that is called in partList
long int MixingManager::Mix(const std::vector<THGlobal::ParticleCoor> &partList, std::function<bool(THGlobal::ParticleCoor,THGlobal::ParticleCoor)> filtr)
{
    long int mixedPairs = 0;
    for (std::size_t it = partList.size(); it > 0; it--)
    {
        for (std::size_t jt = 0; jt < it; jt++)
        {
            if (!filtr(partList[it],partList[jt]))
                break;

            PairKinematics(partList[it],partList[jt],THGlobal::LCMS,mKtMin,mKtMax,0.,1.);
            mixedPairs++;
        }
    }

    return mixedPairs;
}*/