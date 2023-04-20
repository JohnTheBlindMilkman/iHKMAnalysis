#include "FileIO.h"

FileIO::FileIO(TString defIniFile, TString tmpFile) : mPathToFemtoMixer(gSystem->pwd()), mTmpFile(tmpFile), mROOTFileIn("0.root"), mROOTFileOut("femto.root"), mIniFile(defIniFile), mParser(mIniFile)
{
    
}

FileIO::~FileIO()
{

}

void FileIO::ReadIniFile()
{
    TString tmp;
    Configurator tMainConfig;
    tMainConfig = mParser.ReadINI();
    try
    {
        tmp = tMainConfig.GetParameter("PairType");
        if (tmp == "pion-pion")
            mPairType = 0; 
        else if (tmp == "kaon-kaon")
            mPairType = 1; 
        else if (tmp == "proton-proton")
            mPairType = 2; 
        {
            PRINT_MESSAGE("therm2_femto Unknown pair type: " << tmp);
            PRINT_MESSAGE("Please provide the proper pair name in the main INI file.");
            exit(THGlobal::Error::femtoUnknownPairType);
        }

        mTimeCut = tMainConfig.GetParameter("TimeCut").Atof();
        mKtMin = tMainConfig.GetParameter("Ktmin").Atof();
        mKtMax = tMainConfig.GetParameter("Ktmax").Atof();

        mEventsToMix = tMainConfig.GetParameter("EventsToMix").Atoi();
        mEventDir = tMainConfig.GetParameter("InputDir");
        mEventFiles = tMainConfig.GetParameter("EventFiles").Atoi();
        mLogFile = tMainConfig.GetParameter("LogFile");

        (tMainConfig.GetParameter("EnableOnlyPrimordial") == "yes") ? primordial = true : primordial = false;
        (tMainConfig.GetParameter("EnableCoulombCorrection") == "yes") ? coulomb = true : coulomb = false;
        (tMainConfig.GetParameter("EnableSourceHistograms") == "yes") ? writeHist = true : writeHist = false;

    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
}

TChain* FileIO::ReadROOTFile()
{
    TChain *chnEv = new TChain(THGlobal::sTreeName);
    
    for(int i = 0; i < mEventFiles; i++) 
    {
        char Buff[THGlobal::kFileNameMaxChar]; // this could be a TString... - JJ
        sprintf(Buff,"%s%d.root",mEventDir.Data(),i);
        PRINT_DEBUG_1("Adding file: " << Buff);
        chnEv->Add(Buff);
    }

    return chnEv;
}

void FileIO::WriteLogFile(std::vector<TString> fileContent)
{

}

void FileIO::WriteTmpFile(std::vector<TString> fileContent)
{

}

void FileIO::WriteROOTFile(TString fileName, std::vector<std::unique_ptr<TObject>> objList)
{

}

void FileIO::CopyINIFile()
{
    TString  tINI;
    std::ifstream ifs;
    std::ofstream ofs;
    
    tINI = mEventDir + mIniFile;
    ifs.open(mIniFile, std::ios::binary);
    ofs.open(tINI, std::ios::binary);
    if((ifs) && (ofs) && ifs.is_open() && ofs.is_open()) 
    {
        ofs << ifs.rdbuf();
        ifs.close();
        ofs.close();
    } 
    else
        PRINT_MESSAGE("<therm2_femto::main()>\tUnable to copy "<<mIniFile<<" to "<<tINI); 
}