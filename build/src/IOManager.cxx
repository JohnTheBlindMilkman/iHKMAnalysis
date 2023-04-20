#include "IOManager.h"

IOManager::IOManager(TString defIniFile, TString tmpFile) : mFileIO(defIniFile, tmpFile)
{
}

IOManager::~IOManager()
{
}

void IOManager::Init()
{
    mFileIO.ReadIniFile();
    //reader = TTreeReader(mFileIO.ReadROOTFile());
    // create all needed histograms and send them to here
    // pass the reader to Mixing Manager ? 
}

void IOManager::Close()
{
    // save histograms
    // white log file
    // copy ini file
    // run destructors I guess ?
}

std::tuple<TString,int> IOManager::Parse(int argc, char *argv[])
{
    TString iniFile = argv[1];
    int sParentPID = atoi(argv[2]);

    return std::make_tuple(iniFile,sParentPID);
}