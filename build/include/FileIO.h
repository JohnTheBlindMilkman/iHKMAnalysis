#ifndef _FILE_IO_
    #define _FILE_IO_

    #include <iostream>
    #include <fstream>
    #include <TChain.h>
    #include <TFile.h>
    #include <TSystem.h>

    #include "Configurator.h"
    #include "Parser.h"
    #include "THGlobal.h"

    class FileIO
    {
        private:
            void CopyINIFile();

            TFile mInpFile, mOtpFile;
            TString mIniFile, mROOTFileIn, mROOTFileOut, mTmpFile, mLogFile, mPathToFemtoMixer, mEventDir;
            int mPairType, mEventsToMix, mEventFiles;
            float mTimeCut, mKtMin, mKtMax;
            bool primordial,coulomb,writeHist;
            std::map<TString, int> mIniParams; 
            std::map<TString, float> mkTParams;
            Parser mParser;

        public:
            FileIO(TString defIniFile, TString tmpFile);
            ~FileIO();
            void ReadIniFile();
            TChain* ReadROOTFile();
            void WriteTmpFile(std::vector<TString> fileContent);
            void WriteLogFile(std::vector<TString> fileContent);
            void WriteROOTFile(TString fileName, std::vector<std::unique_ptr<TObject>> objList);
            bool GetPrimordial();
            bool GetCoulomb();
            bool GetHistWriteStatus();
            int GetPairType();
            int GetEventsToMix();
            int GetEventFiles();
            float GetTimeCut();
            float GetMinKt();
            float GetMaxKt();
    };

    inline bool FileIO::GetPrimordial() {return primordial;}
    inline bool FileIO::GetCoulomb() {return coulomb;}
    inline bool FileIO::GetHistWriteStatus() {return writeHist;}
    inline int FileIO::GetPairType() {return mPairType;}
    inline int FileIO::GetEventsToMix() {return mEventsToMix;}
    inline int FileIO::GetEventFiles() {return mEventFiles;}
    inline float FileIO::GetTimeCut() {return mTimeCut;}
    inline float FileIO::GetMinKt() {return mKtMin;}
    inline float FileIO::GetMaxKt() {return mKtMax;}

#endif