#ifndef _IO_MANAGER_
    #define _IO_MANAGER_

    #include <TTreeReader.h>
    #include <TTreeReaderArray.h>

    #include "FileIO.h"
    //#include "HistogramIO.h"

    class IOManager
    {
        private:
            FileIO mFileIO;
            //std::vector<std::unique_ptr<TObject>> mObjList;
            //TTreeReader reader;

        public:
            IOManager();
            IOManager(TString defIniFile, TString tmpFile);
            ~IOManager();
            void Init();
            void Close();
            static std::tuple<TString,int> Parse(int argc, char *argv[]);
            void AddHistogram();
            auto GetHistogram();
            auto GetIniParamter(THGlobal::INIParam Epar); //the whole parameter extraction could be a template
    };
    
#endif