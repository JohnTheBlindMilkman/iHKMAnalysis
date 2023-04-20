#ifndef _MIXING_MANAGER_
    #define _MIXING_MANAGER_

    #include <TH1.h>
    #include "FemtoInteraction.h"
    #include "IOManager.h"

    class MixingManager : FemtoInteraction
    {
        private:
            std::unique_ptr<IOManager> mIOManager;
            TTreeReader read;
            TString mIniFile;
            bool mCoul, mPrim, mHist;
            int mPairType, mEvtsToMix, mPpid;
            float mTimeCut, mKtMin, mKtMax;

            bool TrackCut(const THGlobal::ParticleCoor &part);
            bool PairCut(const THGlobal::ParticleCoor &part1, const THGlobal::ParticleCoor &part2);
            std::vector<THGlobal::ParticleCoor> Filter(const std::vector<THGlobal::ParticleCoor> &partList, std::function<bool(THGlobal::ParticleCoor)> filtr);
            std::vector<std::array<double,2>> Mix(const std::vector<THGlobal::ParticleCoor> &partList);
            void AppendToHistogram(const std::vector<std::array<double,2>> &vec, TH1D *hNum, TH1D *hDen, TH1D *hWeight);

        public:
            MixingManager(/* args */);
            MixingManager(int argc, char *argv[]);
            ~MixingManager();
            void Init();
            void Start();
            void Finish();
    };

#endif