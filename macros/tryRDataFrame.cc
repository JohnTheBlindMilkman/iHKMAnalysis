#include <iostream>
#include "TROOT.h"
#include "ROOT/RDataFrame.hxx"

using namespace ROOT;

void tryRDataFrame()
{
    //ROOT::EnableImplicitMT();
    const TString filepath = "/home/jedkol/lustre/hades/user/kjedrzej/iHKM/14p5GeV/0to10cent/0.root";
    const TString treeName = "treefin";

    auto pidCut = [](const RVec<int> &id)
    {
        return VecOps::Map(id,[](const RVec<int> &id){return (id == 2212);});
    };
    /*auto pt = [](RVec<float> px, RVec<float> py){return hypot(px,py);};
    auto rap = [](RVec<float> e, RVec<float> pz){return (0.5 * log((e + pz) / (e - pz)));};
    auto rapCut = [](RVec<double> rap){return (abs(rap) < 0.7);};
    auto peta = [](RVec<float> pt, RVec<float> pz){return -log(tan(atan2(pt,pz)/2.));};
    auto petaCut = [](RVec<double> peta){return (peta < 0.35);};*/

    RDataFrame df(treeName.Data(),filepath.Data(),{"px","py","pz","E","id"});

    // 1. select proton
    // 2. calc rapidity
    // 3. select rap < 0.7
    // 4. calc pt
    // 5. calc peta
    // 6. select peta < 0.35
    // 7. make TH1D
    //auto h1 = df.Filter(pidCut,{"id"}).Define("rap",rap,{"E","pz"}).Filter(rapCut,{"rap"}).Define("pt",pt,{"px","py"}).Define("peta",peta,{"pt","pz"}).Filter(petaCut,{"peta"}).Histo1D("pt");
    auto h1 = df.Filter(pidCut,"id").Histo1D("px");
    h1->Draw();
}