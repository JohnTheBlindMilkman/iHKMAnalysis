#include "simu_lib.h"
using namespace std;

void simu_femto()
{
    const TString filePath = "./";
    const TString treeName = "treefin";
    const int numberFiles = 1;
    const int EventsToMix = 3;
    
    int mixEv; 
    float yi = 0, yj = 0;
    double quantumweight = 1., qscoulombweight = 1., coulombweight = 1.;
    TFile *ofile;
    TH1D *hCF_QS, *hCF_QSC, *hCF_C;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> distr(0,249);

    TH1D *hNumQS = new TH1D("hNumQS","",200,0,1.0);
    TH1D *hNumQSC = new TH1D("hNumQSC","",200,0,1.0);
    TH1D *hNumC = new TH1D("hNumC","",200,0,1.0);
    TH1D *hDen = new TH1D("hDen","",200,0,1.0); // czy mianowniki też niepowinny być oddzielne dla każdego licznika? - JJ
    
    TChain *chMain = new TChain(treeName.Data());
    TChain *chMainMix = new TChain(treeName.Data());

    for (int i = 0; i < numberFiles; i++)
    {
        chMain->Add(Form("%s%d.root",filePath.Data(),i));
        chMainMix->Add(Form("%s%d.root",filePath.Data(),i)); // może Clone() będzie bardziej efektywne? - JJ
    }
    
    TTreeReader ttReader(chMain);
    TTreeReaderValue<int> npart(ttReader,"npart");
    TTreeReaderArray<float> px(ttReader,"px");
    TTreeReaderArray<float> py(ttReader,"py");
    TTreeReaderArray<float> pz(ttReader,"pz");
    TTreeReaderArray<float> E(ttReader,"E");
    TTreeReaderArray<int> id(ttReader,"id");
    TTreeReaderArray<float> t(ttReader,"t");
    TTreeReaderArray<float> x(ttReader,"x");
    TTreeReaderArray<float> y(ttReader,"y");
    TTreeReaderArray<float> z(ttReader,"z");

    TTreeReader ttReaderMix(chMainMix);
    TTreeReaderValue<int> npartMix(ttReaderMix,"npart");
    TTreeReaderArray<float> pxMix(ttReaderMix,"px");
    TTreeReaderArray<float> pyMix(ttReaderMix,"py");
    TTreeReaderArray<float> pzMix(ttReaderMix,"pz");
    TTreeReaderArray<float> EMix(ttReaderMix,"E");
    TTreeReaderArray<int> idMix(ttReaderMix,"id");
    TTreeReaderArray<float> tMix(ttReaderMix,"t");
    TTreeReaderArray<float> xMix(ttReaderMix,"x");
    TTreeReaderArray<float> yMix(ttReaderMix,"y");
    TTreeReaderArray<float> zMix(ttReaderMix,"z");

    while (ttReader.Next())
    {
        cout << "tt " << ttReader.GetCurrentEntry() << "/250" << endl;//, (500 w jednym pliku tym na którym testowałam)"<<endl;
        for (int i = 0; i < *npart; i++)
        {
            yi = 0.5 * TMath::Log((E[i] + pz[i])/(E[i] - pz[i]));  
            if (id[i] == part1pid && fabs(yi) < ABSRAP)
            {
       	    	for (int j = 0; j < *npart; j++) // czy tu nie ma być npartMix? - JJ
            	{
                    yj = 0.5 * TMath::Log((E[j] + pz[j])/(E[j] - pz[j]));
                    if (i != j && id[j] == part2pid && abs(yj) < ABSRAP)
                    {
                        PairKinematics(i, j, px, py, pz, E, t, x, y, z); // czy przepisywanie na zmienne TTreeReaderArray<float> było konieczne? nadal podajesz indeksy, może lepiej przekazać zwykłe floaty? - JJ
                        
                        //quantum statistics only
                        if (twospin == 0) 
                            quantumweight = 1.0+TMath::Cos(-mKO*mRO - mKS*mRS - mKL*mRL + mDE*mDT);
                        
                        else if (twospin == 1) 
                        {
                            if (pcount == 3)
                                quantumweight = 1.0+TMath::Cos(-mKO*mRO - mKS*mRS - mKL*mRL + mDE*mDT);
                            else
                                quantumweight = 1.0-TMath::Cos(-mKO*mRO - mKS*mRS - mKL*mRL + mDE*mDT);
                            
                            pcount++;
                            if (pcount == 4) 
                                pcount=0;
                        }
            
                        hNumQS->Fill(fabs(2.0*mKStarSigned/0.197327),quantumweight); // może bezpieczniej zrobić zmienną z tej liczby? - JJ
                        
                        //quantum statistics and coulomb
                        InitializeGamow();
                        qscoulombweight = GetQuantumCoulomb();
                        hNumQSC->Fill(fabs(2.0*mKStarSigned/0.197327),qscoulombweight);		
                        
                        //coulomb only 	         
                        coulombweight = GetCoulomb();
                        hNumC->Fill(fabs(2.0*mKStarSigned/0.197327),coulombweight);
                        
                        /*
                        //coulomb strong quantum statistics
                        double coulombstrongqsweight = 
                        */
                    }
                }
                
                mixEv = 0;
                while (mixEv < EventsToMix)
                {
                    ttReaderMix.SetEntry(distr(gen));
                    if(ttReader.GetCurrentEntry() != ttReaderMix.GetCurrentEntry())
                    {
                        mixEv++;
                        for (int j = 0; j < *npartMix; j++)
                            if(idMix[j] == part2pid)
                            {
                                PairKinematics2(i, j, px, py, pz, E, pxMix, pyMix, pzMix, EMix);
                                hDen->Fill(fabs(2.0*mKStarSigned/0.197327),1.0);				   
                            }
                    }
                }
            }
        }
    }
    
    ofile = TFile::Open("results.root","recreate");

    hDen->Write();
    hNumQS->Write();
    hNumQSC->Write();
    hNumC->Write();
    
    hCF_QS = getCF(hNumQS, hDen, "cf_qs",0.3,0.5); 
    hCF_QSC = getCF(hNumQSC, hDen, "cf_qsc",0.3,0.5); 
    hCF_C = getCF(hNumC, hDen, "cf_c",0.3,0.5); 
    
    hCF_QS->Write();
    hCF_QSC->Write();
    hCF_C->Write();

    ofile->Close();
}
