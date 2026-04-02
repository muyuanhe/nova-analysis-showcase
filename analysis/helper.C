// Make two simple spectrum plots

#include "CAFAna/Core/Binning.h"
#include "CAFAna/Cuts/Cuts.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Vars/Vars.h"
#include "3FlavorAna/Cuts/NumuCuts.h" // kNumuCC cut
#include "3FlavorAna/Cuts/NueCuts2024.h" // kNue2024FD cut https://github.com/novaexperiment/novasoft/blob/main/3FlavorAna/Cuts/NueCuts2024.h
// reconstructs direction
#include "CosRej/NueCosRej.h" //https://github.com/novaexperiment/novasoft/blob/a86315774e7fc0f3e6e1a39970ba1ec723e0754c/CosRej/NueCosRej.h
#include "TVector3.h"


#include <cmath>
#include "StandardRecord/Proxy/SRProxy.h"
#include "OscLib/OscCalcPMNSOpt.h"

#include "TCanvas.h"
#include "TH2.h"
#include "TLegend.h"
#include "TFile.h"

// Containment cut from 3flavor
#include "3FlavorAna/Cuts/NumuCuts2024.h"
#include "3FlavorAna/Cuts/NueCuts2024.h"

#include "3FlavorAna/Vars/NueEnergy2024.h"
#include "3FlavorAna/Vars/NumuVars.h"

#include "Utilities/rootlogon.C"


TFile *Oscillogram_NSImutomu = TFile::Open("/exp/nova/app/users/muyuanh/OscProb_11_01_2024/OscProb/tutorial/histogram_dir/Oscillograms_mutomu.root");
TFile *Oscillogram_NSImutoe = TFile::Open("/exp/nova/app/users/muyuanh/OscProb_11_01_2024/OscProb/tutorial/histogram_dir/Oscillograms_mutoe.root");
TFile *Oscillogram_NSIetomu = TFile::Open("/exp/nova/app/users/muyuanh/OscProb_11_01_2024/OscProb/tutorial/histogram_dir/Oscillograms_etomu.root");
TFile *Oscillogram_NSIetoe = TFile::Open("/exp/nova/app/users/muyuanh/OscProb_11_01_2024/OscProb/tutorial/histogram_dir/Oscillograms_etoe.root");

using namespace ana;

// Defined later
extern const NuTruthCut kIsTrueNumuCC_NT;
extern const Cut kIsTrueNumuCC;
extern const NuTruthCut kIsTrueNueCC_NT;
extern const Cut kIsTrueNueCC;
// extern const NuTruthCut kNoCut_NT;
extern const Cut kIsTrueNC;
// extern const NuTruthCut kNoCut_NT;
// extern const Cut kNoCut;
// extern const Var kNHit;
extern const Var kNtracks;
extern const MultiVar kProtonProngL;
extern const Var kVertexX;
extern const Var kVertexY;
extern const Var kVertexZ;
extern const Var kVertexXCosmic;
extern const Var kVertexYCosmic;
extern const Var kVertexZCosmic;
extern const MultiVar kRecoShwE;
extern const MultiVar kCosmicRecoShwE;
extern const MultiVar kMultiNHit;
extern const MultiVar kCosmicMultiNHit;
extern const Var knumofpng;
extern const Var kcosmicnumofpng;
extern const Var knueE;
extern const Var kTotalCalE;
extern const Var kCosmicTotalCalE;
extern const Var kCalEPerHit;
extern const Cut numuContainCut;
extern const Cut nueContainCut;
extern const MultiVar kCalEPerHit_total;
extern const MultiVar kCosmicCalEPerHit_total;
extern const Var kCosmicCalEPerHit;
extern const Cut kIsValid;
// NNBar variables 
extern const Var totalNHit;
extern const Var MinStartX;
extern const Var MaxStartX;
extern const Var MinStartY;
extern const Var MaxStartY;
extern const Var MinStartZ;
extern const Var MaxStartZ;
extern const Var Asymmetry;
// NNBar variables not defined yet

// Declare oscillation weights for numu to numu oscillation with various NSI parameters
extern const Weight OscProb_mutomu_epsee;
extern const Weight OscProb_mutomu_epsemu;
extern const Weight OscProb_mutomu_epsetau;
extern const Weight OscProb_mutomu_epsmumu;
extern const Weight OscProb_mutomu_epsmutau;
extern const Weight OscProb_mutomu_epstauetau;

// Declare oscillation weights for numu to nue oscillation with various NSI parameters
extern const Weight OscProb_mutoe_epsee;
extern const Weight OscProb_mutoe_epsemu;
extern const Weight OscProb_mutoe_epsetau;
extern const Weight OscProb_mutoe_epsmumu;
extern const Weight OscProb_mutoe_epsmutau;
extern const Weight OscProb_mutoe_epstauetau;

// Declare oscillation weights for nue to numu oscillation with various NSI parameters
extern const Weight OscProb_etomu_epsee;
extern const Weight OscProb_etomu_epsemu;
extern const Weight OscProb_etomu_epsetau;
extern const Weight OscProb_etomu_epsmumu;
extern const Weight OscProb_etomu_epsmutau;
extern const Weight OscProb_etomu_epstauetau;

// Declare oscillation weights for nue to nue oscillation with various NSI parameters
extern const Weight OscProb_etoe_epsee;
extern const Weight OscProb_etoe_epsemu;
extern const Weight OscProb_etoe_epsetau;
extern const Weight OscProb_etoe_epsmumu;
extern const Weight OscProb_etoe_epsmutau;
extern const Weight OscProb_etoe_epstauetau;

// Declare oscillation spectrums
extern const Var kCosT;
extern const Var kLoE;

extern const Weight OscProb_baseline_mutomu;
extern const Weight OscProb_baseline_mutoe;
extern const Weight OscProb_baseline_etomu;
extern const Weight OscProb_baseline_etoe;

extern const Cut numurecoCut;
extern const Cut nuerecoCut;

extern const Var calcRecMomentum_single_particle;
extern const Var Resolution_TransMomFraction;










const Binning bins = Binning::Simple(100, 0, 1000);
const Binning bins2 = Binning::Simple(10, 0, 10);
const Binning binsxy = Binning::Simple(1600, -800, 800);
const Binning binsz = Binning::Simple(6000, 0, 6000);
const Binning binsE = Binning::Simple(15, 0, 10);
const Binning binsE_nue = Binning::Simple(45, 0, 10);
const Binning binsresE = Binning::Simple(125, -5, 5);
const Binning binsres_nueE = Binning::Simple(50, -2, 2);
const Binning binsnhits = Binning::Simple(40, 0, 200);
const Binning calEPerHitBins = Binning::Simple(100, 0, 0.1);
const Binning bin_cos_theta = Binning::Simple(100, -1, 1);
const Binning bin_theta = Binning::Simple(180, 0, 180);
const Binning bin_LoE = Binning::Simple(50, 0, 10000);
const Binning binstheta = Binning::Simple(10, -1, 1);
const Binning binsCosTheta = Binning::Simple(7, -1, 1);
const Binning binsResolution = Binning::Simple(20, -2, 2);
const Binning binsresolutionvertex = Binning::Simple(24, -1200, 1200);
const Binning binsvtxdist = Binning::Simple(20, 0, 200);

const NuTruthCut kIsTrueNumuCC_NT([](const caf::SRNeutrinoProxy* truth) //true numu CC event
        {
            return (truth->iscc && truth->pdg == 14);
        }); 

  
const Cut kIsTrueNumuCC      = CutFromNuTruthCut(kIsTrueNumuCC_NT);

const NuTruthCut kIsTrueNueCC_NT([](const caf::SRNeutrinoProxy* truth) //true nueCC
        {
        return (truth->iscc && truth->pdg == 12);
        }
        );
const Cut kIsTrueNueCC       = CutFromNuTruthCut(kIsTrueNueCC_NT);  

const NuTruthCut kIsTrueNC_NT([](const caf::SRNeutrinoProxy* truth) {
  return !truth->iscc; // Neutral Current event: iscc is false
});
const Cut kIsTrueNC = CutFromNuTruthCut(kIsTrueNC_NT);


// const NuTruthCut kNoCut_NT([](const caf::SRNeutrinoProxy* truth) //true no cut
//         {
//         return true;
//         });

// const Cut kNoCut      = CutFromNuTruthCut(kNoCut_NT); 

// const Var kNHit([](const caf::SRProxy* sr)
//         {
//         //don't need a check since we aren't looking at any specific part of a vector/array
//         return sr->slc.nhit;
//         });

const Var kNtracks([](const caf::SRProxy* sr)
        {
        return sr->trk.kalman.ntracks;
        });

const MultiVar kProtonProngL([] (const caf::SRProxy* sr)
        {
        std::vector<double> TrkLen;
        for(int i=0; i<sr->vtx.elastic.fuzzyk.png.size(); i++)
          {
            if(sr->vtx.elastic.fuzzyk.png[i].truth.pdg==2212)
              {
                //TrkLen.push_back(sr->vtx.elastic.fuzzyk.png[i].len);
              }
            TrkLen.push_back(sr->vtx.elastic.fuzzyk.png[i].len); //for atmospheric
          }
        return TrkLen;
        });

MultiVarHistAxis ProtonProngLaxis("Proton Prong Length (cm)",bins,kProtonProngL);

const Cut kIsValid([](const caf::SRProxy* sr) // to cut off failed reconstrucion
        {
        return sr->vtx.elastic.IsValid;
        });

const Var kVertexX([](const caf::SRProxy* sr)
        {
        if (sr->vtx.elastic.IsValid){
        return float(sr->vtx.elastic.vtx.x); //for atmospheric
        }
        });

const Var kVertexY([](const caf::SRProxy* sr)
        {
        if (sr->vtx.elastic.IsValid){
        return float(sr->vtx.elastic.vtx.y); //for atmospheric
        }
        });
    
const Var kVertexZ([](const caf::SRProxy* sr)
        {
        if (sr->vtx.elastic.IsValid){
        return float(sr->vtx.elastic.vtx.z);  //for atmospheric
        }
        });


const Var kVertexXCosmic([](const caf::SRProxy* sr)
        {
        if (sr->vtx.elastic.IsValid){
        return float(sr->vtx.elastic.vtx.x);    //for cosmic
        }
        });

const Var kVertexYCosmic([](const caf::SRProxy* sr)
        {
        if (sr->vtx.elastic.IsValid){
        
        return float(sr->vtx.elastic.vtx.y);     //for cosmic
        }
        });
    
const Var kVertexZCosmic([](const caf::SRProxy* sr)
        {
        if (sr->vtx.elastic.IsValid){
          return float(sr->vtx.elastic.vtx.z);     //for cosmic
        }
        });



const MultiVar kRecoShwE([](const caf::SRProxy* sr)
        {
        std::vector<double> ShwEnergy;
        for(int i=0; i<sr->vtx.elastic.fuzzyk.png.size(); i++)
          {
            ShwEnergy.push_back(sr->vtx.elastic.fuzzyk.png[i].shwlid.lidE.shwE);
          }
        return ShwEnergy;
        });
MultiVarHistAxis RecoShwEaxis("Reconstructed shower energy",binsE,kRecoShwE);

const MultiVar kCosmicRecoShwE([](const caf::SRProxy* sr)
        {
        std::vector<double> ShwEnergy;
        for(int i=0; i<sr->vtx.elastic.fuzzyk.png.size(); i++)
          {
            ShwEnergy.push_back(sr->vtx.elastic.fuzzyk.png[i].shwlid.shwE);
          }
        return ShwEnergy;
        });
MultiVarHistAxis RecoCosmicShwEaxis("Reconstructed cosmic shower energy",binsE,kCosmicRecoShwE);


const MultiVar kMultiNHit([](const caf::SRProxy* sr) // number of hits for shower
        {
        std::vector<double> vectorNHits;
        for(int i=0; i<sr->vtx.elastic.fuzzyk.png.size(); i++) //vtx.elastic.fuzzyk.png.nhit
          {
            vectorNHits.push_back(sr->vtx.elastic.fuzzyk.png[i].nhit);
          }
        return vectorNHits;
        });
MultiVarHistAxis NHitaxis("Number of hits",binsnhits,kMultiNHit);

const MultiVar kCosmicMultiNHit([](const caf::SRProxy* sr) // number of hits for cosmic shower
        {
        std::vector<double> vectorNHits;
        for(int i=0; i<sr->vtx.elastic.fuzzyk.png.size(); i++)
          {
            vectorNHits.push_back(sr->vtx.elastic.fuzzyk.png[i].nhit);
          }
        return vectorNHits;
        });
MultiVarHistAxis CosmicNHitaxis("Number of hits",binsnhits,kCosmicMultiNHit);

const Var knumofpng([](const caf::SRProxy* sr)
        {
        return sr->vtx.elastic.fuzzyk.npng; //atm
        // return sr->vtx.elastic.fuzzyk.npng; //cosmic
        });

const Var kcosmicnumofpng([](const caf::SRProxy* sr)
        {
        return sr->vtx.elastic.fuzzyk.npng; //cosmic
        });

const Var knueE([](const caf::SRProxy* sr) //reconstructed energy for neutrino
        {
        return sr->energy.nue.lid.E;
        });

const Var kTotalCalE([](const caf::SRProxy* sr) { //Calorimetric total E
    double totalCalE = 0;
    for (int i = 0; i < sr->vtx.elastic.fuzzyk.png.size(); i++) {
        totalCalE += sr->vtx.elastic.fuzzyk.png[i].calE; // Summing up calorimetric energy of each prong
    }
    return totalCalE;
});

const Var kCosmicTotalCalE([](const caf::SRProxy* sr) { //Calorimetric total E
    double totalCalE = 0;
    for (int i = 0; i < sr->vtx.elastic.fuzzyk.png.size(); i++) {
        totalCalE += sr->vtx.elastic.fuzzyk.png[i].calE; // Summing up calorimetric energy of each prong
    }
    return totalCalE;
});

const Var kCalEPerHit([](const caf::SRProxy* sr) { //Calorimetric total E per hit
    double totalCalE = 0;
    int totalNHits = 0;
    for (int i = 0; i < sr->vtx.elastic.fuzzyk.png.size(); i++) {
        totalCalE += sr->vtx.elastic.fuzzyk.png[i].calE;   // Summing up calorimetric energy
        totalNHits += sr->vtx.elastic.fuzzyk.png[i].nhit;  // Summing up number of hits
    }
    if (totalNHits > 0) {  // Avoid division by zero
        return totalCalE / totalNHits;
    } else {
        return -10.0;  // Or some other suitable value when there are no hits
    }
});

const Var kCosmicCalEPerHit([](const caf::SRProxy* sr) { //Calorimetric total E per hit
    double totalCalE = 0;
    int totalNHits = 0;
    for (int i = 0; i < sr->vtx.elastic.fuzzyk.png.size(); i++) {
        totalCalE += sr->vtx.elastic.fuzzyk.png[i].calE;   // Summing up calorimetric energy
        totalNHits += sr->vtx.elastic.fuzzyk.png[i].nhit;  // Summing up number of hits
    }
    if (totalNHits > 0) {  // Avoid division by zero
        return totalCalE / totalNHits;
    } else {
        return -10.0;  // Or some other suitable value when there are no hits
    }
});

const MultiVar kCalEPerHit_total([](const caf::SRProxy* sr)  //Calorimetric total E per hit, multivar
{
    std::vector<double> vectorCalEperhit;
    for (int i = 0; i < sr->vtx.elastic.fuzzyk.png.size(); i++) 
    {
        vectorCalEperhit.push_back(sr->vtx.elastic.fuzzyk.png[i].calE / sr->vtx.elastic.fuzzyk.png[i].nhit);  // Summing up number of hits
    }
    return vectorCalEperhit;
});
MultiVarHistAxis CalEPerHit_totalaxis("Total Avg. Calorimetric Energy per Hit (GeV/hit)",calEPerHitBins,kCalEPerHit_total);

const MultiVar kCosmicCalEPerHit_total([](const caf::SRProxy* sr)  //Calorimetric total E per hit, multivar
{
    std::vector<double> vectorCosmicCalEperhit;
    for (int i = 0; i < sr->vtx.elastic.fuzzyk.png.size(); i++) 
    {
        vectorCosmicCalEperhit.push_back(sr->vtx.elastic.fuzzyk.png[i].calE / sr->vtx.elastic.fuzzyk.png[i].nhit);  // Summing up number of hits
    }
    return vectorCosmicCalEperhit;
});
MultiVarHistAxis CosmicCalEPerHit_totalaxis("Total Avg. Calorimetric Energy per Hit (GeV/hit)",calEPerHitBins,kCosmicCalEPerHit_total);



const Cut numuContainCut = kNumuContainFD2024;
const Cut nueContainCut = kNue2024ProngContain;

const Cut numuContainCut_modified([](const caf::SRProxy* sr)
        {
        // containment cut
        for( unsigned int i = 0; i < sr->vtx.elastic.fuzzyk.nshwlid; ++i ) {
            TVector3 start = sr->vtx.elastic.fuzzyk.png[i].shwlid.start;
            TVector3 stop  = sr->vtx.elastic.fuzzyk.png[i].shwlid.stop;
            if( std::min( start.X(), stop.X() ) < -784.52 ) return false;
            if( std::max( start.X(), stop.X() ) >  788.39 ) return false;
            if( std::min( start.Y(), stop.Y() ) < -776.78 ) return false;
            if( std::max( start.Y(), stop.Y() ) >  672.29 ) return false;
            if( std::min( start.Z(), stop.Z() ) <   19.8 ) return false;
            if( std::max( start.Z(), stop.Z() ) > 5880.6 ) return false;
            }
        });

const Var totalNHit([](const caf::SRProxy* sr)
        {
        if (sr->vtx.elastic.fuzzyk.npng == 0) return -99999;
        if (sr->vtx.elastic.fuzzyk.nshwlid == 0) return -99999;
        int totalNHitsx = 0;
        int totalNHitsy = 0;
        int totalNHits = 0;
        for (int i = 0; i < sr->vtx.elastic.fuzzyk.png.size(); i++) {
                totalNHitsx += sr->vtx.elastic.fuzzyk.png[i].nhitx;  // Summing up number of hits
                totalNHitsy += sr->vtx.elastic.fuzzyk.png[i].nhity;  // Summing up number of hits
        }
        totalNHits = totalNHitsx + totalNHitsy;
        return (totalNHits);
        });

const Var MinStartX ([](const caf::SRProxy* sr)
        {
        double minStartX = 800;
        for( unsigned int i = 0; i < sr->vtx.elastic.fuzzyk.nshwlid; ++i ) {
            TVector3 start = sr->vtx.elastic.fuzzyk.png[i].shwlid.start;
            TVector3 stop  = sr->vtx.elastic.fuzzyk.png[i].shwlid.stop;
            if (std::min( start.X(), stop.X() ) < minStartX) minStartX = std::min( start.X(), stop.X());
        }      
        return minStartX;
        });
const Var MaxStartX ([](const caf::SRProxy* sr)
        {
        double maxStartX = -800;
        for( unsigned int i = 0; i < sr->vtx.elastic.fuzzyk.nshwlid; ++i ) {
            TVector3 start = sr->vtx.elastic.fuzzyk.png[i].shwlid.start;
            TVector3 stop  = sr->vtx.elastic.fuzzyk.png[i].shwlid.stop;
            if( std::max( start.X(), stop.X() ) > maxStartX) maxStartX = std::max( start.X(), stop.X() );
            }
        return maxStartX;
        });

const Var MinStartY ([](const caf::SRProxy* sr)
        {
        double minStartY = 800;
        for( unsigned int i = 0; i < sr->vtx.elastic.fuzzyk.nshwlid; ++i ) {
            TVector3 start = sr->vtx.elastic.fuzzyk.png[i].shwlid.start;
            TVector3 stop  = sr->vtx.elastic.fuzzyk.png[i].shwlid.stop;
            if( std::min( start.Y(), stop.Y() ) < minStartY) minStartY = std::min( start.Y(), stop.Y() );
            }
        return minStartY;
        });

const Var MaxStartY ([](const caf::SRProxy* sr)
        {
        double maxStartY = -800;
        for( unsigned int i = 0; i < sr->vtx.elastic.fuzzyk.nshwlid; ++i ) {
            TVector3 start = sr->vtx.elastic.fuzzyk.png[i].shwlid.start;
            TVector3 stop  = sr->vtx.elastic.fuzzyk.png[i].shwlid.stop;
            if( std::max( start.Y(), stop.Y() ) > maxStartY) maxStartY = std::max( start.Y(), stop.Y() );
            }
        return maxStartY;
        });

const Var MinStartZ ([](const caf::SRProxy* sr)
        {
        double minStartZ = 800;
        for( unsigned int i = 0; i < sr->vtx.elastic.fuzzyk.nshwlid; ++i ) {
            TVector3 start = sr->vtx.elastic.fuzzyk.png[i].shwlid.start;
            TVector3 stop  = sr->vtx.elastic.fuzzyk.png[i].shwlid.stop;
            if( std::min( start.Z(), stop.Z() ) < minStartZ) minStartZ = std::min( start.Z(), stop.Z() );
            }
        return minStartZ;
        });

const Var MaxStartZ ([](const caf::SRProxy* sr)
        {
        double maxStartZ = -800;
        for( unsigned int i = 0; i < sr->vtx.elastic.fuzzyk.nshwlid; ++i ) {
            TVector3 start = sr->vtx.elastic.fuzzyk.png[i].shwlid.start;
            TVector3 stop  = sr->vtx.elastic.fuzzyk.png[i].shwlid.stop;
            if( std::max( start.Z(), stop.Z() ) > maxStartZ) maxStartZ = std::max( start.Z(), stop.Z() );
            }
        return maxStartZ;
        });

const Var Asymmetry([](const caf::SRProxy* sr)
        {
        if (sr->vtx.elastic.fuzzyk.npng == 0) return -99999;
        if (sr->vtx.elastic.fuzzyk.nshwlid == 0) return -99999;
        int totalNHitsx = 0;
        int totalNHitsy = 0;
        for (int i = 0; i < sr->vtx.elastic.fuzzyk.png.size(); i++) {
                totalNHitsx += sr->vtx.elastic.fuzzyk.png[i].nhitx;  // Summing up number of hits
                totalNHitsy += sr->vtx.elastic.fuzzyk.png[i].nhity;  // Summing up number of hits
        }
        return abs(totalNHitsx - totalNHitsy) / (totalNHitsx + totalNHitsy);
        });
// Function to find the bin and retrieve the weight from the oscillogram

/////////////////////////////////////////// Oscillogram ///////////////////////////////////////////
float findOscBin(const caf::SRProxy* sr, TH2F* oscillogram)
{
    float px = sr->mc.nu[0].p.px;
    float py = sr->mc.nu[0].p.py;
    float pz = sr->mc.nu[0].p.pz;
    float cosT = py / sqrt(px * px + py * py + pz * pz);
    float E = sr->mc.nu[0].E;
    constexpr float DIAMETER = 6386.0f; // Earth's diameter in km
    constexpr float DetDiameter = 6371.0f; // Detector diameter in km
    float L = -1 * DetDiameter * cosT + sqrt(DIAMETER * DIAMETER - DetDiameter * DetDiameter * (1 - cosT * cosT));

    int bin = oscillogram->FindBin(L / E, cosT);
    return oscillogram->GetBinContent(bin); // Retrieve bin content as the weight
}


// Function to create a weight based on an oscillation type and NSI parameter
const Weight CreateOscProbWeight(const TString& oscType, const TString& nsiParam) {
    return [oscType, nsiParam](const caf::SRProxy* sr) -> float {
        // Select the file based on the oscillation type
        TFile* file;
        if (oscType == "mutomu") {
            file = Oscillogram_NSImutomu;
        } else if (oscType == "mutoe") {
            file = Oscillogram_NSImutoe;
        } else if (oscType == "etomu") {
            file = Oscillogram_NSIetomu;
        } else if (oscType == "etoe") {
            file = Oscillogram_NSIetoe;
        } else {
            std::cerr << "Error: Unknown oscillation type " << oscType << std::endl;
            return 1.0;
        }
        TString histName = Form("Oscillogram_%s", nsiParam.Data());
        const auto hist = (TH2F*)file->Get(histName);
        
        if (!hist) {
            std::cerr << "Error: Histogram " << histName << " not found!" << std::endl;
            return 1.0; // Default weight if histogram is not found
        }
        return findOscBin(sr, hist);
    };
}

// Define weights for different oscillation types and NSI parameters
const std::vector<TString> nsiParams = {
    "eps_ee", "eps_emu", "eps_etau", "eps_mumu", "eps_mutau", "eps_tautau"
};

// Define weights for different oscillation scenarios
const Weight OscProb_mutomu_epsee = CreateOscProbWeight("mutomu", "eps_ee");
const Weight OscProb_mutomu_epsemu = CreateOscProbWeight("mutomu", "eps_emu");
const Weight OscProb_mutomu_epsetau = CreateOscProbWeight("mutomu", "eps_etau");
const Weight OscProb_mutomu_epsmumu = CreateOscProbWeight("mutomu", "eps_mumu");
const Weight OscProb_mutomu_epsmutau = CreateOscProbWeight("mutomu", "eps_mutau");
const Weight OscProb_mutomu_epstautau = CreateOscProbWeight("mutomu", "eps_tautau");

const Weight OscProb_mutoe_epsee = CreateOscProbWeight("mutoe", "eps_ee");
const Weight OscProb_mutoe_epsemu = CreateOscProbWeight("mutoe", "eps_emu");
const Weight OscProb_mutoe_epsetau = CreateOscProbWeight("mutoe", "eps_etau");
const Weight OscProb_mutoe_epsmumu = CreateOscProbWeight("mutoe", "eps_mumu");
const Weight OscProb_mutoe_epsmutau = CreateOscProbWeight("mutoe", "eps_mutau");
const Weight OscProb_mutoe_epstautau = CreateOscProbWeight("mutoe", "eps_tautau");

const Weight OscProb_etomu_epsee = CreateOscProbWeight("etomu", "eps_ee");
const Weight OscProb_etomu_epsemu = CreateOscProbWeight("etomu", "eps_emu");
const Weight OscProb_etomu_epsetau = CreateOscProbWeight("etomu", "eps_etau");
const Weight OscProb_etomu_epsmumu = CreateOscProbWeight("etomu", "eps_mumu");
const Weight OscProb_etomu_epsmutau = CreateOscProbWeight("etomu", "eps_mutau");
const Weight OscProb_etomu_epstautau = CreateOscProbWeight("etomu", "eps_tautau");

const Weight OscProb_etoe_epsee = CreateOscProbWeight("etoe", "eps_ee");
const Weight OscProb_etoe_epsemu = CreateOscProbWeight("etoe", "eps_emu");
const Weight OscProb_etoe_epsetau = CreateOscProbWeight("etoe", "eps_etau");
const Weight OscProb_etoe_epsmumu = CreateOscProbWeight("etoe", "eps_mumu");
const Weight OscProb_etoe_epsmutau = CreateOscProbWeight("etoe", "eps_mutau");
const Weight OscProb_etoe_epstautau = CreateOscProbWeight("etoe", "eps_tautau");
/////////////////////////////////////////// cosT / L/E ///////////////////////////////////////////
const Var kCosT([](const caf::SRProxy* sr)
        {
        float px = sr->mc.nu[0].p.px;
        float py = sr->mc.nu[0].p.py;
        float pz = sr->mc.nu[0].p.pz;
        return py / sqrt(px * px + py * py + pz * pz);
        });

const Var kTheta([](const caf::SRProxy* sr)
{
    float px = sr->mc.nu[0].p.px;
    float py = sr->mc.nu[0].p.py;
    float pz = sr->mc.nu[0].p.pz;
    float magnitude = sqrt(px * px + py * py + pz * pz);
    return std::acos(py / magnitude) * 180 / M_PI; // Result in degrees
});

const Var kCosTheta([](const caf::SRProxy* sr)
        {
        float px = sr->mc.nu[0].p.px;
        float py = sr->mc.nu[0].p.py;
        float pz = sr->mc.nu[0].p.pz;
        return py / sqrt(px * px + py * py + pz * pz);
        });

const Var kLoE([](const caf::SRProxy* sr)
        {
        float px = sr->mc.nu[0].p.px;
        float py = sr->mc.nu[0].p.py;
        float pz = sr->mc.nu[0].p.pz;
        float cosT = py / sqrt(px * px + py * py + pz * pz);
        float E = sr->mc.nu[0].E;
        constexpr float DIAMETER = 6386.0f; // Earth's diameter in km
        constexpr float DetDiameter = 6371.0f; // Detector diameter in km
        float L = -1 * DetDiameter * cosT + sqrt(DIAMETER * DIAMETER - DetDiameter * DetDiameter * (1 - cosT * cosT));
        return L / E;
        });

const Var kE([](const caf::SRProxy* sr)
        {
        return sr->mc.nu[0].E;
        });

HistAxis CosTaxis("cos(theta)",bin_cos_theta,kCosT);
HistAxis LoEaxis("L/E",bin_LoE,kLoE);
HistAxis LoEvsCosT("kLoE",bin_LoE,kLoE, "cos(theta)",bin_cos_theta,kCosT);
HistAxis LoEvsTheta("kLoE",bin_LoE,kLoE, "Zenith Angle (degrees)",bin_theta,kTheta);
HistAxis EvsTheta("Energy (GeV)",binsE,kE, "Zenith Angle (degrees)",bin_theta,kTheta);


// /////////////////////////////////////////// For event_rate.C ///////////////////////////////////////////
// float findOscBin(const caf::SRProxy* sr, TH2F* oscillogram)
// {
//     float px = sr->mc.nu[0].p.px;
//     float py = sr->mc.nu[0].p.py;
//     float pz = sr->mc.nu[0].p.pz;
//     float cosT = py / sqrt(px * px + py * py + pz * pz);
//     float E = sr->mc.nu[0].E;
//     constexpr float DIAMETER = 6386.0f; // Earth's diameter in km
//     constexpr float DetDiameter = 6371.0f; // Detector diameter in km
//     float L = -1 * DetDiameter * cosT + sqrt(DIAMETER * DIAMETER - DetDiameter * DetDiameter * (1 - cosT * cosT));

//     int bin = oscillogram->FindBin(L / E, cosT);
//     return oscillogram->GetBinContent(bin); // Retrieve bin content as the weight
// }

// const Weight CreateOscProbWeight(const TString& oscType, const TString& nsiParam, const TString& initialFlavor) {
//     return [oscType, nsiParam, initialFlavor](const caf::SRProxy* sr) -> float {
//         // Check compatibility based on initial flavor
//         bool isNue = (initialFlavor == "nue");
//         bool isNumu = (initialFlavor == "numu");

//         if ((oscType == "mutomu" || oscType == "mutoe") && !isNumu) return 1.0;
//         if ((oscType == "etomu" || oscType == "etoe") && !isNue) return 1.0;

//         // Load the appropriate histogram based on oscillation type
//         TString histName = Form("Oscillogram_%s", nsiParam.Data());
//         TH2F* hist = nullptr;

//         if (oscType == "mutomu") hist = (TH2F*)Oscillogram_NSImutomu->Get(histName);
//         else if (oscType == "mutoe") hist = (TH2F*)Oscillogram_NSImutoe->Get(histName);
//         else if (oscType == "etomu") hist = (TH2F*)Oscillogram_NSIetomu->Get(histName);
//         else if (oscType == "etoe") hist = (TH2F*)Oscillogram_NSIetoe->Get(histName);

//         if (!hist) {
//             std::cerr << "Error: Histogram " << histName << " not found!" << std::endl;
//             return 1.0; // Default weight if histogram is not found
//         }

//         return findOscBin(sr, hist);
//     };
// }

// // const Binning bin_cos_theta = Binning::Simple(100, -1, 1);
// // const Binning bin_LoE = Binning::Simple(100, 0, 10000);

// const Var kCosT([](const caf::SRProxy* sr) {
//     float px = sr->mc.nu[0].p.px;
//     float py = sr->mc.nu[0].p.py;
//     float pz = sr->mc.nu[0].p.pz;
//     return py / sqrt(px * px + py * py + pz * pz);
// });

// const Var kLoE([](const caf::SRProxy* sr) {
//     float px = sr->mc.nu[0].p.px;
//     float py = sr->mc.nu[0].p.py;
//     float pz = sr->mc.nu[0].p.pz;
//     float cosT = py / sqrt(px * px + py * py + pz * pz);
//     float E = sr->mc.nu[0].E;
//     constexpr float DIAMETER = 6386.0f; // Earth's diameter in km
//     constexpr float DetDiameter = 6371.0f; // Detector diameter in km
//     float L = -1 * DetDiameter * cosT + sqrt(DIAMETER * DIAMETER - DetDiameter * DetDiameter * (1 - cosT * cosT));
//     return L / E;
// });

// // Define the 2D histogram axis for L/E vs cos(theta)
// HistAxis LoEvsCosT("L/E", bin_LoE, kLoE, "cos(theta)", bin_cos_theta, kCosT);

// Function to create a weight for the baseline (standard) oscillation

/// Create a weight based on the standrad oscillation  ///
const Weight CreateBaselineOscWeight(const TString& oscType) {
    return [oscType](const caf::SRProxy* sr) -> float {
        TString histName = "Baseline_Oscillogram";
        TFile* file;

        // Select the file based on the oscillation type
        if (oscType == "mutomu") {
            file = Oscillogram_NSImutomu;
        } else if (oscType == "mutoe") {
            file = Oscillogram_NSImutoe;
        } else if (oscType == "etomu") {
            file = Oscillogram_NSIetomu;
        } else if (oscType == "etoe") {
            file = Oscillogram_NSIetoe;
        } else {
            std::cerr << "Error: Unknown oscillation type " << oscType << std::endl;
            return 1.0;
        }

        // Load the baseline histogram
        const auto hist = (TH2F*)file->Get(histName);
        if (!hist) {
            std::cerr << "Error: Histogram " << histName << " not found!" << std::endl;
            return 1.0; // Default weight if histogram is not found
        }

        return findOscBin(sr, hist);
    };
}
// Define weights for the baseline oscillation scenarios without NSI
const Weight OscProb_baseline_mutomu = CreateBaselineOscWeight("mutomu");
const Weight OscProb_baseline_mutoe = CreateBaselineOscWeight("mutoe");
const Weight OscProb_baseline_etomu = CreateBaselineOscWeight("etomu");
const Weight OscProb_baseline_etoe = CreateBaselineOscWeight("etoe");


const Var calcRecMomentum([](const caf::SRProxy* sr)
        {
        // constexpr double c = 2.99792458e8; // worng
        // int c = 1; // GeV/c for mass
        // int counter = 0;
        if (sr->mc.nnu <= 0) return -5.0; // Avoid cosmic events

        double px = 0;
        double py = 0;
        double pz = 0;
        for( unsigned int i = 0; i < sr->vtx.elastic.fuzzyk.npng; ++i ) {
            float dirx = sr->vtx.elastic.fuzzyk.png[i].dir.x;
            float diry = sr->vtx.elastic.fuzzyk.png[i].dir.y;
            float dirz = sr->vtx.elastic.fuzzyk.png[i].dir.z;
            int pdg = sr->vtx.elastic.fuzzyk.png[i].cvnpart.pdgmax;
            float calE = sr->vtx.elastic.fuzzyk.png[i].calE;

            double mass = 0;
            switch (pdg) {
            case 2212: mass = 0.938272; break; // Proton
            case 211:  mass = 0.139570; break; // Pion
            case 13:   mass = 0.105658; break; // Muon
            case 11:   mass = 0.000511; break; // Electron
            case 111:  mass = 0.134977; break; // Pion0
            case 22:   mass = 0; break;       // Photon
            case 2112: mass = 0.939565; break; // Neutron
            case 321:  mass = 0.493677; break; // Kaon
            case 311:  mass = 0.497611; break; // Kaon0
            case 3122: mass = 1.115683; break; // Lambda
            case 3222: mass = 1.18937; break;  // Sigma+
            default:
                std::cerr << "Warning: Unknown PDG " << pdg << ", assuming mass = 0" << std::endl;
                std::cerr << "direction: " << dirx << " " << diry << " " << dirz << std::endl;
                std::cerr << "calE: " << calE << std::endl;
                std::cerr << "truePDG: " << sr->vtx.elastic.fuzzyk.png[i].truth.pdg << std::endl;
                // counter += 1;
                // std::cerr << "counter: " << counter << std::endl;
                break;
            }
            
            // double abs_p = std::sqrt( std::pow(calE,2) - std::pow(mass,2) * std::pow(c,4) ) / c;
            double abs_p = std::sqrt(std::pow(calE, 2) - std::pow(mass, 2));
	        px += abs_p * dirx;
            py += abs_p * diry;
            pz += abs_p * dirz; // only adding 1 slice together
        }
        // double theta = std::acos(py / std::sqrt(std::pow(px, 2) + std::pow(py,2) + std::pow(pz,2))); //in radians
        double norm = std::sqrt(std::pow(px, 2) + std::pow(py, 2) + std::pow(pz, 2));
	if (norm == 0) return -5.0;
    double cosTheta = (norm > 0) ? py / norm : 1.0;
	cosTheta = std::clamp(cosTheta, -1.0, 1.0);  // Prevent NaNs
	// double theta = std::acos(cosTheta);
        return cosTheta;
        });

const Var trueMomentum([](const caf::SRProxy* sr)
    {
    if (sr->mc.nnu <= 0) return -5.0; // Avoid cosmic events
    double px = sr->mc.nu[0].p.px;
    double py = sr->mc.nu[0].p.py;
    double pz = sr->mc.nu[0].p.pz;
    // double theta = std::acos(pz / std::sqrt(std::pow(px, 2) + std::pow(py,2) + std::pow(pz,2))); //in radians
    double norm = std::sqrt(std::pow(px, 2) + std::pow(py, 2) + std::pow(pz, 2));
	if (norm == 0) return -5.0;
	double cosTheta = py / norm;
    //double cosTheta = pz / std::sqrt(std::pow(px, 2) + std::pow(py, 2) + std::pow(pz, 2));
    cosTheta = std::clamp(cosTheta, -1.0, 1.0); // Ensure the value is within acos domain
    // double theta = std::acos(cosTheta);
    return cosTheta;
    });

const Cut kCosmicCut([](const caf::SRProxy* sr) // cosmic cut
        {
            return (sr->mc.nnu > 0);
        }); 

const Cut kHas1Prongs([](const caf::SRProxy* sr) // cosmic cut
        {
            return (sr->vtx.elastic.fuzzyk.npng > 0);
        });
    
const Var resolution([](const caf::SRProxy* sr)
        {
        return (trueMomentum(sr) - calcRecMomentum(sr)) / trueMomentum(sr);
        });



// const Cut numurecoCut = kNumuCC; // from 3FlavorAna/Cuts/numuCuts.h
// However, issue is there's a 5GeV cut
// const Cut kNumuCC = kNumuQuality && kNumuNCRej;
// const Cut kNumuQuality = kNumuBasicQuality && kEbelow5GeV;
const Cut numurecoCut = kNumuBasicQuality && kNumuNCRej; // to get rid of the 5GeV cut


// const Cut nuerecoCut = kNue2024FD; // from 3FlavorAna/Cuts/nueCuts2024.h // can't use it beause of 4GeV limit
// const Cut kNue2024FD = kNue2024CorePresel && kNue2024PID && kNue2024FDNearestSlice;
// const Cut kNue2024CorePresel = kNue2024FDBasicQuality &&  kNue2024CorePart;
// const Cut kNue2024FDBasicQuality =    kApplySecondAnalysisMask && k3flavor2024FDVeto && kNue2024RecoQuality;
// const Cut kNue2024CorePart = kNue2024CoreBasicEventCut && kNue2024ProngContain;
// const Cut kNue2024CoreBasicEventCut = kNueEnergy2024 > 1 && kNueEnergy2024 < 4 && kNHit > 30 && kNHit < 150 && kLongestProng > 100 && kLongestProng < 500;
// const Cut kNue2024ProngContain = kDistAllTop   > 63 && kDistAllBottom > 12 && kDistAllEast  > 12 && kDistAllWest   > 12 && kDistAllFront > 18 && kDistAllBack   > 18;

const Cut kModifiedNue2024CoreBasicEventCut = 
    kNueEnergy2024 > 1 &&  // Keep the lower bound but remove upper bound
    kNHit > 30 && kNHit < 150 &&
    kLongestProng > 100 && kLongestProng < 500;
const Cut kModifiedNue2024CorePart = kModifiedNue2024CoreBasicEventCut && kNue2024ProngContain;
const Cut nuerecoCut = kNue2024FDBasicQuality && kModifiedNue2024CorePart && kNue2024PID && kNue2024FDNearestSlice;

HistAxis ResolutionvsE("Energy (GeV)",binsE,kE, "Resolution",binsResolution,resolution);


const Var calcRecMomentum_single_particle([](const caf::SRProxy* sr)
        {
        // constexpr double c = 2.99792458e8; // worng
        // int c = 1; // GeV/c for mass
        // int counter = 0;
        if (sr->mc.nnu <= 0) return -5.0; // Avoid cosmic events

        double px = 0;
        double py = 0;
        double pz = 0;
        for( unsigned int i = 0; i < sr->vtx.elastic.fuzzyk.npng; ++i ) {
            float dirx = sr->vtx.elastic.fuzzyk.png[i].dir.x;
            float diry = sr->vtx.elastic.fuzzyk.png[i].dir.y;
            float dirz = sr->vtx.elastic.fuzzyk.png[i].dir.z;
            int pdg = sr->vtx.elastic.fuzzyk.png[i].spprongcvnpart5label.pdgmax;
            float calE = sr->vtx.elastic.fuzzyk.png[i].calE;

            double mass = 0;
            switch (pdg) {
            case 2212: mass = 0.938272; break; // Proton
            case 211:  mass = 0.139570; break; // Pion
            case 13:   mass = 0.105658; break; // Muon
            case 11:   mass = 0.000511; break; // Electron
            case 111:  mass = 0.134977; break; // Pion0
            case 22:   mass = 0; break;       // Photon
            case 2112: mass = 0.939565; break; // Neutron
            case 321:  mass = 0.493677; break; // Kaon
            case 311:  mass = 0.497611; break; // Kaon0
            case 3122: mass = 1.115683; break; // Lambda
            case 3222: mass = 1.18937; break;  // Sigma+
            default:
                std::cerr << "Warning: Unknown PDG " << pdg << ", assuming mass = 0" << std::endl;
                std::cerr << "direction: " << dirx << " " << diry << " " << dirz << std::endl;
                std::cerr << "calE: " << calE << std::endl;
                std::cerr << "truePDG: " << sr->vtx.elastic.fuzzyk.png[i].truth.pdg << std::endl;
                // counter += 1;
                // std::cerr << "counter: " << counter << std::endl;
                break;
            }
            
            // double abs_p = std::sqrt( std::pow(calE,2) - std::pow(mass,2) * std::pow(c,4) ) / c;
            double abs_p = std::sqrt(std::pow(calE, 2) - std::pow(mass, 2));
	        px += abs_p * dirx;
            py += abs_p * diry;
            pz += abs_p * dirz; // only adding 1 slice together
        }
        // double theta = std::acos(py / std::sqrt(std::pow(px, 2) + std::pow(py,2) + std::pow(pz,2))); //in radians
        double norm = std::sqrt(std::pow(px, 2) + std::pow(py, 2) + std::pow(pz, 2));
	if (norm == 0) return -5.0;
    double cosTheta = (norm > 0) ? py / norm : 1.0;
	cosTheta = std::clamp(cosTheta, -1.0, 1.0);  // Prevent NaNs
	// double theta = std::acos(cosTheta);
        return cosTheta;
        });

const Var Resolution_sp([](const caf::SRProxy* sr) //wrong
        {
        return trueMomentum(sr) - calcRecMomentum_single_particle(sr);
        });

HistAxis ResolutionvsE_sp("Energy (GeV)",binsE,kE, "Resolution",binsResolution,Resolution_sp);//wrong

// const Cut kNoCVN_sp([](const caf::SRProxy* sr)
//         {
//             for( unsigned int i = 0; i < sr->vtx.elastic.fuzzyk.npng; ++i ) {
//                 if (sr->vtx.elastic.fuzzyk.png[i].spprongcvnpart5label.pdgmax == 0) return false;
//             }
//             return true;
//         });


HistAxis reco_vs_truth("true cos(theta)",bin_cos_theta,trueMomentum, "reco cos(theta)",bin_cos_theta,calcRecMomentum);



const Var TransMomFraction([](const caf::SRProxy* sr)
{
    if (sr->vtx.elastic.fuzzyk.npng == 0) return -5.0;

    TVector3 ret;
    for (const auto& prong : sr->vtx.elastic.fuzzyk.png)
    {
        const double w = prong.calE;
        ret += w * prong.dir; // dir.x, dir.y, dir.z
    }

    ret = ret.Unit();  // Normalize the direction vector

    // Define Zenith Direction (Y-axis)
    const TVector3 ZenithDir(0.0, 1.0, 0.0);

    // // Project onto the Zenith direction
    // const double zenithProj = ret.Dot(ZenithDir);

    // // Compute transverse momentum component
    // return (ret - zenithProj * ZenithDir).Mag();

    // Compute cos(theta_zenith)
    double costheta = ret.Dot(ZenithDir);
    return costheta;
});

const Var MuonAngle([](const caf::SRProxy* sr)
{
    const TVector3 ZenithDir(0.0, 1.0, 0.0);
    // double costheta = sr->trk.kalman.tracks[idxmuonid].dir.Dot(ZenithDir); // idxmuonid not filled
    // double costheta = sr->trk.kalman.tracks[idxlongest].dir.Dot(ZenithDir); //wrong way
    if (sr->trk.kalman.tracks.size() == 0) return -5.0;
    double costheta = sr->trk.kalman.tracks[sr->trk.kalman.idxlongest].dir.Dot(ZenithDir);
    return costheta;
});

HistAxis TransMomFraction_vs_truth("true cos(theta)",bin_cos_theta,trueMomentum, "reco cos(theta)",bin_cos_theta,TransMomFraction);
HistAxis muon_vs_truth("true cos(theta)",bin_cos_theta,trueMomentum, "reco cos(theta)",bin_cos_theta,MuonAngle);

const Var Resolution_TransMomFraction([](const caf::SRProxy* sr)
        {
        return (trueMomentum(sr) - TransMomFraction(sr));
        });

// const Var resolution([](const caf::SRProxy* sr)
//         {
//         return (trueMomentum(sr) - calcRecMomentum(sr)) / trueMomentum(sr);
//         });
// const Cut kHas1Prongs([](const caf::SRProxy* sr) // cosmic cut
//         {
//             return (sr->vtx.elastic.fuzzyk.npng > 0);
//         });
const Var abs_resolution([](const caf::SRProxy* sr)
        {
        return (trueMomentum(sr) - TransMomFraction(sr));
        });

const Cut kOffAxisCut([](const caf::SRProxy* sr)
        {
            return (trueMomentum(sr) > 0.5 && TransMomFraction(sr) < -0.5 or trueMomentum(sr) < -0.5 && TransMomFraction(sr) > 0.5 && TransMomFraction(sr) != -5 && trueMomentum(sr)!= -5);
        });
const Cut nonegative5Cut([](const caf::SRProxy* sr)
        {
            if (sr->vtx.elastic.fuzzyk.npng == 0) return false;
            else return true;
        });

const Var kvertex_diff([](const caf::SRProxy* sr)
        {
        if (sr->vtx.elastic.fuzzyk.npng == 0) return -5.0;
        double diff = 0;
        if (sr->vtx.elastic.IsValid){
        float reco_vertex = float(sr->vtx.elastic.vtx.z); // reco vertex
        float true_vertex = float(sr->mc.nu[0].vtx.z); // true vertex
        diff = true_vertex - reco_vertex;
        }
        return diff;
        });

const Var krecotruedirdiff_y([](const caf::SRProxy* sr)
        {
        if (sr->vtx.elastic.fuzzyk.npng == 0) return -5.0;
        double diff = 0;
        if (sr->vtx.elastic.IsValid){
        float reco_dir = float(TransMomFraction(sr)); // reco dir
        float true_dir = float(trueMomentum(sr)); // true dir
        diff = (reco_dir - true_dir) / true_dir;
        }
        return diff;
        });
// HistAxis kdirdiff_vs_dir("true cos(theta)",bin_cos_theta,trueMomentum, "(reco-true)/true",binsResolution,krecotruedirdiff);
const Var TransMomFraction_CosTheta_x([](const caf::SRProxy* sr)
{
    if (sr->vtx.elastic.fuzzyk.npng == 0) return -5.0; // Handle empty prongs

    TVector3 ret;
    for (const auto& prong : sr->vtx.elastic.fuzzyk.png)
    {
        const double w = prong.calE;
        ret += w * prong.dir; // Energy-weighted direction
    }

    ret = ret.Unit();  // Normalize the direction vector

    // Define X-direction
    const TVector3 XDir(1.0, 0.0, 0.0);

    // Compute cos(theta_x)
    double costheta_x = ret.Dot(XDir);
    return costheta_x;
});

const Var TransMomFraction_CosTheta_z([](const caf::SRProxy* sr)
{
    if (sr->vtx.elastic.fuzzyk.npng == 0) return -5.0; // Handle empty prongs

    TVector3 ret;
    for (const auto& prong : sr->vtx.elastic.fuzzyk.png)
    {
        const double w = prong.calE;
        ret += w * prong.dir; // Energy-weighted direction
    }

    ret = ret.Unit();  // Normalize the direction vector

    // Define Z-direction
    const TVector3 ZDir(0.0, 0.0, 1.0);

    // Compute cos(theta_z)
    double costheta_z = ret.Dot(ZDir);
    return costheta_z;
});

const Var trueMomentum_CosTheta_x([](const caf::SRProxy* sr)
    {
    if (sr->mc.nnu <= 0) return -5.0; // Avoid cosmic events
    double px = sr->mc.nu[0].p.px;
    double py = sr->mc.nu[0].p.py;
    double pz = sr->mc.nu[0].p.pz;
    // double theta = std::acos(pz / std::sqrt(std::pow(px, 2) + std::pow(py,2) + std::pow(pz,2))); //in radians
    double norm = std::sqrt(std::pow(px, 2) + std::pow(py, 2) + std::pow(pz, 2));
	if (norm == 0) return -5.0;
	double cosTheta = py / norm;
    //double cosTheta = pz / std::sqrt(std::pow(px, 2) + std::pow(py, 2) + std::pow(pz, 2));
    cosTheta = std::clamp(cosTheta, -1.0, 1.0); // Ensure the value is within acos domain
    // double theta = std::acos(cosTheta);
    return cosTheta;
    });

    const Var trueMomentum_CosTheta_z([](const caf::SRProxy* sr)
    {
    if (sr->mc.nnu <= 0) return -5.0; // Avoid cosmic events
    double px = sr->mc.nu[0].p.px;
    double py = sr->mc.nu[0].p.py;
    double pz = sr->mc.nu[0].p.pz;
    // double theta = std::acos(pz / std::sqrt(std::pow(px, 2) + std::pow(py,2) + std::pow(pz,2))); //in radians
    double norm = std::sqrt(std::pow(px, 2) + std::pow(py, 2) + std::pow(pz, 2));
	if (norm == 0) return -5.0;
	double cosTheta = py / norm;
    //double cosTheta = pz / std::sqrt(std::pow(px, 2) + std::pow(py, 2) + std::pow(pz, 2));
    cosTheta = std::clamp(cosTheta, -1.0, 1.0); // Ensure the value is within acos domain
    // double theta = std::acos(cosTheta);
    return cosTheta;
    });

const Var krecotruedirdiff_x([](const caf::SRProxy* sr)
        {
        if (sr->vtx.elastic.fuzzyk.npng == 0) return -5.0;
        double diff = 0;
        if (sr->vtx.elastic.IsValid){
        float reco_dir = float(TransMomFraction_CosTheta_x(sr)); // reco dir
        float true_dir = float(trueMomentum_CosTheta_x(sr)); // true dir
        diff = (reco_dir - true_dir) / true_dir;
        }
        return diff;
        });

const Var krecotruedirdiff_z([](const caf::SRProxy* sr)
        {
        if (sr->vtx.elastic.fuzzyk.npng == 0) return -5.0;
        double diff = 0;
        if (sr->vtx.elastic.IsValid){
        float reco_dir = float(TransMomFraction_CosTheta_z(sr)); // reco dir
        float true_dir = float(trueMomentum_CosTheta_z(sr)); // true dir
        diff = (reco_dir - true_dir) / true_dir;
        }
        return diff;
        });

const Var kresolution_x([](const caf::SRProxy* sr)
        {
        return (TransMomFraction_CosTheta_x(sr) - trueMomentum_CosTheta_x(sr));
        });
const Var kresolution_y([](const caf::SRProxy* sr)
        {
        return (TransMomFraction(sr) - trueMomentum(sr));
        });
const Var kresolution_z([](const caf::SRProxy* sr)
        {
        return (TransMomFraction_CosTheta_z(sr) - trueMomentum_CosTheta_z(sr));
        });

const Var kdist_bt_vtx_recn_truth([](const caf::SRProxy* sr)
        {
        if (sr->vtx.elastic.fuzzyk.npng == 0) return -500.0;
        double diff = 0;
        if (sr->vtx.elastic.IsValid){
        float reco_vertex_x = float(sr->vtx.elastic.vtx.x); // reco vertex
        float reco_vertex_y = float(sr->vtx.elastic.vtx.y); // reco vertex
        float reco_vertex_z = float(sr->vtx.elastic.vtx.z); // reco vertex

        float true_vertex_x = float(sr->mc.nu[0].vtx.x); // true vertex
        float true_vertex_y = float(sr->mc.nu[0].vtx.y); // true vertex
        float true_vertex_z = float(sr->mc.nu[0].vtx.z); // true vertex
        diff = sqrt(pow(reco_vertex_x - true_vertex_x, 2) + pow(reco_vertex_y - true_vertex_y, 2) + pow(reco_vertex_z - true_vertex_z, 2));

        }
        return diff;
        });

const Cut kdist_bt_vtx_recn_truthCut([](const caf::SRProxy* sr)
        {
            return kdist_bt_vtx_recn_truth(sr) <= 12;
        });

const Var kRecoNumuE = kNumuE2024; // from "3FlavorAna/Vars/NueEnergy2024.h"
const Var kRecoNueE = kNueEnergy2024_2D3D; // from "3FlavorAna/Vars/NueEnergy2024.h", true energy is kE

HistAxis Numu_recoE_vs_trueE("true Energy (GeV)",binsE,kE, "reco Energy (GeV)",binsE,kRecoNumuE);
HistAxis Nue_recoE_vs_trueE("true Energy (GeV)",binsE,kE, "reco Energy (GeV)",binsE,kRecoNueE);

const Var kResolutionNumuE([](const caf::SRProxy* sr)
        {
        if (sr->vtx.elastic.fuzzyk.npng == 0) return -5.0;
        double diff = 0;
        if (sr->vtx.elastic.IsValid){
        float reco_energy = float(kRecoNumuE(sr)); // reco energy
        float true_energy = float(kE(sr)); // true energy
        diff = (reco_energy - true_energy);
        }
        return diff;
        });
const Var kResolutionNueE([](const caf::SRProxy* sr)
        {
        if (sr->vtx.elastic.fuzzyk.npng == 0) return -5.0;
        double diff = 0;
        if (sr->vtx.elastic.IsValid){
        float reco_energy = float(kRecoNueE(sr)); // reco energy
        float true_energy = float(kE(sr)); // true energy
        diff = (reco_energy - true_energy);
        }
        return diff;
        });

HistAxis Numu_recozenith_vs_angular_resolution("Angular Resolution",binsResolution,Resolution_TransMomFraction, "Zenith Angle (degrees)",bin_theta,kTheta);
HistAxis Numu_truezenith_vs_angular_resolution("Angular Resolution",binsResolution,Resolution_TransMomFraction, "Zenith Angle (degrees)",binstheta,trueMomentum);
HistAxis Numu_recozenith_vs_kdist_bt_vtx_recn_truth("Distance between reco vertex and true vertex (cm)",binsvtxdist,kdist_bt_vtx_recn_truth, "Zenith Angle (degrees)",bin_theta,kTheta);
HistAxis Numu_truezenith_vs_kdist_bt_vtx_recn_truth("Distance between reco vertex and true vertex (cm)",binsvtxdist,kdist_bt_vtx_recn_truth, "Zenith Angle (degrees)",binstheta,trueMomentum);

// truth
// HistAxis true_Costheta_vs_E("True Energy (GeV)",binsE,kE, "Cosine of Zenith Angle",binsCosTheta,kCosTheta);
HistAxis true_Costheta_vs_E("True Energy (GeV)",binsE_nue,kE, "Cosine of Zenith Angle",binsCosTheta,kCosTheta);

//reco numu
HistAxis reco_numu_Costheta_vs_E("Reco Energy (GeV)",binsE,kRecoNumuE, "Cosine of Zenith Angle",binsCosTheta,TransMomFraction);

//reco nue
HistAxis reco_nue_Costheta_vs_E("Reco Energy (GeV)",binsE_nue,kRecoNueE, "Cosine of Zenith Angle",binsCosTheta,TransMomFraction);
