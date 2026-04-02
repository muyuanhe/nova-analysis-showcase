#include "helper.C"
#include <iostream>
#include "NNBar_cut.C"

using namespace ana;

void pointing_nue() {
    // const std::string fname = "muyuanh_reco_triggered_atm_cosmic_oct23_2024";
    const std::string fname = "muyuanh_caf_reco_triggered_atm_cosmic_Feb32025"; // made another file for Feb3, 2025
    // corresponding to pid definition: muyuanh_pid_reco_triggered_atm_cosmic_Feb32025
    SpectrumLoader loader(fname);

    // const Cut interactionType = kIsNumu && !kIsNC; // You can change this to kIsNueCC, kIsNumuCC or kIsNC as needed
    // const Cut interactionType = kIsNumu && !kIsNC; // NueCC event != kIsNueCC
    // const Cut interactionType = kIsNue && !kIsNC; 
    // const Cut interactionType = kIsNC;
    const Cut interactionType = kIsNumu && kIsNue && !kIsNC;

    // Define interaction types for each initial flavor
    const Cut interactionType_nue = kIsNue && !kIsNC;
    const Cut interactionType_numu = kIsNumu && !kIsNC;

    extern const Cut numurecoCut;
    extern const Cut nuerecoCut;

    // Create output file
    // TFile* outputFile = new TFile("numucc_output_pointing.root", "RECREATE"); //numuCC
    TFile* outputFile = new TFile("nuecc_output_pointing.root", "RECREATE"); //nueCC

     // theta of truth and reco
    /*
    // numuCC cut
    Spectrum theta_rec("zenith angle theta", binstheta, loader, calcRecMomentum, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut);
    Spectrum theta_rec_sp("zenith angle theta", binstheta, loader, calcRecMomentum_single_particle, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut);
    Spectrum theta_truth("zenith angle theta", binstheta, loader, trueMomentum, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut);
    Spectrum spec_TransMomFraction("Transverse Mom. Fraction", binstheta, loader, TransMomFraction, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut);
    Spectrum spec_MuonAngle("Muon Angle", binstheta, loader, MuonAngle, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut);

    Spectrum spec_resolution("resolution", binstheta, loader, resolution, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut);
    Spectrum spec_Resolution_sp("resolution", binstheta, loader, Resolution_sp, kCosmicCut && kHas1Prongs  && kIsNumu && numurecoCut); // single particle cvn
    Spectrum spec_abs_resolution("abs_resolution", binstheta, loader, abs_resolution, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut);
    Spectrum spec_Resolution_TransMomFraction("Resolution vs Transverse Mom. Fraction", binstheta, loader, Resolution_TransMomFraction, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut);
    
    //histaxis is loaded differently
    Spectrum spec_ResolutionvsE(loader, ResolutionvsE, kIsNumu && numurecoCut, kNoShift);
    Spectrum spec_ResolutionvsE_sp(loader, ResolutionvsE_sp, kIsNumu && numurecoCut, kNoShift);

    Spectrum spec_reco_vs_truth(loader, reco_vs_truth, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut, kNoShift);
    Spectrum spec_TransMomFraction_vs_truth(loader, TransMomFraction_vs_truth, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut, kNoShift);
    // Spectrum spec_muon_vs_truth(loader, muon_vs_truth, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut, kNoShift);    
    */

    // /*
    // nueCC cut
    Spectrum theta_rec("zenith angle theta", binstheta, loader, calcRecMomentum, kCosmicCut  && kHas1Prongs && kIsNue && nuerecoCut);
    Spectrum theta_rec_sp("zenith angle theta", binstheta, loader, calcRecMomentum_single_particle, kCosmicCut && kHas1Prongs && kIsNue && nuerecoCut);
    Spectrum theta_truth("zenith angle theta", binstheta, loader, trueMomentum, kCosmicCut && kHas1Prongs && kIsNue && nuerecoCut);
    Spectrum spec_TransMomFraction("Transverse Mom. Fraction", binstheta, loader, TransMomFraction, kCosmicCut && kHas1Prongs && kIsNue && nuerecoCut);
    // Spectrum spec_TransMomFraction("Transverse Mom. Fraction", binstheta, loader, TransMomFraction, kCosmicCut && kHas1Prongs && kIsNue);



    Spectrum spec_resolution("resolution", binstheta, loader, resolution, kCosmicCut && kHas1Prongs  && kIsNue && nuerecoCut);
    Spectrum spec_Resolution_sp("resolution", binstheta, loader, Resolution_sp, kCosmicCut && kHas1Prongs  && kIsNue && nuerecoCut); // single particle cvn
    Spectrum spec_abs_resolution("abs_resolution", binsResolution, loader, abs_resolution, kCosmicCut && kHas1Prongs  && kIsNue && nuerecoCut) ;
    Spectrum spec_Resolution_TransMomFraction("Resolution vs Transverse Mom. Fraction", binstheta, loader, Resolution_TransMomFraction, kCosmicCut && kHas1Prongs  && kIsNue && nuerecoCut);
    
    Spectrum spec_kRecoNueE("RecoNueEnergy", binsE, loader, kRecoNueE, kCosmicCut && kHas1Prongs  && kIsNue && nuerecoCut);
    Spectrum spec_kResolutionNueE("ResolutionNueEnergy", binsres_nueE, loader, kResolutionNueE, kCosmicCut && kHas1Prongs  && kIsNue && nuerecoCut);

    
    Spectrum spec_ResolutionvsE(loader, ResolutionvsE, kIsNue && nuerecoCut, kNoShift);
    Spectrum spec_ResolutionvsE_sp(loader, ResolutionvsE_sp, kIsNue && nuerecoCut, kNoShift);
    Spectrum spec_reco_vs_truth(loader, reco_vs_truth, kCosmicCut && kHas1Prongs && kIsNue && nuerecoCut, kNoShift);
    Spectrum spec_TransMomFraction_vs_truth(loader, TransMomFraction_vs_truth, kCosmicCut && kHas1Prongs && kIsNue && nuerecoCut, kNoShift);
    Spectrum spec_Nue_recoE_vs_trueE(loader, Nue_recoE_vs_trueE, kCosmicCut && kHas1Prongs && kIsNue && nuerecoCut, kNoShift);

    Spectrum spec_true_Costheta_vs_E(loader, true_Costheta_vs_E, kCosmicCut && kHas1Prongs && kIsNue && nuerecoCut, kNoShift);
    Spectrum spec_reco_nue_Costheta_vs_E(loader, reco_nue_Costheta_vs_E, kCosmicCut && kHas1Prongs && kIsNue && nuerecoCut, kNoShift);

    // */

    loader.Go();

    const double pot = 18e20;


    // Override livetime for no oscillation spectrum
    theta_rec.OverrideLivetime(pot);
    theta_truth.OverrideLivetime(pot);
    spec_TransMomFraction.OverrideLivetime(pot);
    // spec_MuonAngle.OverrideLivetime(pot);
    spec_resolution.OverrideLivetime(pot);
    spec_abs_resolution.OverrideLivetime(pot);
    spec_ResolutionvsE.OverrideLivetime(pot);
    spec_Resolution_TransMomFraction.OverrideLivetime(pot);
    spec_kRecoNueE.OverrideLivetime(pot);
    spec_kResolutionNueE.OverrideLivetime(pot);

    theta_rec_sp.OverrideLivetime(pot);
    spec_reco_vs_truth.OverrideLivetime(pot);
    spec_TransMomFraction_vs_truth.OverrideLivetime(pot);
    // spec_muon_vs_truth.OverrideLivetime(pot);
    spec_Nue_recoE_vs_trueE.OverrideLivetime(pot);

    spec_true_Costheta_vs_E.OverrideLivetime(pot);
    spec_reco_nue_Costheta_vs_E.OverrideLivetime(pot);



    theta_rec.SaveTo(outputFile, "theta_rec");
    theta_truth.SaveTo(outputFile, "theta_truth");
    spec_TransMomFraction.SaveTo(outputFile, "spec_TransMomFraction");
    // spec_MuonAngle.SaveTo(outputFile, "MuonAngle");
    spec_resolution.SaveTo(outputFile, "spec_resolution");
    spec_abs_resolution.SaveTo(outputFile, "spec_abs_resolution");
    spec_ResolutionvsE.SaveTo(outputFile, "spec_ResolutionvsE");
    spec_Resolution_sp.SaveTo(outputFile, "spec_Resolution_sp");
    spec_ResolutionvsE_sp.SaveTo(outputFile, "spec_ResolutionvsE_sp");
    spec_Resolution_TransMomFraction.SaveTo(outputFile, "spec_Resolution_TransMomFraction");
    spec_kRecoNueE.SaveTo(outputFile, "spec_kRecoNueE");
    spec_kResolutionNueE.SaveTo(outputFile, "spec_kResolutionNueE");

    spec_reco_vs_truth.SaveTo(outputFile, "spec_reco_vs_truth");
    spec_TransMomFraction_vs_truth.SaveTo(outputFile, "spec_TransMomFraction_vs_truth");
    // spec_muon_vs_truth.SaveTo(outputFile, "spec_muon_vs_truth");
    spec_Nue_recoE_vs_trueE.SaveTo(outputFile, "spec_Nue_recoE_vs_trueE");

    spec_true_Costheta_vs_E.SaveTo(outputFile, "spec_true_Costheta_vs_E");
    spec_reco_nue_Costheta_vs_E.SaveTo(outputFile, "spec_reco_nue_Costheta_vs_E");  



    outputFile->Close();

    std::cout << "All spectra saved to nuecc_output_pointing.root for inspection." << std::endl;
}
