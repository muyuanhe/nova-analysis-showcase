#include "helper.C"
#include <iostream>
#include "NNBar_cut.C"

using namespace ana;

void pointing_numu() {
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
    TFile* outputFile = new TFile("numucc_output_pointing_asis.root", "RECREATE"); //numuCC
    // TFile* outputFile = new TFile("numucc_output_pointing_optimal.root", "RECREATE"); //numuCC



    // theta of truth and reco

    // numuCC cut
    Spectrum theta_rec("zenith angle theta", binstheta, loader, calcRecMomentum, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut);
    Spectrum theta_rec_sp("zenith angle theta", binstheta, loader, calcRecMomentum_single_particle, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut);
    Spectrum theta_truth("zenith angle theta", binstheta, loader, trueMomentum, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut);
    Spectrum spec_TransMomFraction("Transverse Mom. Fraction", binstheta, loader, TransMomFraction, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut && nonegative5Cut); //as is
    // Spectrum spec_TransMomFraction("Transverse Mom. Fraction", binstheta, loader, TransMomFraction, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut && nonegative5Cut && kdist_bt_vtx_recn_truthCut); // optimal

    // Spectrum spec_TransMomFraction("Transverse Mom. Fraction", binstheta, loader, TransMomFraction, kCosmicCut && kHas1Prongs && kIsNumu);

    Spectrum spec_MuonAngle("Muon Angle", binstheta, loader, MuonAngle, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut);

    Spectrum spec_resolution("resolution", binstheta, loader, resolution, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut);
    Spectrum spec_Resolution_sp("resolution", binstheta, loader, Resolution_sp, kCosmicCut && kHas1Prongs  && kIsNumu && numurecoCut); // single particle cvn
    Spectrum spec_abs_resolution("abs_resolution", binsResolution, loader, abs_resolution, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut && nonegative5Cut);
    // Spectrum spec_Resolution_TransMomFraction("Resolution vs Transverse Mom. Fraction", binsResolution, loader, Resolution_TransMomFraction, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut && nonegative5Cut && kdist_bt_vtx_recn_truthCut);
    // Spectrum spec_Resolution_TransMomFraction("Resolution vs Transverse Mom. Fraction", binsResolution, loader, Resolution_TransMomFraction, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut && nonegative5Cut ); // as is
    Spectrum spec_Resolution_TransMomFraction("Resolution vs Transverse Mom. Fraction", binsResolution, loader, Resolution_TransMomFraction, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut && nonegative5Cut); // as is
    // Spectrum spec_Resolution_TransMomFraction("Resolution vs Transverse Mom. Fraction", binsResolution, loader, Resolution_TransMomFraction, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut && nonegative5Cut && kdist_bt_vtx_recn_truthCut); // optimal

    
    Spectrum spec_kvertex_diff("vertex_diff", binsresolutionvertex, loader, kvertex_diff, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut);
    
    Spectrum spec_krecotruedirdiff_x("(reco-true)/true cos(theta)", binsResolution, loader, krecotruedirdiff_x, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut); //
    Spectrum spec_krecotruedirdiff_z("(reco-true)/true cos(theta)", binsResolution, loader, krecotruedirdiff_z, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut);
    Spectrum spec_krecotruedirdiff_y("(reco-true)/true cos(theta)",binsResolution, loader, krecotruedirdiff_y, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut); // (reco-true)/true

    Spectrum spec_kresolution_x("resolution", binstheta, loader, kresolution_x, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut);
    Spectrum spec_kresolution_y("resolution", binstheta, loader, kresolution_y, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut);
    Spectrum spec_kresolution_z("resolution", binstheta, loader, kresolution_z, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut);

    Spectrum spec_kdist_bt_vtx_recn_truth("dist_bt_vtx_recn_truth", binsvtxdist, loader, kdist_bt_vtx_recn_truth, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut);
    
    Spectrum spec_kRecoNumuE("RecoNumuEnergy", binsE, loader, kRecoNumuE, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut); // as is
    // Spectrum spec_kRecoNumuE("RecoNumuEnergy", binsE, loader, kRecoNumuE, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut && kdist_bt_vtx_recn_truthCut); // optimal


    Spectrum spec_kResolutionNumuE("ResolutionNumuEnergy", binsresE, loader, kResolutionNumuE, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut); //as is resolution
    // Spectrum spec_kResolutionNumuE("ResolutionNumuEnergy", binsresE, loader, kResolutionNumuE, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut && kdist_bt_vtx_recn_truthCut); //optimal resolution
    
    // Trying to figure why singal cut so low from numureco cut 
    Spectrum s_remid("RemID", Binning::Simple(100, 0, 1), loader, kRemID, kCosmicCut && kIsNumu&& kNumuBasicQuality);
    // How many fail remid > 0.75?




    
    //////// histaxis is loaded differently //////////
    Spectrum spec_ResolutionvsE(loader, ResolutionvsE, kIsNumu && numurecoCut, kNoShift);
    Spectrum spec_ResolutionvsE_sp(loader, ResolutionvsE_sp, kIsNumu && numurecoCut, kNoShift);

    Spectrum spec_reco_vs_truth(loader, reco_vs_truth, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut, kNoShift);
    // Spectrum spec_TransMomFraction_vs_truth(loader, TransMomFraction_vs_truth, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut && nonegative5Cut && kdist_bt_vtx_recn_truthCut, kNoShift); //optimal
    Spectrum spec_TransMomFraction_vs_truth(loader, TransMomFraction_vs_truth, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut && nonegative5Cut, kNoShift); // as is
    Spectrum spec_muon_vs_truth(loader, muon_vs_truth, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut, kNoShift); 
    Spectrum spec_Numu_recoE_vs_trueE(loader, Numu_recoE_vs_trueE, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut, kNoShift);
    // Spectrum spec_kdirdiff_vs_dir(loader, kdirdiff_vs_dir, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut, kNoShift); // (reco dir - true dir) / true dir   vs true dir
    Spectrum spec_Numu_recozenith_vs_angular_resolution(loader, Numu_recozenith_vs_angular_resolution, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut, kNoShift);
    Spectrum spec_Numu_truezenith_vs_angular_resolution(loader, Numu_truezenith_vs_angular_resolution, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut, kNoShift);
    Spectrum spec_Numu_recozenith_vs_kdist_bt_vtx_recn_truth(loader, Numu_recozenith_vs_kdist_bt_vtx_recn_truth, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut, kNoShift);
    Spectrum spec_Numu_truezenith_vs_kdist_bt_vtx_recn_truth(loader, Numu_truezenith_vs_kdist_bt_vtx_recn_truth, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut, kNoShift);
    // as is
    Spectrum spec_true_Costheta_vs_E(loader, true_Costheta_vs_E, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut, kNoShift);
    Spectrum spec_reco_numu_Costheta_vs_E(loader, reco_numu_Costheta_vs_E, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut, kNoShift);
    //optimal
    // Spectrum spec_true_Costheta_vs_E(loader, true_Costheta_vs_E, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut && kdist_bt_vtx_recn_truthCut, kNoShift);
    // Spectrum spec_reco_numu_Costheta_vs_E(loader, reco_numu_Costheta_vs_E, kCosmicCut && kHas1Prongs && kIsNumu && numurecoCut && kdist_bt_vtx_recn_truthCut, kNoShift);


    loader.Go();

    const double pot = 18e20;


    // Override livetime for no oscillation spectrum
    theta_rec.OverrideLivetime(pot);
    theta_truth.OverrideLivetime(pot);
    spec_TransMomFraction.OverrideLivetime(pot);
    spec_MuonAngle.OverrideLivetime(pot);
    spec_resolution.OverrideLivetime(pot);
    spec_abs_resolution.OverrideLivetime(pot);
    spec_ResolutionvsE.OverrideLivetime(pot);
    spec_Resolution_TransMomFraction.OverrideLivetime(pot);
    spec_kvertex_diff.OverrideLivetime(pot);
    spec_krecotruedirdiff_x.OverrideLivetime(pot);
    spec_krecotruedirdiff_z.OverrideLivetime(pot);
    spec_krecotruedirdiff_y.OverrideLivetime(pot);
    spec_kresolution_x.OverrideLivetime(pot);
    spec_kresolution_y.OverrideLivetime(pot);
    spec_kresolution_z.OverrideLivetime(pot);
    spec_kdist_bt_vtx_recn_truth.OverrideLivetime(pot);
    spec_kRecoNumuE.OverrideLivetime(pot);
    spec_kResolutionNumuE.OverrideLivetime(pot);
    s_remid.OverrideLivetime(pot);



    theta_rec_sp.OverrideLivetime(pot);
    spec_reco_vs_truth.OverrideLivetime(pot);
    spec_TransMomFraction_vs_truth.OverrideLivetime(pot);
    spec_muon_vs_truth.OverrideLivetime(pot);
    // spec_kdirdiff_vs_dir.OverrideLivetime(pot);
    spec_Numu_recoE_vs_trueE.OverrideLivetime(pot);
    spec_Numu_recozenith_vs_angular_resolution.OverrideLivetime(pot);
    spec_Numu_truezenith_vs_kdist_bt_vtx_recn_truth.OverrideLivetime(pot);

    spec_true_Costheta_vs_E.OverrideLivetime(pot);
    spec_reco_numu_Costheta_vs_E.OverrideLivetime(pot);





    theta_rec.SaveTo(outputFile, "theta_rec");
    theta_truth.SaveTo(outputFile, "theta_truth");
    spec_TransMomFraction.SaveTo(outputFile, "spec_TransMomFraction");
    spec_MuonAngle.SaveTo(outputFile, "MuonAngle");
    spec_resolution.SaveTo(outputFile, "spec_resolution");
    spec_abs_resolution.SaveTo(outputFile, "spec_abs_resolution");
    spec_ResolutionvsE.SaveTo(outputFile, "spec_ResolutionvsE");
    spec_Resolution_sp.SaveTo(outputFile, "spec_Resolution_sp");
    spec_ResolutionvsE_sp.SaveTo(outputFile, "spec_ResolutionvsE_sp");
    spec_Resolution_TransMomFraction.SaveTo(outputFile, "spec_Resolution_TransMomFraction");
    spec_kvertex_diff.SaveTo(outputFile, "spec_kvertex_diff");
    spec_krecotruedirdiff_x.SaveTo(outputFile, "spec_krecotruedirdiff_x");
    spec_krecotruedirdiff_z.SaveTo(outputFile, "spec_krecotruedirdiff_z");
    spec_krecotruedirdiff_y.SaveTo(outputFile, "spec_krecotruedirdiff_y");
    spec_kresolution_x.SaveTo(outputFile, "spec_kresolution_x");
    spec_kresolution_y.SaveTo(outputFile, "spec_kresolution_y");
    spec_kresolution_z.SaveTo(outputFile, "spec_kresolution_z");
    spec_kdist_bt_vtx_recn_truth.SaveTo(outputFile, "spec_kdist_bt_vtx_recn_truth");
    spec_kRecoNumuE.SaveTo(outputFile, "spec_kRecoNumuE");
    spec_kResolutionNumuE.SaveTo(outputFile, "spec_kResolutionNumuE");
    s_remid.SaveTo(outputFile, "s_remid");



    spec_reco_vs_truth.SaveTo(outputFile, "spec_reco_vs_truth");
    spec_TransMomFraction_vs_truth.SaveTo(outputFile, "spec_TransMomFraction_vs_truth");
    spec_muon_vs_truth.SaveTo(outputFile, "spec_muon_vs_truth");
    // spec_kdirdiff_vs_dir.SaveTo(outputFile, "spec_kdirdiff_vs_dir");
    spec_Numu_recoE_vs_trueE.SaveTo(outputFile, "spec_Numu_recoE_vs_trueE");
    spec_Numu_recozenith_vs_angular_resolution.SaveTo(outputFile, "spec_Numu_recozenith_vs_angular_resolution");
    spec_Numu_truezenith_vs_angular_resolution.SaveTo(outputFile, "spec_Numu_truezenith_vs_angular_resolution");
    spec_Numu_recozenith_vs_kdist_bt_vtx_recn_truth.SaveTo(outputFile, "spec_Numu_recozenith_vs_kdist_bt_vtx_recn_truth");
    spec_Numu_truezenith_vs_kdist_bt_vtx_recn_truth.SaveTo(outputFile, "spec_Numu_truezenith_vs_kdist_bt_vtx_recn_truth");

    spec_true_Costheta_vs_E.SaveTo(outputFile, "spec_true_Costheta_vs_E");
    spec_reco_numu_Costheta_vs_E.SaveTo(outputFile, "spec_reco_numu_Costheta_vs_E");


    outputFile->Close();

    std::cout << "All spectra saved to numucc_output_pointing.root for inspection." << std::endl;
}
