#include "helper.C"
#include <iostream>
#include "NNBar_cut.C"

using namespace ana;

void event_rate() {
    const std::string fname = "muyuanh_reco_triggered_atm_cosmic_oct23_2024";
    SpectrumLoader loader(fname);

    // Define interaction types for each initial flavor
    const Cut interactionType_nue = kIsNue && !kIsNC;
    const Cut interactionType_numu = kIsNumu && !kIsNC;

    // Create output file
    TFile* outputFile = new TFile("event_rate_output_EvsTheta_corrected.root", "RECREATE");

    // No oscillation spectrum for both nue and numu
    Spectrum LoE_vs_CosT_no_osc_nue(loader, EvsTheta, interactionType_nue);
    Spectrum LoE_vs_CosT_no_osc_numu(loader, EvsTheta, interactionType_numu);

    // Define spectra with appropriate interaction types for each oscillation transition and weight
    Spectrum LoE_vs_CosT_osc_mutomu_epsee(loader, EvsTheta, interactionType_numu, kNoShift, OscProb_mutomu_epsee);
    Spectrum LoE_vs_CosT_osc_mutomu_epsemu(loader, EvsTheta, interactionType_numu, kNoShift, OscProb_mutomu_epsemu);
    Spectrum LoE_vs_CosT_osc_mutomu_epsetau(loader, EvsTheta, interactionType_numu, kNoShift, OscProb_mutomu_epsetau);
    Spectrum LoE_vs_CosT_osc_mutomu_epsmumu(loader, EvsTheta, interactionType_numu, kNoShift, OscProb_mutomu_epsmumu);
    Spectrum LoE_vs_CosT_osc_mutomu_epsmutau(loader, EvsTheta, interactionType_numu, kNoShift, OscProb_mutomu_epsmutau);
    Spectrum LoE_vs_CosT_osc_mutomu_epstautau(loader, EvsTheta, interactionType_numu, kNoShift, OscProb_mutomu_epstautau);

    Spectrum LoE_vs_CosT_osc_mutoe_epsee(loader, EvsTheta, interactionType_numu, kNoShift, OscProb_mutoe_epsee);
    Spectrum LoE_vs_CosT_osc_mutoe_epsemu(loader, EvsTheta, interactionType_numu, kNoShift, OscProb_mutoe_epsemu);
    Spectrum LoE_vs_CosT_osc_mutoe_epsetau(loader, EvsTheta, interactionType_numu, kNoShift, OscProb_mutoe_epsetau);
    Spectrum LoE_vs_CosT_osc_mutoe_epsmumu(loader, EvsTheta, interactionType_numu, kNoShift, OscProb_mutoe_epsmumu);
    Spectrum LoE_vs_CosT_osc_mutoe_epsmutau(loader, EvsTheta, interactionType_numu, kNoShift, OscProb_mutoe_epsmutau);
    Spectrum LoE_vs_CosT_osc_mutoe_epstautau(loader, EvsTheta, interactionType_numu, kNoShift, OscProb_mutoe_epstautau);

    Spectrum LoE_vs_CosT_osc_etomu_epsee(loader, EvsTheta, interactionType_nue, kNoShift, OscProb_etomu_epsee);
    Spectrum LoE_vs_CosT_osc_etomu_epsemu(loader, EvsTheta, interactionType_nue, kNoShift, OscProb_etomu_epsemu);
    Spectrum LoE_vs_CosT_osc_etomu_epsetau(loader, EvsTheta, interactionType_nue, kNoShift, OscProb_etomu_epsetau);
    Spectrum LoE_vs_CosT_osc_etomu_epsmumu(loader, EvsTheta, interactionType_nue, kNoShift, OscProb_etomu_epsmumu);
    Spectrum LoE_vs_CosT_osc_etomu_epsmutau(loader, EvsTheta, interactionType_nue, kNoShift, OscProb_etomu_epsmutau);
    Spectrum LoE_vs_CosT_osc_etomu_epstautau(loader, EvsTheta, interactionType_nue, kNoShift, OscProb_etomu_epstautau);

    Spectrum LoE_vs_CosT_osc_etoe_epsee(loader, EvsTheta, interactionType_nue, kNoShift, OscProb_etoe_epsee);
    Spectrum LoE_vs_CosT_osc_etoe_epsemu(loader, EvsTheta, interactionType_nue, kNoShift, OscProb_etoe_epsemu);
    Spectrum LoE_vs_CosT_osc_etoe_epsetau(loader, EvsTheta, interactionType_nue, kNoShift, OscProb_etoe_epsetau);
    Spectrum LoE_vs_CosT_osc_etoe_epsmumu(loader, EvsTheta, interactionType_nue, kNoShift, OscProb_etoe_epsmumu);
    Spectrum LoE_vs_CosT_osc_etoe_epsmutau(loader, EvsTheta, interactionType_nue, kNoShift, OscProb_etoe_epsmutau);
    Spectrum LoE_vs_CosT_osc_etoe_epstautau(loader, EvsTheta, interactionType_nue, kNoShift, OscProb_etoe_epstautau);

    // Define oscillation with SI parameter 
    Spectrum LoE_vs_CosT_osc_mutomu_SI(loader, EvsTheta, interactionType_numu, kNoShift, OscProb_baseline_mutomu);
    Spectrum LoE_vs_CosT_osc_mutoe_SI(loader, EvsTheta, interactionType_numu, kNoShift, OscProb_baseline_mutoe);
    Spectrum LoE_vs_CosT_osc_etomu_SI(loader, EvsTheta, interactionType_nue, kNoShift, OscProb_baseline_etomu);
    Spectrum LoE_vs_CosT_osc_etoe_SI(loader, EvsTheta, interactionType_nue, kNoShift, OscProb_baseline_etoe);

    // Execute loader and set livetime after loader.Go()
    loader.Go();

    // Create combined final-state spectra for numu and nue

    // Separate final states by NSI parameter for numu and nue
    Spectrum final_numu_epsee = LoE_vs_CosT_osc_mutomu_epsee;
    Spectrum final_numu_epsemu = LoE_vs_CosT_osc_mutomu_epsemu;
    Spectrum final_numu_epsetau = LoE_vs_CosT_osc_mutomu_epsetau;
    Spectrum final_numu_epsmumu = LoE_vs_CosT_osc_mutomu_epsmumu;
    Spectrum final_numu_epsmutau = LoE_vs_CosT_osc_mutomu_epsmutau;
    Spectrum final_numu_epstautau = LoE_vs_CosT_osc_mutomu_epstautau;

    // Final states for nue
    Spectrum final_nue_epsee = LoE_vs_CosT_osc_mutoe_epsee;
    Spectrum final_nue_epsemu = LoE_vs_CosT_osc_mutoe_epsemu;
    Spectrum final_nue_epsetau = LoE_vs_CosT_osc_mutoe_epsetau;
    Spectrum final_nue_epsmumu = LoE_vs_CosT_osc_mutoe_epsmumu;
    Spectrum final_nue_epsmutau = LoE_vs_CosT_osc_mutoe_epsmutau;
    Spectrum final_nue_epstautau = LoE_vs_CosT_osc_mutoe_epstautau;

    // Final states for SI
    Spectrum final_numu_SI = LoE_vs_CosT_osc_mutomu_SI; // mu->mu correct
    Spectrum final_nue_SI = LoE_vs_CosT_osc_mutoe_SI; // mu->e incorrect


    






    const double pot = 18e20;


    // Override livetime for no oscillation spectrum
    LoE_vs_CosT_no_osc_nue.OverrideLivetime(pot);
    LoE_vs_CosT_no_osc_numu.OverrideLivetime(pot);

    // Set livetime for each oscillated spectrum explicitly
    LoE_vs_CosT_osc_mutomu_epsee.OverrideLivetime(pot);
    LoE_vs_CosT_osc_mutomu_epsemu.OverrideLivetime(pot);
    LoE_vs_CosT_osc_mutomu_epsetau.OverrideLivetime(pot);
    LoE_vs_CosT_osc_mutomu_epsmumu.OverrideLivetime(pot);
    LoE_vs_CosT_osc_mutomu_epsmutau.OverrideLivetime(pot);
    LoE_vs_CosT_osc_mutomu_epstautau.OverrideLivetime(pot);

    LoE_vs_CosT_osc_mutoe_epsee.OverrideLivetime(pot);
    LoE_vs_CosT_osc_mutoe_epsemu.OverrideLivetime(pot);
    LoE_vs_CosT_osc_mutoe_epsetau.OverrideLivetime(pot);
    LoE_vs_CosT_osc_mutoe_epsmumu.OverrideLivetime(pot);
    LoE_vs_CosT_osc_mutoe_epsmutau.OverrideLivetime(pot);
    LoE_vs_CosT_osc_mutoe_epstautau.OverrideLivetime(pot);

    LoE_vs_CosT_osc_etomu_epsee.OverrideLivetime(pot);
    LoE_vs_CosT_osc_etomu_epsemu.OverrideLivetime(pot);
    LoE_vs_CosT_osc_etomu_epsetau.OverrideLivetime(pot);
    LoE_vs_CosT_osc_etomu_epsmumu.OverrideLivetime(pot);
    LoE_vs_CosT_osc_etomu_epsmutau.OverrideLivetime(pot);
    LoE_vs_CosT_osc_etomu_epstautau.OverrideLivetime(pot);

    LoE_vs_CosT_osc_etoe_epsee.OverrideLivetime(pot);
    LoE_vs_CosT_osc_etoe_epsemu.OverrideLivetime(pot);
    LoE_vs_CosT_osc_etoe_epsetau.OverrideLivetime(pot);
    LoE_vs_CosT_osc_etoe_epsmumu.OverrideLivetime(pot);
    LoE_vs_CosT_osc_etoe_epsmutau.OverrideLivetime(pot);
    LoE_vs_CosT_osc_etoe_epstautau.OverrideLivetime(pot);

    // Set livetime for SI
    LoE_vs_CosT_osc_mutomu_SI.OverrideLivetime(pot);
    LoE_vs_CosT_osc_mutoe_SI.OverrideLivetime(pot);
    LoE_vs_CosT_osc_etomu_SI.OverrideLivetime(pot);
    LoE_vs_CosT_osc_etoe_SI.OverrideLivetime(pot);


    final_numu_epsee.OverrideLivetime(pot);
    final_numu_epsemu.OverrideLivetime(pot);
    final_numu_epsetau.OverrideLivetime(pot);
    final_numu_epsmumu.OverrideLivetime(pot);
    final_numu_epsmutau.OverrideLivetime(pot);
    final_numu_epstautau.OverrideLivetime(pot);

    final_nue_epsee.OverrideLivetime(pot);
    final_nue_epsemu.OverrideLivetime(pot);
    final_nue_epsetau.OverrideLivetime(pot);
    final_nue_epsmumu.OverrideLivetime(pot);
    final_nue_epsmutau.OverrideLivetime(pot);
    final_nue_epstautau.OverrideLivetime(pot);

    final_numu_SI.OverrideLivetime(pot);
    final_nue_SI.OverrideLivetime(pot);


    // Add additional contributions after data loading, can not += before loader.GO(), otherwise nothing to add on to since spectrum is empty
    final_numu_epsee += LoE_vs_CosT_osc_etomu_epsee;
    final_numu_epsemu += LoE_vs_CosT_osc_etomu_epsemu;
    final_numu_epsetau += LoE_vs_CosT_osc_etomu_epsetau;
    final_numu_epsmumu += LoE_vs_CosT_osc_etomu_epsmumu;
    final_numu_epsmutau += LoE_vs_CosT_osc_etomu_epsmutau;
    final_numu_epstautau += LoE_vs_CosT_osc_etomu_epstautau;

    final_nue_epsee += LoE_vs_CosT_osc_etoe_epsee;
    final_nue_epsemu += LoE_vs_CosT_osc_etoe_epsemu;
    final_nue_epsetau += LoE_vs_CosT_osc_etoe_epsetau;
    final_nue_epsmumu += LoE_vs_CosT_osc_etoe_epsmumu;
    final_nue_epsmutau += LoE_vs_CosT_osc_etoe_epsmutau;
    final_nue_epstautau += LoE_vs_CosT_osc_etoe_epstautau;

    // add additional contributions for SI
    final_numu_SI += LoE_vs_CosT_osc_etomu_SI;
    final_nue_SI += LoE_vs_CosT_osc_etoe_SI;

    // Save each spectrum to the output file
    LoE_vs_CosT_no_osc_nue.SaveTo(outputFile, "LoE_vs_CosT_no_osc_nue");
    LoE_vs_CosT_no_osc_numu.SaveTo(outputFile, "LoE_vs_CosT_no_osc_numu");

    LoE_vs_CosT_osc_mutomu_epsee.SaveTo(outputFile, "LoE_vs_CosT_osc_mutomu_epsee");
    LoE_vs_CosT_osc_mutomu_epsemu.SaveTo(outputFile, "LoE_vs_CosT_osc_mutomu_epsemu");
    LoE_vs_CosT_osc_mutomu_epsetau.SaveTo(outputFile, "LoE_vs_CosT_osc_mutomu_epsetau");
    LoE_vs_CosT_osc_mutomu_epsmumu.SaveTo(outputFile, "LoE_vs_CosT_osc_mutomu_epsmumu");
    LoE_vs_CosT_osc_mutomu_epsmutau.SaveTo(outputFile, "LoE_vs_CosT_osc_mutomu_epsmutau");
    LoE_vs_CosT_osc_mutomu_epstautau.SaveTo(outputFile, "LoE_vs_CosT_osc_mutomu_epstautau");

    LoE_vs_CosT_osc_mutoe_epsee.SaveTo(outputFile, "LoE_vs_CosT_osc_mutoe_epsee");
    LoE_vs_CosT_osc_mutoe_epsemu.SaveTo(outputFile, "LoE_vs_CosT_osc_mutoe_epsemu");
    LoE_vs_CosT_osc_mutoe_epsetau.SaveTo(outputFile, "LoE_vs_CosT_osc_mutoe_epsetau");
    LoE_vs_CosT_osc_mutoe_epsmumu.SaveTo(outputFile, "LoE_vs_CosT_osc_mutoe_epsmumu");
    LoE_vs_CosT_osc_mutoe_epsmutau.SaveTo(outputFile, "LoE_vs_CosT_osc_mutoe_epsmutau");
    LoE_vs_CosT_osc_mutoe_epstautau.SaveTo(outputFile, "LoE_vs_CosT_osc_mutoe_epstautau");

    LoE_vs_CosT_osc_etomu_epsee.SaveTo(outputFile, "LoE_vs_CosT_osc_etomu_epsee");
    LoE_vs_CosT_osc_etomu_epsemu.SaveTo(outputFile, "LoE_vs_CosT_osc_etomu_epsemu");
    LoE_vs_CosT_osc_etomu_epsetau.SaveTo(outputFile, "LoE_vs_CosT_osc_etomu_epsetau");
    LoE_vs_CosT_osc_etomu_epsmumu.SaveTo(outputFile, "LoE_vs_CosT_osc_etomu_epsmumu");
    LoE_vs_CosT_osc_etomu_epsmutau.SaveTo(outputFile, "LoE_vs_CosT_osc_etomu_epsmutau");
    LoE_vs_CosT_osc_etomu_epstautau.SaveTo(outputFile, "LoE_vs_CosT_osc_etomu_epstautau");

    LoE_vs_CosT_osc_etoe_epsee.SaveTo(outputFile, "LoE_vs_CosT_osc_etoe_epsee");
    LoE_vs_CosT_osc_etoe_epsemu.SaveTo(outputFile, "LoE_vs_CosT_osc_etoe_epsemu");
    LoE_vs_CosT_osc_etoe_epsetau.SaveTo(outputFile, "LoE_vs_CosT_osc_etoe_epsetau");
    LoE_vs_CosT_osc_etoe_epsmumu.SaveTo(outputFile, "LoE_vs_CosT_osc_etoe_epsmumu");
    LoE_vs_CosT_osc_etoe_epsmutau.SaveTo(outputFile, "LoE_vs_CosT_osc_etoe_epsmutau");
    LoE_vs_CosT_osc_etoe_epstautau.SaveTo(outputFile, "LoE_vs_CosT_osc_etoe_epstautau");

    // Save final state spectra to the output file
    final_numu_epsee.SaveTo(outputFile, "final_numu_epsee");
    final_numu_epsemu.SaveTo(outputFile, "final_numu_epsemu");
    final_numu_epsetau.SaveTo(outputFile, "final_numu_epsetau");
    final_numu_epsmumu.SaveTo(outputFile, "final_numu_epsmumu");
    final_numu_epsmutau.SaveTo(outputFile, "final_numu_epsmutau");
    final_numu_epstautau.SaveTo(outputFile, "final_numu_epstautau");

    final_nue_epsee.SaveTo(outputFile, "final_nue_epsee");
    final_nue_epsemu.SaveTo(outputFile, "final_nue_epsemu");
    final_nue_epsetau.SaveTo(outputFile, "final_nue_epsetau");
    final_nue_epsmumu.SaveTo(outputFile, "final_nue_epsmumu");
    final_nue_epsmutau.SaveTo(outputFile, "final_nue_epsmutau");
    final_nue_epstautau.SaveTo(outputFile, "final_nue_epstautau");

    LoE_vs_CosT_osc_mutomu_SI.SaveTo(outputFile, "LoE_vs_CosT_osc_mutomu_SI");
    LoE_vs_CosT_osc_mutoe_SI.SaveTo(outputFile, "LoE_vs_CosT_osc_mutoe_SI");
    LoE_vs_CosT_osc_etomu_SI.SaveTo(outputFile, "LoE_vs_CosT_osc_etomu_SI");
    LoE_vs_CosT_osc_etoe_SI.SaveTo(outputFile, "LoE_vs_CosT_osc_etoe_SI");

    // Save final state spectra for SI
    final_numu_SI.SaveTo(outputFile, "final_numu_SI");
    final_nue_SI.SaveTo(outputFile, "final_nue_SI");

    // Close the output file
    outputFile->Close();

    std::cout << "All spectra saved to event_rate_output_theta_corrected.root for inspection." << std::endl;
}

    

