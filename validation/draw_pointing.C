// Make two simple spectrum plots

#include "CAFAna/Core/Binning.h"
#include "CAFAna/Cuts/Cuts.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Vars/Vars.h"

#include "StandardRecord/Proxy/SRProxy.h"
#include "OscLib/OscCalcPMNSOpt.h"

#include "TCanvas.h"
#include "TH2.h"
#include "TLegend.h"
#include "TFile.h"

#include "TStyle.h"
using namespace ana;

#include <iostream>
#include "Utilities/rootlogon.C"


void draw_pointing()
{
// TFile *myfile = new TFile("output_pointing.root","READ");
// TFile *myfile = new TFile("numucc_output_pointing.root","READ"); //numuCC
TFile *myfile = new TFile("numucc_output_pointing_optimal.root","READ"); //numuCC optimistic
// TFile *myfile = new TFile("numucc_output_pointing_asis.root","READ"); //numuCC as is
// TFile *myfile = new TFile("nuecc_output_pointing.root","READ"); //nueCC
//gStyle->SetOptStat(1111);
//////////////////////// Load Spectrum from TFile ////////////////////////

std::unique_ptr<Spectrum> theta_rec = Spectrum::LoadFrom(myfile, "theta_rec");
std::unique_ptr<Spectrum> theta_truth = Spectrum::LoadFrom(myfile, "theta_truth");
std::unique_ptr<Spectrum> spec_TransMomFraction = Spectrum::LoadFrom(myfile, "spec_TransMomFraction");
std::unique_ptr<Spectrum> spec_Resolution_TransMomFraction = Spectrum::LoadFrom(myfile, "spec_Resolution_TransMomFraction");

std::unique_ptr<Spectrum> spec_resolution = Spectrum::LoadFrom(myfile, "spec_resolution");
std::unique_ptr<Spectrum> spec_resolutionvsE = Spectrum::LoadFrom(myfile, "spec_ResolutionvsE");

std::unique_ptr<Spectrum> spec_reco_vs_truth = Spectrum::LoadFrom(myfile, "spec_reco_vs_truth");
std::unique_ptr<Spectrum> spec_TransMomFraction_vs_truth = Spectrum::LoadFrom(myfile, "spec_TransMomFraction_vs_truth");
std::unique_ptr<Spectrum> sped_RecoE = Spectrum::LoadFrom(myfile, "spec_kRecoNumuE"); //numu
// std::unique_ptr<Spectrum> sped_RecoE = Spectrum::LoadFrom(myfile, "spec_kRecoNueE"); //nue
std::unique_ptr<Spectrum> spec_resolutionE = Spectrum::LoadFrom(myfile, "spec_kResolutionNumuE"); // resolution vs reco energy numu
// std::unique_ptr<Spectrum> spec_resolutionE = Spectrum::LoadFrom(myfile, "spec_kResolutionNueE"); // nue
// std::unique_ptr<Spectrum> spec_kdirdiff_vs_dir = Spectrum::LoadFrom(myfile, "spec_kdirdiff_vs_dir"); // (reco dir - true dir) / true dir   vs true dir
std::unique_ptr<Spectrum> spec_spec_Numu_recozenith_vs_angular_resolution = Spectrum::LoadFrom(myfile, "spec_Numu_recozenith_vs_angular_resolution");
std::unique_ptr<Spectrum> spec_spec_Numu_truezenith_vs_angular_resolution = Spectrum::LoadFrom(myfile, "spec_Numu_truezenith_vs_angular_resolution");
std::unique_ptr<Spectrum> spec_Numu_recozenith_vs_kdist_bt_vtx_recn_truth = Spectrum::LoadFrom(myfile, "spec_Numu_recozenith_vs_kdist_bt_vtx_recn_truth");
std::unique_ptr<Spectrum> spec_Numu_truezenith_vs_kdist_bt_vtx_recn_truth = Spectrum::LoadFrom(myfile, "spec_Numu_truezenith_vs_kdist_bt_vtx_recn_truth"); 


std::unique_ptr<Spectrum> spec_true_Costheta_vs_E = Spectrum::LoadFrom(myfile, "spec_true_Costheta_vs_E"); //both numu and nue. truth
std::unique_ptr<Spectrum> spec_reco_numu_Costheta_vs_E = Spectrum::LoadFrom(myfile, "spec_reco_numu_Costheta_vs_E"); //numu
// std::unique_ptr<Spectrum> spec_reco_nue_Costheta_vs_E = Spectrum::LoadFrom(myfile, "spec_reco_nue_Costheta_vs_E"); // nue


// Note, only used for numuCC, comment out for nue
// std::unique_ptr<Spectrum> spec_muon_vs_truth = Spectrum::LoadFrom(myfile, "spec_muon_vs_truth");



const double pot = 18e20;

//////////////////////// Override Livetime ////////////////////////

double livetime = 18e20;


theta_rec->OverrideLivetime(livetime);
theta_truth->OverrideLivetime(livetime);
spec_TransMomFraction->OverrideLivetime(livetime);
spec_Resolution_TransMomFraction->OverrideLivetime(livetime);

spec_resolution->OverrideLivetime(livetime);
spec_resolutionvsE->OverrideLivetime(livetime);

spec_reco_vs_truth->OverrideLivetime(livetime);
spec_TransMomFraction_vs_truth->OverrideLivetime(livetime);
sped_RecoE->OverrideLivetime(livetime);
spec_resolutionE->OverrideLivetime(livetime); // resolution vs reco energy
// spec_kdirdiff_vs_dir->OverrideLivetime(livetime); // (reco dir - true dir) / true dir   vs true dir
spec_spec_Numu_recozenith_vs_angular_resolution->OverrideLivetime(livetime);
spec_spec_Numu_truezenith_vs_angular_resolution->OverrideLivetime(livetime);
spec_Numu_recozenith_vs_kdist_bt_vtx_recn_truth->OverrideLivetime(livetime);

spec_true_Costheta_vs_E->OverrideLivetime(livetime); //both numu and nue. truth
spec_reco_numu_Costheta_vs_E->OverrideLivetime(livetime); //numu
// spec_reco_nue_Costheta_vs_E->OverrideLivetime(livetime); // nue



// spec_muon_vs_truth->OverrideLivetime(livetime); // comment out for nuecc


//////////////////////// Convert Spectrum to TH1D ////////////////////////

TH1D *htheta_rec = theta_rec->ToTH1(pot,kLivetime);
TH1D *htheta_truth = theta_truth->ToTH1(pot,kLivetime);
TH1D *hspec_TransMomFraction = spec_TransMomFraction->ToTH1(pot,kLivetime);
TH1D *hspec_Resolution_TransMomFraction = spec_Resolution_TransMomFraction->ToTH1(pot,kLivetime);

// Normalize the histograms
hspec_TransMomFraction->Scale(1.0/hspec_TransMomFraction->Integral()); // Normalize the histogram
if (hspec_TransMomFraction->Integral() > 0) hspec_TransMomFraction->Scale(1.0 / hspec_TransMomFraction->Integral()); // Prevent 0 total counts
hspec_Resolution_TransMomFraction->Scale(1.0/hspec_Resolution_TransMomFraction->Integral()); // Normalize the histogram
if (hspec_Resolution_TransMomFraction->Integral() > 0) hspec_Resolution_TransMomFraction->Scale(1.0 / hspec_Resolution_TransMomFraction->Integral()); // Prevent 0 total counts

TH1D *hspec_resolution = spec_resolution->ToTH1(pot,kLivetime);

TH1D *hspec_RecoE = sped_RecoE->ToTH1(pot,kLivetime); // reco energy spectrum
hspec_RecoE->Rebin(5);// for nue
// hspec_RecoE->Rebin(2);// for numu
// Normalize the histogram
hspec_RecoE->Scale(1.0/hspec_RecoE->Integral()); // Normalize the histogram
if (hspec_RecoE->Integral() > 0) hspec_RecoE->Scale(1.0 / hspec_RecoE->Integral()); // Prevent 0 total counts

TH1D *hspec_resolutionE = spec_resolutionE->ToTH1(pot,kLivetime); // resolution vs reco energy
hspec_resolutionE->Rebin(5); // Rebin for better visualization
// Normalize the histogram
hspec_resolutionE->Scale(1.0/hspec_resolutionE->Integral()); // Normalize the histogram
if (hspec_resolutionE->Integral() > 0) hspec_resolutionE->Scale(1.0 / hspec_resolutionE->Integral()); // Prevent 0 total

// for 2D h1->GetXaxis()->Rebin(5)

TH2 *hspec_resolutionvsE = spec_resolutionvsE->ToTH2(pot,kLivetime);
hspec_resolutionvsE->Rebin(5); // Rebin for better visualization

TH2 *hspec_reco_vs_truth = spec_reco_vs_truth->ToTH2(pot,kLivetime);
TH2 *hspec_TransMomFraction_vs_truth = spec_TransMomFraction_vs_truth->ToTH2(pot,kLivetime);
hspec_TransMomFraction_vs_truth->RebinX(5); // Rebin for better visualization
hspec_TransMomFraction_vs_truth->RebinY(5); // Rebin for better visualization

TH2 *hspec_spec_Numu_recozenith_vs_angular_resolution = spec_spec_Numu_recozenith_vs_angular_resolution->ToTH2(pot,kLivetime);
hspec_spec_Numu_recozenith_vs_angular_resolution->RebinY(10);

TH2 *hspec_spec_Numu_truezenith_vs_angular_resolution = spec_spec_Numu_truezenith_vs_angular_resolution->ToTH2(pot,kLivetime);

TH2 *hspec_Numu_recozenith_vs_kdist_bt_vtx_recn_truth = spec_Numu_recozenith_vs_kdist_bt_vtx_recn_truth->ToTH2(pot,kLivetime);
hspec_Numu_recozenith_vs_kdist_bt_vtx_recn_truth->RebinY(10);

TH2 *hspec_Numu_truezenith_vs_kdist_bt_vtx_recn_truth = spec_Numu_truezenith_vs_kdist_bt_vtx_recn_truth->ToTH2(pot,kLivetime);

TH2 *hspec_true_Costheta_vs_E = spec_true_Costheta_vs_E->ToTH2(pot,kLivetime); //both numu and nue. truth

TH2 *hspec_reco_numu_Costheta_vs_E = spec_reco_numu_Costheta_vs_E->ToTH2(pot,kLivetime); //numu
// TH2 *hspec_reco_nue_Costheta_vs_E = spec_reco_nue_Costheta_vs_E->ToTH2(pot,kLivetime); // nue        




// TH2 *hspec_kdirdiff_vs_dir = spec_kdirdiff_vs_dir->ToTH2(pot,kLivetime); // (reco dir - true dir) / true dir   vs true dir

// TH2 *hspec_muon_vs_truth = spec_muon_vs_truth->ToTH2(pot,kLivetime); // comment out for nue

////////////////////////////// count how many events in the spectrum //////////////////////////////

// Compute the integral (sum of bin contents)
double num_entries = hspec_TransMomFraction->Integral();

std::cout << "Total number of events in the spectrum (including binning effects): " << num_entries << std::endl;
////////////////////////////// end of count how many events in the spectrum //////////////////////////////


//////////////////////// Create Slices ////////////////////////
std::vector<TH1D*> hResSlices_Numu; // Store slices
for (int i = 0; i < 6; i++) {
    // Create projection for each slice
    hResSlices_Numu.push_back(
        hspec_resolutionvsE->ProjectionY(Form("Numu_%i", i), i * 10 + 1, i * 10 + 10)
    );
}

//////////////////////// drawing histograms ////////////////////////

TCanvas *c = new TCanvas("c");
// TCanvas *c = new TCanvas("c", "Slices", 1200, 800);

/*
c->Divide(3, 2); // Divide the canvas into a 3x2 grid

for (int i = 0; i < hResSlices_Numu.size(); i++) {
    c->cd(i + 1); // Go to the ith subpad
    hResSlices_Numu[i]->SetLineColor(i + 1); // Assign different colors
    hResSlices_Numu[i]->SetTitle(Form("Slice %i (Energy Bin %i-%i)", i, i * 10, (i + 1) * 10));
    hResSlices_Numu[i]->GetXaxis()->SetTitle("Resolution");
    hResSlices_Numu[i]->GetYaxis()->SetTitle("Counts");
    hResSlices_Numu[i]->Draw("HIST");
}

c->Print("pointing_resolution_slices_grid.png");
*/


htheta_rec->SetTitle("Reconstructed Zenith Angle");
htheta_rec->Draw("hist");
c->Print("pointing_rec.png");

htheta_truth->SetTitle("True Zenith Angle");
htheta_truth->Draw("hist");
c->Print("pointing_truth.png");

// hspec_TransMomFraction->SetTitle("Transverse Momentum Fraction");
hspec_TransMomFraction->Draw("hist");
hspec_TransMomFraction->GetXaxis()->SetTitle("Zenith Angle (cos(#theta))");
Simulation();
c->Print("pointing_TransMomFraction.png");

// hspec_Resolution_TransMomFraction->SetTitle("Resolution vs Transverse Momentum Fraction");
hspec_Resolution_TransMomFraction->Draw("hist");
hspec_Resolution_TransMomFraction->GetXaxis()->SetTitle("Angular Resolution");
Simulation();
c->Print("pointing_Resolution_TransMomFraction.png");

// hsped_RecoE->SetTitle("Reconstructed Energy Spectrum");
hspec_RecoE->GetXaxis()->SetTitle("Reconstructed Energy (GeV)");
hspec_RecoE->Draw("hist");
Simulation();
c->Print("pointing_RecoE.png");

hspec_resolutionE->GetXaxis()->SetTitle("Reconstructed Energy Resolution");
// hspec_resolutionE->SetTitle("Reconstructed Energy Resolution");
hspec_resolutionE->Draw("hist");
Simulation();
c->Print("pointing_resolutionE.png");

hspec_resolution->SetTitle("Resolution");
hspec_resolution->Draw("hist");
c->Print("pointing_resolution.png");

hspec_resolutionvsE->SetTitle("Resolution vs Energy");
hspec_resolutionvsE->Draw("COLZ");
c->Print("pointing_resolutionvsE.png");

// hspec_reco_vs_truth->SetTitle("Reconstructed vs True Zenith Angle");
hspec_reco_vs_truth->Draw("COLZ");
Simulation();
c->Print("pointing_reco_vs_truth.png");

// hspec_TransMomFraction_vs_truth->SetTitle("Transverse Momentum Fraction vs True Zenith Angle");
hspec_TransMomFraction_vs_truth->Draw("COLZ");
hspec_TransMomFraction_vs_truth->GetXaxis()->SetTitle("True Zenith Angle (cos(#theta))");
hspec_TransMomFraction_vs_truth->GetYaxis()->SetTitle("Reconstructed Zenith Angle (cos(#theta))");
Simulation();
c->Print("pointing_TransMomFraction_vs_truth.png");

// hspec_kdirdiff_vs_dir->SetTitle("kdirdiff vs dir");
// hspec_kdirdiff_vs_dir->Draw("COLZ");
// c->Print("pointing_kdirdiff_vs_dir.png");

// hspec_muon_vs_truth->SetTitle("Muon Angle vs True Zenith Angle"); // comment out for nue
// hspec_muon_vs_truth->Draw("COLZ");
// c->Print("pointing_muon_vs_truth.png");

hspec_spec_Numu_recozenith_vs_angular_resolution->SetTitle("Reconstructed Zenith Angle vs Angular Resolution");
hspec_spec_Numu_recozenith_vs_angular_resolution->Draw("COLZ");
hspec_spec_Numu_recozenith_vs_angular_resolution->GetXaxis()->SetTitle("Angular Resolution (degrees)");
hspec_spec_Numu_recozenith_vs_angular_resolution->GetYaxis()->SetTitle("Reconstructed Zenith Angle (degrees)");
Simulation();
c->Print("pointing_Numu_recozenith_vs_angular_resolution.png");

hspec_spec_Numu_truezenith_vs_angular_resolution->SetTitle("True Zenith Angle vs Angular Resolution");
hspec_spec_Numu_truezenith_vs_angular_resolution->Draw("COLZ");
hspec_spec_Numu_truezenith_vs_angular_resolution->GetXaxis()->SetTitle("Angular Resolution (degrees)");
hspec_spec_Numu_truezenith_vs_angular_resolution->GetYaxis()->SetTitle("True Zenith Angle (degrees)");
Simulation();
c->Print("pointing_Numu_truezenith_vs_angular_resolution.png");

hspec_Numu_recozenith_vs_kdist_bt_vtx_recn_truth->SetTitle("Reconstructed Zenith Angle vs Distance between Reconstructed Vertex and True Vertex");
hspec_Numu_recozenith_vs_kdist_bt_vtx_recn_truth->Draw("COLZ");
hspec_Numu_recozenith_vs_kdist_bt_vtx_recn_truth->GetXaxis()->SetTitle("Distance between Reconstructed Vertex and True Vertex (cm)");
hspec_Numu_recozenith_vs_kdist_bt_vtx_recn_truth->GetYaxis()->SetTitle("Reconstructed Zenith Angle (degrees)");
Simulation();
c->Print("pointing_Numu_recozenith_vs_kdist_bt_vtx_recn_truth.png");

hspec_Numu_truezenith_vs_kdist_bt_vtx_recn_truth->SetTitle("True Zenith Angle vs Distance between Reconstructed Vertex and True Vertex");
hspec_Numu_truezenith_vs_kdist_bt_vtx_recn_truth->Draw("COLZ");
hspec_Numu_truezenith_vs_kdist_bt_vtx_recn_truth->GetXaxis()->SetTitle("Distance between Reconstructed Vertex and True Vertex (cm)");
hspec_Numu_truezenith_vs_kdist_bt_vtx_recn_truth->GetYaxis()->SetTitle("True Zenith Angle (degrees)");
Simulation();
c->Print("pointing_Numu_truezenith_vs_kdist_bt_vtx_recn_truth.png");

TCanvas *c1 = new TCanvas("c_CosthetaVsE");

// --------- True Cos(θ) vs E ---------
hspec_true_Costheta_vs_E->SetTitle("True Cos(#theta) vs Energy");
hspec_true_Costheta_vs_E->GetXaxis()->SetTitle("Energy (GeV)");
hspec_true_Costheta_vs_E->GetYaxis()->SetTitle("True Cos(#theta)");
hspec_true_Costheta_vs_E->Draw("COLZ");
Simulation();  // Optional annotation
c1->Print("pointing_true_Costheta_vs_E.png");

// --------- Reco Cos(θ) vs E ---------
hspec_reco_numu_Costheta_vs_E->SetTitle("Reco Cos(#theta) vs Energy");
hspec_reco_numu_Costheta_vs_E->GetXaxis()->SetTitle("Energy (GeV)");
hspec_reco_numu_Costheta_vs_E->GetYaxis()->SetTitle("Reco Cos(#theta)");
hspec_reco_numu_Costheta_vs_E->Draw("COLZ");
Simulation();  // Optional annotation
c1->Print("pointing_reco_Costheta_vs_E.png");

// --------- Save both histograms into one ROOT file ---------
TFile *fout = new TFile("pointing_numu_optimal_Costheta_vs_E_spectra.root", "RECREATE");
hspec_true_Costheta_vs_E->Write("true_Costheta_vs_E");
hspec_reco_numu_Costheta_vs_E->Write("reco_Costheta_vs_E");
fout->Close();

/*
hspec_reco_nue_Costheta_vs_E->SetTitle("Reco Costheta vs Energy"); // nue
hspec_reco_nue_Costheta_vs_E->Draw("COLZ");
hspec_reco_nue_Costheta_vs_E->GetXaxis()->SetTitle("Energy (GeV)");
hspec_reco_nue_Costheta_vs_E->GetYaxis()->SetTitle("Reco Costheta");
// Simulation();
c->Print("pointing_reco_nue_Costheta_vs_E.png"); // nue
*/


}
