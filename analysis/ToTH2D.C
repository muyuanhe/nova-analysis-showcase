#include <iostream>
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "CAFAna/Core/Spectrum.h"
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

#include "Utilities/rootlogon.C" // For setting up ROOT styles

using namespace ana;
using namespace std;

/*
void ToTH2D()
{
    // Open file and load histograms
    TFile *myfile = new TFile("event_rate_output_corrected.root", "READ"); 
    TH2D* noosc = (TH2D*)Spectrum::LoadFrom(myfile, "LoE_vs_CosT_no_osc").release()->ToTH2(18e20, kLivetime);
    TH2D* yesosc = (TH2D*)Spectrum::LoadFrom(myfile, "LoE_vs_CosT_osc").release()->ToTH2(18e20, kLivetime);
    myfile->Close();
    
    // Create a new histogram for the ratio
    TH2D* ratio = (TH2D*)yesosc->Clone("ratio");
    ratio->SetTitle("Ratio of yesosc to noosc;L/E;cos(theta)");
    ratio->Divide(noosc); // Perform element-wise division

    // Create a canvas and draw the histograms
    TCanvas *c = new TCanvas("c", "c", 800, 600);
    
    // Draw and save the no oscillation histogram
    noosc->Draw("colz");
    c->SaveAs("noosc.png");

    // Draw and save the oscillated histogram
    yesosc->Draw("colz");
    c->SaveAs("yesosc.png");

    // Draw and save the ratio histogram to see the difference
    ratio->Draw("colz");
    c->SaveAs("ratio.png");

    // Optionally, save the ratio histogram to a new file
    TFile *outFile = new TFile("oscillation_comparison.root", "RECREATE");
    ratio->Write();
    outFile->Close();

    // Clean up
    delete c;
    delete myfile;
    delete outFile;
}

*/

/*
void ToTH2D()
{
    // Open the file containing the spectra
    TFile *myfile = new TFile("event_rate_output_theta_corrected.root", "READ");

    // Define the spectrum names for easy handling
    vector<string> spectraNames = {
        "final_numu_epsee", "final_numu_epsemu", "final_numu_epsetau",
        "final_numu_epsmumu", "final_numu_epsmutau", "final_numu_epstautau",
        "final_nue_epsee", "final_nue_epsemu", "final_nue_epsetau",
        "final_nue_epsmumu", "final_nue_epsmutau", "final_nue_epstautau", 
        "final_numu_SI", "final_nue_SI", // Added SI spectra
        "LoE_vs_CosT_no_osc_nue", "LoE_vs_CosT_no_osc_numu"
    };

    const double pot = 18e20;
    
    for(const string& name : spectraNames)
    {
        // Load the spectrum and convert it to TH2D
        auto loadedSpectrum = Spectrum::LoadFrom(myfile, name.c_str()).release();
        TH2D* hist = (TH2D*)loadedSpectrum->ToTH2(pot, kLivetime);

        // Set titles for the histogram and axes
        hist->SetTitle(("2D Histogram for " + name).c_str());
        hist->GetXaxis()->SetTitle("L/E (km/GeV)");
        hist->GetYaxis()->SetTitle("theta");

        // Create a canvas and draw the histogram
        TCanvas *c = new TCanvas(("c_" + name).c_str(), ("2D Histogram for " + name).c_str(), 800, 600);
        hist->Draw("COLZ");

        // Save each 2D histogram as a PNG image
        c->SaveAs(("2D_projection_" + name + ".png").c_str());

        // Save each histogram to a new root file
        TFile *outFile = new TFile(("TH2D_" + name + ".root").c_str(), "RECREATE");
        hist->Write(name.c_str());
        outFile->Close();

        // Clean up
        delete hist;
        delete outFile;
        delete c;
        delete loadedSpectrum;
    }


    // Close the input file
    myfile->Close();
    delete myfile;
}
*/
// the script above wil make each 2D histogram into a TH2D object and save it as a root file


//later part of the script, I'll try to make some histogram divide between NSI and SI to see the difference


/*
void ToTH2D()
{
    // Open the file containing the spectra
    TFile *myfile = new TFile("event_rate_output_corrected.root", "READ");

    // Define the spectrum names for numu and nue, and the SI names for dividing
    vector<string> numuSpectraNames = {
        "final_numu_epsee", "final_numu_epsemu", "final_numu_epsetau",
        "final_numu_epsmumu", "final_numu_epsmutau", "final_numu_epstautau"
    };
    
    vector<string> nueSpectraNames = {
        "final_nue_epsee", "final_nue_epsemu", "final_nue_epsetau",
        "final_nue_epsmumu", "final_nue_epsmutau", "final_nue_epstautau"
    };

    const double pot = 18e20;

    // Load the SI spectra for numu and nue
    auto final_numu_SI = Spectrum::LoadFrom(myfile, "final_numu_SI").release();
    auto final_nue_SI = Spectrum::LoadFrom(myfile, "final_nue_SI").release();
    
    TH2D* hist_numu_SI = (TH2D*)final_numu_SI->ToTH2(pot, kLivetime);
    TH2D* hist_nue_SI = (TH2D*)final_nue_SI->ToTH2(pot, kLivetime);

    for(const string& name : numuSpectraNames)
    {
        auto loadedSpectrum = Spectrum::LoadFrom(myfile, name.c_str()).release();
        TH2D* hist_eps = (TH2D*)loadedSpectrum->ToTH2(pot, kLivetime);

        // Perform the histogram division
        TH2D* hist_ratio = (TH2D*)hist_eps->Clone(("ratio_" + name).c_str());
        hist_ratio->Divide(hist_numu_SI);

        // Set titles and labels for the ratio histogram
        hist_ratio->SetTitle(("Ratio of " + name + " to final_numu_SI").c_str());
        hist_ratio->GetXaxis()->SetTitle("L/E (km/GeV)");
        hist_ratio->GetYaxis()->SetTitle("cos(theta)");
        hist_ratio->SetMaximum(5);

        // Create a canvas and draw the ratio histogram
        TCanvas *c = new TCanvas(("c_ratio_" + name).c_str(), ("Ratio for " + name).c_str(), 800, 600);
        hist_ratio->Draw("COLZ");

        // Save the ratio histogram as a PNG and root file
        c->SaveAs(("Ratio_projection_" + name + ".png").c_str());

        TFile *outFile = new TFile(("Ratio_TH2D_" + name + ".root").c_str(), "RECREATE");
        hist_ratio->Write(("ratio_" + name).c_str());
        outFile->Close();

        // Clean up
        delete hist_eps;
        delete hist_ratio;
        delete outFile;
        delete c;
        delete loadedSpectrum;
    }

    for(const string& name : nueSpectraNames)
    {
        auto loadedSpectrum = Spectrum::LoadFrom(myfile, name.c_str()).release();
        TH2D* hist_eps = (TH2D*)loadedSpectrum->ToTH2(pot, kLivetime);

        // Perform the histogram division
        TH2D* hist_ratio = (TH2D*)hist_eps->Clone(("ratio_" + name).c_str());
        hist_ratio->Divide(hist_nue_SI);

        // Set titles and labels for the ratio histogram
        hist_ratio->SetTitle(("Ratio of " + name + " to final_nue_SI").c_str());
        hist_ratio->GetXaxis()->SetTitle("L/E (km/GeV)");
        hist_ratio->GetYaxis()->SetTitle("cos(theta)");
        hist_ratio->SetMaximum(5);

        // Create a canvas and draw the ratio histogram
        TCanvas *c = new TCanvas(("c_ratio_" + name).c_str(), ("Ratio for " + name).c_str(), 800, 600);
        hist_ratio->Draw("COLZ");

        // Save the ratio histogram as a PNG and root file
        c->SaveAs(("Ratio_projection_" + name + ".png").c_str());

        TFile *outFile = new TFile(("Ratio_TH2D_" + name + ".root").c_str(), "RECREATE");
        hist_ratio->Write(("ratio_" + name).c_str());
        outFile->Close();

        // Clean up
        delete hist_eps;
        delete hist_ratio;
        delete outFile;
        delete c;
        delete loadedSpectrum;
    }

    // Clean up SI histograms
    delete hist_numu_SI;
    delete hist_nue_SI;

    // Close the input file
    myfile->Close();
    delete myfile;
}

*/

/*
// Previous part shows the hist_divide. Below shows relative difference: (NSI-SI)/(NSI+SI)
void ToTH2D()
{
    // Open the file containing the spectra
    TFile *myfile = new TFile("event_rate_output_corrected.root", "READ"); // L/E vs costheta
    // TFile *myfile = new TFile("event_rate_output_theta_corrected.root", "READ"); // L/E vs theta
    // TFile *myfile = new TFile("event_rate_output_EvsTheta_corrected.root", "READ"); // E vs theta

    // Define the spectrum names for numu and nue
    vector<string> numuSpectraNames = {
        "final_numu_epsee", "final_numu_epsemu", "final_numu_epsetau",
        "final_numu_epsmumu", "final_numu_epsmutau", "final_numu_epstautau"
    };

    vector<string> nueSpectraNames = {
        "final_nue_epsee", "final_nue_epsemu", "final_nue_epsetau",
        "final_nue_epsmumu", "final_nue_epsmutau", "final_nue_epstautau"
    };

    const double pot = 18e20;

    // Load the SI spectra for numu and nue
    auto final_numu_SI = Spectrum::LoadFrom(myfile, "final_numu_SI").release();
    auto final_nue_SI = Spectrum::LoadFrom(myfile, "final_nue_SI").release();

    TH2D* hist_numu_SI = (TH2D*)final_numu_SI->ToTH2(pot, kLivetime);
    TH2D* hist_nue_SI = (TH2D*)final_nue_SI->ToTH2(pot, kLivetime);

    for(const string& name : numuSpectraNames)
    {
        auto loadedSpectrum = Spectrum::LoadFrom(myfile, name.c_str()).release();
        TH2D* hist_eps = (TH2D*)loadedSpectrum->ToTH2(pot, kLivetime);

        // Perform the new histogram calculation
        TH2D* hist_new_ratio = (TH2D*)hist_eps->Clone(("new_ratio_" + name).c_str()); // Copy the NSI histogram (numerator)
        hist_new_ratio->Add(hist_numu_SI, -1); // hist_NSI - hist_SI (numerator)
        TH2D* hist_sum = (TH2D*)hist_eps->Clone("hist_sum"); // Copy the NSI histogram (denominator)
        hist_sum->Add(hist_numu_SI); // hist_NSI + hist_SI
        hist_new_ratio->Divide(hist_sum); // (hist_NSI - hist_SI) / (hist_NSI + hist_SI)

        // Set titles and labels for the new ratio histogram
        // hist_new_ratio->SetTitle(("New Ratio (NSI-SI)/(NSI+SI) for " + name).c_str());
        // hist_new_ratio->GetXaxis()->SetTitle("L/E (km/GeV)"); // L/E vs theta
        hist_new_ratio->GetXaxis()->SetTitle("L/E (km/GeV)"); // E vs theta
        hist_new_ratio->GetYaxis()->SetTitle("cos(#theta)"); // L/E vs costheta
        hist_new_ratio->SetMaximum(1); // Normalized range
        hist_new_ratio->SetMinimum(-1);

        // Create a canvas and draw the new ratio histogram
        TCanvas *c = new TCanvas(("c_new_ratio_" + name).c_str(), ("New Ratio for " + name).c_str(), 800, 600);
        hist_new_ratio->Draw("COLZ");
        Simulation();

        // Save the new ratio histogram as a PNG and root file
        c->SaveAs(("New_Ratio_projection_" + name + ".png").c_str());

        TFile *outFile = new TFile(("New_Ratio_TH2D_" + name + ".root").c_str(), "RECREATE");
        hist_new_ratio->Write(("new_ratio_" + name).c_str());
        outFile->Close();

        // Clean up
        delete hist_eps;
        delete hist_new_ratio;
        delete hist_sum;
        delete outFile;
        delete c;
        delete loadedSpectrum;
    }

    for(const string& name : nueSpectraNames)
    {
        auto loadedSpectrum = Spectrum::LoadFrom(myfile, name.c_str()).release();
        TH2D* hist_eps = (TH2D*)loadedSpectrum->ToTH2(pot, kLivetime);

        // Perform the new histogram calculation
        TH2D* hist_new_ratio = (TH2D*)hist_eps->Clone(("new_ratio_" + name).c_str()); // Copy the NSI histogram (numerator)
        hist_new_ratio->Add(hist_nue_SI, -1); // hist_NSI - hist_SI
        TH2D* hist_sum = (TH2D*)hist_eps->Clone("hist_sum"); // Copy the NSI histogram (denominator)
        hist_sum->Add(hist_nue_SI); // hist_NSI + hist_SI
        hist_new_ratio->Divide(hist_sum); // (hist_NSI - hist_SI) / (hist_NSI + hist_SI)

        // Set titles and labels for the new ratio histogram
        // hist_new_ratio->SetTitle(("New Ratio (NSI-SI)/(NSI+SI) for " + name).c_str());
        // hist_new_ratio->GetXaxis()->SetTitle("L/E (km/GeV)");
        hist_new_ratio->GetXaxis()->SetTitle("L/E (km/GeV)");
        hist_new_ratio->GetYaxis()->SetTitle("cos(#theta)");
        hist_new_ratio->SetMaximum(1); // Normalized range
        hist_new_ratio->SetMinimum(-1);

        // Create a canvas and draw the new ratio histogram
        TCanvas *c = new TCanvas(("c_new_ratio_" + name).c_str(), ("New Ratio for " + name).c_str(), 800, 600);
        hist_new_ratio->Draw("COLZ");
        Simulation();

        // Save the new ratio histogram as a PNG and root file
        c->SaveAs(("New_Ratio_projection_" + name + ".png").c_str());

        TFile *outFile = new TFile(("New_Ratio_TH2D_" + name + ".root").c_str(), "RECREATE");
        hist_new_ratio->Write(("new_ratio_" + name).c_str());
        outFile->Close();

        // Clean up
        delete hist_eps;
        delete hist_new_ratio;
        delete hist_sum;
        delete outFile;
        delete c;
        delete loadedSpectrum;
    }

    // Clean up SI histograms
    delete hist_numu_SI;
    delete hist_nue_SI;

    // Close the input file
    myfile->Close();
    delete myfile;
}
*/








// script above worked great and showed some interesting patterns. I want to put them into 1D histograms to see the difference more clearly



// Section below plots NSI and SI events on the same canvas for comparison. It showed interesting features with integrated costheta
// I used MakeSlicesFrom2D.C to make 4 slices of costheta and plot them in 1D histograms to see the difference more clearly

void ToTH2D()
{
    // Open the file containing the spectra
    TFile *myfile = new TFile("event_rate_output_corrected.root", "READ");

    // Define the spectrum names for numu and nue
    vector<string> numuSpectraNames = {
        "final_numu_epsee", "final_numu_epsemu", "final_numu_epsetau",
        "final_numu_epsmumu", "final_numu_epsmutau", "final_numu_epstautau"
    };

    vector<string> nueSpectraNames = {
        "final_nue_epsee", "final_nue_epsemu", "final_nue_epsetau",
        "final_nue_epsmumu", "final_nue_epsmutau", "final_nue_epstautau"
    };

    const double pot = 18e20;

    // Load the SI spectra for numu and nue
    auto final_numu_SI = Spectrum::LoadFrom(myfile, "final_numu_SI").release();
    auto final_nue_SI = Spectrum::LoadFrom(myfile, "final_nue_SI").release();
    
    TH2D* hist_numu_SI_2D = (TH2D*)final_numu_SI->ToTH2(pot, kLivetime);
    TH2D* hist_nue_SI_2D = (TH2D*)final_nue_SI->ToTH2(pot, kLivetime);

    // Project SI histograms to 1D
    TH1D* hist_numu_SI_1D = hist_numu_SI_2D->ProjectionX("hist_numu_SI_1D");
    TH1D* hist_nue_SI_1D = hist_nue_SI_2D->ProjectionX("hist_nue_SI_1D");

    for(const string& name : numuSpectraNames)
    {
        auto loadedSpectrum = Spectrum::LoadFrom(myfile, name.c_str()).release();
        TH2D* hist_eps_2D = (TH2D*)loadedSpectrum->ToTH2(pot, kLivetime);

        // Project NSI histogram to 1D
        TH1D* hist_eps_1D = hist_eps_2D->ProjectionX(("projectionX_" + name).c_str());

        // Create a canvas to plot both SI and NSI
        TCanvas *c = new TCanvas(("c_" + name).c_str(), ("SI vs NSI for " + name).c_str(), 800, 600);

        // Set logarithmic scale for the y-axis
        c->SetLogy();
        
        // Style SI histogram
        hist_numu_SI_1D->SetLineColor(kRed);
        hist_numu_SI_1D->SetLineWidth(2);
        hist_numu_SI_1D->SetTitle(("Number of Events: SI vs NSI for " + name).c_str());
        hist_numu_SI_1D->GetXaxis()->SetTitle("L/E (km/GeV)");
        hist_numu_SI_1D->GetYaxis()->SetTitle("Number of Events");

        // Style NSI histogram
        hist_eps_1D->SetLineColor(kBlue);
        hist_eps_1D->SetLineWidth(2);

        // Draw both histograms on the same canvas
        hist_numu_SI_1D->Draw("E"); // E for error bars
        hist_eps_1D->Draw("E SAME");

        // Add a legend
        auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->AddEntry(hist_numu_SI_1D, "SI", "l");
        legend->AddEntry(hist_eps_1D, ("NSI: " + name).c_str(), "l");
        legend->Draw();

        // Save the plot as a PNG
        c->SaveAs(("SI_vs_NSI_" + name + ".png").c_str());

        // Clean up
        delete hist_eps_2D;
        delete hist_eps_1D;
        delete c;
        delete legend;
        delete loadedSpectrum;
    }

    for(const string& name : nueSpectraNames)
    {
        auto loadedSpectrum = Spectrum::LoadFrom(myfile, name.c_str()).release();
        TH2D* hist_eps_2D = (TH2D*)loadedSpectrum->ToTH2(pot, kLivetime);

        // Project NSI histogram to 1D
        TH1D* hist_eps_1D = hist_eps_2D->ProjectionX(("projectionX_" + name).c_str());

        // Create a canvas to plot both SI and NSI
        TCanvas *c = new TCanvas(("c_" + name).c_str(), ("SI vs NSI for " + name).c_str(), 800, 600);

        // Set logarithmic scale for the y-axis
        c->SetLogy();

        // Style SI histogram
        hist_nue_SI_1D->SetLineColor(kRed);
        hist_nue_SI_1D->SetLineWidth(2);
        hist_nue_SI_1D->SetTitle(("Number of Events: SI vs NSI for " + name).c_str());
        hist_nue_SI_1D->GetXaxis()->SetTitle("L/E (km/GeV)");
        hist_nue_SI_1D->GetYaxis()->SetTitle("Number of Events");

        // Style NSI histogram
        hist_eps_1D->SetLineColor(kBlue);
        hist_eps_1D->SetLineWidth(2);

        // Draw both histograms on the same canvas
        hist_nue_SI_1D->Draw("E");
        hist_eps_1D->Draw("E SAME");

        // Add a legend
        auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->AddEntry(hist_nue_SI_1D, "SI", "l");
        legend->AddEntry(hist_eps_1D, ("NSI: " + name).c_str(), "l");
        legend->Draw();

        // Save the plot as a PNG
        c->SaveAs(("SI_vs_NSI_" + name + ".png").c_str());

        // Clean up
        delete hist_eps_2D;
        delete hist_eps_1D;
        delete c;
        delete legend;
        delete loadedSpectrum;
    }

    // Clean up SI histograms
    delete hist_numu_SI_2D;
    delete hist_nue_SI_2D;
    delete hist_numu_SI_1D;
    delete hist_nue_SI_1D;

    // Close the input file
    myfile->Close();
    delete myfile;
}
