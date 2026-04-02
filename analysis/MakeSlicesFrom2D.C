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
void MakeSlicesFrom2D()
{
    // Open the file containing the spectra
    // TFile *myfile = new TFile("event_rate_output_corrected.root", "READ");
    // TFile *myfile = new TFile("event_rate_output_theta_corrected.root", "READ"); // L/E vs theta
    TFile *myfile = new TFile("event_rate_output_EvsTheta_corrected.root", "READ"); // E vs theta

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

    // int num_y_bins_numu = hist_numu_SI_2D->GetNbinsY();
    // int num_y_bins_nue = hist_nue_SI_2D->GetNbinsY(); //Getting 100 bins, too much
    int num_y_bins_numu = 4;
    int num_y_bins_nue = 4;
    
    std::cout << "Total Y bins: " << hist_numu_SI_2D->GetNbinsY() << std::endl;

    for (const string& name : numuSpectraNames)
    {
        auto loadedSpectrum = Spectrum::LoadFrom(myfile, name.c_str()).release();
        TH2D* hist_eps_2D = (TH2D*)loadedSpectrum->ToTH2(pot, kLivetime);

        for (int i = 1; i <= num_y_bins_numu; i++)
        {
            // Project SI and NSI onto the x-axis for this y-bin
            TH1D* hist_numu_SI_slice = hist_numu_SI_2D->ProjectionX(("numu_SI_SliceY_" + to_string(i)).c_str(), 1+(i-1)*hist_numu_SI_2D->GetNbinsY()/num_y_bins_numu,i*hist_numu_SI_2D->GetNbinsY()/num_y_bins_numu);
            TH1D* hist_eps_slice = hist_eps_2D->ProjectionX(("numu_NSI_SliceY_" + name + "_SliceY_" + to_string(i)).c_str(), 1+(i-1)*hist_numu_SI_2D->GetNbinsY()/num_y_bins_numu,i*hist_numu_SI_2D->GetNbinsY()/num_y_bins_numu);

            // Create a canvas for this slice
            TCanvas* c = new TCanvas(("c_numu_" + name + "_SliceY_" + to_string(i)).c_str(), ("SI vs NSI for Slice " + to_string(i) + " of " + name).c_str(), 800, 600);
            // c->SetLogy();

            // Style the SI histogram
            hist_numu_SI_slice->SetLineColor(kRed);
            hist_numu_SI_slice->SetLineWidth(2);
            hist_numu_SI_slice->SetTitle(("Slice " + to_string(i) + ": SI vs NSI for " + name).c_str());
            // hist_numu_SI_slice->GetXaxis()->SetTitle("L/E (km/GeV)");
            hist_numu_SI_slice->GetXaxis()->SetTitle("E (GeV)");
            hist_numu_SI_slice->GetYaxis()->SetTitle("Number of Events");

            // Style the NSI histogram
            hist_eps_slice->SetLineColor(kBlue);
            hist_eps_slice->SetLineWidth(2);

            // Draw both histograms on the same canvas
            // hist_eps_slice->Draw("E");
            // hist_numu_SI_slice->Draw("E SAME");
            hist_numu_SI_slice->Draw("E");
            hist_eps_slice->Draw("E SAME");

            // Add a legend
            auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);
            legend->AddEntry(hist_numu_SI_slice, "SI", "l");
            legend->AddEntry(hist_eps_slice, ("NSI: " + name).c_str(), "l");
            legend->Draw();

            // Save the plot as a PNG
            c->SaveAs(("SliceY_" + to_string(i) + "_SI_vs_NSI_" + name + ".png").c_str());

            // Clean up
            delete hist_numu_SI_slice;
            delete hist_eps_slice;
            delete c;
            delete legend;
        }

        delete hist_eps_2D;
        delete loadedSpectrum;
    }

    for (const string& name : nueSpectraNames)
    {
        auto loadedSpectrum = Spectrum::LoadFrom(myfile, name.c_str()).release();
        TH2D* hist_eps_2D = (TH2D*)loadedSpectrum->ToTH2(pot, kLivetime);

        for (int i = 1; i <= num_y_bins_nue; i++)
        {
            // Project SI and NSI onto the x-axis for this y-bin
            TH1D* hist_nue_SI_slice = hist_nue_SI_2D->ProjectionX(("nue_SI_SliceY_" + to_string(i)).c_str(), 1+(i-1)*hist_nue_SI_2D->GetNbinsY()/num_y_bins_nue,i*hist_nue_SI_2D->GetNbinsY()/num_y_bins_nue);
            TH1D* hist_eps_slice = hist_eps_2D->ProjectionX(("nue_NSI_SliceY_" + name + "_SliceY_" + to_string(i)).c_str(), 1+(i-1)*hist_nue_SI_2D->GetNbinsY()/num_y_bins_nue,i*hist_nue_SI_2D->GetNbinsY()/num_y_bins_nue);

            // Create a canvas for this slice
            TCanvas* c = new TCanvas(("c_nue_" + name + "_SliceY_" + to_string(i)).c_str(), ("SI vs NSI for Slice " + to_string(i) + " of " + name).c_str(), 800, 600);
            // c->SetLogy();

            // Style the SI histogram
            hist_nue_SI_slice->SetLineColor(kRed);
            hist_nue_SI_slice->SetLineWidth(2);
            hist_nue_SI_slice->SetTitle(("Slice " + to_string(i) + ": SI vs NSI for " + name).c_str());
            // hist_numu_SI_slice->GetXaxis()->SetTitle("L/E (km/GeV)");
            hist_nue_SI_slice->GetXaxis()->SetTitle("E (GeV)");
            hist_nue_SI_slice->GetYaxis()->SetTitle("Number of Events");

            // Style the NSI histogram
            hist_eps_slice->SetLineColor(kBlue);
            hist_eps_slice->SetLineWidth(2);

            // Draw both histograms on the same canvas
            hist_nue_SI_slice->Draw("E");
            hist_eps_slice->Draw("E SAME");

            // Add a legend
            auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);
            legend->AddEntry(hist_nue_SI_slice, "SI", "l");
            legend->AddEntry(hist_eps_slice, ("NSI: " + name).c_str(), "l");
            legend->Draw();

            // Save the plot as a PNG
            c->SaveAs(("SliceY_" + to_string(i) + "_SI_vs_NSI_" + name + ".png").c_str());

            // Clean up
            delete hist_nue_SI_slice;
            delete hist_eps_slice;
            delete c;
            delete legend;
        }

        delete hist_eps_2D;
        delete loadedSpectrum;
    }

    // Clean up SI histograms
    delete hist_numu_SI_2D;
    delete hist_nue_SI_2D;

    // Close the input file
    myfile->Close();
    delete myfile;
}
*/

// To make NOvA format 
/*
void MakeSlicesFrom2D()
{
    gStyle->SetOptStat(0); // Disable stats box

    // Open the file
    TFile *myfile = new TFile("event_rate_output_EvsTheta_corrected.root", "READ");

    vector<string> numuSpectraNames = {
        "final_numu_epsee", "final_numu_epsemu", "final_numu_epsetau",
        "final_numu_epsmumu", "final_numu_epsmutau", "final_numu_epstautau"
    };

    vector<string> nueSpectraNames = {
        "final_nue_epsee", "final_nue_epsemu", "final_nue_epsetau",
        "final_nue_epsmumu", "final_nue_epsmutau", "final_nue_epstautau"
    };

    const double pot = 18e20;

    auto final_numu_SI = Spectrum::LoadFrom(myfile, "final_numu_SI").release();
    auto final_nue_SI = Spectrum::LoadFrom(myfile, "final_nue_SI").release();

    TH2D* hist_numu_SI_2D = (TH2D*)final_numu_SI->ToTH2(pot, kLivetime);
    TH2D* hist_nue_SI_2D = (TH2D*)final_nue_SI->ToTH2(pot, kLivetime);

    int num_y_bins = 4;

    auto drawSlices = [&](TH2D* h2_SI, vector<string>& names, bool isNumu) {
        for (const string& name : names)
        {
            auto loadedSpectrum = Spectrum::LoadFrom(myfile, name.c_str()).release();
            TH2D* h2_NSI = (TH2D*)loadedSpectrum->ToTH2(pot, kLivetime);

            for (int i = 1; i <= num_y_bins; i++)
            {
                int low = 1 + (i - 1) * h2_SI->GetNbinsY() / num_y_bins;
                int high = i * h2_SI->GetNbinsY() / num_y_bins;

                TH1D* h1_SI = h2_SI->ProjectionX(("SI_" + name + "_Y" + to_string(i)).c_str(), low, high);
                TH1D* h1_NSI = h2_NSI->ProjectionX(("NSI_" + name + "_Y" + to_string(i)).c_str(), low, high);

                TCanvas* c = new TCanvas(("c_" + name + "_Y" + to_string(i)).c_str(), "", 800, 600);

                h1_NSI->SetLineColor(kBlue);
                h1_NSI->SetFillColorAlpha(kBlue, 0.3);
                h1_NSI->SetLineWidth(2);
                // h1_NSI->SetTitle(("Slice " + to_string(i) + ": SI vs NSI for " + name).c_str());
                h1_NSI->GetXaxis()->SetTitle("E (GeV)");
                h1_NSI->GetYaxis()->SetTitle("Events / bin");
                h1_NSI->Rebin(5);

                h1_SI->SetLineColor(kRed);
                h1_SI->SetFillColorAlpha(kRed, 0.3);
                h1_SI->SetLineWidth(2);
                h1_SI->Rebin(5);

                h1_NSI->Draw("HIST");
                h1_SI->Draw("HIST SAME");

                h1_NSI->Draw("E SAME");
                h1_SI->Draw("E SAME");

                auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);
                legend->AddEntry(h1_SI, "SI", "f");
                legend->AddEntry(h1_NSI, ("NSI: " + name).c_str(), "f");
                legend->Draw();
                Simulation();

                string flavor = isNumu ? "numu" : "nue";
                c->SaveAs(("SliceY_" + to_string(i) + "_SI_vs_NSI_" + flavor + "_" + name + ".png").c_str());

                delete h1_SI;
                delete h1_NSI;
                delete c;
                delete legend;
            }

            delete h2_NSI;
            delete loadedSpectrum;
        }
    };

    drawSlices(hist_numu_SI_2D, numuSpectraNames, true);
    drawSlices(hist_nue_SI_2D, nueSpectraNames, false);

    delete hist_numu_SI_2D;
    delete hist_nue_SI_2D;
    myfile->Close();
    delete myfile;
}
*/
void MakeSlicesFrom2D()
{
    gStyle->SetOptStat(0); // Disable stats box

    // Open the file
    TFile *myfile = new TFile("event_rate_output_EvsTheta_corrected.root", "READ");

    vector<string> numuSpectraNames = {
        "final_numu_epsee", "final_numu_epsemu", "final_numu_epsetau",
        "final_numu_epsmumu", "final_numu_epsmutau", "final_numu_epstautau"
    };

    vector<string> nueSpectraNames = {
        "final_nue_epsee", "final_nue_epsemu", "final_nue_epsetau",
        "final_nue_epsmumu", "final_nue_epsmutau", "final_nue_epstautau"
    };

    std::map<std::string, std::string> nsiLabelMap = {
        {"final_numu_epsee", "#varepsilon_{ee}"},
        {"final_numu_epsemu", "#varepsilon_{e#mu}"},
        {"final_numu_epsetau", "#varepsilon_{e#tau}"},
        {"final_numu_epsmumu", "#varepsilon_{#mu#mu}"},
        {"final_numu_epsmutau", "#varepsilon_{#mu#tau}"},
        {"final_numu_epstautau", "#varepsilon_{#tau#tau}"},
        {"final_nue_epsee", "#varepsilon_{ee}"},
        {"final_nue_epsemu", "#varepsilon_{e#mu}"},
        {"final_nue_epsetau", "#varepsilon_{e#tau}"},
        {"final_nue_epsmumu", "#varepsilon_{#mu#mu}"},
        {"final_nue_epsmutau", "#varepsilon_{#mu#tau}"},
        {"final_nue_epstautau", "#varepsilon_{#tau#tau}"}
    };

    const double pot = 18e20;

    auto final_numu_SI = Spectrum::LoadFrom(myfile, "final_numu_SI").release();
    auto final_nue_SI = Spectrum::LoadFrom(myfile, "final_nue_SI").release();

    TH2D* hist_numu_SI_2D = (TH2D*)final_numu_SI->ToTH2(pot, kLivetime);
    TH2D* hist_nue_SI_2D = (TH2D*)final_nue_SI->ToTH2(pot, kLivetime);

    int num_y_bins = 4;

    auto drawSlices = [&](TH2D* h2_SI, vector<string>& names, bool isNumu) {
        for (const string& name : names)
        {
            auto loadedSpectrum = Spectrum::LoadFrom(myfile, name.c_str()).release();
            TH2D* h2_NSI = (TH2D*)loadedSpectrum->ToTH2(pot, kLivetime);

            for (int i = 1; i <= num_y_bins; i++)
            {
                int low = 1 + (i - 1) * h2_SI->GetNbinsY() / num_y_bins;
                int high = i * h2_SI->GetNbinsY() / num_y_bins;

                TH1D* h1_SI = h2_SI->ProjectionX(("SI_" + name + "_Y" + to_string(i)).c_str(), low, high);
                TH1D* h1_NSI = h2_NSI->ProjectionX(("NSI_" + name + "_Y" + to_string(i)).c_str(), low, high);

                TCanvas* c = new TCanvas(("c_" + name + "_Y" + to_string(i)).c_str(), "", 800, 600);

                h1_NSI->SetLineColor(kBlue);
                h1_NSI->SetFillColorAlpha(kBlue, 0.3);
                h1_NSI->SetLineWidth(2);
                h1_NSI->GetXaxis()->SetTitle("E (GeV)");
                h1_NSI->GetYaxis()->SetTitle("Events / bin");
                h1_NSI->Rebin(5);

                h1_SI->SetLineColor(kRed);
                h1_SI->SetFillColorAlpha(kRed, 0.3);
                h1_SI->SetLineWidth(2);
                h1_SI->Rebin(5);

                // To normalize them
                h1_SI->Scale(1.0 / h1_SI->Integral());
                h1_NSI->Scale(1.0 / h1_NSI->Integral());  
                // and prevent 0 total counts
                if (h1_SI->Integral() > 0) h1_SI->Scale(1.0 / h1_SI->Integral());
                if (h1_NSI->Integral() > 0) h1_NSI->Scale(1.0 / h1_NSI->Integral());

                h1_NSI->Draw("HIST");
                h1_SI->Draw("HIST SAME");

                h1_NSI->Draw("E SAME");
                h1_SI->Draw("E SAME");

                auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);
                legend->SetTextSize(0.045); // Larger font size
                legend->AddEntry(h1_SI, "#bf{SI}", "f");
                legend->AddEntry(h1_NSI, nsiLabelMap[name].c_str(), "f");
                // legend->AddEntry(h1_NSI, Form("#epsilon_{%s}", latexLabel.c_str()), "f");
                legend->Draw();

                Simulation();

                string flavor = isNumu ? "numu" : "nue";
                c->SaveAs(("SliceY_" + to_string(i) + "_SI_vs_NSI_" + flavor + "_" + name + ".png").c_str());

                delete h1_SI;
                delete h1_NSI;
                delete c;
                delete legend;
            }

            delete h2_NSI;
            delete loadedSpectrum;
        }
    };

    drawSlices(hist_numu_SI_2D, numuSpectraNames, true);
    drawSlices(hist_nue_SI_2D, nueSpectraNames, false);

    delete hist_numu_SI_2D;
    delete hist_nue_SI_2D;
    myfile->Close();
    delete myfile;
}