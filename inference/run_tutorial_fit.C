#include "CAFAna/Analysis/Calcs.h"
#include "CAFAna/Analysis/Exposures.h"
#include "CAFAna/Fit/Fit.h"
#include "CAFAna/Analysis/Plots.h"
#include "3FlavorAna/Plotting/NuePlotStyle.h"
#include "CAFAna/Analysis/Style.h"
#include "CAFAna/Experiment/Dmsq32Constraint.h"
#include "CAFAna/Experiment/MultiExperiment.h"
#include "CAFAna/Experiment/ReactorExperiment.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/FC/FCSurface.h"
#include "CAFAna/Fit/FrequentistSurface.h"
#include "3FlavorAna/Prediction/PredictionSystJoint2018.h"
#include "CAFAna/Prediction/PredictionCombinePeriods.h"
#include "3FlavorAna/Systs/NumuSysts.h"
#include "CAFAna/Vars/FitVars.h"
#include "OscLib/IOscCalc.h"
#include "Utilities/rootlogon.C"

//#include "../joint_fit_2020_loader_tools.h"
#include "3FlavorAna/Ana2020/joint_fit_2020_loader_tools.h"

#include "TCanvas.h"
#include "TBox.h"
#include "TColor.h"
#include "TGraph.h"
#include "TVectorD.h"
#include "TF1.h"
#include "TLegend.h"
#include "TText.h"
#include "TLatex.h"
#include "TPad.h"
#include "TLine.h"
#include "TMarker.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TGaxis.h"

#include <algorithm>
#include <vector>
#include <string>
#include <fstream>

using namespace ana;


void PrettifyTGraph(TGraph* gr, const IFitVar* par, int col, double xmin = 0.0, double xmax = 0.0)
{
  std::string title = "Profiled slice;";
  title += par->LatexName();
  title += ";#chi^{2}";
  gr->SetTitle(title.c_str());
  gr->SetLineColor(col);
  gr->SetLineWidth(3);

  if(xmax != xmin)
    gr->GetXaxis()->SetLimits(xmin, xmax);
}

void run_tutorial_fit(TString options="nueOnly_fake2019_FHCOnly",
                      TString opt="demo1",
                      bool simple = true)
{
  bool corrSysts = false;
  bool nueOnly = options.Contains("nueOnly");
  bool numuOnly = options.Contains("numuOnly");
  bool joint = options.Contains("joint");
  assert (nueOnly || numuOnly || joint);

  bool FHCOnly = options.Contains("FHCOnly");
  bool RHCOnly = options.Contains("RHCOnly");
  bool both = options.Contains("both");
  assert (FHCOnly || RHCOnly || both);

  bool fake2019 = options.Contains("fake2019");
  bool mockData = options.Contains("mockData"); 

  auto suffix = options;
  if(corrSysts) suffix+="_systs";
  else suffix+="_stats";

  TString outdir = "./";
  TString outfilename   (outdir +  "bestfits_"    + suffix);

  //////////////////////////////////////////////////
  // Load Nue and Numu experiments
  //////////////////////////////////////////////////
  //need numu only for prestage seeds

  std::string decomp = ""; // make Pt extrap all the time
  std::vector <ana::predictions> preds;

  // This is a little bit dumb :) 
  if(simple)
  {
    std::vector <ana::predictions> preds_tmp= LoadPredictions (corrSysts, false, decomp, nueOnly || joint, numuOnly || joint, FHCOnly || both, RHCOnly || both);
    preds.push_back(preds_tmp[0]);
    preds.push_back(preds_tmp[1]);
    preds.push_back(preds_tmp[2]);
    preds.push_back(preds_tmp[6]);
  }
  else
    preds = LoadPredictions (corrSysts, false, decomp, nueOnly || joint, numuOnly || joint, FHCOnly || both, RHCOnly || both);

  std::vector <const IExperiment * > expts;
  std::vector <Spectrum * > data;
  auto calc_fake = DefaultOscCalc();
  if(fake2019) SetFakeCalc(calc_fake, 0.565, 2.48e-3, 0);
  else if(!mockData) {std::cerr << "need setting for data\n"; exit(1);}

  // Load the data & create SingleSampleExperiments
  TFile* datafile = new TFile((outdir+ "tutorial_mockdata.root").Data(), "read");
  for(int i = 0; i < int(preds.size()); ++i){
      double POT = preds[i].pot;

      if(mockData)
        data.push_back(LoadFrom<Spectrum>(datafile, preds[i].name).release());
      else
        data.push_back(GetFakeData(preds[i].pred, calc_fake, POT));//, preds[i].cos.first, preds[i].livetime);

      expts.push_back(new SingleSampleExperiment(preds[i].pred, *data[i]));//, *preds[i].cos.first, preds[i].cos.second, true));
  }

  osc::IOscCalcAdjustable* calc = DefaultOscCalc();
  SystShifts auxShifts = SystShifts::Nominal();

	std::vector <const IFitVar*> fitVars = {&kFitDeltaInPiUnits,
		                                      &kFitSinSqTheta23,
						                              &kFitDmSq32Scaled,
						                              &kFitSinSq2Theta13};

  if(opt.Contains("demo1")){
    // Creating a multi-experiment from all the SingleSampleExperiments
    auto exptThis = new MultiExperiment(expts);

    // Setting the fitter and fit!
    MinuitFitter fitter(exptThis, 
                        {
                          &kFitDeltaInPiUnits,
                          &kFitSinSqTheta23,
                          &kFitDmSq32Scaled,
                          &kFitSinSq2Theta13
                        },
                        {});
    fitter.Fit(calc, auxShifts);

    // Let's plot our best fit against the data
    for(int i = 0; i < int(preds.size()); ++i){
      DataMCComparison (*data[i], preds[i].pred, calc);
      gPad->Print(outdir + "demo1_fit_" + suffix + "_" + std::to_string(i) + ".pdf");
    }
  }

  if(opt.Contains("demo2")){
    // Adding the PDG th13 constraint
    expts.push_back(WorldReactorConstraint2019());

    // Creating a multi-experiment from all the SingleSampleExperiments & reactor constraint
    auto exptThis = new MultiExperiment(expts);

    // Create dcp profile by fitting other osc parameter for each dcp grid value
    auto profile_dcp = Profile(exptThis, calc, 
                               &kFitDeltaInPiUnits, 50, 0, 2, -1, 
                               {&kFitSinSqTheta23, &kFitDmSq32Scaled, &kFitSinSq2Theta13}
                               );

    PrettifyTGraph(profile_dcp, &kFitDeltaInPiUnits, kMagenta+1, 0, 2);

    profile_dcp->Draw("al");
    gPad->Print(outdir + "demo2_profile_dcp_notseeded" + suffix + ".pdf");

    // Now redo but with seeding th23 values
    SeedList oscSeeds(
      {
        {&kFitSinSqTheta23, {0.45, 0.5, 0.55}}
      }
    );

    // Create dcp profile by fitting other osc parameter for each dcp grid value
    auto profile_dcp_seeded = Profile(exptThis, calc, 
                               &kFitDeltaInPiUnits, 50, 0, 2, -1, 
                               {&kFitSinSqTheta23, &kFitDmSq32Scaled, &kFitSinSq2Theta13},
                               {},
                               oscSeeds);

    PrettifyTGraph(profile_dcp_seeded, &kFitDeltaInPiUnits, kMagenta+1, 0, 2);

    profile_dcp_seeded->Draw("al");
    gPad->Print(outdir + "demo2_profile_dcp_seeded" + suffix + ".pdf");
  }

  if(opt.Contains("demo3")){
    // Adding the PDG th13 constraint
    expts.push_back(WorldReactorConstraint2019());

    // Creating a multi-experiment from all the SingleSampleExperiments
    auto exptThis = new MultiExperiment(expts);

    // Create a chisq surface for two parameters
    FrequentistSurface surface_dcp_th23(exptThis, calc, 
                                        &kFitDeltaInPiUnits, 30, 0, 2,
                                        &kFitSinSqTheta23, 30, 0.3, 0.7,
                                        {&kFitDmSq32Scaled, &kFitSinSq2Theta13});

    TCanvas c;
    surface_dcp_th23.Draw();
    TH2* surf_1Sigma = Gaussian68Percent2D(surface_dcp_th23);
    TH2* surf_2Sigma = Gaussian2Sigma2D(surface_dcp_th23);
    TH2* surf_3Sigma = Gaussian3Sigma2D(surface_dcp_th23);
    surface_dcp_th23.DrawContour(surf_1Sigma, kSolid, kRed - 4);
    surface_dcp_th23.DrawContour(surf_2Sigma, kSolid, kRed - 7);
    surface_dcp_th23.DrawContour(surf_3Sigma, kSolid, kRed - 9);

    gPad->Print(outdir + "demo3_surface_notseeded_dcp_th23" + suffix + ".pdf");
 
    auto margs = surface_dcp_th23.GetProfiledHists();
    margs[0]->Draw("colz");
    gPad->Print(outdir + "demo3_surface_notseeded_dcp_th23_prof_dmsq32" + suffix + ".pdf");

    margs[1]->Draw("colz");
    gPad->Print(outdir + "demo3_surface_notseeded_dcp_th23_prof_th13" + suffix + ".pdf");
  }

  if(opt.Contains("demo4")){
    // Adding the PDG th13 constraint
    expts.push_back(WorldReactorConstraint2019());

    // Creating a multi-experiment from all the SingleSampleExperiments
    auto exptThis = new MultiExperiment(expts);


    // Generate th23 vs dm32 surface, without the seeds
    FrequentistSurface surface_th23_dm32(exptThis, calc, 
                                        &kFitSinSqTheta23, 30, 0.3, 0.7,
                                        &kFitDmSq32Scaled, 30, 2.2, 2.9,
                                        {&kFitDeltaInPiUnits, &kFitSinSq2Theta13});

    TCanvas c;
    surface_th23_dm32.Draw();
    TH2* surf_1Sigma = Gaussian68Percent2D(surface_th23_dm32);
    TH2* surf_2Sigma = Gaussian2Sigma2D(surface_th23_dm32);
    TH2* surf_3Sigma = Gaussian3Sigma2D(surface_th23_dm32);
    surface_th23_dm32.DrawContour(surf_1Sigma, kSolid, kRed - 4);
    surface_th23_dm32.DrawContour(surf_2Sigma, kSolid, kRed - 7);
    surface_th23_dm32.DrawContour(surf_3Sigma, kSolid, kRed - 9);


    gPad->Print(outdir + "demo4_surface_notseeded_th23_dm32" + suffix + ".pdf");
 
    auto margs = surface_th23_dm32.GetProfiledHists();
    margs[0]->Draw("colz");
    gPad->Print(outdir + "demo4_surface_notseeded_th23_dm32_prof_dcp" + suffix + ".pdf");

    margs[1]->Draw("colz");
    gPad->Print(outdir + "demo4_surface_notseeded_th23_dm32_prof_th13" + suffix + ".pdf");
  }

  if(opt.Contains("demo5")){
    // Adding the PDG th13 constraint
    expts.push_back(WorldReactorConstraint2019());

    // Creating a multi-experiment from all the SingleSampleExperiments
    auto exptThis = new MultiExperiment(expts);


    // Generate th23 vs dm32 surface, with the seeds
    FrequentistSurface surface_th23_dm32(exptThis, calc, 
                                        &kFitSinSqTheta23, 30, 0.3, 0.7,
                                        &kFitDmSq32Scaled, 30, 2.2, 2.9,
                                        {&kFitDeltaInPiUnits, &kFitSinSq2Theta13},
                                        {},
                                        SeedList({{&kFitDeltaInPiUnits, {0.5, 1.0, 1.5, 2.0}}}));

    TCanvas c;
    surface_th23_dm32.Draw();
    TH2* surf_1Sigma = Gaussian68Percent2D(surface_th23_dm32);
    TH2* surf_2Sigma = Gaussian2Sigma2D(surface_th23_dm32);
    TH2* surf_3Sigma = Gaussian3Sigma2D(surface_th23_dm32);
    surface_th23_dm32.DrawContour(surf_1Sigma, kSolid, kRed - 4);
    surface_th23_dm32.DrawContour(surf_2Sigma, kSolid, kRed - 7);
    surface_th23_dm32.DrawContour(surf_3Sigma, kSolid, kRed - 9);


    gPad->Print(outdir + "demo5_surface_seeded_th23_dm32" + suffix + ".pdf");
 
    auto margs = surface_th23_dm32.GetProfiledHists();
    margs[0]->Draw("colz");
    gPad->Print(outdir + "demo5_surface_seeded_th23_dm32_prof_dcp" + suffix + ".pdf");

    margs[1]->Draw("colz");
    gPad->Print(outdir + "demo5_surface_seeded_th23_dm32_prof_th13" + suffix + ".pdf");
  }

}

