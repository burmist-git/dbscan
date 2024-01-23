//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include <time.h>

using namespace std;

Int_t plot_dbscan(){

  TString fileN01;
  //
  TString objectName;
  //fileN01 = "./hist_dbscan.root";
  //fileN01 = "./hist_dbscan_Compound.root"; 
  //fileN01 = "./hist_dbscan_Aggregation.root";
  //fileN01 = "./hist_dbscan_spiral.root";
  fileN01 = "./hist_dbscan_S_sets.root";
  //fileN01 = "./hist_dbscan_worms_2d.root";

  //
  TFile *f01 = new TFile(fileN01.Data());
  //
  TGraph *gr = (TGraph*)f01->Get("gr");
  TGraph *gr_cl00 = (TGraph*)f01->Get("gr_cl00");
  TGraph *gr_cl01 = (TGraph*)f01->Get("gr_cl01");
  TGraph *gr_cl02 = (TGraph*)f01->Get("gr_cl02");
  TGraph *gr_cl03 = (TGraph*)f01->Get("gr_cl03");
  TGraph *gr_cl00_BORDER = (TGraph*)f01->Get("gr_cl00_BORDER");
  TGraph *gr_cl01_BORDER = (TGraph*)f01->Get("gr_cl01_BORDER");
  TGraph *gr_cl02_BORDER = (TGraph*)f01->Get("gr_cl02_BORDER");
  TGraph *gr_cl03_BORDER = (TGraph*)f01->Get("gr_cl03_BORDER");
  TGraph *gr_NOISE = (TGraph*)f01->Get("gr_NOISE");
  //
  TCanvas *c1 = new TCanvas("c1",fileN01.Data(),10,10,800,800);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
  //
  gr->SetMarkerColor(kBlack);
  gr->SetMarkerStyle(24);
  //
  gr_cl00->SetMarkerColor(kRed+2);
  gr_cl00->SetMarkerSize(1);
  gr_cl00->SetLineWidth(2);
  gr_cl00->SetMarkerStyle(2);
  //
  gr_cl00_BORDER->SetMarkerColor(kRed+2);
  gr_cl00_BORDER->SetMarkerSize(1);
  gr_cl00_BORDER->SetLineWidth(2);
  gr_cl00_BORDER->SetMarkerStyle(20);
  //
  gr_cl01->SetMarkerColor(kBlue+2);
  gr_cl01->SetMarkerSize(1);
  gr_cl01->SetLineWidth(2);
  gr_cl01->SetMarkerStyle(2);
  //
  gr_cl01_BORDER->SetMarkerColor(kBlue+2);
  gr_cl01_BORDER->SetMarkerSize(1);
  gr_cl01_BORDER->SetLineWidth(2);
  gr_cl01_BORDER->SetMarkerStyle(20);
  //
  gr_cl02->SetMarkerColor(kMagenta+2);
  gr_cl02->SetMarkerSize(1);
  gr_cl02->SetLineWidth(2);
  gr_cl02->SetMarkerStyle(2);
  //
  gr_cl02_BORDER->SetMarkerColor(kMagenta+2);
  gr_cl02_BORDER->SetMarkerSize(1);
  gr_cl02_BORDER->SetLineWidth(2);
  gr_cl02_BORDER->SetMarkerStyle(20);
  //
  gr_cl03->SetMarkerColor(kGreen+2);
  gr_cl03->SetMarkerSize(1);
  gr_cl03->SetLineWidth(2);
  gr_cl03->SetMarkerStyle(2);
  //
  gr_cl03_BORDER->SetMarkerColor(kGreen+2);
  gr_cl03_BORDER->SetMarkerSize(1);
  gr_cl03_BORDER->SetLineWidth(2);
  gr_cl03_BORDER->SetMarkerStyle(20);
  //
  gr_NOISE->SetMarkerColor(kBlack);
  gr_NOISE->SetMarkerSize(1);
  gr_NOISE->SetLineWidth(2);
  gr_NOISE->SetMarkerStyle(20);
  //
  TMultiGraph *mg = new TMultiGraph();  
  mg->Add(gr);
  mg->Add(gr_cl00);
  mg->Add(gr_cl01);
  mg->Add(gr_cl02);
  mg->Add(gr_cl03);
  mg->Add(gr_cl00_BORDER);
  mg->Add(gr_cl01_BORDER);
  mg->Add(gr_cl02_BORDER);
  mg->Add(gr_cl03_BORDER);
  mg->Add(gr_NOISE);
  mg->Draw("ap");
  //
  TString file_out = fileN01;
  file_out += ".pdf";
  c1->SaveAs(file_out.Data());
  //
  /*
  mg->GetXaxis()->SetTitle("ValueX, Unit");
  mg->GetYaxis()->SetTitle("ValueY, Unit");
  
  TString legInfo;

  for(i = 0;i<nChannels;i++){
    legInfo = "ch ";legInfo += i;
    leg->AddEntry(gr_Arr[i], legInfo.Data(), "apl");
  }
  leg->Draw();
  
  */
  return 0;
}
