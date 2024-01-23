#include <stdio.h>
#include <iostream>
#include "src/dbscan.hh"
#include <fstream>
#include <string>

//root
#include <TGraph.h>
#include <TFile.h>
#include <TRandom3.h>

using namespace std;

void genData(vector<point>& points_v, TString filename);
void genData(vector<point>& points_v);
void plot_and_save(dbscan &ds, TString outrootFile, vector<Double_t> &k_dist_graph);

void genData_spiral(vector<point>& points_v, TString filename){
  ifstream cluster_file (filename.Data());
  Double_t x;
  Double_t y;
  Double_t z;
  if (cluster_file.is_open()){
    while(cluster_file>>x>>y>>z){
      point p;
      p.x = x;
      p.y = y;
      p.z = 0.0;
      p.point_id = (Int_t)points_v.size();
      points_v.push_back(p);
    }
    cluster_file.close();
  }
}

void genData_worms(vector<point>& points_v, TString filename){
  ifstream cluster_file(filename.Data());
  Double_t x;
  Double_t y;
  Int_t counter=0;
  if (cluster_file.is_open()){
    while(cluster_file>>x>>y){
      //if(x>3000 & y>3000){
      if(x>0 & y>0){
	point p;
	p.x = x;
	p.y = y;
	p.z = 0.0;
	p.point_id = (Int_t)points_v.size();
	points_v.push_back(p);      
      }
      counter++;
    }
    cluster_file.close();
  }
}

void genData(vector<point>& points_v, TString filename){
  ifstream cluster_file (filename.Data());
  Double_t x;
  Double_t y;
  Int_t counter=0;
  if (cluster_file.is_open()){
    while(cluster_file>>x>>y){
      //if(counter%5==0){
	point p;
	p.x = x;
	p.y = y;
	p.z = 0.0;
	p.point_id = (Int_t)points_v.size();
	points_v.push_back(p);
	//}
      counter++;
    }
    cluster_file.close();
  }
}

void genData(vector<point>& points_v){
  //
  unsigned int num_points_cluster1 = 100;
  unsigned int num_points_cluster2 = 100;
  unsigned int num_points_noise    = 20;
  //unsigned int num_points = num_points_cluster1 + num_points_cluster2 + num_points_noise;
  //
  Double_t x1 = -1.0;
  Double_t y1 = -1.0;
  Double_t rms1 = 0.1;
  //
  Double_t x2 = 1.0;
  Double_t y2 = 1.0;
  Double_t rms2 = 0.1;
  //
  TRandom3 *rnd = new TRandom3(123124);
  //
  for(unsigned int k = 0;k<num_points_cluster1;k++){
    point p;
    p.x = rnd->Gaus(x1,rms1);
    p.y = rnd->Gaus(y1,rms1);
    p.z = 0.0;
    p.point_id = (Int_t)points_v.size();
    points_v.push_back(p);
    //
    p.x = rnd->Gaus(x2,rms2);
    p.y = rnd->Gaus(y2,rms2);
    p.z = 0.0;
    p.point_id = (Int_t)points_v.size();
    points_v.push_back(p);
    //
    p.x = rnd->Gaus(x2-1,rms2);
    p.y = rnd->Gaus(y2+2,rms2);
    p.z = 0.0;
    p.point_id = (Int_t)points_v.size();
    points_v.push_back(p);
    //
    p.x = rnd->Gaus(x2-0.4,rms2);
    p.y = rnd->Gaus(y2+0.4,rms2);
    p.z = 0.0;
    p.point_id = (Int_t)points_v.size();
    points_v.push_back(p);
  }
  //
  for(unsigned int k = 0;k<num_points_noise;k++){
    point p;
    p.x = rnd->Uniform(-2.0,2.0);
    p.y = rnd->Uniform(-2.0,2.0);
    p.z = 0.0;
    p.point_id = (Int_t)points_v.size();
    points_v.push_back(p);
  }
}

void shuffle(vector<point>& points_v){
  vector<point> points_copy_v;
  for(unsigned int k = 0;k<points_v.size();k++){
    point p;
    p.x = points_v.at(k).x;
    p.y = points_v.at(k).y;
    p.z = points_v.at(k).z;
    p.point_id = points_v.at(k).point_id;
    points_copy_v.push_back(p);
  }
  TRandom3 *rnd = new TRandom3(123123);
  unsigned int k;
  unsigned int counter = 0;
  while(points_copy_v.size() != 0){
    k = rnd->Uniform( 0.0, ((Double_t)points_copy_v.size()-0.00001));
    if(counter<points_v.size()){
      points_v.at(counter).x = points_copy_v.at(k).x;
      points_v.at(counter).y = points_copy_v.at(k).y;
      points_v.at(counter).z = points_copy_v.at(k).z;
      points_v.at(counter).point_id = counter;
    }
    points_copy_v.erase(points_copy_v.begin()+k);
    counter++;
  }
}

int main(){    
  vector<point> points_v;
  vector<Double_t> k_dist_graph;
  dbscan ds;
  unsigned int minPts = 4;
  Double_t eps = 1.4;
  //
  //genData(points_v, "worms_2d.txt");
  //genData(points_v, "S_sets.s1");
  //genData_spiral(points_v, "spiral.txt");
  //genData_spiral(points_v, "Aggregation.txt");
  //genData_spiral(points_v, "Compound.txt");
  //genData(points_v);
  //
  //my
  //unsigned int minPts = 10;
  //Double_t eps = 0.08;
  //S_sets
  //unsigned int minPts = 15;
  //Double_t eps = 30000;
  //Aggregation
  //unsigned int minPts = 7;
  //Double_t eps = 1.7;
  //Compound
  //unsigned int minPts = 4;
  //Double_t eps = 1.4;
  //spiral
  //unsigned int minPts = 4;
  //Double_t eps = 3;
  //spiral
  //unsigned int minPts = 10000;
  //Double_t eps = 0.0001;
  //
  points_v.clear();
  k_dist_graph.clear();
  genData_spiral(points_v, "Compound.txt");
  //shuffle(points_v);
  minPts = 4;
  eps = 1.1;
  ds.run( minPts, eps, points_v);
  ds.print_points_info();
  ds.get_cluster_stats();
  ds.print_cluster_stats();
  k_dist_graph = ds.build_k_dist_graph(4);
  plot_and_save( ds, "hist_dbscan_Compound.root", k_dist_graph);
  ds.clear();
  //
  points_v.clear();
  k_dist_graph.clear();
  genData_spiral(points_v, "Aggregation.txt");
  minPts = 7;
  eps = 1.7;
  ds.run( minPts, eps, points_v);
  //ds.print_points_info();
  ds.get_cluster_stats();
  ds.print_cluster_stats();
  k_dist_graph = ds.build_k_dist_graph(4);
  plot_and_save( ds, "hist_dbscan_Aggregation.root", k_dist_graph);
  ds.clear();  
  //
  points_v.clear();
  k_dist_graph.clear();
  genData_spiral(points_v, "spiral.txt");
  minPts = 4;
  eps = 3;
  ds.run( minPts, eps, points_v);
  //ds.print_points_info();
  ds.get_cluster_stats();
  ds.print_cluster_stats();
  k_dist_graph = ds.build_k_dist_graph(4);
  plot_and_save( ds, "hist_dbscan_spiral.root", k_dist_graph);
  ds.clear();  
  //
  points_v.clear();
  k_dist_graph.clear();
  genData(points_v, "S_sets.s1");
  minPts = 15;
  eps = 30000;
  ds.run( minPts, eps, points_v);
  //ds.print_points_info();
  ds.get_cluster_stats();
  ds.print_cluster_stats();
  //k_dist_graph = ds.build_k_dist_graph(4);
  plot_and_save( ds, "hist_dbscan_S_sets.root", k_dist_graph);
  ds.clear();  
  //
  //points_v.clear();
  //genData_worms(points_v, "worms_2d.txt");
  //minPts = 23;
  //eps = 17;
  //ds.set_points(points_v);
  //ds.run( minPts, eps, points_v);
  //ds.print_points_info();
  //ds.get_cluster_stats();
  //ds.print_cluster_stats();
  //plot_and_save(ds,"hist_dbscan_worms_2d.root");
  //ds.clear();
  //
  return 0;
}

void plot_and_save(dbscan &ds, TString outrootFile, vector<Double_t> &k_dist_graph){
  //
  TGraph *gr = new TGraph();
  gr->SetNameTitle("gr","gr");
  //
  TGraph *gr_cl00 = new TGraph();
  gr_cl00->SetNameTitle("gr_cl00","gr_cl00");
  TGraph *gr_cl01 = new TGraph();
  gr_cl01->SetNameTitle("gr_cl01","gr_cl01");
  TGraph *gr_cl02 = new TGraph();
  gr_cl02->SetNameTitle("gr_cl02","gr_cl02");
  TGraph *gr_cl03 = new TGraph();
  gr_cl03->SetNameTitle("gr_cl03","gr_cl03");
  //
  TGraph *gr_cl00_BORDER = new TGraph();
  gr_cl00_BORDER->SetNameTitle("gr_cl00_BORDER","gr_cl00_BORDER");
  TGraph *gr_cl01_BORDER = new TGraph();
  gr_cl01_BORDER->SetNameTitle("gr_cl01_BORDER","gr_cl01_BORDER");
  TGraph *gr_cl02_BORDER = new TGraph();
  gr_cl02_BORDER->SetNameTitle("gr_cl02_BORDER","gr_cl02_BORDER");
  TGraph *gr_cl03_BORDER = new TGraph();
  gr_cl03_BORDER->SetNameTitle("gr_cl03_BORDER","gr_cl03_BORDER");
  //
  TGraph *gr_k_dist_graph = new TGraph();
  gr_k_dist_graph->SetNameTitle("gr_k_dist_graph","gr_k_dist_graph");
  //
  TGraph *gr_NOISE = new TGraph();
  gr_NOISE->SetNameTitle("gr_NOISE","gr_NOISE");
  //
  for(unsigned int k = 0; k < ds.get_points_v().size(); k++)
    gr->SetPoint( gr->GetN(), ds.get_points_v().at(k).x, ds.get_points_v().at(k).y);
  for(unsigned int k = 0; k < k_dist_graph.size(); k++)
    gr_k_dist_graph->SetPoint( k, k, k_dist_graph.at(k_dist_graph.size()-1-k));
  if(ds.get_nclusters()>0){
    for(unsigned int k = 0; k < ds.get_points_v().size(); k++){
      if(ds.get_points_v().at(k).clusterID == 0){
	gr_cl00->SetPoint( gr_cl00->GetN(), ds.get_points_v().at(k).x, ds.get_points_v().at(k).y);
	if(ds.get_points_v().at(k).point_type == BORDER_POINT)
	  gr_cl00_BORDER->SetPoint( gr_cl00_BORDER->GetN(), ds.get_points_v().at(k).x, ds.get_points_v().at(k).y);
      }
    }
    for(unsigned int k = 0; k < ds.get_points_v().size(); k++){
      if(ds.get_points_v().at(k).clusterID == 1){
	gr_cl01->SetPoint( gr_cl01->GetN(), ds.get_points_v().at(k).x, ds.get_points_v().at(k).y);
	if(ds.get_points_v().at(k).point_type == BORDER_POINT)
	  gr_cl01_BORDER->SetPoint( gr_cl01_BORDER->GetN(), ds.get_points_v().at(k).x, ds.get_points_v().at(k).y);
      }
    }
    for(unsigned int k = 0; k < ds.get_points_v().size(); k++){
      if(ds.get_points_v().at(k).clusterID == 2){
	gr_cl02->SetPoint( gr_cl02->GetN(), ds.get_points_v().at(k).x, ds.get_points_v().at(k).y);
	if(ds.get_points_v().at(k).point_type == BORDER_POINT)
	  gr_cl02_BORDER->SetPoint( gr_cl02_BORDER->GetN(), ds.get_points_v().at(k).x, ds.get_points_v().at(k).y);
      }
    }
    for(unsigned int k = 0; k < ds.get_points_v().size(); k++){
      if(ds.get_points_v().at(k).clusterID == 3){
	gr_cl03->SetPoint( gr_cl03->GetN(), ds.get_points_v().at(k).x, ds.get_points_v().at(k).y);
	if(ds.get_points_v().at(k).point_type == BORDER_POINT)
	  gr_cl03_BORDER->SetPoint( gr_cl03_BORDER->GetN(), ds.get_points_v().at(k).x, ds.get_points_v().at(k).y);
      }
    }
    for(unsigned int k = 0; k < ds.get_points_v().size(); k++)
      if(ds.get_points_v().at(k).point_type == NOISE)
	gr_NOISE->SetPoint( gr_NOISE->GetN(), ds.get_points_v().at(k).x, ds.get_points_v().at(k).y);
  }
  //
  //
  TFile* rootFile = new TFile(outrootFile.Data(), "RECREATE", " Histograms", 1);
  rootFile->cd();
  if (rootFile->IsZombie()){
    cout<<"  ERROR ---> file "<<outrootFile<<" is zombi"<<endl;
    assert(0);
  }
  else
    cout<<"  Output Histos file ---> "<<outrootFile<<endl;
  //
  gr->Write();
  gr_cl00->Write();
  gr_cl01->Write();
  gr_cl02->Write();
  gr_cl03->Write();
  gr_cl00_BORDER->Write();
  gr_cl01_BORDER->Write();
  gr_cl02_BORDER->Write();
  gr_cl03_BORDER->Write();
  gr_NOISE->Write();
  gr_k_dist_graph->Write();
  //
  rootFile->Close();
  //

  //vector<int> rr;
  //rr.push_back(1);
  //rr.push_back(2);
  //rr.push_back(3);
  //while( rr.size() != 0){
  //cout<<rr.at(0)<<endl;
  //rr.erase(rr.begin());
  //}
  //for(unsigned int k = 0;k<rr.size();k++)
  //cout<<rr.at(k)<<endl;
  //rr.erase(rr.begin(),rr.end());
  //rr.erase(rr.begin());
  //rr.erase(rr.end()-1);
  //for(unsigned int k = 0;k<rr.size();k++)
  //cout<<rr.at(k)<<endl;
  //rr.erase(rr.begin());
  //rr.erase(rr.end()-1);
  //for(unsigned int k = 0;k<rr.size();k++)
  //cout<<rr.at(k)<<endl;

  //vector<point> pv;
  //point p;
  //p.x=1;
  //p.y=1;
  //pv.push_back(p);
  //p.x=2;
  //p.y=2;
  //pv.push_back(p);

  //dbscan::print_points_info_st(pv);
  //point p2 = pv.at(0);
  //p2.x=100;
  //p2.y=100;
  //dbscan::print_points_info_st(pv);
}  
