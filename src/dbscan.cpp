//my
#include "dbscan.hh"

//root
#include "TROOT.h"
#include "TMath.h"

//C, C++
#include <iostream>
#include <assert.h>
#include <iomanip>
#include <stdlib.h>

dbscan::dbscan(){
  _minPts = 4;
  _eps = 0.1;
}

dbscan::~dbscan(){;}

Int_t dbscan::get_nclusters() const {
  Int_t nclus = -1;
  for(unsigned int k = 0; k<_points_v.size(); k++){
    if(_points_v.at(k).clusterID > nclus)
      nclus = _points_v.at(k).clusterID;
  }
  return (nclus+1);
}

void dbscan::get_cluster_stats(){
  //
  unsigned int nclus = (unsigned int)get_nclusters();
  if(nclus>0){
    for(unsigned int k = 0; k<nclus; k++){
      cluster_info cl;
      _clusters_v.push_back(cl);
    }
    //
    for(unsigned int k = 0; k<_points_v.size(); k++){
      if(_points_v.at(k).clusterID>=0 && _points_v.at(k).clusterID<(Int_t)nclus){
	_clusters_v.at(_points_v.at(k).clusterID).clusterID = _points_v.at(k).clusterID;
	_clusters_v.at(_points_v.at(k).clusterID).number_of_points++;
	//
	_clusters_v.at(_points_v.at(k).clusterID).mean_x += _points_v.at(k).x;
	_clusters_v.at(_points_v.at(k).clusterID).mean_y += _points_v.at(k).y;
	_clusters_v.at(_points_v.at(k).clusterID).mean_time_ii += _points_v.at(k).time_ii;
	//
	if(_points_v.at(k).point_type == CORE_POINT)
	  _clusters_v.at(_points_v.at(k).clusterID).number_of_CORE_POINT++;
	if(_points_v.at(k).point_type == BORDER_POINT)
	  _clusters_v.at(_points_v.at(k).clusterID).number_of_BORDER_POINT++;
      }
    }
    for(unsigned int k = 0; k < nclus; k++){
      if(_clusters_v.at(k).number_of_points>0){
	_clusters_v.at(k).mean_x /= _clusters_v.at(k).number_of_points;
	_clusters_v.at(k).mean_y /= _clusters_v.at(k).number_of_points;
	_clusters_v.at(k).mean_time_ii /= _clusters_v.at(k).number_of_points;
      }
    }
  }
}

void dbscan::clear(){
  _points_v.clear();
  _clusters_v.clear();
}

Int_t dbscan::get_number_of_NOISE() const{
  Int_t n_NOISE = 0; 
  for(unsigned int k = 0; k<_points_v.size(); k++)
    if(_points_v.at(k).point_type == NOISE)
      n_NOISE++;
  return n_NOISE;
}

void dbscan::print_cluster_stats(vector<cluster_info> clusters_v) {
  cout<<"N_clusters : "<<clusters_v.size()<<endl;
  if(clusters_v.size()>0){
    clusters_v.at(0).print_cluster_info_header();
    for(unsigned int i = 0; i<clusters_v.size(); i++)
      clusters_v.at(i).print_cluster_info();
  }  
}

void dbscan::print_cluster_stats() const {
  cout<<"N_points   : "<<_points_v.size()<<endl
      <<"N_clusters : "<<_clusters_v.size()<<endl
      <<"N_NOISE    : "<<get_number_of_NOISE()<<endl;
  cout<<"_minPts    : "<<_minPts<<endl
      <<"_eps       : "<<_eps<<endl;
  if(_clusters_v.size()>0){
    _clusters_v.at(0).print_cluster_info_header();
    for(unsigned int i = 0; i<_clusters_v.size(); i++)
      _clusters_v.at(i).print_cluster_info();
  }  
}

Int_t dbscan::run(unsigned int minPts, Double_t eps, vector<point> point_v){
  _minPts = minPts;
  _eps = eps;
  _points_v.clear();
  _points_v = point_v;
  //
  int clusterID = 0;
  for(unsigned int i = 0; i<_points_v.size(); i++)
    if((_points_v.at(i).point_type == UNCLASSIFIED) && expandCluster(_points_v.at(i), clusterID))
      clusterID++;
  return 0;
}

void dbscan::set_points(vector<point> point_v){
  _points_v.clear();
  _points_v = point_v;
}

bool dbscan::expandCluster(point &p, int clusterID){
  vector<unsigned int> seeds_v = regionQuery(p);
  if((seeds_v.size()+1)<_minPts){
    p.point_type = NOISE;
    return false;
  }
  else{
    p.point_type = CORE_POINT;
    p.clusterID  = clusterID;
    //
    while( seeds_v.size() != 0){
      point &seeds_p = _points_v.at(seeds_v.at(0));
      vector<unsigned int> neighbours_v = regionQuery(seeds_p);
      if((neighbours_v.size()+1) >= _minPts){
	for(unsigned int k = 0; k < neighbours_v.size(); k++){
	  point &neighbour_p = _points_v.at(neighbours_v.at(k));	  
	  if(neighbour_p.point_type == UNCLASSIFIED ||
	     neighbour_p.point_type == NOISE){
	    if(neighbour_p.point_type == UNCLASSIFIED)
	      seeds_v.push_back(neighbour_p.point_id);
	    neighbour_p.point_type = CORE_POINT;
	    neighbour_p.clusterID  = clusterID;
	  }
	}
      }
      else{
	seeds_p.point_type = BORDER_POINT;
	seeds_p.clusterID  = clusterID;
      }	
      seeds_v.erase(seeds_v.begin());
    }
    return true;
  }
  return false;
}

vector<Double_t> dbscan::build_k_dist_graph(Int_t kk){
  vector<Double_t> k_dist_graph;
  for(unsigned int k = 0;k<_points_v.size();k++){
    point_kdis pkd;
    pkd.p = _points_v.at(k);
    for(unsigned int l = 0;l<_points_v.size();l++)
      pkd.dists.push_back(calculateDistance(_points_v.at(k), _points_v.at(l)));
    point_kdis::bubbleSort(pkd.dists);
    _points_kdis_v.push_back(pkd);    
  }
  for(unsigned int k = 0;k<_points_kdis_v.size();k++){
    if((unsigned int)kk<_points_kdis_v.at(k).dists.size())
      k_dist_graph.push_back(_points_kdis_v.at(k).dists.at(kk));
  }
  point_kdis::bubbleSort(k_dist_graph);
  return k_dist_graph;
}

vector<unsigned int> dbscan::regionQuery(point p){
  vector<unsigned int> neighbours_index;
  for(unsigned int k = 0;k<_points_v.size();k++){
    if(p.point_id != _points_v.at(k).point_id){
      if(calculateDistance( p, _points_v.at(k)) <= _eps){
	neighbours_index.push_back(k);
      }
    }
  }
  return neighbours_index;
}

inline Double_t dbscan::calculateDistance( const point& pointCore, const point& pointTarget){
  return TMath::Sqrt(pow(pointCore.x - pointTarget.x,2) +
		     pow(pointCore.y - pointTarget.y,2) +
		     pow(pointCore.z - pointTarget.z,2));
}

void dbscan::print_points_info_st(const vector<point> point_v){
  std::cout<<setw(12)<<"x"
	   <<setw(12)<<"y"
	   <<setw(12)<<"z"
	   <<setw(12)<<"time_ii"
	   <<setw(12)<<"pixel_id"
	   <<setw(12)<<"clusterID"
	   <<setw(12)<<"point_type"
	   <<setw(12)<<"point_id"
	   <<std::endl;
  for(unsigned int k = 0;k<point_v.size();k++)
    std::cout<<setw(12)<<point_v.at(k).x
	     <<setw(12)<<point_v.at(k).y
	     <<setw(12)<<point_v.at(k).z
      	     <<setw(12)<<point_v.at(k).time_ii
      	     <<setw(12)<<point_v.at(k).pixel_id
	     <<setw(12)<<point_v.at(k).clusterID
	     <<setw(12)<<point_v.at(k).point_type
	     <<setw(12)<<point_v.at(k).point_id
	     <<std::endl;
}

void dbscan::print_points_info(vector<point> point_v){
  _points_v = point_v;
  print_points_info();
}

const void dbscan::print_points_info() const {
  std::cout<<"_minPts      "<<_minPts<<endl
	   <<"_eps         "<<_eps<<endl
	   <<"CORE_POINT   "<<CORE_POINT<<endl
	   <<"BORDER_POINT "<<BORDER_POINT<<endl
	   <<"UNCLASSIFIED "<<UNCLASSIFIED<<endl
	   <<"NOISE        "<<NOISE<<endl;   
  std::cout<<" nclusters : "<<get_nclusters()<<std::endl;  
  std::cout<<setw(12)<<"x"
	   <<setw(12)<<"y"
	   <<setw(12)<<"z"
	   <<setw(12)<<"time_ii"
	   <<setw(12)<<"pixel_id"
	   <<setw(12)<<"clusterID"
	   <<setw(12)<<"point_type"
	   <<setw(12)<<"point_id"
	   <<std::endl;
  for(unsigned int k = 0;k<_points_v.size();k++)
    std::cout<<setw(12)<<_points_v.at(k).x
	     <<setw(12)<<_points_v.at(k).y
	     <<setw(12)<<_points_v.at(k).z
      	     <<setw(12)<<_points_v.at(k).time_ii
      	     <<setw(12)<<_points_v.at(k).pixel_id
	     <<setw(12)<<_points_v.at(k).clusterID
	     <<setw(12)<<_points_v.at(k).point_type
	     <<setw(12)<<_points_v.at(k).point_id
	     <<std::endl;
}
