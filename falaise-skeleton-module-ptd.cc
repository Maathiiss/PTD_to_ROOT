#include <bayeux/dpp/chain_module.h>
#include <falaise/snemo/datamodels/particle_track_data.h>
#include <falaise/snemo/datamodels/tracker_trajectory_solution.h>
#include <falaise/snemo/datamodels/tracker_trajectory_data.h>
#include <falaise/snemo/datamodels/tracker_clustering_data.h>
#include <falaise/snemo/datamodels/tracker_cluster.h>
#include <falaise/snemo/datamodels/geomid_utils.h>  
#include <bayeux/mctools/base_step_hit.h>
#include <bayeux/mctools/simulated_data.h>
#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <string>
#include "TLatex.h"
#include "TVector3.h"
#include "TSystem.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <cmath>

using namespace std;


class falaise_skeleton_module_ptd : public dpp::chain_module
{

public:
  // Constructor
  falaise_skeleton_module_ptd();

  // Destructor
  virtual ~falaise_skeleton_module_ptd();

  // Initialisation function
  virtual void initialize (const datatools::properties &,
                           datatools::service_manager &,
			   dpp::module_handle_dict_type &);

  double compute_ellipse(double y_vertex, double z_vertex, int& source_num) const;
  // Event processing function
  dpp::chain_module::process_status process (datatools::things & event);
  
private:
  int event_number;
  int  ptd_event_counter;
  int cellules_non_associated, cellules_SD_non_associated, number_of_kinks;
  bool ptd_details, has_an_electron_and_positron, has_SD_electron_and_positron, has_SD_two_electrons, opposite_side_e_gamma, same_side_elec, same_cluster_elec, has_a_same, ttd_details, sd_calo_details, sd_tracker_details, hit_the_same_calo_hit;

  double diff_time_elec, angle_3D_between_ep_em, delta_y_elec, delta_z_elec, energy_elec_sum, corrected_energy_elec_sum, time_of_flight_gamma, internal_theoretical_time_diff, external_theoretical_time_diff,t1_th, t2_th, start_run_time, end_run_time, delta_r_calo, kink_angle, one_kink_x, one_kink_y, one_kink_z, total_volume;
  int nb_gamma, nb_elec_ptd_per_event,nb_elec_SD_per_event, nb_wire_hit;
  vector<string>  type_elec, g4_process, material, vertex_type, g4_material;
  vector<int> gamma_type, num_om, num_om_elec, track_number, num_gg;
  vector<double> time_gamma, time_gamma_before, time_gamma_after, angle_SD, volume, total_step_length;
  vector<double> energy_gamma, energy_elec, corrected_energy_elec, time_elec,track_lenght, energy_gamma_after, vertex_SD_x, vertex_SD_y, vertex_SD_z, ellipse_source;
  vector<double> vertex_3D_start_x, vertex_3D_start_y, vertex_3D_start_z,vertex_3D_end_x, vertex_3D_end_y, vertex_3D_end_z, vertex_gamma, kink_x, kink_y, kink_z;
  vector<int> side_elec, cluster_elec_num, calo_num;
  TFile *save_file;
  TTree *tree;
  std::vector<double> y_source_pos = {
    -2087.5, -2087.5, -2087.5, -2087.5, -2087.5, -2087.5, -2087.5,
    -1252.5, -1252.5, -1252.5, -1252.5, -1252.5, -1252.5, -1252.5,
    -417.5,  -417.5,  -417.5,  -417.5,  -417.5,  -417.5,  -417.5,
    417.5,   417.5,   417.5,   417.5,   417.5,   417.5,   417.5,
    1252.5,  1252.5,  1252.5,  1252.5,  1252.5,  1252.5,  1252.5,
    2087.5,  2087.5,  2087.5,  2087.5,  2087.5,  2087.5,  2087.5
  };

  std::vector<double> z_source_pos = {
    1317.5,   882.5,   447.5,     2.5,  -442.5,  -882.5, -1317.5,
    1317.5,   887.5,   447.5,     2.5,  -442.5,  -887.5, -1317.5,
    1312.5,   882.5,   447.5,     2.5,  -442.5,  -882.5, -1317.5,
    1317.5,   887.5,   447.5,     2.5,  -437.5,  -882.5, -1317.5,
    1317.5,   882.5,   442.5,     2.5,  -442.5,  -882.5, -1317.5,
    1317.5,   882.5,   442.5,     2.5,  -442.5,  -882.5, -1317.5
  };

  std::vector<double> source_pos_num = {
    6,  5,  4,  3,   2,  1,  0,
    13, 12, 11, 10,  9,  8,  7,
    20, 19, 18, 17, 16, 15, 14,
    27, 26, 25, 24, 23, 22, 21,
    34, 33, 32, 31, 30, 29, 28,
    41, 40, 39, 38, 37, 36, 35,
  };

  // Macro to register the module
  DPP_MODULE_REGISTRATION_INTERFACE(falaise_skeleton_module_ptd);

};

////////////////////////////////////////////////////////////////////

// Macro to add the module in the global register of data processing modules:
// The module defined by this class 'falaise_skeleton_module_ptd' will be registered
// with the label ID 'FalaiseSkeletonModule_PTD' (to use in pipeline configuration file)
DPP_MODULE_REGISTRATION_IMPLEMENT(falaise_skeleton_module_ptd, "FalaiseSkeletonModule_PTD")


falaise_skeleton_module_ptd::falaise_skeleton_module_ptd()
{
  save_file = new TFile("extracted_data.root", "RECREATE");
  tree = new TTree("Event", "Event information");  
  tree->Branch("event_number", &event_number);
  tree->Branch("has_SD_electron_and_positron", &has_SD_electron_and_positron);
  tree->Branch("has_SD_two_electrons", &has_SD_two_electrons);
  tree->Branch("nb_elec_SD_per_event", &nb_elec_SD_per_event);
  tree->Branch("cellules_SD_non_associated", &cellules_SD_non_associated);
  tree->Branch("angle_SD", &angle_SD);
  tree->Branch("total_step_length", &total_step_length);
  tree->Branch("volume", &volume);
  tree->Branch("total_volume", &total_volume);
  tree->Branch("g4_process", &g4_process);
  tree->Branch("material", &material);
  tree->Branch("g4_material", &g4_material);
  tree->Branch("vertex_type", &vertex_type);
  tree->Branch("track_number", &track_number);
  tree->Branch("nb_wire_hit", &nb_wire_hit);
  tree->Branch("vertex_SD_x", &vertex_SD_x);
  tree->Branch("vertex_SD_y", &vertex_SD_y);
  tree->Branch("vertex_SD_z", &vertex_SD_z);
  tree->Branch("one_kink_x", &one_kink_x);
  tree->Branch("one_kink_y", &one_kink_y);
  tree->Branch("one_kink_z", &one_kink_z);
  tree->Branch("nb_gamma", &nb_gamma);
  tree->Branch("gamma_type", &gamma_type);
  tree->Branch("time_gamma", &time_gamma);
  tree->Branch("energy_gamma", &energy_gamma);
  tree->Branch("time_gamma_before", &time_gamma_before);
  tree->Branch("time_gamma_after", &time_gamma_after);
  tree->Branch("energy_gamma_after", &energy_gamma_after);   
  tree->Branch("vertex_3D_start_x", &vertex_3D_start_x);
  tree->Branch("vertex_3D_start_y", &vertex_3D_start_y);
  tree->Branch("vertex_3D_start_z", &vertex_3D_start_z);
  tree->Branch("vertex_3D_end_x", &vertex_3D_end_x);
  tree->Branch("vertex_3D_end_y", &vertex_3D_end_y);
  tree->Branch("vertex_3D_end_z", &vertex_3D_end_z);
  tree->Branch("delta_r_calo", &delta_r_calo);
  tree->Branch("nb_elec_ptd_per_event", &nb_elec_ptd_per_event);
  tree->Branch("energy_elec", &energy_elec);
  tree->Branch("num_om", &num_om);
  tree->Branch("num_gg", &num_gg);
  tree->Branch("num_om_elec", &num_om_elec);
  tree->Branch("corrected_energy_elec", &corrected_energy_elec);
  tree->Branch("time_elec", &time_elec);
  tree->Branch("side_elec", &side_elec);
  tree->Branch("type_elec", &type_elec);
  tree->Branch("track_lenght", &track_lenght);
  tree->Branch("same_side_elec", &same_side_elec);
  tree->Branch("same_cluster_elec",&same_cluster_elec);
  tree->Branch("cluster_elec_num",&cluster_elec_num);
  tree->Branch("diff_time_elec", &diff_time_elec);
  tree->Branch("energy_elec_sum", &energy_elec_sum);
  tree->Branch("corrected_energy_elec_sum", &corrected_energy_elec_sum);
  tree->Branch("delta_y_elec", &delta_y_elec);
  tree->Branch("delta_z_elec", &delta_z_elec);
  tree->Branch("angle_3D_between_ep_em", &angle_3D_between_ep_em);
  tree->Branch("has_a_same", &has_a_same);
  tree->Branch("hit_the_same_calo_hit", &hit_the_same_calo_hit);
  tree->Branch("has_an_electron_and_positron", &has_an_electron_and_positron);
  tree->Branch("opposite_side_e_gamma",&opposite_side_e_gamma);
  tree->Branch("cellules_non_associated", &cellules_non_associated);
  tree->Branch("time_of_flight_gamma", &time_of_flight_gamma);
  tree->Branch("internal_theoretical_time_diff", &internal_theoretical_time_diff);
  tree->Branch("external_theoretical_time_diff", &external_theoretical_time_diff);
  tree->Branch("t1_th", &t1_th);
  tree->Branch("t2_th", &t2_th);
  tree->Branch("number_of_kinks", &number_of_kinks);
  tree->Branch("kink_x",&kink_x);
  tree->Branch("kink_y",&kink_y);
  tree->Branch("kink_z",&kink_z);
  tree->Branch("kink_angle",&kink_angle);
  tree->Branch("ellipse_source", &ellipse_source);

    
}


falaise_skeleton_module_ptd::~falaise_skeleton_module_ptd()
{
  save_file->cd();
  tree->Write();
  save_file->Close();
  delete save_file;
  std::cout << "falaise_skeleton_module_ptd::~falaise_skeleton_module_ptd() called" << std::endl;
}


void falaise_skeleton_module_ptd::initialize (const datatools::properties & module_properties, datatools::service_manager &, dpp::module_handle_dict_type &)
{
  std::cout << "falaise_skeleton_module_ptd::initialize() called" << std::endl;
  event_number=0;

  if ( module_properties.has_key("ptd_details"))
    ptd_details = module_properties.fetch_boolean("ptd_details");
  else ptd_details = false;

  if ( module_properties.has_key("ttd_details"))
    ttd_details = module_properties.fetch_boolean("ttd_details");
  else ttd_details = false;

  if ( module_properties.has_key("calo_details"))
    sd_calo_details = module_properties.fetch_boolean("calo_details");
  else sd_calo_details = false;

  if ( module_properties.has_key("tracker_details"))
    sd_tracker_details = module_properties.fetch_boolean("tracker_details");
  else sd_tracker_details = false;

  this->_set_initialized(true);
}


double compute_angle(const double* vertex_start_p1, const double* vertex_end_p1,
		     const double* vertex_start_p2, const double* vertex_end_p2)
{// this function compute the angle between two g4 steps
  double v1[3] = {
    vertex_end_p1[0] - vertex_start_p1[0], 
    vertex_end_p1[1] - vertex_start_p1[1],  
    vertex_end_p1[2] - vertex_start_p1[2] 
  };
  double v2[3] = {
    vertex_end_p2[0] - vertex_start_p2[0], 
    vertex_end_p2[1] - vertex_start_p2[1], 
    vertex_end_p2[2] - vertex_start_p2[2]   
  };
  double dot_product = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
  double norm_v1 = std::sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
  double norm_v2 = std::sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
  double angle = std::acos(dot_product / (norm_v1 * norm_v2)) * (180.0 / M_PI);
  return angle;
}


double compute_step_length(const double* vertex_start_p1, const double* vertex_end_p1)
{//this function compute the length of a step
  return std::sqrt(
		   std::pow(vertex_end_p1[0] - vertex_start_p1[0], 2) +
		   std::pow(vertex_end_p1[1] - vertex_start_p1[1], 2) +
		   std::pow(vertex_end_p1[2] - vertex_start_p1[2], 2)
		   );
}


double compute_volume(const double* start_track, const double* end_track,
		      const double* current_track)
{//this function compute the volumecreated by a step with the general track
  double v1[3] = {//general direction of the track
    end_track[0] - start_track[0], 
    end_track[1] - start_track[1],
    end_track[2] - start_track[2]
  };
  double v2[3] = {//current direction with the step
    current_track[0] - start_track[0],
    current_track[1] - start_track[1],
    current_track[2] - start_track[2] 
  };
  double v3[3] = {
    end_track[0] - current_track[0], 
    end_track[1] - current_track[1],
    end_track[2] - current_track[2]
  };
  double cross[3] = {
    v2[1] * v3[2] - v2[2] * v3[1],
    v2[2] * v3[0] - v2[0] * v3[2],
    v2[0] * v3[1] - v2[1] * v3[0]
  }; //cross product
  double combination_product = v1[0] * cross[0] + v1[1] * cross[1] + v1[2] * cross[2];
  double volume = combination_product/6.0; //the volume define by a point with a line is 1/6 of the combination product
  // if(volume==0){
  // }
  
  return volume;
}


double falaise_skeleton_module_ptd::compute_ellipse(double y_vertex, double z_vertex, int& source_num) const
{
  double min_dist2 = std::numeric_limits<double>::max();
  for (size_t i = 0; i < y_source_pos.size(); ++i) {
    double dy = y_source_pos[i] - y_vertex;
    double dz = z_source_pos[i] - z_vertex;
    double dist2 = (dy * dy) / (25.0 * 25.0) + (dz * dz) / (30.0 * 30.0);
    if (dist2 < min_dist2) {
      min_dist2 = dist2;
      source_num = source_pos_num[i];
    }
  }
  return min_dist2;
}



dpp::chain_module::process_status falaise_skeleton_module_ptd::process (datatools::things & event)
{
  // Skip processing if PTD bank is not present
  if (!event.has("PTD"))
    {
      std::cout << "======== no PTD bank in event " << ptd_event_counter++ << " ========" << std::endl;
      return dpp::base_module::PROCESS_SUCCESS;
      event_number++;
    }

  // Retrieve the PTD bank  
  const snemo::datamodel::particle_track_data & PTD = event.get<snemo::datamodel::particle_track_data>("PTD");

  nb_gamma=0;
  nb_wire_hit=0;
  vertex_gamma.clear();
  nb_elec_ptd_per_event=0;
  nb_elec_SD_per_event=0;
  cellules_non_associated=0;
  cellules_SD_non_associated=0;
  gamma_type.clear();
  time_gamma.clear();
  time_gamma_before.clear();
  time_gamma_after.clear();
  energy_gamma_after.clear();
  type_elec.clear();
  side_elec.clear();
  cluster_elec_num.clear();
  energy_elec.clear();
  calo_num.clear();
  num_om.clear();
  num_gg.clear();
  num_om_elec.clear();
  corrected_energy_elec.clear();
  time_elec.clear();
  track_lenght.clear();
  energy_gamma.clear();
  vertex_3D_start_x.clear();
  vertex_3D_start_y.clear();
  vertex_3D_start_z.clear();
  vertex_3D_end_x.clear();
  vertex_3D_end_y.clear();
  vertex_3D_end_z.clear();
  vertex_SD_x.clear();
  vertex_SD_y.clear();
  vertex_SD_z.clear();
  kink_x.clear();
  kink_y.clear();
  kink_z.clear();
  time_of_flight_gamma = 0.0;
  energy_elec_sum=0;
  corrected_energy_elec_sum = 0;
  kink_angle=0.0;
  same_side_elec=0;
  same_cluster_elec = 0;
  diff_time_elec= 0;
  delta_y_elec=0;
  delta_z_elec=0;
  has_a_same=0;
  angle_3D_between_ep_em = 0.0;
  angle_SD.clear();
  total_step_length.clear();
  volume.clear();
  total_volume=0.0;
  g4_process.clear();
  material.clear();
  g4_material.clear();
  vertex_type.clear();
  track_number.clear();
  delta_r_calo=0.0;
  has_an_electron_and_positron = 0;
  has_SD_two_electrons=0;
  has_SD_electron_and_positron=0;
  opposite_side_e_gamma=0;
  internal_theoretical_time_diff = 0.0;
  external_theoretical_time_diff =0.0;
  t1_th=0.0;
  t2_th =0.0;
  start_run_time=0.0;
  end_run_time=0.0;
  one_kink_x=0.0;
  one_kink_y=0.0;
  one_kink_z=0.0;
  number_of_kinks=0;
  hit_the_same_calo_hit=false;
  ellipse_source.clear();
  //SD extraction
  if (event.has("SD")){
    const mctools::simulated_data & SD = event.get<mctools::simulated_data>("SD");
    const mctools::simulated_data::primary_event_type &pevent = SD.get_primary_event();
    const genbb::primary_event::particles_col_type &particles =pevent.get_particles();
    if(particles.size()==2){
      if(particles.front().get_particle_label()=="e-" && particles.back().get_particle_label()=="e-"){
	has_SD_two_electrons = true;
      }
      if(particles.front().get_particle_label()=="e+" && particles.back().get_particle_label()=="e-" || particles.front().get_particle_label()=="e-" && particles.back().get_particle_label()=="e+"){
	has_SD_electron_and_positron = true;
      }
    }
      
    vector<int> old_step_hit;
    vector<double> time_SD;
    if (SD.has_step_hits("__visu.tracks")) {
      for (UInt_t ihit = 0; ihit < SD.get_number_of_step_hits("__visu.tracks"); ihit++){
	auto stepHit = SD.get_step_hit("__visu.tracks", ihit);
	if(stepHit.get_particle_name() != "gamma" && stepHit.get_particle_name() != "alpha"){
	  old_step_hit.push_back(ihit);
	  time_SD.push_back(stepHit.get_time_start());
	}
      }
      std::vector<std::pair<double, int>> paired;
      for (size_t i = 0; i < old_step_hit.size(); i++) {
	paired.push_back({time_SD[i], old_step_hit[i]});
      }
      std::sort(paired.begin(), paired.end());

      std::vector<std::vector<std::array<double, 3>>> vertex_3D_start(30);
      std::vector<std::vector<std::array<double, 3>>> vertex_3D_end(30);
      std::vector<std::vector<double>> time(30);
      std::vector<std::vector<string>> process(30);
      std::vector<std::vector<string>> g4_SD_material(30);
      std::vector<std::vector<string>> SD_material(30);
      std::vector<std::vector<double>> step_length(30);
      int track_max=0;
      for (size_t i = 0; i < paired.size(); i++) {
	auto new_stepHit = SD.get_step_hit("__visu.tracks", paired[i].second);
	int currentTrackID = new_stepHit.get_track_id();
	if(currentTrackID>29){
	  break;
	}
	if(currentTrackID>track_max){
	  track_max=currentTrackID;
	}
	vertex_3D_start[currentTrackID].push_back({new_stepHit.get_position_start().getX(), new_stepHit.get_position_start().getY(), new_stepHit.get_position_start().getZ()});	
	vertex_3D_end[currentTrackID].push_back({new_stepHit.get_position_stop().getX(), new_stepHit.get_position_stop().getY(), new_stepHit.get_position_stop().getZ()});
	time[currentTrackID].push_back(new_stepHit.get_time_start());
	SD_material[currentTrackID].push_back(new_stepHit.get_material_name());
	g4_SD_material[currentTrackID].push_back(new_stepHit.get_g4_volume_name());
	step_length[currentTrackID].push_back(new_stepHit.get_step_length());
	  
	if(new_stepHit.has_creator_process_name()){
	  process[currentTrackID].push_back(new_stepHit.get_creator_process_name());
	}
	else{
	  process[currentTrackID].push_back("no_process");		    
	}
      }    

      set<string> allowed_materials_for_angle = {"basic::mylar", "", "tracking_gas", "bb_source_material.basic"};
      
      vector<int> nb_hits;
      std::vector<vector<string>> vertex;
      std::vector<int> track_id;
      int index_start, index_end; //to compute start and end of the tracks
      index_start = 0;
      index_end = 0;
      for (size_t i = 1; i < track_max+1; i++) {
	if(track_max>30){
	  cout<<"event number weird "<<endl;
	  break;
	}
	vertex.push_back({});
	nb_hits.push_back(vertex_3D_start[i].size());
	for(size_t j=0; j<vertex_3D_start[i].size(); j++){
	  if(j>0){	    
	    if ((abs(vertex_3D_start[i][j-1][0]) > 400) != (abs(vertex_3D_start[i][j][0]) > 400)){ //the track cross the geometrical condition x=400
	      vertex[i-1].push_back("calo");
	      index_end = j;
	    }
	    else if ((abs(vertex_3D_start[i][j-1][0]) > 40) != (abs(vertex_3D_start[i][j][0]) > 40)){ // the track cross the geometrical condition x=40
	      vertex[i-1].push_back("foil");
	      index_start=j;
	    }
	  }	
	}
      }
      
      for(int i=0; i<vertex.size(); i++){//loop on tracks
	bool already_registered=0;
	if(vertex[i].size()<2){
	  if(vertex[i].size()==1 && nb_hits[i]>30){ //track big enough to worried us
	    cellules_SD_non_associated++;
	  }
	  continue;
	}
	double start_kink_step1[3];
	double start_kink_step2[3];
	double end_kink_step1[3];
	double end_kink_step2[3];
	
	for(int j=0; j<vertex[i].size(); j++){//loop on vertices
	  track_number.push_back(i);
	  vertex_type.push_back(vertex[i][j]);
	  if(j>0){
	    if((vertex[i][j-1]=="calo" && vertex[i][j]=="foil") || (vertex[i][j-1]=="foil" && vertex[i][j]=="calo")){
	      nb_elec_SD_per_event++;
	      int index = j+1;
	      if(index!=(vertex[i].size()-1)){
		j++;
		vertex_type.push_back(vertex[i][j-1]);
	      }
	      else{//if you found calo-foil or foil-calo just 1 index before the end of the vector it means that there is one foil or calo that is the last point -> noise
		cellules_SD_non_associated++;
	      }
	    
	      //compute general edge of the track
	      double start_track[3] = {vertex_3D_start[i+1][index_start][0], vertex_3D_start[i+1][index_start][1], vertex_3D_start[i+1][index_start][2]};
	      double end_track[3] = {vertex_3D_start[i+1][index_end][0], vertex_3D_start[i+1][index_end][1], vertex_3D_start[i+1][index_end][2]};
	      
	      //we add everythings if we found 2 vertices in SD
	      if(already_registered==0){
		already_registered=1;
		for(size_t k=0; k<vertex_3D_start[i+1].size(); k++){
		  if(k>0){
		    double start_p1[3] = {vertex_3D_start[i+1][k-1][0], vertex_3D_start[i+1][k-1][1], vertex_3D_start[i+1][k-1][2]};
		    double end_p1[3] = {vertex_3D_end[i+1][k-1][0], vertex_3D_end[i+1][k-1][1], vertex_3D_end[i+1][k-1][2]};
		    double start_p2[3] = {vertex_3D_start[i+1][k][0], vertex_3D_start[i+1][k][1], vertex_3D_start[i+1][k][2]};
		    double end_p2[3] = {vertex_3D_end[i+1][k][0], vertex_3D_end[i+1][k][1], vertex_3D_end[i+1][k][2]};
		    if(start_p2[0] == end_p2[0] && start_p2[1] == end_p2[1] && start_p2[2] == end_p2[2]){
		      //sometimes there are 2 steps with exactly same coordinates (on material change)
		      if (k + 1 < vertex_3D_start[i+1].size()){
			if(SD_material[i+1][k]=="tracking_gas" && (SD_material[i+1][k+1]=="anode"||SD_material[i+1][k+1]=="cathode" || SD_material[i+1][k+1]=="anode_structure")){
			  if(SD_material[i+1][k+1]!="anode_structure"){
			    nb_wire_hit++;
			  }
			  for (int z = 0; z < 3; z++) {						    
			    start_kink_step1[z] = start_p1[z];
			    start_kink_step2[z] = vertex_3D_end[i+1][k+1][z];
			  }
			}
			for (int z = 0; z < 3; z++) {
			  start_p2[z] = vertex_3D_start[i+1][k+1][z];
			  end_p2[z] = vertex_3D_end[i+1][k+1][z];
			}
		      
			k++;		      
		      }	  
		    }	  
		    g4_process.push_back(process[i+1][k]);
		    material.push_back(SD_material[i+1][k]);
		    g4_material.push_back(g4_SD_material[i+1][k]);
		    vertex_SD_x.push_back(vertex_3D_start[i+1][k][0]);
		    //we take the start of the second particle to compare
		    vertex_SD_y.push_back(vertex_3D_start[i+1][k][1]);
		    vertex_SD_z.push_back(vertex_3D_start[i+1][k][2]);
		    if(k>index_start && k<index_end){
		      //we compute the volume only if the vertices are after the end and start point of the track
		      double current_volume = compute_volume(start_track, end_track, end_p2);
		      //end_p2 is the current track index		      
		      volume.push_back(current_volume);
		      total_volume += current_volume;
		    }
		    else{
		      volume.push_back(100);
		    }
		    double first_step_length = compute_step_length(start_p1,end_p1);
		    double second_step_length = compute_step_length(start_p2,end_p2);
		    total_step_length.push_back(first_step_length+second_step_length);
		    
		    string material1 = SD_material[i+1][k];
		    string material2 = SD_material[i+1][k-1];
		    string anode_material = "";
		    if(k<vertex_3D_start[i+1].size()-1){
		      anode_material = SD_material[i+1][k+1];
		    }

		    if(material1 != material2 && allowed_materials_for_angle.count(material1) > 0 && allowed_materials_for_angle.count(material2) > 0){
		      angle_SD.push_back(100);
		    }
		    else{
		      if(material1==anode_material && (anode_material=="anode"||anode_material=="cathode" || anode_material=="anode_structure")){
			//inside cathodic and anodic wires -> angle not computed
			angle_SD.push_back(100);
		      }
		      else if(material1=="tracking_gas" && (anode_material=="anode"||anode_material=="cathode" || anode_material=="anode_structure")){
			//starting point of the kink in wire
			if(anode_material!="anode_structure"){
			  nb_wire_hit++;
			}
			for (int j = 0; j < 3; j++) {
			  start_kink_step1[j] = start_p1[j]; //save value
			  start_kink_step2[j] = end_p1[j];
			}
			angle_SD.push_back(100);	
		      }
		      
		      else if((material1=="anode"||material1=="cathode"||material1=="anode_structure") && anode_material=="tracking_gas"){
			//ending point of the kink in wire
			for (int j = 0; j < 3; j++) {
			  end_kink_step1[j] = start_p2[j]; //save value
			  end_kink_step2[j] = end_p2[j];
			}
			angle_SD.push_back(compute_angle(start_kink_step1, start_kink_step2, end_kink_step1,end_kink_step2));
			//we compute the angle between the first and the last step in the cathode/anode
		      }
		      else{
			if(k==vertex_3D_start[i+1].size()-1 && material1=="anode_structure"){  
			  //if anodic structure = end of vector (particle stop in anode)
			  angle_SD.push_back(100);
			  
			}
			else{
			  angle_SD.push_back(compute_angle(start_p1, end_p1, start_p2,end_p2));
			}
		      }
		    }
		  }
		}
	      }
	    }
	    else{//if foil foil or calo calo -> noise
	      //cout<<"noise "<<event_number<<endl;
	      cellules_SD_non_associated++;
	    }
	  }	  
	}
      }
    }
  }
  
  

  
    
  // PTD extraction
  if (ptd_details)    
    {
      const snemo::datamodel::tracker_trajectory_data & TTD = event.get<snemo::datamodel::tracker_trajectory_data>("TTD");
      //if(TTD.get_solution_id()==0){
      const snemo::datamodel::tracker_trajectory_solution & ttd_solution = TTD.get_default_solution();
      //const snemo::datamodel::tracker_trajectory_solution & ttd_solution = TTD.get_solution(0);
      //cout<<"solution number "<<ttd_solution.get_solution_id()<<endl;
      std::vector<std::vector<int>> same_clusters;
      std::vector<std::vector<double>> kinks_vec;
      int nb_kinks_save=0;
      int count_equal_clusters = 0;
      int nb_clusters=0;
      //check unfitted cluster
      
      if (ttd_details)
	{
	  //const snemo::datamodel::tracker_trajectory_solution & ttd_solution = TTD.get_default_solution();
	  if (ttd_solution.has_unfitted_clusters()){
	    //cout<<"event "<<event_number<<endl;
	    cellules_non_associated++;	
	  }

	}       
      
      
      //check gammas
      vector<int>  particle_type;
      if (PTD.hasIsolatedCalorimeters()) {//gamma start
	const snemo::datamodel::CalorimeterHitHdlCollection &cc_collection = PTD.isolatedCalorimeters();
	for (const auto &it_hit : cc_collection) {
	  //possible de tout récuperer mais pas de vertex -> pas de gammatype  
	  energy_gamma.push_back(it_hit->get_energy());
	  gamma_type.push_back(it_hit->get_geom_id().get(0));
	  vertex_gamma.push_back(it_hit->get_geom_id().get(1));
	  time_gamma.push_back(it_hit->get_time());
	  nb_gamma++;
	  num_om.push_back(snemo::datamodel::om_num(it_hit->get_geom_id()));
	}
      }
     
      //electron start
      //vector<vector<double>> vertex_gamma_tot;
    
      vector<int> cluster_id;
      for (const datatools::handle<snemo::datamodel::particle_track> & particle : PTD.particles())
	{
	  const auto& trajectory_pattern = particle->get_trajectory_handle()->get_pattern();
	  int nb_kinks = trajectory_pattern.number_of_kinks();
	  number_of_kinks = nb_kinks;
	  if(nb_kinks>0){
	    nb_kinks_save = nb_kinks;
	    for (unsigned int i = 0; i < nb_kinks; ++i) {
	      kink_x.push_back(trajectory_pattern.get_kink(i).getX());
	      kink_y.push_back(trajectory_pattern.get_kink(i).getY());
	      kink_z.push_back(trajectory_pattern.get_kink(i).getZ());
	    }
	  }
	
	  
	  //particle->get_trajectory_handle()->get_pattern().get_first();
	    

	  
	  if(particle->get_charge() == snemo::datamodel::particle_track::CHARGE_UNDEFINED || particle->get_charge() == snemo::datamodel::particle_track::CHARGE_POSITIVE || particle->get_charge() == snemo::datamodel::particle_track::CHARGE_NEGATIVE){
	    if(particle->get_charge() == snemo::datamodel::particle_track::CHARGE_POSITIVE){
	      particle_type.push_back(1);
	    }
	    if(particle->get_charge() == snemo::datamodel::particle_track::CHARGE_NEGATIVE){
	      particle_type.push_back(-1);
	    }
	    bool vertex_close_to_the_source = false;
	    bool vertex_associated_to_a_calo = false;
	    double x_calo, y_calo, z_calo, x_foil, y_foil, z_foil;	  
	      
	    for(const datatools::handle<snemo::datamodel::vertex> & vertex : particle->get_vertices()){
	      if(/*vertex->is_on_reference_source_plane() ||*/ vertex->is_on_source_foil()){
		x_foil = vertex->get_spot().get_position().getX();
		y_foil = vertex->get_spot().get_position().getY();
		z_foil = vertex->get_spot().get_position().getZ();
		vertex_close_to_the_source=1;
	      }
	      if(vertex->is_on_main_calorimeter() /*|| vertex->is_on_x_calorimeter()*/ && particle->get_associated_calorimeter_hits().size()==1) //we forced MW only analysis
		{
		  // if(vertex->is_on_main_calorimeter())
		  //   {
		  //     type_elec.push_back("MW");
		  //   }
		  // else
		  //   {
		  //     type_elec.push_back("XW");
		  //   }
		  x_calo = vertex->get_spot().get_position().getX();
		  y_calo = vertex->get_spot().get_position().getY();
		  z_calo = vertex->get_spot().get_position().getZ();		
		  vertex_associated_to_a_calo = true;
		}
	    }
	    int source_num;
	    
	    if(vertex_close_to_the_source && vertex_associated_to_a_calo){
	      ellipse_source.push_back(compute_ellipse(y_foil,z_foil,source_num));
	      
	      const auto& trajectory_pattern_elec = particle->get_trajectory_handle()->get_pattern();
	      int nb_kinks_elec = trajectory_pattern_elec.number_of_kinks();
	      if(nb_kinks_elec==1){
		double x2 = trajectory_pattern_elec.get_kink(0).getX();
		double y2 = trajectory_pattern_elec.get_kink(0).getY();
		double z2 = trajectory_pattern_elec.get_kink(0).getZ();
		one_kink_x = x2;
		one_kink_y= y2;
		one_kink_z=z2;
		double dot_product = (x2 - x_foil) * (x_calo - x2) + (y2 - y_foil) * (y_calo - y2) + (z2 - z_foil) * (z_calo - z2);
		double norm_u = sqrt((x2 - x_foil) * (x2 - x_foil) + (y2 - y_foil) * (y2 - y_foil) + (z2 - z_foil) * (z2 - z_foil));
		double norm_v = sqrt((x_calo - x2) * (x_calo - x2) + (y_calo - y2) * (y_calo - y2) + (z_calo - z2) * (z_calo - z2));
		kink_angle=(180.0/M_PI)*acos(dot_product / (norm_u * norm_v));
	      }
	      else if(nb_kinks_elec>1){
		kink_angle=100;
	      }
	      
	      vertex_3D_start_x.push_back(x_foil);
	      vertex_3D_start_y.push_back(y_foil);
	      vertex_3D_start_z.push_back(z_foil);
	      vertex_3D_end_x.push_back(x_calo);
	      vertex_3D_end_y.push_back(y_calo);
	      vertex_3D_end_z.push_back(z_calo);
	      cluster_id.push_back(particle->get_trajectory().get_cluster().get_hit_id());
	      //cout<<endl;
	      //cout<<"event "<<event_number<<endl;
	      const auto& calorimeter_hits = particle->get_associated_calorimeter_hits();
	      if(calorimeter_hits[0]->has_auxiliaries()){//you can take the index 0 because I verified that calorimeter_hits.size()==1
		//cout<<"has auxiliaries "<<calorimeter_hits[0]->has_auxiliaries()<<endl;
		const auto& aux = calorimeter_hits[0]->get_auxiliaries();
		for (const auto& key : aux.keys()) {
		  const auto& value = aux.get(key);
		  // std::cout << "Key: " << key<<endl;
		  // cout<<"type = "<<value.get_type()<<endl;
		  // cout<<"size = "<<value.size()<<endl;
		  if (key == "Ef_optical_loss") {
		    double real_value = 0.0;
		    value.get_value(real_value);
		    corrected_energy_elec.push_back(real_value);
		    //cout << "Value: " << real_value << endl;
		  }		    
		}
	      }
	      if(nb_kinks>0){//looking if this particle get kinks
		kinks_vec.push_back({trajectory_pattern.get_kink(0).getX(),trajectory_pattern.get_kink(0).getY(),trajectory_pattern.get_kink(0).getZ()});
	      }
	      else{
		kinks_vec.push_back({});
	      }
	      calo_num.push_back(snemo::datamodel::om_num(calorimeter_hits[0]->get_geom_id()));
	      energy_elec.push_back(calorimeter_hits[0]->get_energy());
	      time_elec.push_back(calorimeter_hits[0]->get_time());
	      side_elec.push_back(calorimeter_hits[0]->get_geom_id().get(1));
	      num_om_elec.push_back(snemo::datamodel::om_num(calorimeter_hits[0]->get_geom_id()));
	      num_om.push_back(snemo::datamodel::om_num(calorimeter_hits[0]->get_geom_id()));

	      const auto gg = particle->get_trajectory_handle()->get_cluster().hits();
	      for(const auto hits : gg){
		num_gg.push_back(snemo::datamodel::gg_num(hits->get_geom_id()));
	      }

		
	      if (particle->has_trajectory()) {
		const auto& trajectory = particle->get_trajectory();
		if (trajectory.has_cluster()) {
		  cluster_elec_num.push_back(trajectory.get_cluster().get_hit_id());
		}
	      }
	      //cluster_elec_num.push_back(particle.get_trajectory().get_cluster.get_hit_id());
	      //get from inheritance of the class geomtools::base_hit
	      //nb_elec_ptd_per_event++;
	    }
	    else{
	      cellules_non_associated++;	      
	    }
	  }//end e-	  	  
	  
	  else{
	    cellules_non_associated++;
	  }
	}//end loop on particle
      
      for (size_t i=0; i<cluster_id.size(); i++) {
	for (size_t j=i+1; j<cluster_id.size(); j++) {
          if(cluster_id[i]==cluster_id[j]){
            same_clusters.push_back({cluster_id[i],cluster_id[j]});
            count_equal_clusters++;
            has_a_same=1;
          }
        }
      }

      vector<int> indices_calo_final;
      vector<int> tot_same;
      for(int i=0; i< same_clusters.size(); i++){
        for(int j=0; j< same_clusters[i].size(); j++){
          tot_same.push_back(same_clusters[i][j]);
        }
      }

      for(int j=0; j<cluster_id.size(); j++){
        if(std::find(tot_same.begin(), tot_same.end(), cluster_id[j]) == tot_same.end()){
          indices_calo_final.push_back(j);
        }
      }


      



      if(count_equal_clusters!=0){
        int j_final=0;
        for(int i=0; i<same_clusters.size(); i++){// we do the operation for every same clusters    
          bool cluster_found=false;
          for(int j=0; j<cluster_id.size();j++){
            if(same_clusters[i][0]==cluster_id[j] || same_clusters[i][1]==cluster_id[j]){ //we suppose one cluster can give only 2 solutions                                                                  
              cluster_found=true;
              double min_y_final = 1000;
              double min_z_final = 1000;
              for(int k=0;k<vertex_3D_start_y.size();k++){
                if(k==j){continue;}
                double min_temp_y = abs(vertex_3D_start_y[j]-vertex_3D_start_y[k]);
                double min_temp_z = abs(vertex_3D_start_z[j]-vertex_3D_start_z[k]);
                if(min_y_final+min_z_final>min_temp_y+min_temp_z){
                  min_y_final=min_temp_y;
                  min_z_final = min_temp_z;
                  j_final = j;
                }
              }
            }
          }
          if(cluster_found==true){
            indices_calo_final.push_back(j_final);
          }
        }
      }
      nb_clusters = cluster_id.size();
      nb_elec_ptd_per_event=indices_calo_final.size();

      
      if(nb_elec_ptd_per_event==2){	
	if(num_om_elec[indices_calo_final.at(0)]==num_om_elec[indices_calo_final.at(1)]){
	  hit_the_same_calo_hit=true;
	  //cout<<"event number "<<event_number<<endl;
	}
      }
      if(indices_calo_final.size()==2 /*&& num_om[indices_calo_final.at(0)]!=num_om[indices_calo_final.at(1)]*/){//final search
	//cout<<"event number "<<event_number<<endl;
	if(nb_clusters-count_equal_clusters!=2){
	  cellules_non_associated++;
	}
	//cout<<"taille indice finale : "<<indices_calo_final.size()<<"entrie  = "<<event_number<<endl;
	if(cluster_elec_num[indices_calo_final.at(0)] == cluster_elec_num[indices_calo_final.at(1)]){
	  same_cluster_elec=1;
	}
	if(side_elec[indices_calo_final.at(0)]==side_elec[indices_calo_final.at(1)]){
	  same_side_elec=1;
	}
	if(particle_type.size()>0){	    //for magnetic field
	  if((particle_type[indices_calo_final.at(0)]==1 && particle_type[indices_calo_final.at(1)]==-1) || (particle_type[indices_calo_final.at(0)]==-1 && particle_type[indices_calo_final.at(1)]==1)){
	    has_an_electron_and_positron = 1;
	  }
	}
	delta_r_calo=sqrt(pow(vertex_3D_end_y[indices_calo_final.at(0)]-vertex_3D_end_y[indices_calo_final.at(1)],2)+pow(vertex_3D_end_z[indices_calo_final.at(0)]-vertex_3D_end_z[indices_calo_final.at(1)],2));
	energy_elec_sum=energy_elec[indices_calo_final.at(0)]+energy_elec[indices_calo_final.at(1)];
	if(corrected_energy_elec.size()>0){
	  corrected_energy_elec_sum = corrected_energy_elec[indices_calo_final.at(0)]+corrected_energy_elec[indices_calo_final.at(1)];
	}
		
	delta_y_elec = abs(vertex_3D_start_y[indices_calo_final.at(0)]-vertex_3D_start_y[indices_calo_final.at(1)]);
	delta_z_elec = abs(vertex_3D_start_z[indices_calo_final.at(0)]-vertex_3D_start_z[indices_calo_final.at(1)]);
	//starting point of the two particles
	double x1 = vertex_3D_start_x[indices_calo_final.at(0)];
	double y1 = vertex_3D_start_y[indices_calo_final.at(0)];
	double z1 = vertex_3D_start_z[indices_calo_final.at(0)];
	double x1_ = vertex_3D_start_x[indices_calo_final.at(1)];
	double y1_ = vertex_3D_start_y[indices_calo_final.at(1)];
	double z1_ = vertex_3D_start_z[indices_calo_final.at(1)];	
	//ending point of the two particles
	double x2 = vertex_3D_end_x[indices_calo_final.at(0)];
	double y2 = vertex_3D_end_y[indices_calo_final.at(0)];
	double z2 = vertex_3D_end_z[indices_calo_final.at(0)];
	double x2_ = vertex_3D_end_x[indices_calo_final.at(1)];
	double y2_ = vertex_3D_end_y[indices_calo_final.at(1)];
	double z2_ = vertex_3D_end_z[indices_calo_final.at(1)];	
	if(kinks_vec[indices_calo_final.at(0)].size()>0){
	  x2 = kinks_vec[indices_calo_final.at(0)][0];
	  y2 = kinks_vec[indices_calo_final.at(0)][1];
	  z2 = kinks_vec[indices_calo_final.at(0)][2];
	}
	if(kinks_vec[indices_calo_final.at(1)].size()>0){
	  x2 = kinks_vec[indices_calo_final.at(1)][0];
	  y2 = kinks_vec[indices_calo_final.at(1)][1];
	  z2 = kinks_vec[indices_calo_final.at(1)][2];
	}
	double start_p1[3] = {x1, y1, z1};
	double start_p2[3] = {x1_, y1_, z1_};
	double end_p1[3] = {x2, y2, z2};
	double end_p2[3] = {x2_, y2_, z2_};	
	angle_3D_between_ep_em = abs(compute_angle(start_p1, end_p1, start_p2, end_p2));
	diff_time_elec = time_elec[indices_calo_final.at(0)]-time_elec[indices_calo_final.at(1)];

	
	//electron crossing or pair creation ??
	double new_c = 299.792458;//change c unit to mm/ns
	double track_length_1 = sqrt(pow(x1-x2,2) + pow(y1-y2,2) + pow(z1-z2,2));
	double E_1 = energy_elec[indices_calo_final.at(0)];
	double beta_1 = sqrt(E_1*(E_1+2*0.511))/(E_1+0.511);	
	double t_1 = track_length_1/(beta_1*new_c);
	double track_length_2 = sqrt(pow(x1_-x2_,2) + pow(y1_-y2_,2) + pow(z1_-z2_,2));
	double E_2 = energy_elec[indices_calo_final.at(1)];
	double beta_2 = sqrt(E_2*(E_2+2*0.511))/(E_2+0.511);
	double t_2 = track_length_2/(beta_2*new_c);
	internal_theoretical_time_diff = (t_1-t_2) - diff_time_elec;
	external_theoretical_time_diff = abs(diff_time_elec) - (t_1+t_2) ;
	t1_th = t_1;
	t2_th = t_2;
	track_lenght.push_back(track_length_1);
	track_lenght.push_back(track_length_2);
	// if(internal_theoretical_time_diff<-2){
	// cout<<"entry "<<event_number<<endl;
	// cout<<"t1 th = "<<t_1<<" t2 th "<< t_2<< " diff time elec "<<diff_time_elec<<endl;
	// cout<<"t1 exp = "<<time_elec[indices_calo_final.at(0)]<<" t2 exp "<< time_elec[indices_calo_final.at(1)]<<endl;
 	
	// //Errorbar to compute and add to conclude
	// cout<<"internal " <<internal_theoretical_time_diff<<" external "<<external_theoretical_time_diff<<endl;
	// cout<<endl;
	// }
	

      
	for(int time_gamma_value=0; time_gamma_value<time_gamma.size();time_gamma_value++){
	  time_gamma_after.push_back(time_gamma.at(time_gamma_value)-*std::max_element(time_elec.begin(),time_elec.end()));
	  time_gamma_before.push_back(time_gamma.at(time_gamma_value)-*std::min_element(time_elec.begin(),time_elec.end()));
	  energy_gamma_after.push_back(energy_gamma.at(time_gamma_value));
	  //no problem with time and energy for double case beacuse it is the same
	  if(nb_gamma==1){
	    // double x_foil, y_foil, z_foil;
	    // double x_gamma, y_gamma, z_gamma;
	    // x_gamma = vertex_gamma_tot[0][0];
	    // y_gamma = vertex_gamma_tot[0][1];
	    // z_gamma = vertex_gamma_tot[0][2];
	    // x_foil = vertex_3D_start_x[indices_calo_final.at(0)];
	    // y_foil = vertex_3D_start_y[indices_calo_final.at(0)];
	    // z_foil = vertex_3D_start_z[indices_calo_final.at(0)];
	    // time_of_flight_gamma = sqrt((x_foil-x_gamma)*(x_foil-x_gamma)+(y_foil-y_gamma)*(y_foil-y_gamma)+(z_foil-z_gamma)*(z_foil-z_gamma));
	    if(vertex_gamma[0]!=side_elec[indices_calo_final.at(0)] && same_side_elec==1){
	      opposite_side_e_gamma = 1;
	    }
	  }
	}
      }
      
      tree->Fill();
      event_number++;      
    }
  return dpp::base_module::PROCESS_SUCCESS;
}

