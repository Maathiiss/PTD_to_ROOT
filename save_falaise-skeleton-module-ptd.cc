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

  // Event processing function
  dpp::chain_module::process_status process (datatools::things & event);
  
private:
  int  ptd_event_counter;
  int cellules_non_associated, number_of_kinks;
  bool ptd_details, has_an_electron_and_positron, has_SD_electron_and_positron, has_SD_two_electrons, opposite_side_e_gamma, same_side_elec, same_cluster_elec, has_a_same, ttd_details, sd_calo_details, sd_tracker_details, hit_the_same_calo_hit;
  int event_number;
  double diff_time_elec, angle_2D_between_ep_em, angle_3D_between_ep_em, delta_y_elec, delta_z_elec, energy_elec_sum, corrected_energy_elec_sum, time_of_flight_gamma, internal_theoretical_time_diff, external_theoretical_time_diff,t1_th, t2_th, start_run_time, end_run_time, delta_r_calo;
  int nb_gamma, nb_elec_ptd_per_event, kink_angle;
  vector<string>  type_elec;
  vector<int> gamma_type, num_om;
  vector<double> time_gamma, time_gamma_before, time_gamma_after;
  vector<double> energy_gamma, energy_elec, corrected_energy_elec, time_elec,track_lenght, energy_gamma_after;
  vector<double> vertex_3D_start_x, vertex_3D_start_y, vertex_3D_start_z,vertex_3D_end_x, vertex_3D_end_y, vertex_3D_end_z, vertex_gamma, kink_x, kink_y, kink_z;
  vector<int> side_elec, cluster_elec_num;
  TFile *save_file;
  TTree *tree;

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
  tree->Branch("angle_2D_between_ep_em", &angle_2D_between_ep_em);
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
  vertex_gamma.clear();
  nb_elec_ptd_per_event=0;
  cellules_non_associated=0;
  gamma_type.clear();
  time_gamma.clear();
  time_gamma_before.clear();
  time_gamma_after.clear();
  energy_gamma_after.clear();
  type_elec.clear();
  side_elec.clear();
  cluster_elec_num.clear();
  energy_elec.clear();
  num_om.clear();
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
  angle_2D_between_ep_em = 0.0;
  angle_3D_between_ep_em = 0.0;
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
  number_of_kinks=0;
  hit_the_same_calo_hit=false;
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
  }


  // PTD extraction
  if (ptd_details)    
    {
      
      std::vector<std::vector<int>> same_clusters;
      std::vector<std::vector<double>> kinks_vec;
      int nb_kinks_save=0;
      int count_equal_clusters = 0;
      int nb_clusters=0;
      //check unfitted cluster
      
      const snemo::datamodel::tracker_trajectory_data & TTD = event.get<snemo::datamodel::tracker_trajectory_data>("TTD"); 
      if (ttd_details)
    {
      const snemo::datamodel::tracker_trajectory_solution & ttd_solution = TTD.get_default_solution();
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
	      if(/*vertex->is_on_reference_source_plane()*/ vertex->is_on_source_foil()){
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
	    if(vertex_close_to_the_source && vertex_associated_to_a_calo){
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
		kinks_vec.push_back({trajectory_pattern.get_kink(0).getY(),trajectory_pattern.get_kink(0).getZ()});
	      }
	      else{
		kinks_vec.push_back({});
	      }
	      
	      energy_elec.push_back(calorimeter_hits[0]->get_energy());
              time_elec.push_back(calorimeter_hits[0]->get_time());
	      side_elec.push_back(calorimeter_hits[0]->get_geom_id().get(1));
	      num_om.push_back(snemo::datamodel::om_num(calorimeter_hits[0]->get_geom_id()));
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
      if(nb_kinks_save==1 && nb_elec_ptd_per_event==1){
        double x1 = vertex_3D_start_x[0];
        double y1 = vertex_3D_start_y[0];
        double z1 = vertex_3D_start_z[0];
        double x2 = kink_x[0];
        double y2 = kink_y[0];
        double z2 = kink_z[0];
	double x3 = vertex_3D_end_x[0];
        double y3 = vertex_3D_end_y[0];
        double z3 = vertex_3D_end_z[0];
	double dot_product = (x2 - x1) * (x3 - x2) + (y2 - y1) * (y3 - y2) + (z2 - z1) * (z3 - z2);
	double norm_u = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));
        double norm_v = sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2) + (z3 - z2) * (z3 - z2));
        kink_angle=(180.0/M_PI)*acos(dot_product / (norm_u * norm_v));
      }
      


      




      //check number if electrons clusters are the same
      const snemo::datamodel::tracker_clustering_data & TCD = event.get<snemo::datamodel::tracker_clustering_data>("TCD");
      for (const auto &tcd_solution : TCD.solutions()) {
	const auto &tcd_clusters = tcd_solution->get_clusters();
	nb_clusters = tcd_clusters.size();	  
	for (size_t i : cluster_id) {
	  const auto &tcd_cluster = tcd_clusters[i];	    
	  for (size_t j : cluster_id){
	    if(j<=i){continue;}
	    const auto &tcd_cluster_same = tcd_clusters[j];	   	    
	    if (tcd_cluster->hits() == tcd_cluster_same->hits()) { //same size clusters
	      bool is_equal = true;
	      const auto &hits = tcd_cluster->hits();
	      const auto &hits_same = tcd_cluster_same->hits();		
	      auto hit_it = hits.begin();
	      auto hit_same_it = hits_same.begin();		
	      while (hit_it != hits.end() && hit_same_it != hits_same.end() && is_equal) {
		if ((*hit_it)->get_id() != (*hit_same_it)->get_id()) {
		  is_equal = false;
		  break;
		}
		++hit_it;
		++hit_same_it;
	      }
	      if (is_equal) {
		count_equal_clusters++;
		const int id_cluster = tcd_cluster->get_cluster_id();
		const int id_cluster_same = tcd_cluster_same->get_cluster_id();
		same_clusters.push_back({id_cluster,id_cluster_same});
	      }		
	    }
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
    

      nb_elec_ptd_per_event=indices_calo_final.size();
      if(nb_elec_ptd_per_event==2){	
	if(num_om[indices_calo_final.at(0)]==num_om[indices_calo_final.at(1)]){
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
	  y2 = kinks_vec[indices_calo_final.at(0)][0];
	  z2 = kinks_vec[indices_calo_final.at(0)][1];
	}
	if(kinks_vec[indices_calo_final.at(1)].size()>0){
	  y2 = kinks_vec[indices_calo_final.at(1)][0];
	  z2 = kinks_vec[indices_calo_final.at(1)][1];
	}
	double dot_product = (x2 - x1) * (x2_ - x1_) + (y2 - y1) * (y2_ - y1_);
	double norm_v1 = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
	double norm_v2 = sqrt(pow(x2_ - x1_, 2) + pow(y2_ - y1_, 2));
	//we want the 2d angle
	angle_2D_between_ep_em = (180.0/M_PI)*acos(dot_product/(norm_v1*norm_v2));
	
	angle_3D_between_ep_em = (180.0/M_PI)*acos(((x2-x1)*(x2_-x1_)+(y2-y1)*(y2_-y1_)+(z2-z1)*(z2_-z1_))/(sqrt(pow((x2-x1),2)+pow((y2-y1),2)+pow((z2-z1),2)) * sqrt(pow((x2_-x1_),2)+pow((y2_-y1_),2)+pow((z2_-z1_),2))));
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

