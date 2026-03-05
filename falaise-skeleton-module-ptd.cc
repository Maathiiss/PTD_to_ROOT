#include <bayeux/dpp/chain_module.h>
#include <falaise/snemo/datamodels/particle_track_data.h>
#include <falaise/snemo/datamodels/tracker_trajectory_solution.h>
#include <falaise/snemo/datamodels/tracker_trajectory_data.h>
#include <falaise/snemo/datamodels/tracker_clustering_data.h>
#include <falaise/snemo/datamodels/tracker_cluster.h>
#include <falaise/snemo/datamodels/precalibrated_data.h>
#include <falaise/snemo/datamodels/geomid_utils.h>
#include <falaise/snemo/datamodels/event_header.h>
#include <bayeux/mctools/base_step_hit.h>
#include <bayeux/mctools/simulated_data.h>
#include <falaise/snemo/datamodels/unified_digitized_data.h>
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
#include "TParameter.h"

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
  double extract_unix_start_time(int run_number) const;

  // Event processing function
  dpp::chain_module::process_status process (datatools::things & event);
  
private:
  int event_number;
  int  ptd_event_counter;
  int cellules_non_associated, cellules_SD_non_associated, number_of_kinks, run_number, nb_elec_real, trigger_id, nb_unfitted_cells;
  bool ptd_details, has_an_electron_and_positron, has_SD_electron_and_positron, has_SD_two_electrons, opposite_side_e_gamma, same_side_elec, same_cluster_elec, has_a_same, ttd_details, sd_calo_details, sd_tracker_details, hit_the_same_calo_hit, two_elec_more_350_keV, OM_are_neighbourg, has_kinks;

  double diff_time_elec, angle_3D_between_ep_em, delta_y_elec, delta_z_elec, energy_elec_sum, corrected_energy_elec_sum, time_of_flight_gamma, internal_theoretical_time_diff, external_theoretical_time_diff,t1_th, t2_th, start_run_time, end_run_time, delta_r_calo, kink_angle, one_kink_x, one_kink_y, one_kink_z, total_volume, unix_start_time=0, first_time, last_time, closest_gamma, closest_elec, closest_track, alpha_elec_time_diff, energy_elec_1, energy_elec_2, vertex_3D_projected_x, vertex_3D_projected_y, vertex_3D_projected_z, real_vertex_SD_x, real_vertex_SD_y, real_vertex_SD_z, closest_time_track;
  int nb_gamma, nb_elec_ptd_per_event,nb_elec_SD_per_event, nb_wire_hit;
  vector<string>  type_elec, g4_process, material, vertex_type, g4_material;
  vector<int> gamma_type, num_om, num_om_elec, track_number, num_gg, num_om_elec_f, num_om_gamma, num_om_track, indices_calo_final;
  vector<double> time_gamma, time_gamma_before, time_gamma_after, angle_SD, volume, total_step_length;
  vector<double> energy_gamma, energy_elec, corrected_energy_elec, time_elec, time_elec_alone,track_lenght, energy_gamma_after, vertex_SD_x, vertex_SD_y, vertex_SD_z, ellipse_source, time_track, energy_track, chi2_track, timestamp_elec_alone;
  vector<double> vertex_3D_start_x, vertex_3D_start_y, vertex_3D_start_z,vertex_3D_end_x, vertex_3D_end_y, vertex_3D_end_z, vertex_gamma, kink_x, kink_y, kink_z, vertex_3D_track_y, vertex_3D_track_z, mean_alpha_anodic_time, vertex_ttd_start_x, vertex_ttd_start_y, vertex_ttd_start_z, vertex_ttd_end_x, vertex_ttd_end_y, vertex_ttd_end_z, mean_alpha_anodic_timestamp;
  vector<int> side_elec, cluster_elec_num, calo_num, nb_kink_per_track;
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
  tree->Branch("run_number", &run_number);
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
  tree->Branch("real_vertex_SD_x", &real_vertex_SD_x);
  tree->Branch("real_vertex_SD_y", &real_vertex_SD_y);
  tree->Branch("real_vertex_SD_z", &real_vertex_SD_z);
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
  tree->Branch("vertex_ttd_start_x", &vertex_ttd_start_x);
  tree->Branch("vertex_ttd_start_y", &vertex_ttd_start_y);
  tree->Branch("vertex_ttd_start_z", &vertex_ttd_start_z);
  tree->Branch("vertex_3D_projected_x", &vertex_3D_projected_x);
  tree->Branch("vertex_3D_projected_y", &vertex_3D_projected_y);
  tree->Branch("vertex_3D_projected_z", &vertex_3D_projected_z);
  tree->Branch("vertex_3D_end_x", &vertex_3D_end_x);
  tree->Branch("vertex_3D_end_y", &vertex_3D_end_y);
  tree->Branch("vertex_3D_end_z", &vertex_3D_end_z);
  tree->Branch("vertex_ttd_end_x", &vertex_ttd_end_x);
  tree->Branch("vertex_ttd_end_y", &vertex_ttd_end_y);
  tree->Branch("vertex_ttd_end_z", &vertex_ttd_end_z); 
  tree->Branch("delta_r_calo", &delta_r_calo);
  tree->Branch("nb_elec_ptd_per_event", &nb_elec_ptd_per_event);
  tree->Branch("nb_elec_real", &nb_elec_real);
  tree->Branch("energy_elec", &energy_elec);
  tree->Branch("num_om", &num_om);
  tree->Branch("num_om_gamma", &num_om_gamma);
  tree->Branch("num_gg", &num_gg);
  tree->Branch("num_om_elec", &num_om_elec);
  tree->Branch("num_om_elec_f", &num_om_elec_f);
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
  tree->Branch("energy_elec_1", &energy_elec_1);
  tree->Branch("energy_elec_2", &energy_elec_2);
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
  tree->Branch("unix_start_time", &unix_start_time);
  tree->Branch("closest_gamma", &closest_gamma);
  tree->Branch("closest_elec", &closest_elec);
  tree->Branch("closest_track", &closest_track);
  tree->Branch("closest_time_track", &closest_time_track);
  tree->Branch("OM_are_neighbourg", &OM_are_neighbourg);
  tree->Branch("alpha_elec_time_diff", &alpha_elec_time_diff);
  tree->Branch("two_elec_more_350_keV", &two_elec_more_350_keV);
  tree->Branch("trigger_id", &trigger_id);
  tree->Branch("chi2_track", &chi2_track);
  tree->Branch("nb_unfitted_cells", &nb_unfitted_cells);
  tree->Branch("nb_kink_per_track", &nb_kink_per_track);
  tree->Branch("has_kinks", &has_kinks);

    
}


falaise_skeleton_module_ptd::~falaise_skeleton_module_ptd()
{
  double time = last_time - first_time;
  save_file->cd();
  TParameter<double> param("run_time", time);
  param.Write();
  tree->Write();
  save_file->Close();
  delete save_file;
  std::cout << "falaise_skeleton_module_ptd::~falaise_skeleton_module_ptd() called" << std::endl;
}


void falaise_skeleton_module_ptd::initialize (const datatools::properties & module_properties, datatools::service_manager &, dpp::module_handle_dict_type &)
{
  std::cout << "falaise_skeleton_module_ptd::initialize() called" << std::endl;
  event_number=0;
  first_time=0;
  last_time=0;
  
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

double snemo_run_time (int run_number)
{
  const char *cbd_base_path = "/sps/nemo/snemo/snemo_data/raw_data/CBD";
  std::vector<std::string> log_paths;
  log_paths.push_back(Form("%s/run-%d/snemo_trigger_run-%d.log", cbd_base_path, run_number, run_number));
  for (int crate=6; crate>=0; crate--)
    log_paths.push_back(Form("%s/run-%d/snemo_crate-%d_run-%d.log", cbd_base_path, run_number, crate, run_number));
  int unixtime;
  int run_start=0;
  for (const std::string & log_path : log_paths)
    {
      std::ifstream log_file (log_path);
      if (!log_file.is_open()) continue;
      std::string log_line;
      while (getline(log_file, log_line))
        {
          size_t unixtime_index = log_line.find("run.run_unixtime_ms=");
          if (unixtime_index == std::string::npos)
            continue;
          unixtime = std::stoi(log_line.substr(20));
          if (unixtime > run_start)
            run_start = unixtime;
            }
    }
  return unixtime;
}



double end_trip_time(int run_number) {
  string filename;
  if(run_number<1799){
    filename = "/sps/nemo/snemo/snemo_data/reco_data/UDD_betabeta_v1.list";
  }
  else{
    filename = "/sps/nemo/snemo/snemo_data/reco_data/UDD_betabeta_v2.list";
  }
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        int run; double start,dur,stop; std::string comment;
        if (!(iss >> run >> start >> dur >> stop >> comment)) continue;
        if (run == run_number) return static_cast<int>(stop);
    }
    return -1; // run non trouvé
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


inline void decode_om_type0(int num_om,int& side,int& col,int& row)
{
    side = num_om / 260;
    int reste = num_om % 260;
    col = reste / 13;
    row = reste % 13;
}

bool compute_neighbourg(int num_om_1, int num_om_2)
{
    int side1, col1, row1;
    int side2, col2, row2;
    decode_om_type0(num_om_1, side1, col1, row1);
    decode_om_type0(num_om_2, side2, col2, row2);
    if (side1 != side2) return false;
    int dcol = std::abs(col1 - col2);
    int drow = std::abs(row1 - row2);
    if (dcol <= 1 && drow <= 1 && (dcol + drow) > 0)
        return true;
    return false;
}

 

// Struct to represent a 3D vector 
struct Vec3 {
    double x, y, z;

    Vec3 operator-(const Vec3& v) const { return {x - v.x, y - v.y, z - v.z}; }
    Vec3 operator+(const Vec3& v) const { return {x + v.x, y + v.y, z + v.z}; }
    Vec3 operator*(double k) const { return {x * k, y * k, z * k}; }
};
// Scalar product
double dot(const Vec3& a, const Vec3& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}
// Norme
double norm(const Vec3& v) {
    return sqrt(dot(v, v));
}
//Compute the closest points between two lines 
bool closestPointsBetweenLines(
    const Vec3& A, const Vec3& B,
    const Vec3& A2, const Vec3& B2,
    Vec3& C, Vec3& C2)
{
    Vec3 u = B - A;     // direction de la première droite
    Vec3 v = B2 - A2;   // direction de la deuxième droite
    Vec3 w0 = A - A2;   // vecteur entre les points de départ
    double a = dot(u, u);     // u·u
    double b = dot(u, v);     // u·v
    double c = dot(v, v);     // v·v
    double d = dot(u, w0);    // u·w0
    double e = dot(v, w0);    // v·w0
    double denom = a*c - b*b; // dénominateur
    // Si denom == 0 --> droites parallèles
    if (fabs(denom) < 1e-9) {
        return false; // parallèle ou quasi-parallèle
    }
    double t = (b*e - c*d) / denom;
    double s = (a*e - b*d) / denom;
    C  = A  + u * t;   // point sur D
    C2 = A2 + v * s;   // point sur D'
    return true;
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


double falaise_skeleton_module_ptd::extract_unix_start_time(int run_number) const
{
    std::string cmd = "grep ^" + std::to_string(run_number) + " /sps/nemo/scratch/chauveau/commissioning/software/run-sync-time.txt | awk '{print $2}'";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) {
        std::cerr << "Erreur: impossible d'exécuter grep\n";
        return -1;
    }

    double timestamp = -1;
    fscanf(pipe, "%lf", &timestamp);
    pclose(pipe);
    return timestamp;
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
  chi2_track.clear();
  has_kinks=0;
  nb_kink_per_track.clear();
  indices_calo_final.clear();
  vertex_gamma.clear();
  nb_elec_ptd_per_event=0;
  nb_elec_real = 0;
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
  energy_track.clear();
  calo_num.clear();
  num_om.clear();
  num_om_gamma.clear();
  num_gg.clear();
  num_om_elec.clear();
  num_om_elec_f.clear();
  num_om_track.clear();
  corrected_energy_elec.clear();
  time_elec.clear();
  time_elec_alone.clear();
  timestamp_elec_alone.clear();
  time_track.clear();
  track_lenght.clear();
  energy_gamma.clear();
  mean_alpha_anodic_time.clear();
  mean_alpha_anodic_timestamp.clear();
  vertex_3D_start_x.clear();
  vertex_3D_start_y.clear();
  vertex_3D_start_z.clear();
  vertex_ttd_start_x.clear();
  vertex_ttd_start_y.clear();
  vertex_ttd_start_z.clear();
  vertex_ttd_end_x.clear();
  vertex_ttd_end_y.clear();
  vertex_ttd_end_z.clear();
  vertex_3D_track_y.clear();
  vertex_3D_track_z.clear();
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
  nb_unfitted_cells=0;
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
  real_vertex_SD_x=0.0;
  real_vertex_SD_y=0.0;
  real_vertex_SD_z=0.0;
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
  run_number=0;
  closest_gamma=1e6;
  closest_elec=1e6;
  closest_time_track = 1e6;
  closest_track = 1e6;
  bool projected_vertex = false;
  bool has_tracks=false;
  std::vector<double> x_alpha;
  std::vector<double> y_alpha;
  std::vector<double> z_alpha;
  
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
    real_vertex_SD_x = SD.get_vertex()[0]/CLHEP::mm; 
    real_vertex_SD_y = SD.get_vertex()[1]/CLHEP::mm;
    real_vertex_SD_z = SD.get_vertex()[2]/CLHEP::mm;
  }
 
  
    
  // PTD extraction
  if (ptd_details)    
    {

      std::vector<std::vector<std::pair<double,double>>> gg_elec_pos;
      std::vector<std::vector<std::pair<double,double>>> global_gg_pos;
      int indice_elec=0;
      const snemo::datamodel::event_header & eh = event.get<snemo::datamodel::event_header>("EH");
      run_number = eh.get_id().get_run_number();      
      if(event.has("SD")==0){
	const snemo::datamodel::unified_digitized_data & UDD = event.get<snemo::datamodel::unified_digitized_data>("UDD");
	const auto& trigger_ids = UDD.get_origin_trigger_ids();
	trigger_id = -1;
	if (!trigger_ids.empty()) {
	  trigger_id = *trigger_ids.begin(); 
	}
      }

      if(eh.has_timestamp()){
	double time = (eh.get_timestamp().get_seconds()+eh.get_timestamp().get_picoseconds()/1E12)*1e-9;
	if(time>end_trip_time(run_number) && end_trip_time(run_number)!=0) {
	  return dpp::base_module::PROCESS_SUCCESS;
	}
        last_time = time;
        if(first_time==0){
          first_time = time;
        }
      }
      if(unix_start_time==0){
	if(run_number<1556){
          unix_start_time = snemo_run_time(run_number);
        }
        else{
          unix_start_time=extract_unix_start_time(run_number);
        }
      }

      
      const snemo::datamodel::tracker_trajectory_data & TTD = event.get<snemo::datamodel::tracker_trajectory_data>("TTD");
      //if(TTD.get_solution_id()==0){
      const snemo::datamodel::tracker_trajectory_solution & ttd_solution = TTD.get_default_solution();
      for (const auto & trajectory : ttd_solution.get_trajectories()) {
        if(trajectory->get_fit_infos().has_chi2() && trajectory->get_fit_infos().has_ndof()){
          chi2_track.push_back(trajectory->get_fit_infos().get_chi2() / trajectory->get_fit_infos().get_ndof());
        }
      }
          const snemo::datamodel::tracker_clustering_data & TCD = event.get<snemo::datamodel::tracker_clustering_data>("TCD");
	  const snemo::datamodel::precalibrated_data & pCD = event.get<snemo::datamodel::precalibrated_data>("pCD");
      const snemo::datamodel::tracker_clustering_solution & tcd_solution = TCD.get_default();
      nb_unfitted_cells=tcd_solution.get_unclustered_hits().size();

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
	  time_gamma.push_back(it_hit->get_time()/ CLHEP::ns);
	  nb_gamma++;
	  num_om_gamma.push_back(snemo::datamodel::om_num(it_hit->get_geom_id()));
	}
      }
     
      //electron start
      //vector<vector<double>> vertex_gamma_tot;
    
      vector<int> cluster_id;
      vector<int> cluster_id_tot;
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
	    
	  global_gg_pos.push_back({});
	  auto &gg_vec = global_gg_pos.back();
	  const auto gg = particle->get_trajectory_handle()->get_cluster().hits();
	  for (const auto hits : gg) {
	    if (hits->has_xy()) {
	      gg_vec.push_back({ hits->get_x()/CLHEP::mm, hits->get_y()/CLHEP::mm });
	    }
	  }

	  
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
	       if(vertex->is_on_reference_source_plane()){
                x_foil = vertex->get_spot().get_position().getX();
		y_foil = vertex->get_spot().get_position().getY();
                z_foil = vertex->get_spot().get_position().getZ();
                vertex_close_to_the_source=1;
              }
	       else if(vertex->is_on_source_foil() && vertex_close_to_the_source==0){
		 x_foil = vertex->get_spot().get_position().getX();
                y_foil = vertex->get_spot().get_position().getY();
                z_foil = vertex->get_spot().get_position().getZ();
                vertex_close_to_the_source=1;
              }	    
	       else if(vertex->is_on_main_calorimeter() /*|| vertex->is_on_x_calorimeter()*/) //we forced MW only analysis
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
	    }//end vertices loop

	    int source_num;
	    const auto& calorimeter_hits = particle->get_associated_calorimeter_hits();
	    if(particle->get_associated_calorimeter_hits().size()>1) continue;
	    if(vertex_close_to_the_source && vertex_associated_to_a_calo && particle->get_associated_calorimeter_hits().size() == 1 && calorimeter_hits[0]>0.050){
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
	      vertex_ttd_start_x.push_back(particle->get_trajectory_handle()->get_pattern().get_first()[0]);
	      vertex_ttd_start_y.push_back(particle->get_trajectory_handle()->get_pattern().get_first()[1]);
	      vertex_ttd_start_z.push_back(particle->get_trajectory_handle()->get_pattern().get_first()[2]);
	      vertex_3D_end_x.push_back(x_calo);
	      vertex_3D_end_y.push_back(y_calo);
	      vertex_3D_end_z.push_back(z_calo);
 	      vertex_ttd_end_x.push_back(particle->get_trajectory_handle()->get_pattern().get_last()[0]);
	      vertex_ttd_end_y.push_back(particle->get_trajectory_handle()->get_pattern().get_last()[1]);
	      vertex_ttd_end_z.push_back(particle->get_trajectory_handle()->get_pattern().get_last()[2]);	      
	      cluster_id_tot.push_back(particle->get_trajectory().get_cluster().get_hit_id());
	      cluster_id.push_back(particle->get_trajectory().get_cluster().get_hit_id());
	      if(nb_kinks>0){//looking if this particle get kinks
		kinks_vec.push_back({trajectory_pattern.get_kink(0).getX(),trajectory_pattern.get_kink(0).getY(),trajectory_pattern.get_kink(0).getZ()});
	      }
	      else{
		kinks_vec.push_back({});
	      }
	      calo_num.push_back(snemo::datamodel::om_num(calorimeter_hits[0]->get_geom_id()));
	      const datatools::handle<snemo::datamodel::precalibrated_calorimeter_hit> & pcd_calo_hit = pCD.calorimeter_hits()[calorimeter_hits[0]->get_hit_id()];
	      timestamp_elec_alone.push_back(pcd_calo_hit->get_time()/CLHEP::ns);
	      energy_elec.push_back(calorimeter_hits[0]->get_energy());
	      time_elec.push_back(calorimeter_hits[0]->get_time()/ CLHEP::ns);
	      time_elec_alone.push_back(calorimeter_hits[0]->get_time()/ CLHEP::ns);
	      side_elec.push_back(calorimeter_hits[0]->get_geom_id().get(1));
	      num_om_elec.push_back(snemo::datamodel::om_num(calorimeter_hits[0]->get_geom_id()));
	      num_om.push_back(snemo::datamodel::om_num(calorimeter_hits[0]->get_geom_id()));
  	      gg_elec_pos.push_back({});
	      auto &pos_vec = gg_elec_pos.back();  
	      const auto gg = particle->get_trajectory_handle()->get_cluster().hits();
	      for (const auto hits : gg) {
		if (hits->has_xy()) {				   
		  pos_vec.push_back({ hits->get_x()/CLHEP::mm, hits->get_y()/CLHEP::mm });
		}
	      }
			      
	      if (particle->has_trajectory()) {
		const auto& trajectory = particle->get_trajectory();
		if (trajectory.has_cluster()) {
		  cluster_elec_num.push_back(trajectory.get_cluster().get_hit_id());
		}
	      }
	      //const auto& trajectory_pattern_elec = particle->get_trajectory_handle()->get_pattern();
              int nb_kinks_elec_save = trajectory_pattern_elec.number_of_kinks();
	      nb_kink_per_track.push_back(nb_kinks_elec_save);
	      nb_elec_ptd_per_event++;
	      indices_calo_final.push_back(indice_elec);
	      indice_elec++;
	    }
	    
	    else if (vertex_associated_to_a_calo && vertex_close_to_the_source==0 && particle->get_associated_calorimeter_hits().size()==1 && calorimeter_hits[0]>0.050){//calo search 
	      //const auto& calorimeter_hits = particle->get_associated_calorimeter_hits();
	      // const auto& calorimeter_hits = particle->get_associated_calorimeter_hits();     
	      // if(calorimeter_hits[0]->get_energy()<0.350) continue;	     
	      num_om.push_back(snemo::datamodel::om_num(calorimeter_hits[0]->get_geom_id()));
	      num_om_track.push_back(snemo::datamodel::om_num(calorimeter_hits[0]->get_geom_id()));
              time_elec.push_back(calorimeter_hits[0]->get_time()/ CLHEP::ns);
	      time_track.push_back(calorimeter_hits[0]->get_time()/ CLHEP::ns);
	      energy_track.push_back(calorimeter_hits[0]->get_energy());
	      cluster_id_tot.push_back(particle->get_trajectory().get_cluster().get_hit_id());
	    }
	  
	    else if(particle->get_associated_calorimeter_hits().size()==0 || calorimeter_hits[0]<0.050){//track search
	      vector<int> alpha_gg_num;
	      const auto gg = particle->get_trajectory_handle()->get_cluster().hits();
	      for (const auto hits : gg) {
		alpha_gg_num.push_back(snemo::datamodel::gg_num(hits->get_geom_id()));
	      }
	      double sum_cell = 0;
	      int n_cell = 0;
	      const auto pCD_tracker_hits = pCD.tracker_hits();
	      for(const auto pCD_track_hits : pCD_tracker_hits){
		if(std::find(alpha_gg_num.begin(), alpha_gg_num.end(), snemo::datamodel::gg_num(pCD_track_hits->get_geom_id())) != alpha_gg_num.end() ){
		  sum_cell+=pCD_track_hits->get_anodic_time()/ CLHEP::ns;
		  n_cell++;
		}
	      }
	      double new_mean = (n_cell>0 ? sum_cell/n_cell : 0);
	      mean_alpha_anodic_timestamp.push_back(new_mean);
 	      projected_vertex = vertex_close_to_the_source;
	      if(vertex_close_to_the_source){
		x_alpha.push_back(x_foil);
		y_alpha.push_back(y_foil);
		z_alpha.push_back(z_foil);
	      }
	      else{
		x_alpha.push_back(particle->get_trajectory_handle()->get_pattern().get_first()[0]);
		x_alpha.push_back(particle->get_trajectory_handle()->get_pattern().get_last()[0]);
		y_alpha.push_back(particle->get_trajectory_handle()->get_pattern().get_first()[1]);
		y_alpha.push_back(particle->get_trajectory_handle()->get_pattern().get_last()[1]);
		z_alpha.push_back(particle->get_trajectory_handle()->get_pattern().get_first()[2]);
		z_alpha.push_back(particle->get_trajectory_handle()->get_pattern().get_last()[2]);
		
	      }
	    }//end track search	   
	  }//end e-	  	  
	  
	  else{
	    cellules_non_associated++;
	  }
	}//end loop on particle
    
      
      if (nb_elec_ptd_per_event < 2) {
	event_number++;      
	return dpp::base_module::PROCESS_SUCCESS;//we stop the search if we have less than 2 tracks
      }

      // now we associate tracks by pairs by computing the closest time
      std::vector<std::vector<int>> pairs; // this vector will contain all double beta candidates
      std::vector<int> remaining = indices_calo_final;
      while (remaining.size() >= 2) {
	double best_dist2 = std::numeric_limits<double>::max();
	int best_i = -1;
	int best_j = -1;
	for (int i = 0; i < (int)remaining.size(); ++i) {
	  for (int j = i + 1; j < (int)remaining.size(); ++j) {
            int a = remaining[i];
            int b = remaining[j];
	    double dt = std::abs(time_elec_alone[a] - time_elec_alone[b]);
	    if (dt < best_dist2) {
	      best_dist2 = dt;
	      best_i = i;
	      best_j = j;
	    }
	  }
	}
	// add the 2 tracks found
	pairs.push_back({ remaining[best_i], remaining[best_j] });
	if (best_i > best_j) {
	  std::swap(best_i, best_j);
	}
	remaining.erase(remaining.begin() + best_j);
	remaining.erase(remaining.begin() + best_i);//remove the 2 indexes from the vector
      }

      
      std::vector<double> min_dt_for_pair; 
      min_dt_for_pair.resize(pairs.size(), 1e9);  
      for (size_t p = 0; p < pairs.size(); ++p) {
	int best_a = pairs[p][0];
	int best_b = pairs[p][1];
	double tA = time_elec_alone[best_a];
	double tB = time_elec_alone[best_b];
	double min_dt = 1e9;  // min dt for this pair
	// --- 1. Compare with other electrons ---
	for (size_t idx = 0; idx < cluster_id.size(); ++idx) {
	  if (idx == (size_t)best_a || idx == (size_t)best_b) continue;
          double dtA = std::abs(time_elec_alone[idx] - tA);
	  double dtB = std::abs(time_elec_alone[idx] - tB);
	  min_dt = std::min(min_dt, std::min(dtA, dtB));
	}
	// --- 2. Compare with time_track (no cluster check needed) ---
	for (size_t idx = 0; idx < time_track.size(); ++idx) {
	  double dtA = std::abs(time_track[idx] - tA);
	  double dtB = std::abs(time_track[idx] - tB);
	  min_dt = std::min(min_dt, std::min(dtA, dtB));
	}
	// store result
	min_dt_for_pair[p] = min_dt;
      }

      
    
      //nb_elec_ptd_per_event=indices_calo_final.size();
      for (size_t p = 0; p < pairs.size(); ++p) {
	energy_elec_1 =100;
	energy_elec_2=100;
	two_elec_more_350_keV=false;
	OM_are_neighbourg=false;
	has_kinks=false;
	int best_a = pairs[p][0];
	int best_b = pairs[p][1];
	num_om_elec_f.push_back(num_om_elec[best_a]);
	num_om_elec_f.push_back(num_om_elec[best_b]);

	  if(compute_neighbourg(num_om_elec[best_a],num_om_elec[best_b]) == 1){
	    OM_are_neighbourg = true;
	  }
	  if(nb_kink_per_track[best_a]!=0 || nb_kink_per_track[best_b]!=0){has_kinks=1;}
	  // cout<<event_number<<endl;
	  // cout<<num_om_elec[best_a]<<" "<<num_om_elec[best_b]<<" "<<OM_are_neighbourg<<endl;
	  if(num_om_elec[best_a] == num_om_elec[best_b]) hit_the_same_calo_hit = true;
	//if(energy_elec[best_a] == 0 || energy_elec[best_b] == 0) cellules_non_associated++;
	if(cluster_elec_num[best_a] == cluster_elec_num[best_b]) same_cluster_elec = 1;
	if(side_elec[best_a] == side_elec[best_b]) same_side_elec = 1;
	//spatial computation 
	delta_y_elec = abs(vertex_3D_start_y[best_a] - vertex_3D_start_y[best_b]);
	delta_z_elec = abs(vertex_3D_start_z[best_a] - vertex_3D_start_z[best_b]);
	delta_r_calo = sqrt(pow(vertex_3D_end_y[best_a] - vertex_3D_end_y[best_b],2) +
			    pow(vertex_3D_end_z[best_a] - vertex_3D_end_z[best_b],2));
	energy_elec_sum = energy_elec[best_a] + energy_elec[best_b];
	energy_elec_1 = energy_elec[best_a];
	energy_elec_2 = energy_elec[best_b];
	closest_elec = min_dt_for_pair[p];
	if(energy_elec[best_a] > 0.050 && energy_elec[best_b] > 0.050) two_elec_more_350_keV = true;
	//angle computation
	double start_p1[3] = {vertex_3D_start_x[best_a], vertex_3D_start_y[best_a], vertex_3D_start_z[best_a]};
	double start_p2[3] = {vertex_3D_start_x[best_b], vertex_3D_start_y[best_b], vertex_3D_start_z[best_b]};
	double end_p1[3] = {vertex_3D_end_x[best_a], vertex_3D_end_y[best_a], vertex_3D_end_z[best_a]};
	double end_p2[3] = {vertex_3D_end_x[best_b], vertex_3D_end_y[best_b], vertex_3D_end_z[best_b]};
 	if(!kinks_vec[best_a].empty()) {
	  end_p1[0] = kinks_vec[best_a][0];
	  end_p1[1] = kinks_vec[best_a][1];
	  end_p1[2] = kinks_vec[best_a][2];
	}
	if(!kinks_vec[best_b].empty()) {
	  end_p2[0] = kinks_vec[best_b][0];
	  end_p2[1] = kinks_vec[best_b][1];
	  end_p2[2] = kinks_vec[best_b][2];
	}
	angle_3D_between_ep_em = std::abs(compute_angle(start_p1, end_p1, start_p2, end_p2));

	//projection computation
 	double segA_start[3] = {vertex_ttd_start_x[best_a], vertex_ttd_start_y[best_a], vertex_ttd_start_z[best_a]};
	double segB_start[3] = {vertex_ttd_start_x[best_b], vertex_ttd_start_y[best_b], vertex_ttd_start_z[best_b]};

	double segA_end[3]   = {vertex_ttd_end_x[best_a], vertex_ttd_end_y[best_a], vertex_ttd_end_z[best_a]};
	double segB_end[3]   = {vertex_ttd_end_x[best_b], vertex_ttd_end_y[best_b], vertex_ttd_end_z[best_b]};
	if(!kinks_vec[best_a].empty()) {
	  segA_end[0] = kinks_vec[best_a][0];
	  segA_end[1] = kinks_vec[best_a][1];
	  segA_end[2] = kinks_vec[best_a][2];
	}
	if(!kinks_vec[best_b].empty()) {
	  segB_end[0] = kinks_vec[best_b][0];
	  segB_end[1] = kinks_vec[best_b][1];
	  segB_end[2] = kinks_vec[best_b][2];
	}
	Vec3 A  = {segA_start[0], segA_start[1], segA_start[2]};
	Vec3 B  = {segA_end[0],   segA_end[1],   segA_end[2]};
	Vec3 A2 = {segB_start[0], segB_start[1], segB_start[2]};
	Vec3 B2 = {segB_end[0],   segB_end[1],   segB_end[2]};
	Vec3 C, C2;
	if (closestPointsBetweenLines(A, B, A2, B2, C, C2)) {
	  // Distance minimale entre les droites
	  double distance = norm(C - C2);
	  // Point du milieu du segment CC'
	  Vec3 M = (C + C2) * 0.5; 
	  vertex_3D_projected_x = M.x;
	  vertex_3D_projected_y = M.y;
	  vertex_3D_projected_z = M.z;
	}
	diff_time_elec = time_elec_alone[best_a] - time_elec_alone[best_b];

	//TOF computation
	double new_c = 299.792458; // mm/ns
	double track_length_1 = sqrt(pow(start_p1[0]-end_p1[0],2) + pow(start_p1[1]-end_p1[1],2) + pow(start_p1[2]-end_p1[2],2));
	double track_length_2 = sqrt(pow(start_p2[0]-end_p2[0],2) + pow(start_p2[1]-end_p2[1],2) + pow(start_p2[2]-end_p2[2],2));
	double beta_1 = sqrt(energy_elec[best_a]*(energy_elec[best_a]+2*0.511)) / (energy_elec[best_a]+0.511);
	double beta_2 = sqrt(energy_elec[best_b]*(energy_elec[best_b]+2*0.511)) / (energy_elec[best_b]+0.511);
	double t_1 = track_length_1/(beta_1*new_c);
	double t_2 = track_length_2/(beta_2*new_c);

	internal_theoretical_time_diff = (t_1 - t_2) - diff_time_elec;
	external_theoretical_time_diff = std::abs(diff_time_elec) - (t_1 + t_2);
	t1_th = t_1;
	t2_th = t_2;
	track_lenght.push_back(track_length_1);
	track_lenght.push_back(track_length_2);

	//gamma loop -> compute time differences with all gammas
	std::vector<double> gamma_dt;     
	std::vector<double> gamma_energy; 
	double tA = time_elec_alone[best_a];
	double tB = time_elec_alone[best_b];
	double timestampA = timestamp_elec_alone[best_a];
	double timestampB = timestamp_elec_alone[best_b];
	for (size_t g = 0; g < time_gamma.size(); ++g) {
	  if (energy_gamma[g] < 0.050) continue;
	  double tγ = time_gamma[g];
	  double dtA = std::abs(tγ - tA);
	  double dtB = std::abs(tγ - tB);
	  double dt_min = std::min(dtA, dtB);
	  gamma_dt.push_back(dt_min);
	  gamma_energy.push_back(energy_gamma[g]);
	  if (nb_gamma == 1 &&
	      vertex_gamma[0] != side_elec[best_a] &&
	      same_side_elec == 1)
	    {
	      opposite_side_e_gamma = 1;
	    }
	}
	// --- find closest gamma wrt the e− pair ------------------------------------
	closest_gamma = 1e9;  
	for (double dt : gamma_dt) {
	  if (std::abs(dt) < std::abs(closest_gamma))
	    closest_gamma = dt;
	}

	//remove tracks in time with double beta events
	double closest_alpha_dt_delayed = 1e9;
	for (double t_alpha : mean_alpha_anodic_timestamp) {
	  double dt = std::min(std::abs(t_alpha - timestampA), std::abs(t_alpha - timestampB));
	  if (dt < closest_alpha_dt_delayed)
	    closest_alpha_dt_delayed = dt;
	} 	
	closest_time_track = closest_alpha_dt_delayed * 1e-3;
      
	
	//remove alpha tracks not in time but close in distance
	double closest_alpha_dist = 1e9;
	//if(projected_vertex==1){
	  for (size_t i = 0; i < x_alpha.size(); ++i) {
	    // distance à A
	    double dxA = x_alpha[i] - start_p1[0];
	    double dyA = y_alpha[i] - start_p1[1];
	    double dzA = z_alpha[i] - start_p1[2];
	    double dist1 = std::sqrt(dxA*dxA + dyA*dyA + dzA*dzA);
	    // distance à B
	    double dxB = x_alpha[i] - start_p2[0];
	    double dyB = y_alpha[i] - start_p2[1];
	    double dzB = z_alpha[i] - start_p2[2];
	    double dist2 = std::sqrt(dxB*dxB + dyB*dyB + dzB*dzB);
	    double dist_min = std::min(dist1, dist2);
	    if (dist_min < closest_alpha_dist) {
	      closest_alpha_dist = dist_min;
	    }
	  }
	  closest_track = closest_alpha_dist; // distance minimale en mm
	  //}
      
      	tree->Fill();//tree is fill for every pairs
            
      }// end pair computation
      event_number++;
    }
  return dpp::base_module::PROCESS_SUCCESS;
}

