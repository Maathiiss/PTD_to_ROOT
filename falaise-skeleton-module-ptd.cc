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
#include <cstdint>
#include "TTree.h"
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include "TParameter.h"

using namespace std;

// ============================================================================
// Helper structures and functions
// ============================================================================

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator*(double k) const { return Vec3(x * k, y * k, z * k); }
};

inline double dot(const Vec3& a, const Vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline double norm(const Vec3& v) { return sqrt(dot(v, v)); }

bool closestPointsBetweenLines(const Vec3& A, const Vec3& B, const Vec3& A2, 
                                const Vec3& B2, Vec3& C, Vec3& C2) {
    Vec3 u = B - A, v = B2 - A2, w0 = A - A2;
    double a = dot(u, u), b = dot(u, v), c = dot(v, v);
    double d = dot(u, w0), e = dot(v, w0);
    double denom = a*c - b*b;
    if (fabs(denom) < 1e-9) return false;
    double t = (b*e - c*d) / denom, s = (a*e - b*d) / denom;
    C = A + u * t;
    C2 = A2 + v * s;
    return true;
}

double computeAngle(const double* start1, const double* end1,
                    const double* start2, const double* end2) {
    double v1[3] = {end1[0]-start1[0], end1[1]-start1[1], end1[2]-start1[2]};
    double v2[3] = {end2[0]-start2[0], end2[1]-start2[1], end2[2]-start2[2]};
    double dp = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
    double n1 = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
    double n2 = sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
    if (n1 < 1e-9 || n2 < 1e-9) return 0.0;
    return acos(dp / (n1 * n2)) * (180.0 / M_PI);
}

inline void decodeOM(int om_id, int& side, int& col, int& row) {
    side = om_id / 260;
    int rem = om_id % 260;
    col = rem / 13;
    row = rem % 13;
}

bool areOMNeighbors(int om1, int om2) {
    int s1, c1, r1, s2, c2, r2;
    decodeOM(om1, s1, c1, r1);
    decodeOM(om2, s2, c2, r2);
    if (s1 != s2) return false;
    int dcol = abs(c1 - c2), drow = abs(r1 - r2);
    return (dcol <= 1 && drow <= 1 && (dcol + drow) > 0);
}

double snemo_run_time(int run_number) {
    const char *cbd_base_path = "/sps/nemo/snemo/snemo_data/raw_data/CBD";
    std::vector<std::string> log_paths;
    log_paths.push_back(Form("%s/run-%d/snemo_trigger_run-%d.log", cbd_base_path, run_number, run_number));
    for (int crate=6; crate>=0; crate--)
        log_paths.push_back(Form("%s/run-%d/snemo_crate-%d_run-%d.log", cbd_base_path, run_number, crate, run_number));
    int unixtime, run_start=0;
    for (const std::string & log_path : log_paths) {
        std::ifstream log_file(log_path);
        if (!log_file.is_open()) continue;
        std::string log_line;
        while (getline(log_file, log_line)) {
            size_t idx = log_line.find("run.run_unixtime_ms=");
            if (idx == std::string::npos) continue;
            unixtime = std::stoi(log_line.substr(20));
            if (unixtime > run_start) run_start = unixtime;
        }
    }
    return unixtime;
}

double end_trip_time(int run_number) {
    string filename = (run_number < 1799) ? 
        "/sps/nemo/snemo/snemo_data/reco_data/UDD_betabeta_v1.list" :
        "/sps/nemo/snemo/snemo_data/reco_data/UDD_betabeta_v2.update";
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        int run; double start, dur, stop; std::string comment;
        if (!(iss >> run >> start >> dur >> stop >> comment)) continue;
        if (run == run_number) return static_cast<int>(stop);
    }
    return -1;
}

// Load per-OM energy thresholds from a text file.
// File format: one line per OM → "om_id  threshold_MeV"
// Returns a vector indexed by OM number (size = max_om + 1).
// If the file cannot be opened, fills with fallback_threshold.

inline std::vector<double> load_om_thresholds(const std::string& filepath,double fallback_threshold = 0.350){
    std::vector<double> thresholds(712, fallback_threshold);
    std::ifstream f(filepath);
    if (!f.is_open()) {
        std::cerr << "[WARNING] Could not open threshold file: "
                  << filepath
                  << " — using fallback " << fallback_threshold
                  << " MeV\n";
        return thresholds;
    }
    std::string line;
    // Ignore l'en-tête
    std::getline(f, line);
    while (std::getline(f, line)) {
        if (line.empty())
            continue;
        std::stringstream ss(line);
        std::string om_str, charge_str, energy_str;
        std::getline(ss, om_str, ';');
        std::getline(ss, charge_str, ';');
        std::getline(ss, energy_str, ';');
        if (om_str.empty())
            continue;
        int om_id = std::stoi(om_str);
        // Cas des lignes du type "0;;" ou "12;;"
        if (energy_str.empty())
            continue;
        double thr = std::stod(energy_str);
        if (om_id >= static_cast<int>(thresholds.size()))
            thresholds.resize(om_id + 1, fallback_threshold);
        thresholds[om_id] = thr;
    }
    return thresholds;
}

// ============================================================================
// Main Module Class
// ============================================================================
class falaise_skeleton_module_ptd : public dpp::chain_module {
public:
    falaise_skeleton_module_ptd();
    virtual ~falaise_skeleton_module_ptd();
    virtual void initialize(const datatools::properties&,
                           datatools::service_manager&,
                           dpp::module_handle_dict_type&);
    dpp::chain_module::process_status process(datatools::things& event);

private:
    // ========== CONFIGURATION & STATE ==========
    TFile *save_file;
    TTree *tree;
    int ptd_event_counter;
    std::int64_t event_number;
    bool ptd_details, ttd_details, sd_calo_details, sd_tracker_details;

    // ========== EVENT VARIABLES - Scalars ==========
    int cellules_non_associated, cellules_SD_non_associated, number_of_kinks;
    int run_number, trigger_id, nb_unfitted_cells, _simulation_phase_;
    int nb_gamma, nb_elec_ptd_per_event, nb_elec_SD_per_event, nb_wire_hit;
    
    bool has_an_electron_and_positron, has_SD_electron_and_positron;
    bool has_SD_two_electrons, opposite_side_e_gamma, same_side_elec;
    bool same_cluster_elec, has_a_same, hit_the_same_calo_hit;
    bool two_elec_more_350_keV, OM_are_neighbourg, has_kinks;

    double diff_time_elec, angle_3D_between_ep_em, delta_y_elec, delta_z_elec;
    double energy_elec_sum, corrected_energy_elec_sum, time_of_flight_gamma;
    double internal_theoretical_time_diff, external_theoretical_time_diff;
    double t1_th, t2_th, start_run_time, end_run_time, delta_r_calo;
    double one_kink_x, one_kink_y, one_kink_z, total_volume;
    double unix_start_time, first_time, last_time;
    double closest_gamma, closest_elec, closest_track, alpha_elec_time_diff;
    double energy_elec_1, energy_elec_2, current_timestamp;
    double vertex_3D_projected_x, vertex_3D_projected_y, vertex_3D_projected_z;
    double real_vertex_SD_x, real_vertex_SD_y, real_vertex_SD_z;
    double closest_time_track;

    // ========== EVENT VARIABLES - Vectors ==========
    vector<string> type_elec, g4_process, material, vertex_type, g4_material;
    vector<int> gamma_type, num_om, num_om_elec, track_number, num_gg;
    vector<int> num_om_elec_f, num_om_track, indices_calo_final;
    vector<double> time_gamma, time_gamma_before, time_gamma_after, angle_SD;
    vector<double> volume, total_step_length, energy_gamma, energy_elec, num_om_gamma;
    vector<double> corrected_energy_elec, time_elec, time_elec_alone, individual_energy_elec;
    vector<double> track_lenght, energy_gamma_after, kink_angle;
    vector<double> vertex_SD_x, vertex_SD_y, vertex_SD_z;
    vector<double> time_track, energy_track, chi2_track;
  vector<double> timestamp_elec_alone;
    vector<double> vertex_3D_start_x, vertex_3D_start_y, vertex_3D_start_z, track_length_vector;
    vector<double> vertex_3D_end_x, vertex_3D_end_y, vertex_3D_end_z;
    vector<double> vertex_gamma, kink_x, kink_y, kink_z;
    vector<double> vertex_3D_track_y, vertex_3D_track_z;
    vector<double> mean_alpha_anodic_time, mean_alpha_anodic_timestamp;
    vector<double> vertex_ttd_start_x, vertex_ttd_start_y, vertex_ttd_start_z;
    vector<double> vertex_ttd_end_x, vertex_ttd_end_y, vertex_ttd_end_z;
    vector<int> side_elec, cluster_elec_num, calo_num, nb_kink_per_track;

    // ========== LOOKUP TABLES ==========
    std::vector<double> y_source_pos;
    std::vector<double> z_source_pos;
    std::vector<double> source_pos_num;

    // ========== PER-PHASE, PER-OM ENERGY THRESHOLDS ==========
    // Loaded once at initialize() from energy_threshold_phase_N.txt
    // Indexed by OM number; one vector per phase (0-3).
    std::string threshold_dir;                        // directory containing the files
    std::vector<std::vector<double>> om_thresholds;   // om_thresholds[phase][om_id]

    // ========== PRIVATE METHODS ==========
    void initSourcePositions();
    void resetEventVariables();
    void processSimulatedData(datatools::things& event);
    bool processEventHeader(datatools::things& event);
    void processTrackerData(datatools::things& event);
    void processGammas(const snemo::datamodel::particle_track_data& PTD);
    void processElectrons(const snemo::datamodel::particle_track_data& PTD,
                         const snemo::datamodel::precalibrated_data& pCD,
                         std::vector<std::vector<double>>& kinks_vec,
                         std::vector<double>& x_alpha,
                         std::vector<double>& y_alpha,
                         std::vector<double>& z_alpha,
                         std::vector<int>& indices_for_pairs,
                         std::vector<int>& cluster_id);
  double compute_ellipse(double y_vertex, double z_vertex, int& source_num, double &dy, double &dz) const;    
    void processElectronVertex(
        const datatools::handle<snemo::datamodel::particle_track>& particle,
        const snemo::datamodel::precalibrated_data& pCD,
        std::vector<std::vector<double>>& kinks_vec,
        std::vector<double>& x_alpha, std::vector<double>& y_alpha, std::vector<double>& z_alpha,
        int& indice_elec, std::vector<int>& indices_for_pairs, std::vector<int>& cluster_id_tot,
        bool& found_good_electron);
    
    // Specialized helpers for electron vertex processing
    void extractElectronKinks(const datatools::handle<snemo::datamodel::particle_track>& particle);
    void extractElectronVertices(const datatools::handle<snemo::datamodel::particle_track>& particle,
                                bool& vertex_close_to_source, bool& vertex_associated_to_calo,
                                double& x_foil, double& y_foil, double& z_foil,
                                double& x_calo, double& y_calo, double& z_calo);
    void processGoodElectron(
        const datatools::handle<snemo::datamodel::particle_track>& particle,
        const snemo::datamodel::precalibrated_data& pCD,
        double x_foil, double y_foil, double z_foil,
        double x_calo, double y_calo, double z_calo,
        std::vector<std::vector<double>>& kinks_vec,
        int& indice_elec, std::vector<int>& indices_for_pairs);
    void processTrackOnlyElectron(
        const datatools::handle<snemo::datamodel::particle_track>& particle);
    void processAlphaTrack(
        const datatools::handle<snemo::datamodel::particle_track>& particle,
        const snemo::datamodel::precalibrated_data& pCD,
        bool vertex_close_to_source,
        double x_foil, double y_foil, double z_foil,
        std::vector<double>& x_alpha, std::vector<double>& y_alpha, std::vector<double>& z_alpha);
    void computeKinkAngle(double x2, double y2, double z2,
                         double x_foil, double y_foil, double z_foil,
                         double x_calo, double y_calo, double z_calo,
                         int nb_kinks_elec);
    
    void processPairs(const std::vector<int>& indices_for_pairs,
                     const std::vector<int>& cluster_id,
                     std::vector<std::vector<double>>& kinks_vec,
                     const std::vector<double>& x_alpha,
                     const std::vector<double>& y_alpha,
                     const std::vector<double>& z_alpha);
    
    void processSinglePair(int best_a, int best_b,
                          std::vector<std::vector<double>>& kinks_vec,
                          const std::vector<double>& x_alpha,
                          const std::vector<double>& y_alpha,
                          const std::vector<double>& z_alpha,
                          double min_dt_pair);
    
    void analyzeGammasForPair(int best_a, int best_b,
                             double tA, double tB,
                             double timestampA, double timestampB);
    
    void analyzeAlphasForPair(const double* start_p1, const double* start_p2,
                             double timestampA, double timestampB,
                             const std::vector<double>& x_alpha,
                             const std::vector<double>& y_alpha,
                             const std::vector<double>& z_alpha);

    void computePairGeometry(int best_a, int best_b,
                            const std::vector<std::vector<double>>& kinks_vec);
    
    void computePairTiming(int best_a, int best_b,
                          const double* start_p1, const double* end_p1,
                          const double* start_p2, const double* end_p2);
    
    double extract_unix_start_time(int run_number) const;
    int get_phase(int run_number);
  
    DPP_MODULE_REGISTRATION_INTERFACE(falaise_skeleton_module_ptd);
};

// ============================================================================
// Module Registration
// ============================================================================
DPP_MODULE_REGISTRATION_IMPLEMENT(falaise_skeleton_module_ptd, "FalaiseSkeletonModule_PTD")

// ============================================================================
// Constructor
// ============================================================================
falaise_skeleton_module_ptd::falaise_skeleton_module_ptd() : dpp::chain_module() {
    save_file = new TFile("extracted_data.root", "RECREATE");
    tree = new TTree("Event", "Event information");
    initSourcePositions();

    // Register all ROOT branches
    tree->Branch("run_number", &run_number);
    tree->Branch("event_number", &event_number);
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
    tree->Branch("same_cluster_elec", &same_cluster_elec);
    tree->Branch("cluster_elec_num", &cluster_elec_num);
    tree->Branch("diff_time_elec", &diff_time_elec);
    tree->Branch("energy_elec_sum", &energy_elec_sum);
    tree->Branch("energy_elec_1", &energy_elec_1);
    tree->Branch("energy_elec_2", &energy_elec_2);
    tree->Branch("individual_energy_elec", &individual_energy_elec);
    tree->Branch("delta_y_elec", &delta_y_elec);
    tree->Branch("delta_z_elec", &delta_z_elec);
    tree->Branch("angle_3D_between_ep_em", &angle_3D_between_ep_em);
    tree->Branch("has_a_same", &has_a_same);
    tree->Branch("hit_the_same_calo_hit", &hit_the_same_calo_hit);
    tree->Branch("has_an_electron_and_positron", &has_an_electron_and_positron);
    tree->Branch("opposite_side_e_gamma", &opposite_side_e_gamma);
    tree->Branch("cellules_non_associated", &cellules_non_associated);
    tree->Branch("time_of_flight_gamma", &time_of_flight_gamma);
    tree->Branch("internal_theoretical_time_diff", &internal_theoretical_time_diff);
    tree->Branch("external_theoretical_time_diff", &external_theoretical_time_diff);
    tree->Branch("t1_th", &t1_th);
    tree->Branch("t2_th", &t2_th);
    tree->Branch("number_of_kinks", &number_of_kinks);
    tree->Branch("kink_x", &kink_x);
    tree->Branch("kink_y", &kink_y);
    tree->Branch("kink_z", &kink_z);
    tree->Branch("kink_angle", &kink_angle);
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
    tree->Branch("current_timestamp", &current_timestamp);
}

// ============================================================================
// Destructor
// ============================================================================
falaise_skeleton_module_ptd::~falaise_skeleton_module_ptd() {
    double time = last_time - first_time;
    std::int64_t last_event_number = event_number;
    save_file->cd();
    TParameter<double> param("run_time", time);
    param.Write();
    TParameter<Long64_t> param_num("last_event_number", last_event_number);
    param_num.Write();
    tree->Write();
    save_file->Close();
    delete save_file;
    std::cout << "Module cleaned up" << std::endl;
}

// ============================================================================
// Initialize source positions
// ============================================================================
void falaise_skeleton_module_ptd::initSourcePositions() {
    y_source_pos = {
        -2087.5, -2087.5, -2087.5, -2087.5, -2087.5, -2087.5, -2087.5,
        -1252.5, -1252.5, -1252.5, -1252.5, -1252.5, -1252.5, -1252.5,
         -417.5,  -417.5,  -417.5,  -417.5,  -417.5,  -417.5,  -417.5,
         417.5,   417.5,   417.5,   417.5,   417.5,   417.5,   417.5,
        1252.5,  1252.5,  1252.5,  1252.5,  1252.5,  1252.5,  1252.5,
        2087.5,  2087.5,  2087.5,  2087.5,  2087.5,  2087.5,  2087.5
    };
    z_source_pos = {
        1317.5,   882.5,   447.5,     2.5,  -442.5,  -882.5, -1317.5,
        1317.5,   887.5,   447.5,     2.5,  -442.5,  -887.5, -1317.5,
        1312.5,   882.5,   447.5,     2.5,  -442.5,  -882.5, -1317.5,
        1317.5,   887.5,   447.5,     2.5,  -437.5,  -882.5, -1317.5,
        1317.5,   882.5,   442.5,     2.5,  -442.5,  -882.5, -1317.5,
        1317.5,   882.5,   442.5,     2.5,  -442.5,  -882.5, -1317.5
    };
    source_pos_num = {
        6,  5,  4,  3,   2,  1,  0,
        13, 12, 11, 10,  9,  8,  7,
        20, 19, 18, 17, 16, 15, 14,
        27, 26, 25, 24, 23, 22, 21,
        34, 33, 32, 31, 30, 29, 28,
        41, 40, 39, 38, 37, 36, 35
    };
}



double falaise_skeleton_module_ptd::compute_ellipse(double y_vertex, double z_vertex, int& source_num, double &best_dy, double &best_dz) const
{
  best_dy = 0.0;
  best_dz = 0.0;

  double min_dist2 = std::numeric_limits<double>::max();

  for (size_t i = 0; i < y_source_pos.size(); ++i) {
    double dy = y_source_pos[i] - y_vertex;
    double dz = z_source_pos[i] - z_vertex;
    double dist2 = (dy * dy) / (25.0 * 25.0) + (dz * dz) / (30.0 * 30.0);

    if (dist2 < min_dist2) {
      min_dist2 = dist2;
      source_num = source_pos_num[i];
      best_dy = dy;
      best_dz = dz;
    }
  }

  return min_dist2; // distance minimale                                                                 
}



// ============================================================================
// Reset all event variables
// ============================================================================
void falaise_skeleton_module_ptd::resetEventVariables() {
    nb_gamma = 0;
    nb_elec_ptd_per_event = 0;
    nb_elec_SD_per_event = 0;
    nb_wire_hit = 0;
    cellules_non_associated = 0;
    cellules_SD_non_associated = 0;
    number_of_kinks = 0;
    nb_unfitted_cells = 0;
    gamma_type.clear();
    time_gamma.clear();
    time_gamma_before.clear();
    time_gamma_after.clear();
    energy_gamma_after.clear();
    energy_gamma.clear();
    type_elec.clear();
    side_elec.clear();
    cluster_elec_num.clear();
    energy_elec.clear();
    individual_energy_elec.clear();
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
    mean_alpha_anodic_time.clear();
    mean_alpha_anodic_timestamp.clear();
    vertex_3D_start_x.clear();
    vertex_3D_start_y.clear();
    vertex_3D_start_z.clear();
    track_length_vector.clear();
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
    chi2_track.clear();
    nb_kink_per_track.clear();
    
    has_an_electron_and_positron = false;
    has_SD_electron_and_positron = false;
    has_SD_two_electrons = false;
    opposite_side_e_gamma = false;
    same_side_elec = false;
    same_cluster_elec = false;
    has_a_same = false;
    hit_the_same_calo_hit = false;
    two_elec_more_350_keV = false;
    OM_are_neighbourg = false;
    has_kinks = false;

    diff_time_elec = 0;
    angle_3D_between_ep_em = 0.0;
    delta_y_elec = 0;
    delta_z_elec = 0;
    energy_elec_sum = 0;
    corrected_energy_elec_sum = 0;
    time_of_flight_gamma = 0.0;
    internal_theoretical_time_diff = 0.0;
    external_theoretical_time_diff = 0.0;
    t1_th = 0.0;
    t2_th = 0.0;
    start_run_time = 0.0;
    end_run_time = 0.0;
    delta_r_calo = 0.0;
    kink_angle.clear();
    one_kink_x = 0.0;
    one_kink_y = 0.0;
    one_kink_z = 0.0;
    total_volume = 0.0;
    run_number = 0;
    trigger_id = 0;
    closest_gamma = 1e6;
    closest_elec = 1e6;
    closest_time_track = 1e6;
    closest_track = 1e6;
    energy_elec_1 = 100;
    energy_elec_2 = 100;
    vertex_3D_projected_x = 0;
    vertex_3D_projected_y = 0;
    vertex_3D_projected_z = 0;
    real_vertex_SD_x = 0.0;
    real_vertex_SD_y = 0.0;
    real_vertex_SD_z = 0.0;
    alpha_elec_time_diff = 0;
}

// ============================================================================
// Initialize module
// ============================================================================
void falaise_skeleton_module_ptd::initialize(const datatools::properties& module_properties,
                                             datatools::service_manager&,
                                             dpp::module_handle_dict_type&) {
    std::cout << "Initializing FalaiseSkeletonModule_PTD" << std::endl;
    event_number = 0;
    ptd_event_counter = 0;
    first_time = 0;
    last_time = 0;
    unix_start_time = 0;

    ptd_details = module_properties.has_key("ptd_details") ?
                  module_properties.fetch_boolean("ptd_details") : false;
    ttd_details = module_properties.has_key("ttd_details") ?
                  module_properties.fetch_boolean("ttd_details") : false;
    sd_calo_details = module_properties.has_key("calo_details") ?
                      module_properties.fetch_boolean("calo_details") : false;
    sd_tracker_details = module_properties.has_key("tracker_details") ?
      module_properties.fetch_boolean("tracker_details") : false;

    if (module_properties.has_key("simulation_phase")) {
      _simulation_phase_ = module_properties.fetch_integer("simulation_phase");
    }
    // Load per-phase, per-OM energy thresholds
    threshold_dir = "/sps/nemo/scratch/margnes/data/seuil_max";
    om_thresholds.resize(5);
    for (int phase = 0; phase <= 4; ++phase) {
      double threshold_value = 0.35;
      if(phase != 0){
	threshold_value=0.05;
      }
        std::string fname = threshold_dir + "/phase" + std::to_string(phase) + ".txt";
        om_thresholds[phase] = load_om_thresholds(fname,threshold_value);
        std::cout << "Loaded thresholds for phase " << phase
                  << " from " << fname
                  << " (" << om_thresholds[phase].size() << " OMs)\n";
    }

    this->_set_initialized(true);
}

// ============================================================================
// Process Simulated Data (SD)
// ============================================================================
void falaise_skeleton_module_ptd::processSimulatedData(datatools::things& event) {
    if (!event.has("SD")) return;

    const mctools::simulated_data& SD = event.get<mctools::simulated_data>("SD");
    const mctools::simulated_data::primary_event_type& pevent = SD.get_primary_event();
    const genbb::primary_event::particles_col_type& particles = pevent.get_particles();

    if (particles.size() == 2) {
        if (particles.front().get_particle_label() == "e-" &&
            particles.back().get_particle_label() == "e-") {
            has_SD_two_electrons = true;
        }
        if ((particles.front().get_particle_label() == "e+" &&
             particles.back().get_particle_label() == "e-") ||
            (particles.front().get_particle_label() == "e-" &&
             particles.back().get_particle_label() == "e+")) {
            has_SD_electron_and_positron = true;
        }
    }
    real_vertex_SD_x = SD.get_vertex()[0] / CLHEP::mm;
    real_vertex_SD_y = SD.get_vertex()[1] / CLHEP::mm;
    real_vertex_SD_z = SD.get_vertex()[2] / CLHEP::mm;
}

// ============================================================================
// Process Event Header
// ============================================================================
bool falaise_skeleton_module_ptd::processEventHeader(datatools::things& event) {
    const snemo::datamodel::event_header& eh = event.get<snemo::datamodel::event_header>("EH");
    run_number = eh.get_id().get_run_number();

    if (event.has("SD") == 0) {
        const snemo::datamodel::unified_digitized_data& UDD =
            event.get<snemo::datamodel::unified_digitized_data>("UDD");
        const auto& trigger_ids = UDD.get_origin_trigger_ids();
        trigger_id = -1;
        if (!trigger_ids.empty()) {
            trigger_id = *trigger_ids.begin();
        }
    }

    // Charger unix_start_time AVANT le test temporel
    if (unix_start_time == 0) {
        if (run_number < 1556) {
            unix_start_time = snemo_run_time(run_number);
        } else {
            unix_start_time = extract_unix_start_time(run_number);
        }
    }

    if (eh.has_timestamp()) {
        double time = eh.get_timestamp().get_seconds() +
                      eh.get_timestamp().get_picoseconds() / 1E12;
	current_timestamp = time;
        double trip_end = end_trip_time(run_number);
        //cout << event_number << " time=" << time
	//   << " limit=" << unix_start_time + trip_end
	//   << " unix_start=" << unix_start_time
	//   << " trip_end=" << trip_end
	//   <<" diff "<<time - (unix_start_time + trip_end)
	//   << endl;
        if (time > unix_start_time + trip_end && trip_end > 0) {
	  //cout<<"SHOULD BREAK "<<endl;
            event_number++;
	    return false;
        }
        last_time = time;
        if (first_time == 0) first_time = time;
    }
    return true;
}

// ============================================================================
// Process Tracker Data
// ============================================================================
void falaise_skeleton_module_ptd::processTrackerData(datatools::things& event) {
    const snemo::datamodel::tracker_trajectory_data& TTD =
        event.get<snemo::datamodel::tracker_trajectory_data>("TTD");
    const snemo::datamodel::tracker_trajectory_solution& ttd_solution =
        TTD.get_default_solution();

    for (const auto& trajectory : ttd_solution.get_trajectories()) {
        if (trajectory->get_fit_infos().has_chi2() &&
            trajectory->get_fit_infos().has_ndof()) {
            chi2_track.push_back(trajectory->get_fit_infos().get_chi2() /
                                trajectory->get_fit_infos().get_ndof());
        }
    }

    const snemo::datamodel::tracker_clustering_data& TCD =
        event.get<snemo::datamodel::tracker_clustering_data>("TCD");
    const snemo::datamodel::tracker_clustering_solution& tcd_solution = TCD.get_default();
    nb_unfitted_cells = tcd_solution.get_unclustered_hits().size();

    if (ttd_details && ttd_solution.has_unfitted_clusters()) {
        cellules_non_associated++;
    }
}

// ============================================================================
// Process Gammas
// ============================================================================
void falaise_skeleton_module_ptd::processGammas(const snemo::datamodel::particle_track_data& PTD) {
    if (!PTD.hasIsolatedCalorimeters()) return;
    const snemo::datamodel::CalorimeterHitHdlCollection& cc_collection = PTD.isolatedCalorimeters();
    for (const auto& it_hit : cc_collection) {
      double energy_threshold = 0.0; // fallback                                                
      int om_id = snemo::datamodel::om_num(it_hit->get_geom_id());
      int phase  = get_phase(run_number);
      if (phase >= 0 && phase <= 4 && om_id >= 0) {
	// if threshold is wanted :                                                        
	energy_threshold = om_thresholds[phase][om_id];                                  
	//cout<<event_number<<" "<<energy_threshold<<" "<<phase<<" "<<om_id<<endl;                                                              
	// if no threshold is wanted :                                                     
	//if(phase==0){
	//energy_threshold = 0.35;
	//}
	//else{
	//energy_threshold = 0.05;
	//}
      }
      //cout<<phase<<" "<<energy_threshold<<endl;
      if(it_hit->get_energy()<energy_threshold) continue;
        energy_gamma.push_back(it_hit->get_energy());
        gamma_type.push_back(it_hit->get_geom_id().get(0));
        vertex_gamma.push_back(it_hit->get_geom_id().get(1));
        time_gamma.push_back(it_hit->get_time() / CLHEP::ns);
        nb_gamma++;
        num_om_gamma.push_back(snemo::datamodel::om_num(it_hit->get_geom_id()));
    }
}

// ============================================================================
// Process Electrons 
// ============================================================================
void falaise_skeleton_module_ptd::processElectrons(
    const snemo::datamodel::particle_track_data& PTD,
    const snemo::datamodel::precalibrated_data& pCD,
    std::vector<std::vector<double>>& kinks_vec,
    std::vector<double>& x_alpha, std::vector<double>& y_alpha, std::vector<double>& z_alpha,
    std::vector<int>& indices_for_pairs, std::vector<int>& cluster_id) {

    int indice_elec = 0;
    std::vector<int> cluster_id_tot;

    for (const datatools::handle<snemo::datamodel::particle_track>& particle : PTD.particles()) {
        bool found_good_electron = false;
        
        processElectronVertex(particle, pCD, kinks_vec, x_alpha, y_alpha, z_alpha,
                             indice_elec, indices_for_pairs, cluster_id_tot, found_good_electron);

        if (!found_good_electron) {
            cellules_non_associated++;
        }
    }
}

// ============================================================================
// Process Single Electron Vertex - Dispatches to specialized handlers
// ============================================================================
void falaise_skeleton_module_ptd::processElectronVertex(
    const datatools::handle<snemo::datamodel::particle_track>& particle,
    const snemo::datamodel::precalibrated_data& pCD,
    std::vector<std::vector<double>>& kinks_vec,
    std::vector<double>& x_alpha, std::vector<double>& y_alpha, std::vector<double>& z_alpha,
    int& indice_elec, std::vector<int>& indices_for_pairs, std::vector<int>& cluster_id_tot,
    bool& found_good_electron) {

    // Step 1: Extract kinks from trajectory
    extractElectronKinks(particle);

    // Step 2: Check calorimeter association
    const auto& calorimeter_hits = particle->get_associated_calorimeter_hits();
    if (calorimeter_hits.size() > 1) return;

    // Step 3: Extract vertex positions
    bool vertex_close_to_source = false;
    bool vertex_associated_to_calo = false;
    double x_foil = 0, y_foil = 0, z_foil = 0;
    double x_calo = 0, y_calo = 0, z_calo = 0;
    
    extractElectronVertices(particle, vertex_close_to_source, vertex_associated_to_calo,
                           x_foil, y_foil, z_foil, x_calo, y_calo, z_calo);

    // Step 4: Get the OM-specific energy threshold for the current run/phase
    double energy_threshold = 0.0; // fallback
    if (calorimeter_hits.size() == 1) {
        int om_id = snemo::datamodel::om_num(calorimeter_hits[0]->get_geom_id());
        int phase  = get_phase(run_number);
        if (phase >= 0 && phase <= 4 && om_id >= 0) {
	  // if threshold is wanted :
	  energy_threshold = om_thresholds[phase][om_id];
	  //cout<<event_number<<" "<<energy_threshold<<" "<<phase<<" "<<om_id<<" "<<calorimeter_hits[0]->get_energy()<<endl;
	  // if no threshold is wanted :
	  //if(phase==0){
	  //  energy_threshold = 0.35;
	  //}
	  //else{
	  //  energy_threshold = 0.05;
	  //}
	}
	//cout<<phase<<" "<<energy_threshold<<endl;
    }
    
    // Step 5: Dispatch to specialized handler based on vertex type
    if (vertex_close_to_source && vertex_associated_to_calo &&
        calorimeter_hits.size() == 1 && calorimeter_hits[0]->get_energy() > energy_threshold) {
        // Case 1: Good electron (source + calo)
        found_good_electron = true;
        processGoodElectron(particle, pCD, x_foil, y_foil, z_foil, x_calo, y_calo, z_calo,
                           kinks_vec, indice_elec, indices_for_pairs);

    } else if (vertex_associated_to_calo && !vertex_close_to_source &&
               calorimeter_hits.size() == 1 && calorimeter_hits[0]->get_energy() > energy_threshold) {
        // Case 2: Track only (no source vertex)
        processTrackOnlyElectron(particle);

    } else if (calorimeter_hits.size() == 0 || calorimeter_hits[0]->get_energy() < energy_threshold) {
        // Case 3: Alpha track (no calo hit)
        processAlphaTrack(particle, pCD, vertex_close_to_source,
                         x_foil, y_foil, z_foil, x_alpha, y_alpha, z_alpha);
    }
}

// ============================================================================
// Extract electron kinks from trajectory
// ============================================================================
void falaise_skeleton_module_ptd::extractElectronKinks(
    const datatools::handle<snemo::datamodel::particle_track>& particle) {
    
    const auto& trajectory_pattern = particle->get_trajectory_handle()->get_pattern();
    int nb_kinks = trajectory_pattern.number_of_kinks();
    number_of_kinks = nb_kinks;

    if (nb_kinks > 0) {
        for (unsigned int i = 0; i < nb_kinks; ++i) {
            kink_x.push_back(trajectory_pattern.get_kink(i).getX());
            kink_y.push_back(trajectory_pattern.get_kink(i).getY());
            kink_z.push_back(trajectory_pattern.get_kink(i).getZ());
        }
    }
}

// ============================================================================
// Extract electron vertex positions (foil and calorimeter)
// ============================================================================
void falaise_skeleton_module_ptd::extractElectronVertices(
    const datatools::handle<snemo::datamodel::particle_track>& particle,
    bool& vertex_close_to_source, bool& vertex_associated_to_calo,
    double& x_foil, double& y_foil, double& z_foil,
    double& x_calo, double& y_calo, double& z_calo) {

    vertex_close_to_source = false;
    vertex_associated_to_calo = false;

    for (const datatools::handle<snemo::datamodel::vertex>& vertex :
         particle->get_vertices()) {
        if (vertex->is_on_reference_source_plane()) {
            x_foil = vertex->get_spot().get_position().getX();
            y_foil = vertex->get_spot().get_position().getY();
            z_foil = vertex->get_spot().get_position().getZ();
            vertex_close_to_source = true;
        } else if (vertex->is_on_source_foil() && !vertex_close_to_source) {
            x_foil = vertex->get_spot().get_position().getX();
            y_foil = vertex->get_spot().get_position().getY();
            z_foil = vertex->get_spot().get_position().getZ();
            vertex_close_to_source = true;
        } else if (vertex->is_on_main_calorimeter()) {
            x_calo = vertex->get_spot().get_position().getX();
            y_calo = vertex->get_spot().get_position().getY();
            z_calo = vertex->get_spot().get_position().getZ();
            vertex_associated_to_calo = true;
        }
    }
}

// ============================================================================
// Compute kink angle at the kink point
// ============================================================================
void falaise_skeleton_module_ptd::computeKinkAngle(
    double x2, double y2, double z2,
    double x_foil, double y_foil, double z_foil,
    double x_calo, double y_calo, double z_calo,
    int nb_kinks_elec) {
    
    for (int i = 0; i < nb_kinks_elec; ++i) {
        double dot_prod = (x2 - x_foil) * (x_calo - x2) +
                         (y2 - y_foil) * (y_calo - y2) +
                         (z2 - z_foil) * (z_calo - z2);
        double norm_u = sqrt((x2 - x_foil) * (x2 - x_foil) +
                            (y2 - y_foil) * (y2 - y_foil) +
                            (z2 - z_foil) * (z2 - z_foil));
        double norm_v = sqrt((x_calo - x2) * (x_calo - x2) +
                            (y_calo - y2) * (y_calo - y2) +
                            (z_calo - z2) * (z_calo - z2));
        
        if (norm_u > 1e-9 && norm_v > 1e-9) {
            kink_angle.push_back((180.0 / M_PI) * acos(dot_prod / (norm_u * norm_v)));
        }
    }
}


// ============================================================================
// Process good electron (with source and calorimeter vertices)
// ============================================================================
void falaise_skeleton_module_ptd::processGoodElectron(
    const datatools::handle<snemo::datamodel::particle_track>& particle,
    const snemo::datamodel::precalibrated_data& pCD,
    double x_foil, double y_foil, double z_foil,
    double x_calo, double y_calo, double z_calo,
    std::vector<std::vector<double>>& kinks_vec,
    int& indice_elec, std::vector<int>& indices_for_pairs) {

    int source_num;
    const auto& traj_pattern_elec = particle->get_trajectory_handle()->get_pattern();
    const auto& calorimeter_hits = particle->get_associated_calorimeter_hits();

    int nb_kinks_elec = traj_pattern_elec.number_of_kinks();
    if (nb_kinks_elec > 0) {
      for (int i = 0; i < nb_kinks_elec; ++i) {
        double x2 = traj_pattern_elec.get_kink(i).getX();
        double y2 = traj_pattern_elec.get_kink(i).getY();
        double z2 = traj_pattern_elec.get_kink(i).getZ();
        computeKinkAngle(x2, y2, z2, x_foil, y_foil, z_foil, x_calo, y_calo, z_calo, nb_kinks_elec);
      }
    }
    track_length_vector.push_back(particle->get_trajectory().get_pattern().get_shape().get_length());
    vertex_3D_start_x.push_back(x_foil);
    vertex_3D_start_y.push_back(y_foil);
    vertex_3D_start_z.push_back(z_foil);
    vertex_ttd_start_x.push_back(traj_pattern_elec.get_first()[0]);
    vertex_ttd_start_y.push_back(traj_pattern_elec.get_first()[1]);
    vertex_ttd_start_z.push_back(traj_pattern_elec.get_first()[2]);
    vertex_3D_end_x.push_back(x_calo);
    vertex_3D_end_y.push_back(y_calo);
    vertex_3D_end_z.push_back(z_calo);
    vertex_ttd_end_x.push_back(traj_pattern_elec.get_last()[0]);
    vertex_ttd_end_y.push_back(traj_pattern_elec.get_last()[1]);
    vertex_ttd_end_z.push_back(traj_pattern_elec.get_last()[2]);

    // Store kinks vector (with index alignment guarantee)
    if (number_of_kinks > 0) {
        std::vector<double> kink_point = {
            traj_pattern_elec.get_kink(0).getX(),
            traj_pattern_elec.get_kink(0).getY(),
            traj_pattern_elec.get_kink(0).getZ()
        };
        kinks_vec.push_back(kink_point);
    } else {
        kinks_vec.push_back(std::vector<double>());
    }

    // Store calorimeter information
    int number = snemo::datamodel::om_num(calorimeter_hits[0]->get_geom_id());
    const auto pCD_calo = pCD.calorimeter_hits();
    for(const auto & pCD_calo_hit : pCD_calo){
      if(number==snemo::datamodel::om_num(pCD_calo_hit->get_geom_id())){
	calo_num.push_back(snemo::datamodel::om_num(pCD_calo_hit->get_geom_id()));
	timestamp_elec_alone.push_back(pCD_calo_hit->get_time() / CLHEP::ns);
	break;
      }
    }
    energy_elec.push_back(calorimeter_hits[0]->get_energy());
    time_elec.push_back(calorimeter_hits[0]->get_time() / CLHEP::ns);
    time_elec_alone.push_back(calorimeter_hits[0]->get_time() / CLHEP::ns);
    side_elec.push_back(calorimeter_hits[0]->get_geom_id().get(1));
    num_om_elec.push_back(snemo::datamodel::om_num(calorimeter_hits[0]->get_geom_id()));
    num_om.push_back(snemo::datamodel::om_num(calorimeter_hits[0]->get_geom_id()));

    // Store cluster information
    if (particle->has_trajectory()) {
        const auto& trajectory = particle->get_trajectory();
        if (trajectory.has_cluster()) {
            cluster_elec_num.push_back(trajectory.get_cluster().get_hit_id());
        }
    }

    // Register this as a good electron
    nb_kink_per_track.push_back(nb_kinks_elec);
    nb_elec_ptd_per_event++;
    indices_for_pairs.push_back(indice_elec);
    indice_elec++;
}

// ============================================================================
// Process track-only electron (no source vertex)
// ============================================================================
void falaise_skeleton_module_ptd::processTrackOnlyElectron(
    const datatools::handle<snemo::datamodel::particle_track>& particle) {

    const auto& calorimeter_hits = particle->get_associated_calorimeter_hits();
    
    num_om.push_back(snemo::datamodel::om_num(calorimeter_hits[0]->get_geom_id()));
    num_om_track.push_back(snemo::datamodel::om_num(calorimeter_hits[0]->get_geom_id()));
    time_elec.push_back(calorimeter_hits[0]->get_time() / CLHEP::ns);
    time_track.push_back(calorimeter_hits[0]->get_time() / CLHEP::ns);
    energy_track.push_back(calorimeter_hits[0]->get_energy());
}

// ============================================================================
// Process alpha track (no good calorimeter hit)
// ============================================================================
void falaise_skeleton_module_ptd::processAlphaTrack(
    const datatools::handle<snemo::datamodel::particle_track>& particle,
    const snemo::datamodel::precalibrated_data& pCD,
    bool vertex_close_to_source,
    double x_foil, double y_foil, double z_foil,
    std::vector<double>& x_alpha, std::vector<double>& y_alpha, std::vector<double>& z_alpha) {

    // Collect geiger cell IDs for this alpha track
    std::vector<int> alpha_gg_num;
    const auto gg = particle->get_trajectory_handle()->get_cluster().hits();
    for (const auto hits : gg) {
        alpha_gg_num.push_back(snemo::datamodel::gg_num(hits->get_geom_id()));
    }

    // Compute mean anodic timing
    double sum_cell = 0;
    int n_cell = 0;
    const auto pCD_tracker_hits = pCD.tracker_hits();
    for (const auto pCD_track_hits : pCD_tracker_hits) {
        if (std::find(alpha_gg_num.begin(), alpha_gg_num.end(),
                     snemo::datamodel::gg_num(pCD_track_hits->get_geom_id())) !=
            alpha_gg_num.end()) {
            sum_cell += pCD_track_hits->get_anodic_time() / CLHEP::ns;
            n_cell++;
        }
    }
    double new_mean = (n_cell > 0 ? sum_cell / n_cell : 0);
    mean_alpha_anodic_timestamp.push_back(new_mean);

    // Store alpha vertex position
    if (vertex_close_to_source) {
        x_alpha.push_back(x_foil);
        y_alpha.push_back(y_foil);
        z_alpha.push_back(z_foil);
    } else {
        const auto& traj_pattern = particle->get_trajectory_handle()->get_pattern();
        x_alpha.push_back(traj_pattern.get_first()[0]);
        x_alpha.push_back(traj_pattern.get_last()[0]);
        y_alpha.push_back(traj_pattern.get_first()[1]);
        y_alpha.push_back(traj_pattern.get_last()[1]);
        z_alpha.push_back(traj_pattern.get_first()[2]);
        z_alpha.push_back(traj_pattern.get_last()[2]);
    }
}

// ============================================================================
// Process Pairs
// ============================================================================
void falaise_skeleton_module_ptd::processPairs(
    const std::vector<int>& indices_for_pairs,
    const std::vector<int>& cluster_id,
    std::vector<std::vector<double>>& kinks_vec,
    const std::vector<double>& x_alpha,
    const std::vector<double>& y_alpha,
    const std::vector<double>& z_alpha) {

    if (nb_elec_ptd_per_event < 2) return;

    // Find pairs by closest timing
    std::vector<std::vector<int>> pairs;
    std::vector<int> remaining = indices_for_pairs;

    while (remaining.size() >= 2) {
        double best_dist2 = std::numeric_limits<double>::max();
        int best_i = -1, best_j = -1;

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

        if (best_i < 0 || best_j < 0) break;
        pairs.push_back({remaining[best_i], remaining[best_j]});
        if (best_i > best_j) std::swap(best_i, best_j);
        remaining.erase(remaining.begin() + best_j);
        remaining.erase(remaining.begin() + best_i);
    }

    // Compute min timing to other electrons/tracks
    std::vector<double> min_dt_for_pair;
    min_dt_for_pair.resize(pairs.size(), 1e9);

    for (size_t p = 0; p < pairs.size(); ++p) {
        int best_a = pairs[p][0];
        int best_b = pairs[p][1];
        double tA = time_elec_alone[best_a];
        double tB = time_elec_alone[best_b];
        double min_dt = 1e9;

        for (size_t idx = 0; idx < cluster_id.size(); ++idx) {
            if (idx == (size_t)best_a || idx == (size_t)best_b) continue;
            double dtA = std::abs(time_elec_alone[idx] - tA);
            double dtB = std::abs(time_elec_alone[idx] - tB);
            min_dt = std::min(min_dt, std::min(dtA, dtB));
        }

        for (size_t idx = 0; idx < time_track.size(); ++idx) {
            double dtA = std::abs(time_track[idx] - tA);
            double dtB = std::abs(time_track[idx] - tB);
            min_dt = std::min(min_dt, std::min(dtA, dtB));
        }

        min_dt_for_pair[p] = min_dt;
    }

    // Process each pair
    for (size_t p = 0; p < pairs.size(); ++p) {
        processSinglePair(pairs[p][0], pairs[p][1], kinks_vec, x_alpha, y_alpha, z_alpha,
                         min_dt_for_pair[p]);
    }
}

// ============================================================================
// Process Single Pair
// ============================================================================
void falaise_skeleton_module_ptd::processSinglePair(
    int best_a, int best_b,
    //best a and best b are the indexes of double beta electrons candidates
    std::vector<std::vector<double>>& kinks_vec,
    const std::vector<double>& x_alpha,
    const std::vector<double>& y_alpha,
    const std::vector<double>& z_alpha,
    double min_dt_pair) {

    // Reset pair variables
    energy_elec_1 = 100;
    energy_elec_2 = 100;
    two_elec_more_350_keV = false;
    OM_are_neighbourg = false;
    has_kinks = false;

    // Bounds checking
    if (best_a >= (int)energy_elec.size() || best_b >= (int)energy_elec.size()) {
        std::cerr << "ERROR: processSinglePair - invalid electron indices!" << std::endl;
        std::cerr << "best_a=" << best_a << " best_b=" << best_b;
        std::cerr << " energy_elec.size()=" << energy_elec.size() << std::endl;
        return;
    }

    num_om_elec_f.push_back(num_om_elec[best_a]);
    num_om_elec_f.push_back(num_om_elec[best_b]);

    // Check neighbor and kink status
    if (areOMNeighbors(num_om_elec[best_a], num_om_elec[best_b])) {
        OM_are_neighbourg = true;
    }
    if (nb_kink_per_track[best_a] != 0 || nb_kink_per_track[best_b] != 0) {
        has_kinks = true;
    }

    // Geometry computation
    computePairGeometry(best_a, best_b, kinks_vec);

    // Energy and timing
    energy_elec_sum = energy_elec[best_a] + energy_elec[best_b];
    energy_elec_1 = energy_elec[best_a];
    energy_elec_2 = energy_elec[best_b];
    individual_energy_elec.push_back(energy_elec[best_a]);
    individual_energy_elec.push_back(energy_elec[best_b]);
    closest_elec = min_dt_pair;

    if (energy_elec[best_a] > 0.050 && energy_elec[best_b] > 0.050) {
        two_elec_more_350_keV = true;
    }

    // Compute angles and TOF
    double start_p1[3] = {vertex_3D_start_x[best_a], vertex_3D_start_y[best_a],
                         vertex_3D_start_z[best_a]};
    double start_p2[3] = {vertex_3D_start_x[best_b], vertex_3D_start_y[best_b],
                         vertex_3D_start_z[best_b]};
    double end_p1[3] = {vertex_3D_end_x[best_a], vertex_3D_end_y[best_a],
                       vertex_3D_end_z[best_a]};
    double end_p2[3] = {vertex_3D_end_x[best_b], vertex_3D_end_y[best_b],
                       vertex_3D_end_z[best_b]};

    if (best_a < (int)kinks_vec.size() && !kinks_vec[best_a].empty()) {
        end_p1[0] = kinks_vec[best_a][0];
        end_p1[1] = kinks_vec[best_a][1];
        end_p1[2] = kinks_vec[best_a][2];
    }
    if (best_b < (int)kinks_vec.size() && !kinks_vec[best_b].empty()) {
        end_p2[0] = kinks_vec[best_b][0];
        end_p2[1] = kinks_vec[best_b][1];
        end_p2[2] = kinks_vec[best_b][2];
    }

    angle_3D_between_ep_em = std::abs(computeAngle(start_p1, end_p1, start_p2, end_p2));
    diff_time_elec = time_elec_alone[best_a] - time_elec_alone[best_b];
    computePairTiming(best_a, best_b, start_p1, end_p1, start_p2, end_p2);

    // Gamma and alpha analysis
    analyzeGammasForPair(best_a, best_b, 
                        time_elec_alone[best_a], time_elec_alone[best_b],
                        timestamp_elec_alone[best_a], timestamp_elec_alone[best_b]);

    analyzeAlphasForPair(start_p1, start_p2,
                        timestamp_elec_alone[best_a], timestamp_elec_alone[best_b],
                        x_alpha, y_alpha, z_alpha);

    tree->Fill();
}

// ============================================================================
// Compute pair geometry
// ============================================================================
void falaise_skeleton_module_ptd::computePairGeometry(
    int best_a, int best_b,
    const std::vector<std::vector<double>>& kinks_vec) {

    // Bounds checking
    if (best_a >= (int)num_om_elec.size() || best_b >= (int)num_om_elec.size()) {
        return;
    }

    // Check if same calo, same cluster, same side
    if (num_om_elec[best_a] == num_om_elec[best_b]) {
        hit_the_same_calo_hit = true;
    }
    if (cluster_elec_num[best_a] == cluster_elec_num[best_b]) {
        same_cluster_elec = true;
    }
    if (side_elec[best_a] == side_elec[best_b]) {
        same_side_elec = true;
    }

    // Compute spatial separation
    delta_y_elec = abs(vertex_3D_start_y[best_a] - vertex_3D_start_y[best_b]);
    delta_z_elec = abs(vertex_3D_start_z[best_a] - vertex_3D_start_z[best_b]);
    delta_r_calo = sqrt(pow(vertex_3D_end_y[best_a] - vertex_3D_end_y[best_b], 2) +
                       pow(vertex_3D_end_z[best_a] - vertex_3D_end_z[best_b], 2));

    // Projection onto closest point between two lines
    double segA_start[3] = {vertex_ttd_start_x[best_a], vertex_ttd_start_y[best_a],
                           vertex_ttd_start_z[best_a]};
    double segB_start[3] = {vertex_ttd_start_x[best_b], vertex_ttd_start_y[best_b],
                           vertex_ttd_start_z[best_b]};
    double segA_end[3] = {vertex_ttd_end_x[best_a], vertex_ttd_end_y[best_a],
                         vertex_ttd_end_z[best_a]};
    double segB_end[3] = {vertex_ttd_end_x[best_b], vertex_ttd_end_y[best_b],
                         vertex_ttd_end_z[best_b]};

    // Bounds checking for kinks_vec
    if (best_a < (int)kinks_vec.size() && !kinks_vec[best_a].empty()) {
        segA_end[0] = kinks_vec[best_a][0];
        segA_end[1] = kinks_vec[best_a][1];
        segA_end[2] = kinks_vec[best_a][2];
    }
    if (best_b < (int)kinks_vec.size() && !kinks_vec[best_b].empty()) {
        segB_end[0] = kinks_vec[best_b][0];
        segB_end[1] = kinks_vec[best_b][1];
        segB_end[2] = kinks_vec[best_b][2];
    }

    Vec3 A(segA_start[0], segA_start[1], segA_start[2]);
    Vec3 B(segA_end[0], segA_end[1], segA_end[2]);
    Vec3 A2(segB_start[0], segB_start[1], segB_start[2]);
    Vec3 B2(segB_end[0], segB_end[1], segB_end[2]);
    Vec3 C, C2;

    if (closestPointsBetweenLines(A, B, A2, B2, C, C2)) {
        Vec3 M = (C + C2) * 0.5;
        vertex_3D_projected_x = M.x;
        vertex_3D_projected_y = M.y;
        vertex_3D_projected_z = M.z;
    }
}

// ============================================================================
// Compute pair timing
// ============================================================================
void falaise_skeleton_module_ptd::computePairTiming(
    int best_a, int best_b,
    const double* start_p1, const double* end_p1,
    const double* start_p2, const double* end_p2) {

    // Bounds checking
    if (best_a >= (int)energy_elec.size() || best_b >= (int)energy_elec.size()) {
        return;
    }

    double new_c = 299.792458;  // mm/ns
    double track_length_1 = track_length_vector[best_a];
    double track_length_2 = track_length_vector[best_b];


    double beta_1 = sqrt(energy_elec[best_a] * (energy_elec[best_a] + 2 * 0.511)) /
                   (energy_elec[best_a] + 0.511);
    double beta_2 = sqrt(energy_elec[best_b] * (energy_elec[best_b] + 2 * 0.511)) /
                   (energy_elec[best_b] + 0.511);

    t1_th = track_length_1 / (beta_1 * new_c);
    t2_th = track_length_2 / (beta_2 * new_c);

    internal_theoretical_time_diff = (t1_th - t2_th) - diff_time_elec;
    external_theoretical_time_diff = std::abs(diff_time_elec) - (t1_th + t2_th);

    track_lenght.push_back(track_length_1);
    track_lenght.push_back(track_length_2);
}

// ============================================================================
// Analyze gammas for pair
// ============================================================================
void falaise_skeleton_module_ptd::analyzeGammasForPair(
    int best_a, int best_b,
    double tA, double tB,
    double timestampA, double timestampB) {

    std::vector<double> gamma_dt;
    std::vector<double> gamma_energy;

    for (size_t g = 0; g < time_gamma.size(); ++g) {
        double tgamma = time_gamma[g];
        double dtA = std::abs(tgamma - tA);
        double dtB = std::abs(tgamma - tB);
        double dt_min = std::min(dtA, dtB);
        gamma_dt.push_back(dt_min);
        gamma_energy.push_back(energy_gamma[g]);

        if (nb_gamma == 1 && vertex_gamma[0] != side_elec[best_a] && same_side_elec == 1) {
            opposite_side_e_gamma = true;
        }
    }

    closest_gamma = 1e9;
    for (double dt : gamma_dt) {
        if (std::abs(dt) < std::abs(closest_gamma))
            closest_gamma = dt;
    }
}

// ============================================================================
// Analyze alphas for pair
// ============================================================================
void falaise_skeleton_module_ptd::analyzeAlphasForPair(
    const double* start_p1, const double* start_p2,
    double timestampA, double timestampB,
    const std::vector<double>& x_alpha,
    const std::vector<double>& y_alpha,
    const std::vector<double>& z_alpha) {

    // Time-based alpha veto
    double closest_alpha_dt_delayed = 1e9;
    for (double t_alpha : mean_alpha_anodic_timestamp) {
        double dt = std::min(std::abs(t_alpha - timestampA), std::abs(t_alpha - timestampB));
        if (dt < closest_alpha_dt_delayed)
            closest_alpha_dt_delayed = dt;
    }
    closest_time_track = closest_alpha_dt_delayed * 1e-3;

    // Distance-based alpha veto
    double closest_alpha_dist = 1e9;
    for (size_t i = 0; i < x_alpha.size(); ++i) {
        double dxA = x_alpha[i] - start_p1[0];
        double dyA = y_alpha[i] - start_p1[1];
        double dzA = z_alpha[i] - start_p1[2];
        double dist1 = std::sqrt(dxA * dxA + dyA * dyA + dzA * dzA);

        double dxB = x_alpha[i] - start_p2[0];
        double dyB = y_alpha[i] - start_p2[1];
        double dzB = z_alpha[i] - start_p2[2];
        double dist2 = std::sqrt(dxB * dxB + dyB * dyB + dzB * dzB);
        
        double dist_min = std::min(dist1, dist2);
        if (dist_min < closest_alpha_dist) {
            closest_alpha_dist = dist_min;
        }
    }
    closest_track = closest_alpha_dist;
}

// ============================================================================
// Main process function
// ============================================================================
dpp::chain_module::process_status falaise_skeleton_module_ptd::process(datatools::things& event) {
    resetEventVariables();

    // Skip if no PTD bank
    if (!event.has("PTD")) {
        std::cout << "No PTD bank in event " << ptd_event_counter++ << std::endl;
        event_number++;
        return dpp::base_module::PROCESS_SUCCESS;
    }

    const snemo::datamodel::particle_track_data& PTD =
        event.get<snemo::datamodel::particle_track_data>("PTD");

    // ===== STEP 1: Process simulated data =====
    processSimulatedData(event);

    if (!ptd_details) {
        event_number++;
        return dpp::base_module::PROCESS_SUCCESS;
    }

    // ===== STEP 2: Process event header =====
    if (!processEventHeader(event)) {
        return dpp::base_module::PROCESS_SUCCESS;
    }

    // ===== STEP 3: Process tracker data =====
    processTrackerData(event);

    // ===== STEP 4: Process gammas =====
    processGammas(PTD);

    // ===== STEP 5: Process electrons =====
    std::vector<std::vector<double>> kinks_vec;
    std::vector<double> x_alpha, y_alpha, z_alpha;
    std::vector<int> indices_for_pairs;
    std::vector<int> cluster_id;

    const snemo::datamodel::precalibrated_data& pCD =
        event.get<snemo::datamodel::precalibrated_data>("pCD");

    processElectrons(PTD, pCD, kinks_vec, x_alpha, y_alpha, z_alpha,
                    indices_for_pairs, cluster_id);

    // ===== STEP 6: Process pairs =====
    processPairs(indices_for_pairs, cluster_id, kinks_vec, x_alpha, y_alpha, z_alpha);

    event_number++;
    return dpp::base_module::PROCESS_SUCCESS;
}

double falaise_skeleton_module_ptd::extract_unix_start_time(int run_number) const {
    std::string cmd = "grep ^" + std::to_string(run_number) +
                     " /sps/nemo/scratch/chauveau/commissioning/software/run-sync-time.txt | awk '{print $2}'";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) {
        std::cerr << "Error: cannot execute grep\n";
        return -1;
    }
    double timestamp = -1;
    fscanf(pipe, "%lf", &timestamp);
    pclose(pipe);
    return timestamp;
}

// ============================================================================
// Phase determination based on run number (from data-taking table)
// Phase 0: runs 1546-1798  (high threshold 300 keV + 7% tracker missing)
// Phase 1: runs 2011-2183  (high threshold 30 keV  + 7% tracker missing)
// Phase 2: runs 2683-3467  (high threshold 30 keV  + full tracker)
// Phase 3: runs 3470+      (phase 2 + low Rn activity)
// ============================================================================
int falaise_skeleton_module_ptd::get_phase(int run_number) {
    if(run_number == -2) return _simulation_phase_;
    if (run_number >= 1546 && run_number <= 1798) return 0;
    if (run_number >= 2011 && run_number <= 2183) return 1;
    if (run_number >= 2683 && run_number <= 2869) return 2;
    if (run_number >= 2869 && run_number <= 3467) return 3;
    if (run_number >= 3470) return 4;
    return -1;
}
