#include <bayeux/dpp/chain_module.h>
#include <bayeux/mctools/base_step_hit.h>
#include <bayeux/mctools/simulated_data.h>

#include <falaise/snemo/datamodels/gg_track_utils.h>

class falaise_skeleton_module_sd : public dpp::chain_module
{

public:
  // Constructor
  falaise_skeleton_module_sd();

  // Destructor
  virtual ~falaise_skeleton_module_sd();

  // Initialisation function
  virtual void initialize (const datatools::properties &,
                           datatools::service_manager &,
			   dpp::module_handle_dict_type &);

  // Event processing function
  dpp::chain_module::process_status process (datatools::things & event);
  
private:
  int  sd_event_counter;
  bool sd_calo_details;
  bool sd_tracker_details;

  // Macro to register the module
  DPP_MODULE_REGISTRATION_INTERFACE(falaise_skeleton_module_sd);

};

////////////////////////////////////////////////////////////////////

// Macro to add the module in the global register of data processing modules:
// The module defined by this class 'falaise_skeleton_module_sd' will be registered
// with the label ID 'FalaiseSkeletonModule_SD' (to use in pipeline configuration file)
DPP_MODULE_REGISTRATION_IMPLEMENT(falaise_skeleton_module_sd, "FalaiseSkeletonModule_SD")


falaise_skeleton_module_sd::falaise_skeleton_module_sd()
{
  std::cout << "falaise_skeleton_module_sd::falaise_skeleton_module_sd() called" << std::endl;
}


falaise_skeleton_module_sd::~falaise_skeleton_module_sd()
{
  std::cout << "falaise_skeleton_module_sd::~falaise_skeleton_module_sd() called" << std::endl;
}


void falaise_skeleton_module_sd::initialize (const datatools::properties & module_properties, datatools::service_manager &, dpp::module_handle_dict_type &)
{
  std::cout << "falaise_skeleton_module_sd::initialize() called" << std::endl;

  sd_event_counter = 0;

  if ( module_properties.has_key("calo_details"))
    sd_calo_details = module_properties.fetch_boolean("calo_details");
  else sd_calo_details = false;

  if ( module_properties.has_key("tracker_details"))
    sd_tracker_details = module_properties.fetch_boolean("tracker_details");
  else sd_tracker_details = false;

  this->_set_initialized(true);
}

dpp::chain_module::process_status falaise_skeleton_module_sd::process (datatools::things & event)
{
  // Skip processing if SD bank is not present
  if (!event.has("SD"))
    {
      std::cout << "======== no SD bank in event " << sd_event_counter++ << " ========" << std::endl;
      return dpp::base_module::PROCESS_SUCCESS;
    }

  // Retrieve the SD bank  
  const mctools::simulated_data & SD = event.get<mctools::simulated_data>("SD");
  std::cout << "======== SD bank of event " << sd_event_counter++ << " ========" << std::endl;

  // Calorimeter step hit categories (3 differents labels)
  const std::vector<std::string> calo_hit_categories = {"calo", "xcalo", "gveto"};
  
  for (std::string calo_hit_category : calo_hit_categories)
    {
      size_t nb_calo_step_hits = SD.has_step_hits(calo_hit_category) ? SD.get_step_hits(calo_hit_category).size() : 0;
      std::cout << "=> " << nb_calo_step_hits << " '" << calo_hit_category << "' step hit(s)" << std::endl;

      if (sd_calo_details && (nb_calo_step_hits > 0))
	{
	  // Browse the calorimeter step hits
	  for (const auto & a_step_hit : SD.get_step_hits(calo_hit_category))
	    {
	      const geomtools::geom_id & calo_geom_id = a_step_hit->get_geom_id();
	      const double time_ns  = a_step_hit->get_time_start() / CLHEP::ns;
	      const double edep_mev = a_step_hit->get_energy_deposit() / CLHEP::MeV; 

	      std::cout << " - " << calo_geom_id << "   time = " << time_ns << " ns   edep = " << edep_mev << " MeV" << std::endl;
	    }
	}
    }
  
  // Tracker step hit category (1 unique label)
  size_t nb_tracker_step_hits = SD.has_step_hits("gg") ? SD.get_step_hits("gg").size() : 0;
  std::cout << "=> " << nb_tracker_step_hits << " 'gg' step hit(s)" << std::endl;

  if (sd_tracker_details && (nb_tracker_step_hits > 0))
    {
      // Browse the tracker step hits
      for (const auto & a_step_hit : SD.get_step_hits("gg"))
	{
	  const geomtools::geom_id & gg_geom_id = a_step_hit->get_geom_id();

	  // Minimal ionisation approach distance to anode wire and position is stored in auxiliaries
	  const datatools::properties & step_hit_aux = a_step_hit->get_auxiliaries();
	  const double true_radius_cm = step_hit_aux.fetch_real(snemo::datamodel::gg_track::minimum_approach_distance_key())/CLHEP::cm;
	  std::vector<double> minimum_approach_position;
	  step_hit_aux.fetch(snemo::datamodel::gg_track::minimum_approach_position_key(), minimum_approach_position);
	  const double true_height_m = minimum_approach_position.at(2)/CLHEP::m;

	  std::cout << " - " << gg_geom_id << "   radius = " << true_radius_cm << " cm   height = " << true_height_m << " m" << std::endl;
	}
    }

  return dpp::base_module::PROCESS_SUCCESS;
}

