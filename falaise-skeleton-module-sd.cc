#include <bayeux/dpp/chain_module.h>
#include <bayeux/mctools/base_step_hit.h>
#include <bayeux/mctools/simulated_data.h>

#include <falaise/snemo/datamodels/gg_track_utils.h>

class falaise_skeleton_module_sd : public dpp::chain_module
{

public:

  falaise_skeleton_module_sd();

  virtual ~falaise_skeleton_module_sd();
    
  virtual void initialize (const datatools::properties &,
                           datatools::service_manager &,
			   dpp::module_handle_dict_type &);

  dpp::chain_module::process_status process (datatools::things & event);

  virtual void finalize ();
  
private:
  int sd_event_counter;

  // Macro to register the module
  DPP_MODULE_REGISTRATION_INTERFACE(falaise_skeleton_module_sd);

};

////////////////////////////////////////////////////////////////////

// Macro to register the module in the global register of data processing modules:
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


void falaise_skeleton_module_sd::initialize (const datatools::properties &, datatools::service_manager &, dpp::module_handle_dict_type &)
{
  sd_event_counter = 0;

  this->_set_initialized(true);
}

dpp::chain_module::process_status falaise_skeleton_module_sd::process (datatools::things & event)
{
  std::cout << "======== SD bank of event " << sd_event_counter++ << " ========" << std::endl;

  // Retrieve the SD bank  
  const mctools::simulated_data & SD = event.get<mctools::simulated_data>("SD");

  // Calorimeter step hit categories (3 differents labels)
  const std::vector<std::string> calo_hit_categories = {"calo", "xcalo", "gveto"};
  
  for (std::string calo_hit_category : calo_hit_categories)
    {
      if (SD.has_step_hits(calo_hit_category))
	{
    	  std::cout << "=> " << SD.get_step_hits(calo_hit_category).size() << " '" << calo_hit_category << "' step hit(s):" << std::endl;

	  // browse the calorimeter step hits (mctools::base_step_hit)
	  for (const auto & a_step_hit : SD.get_step_hits(calo_hit_category))
	    {
	      const geomtools::geom_id & calo_geom_id = a_step_hit->get_geom_id();
	      const double time_ns = a_step_hit->get_time_start() / CLHEP::ns;
	      const double edep_mev = a_step_hit->get_energy_deposit() / CLHEP::MeV; 
	      std::cout << " - " << calo_geom_id << "   time = " << time_ns << " ns   edep = " << edep_mev << " MeV" << std::endl;
	    }
	}
      else
	{
	  std::cout << "=> no '" << calo_hit_category << "' step hit"  << std::endl;
	  continue;
	}
    }
  
  // Tracker step hit category (1 unique label)
  if (SD.has_step_hits("gg"))
    {
      std::cout << "=> " << SD.get_step_hits("gg").size() << " 'gg' step hit(s):" << std::endl;

      for (const auto & a_step_hit : SD.get_step_hits("gg"))
	{
	  const geomtools::geom_id & gg_geom_id = a_step_hit->get_geom_id();

	  const datatools::properties & step_hit_aux = a_step_hit->get_auxiliaries();
	  const double true_radius_cm = step_hit_aux.fetch_real(snemo::datamodel::gg_track::minimum_approach_distance_key())/CLHEP::cm;

	  std::vector<double> minimum_approach_position;
	  step_hit_aux.fetch(snemo::datamodel::gg_track::minimum_approach_position_key(), minimum_approach_position);
	  const double & true_height_m = minimum_approach_position.at(2)/CLHEP::m;

	  std::cout << " - " << gg_geom_id << "   radius = " << true_radius_cm << " cm   height = " << true_height_m << " m" << std::endl;
	}
    }
  else
    {
      std::cout << "=> no 'gg' step hit"  << std::endl;
    }

  return dpp::base_module::PROCESS_SUCCESS;
}


void falaise_skeleton_module_sd::finalize()
{
  std::cout << "falaise_skeleton_module_sd::finalize()" << std::endl;
}
