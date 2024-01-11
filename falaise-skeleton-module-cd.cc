#include <bayeux/dpp/chain_module.h>
#include <falaise/snemo/datamodels/calibrated_data.h>

#include "CLHEP/Units/SystemOfUnits.h"

class falaise_skeleton_module_cd : public dpp::chain_module
{

public:
  // Constructor
  falaise_skeleton_module_cd();

  // Destructor
  virtual ~falaise_skeleton_module_cd();

  // Initialisation function
  virtual void initialize (const datatools::properties &,
                           datatools::service_manager &,
			   dpp::module_handle_dict_type &);

  // Event processing function
  dpp::chain_module::process_status process (datatools::things & event);
  
private:
  int  cd_event_counter;
  bool cd_calo_details;
  bool cd_tracker_details;

  // Macro to register the module
  DPP_MODULE_REGISTRATION_INTERFACE(falaise_skeleton_module_cd);

};

////////////////////////////////////////////////////////////////////

// Macro to add the module in the global register of data processing modules:
// The module defined by this class 'falaise_skeleton_module_cd' will be registered
// with the label ID 'FalaiseSkeletonModule_CD' (to use in pipeline configuration file)
DPP_MODULE_REGISTRATION_IMPLEMENT(falaise_skeleton_module_cd, "FalaiseSkeletonModule_CD")


falaise_skeleton_module_cd::falaise_skeleton_module_cd()
{
  std::cout << "falaise_skeleton_module_cd::falaise_skeleton_module_cd() called" << std::endl;
}


falaise_skeleton_module_cd::~falaise_skeleton_module_cd()
{
  std::cout << "falaise_skeleton_module_cd::~falaise_skeleton_module_cd() called" << std::endl;
}


void falaise_skeleton_module_cd::initialize (const datatools::properties & module_properties, datatools::service_manager &, dpp::module_handle_dict_type &)
{
  std::cout << "falaise_skeleton_module_cd::initialize() called" << std::endl;

  cd_event_counter = 0;

  if ( module_properties.has_key("calo_details"))
    cd_calo_details = module_properties.fetch_boolean("calo_details");
  else cd_calo_details = false;

  if ( module_properties.has_key("tracker_details"))
    cd_tracker_details = module_properties.fetch_boolean("tracker_details");
  else cd_tracker_details = false;

  this->_set_initialized(true);
}

dpp::chain_module::process_status falaise_skeleton_module_cd::process (datatools::things & event)
{
  // Skip processing if CD bank is not present
  if (!event.has("CD"))
    {
      std::cout << "======== no CD bank in event " << cd_event_counter++ << " ========" << std::endl;
      return dpp::base_module::PROCESS_SUCCESS;
    }

  // Retrieve the CD bank  
  const snemo::datamodel::calibrated_data & CD = event.get<snemo::datamodel::calibrated_data>("CD");
  std::cout << "======== CD bank of event " << cd_event_counter++ << " ========" << std::endl;

  // Browse calibrated calorimeter hits
  std::cout << "=> " << CD.calorimeter_hits().size() << " calo hit(s)" << std::endl;
  if (cd_calo_details)
    {
      for (const auto & calo_hit : CD.calorimeter_hits())
	{
	  const geomtools::geom_id & calo_geom_id = calo_hit->get_geom_id();
	  const double calo_time_ns  = calo_hit->get_time() / CLHEP::ns;
	  const double calo_energy_mev = calo_hit->get_energy() / CLHEP::MeV; 
	  
	  std::cout << " - " << calo_geom_id << "   time = " << calo_time_ns << " ns   edep = " << calo_energy_mev << " MeV" << std::endl;
	}
    }

  // Browse calibrated tracker hits
  std::cout << "=> " << CD.tracker_hits().size() << " tracker hit(s)" << std::endl;
  if (cd_tracker_details)
    {
      for (const auto & tracker_hit : CD.tracker_hits())
	{
	  const geomtools::geom_id & tracker_geom_id = tracker_hit->get_geom_id();
	  const double tracker_radius_cm = tracker_hit->get_r() / CLHEP::cm;
	  const double tracker_height_m  = tracker_hit->get_z() / CLHEP::m;
	  std::cout << " - " << tracker_geom_id << "   r = " << tracker_radius_cm << " cm   z = " << tracker_height_m << " m" << std::endl;
	}
    }

  return dpp::base_module::PROCESS_SUCCESS;
}

