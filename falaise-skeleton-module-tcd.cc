#include <bayeux/dpp/chain_module.h>
#include <falaise/snemo/datamodels/tracker_clustering_data.h>

class falaise_skeleton_module_tcd : public dpp::chain_module
{

public:
  // Constructor
  falaise_skeleton_module_tcd();

  // Destructor
  virtual ~falaise_skeleton_module_tcd();

  // Initialisation function
  virtual void initialize (const datatools::properties &,
                           datatools::service_manager &,
			   dpp::module_handle_dict_type &);

  // Event processing function
  dpp::chain_module::process_status process (datatools::things & event);
  
private:
  int  tcd_event_counter;
  bool tcd_details;

  // Macro to register the module
  DPP_MODULE_REGISTRATION_INTERFACE(falaise_skeleton_module_tcd);

};

////////////////////////////////////////////////////////////////////

// Macro to add the module in the global register of data processing modules:
// The module defined by this class 'falaise_skeleton_module_tcd' will be registered
// with the label ID 'FalaiseSkeletonModule_TCD' (to use in pipeline configuration file)
DPP_MODULE_REGISTRATION_IMPLEMENT(falaise_skeleton_module_tcd, "FalaiseSkeletonModule_TCD")


falaise_skeleton_module_tcd::falaise_skeleton_module_tcd()
{
  std::cout << "falaise_skeleton_module_tcd::falaise_skeleton_module_tcd() called" << std::endl;
}


falaise_skeleton_module_tcd::~falaise_skeleton_module_tcd()
{
  std::cout << "falaise_skeleton_module_tcd::~falaise_skeleton_module_tcd() called" << std::endl;
}


void falaise_skeleton_module_tcd::initialize (const datatools::properties & module_properties, datatools::service_manager &, dpp::module_handle_dict_type &)
{
  std::cout << "falaise_skeleton_module_tcd::initialize() called" << std::endl;

  tcd_event_counter = 0;

  if ( module_properties.has_key("tcd_details"))
    tcd_details = module_properties.fetch_boolean("tcd_details");
  else tcd_details = false;

  this->_set_initialized(true);
}

dpp::chain_module::process_status falaise_skeleton_module_tcd::process (datatools::things & event)
{
  // Skip processing if TCD bank is not present
  if (!event.has("TCD"))
    {
      std::cout << "======== no TCD bank in event " << tcd_event_counter++ << " ========" << std::endl;
      return dpp::base_module::PROCESS_SUCCESS;
    }

  // Retrieve the TCD bank  
  const snemo::datamodel::tracker_clustering_data & TCD = event.get<snemo::datamodel::tracker_clustering_data>("TCD");
  std::cout << "======== TCD bank of event " << tcd_event_counter++ << " ========" << std::endl;

  // Browse calibrated calorimeter hits
  std::cout << "=> " << TCD.solutions().size() << " clustering solution(s)" << std::endl;

  if (tcd_details)
    {
      for (const auto & tcd_solution : TCD.solutions())
	{
	  const int & tcd_id = tcd_solution->get_solution_id();
	  const auto & tcd_clusters = tcd_solution->get_clusters();

	  // const auto & tcd_unclustered_hits = tcd_solution->get_unclustered_hits();
	  std::cout << " - solution #" << tcd_id << " with " << tcd_clusters.size() << " clusters " << std::endl;

	  // for (const auto & tcd_cluster : tcd_clusters)
	  //   {
	  //     const int   & cluster_id = tcd_cluster->get_cluster_id();

	  //     const bool  & is_prompt  = tcd_cluster->is_prompt();
	  //     const bool  & is_delayed = tcd_cluster->is_delayed();

	  //     for (const auto & tracker_hit : tcd_cluster->hits())
	  // 	{
	  // 	  const geomtools::geom_id & tracker_geom_id = tracker_hit->get_geom_id();
	  // 	  const double tracker_radius_cm = tracker_hit->get_r() / CLHEP::cm;
	  // 	  const double tracker_height_m  = tracker_hit->get_z() / CLHEP::m;
	  // 	}
	  //   }
	}
    }

  return dpp::base_module::PROCESS_SUCCESS;
}

