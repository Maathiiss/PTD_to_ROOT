#include <bayeux/dpp/chain_module.h>
#include <falaise/snemo/datamodels/tracker_trajectory_data.h>
#include <falaise/snemo/datamodels/line_trajectory_pattern.h>
#include <falaise/snemo/datamodels/helix_trajectory_pattern.h>

class falaise_skeleton_module_ttd : public dpp::chain_module
{

public:
  // Constructor
  falaise_skeleton_module_ttd();

  // Destructor
  virtual ~falaise_skeleton_module_ttd();

  // Initialisation function
  virtual void initialize (const datatools::properties &,
                           datatools::service_manager &,
			   dpp::module_handle_dict_type &);

  // Event processing function
  dpp::chain_module::process_status process (datatools::things & event);
  
private:
  int  ttd_event_counter;
  bool ttd_details;

  // Macro to register the module
  DPP_MODULE_REGISTRATION_INTERFACE(falaise_skeleton_module_ttd);

};

////////////////////////////////////////////////////////////////////

// Macro to add the module in the global register of data processing modules:
// The module defined by this class 'falaise_skeleton_module_ttd' will be registered
// with the label ID 'FalaiseSkeletonModule_TTD' (to use in pipeline configuration file)
DPP_MODULE_REGISTRATION_IMPLEMENT(falaise_skeleton_module_ttd, "FalaiseSkeletonModule_TTD")


falaise_skeleton_module_ttd::falaise_skeleton_module_ttd()
{
  std::cout << "falaise_skeleton_module_ttd::falaise_skeleton_module_ttd() called" << std::endl;
}


falaise_skeleton_module_ttd::~falaise_skeleton_module_ttd()
{
  std::cout << "falaise_skeleton_module_ttd::~falaise_skeleton_module_ttd() called" << std::endl;
}


void falaise_skeleton_module_ttd::initialize (const datatools::properties & module_properties, datatools::service_manager &, dpp::module_handle_dict_type &)
{
  std::cout << "falaise_skeleton_module_ttd::initialize() called" << std::endl;

  ttd_event_counter = 0;

  if ( module_properties.has_key("ttd_details"))
    ttd_details = module_properties.fetch_boolean("ttd_details");
  else ttd_details = false;

  this->_set_initialized(true);
}

dpp::chain_module::process_status falaise_skeleton_module_ttd::process (datatools::things & event)
{
  // Skip processing if TTD bank is not present
  if (!event.has("TTD"))
    {
      std::cout << "======== no TTD bank in event " << ttd_event_counter++ << " ========" << std::endl;
      return dpp::base_module::PROCESS_SUCCESS;
    }

  // Retrieve the TTD bank  
  const snemo::datamodel::tracker_trajectory_data & TTD = event.get<snemo::datamodel::tracker_trajectory_data>("TTD");
  std::cout << "======== TTD bank of event " << ttd_event_counter++ << " ========" << std::endl;

  // Browse calibrated calorimeter hits
  std::cout << "=> " << TTD.get_solutions().size() << " trajectory solution(s)" << std::endl;

  if (ttd_details)
    {
      const snemo::datamodel::tracker_trajectory_solution & ttd_solution = TTD.get_default_solution();

      for (const datatools::handle<snemo::datamodel::tracker_trajectory> & ttd_trajectory : ttd_solution.get_trajectories())
	{
	  if (! ttd_trajectory->get_fit_infos().is_best())
	    continue;

	  const snemo::datamodel::base_trajectory_pattern & tracker_pattern = ttd_trajectory->get_pattern();
	  const std::string & pattern_id = tracker_pattern.get_pattern_id();

	  if (pattern_id == snemo::datamodel::line_trajectory_pattern::pattern_id())
	    {
	      const auto & line_trajectory = dynamic_cast<const snemo::datamodel::line_trajectory_pattern &> (tracker_pattern);
	      std::cout << " - line trajectory" << std::endl;
	    }

	  else if (pattern_id == snemo::datamodel::helix_trajectory_pattern::pattern_id())
	    {
	      const auto & helix_trajectory = dynamic_cast<const snemo::datamodel::helix_trajectory_pattern &> (tracker_pattern);
	      std::cout << " - helix_trajectory" << std::endl;
	    }
	}
    }


  return dpp::base_module::PROCESS_SUCCESS;
}

