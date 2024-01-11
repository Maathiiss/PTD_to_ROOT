#include <bayeux/dpp/chain_module.h>
#include <falaise/snemo/datamodels/particle_track_data.h>

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
  bool ptd_details;

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
  std::cout << "falaise_skeleton_module_ptd::falaise_skeleton_module_ptd() called" << std::endl;
}


falaise_skeleton_module_ptd::~falaise_skeleton_module_ptd()
{
  std::cout << "falaise_skeleton_module_ptd::~falaise_skeleton_module_ptd() called" << std::endl;
}


void falaise_skeleton_module_ptd::initialize (const datatools::properties & module_properties, datatools::service_manager &, dpp::module_handle_dict_type &)
{
  std::cout << "falaise_skeleton_module_ptd::initialize() called" << std::endl;

  ptd_event_counter = 0;

  if ( module_properties.has_key("ptd_details"))
    ptd_details = module_properties.fetch_boolean("ptd_details");
  else ptd_details = false;

  this->_set_initialized(true);
}

dpp::chain_module::process_status falaise_skeleton_module_ptd::process (datatools::things & event)
{
  // Skip processing if PTD bank is not present
  if (!event.has("PTD"))
    {
      std::cout << "======== no PTD bank in event " << ptd_event_counter++ << " ========" << std::endl;
      return dpp::base_module::PROCESS_SUCCESS;
    }

  // Retrieve the PTD bank  
  const snemo::datamodel::particle_track_data & PTD = event.get<snemo::datamodel::particle_track_data>("PTD");
  std::cout << "======== PTD bank of event " << ptd_event_counter++ << " ========" << std::endl;

  // Browse calibrated calorimeter hits
  std::cout << "=> " << PTD.particles().size() << " particle(s)" << std::endl;

  if (ptd_details)
    {
      for (const datatools::handle<snemo::datamodel::particle_track> & particle : PTD.particles())
	{
	  std::cout << " - particle with charge = ";

	  switch (particle->get_charge())
	    {
	    case snemo::datamodel::particle_track::CHARGE_UNDEFINED:
	      std::cout << "? " ; break;
	    case snemo::datamodel::particle_track::CHARGE_NEUTRAL:
	      std::cout << "0 " ; break;
	    case snemo::datamodel::particle_track::CHARGE_POSITIVE:
	      std::cout << "+1 " ; break;
	    case snemo::datamodel::particle_track::CHARGE_NEGATIVE:
	      std::cout << "-1 " ; break;
	    }

	  const std::vector<datatools::handle<snemo::datamodel::vertex>> & particle_vertices = particle->get_vertices();
	  std::cout << "and " << particle_vertices.size() << " vertice(s)" << std::endl;
	}
    }

  return dpp::base_module::PROCESS_SUCCESS;
}

