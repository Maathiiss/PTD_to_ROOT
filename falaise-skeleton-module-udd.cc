#include <bayeux/dpp/chain_module.h>
#include <falaise/snemo/datamodels/unified_digitized_data.h>

#include "CLHEP/Units/SystemOfUnits.h"

const double CALO_ADC2MV = 2500./4096; // in mV/tick
const double CALO_SAMPLING_PERIOD_NS = 0.390625; // in ns
const double CALO_TDC2SEC = 6.25E-9; // in seconds
const double TRACKER_TDC2SEC = 12.5E-9; // in seconds

class falaise_skeleton_module_udd : public dpp::chain_module
{

public:
  // Constructor
  falaise_skeleton_module_udd();

  // Destructor
  virtual ~falaise_skeleton_module_udd();

  // Initialisation function
  virtual void initialize (const datatools::properties &,
                           datatools::service_manager &,
			   dpp::module_handle_dict_type &);

  // Event processing function
  dpp::chain_module::process_status process (datatools::things & event);
  
private:
  int udd_event_counter;
  bool udd_calo_details;
  bool udd_tracker_details;

  // Macro to register the module
  DPP_MODULE_REGISTRATION_INTERFACE(falaise_skeleton_module_udd);

};

////////////////////////////////////////////////////////////////////

// Macro to add the module in the global register of data processing modules:
// The module defined by this class 'falaise_skeleton_module_udd' will be registered
// with the label ID 'FalaiseSkeletonModule_UDD' (to use in pipeline configuration file)
DPP_MODULE_REGISTRATION_IMPLEMENT(falaise_skeleton_module_udd, "FalaiseSkeletonModule_UDD")


falaise_skeleton_module_udd::falaise_skeleton_module_udd()
{
  std::cout << "falaise_skeleton_module_udd::falaise_skeleton_module_udd() called" << std::endl;
}


falaise_skeleton_module_udd::~falaise_skeleton_module_udd()
{
  std::cout << "falaise_skeleton_module_udd::~falaise_skeleton_module_udd() called" << std::endl;
}


void falaise_skeleton_module_udd::initialize (const datatools::properties & module_properties, datatools::service_manager &, dpp::module_handle_dict_type &)
{
  std::cout << "falaise_skeleton_module_udd::initialize() called" << std::endl;

  udd_event_counter = 0;

  if ( module_properties.has_key("calo_details"))
    udd_calo_details = module_properties.fetch_boolean("calo_details");
  else udd_calo_details = false;

  if ( module_properties.has_key("tracker_details"))
    udd_tracker_details = module_properties.fetch_boolean("tracker_details");
  else udd_tracker_details = false;

  this->_set_initialized(true);
}

dpp::chain_module::process_status falaise_skeleton_module_udd::process (datatools::things & event)
{
  // Skip processing if UDD bank is not present
  if (!event.has("UDD"))
    {
      std::cout << "======== no UDD bank in event " << udd_event_counter++ << " ========" << std::endl;
      return dpp::base_module::PROCESS_SUCCESS;
    }

  // Retrieve the UDD bank  
  const snemo::datamodel::unified_digitized_data & UDD = event.get<snemo::datamodel::unified_digitized_data>("UDD");

  const int32_t & run_id = UDD.get_run_id();
  const int32_t & event_id = UDD.get_event_id();
  std::cout << "======== UDD bank of run " << run_id << " event " << event_id << " ========" << std::endl;

  // Browse UDD calorimeter hits
  const auto & udd_calo_hits = UDD.get_calorimeter_hits();
  std::cout << "=> " << udd_calo_hits.size() << " calo hit(s)" << std::endl;

  if (udd_calo_details)
    {
      for (const auto & udd_calo_hit : udd_calo_hits)
	{
	  const geomtools::geom_id & udd_calo_geom_id = udd_calo_hit->get_geom_id();

	  // const int64_t & udd_calo_tdc = udd_calo_hit.get_timestamp();

	  // udd_calo_hit.is_high_threshold_only();
	  // udd_calo_hit.is_low_threshold_only();


	  // // Waveform
	  // const std::vector<int16_t> & udd_calo_waveform = udd_calo_hit->get_waveform();
	  // const uint16_t udd_calo_fcr = udd_calo_hit->get_fcr(); // first read cell

	  // // LT counter statistics
	  // const uint16_t & udd_calo_lt_trigger_counter = udd_calo_hit->get_lt_trigger_counter();
	  // const uint32_t & udd_calo_lt_time_counter = udd_calo_hit->get_lt_time_counter();

	  // FW measurements
	  const float udd_calo_fwmeas_baseline_mv = (udd_calo_hit->get_fwmeas_baseline()/16.0) * CALO_ADC2MV;
	  const float udd_calo_fwmeas_amplitude_mv = (udd_calo_hit->get_fwmeas_peak_amplitude()/8.0) * CALO_ADC2MV;
	  const float udd_calo_fwmeas_charge_nvs = udd_calo_hit->get_fwmeas_charge() * 1e-3 * CALO_ADC2MV * CALO_SAMPLING_PERIOD_NS;
	  const float udd_calo_time_cfd_ns = (udd_calo_hit->get_fwmeas_rising_cell()/256.0) * CALO_SAMPLING_PERIOD_NS;

	  std::cout << " - " << udd_calo_geom_id
		    << " baseline = " << udd_calo_fwmeas_baseline_mv << " mV, "
		    << " amplitude = " << udd_calo_fwmeas_amplitude_mv << " mV, "
		    << " charge = " << udd_calo_fwmeas_charge_nvs << " nV.s, "
		    << " time = " << udd_calo_time_cfd_ns << " ns" << std::endl;
	}
    }


  // Browse UDD tracker hits
  const auto & udd_tracker_hits = UDD.get_tracker_hits();
  std::cout << "=> " << udd_tracker_hits.size() << " tracker hit(s)" << std::endl;

  if (udd_tracker_details)
    {
      for (const auto & udd_tracker_hit : udd_tracker_hits)
	{
	  const geomtools::geom_id & udd_tracker_geom_id = udd_tracker_hit->get_geom_id();

	  // only consider the front times here
	  const auto & udd_tracker_front_times = udd_tracker_hit->get_times().front();

	  double udd_tracker_anode_time = -1;
	  double udd_tracker_bottom_cathode_time = -1;
	  double udd_tracker_top_cathode_time = -1;

	  const uint16_t R0 = snemo::datamodel::tracker_digitized_hit::ANODE_R0;

	  if (udd_tracker_front_times.has_anode_time(R0))
	    udd_tracker_anode_time = udd_tracker_front_times.get_anode_time(R0) * TRACKER_TDC2SEC;

	  if (udd_tracker_front_times.has_bottom_cathode_time())
	    udd_tracker_bottom_cathode_time = udd_tracker_front_times.get_bottom_cathode_time();

	  if (udd_tracker_front_times.has_top_cathode_time())
	    udd_tracker_top_cathode_time = udd_tracker_front_times.get_top_cathode_time();

	  std::cout << " - " << udd_tracker_geom_id << " "
		    << (udd_tracker_front_times.has_anode_time(R0) ? "with R0 " : "without R0 ")
		    << (udd_tracker_front_times.has_bottom_cathode_time() ? "with R5 " : "without R5 ")
		    << (udd_tracker_front_times.has_top_cathode_time() ? "with R6" : "without R6")
		    << std::endl;
	}
    }

  return dpp::base_module::PROCESS_SUCCESS;
}

