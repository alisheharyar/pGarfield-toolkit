// Boost library to parse the command line parameters and INI file
#include <boost/program_options.hpp>
namespace po = boost::program_options;

std::string card_filename;

// parameters from the card file
// Dimensions
double param_dim_pitch;
double param_dim_kapton;
double param_dim_metal;
double param_dim_outdia;
double param_dim_middia;

// Gas
std::string param_gas1_name;
double param_gas1_comp;
std::string param_gas2_name;
double param_gas2_comp;
std::string param_gas3_name;
double param_gas3_comp;
double param_gas_temperature;
double param_gas_pressure;
double param_gas_penning;
std::string param_ion_filepath;

// Sensor parameters
double param_sensor_x1, param_sensor_y1, param_sensor_z1, param_sensor_x2, param_sensor_y2, param_sensor_z2;

// Track
double param_track_momentum;
std::string param_track_particle;

// Misc
std::string param_mesh_dir;
double param_num_events;

// Optimizations
bool param_opt_cache_boundingboxes;
bool param_opt_search_tetratree;
bool param_opt_search_neighbors;

bool read_card(int argc, char * argv[]) {
  po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
            ("card", po::value<std::string>(), "set card file")
        ;

  po::variables_map vm;        
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    

  if (vm.count("help")) {
      std::cout << desc << "\n";
      return false;
  }

  if (vm.count("card")) {
      std::cout << "Card file was set to " << vm["card"].as<std::string>() << ".\n";
      card_filename = vm["card"].as<std::string>();
            //<< vm["card"].as<double>() << ".\n";
  } else {
      std::cout << "Card file was not set.\n";
      return false;
  }


  // Setup options.
  po::options_description params("Parameters");
  params.add_options()
    ("Dimensions.Pitch", po::value< double >( &param_dim_pitch ), "dim pitch" )
    ("Dimensions.Kapton", po::value< double >( &param_dim_kapton ), "dim kapton" )
    ("Dimensions.Metal", po::value< double >( &param_dim_metal ), "dim metal" )
    ("Dimensions.Outdia", po::value< double >( &param_dim_outdia ), "dim outdia" )
    ("Dimensions.Middia", po::value< double >( &param_dim_middia ), "dim middia" )
    ("Gas.Name1", po::value< std::string >( &param_gas1_name ), "gas1 name" )
    ("Gas.Comp1", po::value< double >( &param_gas1_comp ), "gas1 comp" )
    ("Gas.Name2", po::value< std::string >( &param_gas2_name ), "gas2 name" )
    ("Gas.Comp2", po::value< double >( &param_gas2_comp ), "gas2 comp" )
    ("Gas.Name3", po::value< std::string >( &param_gas3_name ), "gas3 name" )
    ("Gas.Comp3", po::value< double >( &param_gas3_comp ), "gas3 comp" )
    ("Gas.Temperature", po::value< double >( &param_gas_temperature ), "gas temperature" )
    ("Gas.Pressure", po::value< double >( &param_gas_pressure ), "gas pressure" )
    ("Gas.Penning", po::value< double >( &param_gas_penning ), "gas penning" )
    ("Gas.IonMobilityFile", po::value< std::string >( &param_ion_filepath ), "gas ion mobility file" )
    ("Sensor.X1", po::value< double >( &param_sensor_x1 ), "sensor x1" )
    ("Sensor.Y1", po::value< double >( &param_sensor_y1 ), "sensor y1" )
    ("Sensor.Z1", po::value< double >( &param_sensor_z1 ), "sensor z1" )
    ("Sensor.X2", po::value< double >( &param_sensor_x2 ), "sensor x2" )
    ("Sensor.Y2", po::value< double >( &param_sensor_y2 ), "sensor y2" )
    ("Sensor.Z2", po::value< double >( &param_sensor_z2 ), "sensor z2" )
    ("Track.Momentum", po::value< double >( &param_track_momentum ), "track momentum" )
    ("Track.Particle", po::value< std::string >( &param_track_particle ), "track particle" )
    
    ("Misc.MeshDirectory", po::value< std::string >( &param_mesh_dir ), "mesh directory" )
    ("Misc.NumEvents", po::value< double >( &param_num_events ), "number of events" )
	("Optimizations.CacheBoundingBoxes", po::value< bool >( &param_opt_cache_boundingboxes ), "cache bounding boxes" )
	("Optimizations.SearchWithTree", po::value< bool >( &param_opt_search_tetratree ), "search with tetrahedral tree" )
	("Optimizations.SearchNeighbors", po::value< bool >( &param_opt_search_neighbors ), "search through neighbors" )
    ;

  // Load setting file.
  po::variables_map vm2;
  std::ifstream settings_file( card_filename.c_str() , std::ifstream::in );
  po::store( po::parse_config_file( settings_file , params ), vm2 );
  settings_file.close();
  po::notify( vm2 );

  std::cout << std::endl;
  std::cout << "Parameters from the card file: \n";
  std::cout << "Gas.gas1_name: " << param_gas1_name << std::endl;
  std::cout << "Gas.gas1_comp: " << param_gas1_comp << std::endl;
  std::cout << "Gas.gas2_name: " << param_gas2_name << std::endl;
  std::cout << "Gas.gas2_comp: " << param_gas2_comp << std::endl;
  std::cout << "Gas.gas3_name: " << param_gas3_name << std::endl;
  std::cout << "Gas.gas3_comp: " << param_gas3_comp << std::endl;
  std::cout << "Gas.temperature: " << param_gas_temperature << std::endl;
  std::cout << "Gas.Pressure: " << param_gas_pressure << std::endl;
  std::cout << "Gas.penning: " << param_gas_penning << std::endl;
  std::cout << "Gas.IonMobilityFile: " << param_ion_filepath << std::endl;
  std::cout << "Sensor.x1: " << param_sensor_x1 << std::endl;
  std::cout << "Sensor.y1: " << param_sensor_y1 << std::endl;
  std::cout << "Sensor.z1: " << param_sensor_z1 << std::endl;
  std::cout << "Sensor.x2: " << param_sensor_x2 << std::endl;
  std::cout << "Sensor.y2: " << param_sensor_y2 << std::endl;
  std::cout << "Sensor.z2: " << param_sensor_z2 << std::endl;
  std::cout << "Misc.MeshDirectory: " << param_mesh_dir << std::endl;
  std::cout << "Misc.NumEvents: " << param_num_events << std::endl;
  std::cout << "Dimsensions.Pitch: " << param_dim_pitch << std::endl;
  std::cout << "Dimsensions.Kapton: " << param_dim_kapton << std::endl;
  std::cout << "Dimsensions.Metal: " << param_dim_metal << std::endl;
  std::cout << "Dimsensions.Outdia: " << param_dim_outdia << std::endl;
  std::cout << "Dimsensions.Middia: " << param_dim_middia << std::endl;
  std::cout << "Optimizations.CacheBoundingBoxes: " << param_opt_cache_boundingboxes << std::endl;
  std::cout << "Optimizations.SearchWithTree: " << param_opt_search_tetratree << std::endl;
  std::cout << "Optimizations.SearchNeighbors: " << param_opt_search_neighbors << std::endl;
  std::cout << "Track.Momentum: " << param_track_momentum << std::endl;
  std::cout << "Track.Particle: " << param_track_particle << std::endl;
  std::cout << std::endl;

  return true;
}
