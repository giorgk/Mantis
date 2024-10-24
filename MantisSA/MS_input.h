//
// Created by giorg on 10/16/2024.
//

#ifndef MANTISSA_MS_INPUT_H
#define MANTISSA_MS_INPUT_H

#include <boost/mpi.hpp>
#include <boost/program_options.hpp>
#include "SWAT_data.h"
#include "MS_structures.h"

namespace po = boost::program_options;

namespace MS {

    class UserInput{
    public:
        UserInput(boost::mpi::communicator& world_in);
        bool read(int argc, char* argv[]);

        int NsimYears;
        int porosity;
        int Nwells;
        double wc;
        double minDepth;
        double minRch;
        std::string configFile;
        std::string swat_input_file;
        std::string npsat_VI_file;
        std::string npsat_VD_file;
        std::string init_salt_VI_file;
        std::string init_salt_VD_file;
        std::string cell_well_file;
        std::string depth_input_file;
        std::string depth_name;
        std::string rch_input_file;

        RasterOptions rasteroptions;

        std::string Version;

    private:
        boost::mpi::communicator world;
    };

    UserInput::UserInput(boost::mpi::communicator &world_in)
        :
            world(world_in)
    {
        Version = "0.0.01";
    }

    bool UserInput::read(int argc, char **argv) {
        po::options_description commandLineOptions("Command line options");
        commandLineOptions.add_options()
            ("version,v", "print version information")
            ("help,h", "Get a list of options in the configuration file")
            ("config,c", po::value<std::string >(), "Set configuration file")
            ;

        po::variables_map vm_cmd;
        po::store(po::parse_command_line(argc, argv, commandLineOptions), vm_cmd);

        if (vm_cmd.empty()){
            if (world.rank() == 0) {
                std::cout << " To run MantisSA specify the configuration file as" << std::endl;
                std::cout << "-c config" << std::endl << std::endl;;
                std::cout << "Other command line options are:" << std::endl;
                std::cout << commandLineOptions << std::endl;
            }
            return false;
        }

        if (vm_cmd.count("version")) {
            if (world.rank() == 0) {
                std::cout << "|------------------|" << std::endl;
                std::cout << "|     MantisSA     |" << std::endl;
                std::cout << "| Version : " << Version <<" |" << std::endl;
                std::cout << "|  gwt.ucdavis.edu |" << std::endl;
                std::cout << "|------------------|" << std::endl;
            }
            return false;
        }

        // Configuration file options
        po::options_description config_options("Configuration file options");
        config_options.add_options()
            ("Raster.Ncells", po::value<int>()->default_value(0), "Number of cells with index")
            ("Raster.Nrows", po::value<int>()->default_value(0), "Number of Raster rows")
            ("Raster.Ncols", po::value<int>()->default_value(0), "Number of Raster Columns")
            ("Raster.Xorig", po::value<double>()->default_value(0), "X coordinate of left lower raster corner")
            ("Raster.Yorig", po::value<double>()->default_value(0), "Y coordinate of left lower raster corner")
            ("Raster.CellSize", po::value<double>()->default_value(0), "Cell size of raster")
            ("Raster.File", po::value<std::string>(), "Name of the file with the raster values")

            ("Simulation.Nyears", po::value<int>()->default_value(50), "Number of simulation years")
            ("Simulation.Porosity", po::value<int>()->default_value(2), "Porosity: 1-6")

            ("Simulation.SWAT_Data", po::value<std::string>(), "Swat input file")
            ("Simulation.NPSAT_VI", po::value<std::string>(), "NPSAT input file for VI")
            ("Simulation.NPSAT_VD", po::value<std::string>(), "NPSAT input file for VI")
            ("Simulation.InitSaltVI", po::value<std::string>(), "Initial salt concentration of VI wells")
            ("Simulation.InitSaltVD", po::value<std::string>(), "Initial salt concentration of VI wells")
            ("Simulation.DistribPump", po::value<std::string>(), "Distribution of pumping to Land")

            ("UNSAT.Depth_file", po::value<std::string>(), "Depth file name")
            ("UNSAT.Depth_name", po::value<std::string>(), "Depth scenario name")
            ("UNSAT.Rch_file", po::value<std::string>(), "Recharge file name")
            ("UNSAT.WC", po::value<double>()->default_value(0.2), "Water content")
            ("UNSAT.minDepth", po::value<double>()->default_value(1.0), "Minimum depth")
            ("UNSAT.minRch", po::value<double>()->default_value(10.0), "Minimum recharge")

            ("Other.Version", po::value<std::string>(), "version number")
        ;

        po::variables_map vm_cfg;
        if (vm_cmd.count("config")){
            try {
                configFile = vm_cmd["config"].as<std::string>().c_str();
                std::cout << "--> Configuration file: " << vm_cmd["config"].as<std::string>().c_str() << std::endl;
                po::store(po::parse_config_file<char>(vm_cmd["config"].as<std::string>().c_str(), config_options), vm_cfg);

                if (!vm_cfg.count("Other.Version")){
                    std::cout << "Version is not specified!. The current version is " << Version << std::endl;
                    return false;
                }

                std::string user_version = vm_cfg["Other.Version"].as<std::string>();
                if (Version.compare(user_version) != 0){
                    std::cout << "Mismatch between code version (" << Version << ") and Input file version (" << user_version << ")" << std::endl;
                    return false;
                }
                swat_input_file = vm_cfg["Simulation.SWAT_Data"].as<std::string>();
                npsat_VI_file = vm_cfg["Simulation.NPSAT_VI"].as<std::string>();
                npsat_VD_file = vm_cfg["Simulation.NPSAT_VD"].as<std::string>();
                init_salt_VI_file = vm_cfg["Simulation.InitSaltVI"].as<std::string>();
                init_salt_VD_file = vm_cfg["Simulation.InitSaltVD"].as<std::string>();
                cell_well_file = vm_cfg["Simulation.DistribPump"].as<std::string>();

                {// Read Raster
                    rasteroptions.Ncells = vm_cfg["Raster.Ncells"].as<int>();
                    rasteroptions.Nrows = vm_cfg["Raster.Nrows"].as<int>();
                    rasteroptions.Ncols = vm_cfg["Raster.Ncols"].as<int>();
                    rasteroptions.Xorig = vm_cfg["Raster.Xorig"].as<double>();
                    rasteroptions.Yorig = vm_cfg["Raster.Yorig"].as<double>();
                    rasteroptions.CellSize = vm_cfg["Raster.CellSize"].as<double>();
                    rasteroptions.File = vm_cfg["Raster.File"].as<std::string>();
                }

                NsimYears = vm_cfg["Simulation.Nyears"].as<int>();
                porosity = vm_cfg["Simulation.Porosity"].as<int>();

                depth_input_file = vm_cfg["UNSAT.Depth_file"].as<std::string>();
                depth_name = vm_cfg["UNSAT.Depth_name"].as<std::string>();
                rch_input_file = vm_cfg["UNSAT.Rch_file"].as<std::string>();
                wc = vm_cfg["UNSAT.WC"].as<double>();
                minDepth = vm_cfg["UNSAT.minDepth"].as<double>();
                minRch = vm_cfg["UNSAT.minRch"].as<double>();


            }
            catch (std::exception& E)
            {
                std::cout << E.what() << std::endl;
                return false;
            }
        }
        return true;
    }

}

#endif //MANTISSA_MS_INPUT_H
