//
// Created by giorg on 10/16/2024.
//

#ifndef MANTISSA_MS_INPUT_H
#define MANTISSA_MS_INPUT_H

#include <boost/mpi.hpp>
#include <boost/program_options.hpp>
#include "MS_structures.h"
#include "SWAT_data.h"

namespace po = boost::program_options;

namespace MS {



    class UserInput{
    public:
        UserInput(boost::mpi::communicator& world_in);
        bool read(int argc, char* argv[]);

        int NsimYears;
        int NswatYears;
        int porosity;
        int Nwells;
        int dbg_id;
        int nBuffer;
        int SWAT_data_version;
        bool doDebug = false;
        bool EnableFeedback = false;
        bool bUseInitConc4OutofArea = false;
        //int nurfsVI;
        //int nurfsVD;
        //int nSelectWells;
        double wc;
        double minDepth;
        double minRch;
        double maxConc;
        double maxAge;
        double SurfConcValue;
        double riverConcValue;
        double OutofAreaConc;
        std::string configFile;
        std::string swat_input_file;
        std::string HRU_raster_file;
        std::string hru_idx_file;
        std::string npsat_VI_file;
        std::string npsat_VD_file;
        std::string init_salt_VI_file;
        std::string init_salt_VD_file;
        std::string cell_well_file;
        std::string depth_input_file;
        std::string depth_name;
        std::string rch_input_file;
        std::string SelectedWells_file;
        std::string SelectedWellsGroupFile;
        std::string dbg_file;

        RasterOptions rasteroptions;

        std::string outfile;
        //std::string outfileVD;
        //std::string outfileVIdetail;
        //std::string outfileVDdetail;
        //std::string outfileVImfeed;
        //std::string outfileVDmfeed;
        //std::string outfileVIurfs;
        //std::string outfileVDurfs;
        //std::string outfileVIbtcs;
        //std::string outfileVDbtcs;
        std::string Version;
        bool printURFs;
        bool printLoad;
        bool printBTCs;
        bool bUseInitConcVI;
        bool bUseInitConcVD;

    private:
        boost::mpi::communicator world;
    };

    UserInput::UserInput(boost::mpi::communicator &world_in)
        :
            world(world_in)
    {
        Version = "0.0.27";
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
            ("Simulation.NSwatYears", po::value<int>()->default_value(30), "Number of simulation year in swat data")
            ("Simulation.Porosity", po::value<int>()->default_value(2), "Porosity: 1-6")
            ("Simulation.MaxConc", po::value<double>()->default_value(10000), "Maximum loading concentration (before convolution)")
            ("Simulation.SurfConcValue", po::value<double>()->default_value(-9), "Surface water salt concentration. Negative deactivates")
            ("Simulation.RiverConcValue", po::value<double>()->default_value(-9), "River salt concentration. Negative deactivates")
            ("Simulation.OutofAreaConc", po::value<double>()->default_value(-9), "Concentration for streamlines outside the study area. Negative skip the streamlines")
            ("Simulation.OutofAreaUseInitConc", po::value<int>()->default_value(1), "Use the initial concentration if outside of area.")
            ("Simulation.MaxAge", po::value<double>()->default_value(-9), "River salt concentration. Negative deactivates")
            ("Simulation.EnableFeedback", po::value<int>()->default_value(0), "Set this to non zero to run with feedback")
            ("Simulation.nBuffer", po::value<int>()->default_value(1), "Number of buffer zones around each streamline")

            ("Simulation.SWAT_data_version", po::value<int>()->default_value(1), "Swat data version")
            ("Simulation.SWAT_Data", po::value<std::string>(), "Swat input file")
            ("Simulation.HRU_Raster", po::value<std::string>(), "HRU raster file")
            ("Simulation.SWAT_HRUs", po::value<std::string>(), "Swat HRU index map")
            ("Simulation.NPSAT_VI", po::value<std::string>(), "NPSAT input file for VI")
            //("Simulation.NURFS_VI", po::value<int>(), "Number of VI streamlines")
            //("Simulation.NURFS_VD", po::value<int>(), "Number of VD streamlines")
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

            ("Other.OutFile", po::value<std::string>(), "Output filename")
            ("Other.SelectedWells", po::value<std::string>(), "Selected wells for detailed output")
            ("Other.SelectedWellsGroups", po::value<std::string>(), "A file with group Ids and names")
            //("Other.NselectWells", po::value<int>(), "Number of Selected wells")
            ("Other.dbg_File", po::value<std::string>(), "Debugging output filename")
            ("Other.dbg_ids", po::value<int>(), "well ids for debugging")
            ("Other.PrintMatrices", po::value<int>(), "Print Matrices to test communication works")
            ("Other.printLoad", po::value<int>(), "prints details of loading functions")
            ("Other.printURFs", po::value<int>(), "prints the urfs of the selected wells")
            ("Other.printBTCs", po::value<int>(), "prints the BTCs for each streamline of the selected wells")
            ("Other.Version", po::value<std::string>(), "version number")
        ;

        if (vm_cmd.count("help")) {
            if (world.rank() == 0) {
                std::cout << " To run ICHNOS specify the configuration file as" << std::endl;
                std::cout << "-c config" << std::endl << std::endl;;
                std::cout << "Other command line options are:" << std::endl;
                std::cout << commandLineOptions << std::endl;

                std::cout << "ICHNOS configuration file options:" << std::endl;
                std::cout << "The options without default values are mandatory" << std::endl;
                std::cout << "(All options are case sensitive)" << std::endl;
                std::cout << "------------------------------" << std::endl;
                std::cout << config_options << std::endl;
            }
            return false;
        }

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
                HRU_raster_file = vm_cfg["Simulation.HRU_Raster"].as<std::string>();
                SWAT_data_version = vm_cfg["Simulation.SWAT_data_version"].as<int>();
                swat_input_file = vm_cfg["Simulation.SWAT_Data"].as<std::string>();
                hru_idx_file = vm_cfg["Simulation.SWAT_HRUs"].as<std::string>();
                npsat_VI_file = vm_cfg["Simulation.NPSAT_VI"].as<std::string>();
                //nurfsVI = vm_cfg["Simulation.NURFS_VI"].as<int>();
                npsat_VD_file = vm_cfg["Simulation.NPSAT_VD"].as<std::string>();
                //nurfsVD = vm_cfg["Simulation.NURFS_VD"].as<int>();
                init_salt_VI_file = vm_cfg["Simulation.InitSaltVI"].as<std::string>();
                if (init_salt_VI_file.empty()) {
                    bUseInitConcVI = false;
                }
                else{
                    bUseInitConcVI = true;
                }
                init_salt_VD_file = vm_cfg["Simulation.InitSaltVD"].as<std::string>();
                if (init_salt_VD_file.empty()) {
                    bUseInitConcVD = false;
                }
                else{
                    bUseInitConcVD = true;
                }
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
                NswatYears = vm_cfg["Simulation.NSwatYears"].as<int>();
                porosity = vm_cfg["Simulation.Porosity"].as<int>();
                maxConc = vm_cfg["Simulation.MaxConc"].as<double>();
                maxAge = vm_cfg["Simulation.MaxAge"].as<double>();
                SurfConcValue = vm_cfg["Simulation.SurfConcValue"].as<double>();
                riverConcValue = vm_cfg["Simulation.RiverConcValue"].as<double>();
                OutofAreaConc = vm_cfg["Simulation.OutofAreaConc"].as<double>();
                bUseInitConc4OutofArea = vm_cfg["Simulation.OutofAreaUseInitConc"].as<int>()!= 0;
                EnableFeedback = vm_cfg["Simulation.EnableFeedback"].as<int>() != 0;
                nBuffer = vm_cfg["Simulation.nBuffer"].as<int>();

                depth_input_file = vm_cfg["UNSAT.Depth_file"].as<std::string>();
                depth_name = vm_cfg["UNSAT.Depth_name"].as<std::string>();
                rch_input_file = vm_cfg["UNSAT.Rch_file"].as<std::string>();
                wc = vm_cfg["UNSAT.WC"].as<double>();
                minDepth = vm_cfg["UNSAT.minDepth"].as<double>();
                minRch = vm_cfg["UNSAT.minRch"].as<double>();

                outfile = vm_cfg["Other.OutFile"].as<std::string>();

                //outfileVI = mainOutfile + "VI.dat";
                //outfileVD = mainOutfile + "VD.dat";

                SelectedWells_file = vm_cfg["Other.SelectedWells"].as<std::string>();
                SelectedWellsGroupFile = vm_cfg["Other.SelectedWellsGroups"].as<std::string>();
                //nSelectWells = vm_cfg["Other.NselectWells"].as<int>();

                //outfileVIdetail = mainOutfile + "VI_lf.dat";
                //outfileVDdetail = mainOutfile + "VD_lf.dat";
                //outfileVImfeed = mainOutfile + "VI_mf.dat";
                //outfileVDmfeed = mainOutfile + "VD_mf.dat";
                //outfileVIurfs = mainOutfile + "VI_urf.dat";
                //outfileVDurfs = mainOutfile + "VD_urf.dat";
                //outfileVIbtcs = mainOutfile + "VI_btc.dat";
                //outfileVDbtcs = mainOutfile + "VD_btc.dat";
                dbg_file = vm_cfg["Other.dbg_File"].as<std::string>();
                dbg_id = vm_cfg["Other.dbg_ids"].as<int>();
                PrintMatrices = vm_cfg["Other.PrintMatrices"].as<int>() != 0;
                doDebug = !dbg_file.empty();
                printLoad = vm_cfg["Other.printLoad"].as<int>() != 0;
                printURFs = vm_cfg["Other.printURFs"].as<int>() != 0;
                printBTCs = vm_cfg["Other.printBTCs"].as<int>() != 0;
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
