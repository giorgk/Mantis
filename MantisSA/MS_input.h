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

        std::string configFile;
        int dbg_id;
        bool doDebug = false;
        std::string dbg_file;

        RasterOptions rasterOptions;
        SwatOptions swatOptions;
        HistoricOptions historicOptions;
        NpsatOptions npsatOptions;
        RiverOptions riverOptions;
        SimOptions simOptions;
        UnsatOptions unsatOptions;
        OutputOptions outputOptions;

        std::string Version;

    private:
        boost::mpi::communicator world;
    };

    UserInput::UserInput(boost::mpi::communicator &world_in)
        :
            world(world_in)
    {
        Version = "0.0.30";
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
            // [Raster]
            ("Raster.Ncells", po::value<int>()->default_value(0), "Number of cells with index")
            ("Raster.Nrows", po::value<int>()->default_value(0), "Number of Raster rows")
            ("Raster.Ncols", po::value<int>()->default_value(0), "Number of Raster Columns")
            ("Raster.Xorig", po::value<double>()->default_value(0), "X coordinate of left lower raster corner")
            ("Raster.Yorig", po::value<double>()->default_value(0), "Y coordinate of left lower raster corner")
            ("Raster.CellSize", po::value<double>()->default_value(0), "Cell size of raster")
            ("Raster.File", po::value<std::string>(), "Name of the file with the raster values")

            // [SWAT]
            ("SWAT.Nyears", po::value<int>()->default_value(30), "Number of simulation year in swat data")
            ("SWAT.StartYear", po::value<int>()->default_value(1995), "Number of simulation year in swat data")
            ("SWAT.Data_version", po::value<int>()->default_value(1), "Swat data version")
            ("SWAT.HRU_Raster", po::value<std::string>(), "HRU raster file")
            ("SWAT.HRU_index", po::value<std::string>(), "Swat HRU index map")
            ("SWAT.Data", po::value<std::string>(), "Swat input file")

            // [Historic]
            ("Historic.StartYear", po::value<int>()->default_value(1945), "Start year of historic loading")
            ("Historic.EndYear", po::value<int>()->default_value(2005), "End year of historic loading")
            ("Historic.Interval", po::value<int>()->default_value(2005), "Interval of historic loading")
            ("Historic.Prefix_name", po::value<std::string>(), "Prefix of historic loading")
            ("Historic.Ext", po::value<std::string>(), "Extension of historic loading without the dot e.g h5 or dat")
            ("Historic.BlendStart", po::value<int>()->default_value(1945), "Start year of historic loading")
            ("Historic.BlendEnd", po::value<int>()->default_value(2005), "End year of historic loading")

            // [NPSAT]
            ("NPSAT.Data_version", po::value<int>()->default_value(1), "NPSAT data version")
            ("NPSAT.VIwells", po::value<std::string>(), "NPSAT input file for VI")
            ("NPSAT.VDwells", po::value<std::string>(), "NPSAT input file for VI")
            ("NPSAT.InitSaltVI", po::value<std::string>(), "Initial salt concentration of VI wells")
            ("NPSAT.InitSaltVD", po::value<std::string>(), "Initial salt concentration of VI wells")
            ("NPSAT.DistribPump", po::value<std::string>(), "Distribution of pumping to Land")

            // [Simulation]
            ("Simulation.Nyears", po::value<int>()->default_value(50), "Number of simulation years")
            ("Simulation.StartYear", po::value<int>()->default_value(1945), "Starting year of simulation")
            ("Simulation.Porosity", po::value<int>()->default_value(2), "Porosity: 1-6")
            ("Simulation.MaxConc", po::value<double>()->default_value(10000), "Maximum loading concentration (before convolution)")
            ("Simulation.MaxAge", po::value<double>()->default_value(-9), "River salt concentration. Negative deactivates")
            ("Simulation.SurfConcValue", po::value<double>()->default_value(-9), "Surface water salt concentration. Negative deactivates")
            ("Simulation.ConcValue", po::value<double>()->default_value(-9), "River salt concentration. Negative deactivates")
            ("Simulation.StartDist", po::value<double>()->default_value(100), "Up to this distance the river influence is 1")
            ("Simulation.EndDist", po::value<double>()->default_value(300), "After this distance the river influence is 0")
            ("Simulation.OutofAreaConc", po::value<double>()->default_value(-9), "Concentration for streamlines outside the study area. Negative skip the streamlines")
            ("Simulation.OutofAreaUseInitConc", po::value<int>()->default_value(1), "Use the initial concentration if outside of area.")
            ("Simulation.EnableFeedback", po::value<int>()->default_value(0), "Set this to non zero to run with feedback")
            ("Simulation.nBuffer", po::value<int>()->default_value(1), "Number of buffer zones around each streamline")

            // [UNSAT]
            ("UNSAT.Depth_file", po::value<std::string>(), "Depth file name")
            ("UNSAT.Depth_name", po::value<std::string>(), "Depth scenario name")
            ("UNSAT.Rch_file", po::value<std::string>(), "Recharge file name")
            ("UNSAT.WC", po::value<double>()->default_value(0.2), "Water content")
            ("UNSAT.minDepth", po::value<double>()->default_value(1.0), "Minimum depth")
            ("UNSAT.minRch", po::value<double>()->default_value(10.0), "Minimum recharge")

            // [Output]
            ("Output.FilePrefix", po::value<std::string>(), "Output filename")
            ("Output.SelectedWells", po::value<std::string>(), "Selected wells for detailed output")
            ("Output.SelectedWellsGroups", po::value<std::string>(), "A file with group Ids and names")
            ("Output.printLoad", po::value<int>(), "prints details of loading functions")
            ("Output.printURFs", po::value<int>(), "prints the urfs of the selected wells")
            ("Output.printBTCs", po::value<int>(), "prints the BTCs for each streamline of the selected wells")

            // [Other]
            ("Other.dbg_File", po::value<std::string>(), "Debugging output filename")
            ("Other.dbg_ids", po::value<int>(), "well ids for debugging")
            ("Other.PrintMatrices", po::value<int>(), "Print Matrices to test communication works")
            ("Other.Version", po::value<std::string>(), "version number")
        ;

        if (vm_cmd.count("help")) {
            if (world.rank() == 0) {
                std::cout << " To run MantisSA specify the configuration file as" << std::endl;
                std::cout << "-c config" << std::endl << std::endl;;
                std::cout << "Other command line options are:" << std::endl;
                std::cout << commandLineOptions << std::endl;

                std::cout << "MantisSA configuration file options:" << std::endl;
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

                {// [Raster]
                    rasterOptions.Ncells = vm_cfg["Raster.Ncells"].as<int>();
                    rasterOptions.Nrows = vm_cfg["Raster.Nrows"].as<int>();
                    rasterOptions.Ncols = vm_cfg["Raster.Ncols"].as<int>();
                    rasterOptions.Xorig = vm_cfg["Raster.Xorig"].as<double>();
                    rasterOptions.Yorig = vm_cfg["Raster.Yorig"].as<double>();
                    rasterOptions.CellSize = vm_cfg["Raster.CellSize"].as<double>();
                    rasterOptions.File = vm_cfg["Raster.File"].as<std::string>();
                }

                {// [SWAT]
                    swatOptions.Nyears = vm_cfg["SWAT.Nyears"].as<int>();
                    swatOptions.StartYear = vm_cfg["SWAT.StartYear"].as<int>();
                    swatOptions.version = vm_cfg["SWAT.Data_version"].as<int>();
                    swatOptions.HRU_raster_file = vm_cfg["SWAT.HRU_Raster"].as<std::string>();
                    swatOptions.HRU_index_file = vm_cfg["SWAT.HRU_index"].as<std::string>();
                    swatOptions.Data_file = vm_cfg["SWAT.Data"].as<std::string>();
                }

                {// [Historic]
                    historicOptions.StartYear = vm_cfg["Historic.StartYear"].as<int>();
                    historicOptions.EndYear = vm_cfg["Historic.EndYear"].as<int>();
                    historicOptions.Interval = vm_cfg["Historic.Interval"].as<int>();
                    historicOptions.filename = vm_cfg["Historic.Prefix_name"].as<std::string>();
                    historicOptions.ext = vm_cfg["Historic.Ext"].as<std::string>();
                    historicOptions.BlendStart = vm_cfg["Historic.BlendStart"].as<int>();
                    historicOptions.BlendEnd = vm_cfg["Historic.BlendEnd"].as<int>();
                }

                {// [NPSAT]
                    npsatOptions.version = vm_cfg["NPSAT.Data_version"].as<int>();
                    npsatOptions.VIdataFile = vm_cfg["NPSAT.VIwells"].as<std::string>();
                    npsatOptions.VDdataFile = vm_cfg["NPSAT.VDwells"].as<std::string>();
                    npsatOptions.InitSaltVIFile = vm_cfg["NPSAT.InitSaltVI"].as<std::string>();
                    npsatOptions.bUseInitConcVI =  !npsatOptions.InitSaltVIFile.empty();
                    npsatOptions.InitSaltVDFile = vm_cfg["NPSAT.InitSaltVD"].as<std::string>();
                    npsatOptions.bUseInitConcVD =  !npsatOptions.InitSaltVDFile.empty();
                    npsatOptions.DistribuPumpFile = vm_cfg["NPSAT.DistribPump"].as<std::string>();
                }

                {//[Simulation]
                    simOptions.StartYear = vm_cfg["Simulation.StartYear"].as<int>();
                    simOptions.Nyears = vm_cfg["Simulation.Nyears"].as<int>();
                    simOptions.Porosity = vm_cfg["Simulation.Porosity"].as<int>();
                    simOptions.MaxConc = vm_cfg["Simulation.MaxConc"].as<double>();
                    simOptions.MaxAge = vm_cfg["Simulation.MaxAge"].as<double>();
                    simOptions.SurfConcValue = vm_cfg["Simulation.SurfConcValue"].as<double>();
                    riverOptions.ConcValue = vm_cfg["Simulation.ConcValue"].as<double>();
                    riverOptions.StartDist = vm_cfg["Simulation.StartDist"].as<double>();
                    riverOptions.EndDist = vm_cfg["Simulation.EndDist"].as<double>();
                    simOptions.OutofAreaConc = vm_cfg["Simulation.OutofAreaConc"].as<double>();
                    simOptions.OutofAreaUseInitConc = vm_cfg["Simulation.OutofAreaUseInitConc"].as<int>()!= 0;
                    simOptions.EnableFeedback = vm_cfg["Simulation.EnableFeedback"].as<int>() != 0;
                    simOptions.nBuffer = vm_cfg["Simulation.nBuffer"].as<int>();

                }

                {// [UNSAT]
                    unsatOptions.Depth_file = vm_cfg["UNSAT.Depth_file"].as<std::string>();
                    unsatOptions.Depth_name = vm_cfg["UNSAT.Depth_name"].as<std::string>();
                    unsatOptions.Rch_file = vm_cfg["UNSAT.Rch_file"].as<std::string>();
                    unsatOptions.wc = vm_cfg["UNSAT.WC"].as<double>();
                    unsatOptions.minDepth = vm_cfg["UNSAT.minDepth"].as<double>();
                    unsatOptions.minRch = vm_cfg["UNSAT.minRch"].as<double>();
                }

                {// [Output]
                    outputOptions.OutFile = vm_cfg["Other.OutFile"].as<std::string>();
                    outputOptions.SelectedWells = vm_cfg["Other.SelectedWells"].as<std::string>();
                    outputOptions.SelectedWellGroups = vm_cfg["Other.SelectedWellsGroups"].as<std::string>();
                    outputOptions.printSelectedWells = !outputOptions.SelectedWells.empty() && !outputOptions.SelectedWellGroups.empty();
                    outputOptions.printBTCs = vm_cfg["Other.printBTCs"].as<int>() != 0;
                    outputOptions.printLoad = vm_cfg["Other.printLoad"].as<int>() != 0;
                    outputOptions.printURFs = vm_cfg["Other.printURFs"].as<int>() != 0;
                }

                dbg_file = vm_cfg["Other.dbg_File"].as<std::string>();
                dbg_id = vm_cfg["Other.dbg_ids"].as<int>();
                PrintMatrices = vm_cfg["Other.PrintMatrices"].as<int>() != 0;
                doDebug = !dbg_file.empty();
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
