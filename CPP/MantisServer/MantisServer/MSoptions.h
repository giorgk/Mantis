#pragma once

#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>


namespace po = boost::program_options;

namespace mantisServer {

	template <typename T> 
	bool get_option(std::string optionName, po::variables_map &vm, T &out) {

		if (vm.count(optionName)) {
			out = vm[optionName].as<T>();
			return true;
		}
		else {
			std::cout << optionName << " option is missing" << std::endl;
			return false;
		}
	}

	struct options {
        std::string RegionsListFile;
        int Nregions = 0;
		//! The number of pixels in \p LUfile and \p NGWfile files
		//int Npixels;

		//! The number of rows of the Land use and Nitrate groundwater loading maps
		//int Nrow;

		//! The number of columns of the Land use and Nitrate groundwater loading maps
		//int Ncol;

		//! GNLM_LUfile holds name file for the Land use file.
		/*!
		This file has the following format:

		ID linear_index YR45 YR60 YR75 YR90 YR05

		where ID, is an id with each pixel linear index is the linear index of the pixel position.  
		YRxx are the land use codes. The file  contain the pixels
		that are non zero either in land use.
		*/
		//std::string NO3LoadFile;

		//std::string CVrasterFile;

		//! GNLM_NGWfile holds the Nitrate groundwater load according to GNLM.
		/*!
		This file has the following format:

		ID linear index NGW45 NGW60 NGW75 NGW90 NGW05 NGW20 NGW35 NGW50

		where NGWxx are the nitrate groundwater load. The pixels must be in the same order as in the
		\p LUfile. 
		*/
		//std::string GNLM_NGWfile;

		//std::string SWAT_Mainfile;

		//std::string UNSATfile;

		//! MAPSfile is the name file for the background maps
		/*!
		The format of the file is....
		*/
		//std::string MAPSfile;

		//! Port is the port number
		/*!
		This is the port that the clients use to communicate with the server
		By defaults it gets 1234
		*/
		int port;
        std::string ip;

		//std::string WELLfile;
		//std::string URFfile;
		//std::string RCHfile;
		//std::string RFURFfile;
		//bool bReadRFURF;

		//int yearInterval;
		//int startYear;
		//int nSimulationYears;

		int nThreads;

		bool testMode;

        /**
         * Indicates whether the paths inside the configuration file are absolute or relative to the main path
         */
		//bool bAbsolutePaths = true;

		/**
		 * This is the main path where all input files are relative to.
		 */
		std::string mainPath;
		/**
		 * The name of the configuration file
		 */
		std::string configFile;

		std::string DebugPrefix;

        /// This is the string input file name
        std::string logFile;
        bool bUseLogFile = false;
        int nLogReset = 10;
        int nTimesPrinted = 0;

		int RFmem;
        std::string version = "2.0.03";

	};

	/**
	 * This is the function that parses the input arguments
	 *
	 * There are options that can be defined at command line or inside the configuration file
	 *
	 * Command line options are the following:
	 * - v or version Prints the version and exists
	 * - h or help Prints a help message and exits
	 * - c or config. A file name must be entered after this option with the name of the configuration file
	 * - p or path. A path must be enter after this option. \n
	 *              This path indicates that the configuration file and all the files inside the configuration file
	 *              are set relative to this path. If the paths of the configuration file and the names of the files inside
	 *              the configuration file are absolute then this path must be empty.
	 * - t This runs mantis in test mode. The output messages are just random numbers which have the exact format
	 *     the actual runs.
	 *
	 *
	 *
	 * @param argc
	 * @param argv
	 * @param opt
	 * @return true if parsing was successful.
	 */
	bool readInputParameters(int argc, char *argv[], options& opt) {
		// user command line options
		po::options_description commandLineOptions("Command line options");
		commandLineOptions.add_options()
			("version,v", "print version information")
			("help,h", "Get a list of options in the configuration file")
			("config,c", po::value<std::string >(), "Set configuration file")
			("test,t", "Run Server in test mode [a config file is required]")
			;

		po::variables_map vm_cmd;
		po::store(po::parse_command_line(argc, argv, commandLineOptions), vm_cmd);

		if (vm_cmd.size() == 0) {
			std::cout << " To run MantisServer specify the configuration file and run as" << std::endl;
			std::cout << "mantisServer -c configfilename" << std::endl << std::endl;
			std::cout << "Other command line options are:" << std::endl;
			std::cout << commandLineOptions << std::endl;
			return false;
		}

		if (vm_cmd.count("version")) {
			std::cout << "|------------------|" << std::endl;
			std::cout << "|  Mantis Server   |" << std::endl;
			std::cout << "| Version : " << opt.version << " |" << std::endl;
			std::cout << "|    by  giorgk    |" << std::endl;
			std::cout << "|------------------|" << std::endl;
			return false;
		}

		// Configuration file options
		po::options_description config_options("Configuration file options");
		config_options.add_options()
            //Regions
            ("Regions.File", "The file with a list of regions")
		    // CV_Raster
            //("CV_Raster.Ncells", po::value<int>(), "Number of active pixels in background raster")
            //("CV_Raster.Nrows", po::value<int>(), "Number of rows of the background raster")
            //("CV_Raster.Ncols", po::value<int>(), "Number of columns of the background raster")
            //("CV_Raster.Raster", "The file with the raster indices that point to Nidx")

            // Data
			//("Data.MAPS", "A List of the background maps with their subregions")
			//("Data.NO3", "A file with a list of files that contain the Nitrate loading base scenarios")
            //("Data.UNSAT", "A file that contains the travel time for each LU pixel")
			//("Data.WELLS", "A list of files with the well info for each scenario")
			//("Data.URFS", "A list of files with the URF information")
			//("Data.RFURF", "A list of files with Region and Flow specific URFs")
            //("Data.RCH", "A file that contains the recharge values in mm/year for each LU pixel")
            //("Data.Path", "If this is not empty all data all data must be relative to this Path.")

			// ServerOptions
            ("ServerOptions.PORT", po::value<int>()->default_value(1234), "Port number")
            ("ServerOptions.IP", "IP address as string. Leave empty for default 127.0.0.1")
            ("ServerOptions.NTHREADS", po::value<int>()->default_value(6), "Number of threads to use by server")
            ("ServerOptions.RFmem", po::value<int>()->default_value(5), "RF Memory depth")
            ("ServerOptions.logFile", "If it is empty the output will be printed in the terminal. You can include a full path")
            ("ServerOptions.logClearFreq", po::value<int>()->default_value(50), "Frequency to clear log file")
            ("ServerOptions.DebugPrefix", "Any output file will start with this prefix. You can include a full path")
			;

		if (vm_cmd.count("help")) {
			std::cout << "------------------------------" << std::endl;
			std::cout << commandLineOptions << std::endl;
			std::cout << "------------------------------" << std::endl;
			std::cout << "Copy paste the following list of options" << std::endl;
			std::cout << "into a configuration file." << std::endl;
			std::cout << "The options without default values are mandatory" << std::endl;
			std::cout << "All options are case sensitive" << std::endl;
			std::cout << "------------------------------" << std::endl;
			std::cout << config_options << std::endl;
			return false;
		}

		po::variables_map vm_cfg;

		if (vm_cmd.count("config")) {
            opt.configFile = vm_cmd["config"].as<std::string>();
			if (vm_cmd.count("test")) {
			    std::cout << "Mantis is set to run in test mode!. The results are random numbers!!" << std::endl;
				opt.testMode = true;
			}
			else {
				opt.testMode = false;
			}

			std::cout << opt.configFile << std::endl;
			po::store(po::parse_config_file<char>(opt.configFile.c_str(), config_options), vm_cfg);
            //Regions
            if (!get_option<std::string>("Regions.File", vm_cfg, opt.RegionsListFile)) return false;

            // CV_Raster
			//opt.Npixels = vm_cfg["CV_Raster.Ncells"].as<int>();
            //opt.Nrow = vm_cfg["CV_Raster.Nrows"].as<int>();
            //opt.Ncol = vm_cfg["CV_Raster.Ncols"].as<int>();
            //if (!get_option<std::string>("CV_Raster.Raster", vm_cfg, opt.CVrasterFile)) return false;

            // Data
			//if (!get_option<std::string>("Data.MAPS", vm_cfg, opt.MAPSfile)) return false;
			//if (!get_option<std::string>("Data.NO3", vm_cfg, opt.NO3LoadFile)) return false;
			//if (!get_option<std::string>("Data.UNSAT", vm_cfg, opt.UNSATfile)) return false;
			//if (!get_option<std::string>("Data.WELLS", vm_cfg, opt.WELLfile)) return false;
			//if (!get_option<std::string>("Data.URFS", vm_cfg, opt.URFfile)) return false;
            //if (!get_option<std::string>("Data.RCH", vm_cfg, opt.RCHfile)) return false;

            //if (vm_cfg.count("Data.RFURF")){
			//	if (get_option<std::string>("Data.RFURF", vm_cfg, opt.RFURFfile)){
			//		opt.bReadRFURF = true;
			//	}
			//	else{
			//		opt.bReadRFURF = false;
			//	}
            //}

            //opt.mainPath = "";
			//if (vm_cfg.count("Data.Path")){
			//    if (get_option<std::string>("Data.Path", vm_cfg, opt.mainPath)){
			//    	if (opt.mainPath.empty()){
			//			opt.bAbsolutePaths = true;
			//    	}
			//    	else{
			//			opt.bAbsolutePaths = false;
            //   	}
			//    }
			//}
			//else{
			//    opt.bAbsolutePaths = true;
			//}

            // ServerOptions
			opt.port = vm_cfg["ServerOptions.PORT"].as<int>();
            if (!get_option<std::string>("ServerOptions.IP", vm_cfg, opt.ip)){
                opt.ip = "127.0.0.1";
            }

			opt.nThreads = vm_cfg["ServerOptions.NTHREADS"].as<int>();
            opt.RFmem = vm_cfg["ServerOptions.RFmem"].as<int>();

			if (vm_cfg.count("ServerOptions.DebugPrefix")) {
				if (!get_option<std::string>("ServerOptions.DebugPrefix", vm_cfg, opt.DebugPrefix))
                    opt.DebugPrefix = "";
			}
			else{
                opt.DebugPrefix = "";
			}

            opt.nLogReset = vm_cfg["ServerOptions.logClearFreq"].as<int>();
            if (vm_cfg.count("ServerOptions.logFile")){
                if (get_option<std::string>("ServerOptions.logFile", vm_cfg, opt.logFile)){
                    if (!opt.logFile.empty()){
                        opt.bUseLogFile = true;
                        logStream.open(opt.logFile.c_str(), std::ios::out);
                        std::cout.rdbuf(logStream.rdbuf());
                    }
                }
            }

            return true;
		}
        std::cout << "Mantis received wrong input options" << std::endl;
		std::cout << "Run mantisServer -h for additional help" << std::endl;
        return false;
	}
}
