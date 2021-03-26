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
		//! The number of pixels in \p LUfile and \p NGWfile files
		int gnlmNpixels;

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
		std::string NO3LoadFile;

		//! GNLM_NGWfile holds the Nitrate groundwater load according to GNLM.
		/*!
		This file has the following format:

		ID linear index NGW45 NGW60 NGW75 NGW90 NGW05 NGW20 NGW35 NGW50

		where NGWxx are the nitrate groundwater load. The pixels must be in the same order as in the
		\p LUfile. 
		*/
		//std::string GNLM_NGWfile;

		//std::string SWAT_Mainfile;

		std::string UNSATfile;

		//! MAPSfile is the name file for the background maps
		/*!
		The format of the file is....
		*/
		std::string MAPSfile;

		//! Port is the port number
		/*!
		This is the port that the clients use to communicate with the server
		By defaults it gets 1234
		*/
		int port;

		std::string WELLfile;
		std::string URFfile;

		//int yearInterval;
		//int startYear;
		//int nSimulationYears;

		int nThreads;

		bool testMode;

        /**
         * Indicates whether the paths inside the configuration file are absolute or relative to the main path
         */
		bool bAbsolutePaths = true;

		/**
		 * This is the main path where all input files are relative to.
		 */
		std::string mainPath;
		/**
		 * The name of the configuration file
		 */
		std::string configFile;

		std::string DebugFolder;

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
			std::cout << "mantisServer -c configfilename" << std::endl << std::endl;;
			std::cout << "Other command line options are:" << std::endl;
			std::cout << commandLineOptions << std::endl;
			return false;
		}

		if (vm_cmd.count("version")) {
			std::cout << "|------------------|" << std::endl;
			std::cout << "|  Mantis Server   |" << std::endl;
			std::cout << "| Version : 1.4.11 |" << std::endl;
			std::cout << "|    by  giorgk    |" << std::endl;
			std::cout << "|------------------|" << std::endl;
			return false;
		}

		// Configuration file options
		po::options_description config_options("Configuration file options");
		config_options.add_options()
			("GNLM_Npixels", po::value<int>(), "Number of pixel in GNLM load. It is used for allocation")
			("MAPS", "Geometry of background maps")
			("NO3_LOAD", "A File with a list of files that contain the loading information")
			("WELLS", "A list of files with the well info for each scenario")
			("URFS", "A list of files with the URF information")
			("UNSAT", "A file that containts the travel time for each LU pixel")
            ("DataPath", "If this is not empty all data all data must be relative to DataPath.")
			("PORT", po::value<int>()->default_value(1234), "Port number")
			("NTHREADS", po::value<int>()->default_value(6), "Number of threads to use by server")
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

			bool tf;
			std::cout << opt.configFile << std::endl;
			po::store(po::parse_config_file<char>(opt.configFile.c_str(), config_options), vm_cfg);
			// read mandatory options
			opt.gnlmNpixels = vm_cfg["GNLM_Npixels"].as<int>();
			if (!get_option<std::string>("MAPS", vm_cfg, opt.MAPSfile)) return false;
			if (!get_option<std::string>("NO3_LOAD", vm_cfg, opt.NO3LoadFile)) return false;
			if (!get_option<std::string>("UNSAT", vm_cfg, opt.UNSATfile)) return false;
			if (!get_option<std::string>("WELLS", vm_cfg, opt.WELLfile)) return false;
			if (!get_option<std::string>("URFS", vm_cfg, opt.URFfile)) return false;

            opt.mainPath = "";
			if (vm_cfg.count("DataPath")){
			    if (!get_option<std::string>("DataPath", vm_cfg, opt.mainPath)){
                    opt.bAbsolutePaths = true;
			    }
			    else{
                    opt.bAbsolutePaths = false;
			    }
			}
			else{
			    opt.bAbsolutePaths = true;
			}

			// read optional options
			opt.port = vm_cfg["PORT"].as<int>();

			opt.nThreads = vm_cfg["NTHREADS"].as<int>();

			if (vm_cfg.count("DEBUG_DIR")) {
				if (!get_option<std::string>("DEBUG_DIR", vm_cfg, opt.DebugFolder));
                    opt.DebugFolder = "";
			}
			else{
                opt.DebugFolder = "";
			}
            return true;
		}
        std::cout << "Mantis received wrong input options" << std::endl;
		std::cout << "Run mantisServer -h for additional help" << std::endl;
        return false;
	}
}
