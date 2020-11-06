#pragma once

#include <iostream>
#include <string>
#include <boost/program_options.hpp>


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
	};

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
			std::cout << "| Version : 1.4.00 |" << std::endl;
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
			if (vm_cmd.count("test")) {
				opt.testMode = true;
			}
			else {
				opt.testMode = false;
			}
			bool tf;
			std::cout << vm_cmd["config"].as<std::string>().c_str() << std::endl;
			po::store(po::parse_config_file<char>(vm_cmd["config"].as<std::string>().c_str(), config_options), vm_cfg);
			// read mandatory options
			opt.gnlmNpixels = vm_cfg["GNLM_Npixels"].as<int>();
			tf = get_option<std::string>("MAPS", vm_cfg, opt.MAPSfile);
			tf = get_option<std::string>("NO3_LOAD", vm_cfg, opt.NO3LoadFile);
			tf = get_option<std::string>("UNSAT", vm_cfg, opt.UNSATfile);
			tf = get_option<std::string>("WELLS", vm_cfg, opt.WELLfile);
			tf = get_option<std::string>("URFS", vm_cfg, opt.URFfile);

			// read optional options
			opt.port = vm_cfg["PORT"].as<int>();
			//opt.Nrow = vm_cfg["Nrow"].as<int>();
			//opt.Ncol= vm_cfg["Ncol"].as<int>();
			//opt.yearInterval = vm_cfg["YRINTERVAL"].as<int>();
			//opt.startYear = vm_cfg["StartYR"].as<int>();
			//opt.nSimulationYears = vm_cfg["NYRS"].as<int>();
			opt.nThreads = vm_cfg["NTHREADS"].as<int>();
		}
		return true;
	}

}
