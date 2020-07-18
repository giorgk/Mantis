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
		int Npixels;

		//! The number of rows of the Land use and Nitrate groundwater loading maps
		int Nrow;

		//! The number of columns of the Land use and Nitrate groundwater loading maps
		int Ncol;

		//! LUfile holds name file for the Land use file.
		/*!
		This file has the following format:

		I J YR45 YR60 YR75 YR90 YR05

		where I, J are the row and columns of the pixel and 
		YRxx are the land use codes. The file should contain the pixels
		that are non zero either in land used or nitrate groundwater load.
		*/
		std::string LUfile;

		//! NGWfile holds name file for the Nitrate groundwater load.
		/*!
		This file has the following format:

		NGW45 NGW60 NGW75 NGW90 NGW05 NGW20 NGW35 NGW50

		where NGWxx are the nitrate groundwater load. The pixels must be in the same order as in the
		\p LUfile. 
		*/
		std::string NGWfile;

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

		//! the number of available scenarios
		int Nscenarios;

		std::string WELLfile;
		std::string URFfile;

		int yearInterval;
		int startYear;
		int nSimulationYears;

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
			std::cout << "| Version : 1.3.00 |" << std::endl;
			std::cout << "|    by  giorgk    |" << std::endl;
			std::cout << "|------------------|" << std::endl;
			return false;
		}

		// Configuration file options
		po::options_description config_options("Configuration file options");
		config_options.add_options()
			("Npixels", po::value<int>(), "LU and NGW pixels in files")
			("MAPS", "Geometry of background maps")
			("LU", "Land Use file")
			("NGW", "Nitrate groundwater loading file")
			("NSCEN", po::value<int>()->default_value(1), "Number of scenarios")
			("WELLS", "A list of files with the well info for each scenario")
			("URFS", "A list of files with the URF information")
			("UNSAT", "A file that containts the travel time for each LU pixel")
			("PORT", po::value<int>()->default_value(1234), "Port number")
			("Nrow", po::value<int>()->default_value(12863), "Number of rows")
			("Ncol", po::value<int>()->default_value(7046), "Number of columns")
			("YRINTERVAL", po::value<int>()->default_value(15), "Number of years between LU and NGW")
			// These two should be defined on the clinet side not the server side.
			("StartYR", po::value<int>()->default_value(1945), "The starting year of the simulation")
			//  Especially this one
			("NYRS", po::value<int>()->default_value(150), "Number of years to simulate")
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
			//opt.Npixels = vm_cfg["Npixels"].as<int>();
			tf = get_option<int>("Npixels", vm_cfg, opt.Npixels);
			tf = get_option<std::string>("MAPS", vm_cfg, opt.MAPSfile);
			tf = get_option<std::string>("LU", vm_cfg, opt.LUfile);
			tf = get_option<std::string>("NGW", vm_cfg, opt.NGWfile);
			tf = get_option<std::string>("UNSAT", vm_cfg, opt.UNSATfile);
			tf = get_option<std::string>("WELLS", vm_cfg, opt.WELLfile);
			tf = get_option<std::string>("URFS", vm_cfg, opt.URFfile);

			// read optional options
			opt.Nscenarios = vm_cfg["NSCEN"].as<int>();
			opt.port = vm_cfg["PORT"].as<int>();
			opt.Nrow = vm_cfg["Nrow"].as<int>();
			opt.Ncol= vm_cfg["Ncol"].as<int>();
			opt.yearInterval = vm_cfg["YRINTERVAL"].as<int>();
			opt.startYear = vm_cfg["StartYR"].as<int>();
			opt.nSimulationYears = vm_cfg["NYRS"].as<int>();
			opt.nThreads = vm_cfg["NTHREADS"].as<int>();
		}
		return true;
	}

}
