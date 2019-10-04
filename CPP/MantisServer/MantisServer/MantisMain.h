#pragma once

#include <sstream>
#include <fstream>

#define CGAL_HEADER_ONLY 1

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>


#include "MSoptions.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel ine_Kernel;
typedef ine_Kernel::Point_2 ine_Point2;
typedef CGAL::Polygon_2< ine_Kernel> ine_Poly_2;

namespace mantisServer {

	class Polyregion {
	public:
		// This is a list of boost polygons that define a unit analysis such as a Basin, a farm, a county etc.
		// The unit may consists of more than one polygons
		std::vector<ine_Poly_2> polys;
		// The well ids that are contained by the polygon. This is actually a map so that it can containt mutliple well sets where the name of the set is the key
		std::map<std::string, std::vector<int> > wellids;
	};


	class Mantis {
	public:
		Mantis(mantisServer::options options_in);

		bool parse_input(std::string& msg);
		bool readInputs();

	private:
		mantisServer::options options;

		//! LU is a 3D array of Nrow x Ncol x Years (5)
		std::vector< std::vector< std::vector<int> > > LU;
		//! NGW is a 3D array or Nrow x Ncol x Years (8)
		std::vector< std::vector< std::vector<float> > > NGW;
		//! a list of background maps
		std::map<int, std::map<int, Polyregion> > MAPList;

		//! mapID is the selected background map
		int mapID;
		//! scenarioName is the name of the selected scenario
		std::string scenarioName;
		//! is a list of region to compute the breakthrough curves
		std::vector<int> regionIDs;
		//! LoadReduction is a mpa the sets the nitrate loading reduction for selected land use categories
		std::map<int, double> LoadReduction;

		bool readBackgroundMaps();
	};



	Mantis::Mantis(mantisServer::options options_in)
		:
		options(options_in)
	{
		LU.resize(options.Nrow, std::vector< std::vector<int> >(options.Ncol, std::vector<int>(5, 0)));
		NGW.resize(options.Nrow, std::vector< std::vector<float> >(options.Ncol, std::vector<float>(8, 0)));
	}

	bool Mantis::parse_input(std::string& msg) {
		// convert the string to a stringstream
		std::stringstream ss;
		ss << msg;

		// Get the name of the scenario
		ss >> scenarioName;
		// Get the selected background map id
		ss >> mapID;
		// Get the number of selected regions
		int n, id;
		ss >> n;
		// Get the region IDs
		regionIDs.clear();
		for (int i = 0; i < n; ++i) {
			ss >> id;
			regionIDs.push_back(id);
		}
		// Get the number of crop categories for reduction
		int Ncat;
		ss >> Ncat;
		double perc;
		for (int i = 0; i < Ncat; ++i) {
			ss >> id;
			ss >> perc;
			LoadReduction.insert(std::pair<int, double>(id, perc));
		}
		return true;
	}

	bool Mantis::readBackgroundMaps() {
		std::ifstream MAPSdatafile;
		MAPSdatafile.open(options.MAPSfile);
		if (!MAPSdatafile.is_open()) {
			std::cout << "Cant open file: " << options.MAPSfile << std::endl;
			return false;
		}
		else {
			int Nmaps;
			MAPSdatafile >> Nmaps; // This is the number of background maps e.g. All area, Basins, counties, farms
			for (int imap = 0; imap < Nmaps; ++imap) {
				int Nregions; // This is the number of regions each background maps consists of 1 for All , 3 for Basins, 21 for farms
				MAPSdatafile >> Nregions;
				std::map<int, Polyregion> RegionMap;
				for (int irg = 0; irg < Nregions; ++irg) {
					int Npoly; // This is the number of polygons for each region. For example TLB region may consists of more than on polygon
					Polyregion polyRegion;
					MAPSdatafile >> Npoly;
					for (int ipl = 0; ipl < Npoly; ++ipl) {
						int Nvert; // Number of vertices per polygon
						MAPSdatafile >> Nvert;
						double xm, ym;
						std::vector<ine_Point2> poly_pnts;
						for (int ivrt = 0; ivrt < Nvert; ++ivrt) {
							MAPSdatafile >> xm;
							MAPSdatafile >> ym;
							poly_pnts.push_back(ine_Point2(xm, ym));
						}
						ine_Poly_2 p(poly_pnts.begin(), poly_pnts.end());
						polyRegion.polys.push_back(p);
					}
					RegionMap[irg] = polyRegion;
				}
				MAPList[imap] = RegionMap;
			}
		}
	}

	bool Mantis::readInputs() {
		bool tf = readBackgroundMaps();
		return tf;
	}

}
