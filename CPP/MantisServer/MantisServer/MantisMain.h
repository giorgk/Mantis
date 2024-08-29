#ifndef MANTISSERVER_MANTISMAIN_H
#define MANTISSERVER_MANTISMAIN_H

#include <sstream>
#include <string>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <chrono>
#include <thread>
#include <boost/bimap.hpp>
#include <cstdlib>
//#include <highfive/H5File.hpp>

//


#include "MSoptions.h"
#include "MShelper.h"
#include "Rch.h"
#include "runtimeWells.h"
#include "Scenario.h"


/**
 * @brief All the classes and structures of Mantis are defined under the mantisServer namespace
 * 
 */
namespace mantisServer {

    //! Set the square pi as constant
    //const double sqrt2pi = std::sqrt(2*std::acos(-1));
    //! Set the pi as constant
    //const double pi = std::atan(1)*4;
    //std::vector<double> dummyVector(1,0);








	/*! \class YearIndex
		\brief A map between the years and the indices in the arrays.
	
		This is a wrap for a bidirectional map between the years in 
		XXXX format and a zero based indexing
	*/
	class YearIndex {
	public:
		/// The constructor does nothing other than making sure the map is empty.
		YearIndex();
		/// Returns the index for a given year. If the year is not in the map it returns a negative -9.
		int get_index(int year);
		/// Returns the year for a given index. If the index is not in the map it returns a negative -9.
		int get_year(int index);
		//! This function constructs or reconstructs the bidirectional map.
		/*!
		\param startYear is the starting year in the map. The format is YYYY. 
		the startYear will be associated with index 0.
		\param nYears The number of years to insert in the map.
		*/
		void reset(int startYear, int nYears);
	private:
		/// This is the method that actually creates the map. It is called from the \ref reset.
		void build_index(int startYear, int nYears);
		/// This is the actual bidirectional container.
		boost::bimap<int, int> yearIndexMap;
	};
	YearIndex::YearIndex() {
		yearIndexMap.clear();
	}
	void YearIndex::build_index(int startYear, int nYears) {
		int index = 0;
		int yr = startYear;
		while (true) {
			yearIndexMap.insert(boost::bimap<int, int>::value_type(yr, index));
			index = index + 1;
			yr = yr + 1;
			if (index > nYears)
				break;
		}
	}
	void YearIndex::reset(int startYear, int nYears) {
		yearIndexMap.clear();
		build_index(startYear, nYears);
	}

	int YearIndex::get_index(int year) {
		boost::bimap<int, int>::left_map::const_iterator it = yearIndexMap.left.find(year);
		if (it == yearIndexMap.left.end())
			return -9;
		else
			return it->get_right();
	}
	int YearIndex::get_year(int index) {
		boost::bimap<int, int>::right_map::const_iterator it = yearIndexMap.right.find(index);
		if (it == yearIndexMap.right.end())
			return -9;
		else
			return it->get_left();
	}









	/**
	 * @brief Stores data for each streamline. 
	 * 
	 * Although this is a class, it is used more like struct container.
	 * 
	 */

	/*class streamlineClass {
	public:
		*//*! \brief streamlineClass constructor expects the parameters that define a streamline
		\param row_ind is the index in GNLM loading where this streamline starts from near the land surface.
		\param col_ind is the index in SWAT loading where this streamline starts from near the land surface.
		\param w_in is the weight of this streamline. This is proportional to the velocity at the well side of the streamline.
		\param rch_in is the groundwater recharge rate in m/day according to the flow model
		\param type_in is the type of the unit response function
		\param paramA this is either the mean value or the streamline length.
		\param paramB this is either the standard deviation or the velocity.
		\param paramC if the type is both this is mean, while A and B are length and velocity
		\param paramD if the type is both this is standard deviation
		*//*
		streamlineClass(int row_ind, int col_ind, double w_in, int npxl_in, URFTYPE type_in, int Riv,
                  double paramA, double paramB, double paramC = 0, double paramD = 0);
		
		//! the row number of the pixel where this streamline starts from near the land surface.
		int row;
		
		//! the column number of the pixel where this streamline starts from near the land surface.
		int col;
		
		//! the mean value of the fitted unit response function.
		double mu;
		
		//! the standard deviation of the fitted unit response function.
		double std;
		
		//!the weight of this streamline. This is proportional to the velocity at the well side of the streamline.
		double w;

		//! The streamline length
		double sl;

		//! the mean velocity along the streamline
		double vel;

		//! This is the groundwater recharge according to the flow model
		//double gwrch;

		//! This is true for the streamlines that start from rivers
		bool inRiver;

		//! The type of streamline
		URFTYPE type;

		int Npxl;

		std::vector<cell> SourceArea;


	};
	streamlineClass::streamlineClass(int row_ind, int col_ind, double w_in, int npxl_in, URFTYPE type_in, int Riv,
                                  double paramA, double paramB, double paramC, double paramD) {
		row = row_ind;
		col = col_ind;
		w = w_in;
		Npxl = npxl_in;
		type = type_in;
		inRiver = Riv == 1;
		switch (type)
		{
		case URFTYPE::LGNRM:
			mu = paramA;
			std = paramB;
			sl = paramC;
			vel = paramD;
			break;
		case URFTYPE::ADE:
			sl = paramA;
			vel = paramB;
			mu = paramC;
			std = paramD;
			break;
		case URFTYPE::BOTH:
			sl = paramA;
			vel = paramB;
			mu = paramC;
			std = paramD;
			break;
		default:
			break;
		}
	}*/

	/*! \class wellClass
		\brief Contains the streamline information for a well.

		Each well is associated with multiple streamlines.
	*//*
	class wellClass {
	public:
		*//**
		 * @brief inserts information about a streamline that reaches the well.
		 * 
		 * @param Sid is the streamline id.
		 * @param gnlm_ind see streamlineClass::streamlineClass for the rest of the parameters
		 * @param swat_ind 
		 * @param w 
		 * @param npxl
		 * @param type
		 * @param riv
		 * @param paramA 
		 * @param paramB 
		 * @param paramC 
		 * @param paramD 
		 *//*
		void addStreamline(int Sid, int row_ind, int col_ind, double w, int npxl, URFTYPE type, int riv,
			double paramA, double paramB, double paramC = 0, double paramD = 0);

		//! streamlines is a map where the key is the streamline id and the value is an object of type streamlineClass::streamlineClass.
		std::map<int, streamlineClass> streamlines;
		double xcoord;
		double ycoord;
		double depth;
		double screenLength;
		double pumpingRate;
		double ratio;
		double angle;
		void setAdditionalData(double x, double y, double d, double s, double q, double r, double a);
	};

	void wellClass::addStreamline(int Sid, int row_ind, int col_ind, double w, int npxl, URFTYPE type, int riv,
		double paramA, double paramB, double paramC, double paramD) {
		streamlines.insert(std::pair<int, streamlineClass>(Sid, streamlineClass( row_ind, col_ind, w, npxl, type, riv, paramA, paramB, paramC, paramD)));
	}

	void wellClass::setAdditionalData(double x, double y, double d, double s, double q, double r, double a){
	    xcoord = x;
	    ycoord = y;
	    depth = d;
	    screenLength = s;
	    pumpingRate = q;
	    ratio = r;
	    angle = a;
	}

	class wellCollection{
	public:
        wellCollection(){}
        wellCollection(std::string name, bool tf)
            :
            Name(name),
            calcSourceArea(tf)
            {}

	    std::string Name;
	    std::map<int, wellClass> Wells;
	    bool calcSourceArea = false;
	    std::string rch_map;
	    void addWell(int Eid, wellClass w){
            Wells.insert(std::pair<int, wellClass>(Eid, w));
	    }
	};*/




	/**
	 * @brief Contains a list of wells that included in the unit area.
	 * 
	 * The unit area is the smallest division of a given background map.
	 * For example if the background map is the Basins map then the Sacramento valley is a unit area. \n 
	 * If the background map contains the CVHM farms then each farm is a Polyregion.
	 *
	 * This class contains the ids of the wells that each Polyregion contains.
	 *//*
	class Polyregion {
	public:
		//! The well ids that are contained by the polygon. 
		//! This is actually a map so that it can contain mutliple well sets where the name of the set is the key
		std::map<std::string, std::vector<int> > wellids;
	};*/

	/*! \class Mantis
		\brief This is the main class of the server.

		The workflow for the Mantis class is the following:

		- \b Construction \n
			To declare a Mantis object a mantisServer::options object is needed.
		- \b initialization \n
			The initialization corresponds to calling the \ref readInputs method.
			This reads the input files and does all the prepocessing steps needed for the simulation.
		- <b>Reading client input</b> \n
			Before doing any simulation Mantis nees to read the client message using the Mantis::parse_incoming_msg method
		- \b Simulate \n
			If the previous step was successfull Mantis can start the simulation.
			However prior to that we should get rid of previous replies by calling Mantis::resetReply. \n
			Then Mantis can safely proceed with the simulation by calling Mantis::simulate_with_threads. 
			Note that there is a non threaded method simulate which may not work at the moment.
		- <b>Prepare output</b> \n
			Last we convert the output to a stream of text using Mantis::makeReply method.
		.

	*/
	class Mantis {
	public:
		Mantis(mantisServer::options options_in);

		//void simulate(std::string &msg, std::string &outmsg);

		void simulate_with_threads(int id);//int id, std::string &msg, std::string &outmsg

		void simulate_RF_wells(int id);

		//bool simulate_streamline(int Eid, int Sid, int unsat_idx, int &I, int &J,
        //                         std::vector<cell> &sourceArea, bool &inRiv,
        //                         double &m, double &s, double &w,
        //                         std::vector<double> &weightedBTC,
        //                         std::ofstream &lf_file, std::ofstream &urf_file,
        //                         std::ofstream &btc_file);
		
		bool readInputs();
		/**
		The incoming message is a string where the defferent pieces of information are separated by a single space.
		The end of the message is denoted by the end line character "\n"
		The format of the message is the following. (Note that the entire message is a single line)
		Nyears: Number of years to simulate e.g 150 200 500 etc. (This must be integer)
		Year to start reduction: 4 digit Integer XXXX e.g 2020, 2025, 2034
		Unsaturated zone mobile water content: This is a double in m^3/m^3. According to Stathi, candidate values are 0.05,0.1 ,0.15 and 0.20 
		Scenario Name: This is a code name that is printed in the files that are read by readWellSet and readURFs
		MapID: The id of the background map
		Nregions: The number of regions to include in the simulation
		Region ids: These are the modelArea ids. The message should containt exactly Nregions ids.
		Ncategories: The number of crops to reduce the loading. 
		Then Repeat Ncategories times the following line
		CropID LoadPercentage (Integer, double)
		CropID is the crop id as listed in the LU file.
		LoadPercentage is number between 0 and 1 which correspond to how much we should reduce the loading
		Next Print the following as it appears
		ENDofMSG
		Finish the message always with an \n 
		*/
		bool parse_incoming_msg(std::string &msg, std::string &outmsg);
		bool validate_msg(std::string& outmsg);
		
		void resetReply();

		/**
		Prepares the following message to send.
		1 or 0. If the computation was successfull the message will start with 1.
		If the message starts with 0 then an error message follows.
		In case of 1:
		nBTC The number of breakthrough curves
		*/
		void makeReply(std::string &outmsg);

		void postReplyActions();

        void resetLogfile();

	private:
		mantisServer::options options;

        std::map<std::string, Region> RegionList;

		/*!
		 * CVrasterClass is a 2D raster map that contains the indices of the linear vector.
		 * All the info that is formatted as 1D or 2D arrays with length 1 - Ncells
		 * must use the CVrasterClass indices to indentify the associated record.
		 * For example LU[0][CVrasterClass[col][row]].
		 * Note however that CVrasterClass contains negative numbers. which correspond to the
		 * cells that have no associated info. These are outside the study area.
		 */
		//std::vector<std::vector<int>> CVraster;
		//CVrasterClass cvraster;

		//std::map<std::string, NLoad> NGWLoading;

		//! UNSAT is a linear data structure that holds the depth/recharge values
        //UNSATdataClass unsat;
        //LinearData unsat;
        //! rch is a linear data structure that holds the groundwater recharge in mm/year
        //LinearData rch;
        //RechargeScenarioList rchList;

		//std::map<std::string, int > UNSATscenarios;
        //std::vector<std::vector<double>> UNSATData;
		//! a list of background maps
		//std::map<std::string, std::map<std::string, Polyregion> > MAPList;
		//td::vector<std::string> backgroundMapNames;

		//! A map of well ids and wellClass
		//std::map<std::string, wellCollection > Wellmap;

		//std::map<std::string, runtimeURFSet> RegionFlowURFS;
        //runtimeURFSet temporaryURFs;

		Scenario scenario;

		std::vector<std::string> replymsg;
		std::vector<int> replyLength;

		/**
		First line:
		Nmaps : The Number of background maps.
		Background maps can be the outline polygon, the subbasins, Counties, Subregions etc.
		Repeat Nmaps times the following
		Nregions : The number of regions this map containts, 1 for entire CV, 3 for subbasins etc.
		For each subregion repeat:
		Npoly: The number of polygons that the subregion consist from.
		For each Npoly repeat:
		Nverts : The number of polygon corners of the subregion polygon
		Repeat Nverts the following
		xcoord ycoord : the coordinates of the polygon order
		*/
		//bool readBackgroundMaps();

		/**
		 * @brief Reads a file with information about a well set
		 * 
		 * The format of the file is the following: \n
		 * NWELLS SETNAME \n
		 * where NWELLS is the number of wells in the file and SETNAME is 
		 * a code name. The code name must be identical to the one used in the URFS
		 * 
		 * Next repeate NWELLS times the following line \n 
		 * EID XW YW \n
		 * 
		 * where EID is the entity id (well id), XW and YW are the coordinates of the wells.
		 * 
		 * @param filename This is the name of the file
		 * @return true if the reading was succesfull
		 * @return false if reading failed for any reason
		 */
		//bool readWellSet(std::string filename);


		//bool readMultipleSets(std::string filename, bool isWell);

		/**
		 * @brief Reads the parameteres of the Unit Respons 
		 * Functions of a given well set
		 * 
		 * The first line containts the following labels: \n
		 * NURFS SETNAME URFTYPE \n
		 * where: \n
		 * NURFS is the number of urfs that follow, \n
		 * SETNAME is the name of the set. This is a string which must be 
		 * identical to the one in the wellfile Mantis::readWellSet(std::string filename) and \n 
		 * URFTYPE is a string from the mantisServer::URFTYPE that specifies the types of parameters to read
		 * 
		 * Next repeate NURFS times the following line: \n
		 * EID, SID, ROW, COL, WEIGHT, {parameters} \n
		 * where \n
		 * EID is the entity id (the well id of the well file) \n
		 * SID is the streamline id. \n
		 * ROW and COL is the pixel row and column where this streamline exits \n
		 * WEIGHT is the weight of this urf. This is usually proportional to velocity at the well side \n
		 * Velocity at the land side. This is used for the travel time calculation 
		 * {parameters} the parameters depend on the URFTYPE. \n
		 * For mantisServer::URFTYPE::LGNRM the paramters are MEAN STD (mean and standard deviation). \n
		 * For mantisServer::URFTYPE::ADE the parameters are SL VEL (streamline length and velocity). \n
		 * While for mantisServer::URFTYPE::BOTH the parameteres are SL VEL MEAN STD
		 * 
		 * @param filename The name of the file
		 * @return true if there were no errors during the reading
		 * @return false if the reading failed
		 */

		//bool readURFs(std::string filename);

		//bool readLU_NGW();

		//bool readUNSAT();

		//bool readRCH();

		//bool readRFURF();

		//bool readCVraster();

        bool readRegionList();
		
		//bool buildLoadingFunction(Scenario &scenario, std::vector<double> &LF, int row, int col, double rch);
        bool buildLoadingFunction(Scenario &scenario, std::vector<double> &LF, std::vector<cell> cells);

		//void identifySurroundingPixel(int Npixels, int row, int col, std::vector<int>& lin_inds);
		//bool CalculateSourceArea();

		void getStartEndWells(int id, int Nwells, unsigned  int &startWell, unsigned int &endWell);

		//void manageRFSets();

        void test_simulation(int threadId, int nThreads);
	};



	Mantis::Mantis(mantisServer::options options_in)
		:
		options(options_in)
	{
        scenario.debugPath = options.DebugPrefix;
    }

	bool Mantis::validate_msg(std::string& outmsg) {
        if (options.testMode){
            return true;
        }

		if ((scenario.endSimulationYear <= scenario.startSimulationYear) || (scenario.endSimulationYear > 2500)){
            scenario.startSimulationYear = 1945;
            scenario.endSimulationYear = 2100;
        }
		if (scenario.startReductionYear < scenario.startSimulationYear || scenario.startReductionYear > scenario.endSimulationYear)
			scenario.startReductionYear = scenario.endSimulationYear; //This desables the reduction
		if (scenario.startReductionYear >= scenario.endReductionYear)
            scenario.endReductionYear = scenario.startReductionYear + 1;
		if (scenario.endReductionYear > scenario.endSimulationYear)
			scenario.endReductionYear = scenario.endSimulationYear;
        bool tf;

        {// Validate Region
            std::map<std::string, Region>::iterator regionit = RegionList.find(scenario.modelArea);
            if (regionit == RegionList.end()){
                outmsg += "0 ERROR: The Region [";
                outmsg += scenario.modelArea;
                outmsg += "] could not be found";
                return false;
            }
            else{
                tf = regionit->second.validateScenario(outmsg, scenario);
                if (!tf){
                    return false;
                }
            }
        }

        if (scenario.unsatZoneMobileWaterContent <= 0) {
            scenario.unsatZoneMobileWaterContent = 0.0;
        }

        if (scenario.userSuppliedConstRed){
            scenario.LoadReductionMap.clear();
            scenario.globalReduction = scenario.constReduction;
        }

        if (scenario.debugID.empty()){
            scenario.printLF = false;
            scenario.printURF = false;
            scenario.printBTC = false;
            scenario.printWellBTC = false;
            scenario.printAdditionalInfo = false;
        }
        else{
            if (scenario.printLF | scenario.printURF | scenario.printBTC | scenario.printWellBTC){
                scenario.printAdditionalInfo = true;
            }
        }

        /*if (scenario.wellType.compare("VM") == 0){
            if (scenario.mapID.compare("Townships") != 0 ){
                outmsg += "0 ERROR: You can simulate VIRTUAL Monitoring wells with Townships only";
                return false;
            }
            scenario.bUseMonitoringWells = true;
            // Load data if needed
            for (int i = 0; i < static_cast<int>(scenario.regionIDs.size()); ++i){
                std::map<std::string, runtimeURFSet>::iterator it;
                std::string scen_name = "TWN_" + scenario.flowScen + "_" + scenario.regionIDs[i];
                std::string scen_path;
                if (!options.bAbsolutePaths){
                    scen_path = options.mainPath + "Townships/" + scenario.flowScen + "/" + scen_name;
                }
                else{
                    scen_path = "Townships/" + scenario.flowScen + "/" + scen_name;
                }
                it = RegionFlowURFS.find(scen_name);
                if (it == RegionFlowURFS.end()){
                    // If the scenario does not exist loadit
                    runtimeURFSet rfset;
                    std::vector<std::string> files;
                    files.push_back(scen_path);
                    bool tf = rfset.readWellSet(scen_name, files);
                    if (!tf){
                        outmsg += "0 ERROR: while reading the files of [";
                        outmsg += scen_name;
                        outmsg += "]";
                        return false;
                    }
                    else {
                        rfset.setlife(options.RFmem);
                        RegionFlowURFS.insert(std::pair<std::string, runtimeURFSet>(scen_name, rfset));
                    }
                }
            }
            scenario.printWellIds = false;
            scenario.printAdditionalInfo = false;
        }
        else{
            std::map<std::string, std::map<std::string, Polyregion> >::iterator mapit = MAPList.find(scenario.mapID);
            if (mapit == MAPList.end()) {
                outmsg += "0 ERROR: The Background map with id [";
                outmsg += scenario.mapID;
                outmsg += "] could not be found";
                return false;
            }

            scenario.flowWellScen = scenario.flowScen + '_' + scenario.wellType;

            int Nwells = 0;
            std::map<std::string, Polyregion>::iterator regit;
            std::map<std::string, std::vector<int> >::iterator wellscenit;
            for (int i = 0; i < static_cast<int>(scenario.regionIDs.size()); ++i) {
                regit = mapit->second.find(scenario.regionIDs[i]);
                if (regit == mapit->second.end()) {
                    outmsg += "0 ERROR: The Region with id [";
                    outmsg += scenario.regionIDs[i];
                    outmsg += "] could not be found in the map with id [";
                    outmsg += scenario.mapID;
                    outmsg += "]";
                    return false;
                }

                wellscenit = regit->second.wellids.find(scenario.flowWellScen);
                if (wellscenit != regit->second.wellids.end()) {
                    Nwells = wellscenit->second.size();
                }
            }
            if (Nwells == 0) {
                outmsg += "0 ERROR: There no wells in the selected regions for the combination of flow scenario [";
                outmsg += scenario.flowScen;
                outmsg += "] and well type [";
                outmsg += scenario.wellType;
                outmsg += "]";
                return false;
            }

            std::map<std::string, wellCollection >::iterator wellscenNameit = Wellmap.find(scenario.flowWellScen);
            if (wellscenNameit == Wellmap.end()) {
                outmsg += "0 ERROR: There are no wells and urfs for the combination of flow scenario [";
                outmsg += scenario.flowScen;
                outmsg += "] and well type [";
                outmsg += scenario.wellType;
                outmsg += "]";
                return false;
            }
        }




		if (unsat.ScenarioIndex(scenario.unsatScenario) == -1) {
			outmsg += "0 ERROR: There is no unsaturated scenario with name: [";
			outmsg += scenario.unsatScenario;
            outmsg += "]";
			return false;
		}
        else{
            scenario.unsatScenarioID = unsat.ScenarioIndex(scenario.unsatScenario);
        }

        if (scenario.buseLoadTransition){
            std::map<std::string, NLoad>::iterator baseloadit = NGWLoading.find(scenario.LoadTransitionName);
            if (baseloadit == NGWLoading.end()){
                outmsg += "0 ERROR: There is no loading scenario with name: [";
                outmsg += scenario.LoadTransitionName;
                outmsg += "]";
                return false;
            }
            if (scenario.LoadTransitionStart < 1945){scenario.LoadTransitionStart = std::min(1945, scenario.endSimulationYear-1);}
            if (scenario.LoadTransitionEnd > scenario.endSimulationYear){scenario.LoadTransitionStart = scenario.endSimulationYear;}
            if (scenario.LoadTransitionEnd <= scenario.LoadTransitionStart){scenario.LoadTransitionEnd = scenario.LoadTransitionStart+1;}
        }

        std::map<std::string, NLoad>::iterator loadit = NGWLoading.find(scenario.loadScen);
        if (loadit == NGWLoading.end()) {
            outmsg += "0 ERROR: There is no loading scenario with name: [";
            outmsg += scenario.loadScen;
            outmsg += "]";
            return false;
        }
        else{
            if (loadit->second.getLtype() == LoadType::RASTER){
                bool testforSubScen = false;
                if (!scenario.modifierName.empty()){
                    scenario.userRasterLoad.setNoDataValue(0.0);
                    bool tf = scenario.userRasterLoad.readData(scenario.modifierName, options.Npixels);
                    if (!tf){
                        outmsg += "0 ERROR: While reading loading [";
                        outmsg += scenario.modifierName;
                        outmsg += "]";
                        return false;
                    }
                    if (scenario.modifierType == RasterOperation::Multiply){
                        testforSubScen = true;
                    }
                    else if (scenario.modifierType == RasterOperation::Replace){
                        testforSubScen = false;
                    }
                    else{
                        outmsg += "0 ERROR: The modifier type is UNKNOWN. Valid options are (Replace or Multiply)";
                    }
                }
                else{
                    testforSubScen = true;
                    scenario.modifierType = RasterOperation::DONTUSE;
                }

                if (testforSubScen){
                    scenario.loadSubScenID = loadit->second.getScenarioID(scenario.loadSubScen);
                    if (scenario.loadSubScenID == -1){
                        outmsg += "0 ERROR: The subScenario: (";
                        outmsg += scenario.loadSubScen;
                        outmsg += ") is not a valid option for the loading scenario ";
                        outmsg += scenario.loadScen;
                        return false;
                    }
                }
            }
        }

		if (scenario.unsatZoneMobileWaterContent <= 0) {
			scenario.unsatZoneMobileWaterContent = 0.0;
		}

		if (scenario.bUseFlowRch){
            scenario.rchName = scenario.flowScen.substr(0,scenario.flowScen.size()-3);
		}
        scenario.rchScenID = rch.ScenarioIndex(scenario.rchName);
        if (scenario.rchScenID == -1){
            outmsg += "0 ERROR: Cannot find the Recharge map : [";
            outmsg += scenario.rchName;
            outmsg += "]";
            return false;
        }

        if (scenario.userSuppliedConstRed){
            scenario.LoadReductionMap.clear();
            scenario.globalReduction = scenario.constReduction;
        }*/

		return true;
	}

	bool Mantis::parse_incoming_msg(std::string &msg, std::string &outmsg) {
        return scenario.parse_incoming_msg(msg, outmsg);
        /*
		bool out = false;
		scenario.clear();
		// convert the string to a stringstream
		std::istringstream ss(msg);
		//ss << msg;
		
		int counter = 0;
		while (true) {
			std::string test, tmp;
			counter++;
			ss >> test;

			if (test == "endSimYear") {
				ss >> scenario.endSimulationYear;
				continue;
			}

			if (test == "startRed" ) {
				ss >> scenario.startReductionYear;
				continue;
			}

			if (test == "endRed") {
				ss >> scenario.endReductionYear;
				continue;
			}

			if (test == "flowScen") {
				ss >> scenario.flowScen;
				continue;
			}

			if (test == "wellType"){
			    ss >> scenario.wellType;
			    continue;
			}

			if (test == "loadScen") {
				ss >> scenario.loadScen;
				continue;
			}

            if (test == "loadSubScen") {
                ss >> scenario.loadSubScen;
                continue;
            }

            if (test == "modifierName") {
                ss >> scenario.modifierName;
                continue;
            }

            if (test == "modifierType") {
                ss >> tmp;
                scenario.modifierType = string2RasterOperation(tmp);
                continue;
            }

            if (test == "modifierUnit") {
                ss >> tmp;
                scenario.modifierUnit = string2LoadUnits(tmp);
                continue;
            }

			if (test == "unsatScen") {
				ss >> scenario.unsatScenario;
				continue;
			}

			if (test == "unsatWC") {
				ss >> scenario.unsatZoneMobileWaterContent;
				continue;
			}


			if (test == "bMap") {
				ss >> scenario.mapID;
				continue;
			}

			if (test == "Nregions") {
                int Nregions;
                ss >> Nregions;
                for (int i = 0; i < Nregions; i++) {
                    ss >> test;
                    scenario.regionIDs.push_back(test);
                }
                continue;
            }

			if (test == "RadSelect"){
                scenario.bNarrowSelection = true;
                scenario.useRadSelect = true;
                double cx, cy, r;
                ss >> cx;
                ss >> cy;
                ss >> r;
                scenario.RadSelect.setData(cx,cy,r);
                continue;
            }

            if (test == "RectSelect"){
                scenario.bNarrowSelection = true;
                scenario.useRectSelect = true;
                double xmin, ymin, xmax, ymax;
                ss >> xmin;
                ss >> ymin;
                ss >> xmax;
                ss >> ymax;
                scenario.RectSelect.setData(xmin, ymin, xmax, ymax);
                continue;
            }

            if (test == "DepthRange"){
                scenario.bNarrowSelection = true;
                scenario.useDepthRange = true;
                double dmin, dmax;
                ss >> dmin;
                ss >> dmax;
                scenario.DepthRange.setData(dmin, dmax);
                continue;
            }

            if (test == "ScreenLenRange"){
                scenario.bNarrowSelection = true;
                scenario.useScreenLenghtRange = true;
                double slmin, slmax;
                ss >> slmin;
                ss >> slmax;
                scenario.ScreenLengthRange.setData(slmin, slmax);
                continue;
            }

            if (test == "SourceArea"){
                int nPix, minPix, maxPix;
                double percPix;
                ss >> nPix;
                ss >> minPix;
                ss >> maxPix;
                ss >> percPix;
                scenario.SourceArea.setParameters(nPix, minPix, maxPix, percPix);
                continue;
            }

			if (test == "Ncrops") {
				int Ncrops, cropid;
				double perc;
				ss >> Ncrops;
				for (int i = 0; i < Ncrops; i++) {
					ss >> cropid;
					ss >> perc;
					if (cropid == -9){
					    scenario.globalReduction = perc;
					}
					else{
                        scenario.LoadReductionMap.insert(std::pair<int, double>(cropid, perc));
					}
				}
				continue;
			}

            if (test ==  "DebugID") {
                ss >> scenario.debugID;
                if (!scenario.debugID.empty())
                    scenario.printAdditionalInfo = true;
                continue;
            }

            if (test == "rchMap"){
                scenario.bUseFlowRch = false;
                ss >> scenario.rchName;
                continue;
            }

            if (test == "minRch") {
                ss >> scenario.minRecharge;
                continue;
            }

            if (test == "maxConc") {
                ss >> scenario.maxConc;
                continue;
            }


            if (test == "RFset"){
                scenario.bUseRuntimeRFsets = true;
                scenario.bUseInitRFsets = false;
                int Nsets;
                ss >> Nsets;
                for (int i = 0; i < Nsets; ++i){
                    ss >> test;
                    scenario.RFSets.push_back(test);
                }
            }


            if (test == "constRed"){
                scenario.userSuppliedConstRed = true;
                ss >> scenario.constReduction;
                continue;
            }


            if (test == "loadTrans"){
                ss >> scenario.LoadTransitionName;
                ss >> scenario.LoadTransitionStart;
                ss >> scenario.LoadTransitionEnd;
                if (scenario.LoadTransitionName.compare("NONE") == 0){
                    scenario.buseLoadTransition = false;
                }
                else{
                    scenario.buseLoadTransition = true;
                }

                continue;
            }

            if (test ==  "PixelRadius") {
                ss >> scenario.PixelRadius;
                continue;
            }

            if (test == "getids"){
                int tmp;
                ss >> tmp;
                scenario.printWellIds = tmp != 0;
                continue;
            }

			if (test == "ENDofMSG") {
				out = true;
				break;
			}
			if (test.empty()) {
				outmsg += "0 ERROR: Empty message was found\n";
				out = false;
				break;
			}

			if (counter > 200) {
				outmsg += "0 ERROR: After 200 iterations I cant find ENDofMSG flag\n";
				out = false;
				break;
			}

            outmsg += "0 ERROR: UNKNOWN option [" + test + "]\n";
            out = false;
            break;
		}
		return out;
        */
	}

	void Mantis::getStartEndWells(int id, int Nwells, unsigned int &startWell, unsigned int &endWell) {
        if (Nwells <= options.nThreads){
            if (id > 0){
                startWell = 0;
                endWell = 0;
            }
            else{
                startWell = 0;
                endWell = Nwells;
            }
        }
        else{
            int nwells2calc = Nwells / options.nThreads;
            startWell = id * nwells2calc;
            endWell = (id + 1)*nwells2calc;
            if (id == options.nThreads - 1)
                endWell = Nwells;
        }
	}

	/*bool Mantis::readBackgroundMaps() {
		//auto start = std::chrono::high_resolution_clock::now();
		std::ifstream MAPSdatafile;
		if (!options.bAbsolutePaths)
            options.MAPSfile = options.mainPath + options.MAPSfile;

        MAPSdatafile.open(options.MAPSfile);
		if (!MAPSdatafile.is_open()) {
			std::cout << "Cant open file: " << options.MAPSfile << std::endl;
			return false;
		}
		else {
            Polyregion polyreg;
			std::cout << "Reading " << options.MAPSfile << std::endl;
			std::string line;
			int Nmaps;
			{ // Get the number of background maps e.g. All area, Basins, counties, farms
				getline(MAPSdatafile, line);
				std::istringstream inp(line.c_str());
				inp >> Nmaps;
			}
			//MAPSdatafile >> Nmaps; 
			for (int imap = 0; imap < Nmaps; ++imap) {
				std::string MapKey;
				int Nregions; // This is the number of regions each background maps consists of 1 for All , 3 for Basins, 21 for farms
				{ // Get the  Background map key and the number of Subregions this background map is divided
					std::getline(MAPSdatafile, line);
					std::istringstream inp(line.c_str());
					inp >> MapKey;
					inp >> Nregions;
				}
				std::map<std::string, Polyregion> RegionMap;
				for (int irg = 0; irg < Nregions; ++irg) {
					std::string RegionKey;
					//int Npoly; //number of polygons for this modelArea

					{// Get the Key for the Subregion of the background map and the number of polygons it consists from
						std::getline(MAPSdatafile, line);
						std::istringstream inp(line.c_str());
						inp >> RegionKey;
					}
                    RegionMap.insert(std::pair<std::string, Polyregion>(RegionKey, polyreg));

				}
				MAPList.insert(std::pair<std::string, std::map<std::string, Polyregion>>(MapKey, RegionMap));
                backgroundMapNames.push_back(MapKey);

			}
		}
		MAPSdatafile.close();
		return true;
	}*/

	bool Mantis::readInputs() {

        if (options.testMode){
            return true;
        }

        bool tf = readRegionList();
        if (!tf) { std::cout << "Error reading Region list file" << std::endl; return false; }


        //bool tf = readCVraster();
        //if (!tf) { std::cout << "Error reading Raster file" << std::endl; return false; }

		//tf = readBackgroundMaps();
		//if (!tf) { std::cout << "Error reading Background Maps" << std::endl; return false; }

        //tf = readUNSAT();
        //if (!tf) { std::cout << "Error reading UNSAT data" << std::endl; return false; }

        //tf = readRCH();
        //if (!tf) { std::cout << "Error reading Recharge data" << std::endl; return false; }

        //tf = readMultipleSets(options.WELLfile, true);
        //if (!tf) { std::cout << "Error reading Wells" << std::endl; return false; }

        //tf = readMultipleSets(options.URFfile, false);
        //if (!tf) { std::cout << "Error reading URFs" << std::endl; return false; }

        //tf = readLU_NGW();
        //if (!tf) { std::cout << "Error reading LU or NGW" << std::endl; return false; }

        //tf = CalculateSourceArea();
        //if (!tf) { std::cout << "Error while calculating Source areas" << std::endl; return false; }

        //if (options.bReadRFURF){
        //    bool tf = readRFURF();
        //    if (!tf) { std::cout << "Error reading Region-Flow specific file" << std::endl; return false; }
        //}

        return tf;
	}

	/*bool Mantis::readRFURF() {
        std::ifstream RFMasterfile;
        std::string filename;
        if (!options.bAbsolutePaths){
            filename = options.mainPath + options.RFURFfile;
        }
        else{
            filename = options.RFURFfile;
        }
        RFMasterfile.open(filename);
        if (!RFMasterfile.is_open()) {
            std::cout << "Cant open file: " << filename << std::endl;
            return false;
        }
        else{
            std::string line, nameset, tmpfile;
            int nfiles;
            while (getline(RFMasterfile, line)){
                std::istringstream inp(line.c_str());
                inp >> nameset;

                if (nameset.empty())
                    continue;
                if (nameset.front() == '#')
                    continue;

                inp >> nfiles;
                std::vector<std::string> files;
                for (int i = 0; i < nfiles; ++i){
                    getline(RFMasterfile, line);
                    std::istringstream inp1(line.c_str());
                    inp1 >> tmpfile;
                    if (!options.bAbsolutePaths)
                        tmpfile = options.mainPath + tmpfile;
                    files.push_back(tmpfile);
                }
                std::map<std::string, runtimeURFSet>::iterator it;
                it = RegionFlowURFS.find(nameset);
                if (it == RegionFlowURFS.end()){
                    runtimeURFSet rfset;
                    bool tf = rfset.readWellSet(nameset, files);
                    if (tf){
                        rfset.setlife(options.RFmem);
                        RegionFlowURFS.insert(std::pair<std::string, runtimeURFSet>(nameset, rfset));
                    }
                    else
                        return false;
                }
            }

            //std::map<std::string, runtimeURFSet>::iterator it1;
            //for (it1 = RegionFlowURFS.begin(); it1 != RegionFlowURFS.end(); ++it1){
            //    std::cout <<it1->first << ":" << it1->second.getNwells() << ", " << it1->second.getNurfs() << std::endl;
            //}
        }

        return true;
	}*/

    bool Mantis::readRegionList(){
        std::ifstream RegionListFStream;
        RegionListFStream.open(options.RegionsListFile);
        if (!RegionListFStream.is_open()){
            std::cout << "Cant open file: " << options.RegionsListFile << std::endl;
            return false;
        }
        else{
            std::string line, filename, RegionName;
            while (getline(RegionListFStream, line)){
                {
                    std::string teststring;
                    if (line.empty())
                        continue;
                    if (line.front() == '#')
                        continue;
                    std::istringstream inptest(line.c_str());
                    inptest >> teststring;
                    if (teststring.empty())
                        continue;
                }


                std::istringstream inp(line.c_str());
                inp >> RegionName;
                inp >> filename;
                RegionList.insert(std::pair<std::string, Region>(RegionName, Region()));
                std::map<std::string, Region>::iterator it = RegionList.find(RegionName);
                if (it != RegionList.end()){
                    std::cout << "||============= " << RegionName << " (" << filename << ") ============||" << std::endl;
                    it->second.readRegionData(filename);
                }
            }
        }
        return true;
    }

	/*bool Mantis::readMultipleSets(std::string filename, bool isWell) {
		std::ifstream WellMasterfile;

		if (!options.bAbsolutePaths)
		    filename = options.mainPath + filename;

		WellMasterfile.open(filename);
		if (!WellMasterfile.is_open()) {
			std::cout << "Cant open file: " << filename << std::endl;
			return false;
		}
		else {
			std::string line, filename1;
			while (getline(WellMasterfile, line)) {
                std::istringstream inp(line.c_str());
                inp >> filename1;
			    if (filename1.empty())
			        continue;
			    if (filename1.front() == '#')
			        continue;

                if (!options.bAbsolutePaths)
                    filename1 = options.mainPath + filename1;
				bool tf;
				if (isWell)
                    tf = readWellSet(filename1);
				else
                    tf = readURFs(filename1);

				if (!tf) {
					std::cout << "An error occurred while reading " << filename1 << std::endl;
					return false;
				}
			}
		}
		return true;
	}*/

	/*bool Mantis::readWellSet(std::string filename) {
		auto start = std::chrono::high_resolution_clock::now();
		std::ifstream Welldatafile;

        std::map<std::string, std::map<std::string, Polyregion> >::iterator mapIt;
        std::map<std::string, Polyregion>::iterator regIt;
        std::map<std::string, std::vector<int>>::iterator scenIt;

		Welldatafile.open(filename);
		if (!Welldatafile.is_open()) {
			std::cout << "Cant open file: " << filename << std::endl;
			return false;
		}
		else {
			std::cout << "Reading " << filename << std::endl;
			int Nwells, Eid;
			std::string setName, line, rch_scen;
			int calcSource;
            getline(Welldatafile, line);
            std::istringstream inp(line.c_str());
			inp >> Nwells;
			inp >> setName;
			inp >> calcSource;
			inp >> rch_scen;
			wellCollection wCollection(setName, calcSource == 1);
            wCollection.rch_map = rch_scen;
            Wellmap.insert(std::pair<std::string, wellCollection>(setName,wCollection));
			double xw, yw, D, SL, Q, ratio, angle;
			std::string regionCode;
			for (int i = 0; i < Nwells; ++i) {
                getline(Welldatafile, line);
                std::istringstream inp1(line.c_str());
				inp1 >> Eid;
				//std::cout << Eid << std::endl;
				inp1 >> xw;
				inp1 >> yw;
				inp1 >> D;
				inp1 >> SL;
				inp1 >> Q;
				inp1 >> ratio;
				inp1 >> angle;
				for (unsigned int j = 0; j < backgroundMapNames.size(); ++j){
				    inp1 >> regionCode;
                    mapIt = MAPList.find(backgroundMapNames[j]);
                    if (mapIt != MAPList.end()){
                        regIt = mapIt->second.find(regionCode);
                        if (regIt != mapIt->second.end()){
                            scenIt = regIt->second.wellids.find(setName);
                            if (scenIt == regIt->second.wellids.end()){
                                std::vector<int> tmpid;
                                tmpid.push_back(Eid);
                                regIt->second.wellids.insert(std::pair<std::string, std::vector<int>>(setName, tmpid));
                            }
                            else{
                                scenIt->second.push_back(Eid);
                            }
                        }
                        else{
                            std::cout << "Cannot find modelArea " <<  regionCode << " in background map " << backgroundMapNames[j] << std::endl;
                        }
                    }
                    else{
                        std::cout << " Cannot find background map with name : " << backgroundMapNames[j] << std::endl;
                    }
				}
				//assign_point_in_sets(xw, yw, Eid, setName);
                wellClass w;
				w.setAdditionalData(xw, yw, D, SL, Q, ratio, angle);
				Wellmap[setName].addWell(Eid, w);
			}
		}
		Welldatafile.close();
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		std::cout << "Read Wells in " << elapsed.count() << std::endl;
		return true;
	}*/

	/*bool Mantis::readURFs(std::string filename) {
        auto start = std::chrono::high_resolution_clock::now();
        std::cout << "Reading " << filename << std::endl;
#if _USEHF>0
        std::string ext = getExtension(filename);
        if (ext.compare("h5") == 0){
            const std::string NamesNameSet("Names");
            const std::string IntsNameSet("ESIJRN");
            const std::string FloatNameSet("MSW");
            HighFive::File HDFfile(filename, HighFive::File::ReadOnly);
            HighFive::DataSet datasetNames = HDFfile.getDataSet(NamesNameSet);
            HighFive::DataSet datasetInts = HDFfile.getDataSet(IntsNameSet);
            HighFive::DataSet datasetFloat = HDFfile.getDataSet(FloatNameSet);
            std::vector<std::string> names;
            std::vector<std::vector<int>> IDS;
            std::vector<std::vector<double>> DATA;
            datasetNames.read(names);
            datasetInts.read(IDS);
            datasetFloat.read(DATA);
            if (names.size() != 2){
                std::cout << "2 names are needed for the URFset. " << names.size() << " provided" << std::endl;
                return false;
            }
            if (IDS[0].size() != DATA[0].size()){
                std::cout << "The rows of integer and float data do not match" << std::endl;
                return false;
            }
            if (IDS.size() != 6 || DATA.size() != 3){
                std::cout << "Incorrect number of columns" << std::endl;
                std::cout << "The size of integers must be 6 and for the floats 3" << std::endl;
                return false;
            }

            URFTYPE urfType;
            if (names[1].compare("LGNRM") == 0)
                urfType = URFTYPE::LGNRM;
            else if (names[1].compare("ADE") == 0)
                urfType = URFTYPE::ADE;
            else if (names[1].compare("BOTH") == 0)
                urfType = URFTYPE::BOTH;
            else{
                std::cout << "The URF type " << names[1] << " is not valid" << std::endl;
                return false;
            }

            std::map<std::string, wellCollection >::iterator scenit;
            scenit = Wellmap.find(names[0]);
            if (scenit == Wellmap.end()) {
                std::cout << "The URF Scenario " << names[0] << " is not defined for the wells" << std::endl;
                return false;
            }
            else{
                std::map<int, wellClass>::iterator wellmapit;
                int eid, sid, r, c, riv, npxl;
                double m, s, w;
                for (unsigned int i = 0; i < IDS[0].size(); ++i){
                    wellmapit = scenit->second.Wells.find(IDS[0][i]);
                    if (wellmapit != scenit->second.Wells.end()){
                        wellmapit->second.addStreamline(IDS[1][i],
                                                        IDS[2][i],
                                                        IDS[3][i],
                                                        DATA[2][i],
                                                        IDS[5][i],
                                                        urfType,
                                                        IDS[4][i],
                                                        DATA[0][i],
                                                        DATA[1][i]);
                    }
                }

            }
            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;
            std::cout << "Read URFS in " << elapsed.count() << std::endl;
            return true;
        }

#endif

        std::ifstream ifile;
        ifile.open(filename);
        if (!ifile.is_open()){
            std::cout << "Cant open file: " << filename << std::endl;
            return false;
        }
        else{
            std::cout << "Reading " << filename << std::endl;
            std::string line;
            int Nurfs;
            URFTYPE urfType;
            std::string name, type;
            {// Read the header
                getline(ifile, line);
                std::istringstream inp(line.c_str());
                inp >> Nurfs;
                inp >> name;
                inp >> type;

                if (type.compare("LGNRM") == 0)
                    urfType = URFTYPE::LGNRM;
                else if (type.compare("ADE") == 0)
                    urfType = URFTYPE::ADE;
                else if (type.compare("BOTH") == 0)
                    urfType = URFTYPE::BOTH;
                else{
                    std::cout << "The URF type " << type << " is not valid" << std::endl;
                    return false;
                }
            }
            std::map<std::string, wellCollection >::iterator scenit;
            scenit = Wellmap.find(name);
            if (scenit == Wellmap.end()) {
                std::cout << "The URF Scenario " << name << " is not defined for the wells" << std::endl;
                return false;
            }
            else{// Read the data
                std::map<int, wellClass>::iterator wellmapit;
                int eid, sid, r, c, riv, npxl;
                double m, s, w;
                for (int i = 0; i < Nurfs; ++i){
                    getline(ifile, line);
                    std::istringstream inp(line.c_str());
                    inp >> eid;
                    inp >> sid;
                    inp >> r;
                    inp >> c;
                    inp >> riv;
                    inp >> m;
                    inp >> s;
                    inp >> w;
                    inp >> npxl;
                    wellmapit = scenit->second.Wells.find(eid);
                    if (wellmapit != scenit->second.Wells.end()){
                        wellmapit->second.addStreamline(sid,r,c,w,npxl,urfType,riv,m,s);
                    }
                }
            }
        }

        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "Read URFS in " << elapsed.count() << std::endl;
		return true;
	}*/

	/*
	void Mantis::identifySurroundingPixel(int Npixels, int row, int col, std::vector<int>& lin_inds) {
	    if (Npixels == 0){
	        if (row >= 0 && row < CVraster[0].size() && col >= 0 && col < CVraster.size()){
	            if (CVraster[col][row] >= 0)
                    lin_inds.push_back(CVraster[col][row]);
	        }
	    }
	    else{
            for (int i = row - Npixels; i <= row + Npixels; ++i){
                if (i < 0 || i >= CVraster[0].size())
                    continue;
                for (int j = col - Npixels; j <= col + Npixels; ++j){
                    if (j < 0 || j >= CVraster.size())
                        continue;
                    if (CVraster[j][i] >= 0)
                        lin_inds.push_back(CVraster[j][i]);
                }
            }
	    }
	}
    */
	bool Mantis::buildLoadingFunction(Scenario &scenario, std::vector<double> &LF, std::vector<cell> cells) {
		bool out = false;
        /*int nCells = scenario.SourceArea.getNpixels(cells.size());
        std::vector<int> lin_idx_vec;

        std::map<std::string, NLoad>::iterator loadit, baseloadit;

		loadit = NGWLoading.find(scenario.loadScen);
        std::vector<double> rch_val;

        if (scenario.buseLoadTransition){
            baseloadit = NGWLoading.find(scenario.LoadTransitionName);
        }

        // Find the linear indices of the source area cells and the recharge value
        for (int i = 0; i < nCells; ++i){
            lin_idx_vec.push_back(cvraster.IJ(cells[i].row, cells[i].col));
            double rtmp = rch.getValue(scenario.rchScenID,lin_idx_vec.back());
            if (rtmp < scenario.minRecharge){
                rtmp = scenario.minRecharge;
            }
            rch_val.push_back(rtmp);
        }
        if (lin_idx_vec.empty()){
            return out;
        }

        std::vector<double> preLF, postLF;
        if (scenario.buseLoadTransition){
            out = baseloadit->second.buildLoadingFunction(lin_idx_vec,
                                                    1945,
                                                    scenario.LoadTransitionEnd,
                                                    preLF,
                                                    scenario,
                                                    rch_val);

            out = loadit->second.buildLoadingFunction(lin_idx_vec,
                                                      scenario.LoadTransitionStart,
                                                      scenario.endSimulationYear,
                                                      postLF,
                                                      scenario,
                                                      rch_val);

            int Nyears = scenario.endSimulationYear - 1945;
            int ia = scenario.LoadTransitionStart - 1945;
            int ib = scenario.LoadTransitionEnd - 1945;
            double da = static_cast<double>(ia);
            double dabRange = static_cast<double>(ib) - da;
            LF.clear();
            LF.resize(Nyears, 0.0);
            for (int iyr = 0; iyr < Nyears; iyr++){
                if (iyr <= ia){
                    LF[iyr] = preLF[iyr];
                }
                else if (iyr >= ib){
                    LF[iyr] = postLF[iyr - ia];
                }
                else{
                    double t = (static_cast<double>(iyr) - ia) / dabRange;
                    LF[iyr] = preLF[iyr]*(1-t) + t*postLF[iyr - ia];
                }

            }

        }
        else{
            out = loadit->second.buildLoadingFunction(lin_idx_vec,
                                                          1945,
                                                          scenario.endSimulationYear,
                                                          LF,
                                                          scenario,
                                                          rch_val);
        }*/

		return out;
	}


	/*bool Mantis::readCVraster() {
        if (!options.bAbsolutePaths)
            options.CVrasterFile = options.mainPath + options.CVrasterFile;
        return cvraster.readData(options.CVrasterFile, options.Nrow, options.Ncol, options.Npixels);
        *//*
	    auto start = std::chrono::high_resolution_clock::now();
        const std::string NameSet("Raster");


        hf::File HDFfile(options.CVrasterFile, hf::File::ReadOnly);
        hf::DataSet dataset = HDFfile.getDataSet(NameSet);
        dataset.read(CVraster);
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "Read CV active raster in " << elapsed.count() << std::endl;
        return true;
        *//*
	}*/


	/*bool Mantis::readUNSAT() {
        if (!options.bAbsolutePaths)
            options.UNSATfile = options.mainPath + options.UNSATfile;
        unsat.setNoDataValue(0.0);
	    return unsat.readData(options.UNSATfile, options.Npixels);
	    *//*
        auto start = std::chrono::high_resolution_clock::now();
        const std::string NamesNameSet("Names");
        const std::string DataNameSet("Data");



        hf::File HDFfile(options.UNSATfile, hf::File::ReadOnly);

        hf::DataSet datasetNames = HDFfile.getDataSet(NamesNameSet);
        hf::DataSet datasetData = HDFfile.getDataSet(DataNameSet);
        std::vector<std::string> names;
        datasetNames.read(names);
        datasetData.read(UNSATData);
        for (unsigned int i = 0; i < names.size(); ++i){
            UNSATscenarios.insert(std::pair<std::string, int>(names[i],i));
        }
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "Read Unsaturated Scenarios in " << elapsed.count() << std::endl;
        return true;
        *//*
	}*/

	/*bool Mantis::readRCH() {
        bool tf1 = rchList.readData(options.mainPath,options.RCHfile, options.Npixels);
        return false;
        if (!options.bAbsolutePaths)
            options.RCHfile = options.mainPath + options.RCHfile;
        rch.setNoDataValue(0.0);
        bool tf = rch.readData(options.RCHfile, options.Npixels);
        //if (tf){
        //    rch.multiply(1.0/365000.0);
        //}
        return tf;
	}*/

	/*bool Mantis::readLU_NGW() {
		auto start = std::chrono::high_resolution_clock::now();

		std::ifstream no3MainFile;

		if (!options.bAbsolutePaths)
		    options.NO3LoadFile = options.mainPath + options.NO3LoadFile;

		no3MainFile.open(options.NO3LoadFile);
		if (!no3MainFile.is_open()) {
			std::cout << "Cant open file: " << options.NO3LoadFile << std::endl;
			return false;
		}
		else {
			std::string line;
			while (getline(no3MainFile, line)) {
				std::string Ltype, Lname, Lfile, Lunit;
				LoadUnits loadunit;
				std::istringstream inp(line.c_str());
                inp >> Ltype;

				if (Ltype.empty())
				    continue;

				if (Ltype.front() == '#')
                    continue;

				inp >> Lname;
				inp >> Lunit;
                loadunit = string2LoadUnits(Lunit);
                if (loadunit == LoadUnits::UNKNOWN){
                    std::cout << Lunit << " is not recognized loading Unit" << std::endl;
                    return false;
                }

				inp >> Lfile;

                if (!options.bAbsolutePaths)
                    Lfile = options.mainPath + Lfile;

				if (Ltype.compare("GNLM") == 0) {
					NLoad NL;
                    NGWLoading.insert(std::pair<std::string, NLoad>(Lname, NL));
                    bool tf = NGWLoading[Lname].readData(Lfile, LoadType::GNLM, loadunit);
                    if (!tf){
                        return false;
                    }

				}
				else if (Ltype.compare("SWAT") == 0) {
					NLoad NL;
                    NGWLoading.insert(std::pair<std::string, NLoad>(Lname, NL));
                    bool tf = NGWLoading[Lname].readData(Lfile, LoadType::SWAT, loadunit);
                    if (!tf){
                        return false;
                    }
				}
				else if (Ltype.compare("RASTER") == 0){
                    NLoad NL;
                    NGWLoading.insert(std::pair<std::string, NLoad>(Lname, NL));
                    bool tf = NGWLoading[Lname].readData(Lfile, LoadType::RASTER, loadunit, options.Npixels);
                    if (!tf){
                        return false;
                    }
				}
				else {
					std::cout << "Unknown loading type " << Ltype << std::endl;
					return false;
				}
			}
		}

		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		std::cout << "Read LU and NGW in " << elapsed.count() << std::endl;
		return true;
	}*/

	void Mantis::simulate_RF_wells(int id) {
        std::ofstream urf_file;
        std::ofstream lf_file;
        std::ofstream btc_file;
        std::ofstream well_btc_file;
        if (scenario.printAdditionalInfo) {
            std::string root_name;
            if (options.DebugPrefix.empty()){
                root_name = scenario.debugID + "_" + num2Padstr(id, 2);
            }
            else{
                root_name = options.DebugPrefix + "_" + scenario.debugID + "_" + num2Padstr(id, 2);
            }
            std::string urf_file_name = root_name + "_urf.dat";
            urf_file.open(urf_file_name.c_str());
            std::string lf_file_name = root_name + "_lf.dat";
            lf_file.open(lf_file_name.c_str());
            std::string btc_file_name = root_name + "_btc.dat";
            btc_file.open(btc_file_name.c_str());
            std::string well_btc_file_name = root_name + "_well_btc.dat";
            well_btc_file.open(well_btc_file_name.c_str());
        }

        /*int unsat_idx = unsat.ScenarioIndex(scenario.unsatScenario);
        int cntBTC = 0;
        std::map<std::string, runtimeURFSet>::iterator it;
        for (int irg = 0; irg < static_cast<int>(scenario.regionIDs.size()); ++irg){
            std::string scen_name = "TWN_" + scenario.flowScen + "_" + scenario.regionIDs[irg];
            it = RegionFlowURFS.find(scen_name);
            if (it != RegionFlowURFS.end()){
                it->second.increaseLife();
                int Nwells = it->second.getNwells();
                unsigned int startWell = 0;
                unsigned int endWell = 0;
                getStartEndWells(id, Nwells, startWell, endWell);
                std::cout << "Thread " << id << " will simulate from [" << startWell << " to " << endWell << ")" << std::endl;

                int NsimulationYears = scenario.endSimulationYear - 1945;
                for (unsigned int iw = startWell; iw < endWell; ++iw){
                    double xw = 0;
                    double yw = 0;
                    bool valid_coords = it->second.getWellCoords(iw,xw,yw);
                    if (scenario.bNarrowSelection == true){
                        if (scenario.useRadSelect && valid_coords){
                            if (!scenario.RadSelect.isPointIn(xw,yw)){
                                continue;
                            }
                        }
                        if (scenario.useRectSelect){
                            if (!scenario.RectSelect.isPointIn(xw,yw)) {
                                continue;
                            }
                        }
                    }



                    std::vector<double> weightedBTC(NsimulationYears, 0);
                    double sumW = 0;
                    int nStreamlines = 0;
                    double m, s, w, d;
                    int urfI, urfJ, inRiv;
                    bool tf;

                    for(unsigned int istrml = 0; istrml < it->second.getNsid(iw); ++istrml){
                        //std::cout << istrml << std::endl;
                        bool valid_param = it->second.getParam(iw, istrml, m, s, w, d, urfI, urfJ, inRiv);
                        if (!valid_param)
                            continue;
                        if (scenario.bNarrowSelection == true){
                            if (scenario.useDepthRange){
                                if (!scenario.DepthRange.isInRange(d)){
                                    continue;
                                }
                            }
                        }

                        std::vector<cell> SourceArea;
                        SourceArea.push_back(cell(urfI, urfJ));
                        bool riv = inRiv == 1;
                        tf = simulate_streamline(iw, istrml, unsat_idx,urfI, urfJ, SourceArea, riv,
                                                 m, s, w ,weightedBTC,
                                                 lf_file, urf_file, btc_file);

                        if (tf){
                            sumW += w;
                            nStreamlines++;
                        }
                    }

                    if (nStreamlines == 0) {
                        continue;
                    }
                    //average streamlines
                    for (int iwbtc = 0; iwbtc < NsimulationYears; ++iwbtc){
                        weightedBTC[iwbtc] = weightedBTC[iwbtc] / sumW;
                        replymsg[id] += std::to_string(static_cast<float>(weightedBTC[iwbtc]));
                        replymsg[id] += " ";
                    }
                    if (scenario.printAdditionalInfo){
                        well_btc_file << iw << " " << nStreamlines << " ";
                        for (int iwbtc = 0; iwbtc < NsimulationYears; ++iwbtc)
                            well_btc_file << std::scientific << std::setprecision(10) << weightedBTC[iwbtc] << " ";
                        well_btc_file << std::endl;
                    }
                    cntBTC++;
                }
            }
        }*/

        /*std::cout << " \tThread " << id << " simulated " << cntBTC << " BTCs" << std::endl;
        replyLength[id] = cntBTC;
        if (scenario.printAdditionalInfo) {
            urf_file.close();
            lf_file.close();
            btc_file.close();
            well_btc_file.close();
        }*/
	}


	void Mantis::simulate_with_threads(int id) {//, , std::string &outmsg
        if (options.testMode){
            test_simulation(id, options.nThreads);
            return;
        }

        std::map<std::string, Region>::iterator regionit;
        regionit = RegionList.find(scenario.modelArea);
        if (regionit != RegionList.end()){
            int nBTC = 0;
            regionit->second.runSimulation(id, options.nThreads, scenario, replymsg, nBTC);
            replyLength[id] = nBTC;
            return;
        }

	    if (scenario.bUseMonitoringWells){
            simulate_RF_wells(id);
            return;
	    }

		std::ofstream urf_file;
		std::ofstream lf_file;
		std::ofstream btc_file;
		std::ofstream well_btc_file;
		if (scenario.printAdditionalInfo) {
            std::string root_name;
            if (options.DebugPrefix.empty()){
                root_name = scenario.debugID + "_" + num2Padstr(id, 2);
            }
            else{
                root_name = options.DebugPrefix + "_" + scenario.debugID + "_" + num2Padstr(id, 2);
            }
			std::string urf_file_name = root_name + "_urf.dat";
			urf_file.open(urf_file_name.c_str());
			std::string lf_file_name = root_name + "_lf.dat";
			lf_file.open(lf_file_name.c_str());
			std::string btc_file_name = root_name + "_btc.dat";
			btc_file.open(btc_file_name.c_str());
			std::string well_btc_file_name = root_name + "_well_btc.dat";
			well_btc_file.open(well_btc_file_name.c_str());
		}


		if (scenario.printAdditionalInfo) {
			urf_file.close();
			lf_file.close();
			btc_file.close();
			well_btc_file.close();
		}

	}

	void Mantis::resetReply() {
		replyLength.clear();
		replyLength.resize(options.nThreads, 0);

		replymsg.clear();
		replymsg.resize(options.nThreads, "");
	}

	void Mantis::makeReply(std::string &outmsg) {
		outmsg.clear();
		
		int nBTC = 0;
		for (int i = 0; i < static_cast<int>(replyLength.size()); ++i) {
			nBTC += replyLength[i];
		}
		if (nBTC == 0) {
			outmsg += "0 ERROR: There are no BTCs to send ENDofMSG\n";
			return;
		}

		outmsg += "1 ";
		outmsg += std::to_string(nBTC);
		int Nyears = scenario.endSimulationYear - scenario.startSimulationYear;
		if (scenario.printWellIds){
		    Nyears = Nyears + 1;
		}
		outmsg += " " + std::to_string(Nyears);

		for (int i = 0; i < static_cast<int>(replymsg.size()); ++i) {
			outmsg += " ";
			outmsg += replymsg[i];
		}
		// This is the character that indicates the end of the message
		outmsg += " ENDofMSG\n";

		//std::cout << outmsg << std::endl;
		std::cout << nBTC << " BTCs will be sent" << std::endl;
	}



	void Mantis::postReplyActions() {
        options.nTimesPrinted = options.nTimesPrinted + 1;
        //manageRFSets();
        resetLogfile();
	}

	/*void Mantis::manageRFSets() {
        std::map<std::string, runtimeURFSet>::iterator it;
        std::vector<std::string> sets4delete;
        for (it = RegionFlowURFS.begin(); it != RegionFlowURFS.end(); ++it){
            it->second.reduceLife();
            if (it->second.getlife() < 0){
                sets4delete.push_back(it->first);
            }
        }
        for (unsigned int i = 0; i < sets4delete.size(); ++i){
            RegionFlowURFS.erase(sets4delete[i]);
        }
	}*/

    void Mantis::resetLogfile() {
        if (options.bUseLogFile){
            if (options.nTimesPrinted > options.nLogReset){
                options.nTimesPrinted = 0;
                logStream.close();
                logStream.open(options.logFile.c_str(), std::ios::out);
                std::cout.rdbuf(logStream.rdbuf());
            }
        }
    }

    void Mantis::test_simulation(int threadId, int nThreads) {
        boost::random::uniform_int_distribution<> dist(50, 200);

        int Nwells = dist(gen);
        unsigned int Nyears = 150;
        int nBTCs = 0;
        if (options.testMode){
            int startId, endId;
            getStartEndIndices(threadId, nThreads, Nwells, startId, endId);
            boost::random::uniform_real_distribution<> distf(0.0, 100.0);
            for (int iw = startId; iw < endId; ++iw){
                for (unsigned int i = 0; i < Nyears; ++i){
                    double v = distf(gen);
                    replymsg[threadId] += std::to_string(static_cast<float>(v));
                    replymsg[threadId] += " ";
                }
                nBTCs = nBTCs + 1;
            }
            replyLength[threadId] = nBTCs;
        }
    }
}

#endif //MANTISSERVER_MANTISMAIN_H