#pragma once

#include <sstream>
#include <fstream>
#include <math.h>
#include <chrono>
#include <thread>
#include <boost/bimap.hpp>
#include <cstdlib>

#define CGAL_HEADER_ONLY 1

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>


#include "MSoptions.h"
//#include "MShelper.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel ine_Kernel;
typedef ine_Kernel::Point_2 ine_Point2;
typedef CGAL::Polygon_2< ine_Kernel> ine_Poly_2;

const double sqrt2pi = std::sqrt(2*std::acos(-1));
const double pi = std::atan(1)*4;
/**
 * @brief All the classes and structures of Mantis are defined under the mantisServer namespace
 * 
 */
namespace mantisServer {

	template<typename T>
	void printVector(std::vector<T>& v, std::string varname) {
		std::cout << std::endl;
		std::cout << varname << " = [";
		for (unsigned int i = 0; i < v.size(); ++i) {
			std::cout << v[i] << " ";
		}
		std::cout << "];" << std::endl;
		std::cout << std::endl;
	}

	/**
	 * @brief Enumeration for the type of the Unit Response function.
	 * 
	 * Since it is very inefficient to hold all 200~300 double values we will keep only 
	 * a number of parameters.
	 * 
	 * If the URF is represented as lognormal distribution the parameteres are the mean and the
	 * standard deviation.
	 * 
	 * We can also store the length and velocity and compute at runtime the analytical 
	 * ADE. However this contains evaluations of the error function and exponential function
	 * multiple times and may be not very efficient
	 * 
	 * Both is used as test where the parameters for both representations are stored and
	 * we can compare the efficiency and results 
	 * 
	 */
	enum class URFTYPE{
		LGNRM, /**< The URF is represented as lognormal distribution with mean and standard deviation parameters*/
		ADE, /**< The URF is represented as analitical ADE with Length and velocity parameters*/
		BOTH /**< The URF has parameteres for both of the above definitions. It is used for testing */
	};

	struct ADEoptions{
		double alpha = 0.32;
		double beta = 0.83;
		double R = 1;
		double lambda = 0.0;
	};


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
			index++;
			yr++;
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

	/*! \class URF
		\brief The URF class contains functionality releated to the Unit Respons Functions.

		At the moment of writing a URF is a function of the following form: \n

		\htmlonly
		<a href="https://www.codecogs.com/eqnedit.php?latex=\huge&space;URF(t)&space;=&space;\frac{1}{t&space;\cdot&space;s&space;\sqrt{2\pi}}\cdot&space;e^{-\frac{(\ln{(t)}-m)^2}{2\cdot&space;s^2}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\huge&space;URF(t)&space;=&space;\frac{1}{t&space;\cdot&space;s&space;\sqrt{2\pi}}\cdot&space;e^{-\frac{(\ln{(t)}-m)^2}{2\cdot&space;s^2}}" title="\huge URF(t) = \frac{1}{t \cdot s \sqrt{2\pi}}\cdot e^{-\frac{(\ln{(t)}-m)^2}{2\cdot s^2}}" /></a>
		\endhtmlonly

		\n
		where \a m is the mean and \a s is the standard deviation. These parameters calculated by fitting the numericaly generated URF from the 
		<a href="https://gwt.ucdavis.edu/research-tools-and-applications/npsat-engine">NPSAT engine</a>. Here we store only the two values.
		This has two advantages. First for each urf we store just 2 double values instead ~200. 
		Secondly we can adjust the simulation time at runtime as we can generate urfs for as longer periods. 
		Nevertheless the fitted parameters correspond to the first 200 years and extending the simulation to let's say 1000 years may not be appropriate.

		The units of the time parameter \a are in years. URF(20) is the value after 20 years of simulation.  
	*/
	class URF {
	public:
		/*! 
		\brief The constructor expects three inputs.

		The constructor expands the Unit respons Function using the input parameters for Nyrs.

		\param Nyrs is the number of years to expand the URF
		\param m_in is the mean value of the fitted equation
		\param s_in is the standard deviation
		*/

		/**
		 * @brief Construct a new URF object
		 * 
		 * The constructor expands the Unit respons Function using the input parameters for Nyrs.
		 * 
		 * @param Nyrs is the number of years to expand the URF
		 * @param paramA This is the first parameter. 
		 * @param paramB  This is the second parameter.
		 * @param type If the type is URFTYPE::LGNRM paramA is the mean and paramB is the standard deviation
		 * If the type is URFTYPE::ADE, paramA is the screen length and paramB is the velocity.
		 */
		URF(int Nyrs, double paramA, double paramB, URFTYPE type, ADEoptions ade_opt = ADEoptions());

		/// calc_conc returns the concentration for a given year. This is actuall used internaly by \ref calc_urf method.
		double calc_conc(double t, ADEoptions ade_opt = ADEoptions());
		/*
		\brief convolute executes the convulotion of the urf with the loading function.
		\param LF is the loading function. The size of the loading function must be equal to \ref urf.
		\param BTC is the output of the convolution. The size of BTC must be equal to \LF. 
		Since this function will be called million times the function doesn't perform and check or allocation for the BTC.

		To improve performance the function does not even loop for loading values less than \ref zeroLoading.
		*/
		void convolute(std::vector<double> &LF, std::vector<double> &BTC);
		
	private:
		/// The mean fitted value.
		double paramA;
		///The standard deviation fitted value
		double paramB;
		/// This is a vector where the expanded values of the urfs will be stored.
		std::vector<double> urf;
		/// calc_urf expands the urf. This takes place during contruction.
		void calc_urf(ADEoptions ade_opt = ADEoptions());
		/// A hard coded threshold. Loading values lower than this threshold are treated as zero 
		/// therefore the convolution can skip an entire loop associated with the particular loading value.
		double zeroLoading = 0.00000000001;
		URFTYPE type;
	};

	URF::URF(int Nyrs, double paramA_in, double paramB_in, URFTYPE type_in, ADEoptions ade_opt)
		:
		paramA(paramA_in),
		paramB(paramB_in),
		type(type_in)
	{
		urf.resize(Nyrs, 0);
		calc_urf(ade_opt);
		//printVector<double>(urf, "URF");
	}

	double URF::calc_conc(double t, ADEoptions ade_opt) {
		double out = 0.0;
		switch (type)
		{
		case URFTYPE::LGNRM: 
		{
			//(1 / (x*b*sqrt(2 * pi)))*exp((-(log(x) - a) ^ 2) / (2 * b ^ 2))
			//std::cout << "t: (m,s)" << t << " (" << m << "," << s << ")" << std::endl;
			double p1 = (1 / (t*paramB*sqrt2pi));
			//std::cout << "p1: " << p1 << std::endl;
			//std::cout << "log(t): " << log(t) << std::endl;
			double p2 = (log(t) - paramA);
			//std::cout << "p2: " << p2 << std::endl;
			p2 = -p2 * p2;
			//std::cout << "p2: " << p2 << std::endl;
			p2 = p2 / (2 * paramB*paramB);
			//std::cout << "p2: " << p2 << std::endl;
			//std::cout << "return: " << p1 * exp(p2) << std::endl;
			out = p1*exp(p2);
		}
			break;
		case URFTYPE::ADE:
		{
			double D = std::pow(ade_opt.alpha*paramA, ade_opt.beta);
			double DRt = 2*std::sqrt(D*ade_opt.R*t);
			if (ade_opt.lambda < 0.000001){
				out = 0.5*std::erfc(( ade_opt.R * paramA - paramB * t )/DRt) + 
						std::sqrt(t*paramB*paramB/pi*D*ade_opt.R) * 
						exp(-1*(ade_opt.R*paramA - paramB*t)*(ade_opt.R*paramA - paramB*t)/(4*D*ade_opt.R*t)) - 
					0.5 * (1+ paramB*paramA/D + t*paramB*paramB/D*ade_opt.R ) * 
						exp(paramB*paramA/D)*
						std::erfc((ade_opt.R*paramA + paramB*t )/DRt);
			}
			else{
				std::cout << "Not implemented for lambda > 0" << std::endl;
			}
		}
			break;
		default:
			out = std::numeric_limits<double>::quiet_NaN();
			break;
		}
		return out;
	}

	void URF::calc_urf(ADEoptions ade_opt) {
		for (int i = 0; i < static_cast<int>(urf.size()); ++i)
			urf[i] = calc_conc(static_cast<double>(i + 1), ade_opt);
	}

	void URF::convolute(std::vector<double> &LF, std::vector<double> &BTC) {
		//std::cout << "LF size: " << LF.size() << std::endl;
		//BTC.resize(LF.size(), 0.0);
		//std::cout << "BTC size: " << BTC.size() << std::endl;
		int shift = 0;
		for (int i = 0; i < static_cast<int>(LF.size()); ++i) {
			if (std::abs(LF[i] - 0) < zeroLoading)
				continue;
			for (int k = shift; k < static_cast<int>(LF.size()); ++k) {
				//std::cout << k - shift << " : " << urf[k - shift] << " " << LF[i] << std::endl;
				BTC[k] = BTC[k] + urf[k - shift] * LF[i];
			}
			shift++;
		}
	}


	/**
	 * @brief is a struct variable that contains all the information needed for each simulation scenario.
	 * 
	 * This struct is a convinent way to transfer around all inputs at once.
	 */
	struct Scenario {
		//! mapID is the id of selected background map.
		int mapID;
		//! scenarioName is the code name of the selected scenario.
		std::string name;
		//! regionIDs is a list of regions to compute the breakthrough curves
		std::vector<int> regionIDs;
		//! LoadReductionMap is a map the sets the nitrate loading reduction for selected land use categories.
		//! The key value of the map is the land use category and the value is the percentage of loading reduction.
		std::map<int, double> LoadReductionMap;
		//! ReductionYear is the year to start the reduction. 
		//! The implementation of these is not fully thought. The best is to choose a starting year that corresponds to
		//! any of the default 15 year time incements.
		int ReductionYear;
		//! This is the number of years to simulate starting from 1945
		int NsimulationYears;
		// This is the unsaturated zone mobile water content (m3/m3)
		// Typical values are 0.05,0.1 ,0.15 and 0.20 
		double unsatZoneMobileWaterContent;
		// Once the Nsimulation year has set this is populated with the actual 4digit years.
		// This is used in the loading function building method 
		//std::vector<int> SimulationYears;
		/**
		 * @brief clear is making sure that the scenario has no data from a previous run.
		 * 
		 */
		void clear() {
			mapID = -9;
			name = "";
			regionIDs.clear();
			LoadReductionMap.clear();
			//SimulationYears.clear();
		}
	};

	/*! \class streamlineClass
		\brief Stores data for each streamline. 

		Although this is a class, it is used more like struct container. 
		I guess when I first writing this I was expecting more functionality here.
	*/

	/**
	 * @brief Stores data for each streamline. 
	 * 
	 * Although this is a class, it is used more like struct container.
	 * 
	 */
	class streamlineClass {
	public:
		/*! \brief streamlineClass constructor expects the parameters that define a streamline
		\param r is the row number of the pixel where this streamline starts from near the land surface.
		\param c is the column number of the pixel where this streamline starts from near the land surface.
		\param w_in is the weight of this streamline. This is proportional to the velocity at the well side of the streamline.
		\param type_in is the type of the unit respons function
		\param paramA this is either the mean value or the streamline length.
		\param paramB this is either the standard deviation or the velocity.
		\param paramC if the type is both this is mean, while A and B are length and velocity
		\param paramD if the type is both this is standard deviation
		*/
		streamlineClass(int r, int c, double w_in, URFTYPE type_in, double paramA, double paramB, double paramC = 0, double paramD = 0);
		
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

		//! The type of streamline
		URFTYPE type;


	};
	streamlineClass::streamlineClass(int r, int c, double w_in, URFTYPE type_in, double paramA, double paramB, double paramC, double paramD) {
		row = r;
		col = c;
		w = w_in;
		type = type_in;
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
	}

	/*! \class wellClass
		\brief Contains the streamline information for a well.

		Each well is associated with multiple streamlines.
	*/
	class wellClass {
	public:
		/**
		 * @brief inserts information about a streamline that reaches the well.
		 * 
		 * @param Sid is the streamline id.
		 * @param r see streamlineClass::streamlineClass for the rest of the parameters
		 * @param c 
		 * @param w 
		 * @param type
		 * @param paramA 
		 * @param paramB 
		 * @param paramC 
		 * @param paramD 
		 */
		void addStreamline(int Sid, int r, int c, double w, URFTYPE type,
			double paramA, double paramB, double paramC = 0, double paramD = 0);

		//! streamlines is a map where the key is the streamline id and the value is an obect of type streamlineClass::streamlineClass.
		std::map<int, streamlineClass> streamlines;
	};

	void wellClass::addStreamline(int Sid, int r, int c, double w, URFTYPE type,
		double paramA, double paramB, double paramC, double paramD) {
		streamlines.insert(std::pair<int, streamlineClass>(Sid, streamlineClass(r, c, w, type, paramA, paramB, paramC, paramD)));
	}

	/*! \class Polyregion
		\brief 

		
	*/

	/**
	 * @brief Contains spatial information of a unit area.
	 * 
	 * A unit area is the smallest division of a given background map. 
	 * For example if the background map is the Basins map then the Sacramento valley is a unit area. \n 
	 * If the background map contains the CVHM farms then each farm is a Polyregion.
	 * 
	 * A Polyregion may consist of one or more polygons. The holes in the geometry are not taken into acount. \n
	 * If there is a hole there should not be any wells in the hole from the first place 
	 * 
	 * In addition to the geometry info this class contains the ids of the wells that each Polyregion contains.
	 * The spatial queries take place during server initialization therefore  during simulation time all that 
	 * is needed is searching through std::map objects.
	 */
	class Polyregion {
	public:
		//! This is a list of cgal polygons that define the unit analysis such as a Basin, a farm, a county etc.
		//! The unit may consists of more polygons that's why we use a vector here
		std::vector<ine_Poly_2> polys;

		//! The well ids that are contained by the polygon. 
		//! This is actually a map so that it can contain mutliple well sets where the name of the set is the key
		std::map<std::string, std::vector<int> > wellids;
	};

	/*! \class Mantis
		\brief This is the main class of the server.

		The workflow for the Mantis class is the following:

		- \b Construction \n
			To declare a Mantis object a mantisServer::options object is needed.
		- \b initialization \n
			The initialization corresopns to calling the \ref readInputs method. 
			This reads the input files and does all the prepocessing steps needed for the simulation.
		- <b>Reading client input</b> \n
			Before doing any simulation Mantis nees to read the client message using the Mantis::parse_incoming_msg method
		- \b Simulate \n
			If the previous step was succesfull Mantis can start the simulation. 
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

		void simulate(std::string &msg, std::string &outmsg);

		void simulate_with_threads(int id);//int id, std::string &msg, std::string &outmsg
		
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
		Region ids: These are the region ids. The message should containt exactly Nregions ids.
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
		
		void resetReply();

		/**
		Prepares the following message to send.
		1 or 0. If the computation was successfull the message will start with 1.
		If the message starts with 0 then an error message follows.
		In case of 1:
		nBTC The number of breakthrough curves
		*/
		void makeReply(std::string &outmsg);


	private:
		mantisServer::options options;

		//! LU is a 3D array of Nrow x Ncol x Years (5)
		std::vector< std::vector< std::vector<int> > > LU;
		//! NGW is a 3D array of Nrow x Ncol x Years (8)
		std::vector< std::vector< std::vector<float> > > NGW;
		//! UNSAT is a 2D array of Nrow x Ncol. If we decide to include more than one travel time map then this should be a 3D array
		std::vector< std::vector<float> > UNSAT;
		//! a list of background maps
		std::map<int, std::map<int, Polyregion> > MAPList;
		//! A map of well ids and wellClass
		std::map<std::string, std::map<int, wellClass> > Wellmap;

		std::vector<double> yearDistributor;
		YearIndex yearIndex;

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
		bool readBackgroundMaps();

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
		bool readWellSet(std::string filename);


		bool readMultipleSets(std::string filename, bool isWell);

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

		bool readURFs(std::string filename);

		bool readLU_NGW();

		void assign_point_in_sets(double x, double y, int wellid, std::string setname);
		
		void buildLoadingFunction(Scenario &scenario, std::vector<double> &LF, int row, int col);
		void distributeYears();
		float unsatTravelTime(int row, int col);
	};



	Mantis::Mantis(mantisServer::options options_in)
		:
		options(options_in)
	{
		//LU.resize(options.Nrow, std::vector< std::vector<int> >(options.Ncol, std::vector<int>(5, 0)));
		//NGW.resize(options.Nrow, std::vector< std::vector<float> >(options.Ncol, std::vector<float>(8, 0)));
		distributeYears();
		yearIndex.reset(options.startYear, options.nSimulationYears);
	}

	bool Mantis::parse_incoming_msg(std::string &msg, std::string &outmsg) {
		scenario.clear();
		// convert the string to a stringstream
		std::stringstream ss;
		ss << msg;

		// ------Start reading the incoming message -------
		// Nyears
		ss >> scenario.NsimulationYears;
		if (scenario.NsimulationYears > 500) {
			outmsg += "0 ERROR: The total number of simulation years must be less than 500\n";
			return false;
		}

		// Populate the Simulation years
		//for (int i = 0; i < scenario.NsimulationYears; i++){
		//	scenario.SimulationYears.push_back(options.startYear + i);
		//}
		yearIndex.reset(options.startYear, scenario.NsimulationYears);

		
		// Year to start reduction
		ss >> scenario.ReductionYear;
		
		// Unsaturated zone mobile water content
		ss >> scenario.unsatZoneMobileWaterContent;
		
		// Get the name of the scenario
		ss >> scenario.name;
		std::map<std::string, std::vector<int> >::iterator wellscenit;

		// Get the selected background map id
		ss >> scenario.mapID;

		std::map<int, std::map<int, Polyregion> >::iterator mapit = MAPList.find(scenario.mapID);
		if (mapit == MAPList.end()) {
			outmsg += "0 ERROR: The Background map with id [";
			outmsg += scenario.mapID;
			outmsg += "] could not be found\n";
			return false;
		}

		// Get the number of selected regions
		int Nregions, RegionID;
		ss >> Nregions;
		std::map<int, Polyregion>::iterator regit;
		// Get the region IDs
		for (int i = 0; i < Nregions; ++i) {
			ss >> RegionID;
			regit = mapit->second.find(RegionID);
			if (regit == mapit->second.end()) {
				outmsg += "0 ERROR: The Region with id [";
				outmsg += RegionID;
				outmsg += "] could not be found in the map with id [";
				outmsg += scenario.mapID;
				outmsg += "]\n";
				return false;
			}
			scenario.regionIDs.push_back(RegionID);

			wellscenit = regit->second.wellids.find(scenario.name);
			if (wellscenit == regit->second.wellids.end()) {
				outmsg += "0 ERROR: There is no scenario with name: ";
				outmsg += scenario.name;
				outmsg += "\n";
				return false;
			}
		}
		// Get the number of crop categories for reduction
		int Ncat;
		ss >> Ncat;
		double perc;
		for (int i = 0; i < Ncat; ++i) {
			ss >> RegionID;
			ss >> perc;
			scenario.LoadReductionMap.insert(std::pair<int, double>(RegionID, perc));
		}

		std::string endofmgs;
		ss >> endofmgs;
		if (endofmgs.compare("ENDofMSG") == 0)
			return true;
		else {
			outmsg += "0 ERROR: The keyword ENDofMSG was not identified. Instead I got";
			outmsg += endofmgs;
			outmsg += "\n";
			return false;
		}
	}

	bool Mantis::readBackgroundMaps() {
		auto start = std::chrono::high_resolution_clock::now();
		std::ifstream MAPSdatafile;
		MAPSdatafile.open(options.MAPSfile);
		if (!MAPSdatafile.is_open()) {
			std::cout << "Cant open file: " << options.MAPSfile << std::endl;
			return false;
		}
		else {
			std::cout << "Reading " << options.MAPSfile << std::endl;
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
					RegionMap[irg+1] = polyRegion;
				}
				MAPList[imap+1] = RegionMap;
			}
		}
		MAPSdatafile.close();
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		std::cout << "Read Maps in " << elapsed.count() << std::endl;
		return true;
	}

	bool Mantis::readInputs() {
		bool tf = readBackgroundMaps();
		if (!tf) { std::cout << "Error reading Background Maps" << std::endl; return false; }
		

		tf = readMultipleSets(options.WELLfile, true);
		if (!tf) { std::cout << "Error reading Wells" << std::endl; return false; }

		tf = readMultipleSets(options.URFfile, false);
		if (!tf) { std::cout << "Error reading URFs" << std::endl; return false; }

		if (!options.testMode) {
			tf = readLU_NGW();
			if (!tf) { std::cout << "Error reading LU or NGW" << std::endl; return false; }
		}
		return tf;
	}

	void Mantis::assign_point_in_sets(double x, double y, int wellid, std::string setname) {
		std::map<int, std::map<int, Polyregion> >::iterator mapit;
		std::map<int, Polyregion>::iterator regit;
		std::vector<ine_Poly_2>::iterator polyit;

		ine_Point2 testPoint(x, y);
		for (mapit = MAPList.begin(); mapit != MAPList.end(); ++mapit) {
			bool found = false;
			for (regit = mapit->second.begin(); regit != mapit->second.end(); ++regit) {
				for (polyit = regit->second.polys.begin(); polyit != regit->second.polys.end(); ++polyit) {
					switch (polyit->bounded_side(testPoint)) {
					case CGAL::ON_BOUNDED_SIDE:
						regit->second.wellids[setname].push_back(wellid);
						found = true;
						break;
					}
				}
				if (found)
					break;
			}
		}
	}

	bool Mantis::readMultipleSets(std::string filename, bool isWell) {
		std::ifstream WellMasterfile;
		WellMasterfile.open(filename);
		if (!WellMasterfile.is_open()) {
			std::cout << "Cant open file: " << filename << std::endl;
			return false;
		}
		else {
			std::string line;
			while (getline(WellMasterfile, line)) {
				bool tf;
				if (isWell)
					tf = readWellSet(line);
				else
					tf = readURFs(line);
				if (!tf) {
					std::cout << "An error occured while reading " << line << std::endl;
					return false;
				}
			}
		}
		return true;
	}

	bool Mantis::readWellSet(std::string filename) {
		auto start = std::chrono::high_resolution_clock::now();
		std::ifstream Welldatafile;

		Welldatafile.open(filename);
		if (!Welldatafile.is_open()) {
			std::cout << "Cant open file: " << filename << std::endl;
			return false;
		}
		else {
			std::cout << "Reading " << filename << std::endl;
			int Nwells, Eid;
			std::string setName;
			Welldatafile >> Nwells;
			Welldatafile >> setName;
			double xw, yw;
			for (int i = 0; i < Nwells; ++i) {
				Welldatafile >> Eid;
				//std::cout << Eid << std::endl;
				Welldatafile >> xw;
				Welldatafile >> yw;
				assign_point_in_sets(xw, yw, Eid, setName);
				Wellmap[setName].insert(std::pair<int, wellClass>(Eid, wellClass()));
			}
		}
		Welldatafile.close();
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		std::cout << "Read Wells in " << elapsed.count() << std::endl;
		return true;
	}

	bool Mantis::readURFs(std::string filename) {
		auto start = std::chrono::high_resolution_clock::now();
		std::ifstream URFdatafile;
		URFdatafile.open(filename);
		if (!URFdatafile.is_open()) {
			std::cout << "Cant open file: " << filename << std::endl;
			return false;
		}
		else {
			std::cout << "Reading " << filename << std::endl;
			std::map<std::string, std::map<int, wellClass> >::iterator scenit;
			std::map<int, wellClass>::iterator wellmapit;
			std::string setName, urftypestring;
			URFTYPE urftype;
			int Nurfs;
			URFdatafile >> Nurfs;
			URFdatafile >> setName;
			scenit = Wellmap.find(setName);
			{// Identify urf type
				URFdatafile >> urftypestring;
				if (urftypestring.compare("LGNRM") == 0)
					urftype = URFTYPE::LGNRM;
				else if (urftypestring.compare("ADE") == 0)
					urftype = URFTYPE::ADE;
				else if (urftypestring.compare("BOTH") == 0)
					urftype = URFTYPE::BOTH;
				else{
					std::cout << "The URF type " << urftypestring << " is not valid" << std::endl;
					return false;
				}
			}

			if (scenit == Wellmap.end()) {
				std::cout << "The URF Scenario " << setName << " is not defined for the wells" << std::endl;
				return false;
			}
			else {
				int Eid, Sid, ROW, COL;
				double paramA, paramB, w, paramC, paramD;
				paramC = 0.0;
				paramD = 0.0;
				for (int i = 0; i < Nurfs; ++i) {
					URFdatafile >> Eid;
					URFdatafile >> Sid;
					URFdatafile >> ROW;
					URFdatafile >> COL;
					URFdatafile >> w;
					URFdatafile >> paramA;
					URFdatafile >> paramB;
					if (urftype == URFTYPE::BOTH){
						URFdatafile >> paramC;
						URFdatafile >> paramD;
					}

					//std::cout << Eid << ", " << Sid << ", " << ROW << ", " << COL << ", " << w << ", " << paramA << ", " << paramB << std::endl;
					//if (Eid > 699) {
					//	bool breakHere = true;
					//}
					//Eid++;// This is used since the particle tracking code writes the entities with zero based numbering
					wellmapit = scenit->second.find(Eid);
					if (wellmapit != scenit->second.end()) {
						wellmapit->second.addStreamline(Sid, ROW, COL, w, urftype, paramA, paramB, paramC, paramD);
					}
					else {
						std::cout << "I cant find well with ID " << Eid << " in the set with name: " << setName << std::endl;
					}
				}
			}
		}
		URFdatafile.close();
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		std::cout << "Read URFS in " << elapsed.count() << std::endl;
		return true;
	}

	void Mantis::simulate(std::string &msg, std::string &outmsg) {
		auto start = std::chrono::high_resolution_clock::now();
		outmsg.clear();
		//Scenario scenario;
		//parse_incoming_msg(msg, scenario);
		
		// Find the selected background map
		std::map<int, std::map<int, Polyregion> >::iterator mapit = MAPList.find(scenario.mapID);
		if (mapit == MAPList.end()) {
			outmsg += "ERROR: The Background map with id [";
			outmsg += scenario.mapID;
			outmsg += "] could not be found";
			return;
		}

		std::map<int, Polyregion>::iterator regit;
		std::map<std::string, std::vector<int> >::iterator wellscenit;
		// precheck that the selected regions exist in the selected map
		for (int i = 0; i < static_cast<int>(scenario.regionIDs.size()); ++i) {
			regit = mapit->second.find(scenario.regionIDs[i]);
			if (regit == mapit->second.end()) {
				outmsg += "ERROR: The Region with id [";
				outmsg += scenario.regionIDs[i];
				outmsg += "] could not be found in the map with id [";
				outmsg += scenario.mapID;
				outmsg += "]";
				return;
			}

			wellscenit = regit->second.wellids.find(scenario.name);
			if (wellscenit == regit->second.wellids.end()) {
				outmsg += "ERROR: There is no scenario with name: ";
				outmsg += scenario.name;
				return;
			}

		}

		std::map<std::string, std::map<int, wellClass> >::iterator wellscenNameit = Wellmap.find(scenario.name);
		if (wellscenNameit == Wellmap.end()) {
			outmsg += "ERROR: There are no wells and urfs for the scenario with name: ";
			outmsg += scenario.name;
			return;
		}
		
		// If we have reached here we know that the map and region ids are valid as well as the selected scenario
		
		// For the time being I'm going to hard code these variables
		URFTYPE urftp = URFTYPE::LGNRM;
		ADEoptions ade_opt;

		std::map<int, streamlineClass>::iterator strmlnit;
		std::map<int, wellClass>::iterator wellit;
		for (int irg = 0; irg < static_cast<int>(scenario.regionIDs.size()); ++irg) {
			regit = mapit->second.find(scenario.regionIDs[irg]);
			wellscenit = regit->second.wellids.find(scenario.name);

			// Number of wells in the selected region
			int Nwells = static_cast<int>(wellscenit->second.size());
			std::cout << Nwells << " wells in Region " << scenario.regionIDs[irg] << std::endl;
			int wellid;
			std::vector<double> LF;
			for (int iw = 0; iw < Nwells; ++iw) {
				wellid = wellscenit->second[iw];
				wellit = wellscenNameit->second.find(wellid);
				if (wellit != wellscenNameit->second.end()) {
					std::vector<double> weightBTC(options.nSimulationYears, 0);
					double sumW = 0;
					if (wellit->second.streamlines.size() > 0) {
						//int temp_count = 1;
						for (strmlnit = wellit->second.streamlines.begin(); strmlnit != wellit->second.streamlines.end(); ++strmlnit) {
							//std::cout << "streamline " << temp_count << std::endl;
							std::vector<double> BTC(options.nSimulationYears, 0);
							if (urftp == URFTYPE::LGNRM){
								URF urf(options.nSimulationYears, strmlnit->second.mu, strmlnit->second.std, urftp);
								buildLoadingFunction(scenario, LF, strmlnit->second.row-1, strmlnit->second.col-1);
								urf.convolute(LF, BTC);
							}
							else if (urftp == URFTYPE::ADE){
								URF urf(options.nSimulationYears, strmlnit->second.sl, strmlnit->second.vel, urftp,ade_opt);
								buildLoadingFunction(scenario, LF, strmlnit->second.row-1, strmlnit->second.col-1);
								urf.convolute(LF, BTC);
							}
							// sum BTC
							for (int ibtc = 0; ibtc < options.nSimulationYears; ++ibtc) {
								weightBTC[ibtc] = weightBTC[ibtc] + BTC[ibtc] * strmlnit->second.w;
							}
							sumW += strmlnit->second.w;
							//temp_count++;
						}
						for (int iwbtc = 0; iwbtc < options.nSimulationYears; ++iwbtc) {
							weightBTC[iwbtc] = weightBTC[iwbtc]/ sumW;
							//std::cout << weightBTC[i] << std::endl;
							outmsg += std::to_string(static_cast<float>(weightBTC[iwbtc]));
							outmsg += " ";
						}
					}
				}
			}
		}
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		std::cout << "Simulation executed in " << elapsed.count() << std::endl;
	}

	void Mantis::buildLoadingFunction(Scenario &scenario, std::vector<double> &LF, int row, int col) {
		// NEW ATTEMPT
		{
			LF.resize(scenario.NsimulationYears, 0);
			// starting and ending values of the year interval
			int YearStart = options.startYear;
			int YearEnd = options.startYear + options.yearInterval;
			int LUsi = 0; //LU starting index
			int LUei = 1; //LU  ending index
			int NGWsi = 0; //NGW starting index
			int NGWei = 1; //NGW ending index
			int LUsc, LUec; // LU starting-ending codes
			double NGWsv, NGWev, Rs, Re; // NGW starting-ending values, Reduction starting and ending value
			Rs = 1; Re = 1;
			// find LU codes and NGW values
			NGWsv = static_cast<double>(NGW[row][col][NGWsi]);
			NGWev = static_cast<double>(NGW[row][col][NGWei]);
			LUsc = LU[row][col][LUsi];
			LUec = LU[row][col][LUei];

			int count_yearly_intervals = 0;
			std::map<int, double>::iterator redit;
			int currentYear;

			for (int i = 0; i < scenario.NsimulationYears; ++i) {
				currentYear = yearIndex.get_year(i);
				if (count_yearly_intervals >= options.yearInterval) {
					YearStart += options.yearInterval;
					YearEnd += options.yearInterval;
					NGWsi++;
					NGWei++;
					if (NGWsi >= static_cast<int>(NGW[0][0].size())) {
						NGWsi = static_cast<int>(NGW[0][0].size()) - 1;
					}
					if (NGWei >= static_cast<int>(NGW[0][0].size())) {
						NGWei = static_cast<int>(NGW[0][0].size()) - 1;
					}
					LUsi++;
					LUei++;
					if (LUsi >= static_cast<int>(LU[0][0].size())) {
						LUsi = static_cast<int>(LU[0][0].size()) - 1;
					}
					if (LUei >= static_cast<int>(LU[0][0].size())) {
						LUei = static_cast<int>(LU[0][0].size()) - 1;
					}

					NGWsv = static_cast<double>(NGW[row][col][NGWsi]);
					NGWev = static_cast<double>(NGW[row][col][NGWei]);
					LUsc = LU[row][col][LUsi];
					LUec = LU[row][col][LUei];
					count_yearly_intervals = 0;
				}

				if (currentYear >= scenario.ReductionYear) {
					// look into the land use codes to see if any reduction has been set
					
					redit = scenario.LoadReductionMap.find(LUsc);
					if (redit == scenario.LoadReductionMap.end())
						Rs = 1;
					else
						Rs = redit->second;
					

					redit = scenario.LoadReductionMap.find(LUec);
					if (redit == scenario.LoadReductionMap.end())
						Re = 1;
					else
						Re = redit->second;
					//if (Rs < 1 || Re < 1)
					//	std::cout << Rs << " - " << Re << std::endl;
				}

				LF[i] = (Rs * NGWsv) * (1 - yearDistributor[count_yearly_intervals]) + (Re * NGWev) * yearDistributor[count_yearly_intervals];

				count_yearly_intervals++;
			}
		}

		/*
		// OLD WAY
		{
			LF.resize(scenario.NsimulationYears, 0);
			// get the id of the year when the reduction starts
			int idRed = yearIndex.get_index(scenario.ReductionYear);
			if (idRed < 0)
				idRed = options.nSimulationYears + 100;

			// starting ending values of the year interval
			int LUsi = 0; //LU starting index
			int LUei = 1; //LU  ending index
			int NGWsi = 0; //NGW starting index
			int NGWei = 1; //NGW ending index
			int LUsc, LUec; // LU starting-ending codes
			double NGWsv, NGWev, Rs, Re; // NGW starting-ending values
			Rs = 1; Re = 1;
			// find LU codes and NGW values
			NGWsv = static_cast<double>(NGW[row][col][NGWsi]);
			NGWev = static_cast<double>(NGW[row][col][NGWei]);
			LUsc = LU[row][col][LUsi];
			LUec = LU[row][col][LUei];

			int count_yearly_intervals = 0;
			std::map<int, double>::iterator redit;

			for (int i = 0; i < static_cast<int>(LF.size()); ++i) {
				if (count_yearly_intervals >= options.yearInterval) {
					NGWsi++;
					NGWei++;
					if (NGWsi >= static_cast<int>(NGW[0][0].size())) {
						NGWsi = static_cast<int>(NGW[0][0].size()) - 1;
					}
					if (NGWei >= static_cast<int>(NGW[0][0].size())) {
						NGWei = static_cast<int>(NGW[0][0].size()) - 1;
					}
					LUsi++;
					LUei++;
					if (LUsi >= static_cast<int>(LU[0][0].size())) {
						LUsi = static_cast<int>(LU[0][0].size()) - 1;
					}
					if (LUei >= static_cast<int>(LU[0][0].size())) {
						LUei = static_cast<int>(LU[0][0].size()) - 1;
					}

					NGWsv = static_cast<double>(NGW[row][col][NGWsi]);
					NGWev = static_cast<double>(NGW[row][col][NGWei]);
					LUsc = LU[row][col][LUsi];
					LUec = LU[row][col][LUei];
					count_yearly_intervals = 0;
				}

				if (i >= idRed) {
					// Then we should look into the land use codes to see if any reduction has been set
					redit = scenario.LoadReductionMap.find(LUsc);
					if (redit == scenario.LoadReductionMap.end())
						Rs = 1;
					else
						Rs = redit->second;

					redit = scenario.LoadReductionMap.find(LUec);
					if (redit == scenario.LoadReductionMap.end())
						Re = 1;
					else
						Re = redit->second;
				}

				LF[i] = (Rs*NGWsv) * (1 - yearDistributor[count_yearly_intervals]) + (Re*NGWev) * yearDistributor[count_yearly_intervals];
				count_yearly_intervals++;
			}
		}// OLD WAY
		*/
	}

	float Mantis::unsatTravelTime(int row, int col) {
		return UNSAT[row][col];
	}

	void Mantis::distributeYears() {
		yearDistributor.clear();
		double x = 0;
		double Nx = static_cast<double>(options.yearInterval);
		for (int i = 0; i < options.yearInterval; ++i) {
			yearDistributor.push_back(x / Nx);
			x += 1.0;
		}
	}

	bool Mantis::readLU_NGW() {
		auto start = std::chrono::high_resolution_clock::now();
		// allocate memory
		//12863 x 7046 x 5 
		LU.resize(options.Nrow, std::vector< std::vector<int> >(options.Ncol, std::vector<int>(5, 0)));
		NGW.resize(options.Nrow, std::vector< std::vector<float> >(options.Ncol, std::vector<float>(8, 0)));
		UNSAT.resize(options.Nrow, std::vector<float>(options.Ncol, 0));
		std::ifstream LUdatafile, NGWdatafile, UNSATdatafile;
		
		LUdatafile.open(options.LUfile);
		if (!LUdatafile.is_open()) {
			std::cout << "Cant open file: " << options.LUfile << std::endl;
			return false;
		}

		NGWdatafile.open(options.NGWfile);
		if (!NGWdatafile.is_open()) {
			std::cout << "Cant open file: " << options.NGWfile << std::endl;
			return false;
		}

		UNSATdatafile.open(options.UNSATfile);
		if (!UNSATdatafile.is_open()) {
			std::cout << "Cant open file: " << options.UNSATfile << std::endl;
			return false;
		}
		std::cout << "Reading " << options.LUfile << std::endl;
		std::cout << "Reading " << options.NGWfile << std::endl;
		std::cout << "Reading " << options.UNSATfile << std::endl;

		int data, I, J;
		float d;

		for (int i = 0; i < options.Npixels; ++i) {
			LUdatafile >> data;
			I = data - 1;
			LUdatafile >> data;
			J = data - 1;
			for (unsigned int k = 0; k < 5; ++k) {
				LUdatafile >> data;
				LU[I][J][k] = data;
			}

			for (unsigned int k = 0; k < 8; ++k) {
				NGWdatafile >> d;
				NGW[I][J][k] = d;
			}
			UNSATdatafile >> d;
			UNSAT[I][J] = d;
		}
		LUdatafile.close();
		NGWdatafile.close();
		UNSATdatafile.close();
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		std::cout << "Read LU and NGW in " << elapsed.count() << std::endl;
		return true;
	}

	void Mantis::simulate_with_threads(int id) {//, , std::string &outmsg

		// Get an iterator to the selected map
		std::map<int, std::map<int, Polyregion> >::iterator mapit = MAPList.find(scenario.mapID);
		// Get an iterator to the list of wells for the selected scenario
		std::map<std::string, std::map<int, wellClass> >::iterator wellscenNameit = Wellmap.find(scenario.name);

		std::map<int, Polyregion>::iterator regit;
		std::map<std::string, std::vector<int> >::iterator wellscenit;

		

		//for (int i = startWell; i < endWell; ++i) {
		//	std::cout << id << ":" << i << std::endl;
		//}

		std::map<int, streamlineClass>::iterator strmlnit;
		std::map<int, wellClass>::iterator wellit;
		int cntBTC = 0;
		for (int irg = 0; irg < static_cast<int>(scenario.regionIDs.size()); ++irg) {
			regit = mapit->second.find(scenario.regionIDs[irg]);
			wellscenit = regit->second.wellids.find(scenario.name);

			// Number of wells in the selected region
			int Nwells = static_cast<int>(wellscenit->second.size());
			int startWell, endWell;

			if (Nwells < options.nThreads) {
				if (id > 0)
					return;
				else {
					startWell = 0;
					endWell = Nwells;
				}
			}
			else {
				int nwells2calc = Nwells / options.nThreads;
				startWell = id * nwells2calc;
				endWell = (id + 1)*nwells2calc;
				if (id == options.nThreads - 1)
					endWell = Nwells;
			}
			
			int NsimulationYears = scenario.NsimulationYears;
			
			

			std::cout << id << " will simulate from [" << startWell << " to " << endWell << ")" << std::endl;
			int wellid;
			std::vector<double> LF;
			
			for (int iw = startWell; iw < endWell; ++iw) {
				wellid = wellscenit->second[iw];
				//std::cout << wellid << std::endl;
				wellit = wellscenNameit->second.find(wellid);
				if (wellit != wellscenNameit->second.end()) {
					std::vector<double> weightBTC(NsimulationYears, 0);
					double sumW = 0;
					if (wellit->second.streamlines.size() > 0) {
						if (!options.testMode) {
							int cnt_strmlines = 0;
							for (strmlnit = wellit->second.streamlines.begin(); strmlnit != wellit->second.streamlines.end(); ++strmlnit) {
								//std::cout << cnt_strmlines++ << std::endl;
								// do convolution only if the source of water is not river. When mu and std are 0 then the source area is river
									std::vector<double> BTC(NsimulationYears, 0);
								if (std::abs(strmlnit->second.mu - 0) > 0.00000001) {
									if (strmlnit->second.type == URFTYPE::LGNRM){
										URF urf(NsimulationYears, strmlnit->second.mu, strmlnit->second.std, strmlnit->second.type);
										buildLoadingFunction(scenario, LF, strmlnit->second.row - 1, strmlnit->second.col - 1);
										//printVector<double>(LF, "LF");
										urf.convolute(LF, BTC);
										//printVector<double>(BTC, "BTC");
									}
									else if (strmlnit->second.type == URFTYPE::ADE){
										URF urf(NsimulationYears, strmlnit->second.mu, strmlnit->second.std, strmlnit->second.type, ADEoptions());
										buildLoadingFunction(scenario, LF, strmlnit->second.row - 1, strmlnit->second.col - 1);
										urf.convolute(LF, BTC);
									}

									// Find the travel time in the unsaturated zone
									float tau = unsatTravelTime(strmlnit->second.row-1, strmlnit->second.col-1);
									//std::cout << tau << std::endl;
									tau = std::floor(tau * scenario.unsatZoneMobileWaterContent);
									if (tau < 0)
										tau = 0.0;
									int intTau = static_cast<int>(tau);

									if (intTau < NsimulationYears) {
										int ibtc = 0;
										for (int ii = intTau; ii < NsimulationYears; ++ii) {
											//std::cout << ii << std::endl;
											weightBTC[ii] = weightBTC[ii] + BTC[ibtc] * strmlnit->second.w;
											ibtc++;
										}
									}

									//for (int ibtc = 0; ibtc < NsimulationYears; ++ibtc) {
									//	weightBTC[ibtc] = weightBTC[ibtc] + BTC[ibtc] * strmlnit->second.w;
									//}
								}
								sumW += strmlnit->second.w;
							}
							//average streamlines
							for (int iwbtc = 0; iwbtc < NsimulationYears; ++iwbtc) {
								weightBTC[iwbtc] = weightBTC[iwbtc] / sumW;
								//std::cout << weightBTC[i] << std::endl;
								replymsg[id] += std::to_string(static_cast<float>(weightBTC[iwbtc]));
								replymsg[id] += " ";
							}
							//printVector<double>(weightBTC, "wBTC");
						}
						else {
							for (int iwbtc = 0; iwbtc < options.nSimulationYears; ++iwbtc) {
								replymsg[id] += std::to_string(static_cast<float>(std::rand() % 100) / 10);
								replymsg[id] += " ";
							}
						}
						cntBTC++;
					}
				}
			}
		}// loop through regions
		replyLength[id] = cntBTC;
	}

	void Mantis::resetReply() {
		replyLength.clear();
		replyLength.resize(options.nThreads, 0);

		replymsg.clear();
		replymsg.resize(options.nThreads, "");
	}

	void Mantis::makeReply(std::string &outmsg) {
		outmsg.clear();
		outmsg += "1 ";
		int nBTC = 0;
		for (int i = 0; i < static_cast<int>(replyLength.size()); ++i) {
			nBTC += replyLength[i];
		}

		outmsg += std::to_string(nBTC);
		for (int i = 0; i < static_cast<int>(replymsg.size()); ++i) {
			outmsg += " ";
			outmsg += replymsg[i];
		}
		// This is the character that indicates the end of the message
		outmsg += " ENDofMSG\n";

		std::cout << nBTC << "BTCs will be sent" << std::endl;
	}
}
