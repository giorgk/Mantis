#pragma once

#include <sstream>
#include <fstream>
#include <math.h>
#include <chrono>
#include <thread>
#include <boost/bimap.hpp>
#include <cstdlib>
#include <highfive/H5File.hpp>

namespace hf = HighFive;

//#define CGAL_HEADER_ONLY 1

//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Polygon_2.h>
//#include <CGAL/Polygon_2_algorithms.h>


#include "MSoptions.h"
//#include "MShelper.h"

//typedef CGAL::Exact_predicates_inexact_constructions_kernel ine_Kernel;
//typedef ine_Kernel::Point_2 ine_Point2;
//typedef CGAL::Polygon_2< ine_Kernel> ine_Poly_2;

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

	enum class LoadType {
		GNLM,
		SWAT
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
		std::string mapID;
		//! scenarioName is the code name of the selected scenario.
		std::string flowScen;
		//! The loading scenario
		std::string loadScen;
		//! regionIDs is a list of regions to compute the breakthrough curves
		std::vector<std::string> regionIDs;
		//! LoadReductionMap is a map the sets the nitrate loading reduction for selected land use categories.
		//! The key value of the map is the land use category and the value is the percentage of loading reduction.
		std::map<int, double> LoadReductionMap;

		//! If the input message contains crop id with -9 then this is applied to all crops except to the ones defined
		//! in the LoadReductionMap
		double globalMod;
		//! ReductionYear is the year to start the reduction. 
		//! The implementation of these is not fully thought. The best is to choose a starting year that corresponds to
		//! any of the default 15 year time increments.
		int startReductionYear;
		//! The end reduction is the year when the loading has reduced to the desired amount
		int endReductionYear;
		//! This is the number of years to simulate starting from 1945
		int endSimulationYear;

		std::string unsatScenario;
		// This is the unsaturated zone mobile water content (m3/m3)
		// Typical values are 0.05,0.1 ,0.15 and 0.20 
		double unsatZoneMobileWaterContent;

		bool printAdditionalInfo = false;
		int debugID;
		// Once the Nsimulation year has set this is populated with the actual 4digit years.
		// This is used in the loading function building method 
		//std::vector<int> SimulationYears;
		/**
		 * @brief clear is making sure that the scenario has no data from a previous run.
		 * 
		 */
		void clear() {
			mapID = "";
			flowScen = "";
			loadScen = "";
			regionIDs.clear();
			LoadReductionMap.clear();
			printAdditionalInfo = false;
			globalMod = 1.0;
            unsatZoneMobileWaterContent = 0.0;
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
		\param row_ind is the index in GNLM loading where this streamline starts from near the land surface.
		\param col_ind is the index in SWAT loading where this streamline starts from near the land surface.
		\param w_in is the weight of this streamline. This is proportional to the velocity at the well side of the streamline.
		\param rch_in is the groundwater recharge rate in m/day according to the flow model
		\param type_in is the type of the unit response function
		\param paramA this is either the mean value or the streamline length.
		\param paramB this is either the standard deviation or the velocity.
		\param paramC if the type is both this is mean, while A and B are length and velocity
		\param paramD if the type is both this is standard deviation
		*/
		streamlineClass(int row_ind, int col_ind, double w_in, double rch_in, URFTYPE type_in, int Riv,
                  double paramA, double paramB, double paramC = 0, double paramD = 0);
		
		//! The index of GNLM loading
		//int gnlm_index;
		
		//! The index of SWAT loading
		//int swat_index;
		
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
		double gwrch;

		bool inRiver;

		//! The type of streamline
		URFTYPE type;


	};
	streamlineClass::streamlineClass(int row_ind, int col_ind, double w_in, double rch_in, URFTYPE type_in, int Riv,
                                  double paramA, double paramB, double paramC, double paramD) {
		row = row_ind;
		col = col_ind;
		w = w_in;
		gwrch = rch_in;
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
		 * @param gnlm_ind see streamlineClass::streamlineClass for the rest of the parameters
		 * @param swat_ind 
		 * @param w 
		 * @param rch
		 * @param type
		 * @param riv
		 * @param paramA 
		 * @param paramB 
		 * @param paramC 
		 * @param paramD 
		 */
		void addStreamline(int Sid, int row_ind, int col_ind, double w, double rch, URFTYPE type, int riv,
			double paramA, double paramB, double paramC = 0, double paramD = 0);

		//! streamlines is a map where the key is the streamline id and the value is an object of type streamlineClass::streamlineClass.
		std::map<int, streamlineClass> streamlines;
	};

	void wellClass::addStreamline(int Sid, int row_ind, int col_ind, double w, double rch, URFTYPE type, int riv,
		double paramA, double paramB, double paramC, double paramD) {
		streamlines.insert(std::pair<int, streamlineClass>(Sid, streamlineClass( row_ind, col_ind, w, rch, type, riv, paramA, paramB, paramC, paramD)));
	}

	class NLoad {
	public:
		NLoad() {}
		bool readData(std::string filename, LoadType ltype);
		float getNload(int index,int iyr);
		void getNload(int index, int iyr, double &N1, double &N2, double &u);
		void getLU(int index, int iyr, int &LUcS, int &LUcE, double &perc);
		int getLU(int index, int iyr);
		void buildLoadingFunction(int index, int endYear, std::vector<double>& LF, Scenario& scenario, double mult);
		LoadType getLtype() {
			return loadType;
		}
		bool isValidIndex(int index);
	private:
		LoadType loadType;
		std::vector<std::vector<double> > Ndata;
		std::vector<std::vector<int> > LU;
		std::vector<int> Nidx;
	};

	bool NLoad::isValidIndex(int index) {
		if ((index >= 0) | (index < Ndata[0].size()))
			return true;
		else
			return false;
	}

	int NLoad::getLU(int index, int iyr) {
	    if ((iyr >=0 & iyr < LU.size()) & (index >=0 & index < LU[0].size()))
	        return LU[iyr][index];
		/*switch (loadType)
		{
		case mantisServer::LoadType::GNLM:
		{
			std::cout << "You can't call this method for GNLM loading type" << std::endl;
			return 0;
			break;
		}
		case mantisServer::LoadType::SWAT:
		{
			// Swat loading does not change over time
			return LU[index][0];
			break;
		}
		default:
		{
			return 0;
			break;
		}
		}*/
	}
	void NLoad::getLU(int index, int iyr, int& LUcS, int& LUcE, double& perc) {
		switch (loadType)
		{
		case mantisServer::LoadType::GNLM:
		{
			if (iyr < 1945) {
				LUcS = getLU(index,0);
				LUcE = getLU(index,0);
				perc = 1.0;
			}
			else if (iyr >= 2005) {
				LUcS = getLU(index,4);
				LUcE = getLU(index,4);
				perc = 0.0;
			}
			else {
				int isN = static_cast<int>( std::floor((iyr - 1945) / 15) );
				int ieN = isN + 1;
				LUcS = getLU(index,isN);
				LUcE = getLU(index,ieN);
				perc = static_cast<double>((iyr - 1945) % 15) / 15.0;
			}
			break;
		}
		case mantisServer::LoadType::SWAT:
		{
			std::cout << "You cant call this method for SWAT loading type" << std::endl;
			break;
		}
		default:
			break;
		}
	}

	void NLoad::getNload(int index, int iyr, double& N1, double& N2, double& u) {
		switch (loadType)
		{
		case mantisServer::LoadType::GNLM:
		{
			//std::cout << " index=" << index << " iyr=" << iyr;
			if (iyr < 1945) {
				N1 = Ndata[index][0];
				N2 = Ndata[index][0];
				u = 1.0;
			}
			else if (iyr >= 2050) {
				N1 = Ndata[index][7];
				N2 = Ndata[index][7];
				u = 0.0;
			}
			else {

				int isN = static_cast<int>(std::floor((iyr - 1945) / 15)); //index of starting year
				int ieN = isN + 1; // index of starting year
				//std::cout << " isN=" << isN << " ieN=" << ieN;
				u = static_cast<double>((iyr - 1945) % 15) / 15.f;
				N1 = Ndata[index][isN];
				N2 = Ndata[index][ieN];
				//std::cout << " u=" << u;
				//std::cout << " N1=" << N1 << " N2=" << N2;
			}
			break;
		}
		case mantisServer::LoadType::SWAT:
		{
			std::cout << "You can't call this method for SWAT loading type" << std::endl;
			break;
		}
		default:
			break;
		}
	}

	float NLoad::getNload(int index, int iyr) {
		float value = 0.0;
		switch (loadType)
		{
		case mantisServer::LoadType::GNLM:
		{
			std::cout << "You can't call this method for GNLM loading type" << std::endl;
			break;
		}
		case mantisServer::LoadType::SWAT:
		{
			int load_index = (iyr - 1940) % 25;
			value = Ndata[index][load_index];
			break;
		}
		default:
			break;
		}
		return value;
	}

	bool NLoad::readData(std::string filename, LoadType ltype) {
        const std::string LUNameSet("LU");
        const std::string NidxNameSet("Nidx");
        const std::string NloadNameSet("NLoad");

        hf::File HDFNfile(filename, hf::File::ReadOnly);
        hf::DataSet datasetLU = HDFNfile.getDataSet(LUNameSet);
        hf::DataSet datasetNidx = HDFNfile.getDataSet(NidxNameSet);
        hf::DataSet datasetNLoad = HDFNfile.getDataSet(NloadNameSet);
        datasetLU.read(LU);
        datasetNidx.read(Nidx);
        datasetNLoad.read(Ndata);
        loadType = ltype;

        /*

		std::ifstream ifile;
		ifile.open(filename);
		if (!ifile.is_open()) {
			std::cout << "Cant open file: " << filename << std::endl;
			return false;
		}
		else {
			std::string line;
			int Nlu = 0;
			int Nl = 0;
			
			switch (ltype)
			{
			case mantisServer::LoadType::GNLM:
				Nlu = 5;
				Nl = 8;
				break;
			case mantisServer::LoadType::SWAT:
				Nlu = 1;
				Nl = 25;
				break;
			default:
				break;
			}
			std::vector<int> lu(Nlu);
			std::vector<float> NGWload(Nl);
			int code;
			while (getline(ifile, line)) {
				std::istringstream inp(line.c_str());
				inp >> code;
				for (int i = 0; i < Nlu; i++) {
					inp >> lu[i];
				}
				for (int i = 0; i < Nl; i++) {
					inp >> NGWload[i];
				}
				Ndata.push_back(NGWload);
				LU.push_back(lu);
			}
			ifile.close();
		}
		*/
		return true;
	}

	void NLoad::buildLoadingFunction(int index, int endYear, std::vector<double>& LF, Scenario& scenario, double mult) {
		
		int startYear = 1945;
		int istartReduction = scenario.startReductionYear - startYear;
		int iendReduction = scenario.endReductionYear - startYear;
		int Nyears = endYear - startYear;
		double adoptionCoeff = 0;
		LF.resize(Nyears, 0.0);
		double percReduction = 1.0;
		std::map<int, double>::iterator it;
		if (loadType == LoadType::SWAT) {
			int lucode = getLU(index, 0);
			it = scenario.LoadReductionMap.find(lucode);
			if (it != scenario.LoadReductionMap.end()) {
				percReduction = it->second;
				//std::cout << "pR=" << percReduction << std::endl;
			}
		}

		for (int i = 0; i < Nyears; i++) {
			//std::cout << "i=" << i << " Y=" << i + startYear << std::endl;
			if ((i >= istartReduction) & (i <= iendReduction))
				adoptionCoeff = (static_cast<double>(i) - static_cast<double>(istartReduction)) / (static_cast<double>(iendReduction) - static_cast<double>(istartReduction));
			else if (i > iendReduction)
				adoptionCoeff = 1.0;

			//std::cout << "a=" << adoptionCoeff;

			if (loadType == LoadType::GNLM) {
				int lus = 0;
				int lue = 0;
				double prc = 0.0;
				double rs = 1.0;
				double re = 1.0;
				getLU(index, i + startYear, lus, lue, prc);
				//std::cout << " lus=" << lus << " lue=" << lue << " prc=" << prc;
				if (adoptionCoeff > 0) {
					it = scenario.LoadReductionMap.find(lus);
					if (it != scenario.LoadReductionMap.end()) {
						rs = it->second;
						//std::cout << " rs=" << rs;
					}
					it = scenario.LoadReductionMap.find(lue);
					if (it != scenario.LoadReductionMap.end()) {
						re = it->second;
						//std::cout << " rs=" << rs;
					}
				}
				double N1 = 0, N2 = 0, u = 0;
				getNload(index, i + startYear, N1, N2, u);

				double Nbase = N1* (1 - u) + N2* u;
				//std::cout << " Nbase=" << Nbase;
				
				if ((adoptionCoeff > 0) & ((std::abs(1 - rs) > 0.000000001) | (std::abs(1 - re) > 0.000000001))) {
					double Nred = (N1 * rs) * (1 - u) + (N2 * re) * u;
					//std::cout << " Nred=" << Nred;
					LF[i] = (Nbase * (1 - adoptionCoeff) + Nred * adoptionCoeff) * mult;
				}
				else {
					LF[i] = Nbase * mult;
				}
				//std::cout << " LF[i]=" << LF[i] << std::endl;
			}
			else if (loadType == LoadType::SWAT) {
				double Nbase = static_cast<double>(getNload(index, i + startYear));
				//std::cout << "Index " << index;
				//std::cout << " Nbase " << Nbase;
				double Nred = percReduction * Nbase;
				//std::cout << " Nred " << Nred;
				LF[i] = (Nbase * (1 - adoptionCoeff) + Nred * adoptionCoeff) * mult;
				//std::cout << " LF[i] " << LF[i] << std::endl;
			}
		}
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
		//std::vector<ine_Poly_2> polys;

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

		//void simulate(std::string &msg, std::string &outmsg);

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


	private:
		mantisServer::options options;

		/*!
		 * CVraster is a 2D raster map that contains the indices of the linear vector.
		 * All the info that is formatted as 1D or 2D arrays with length 1 - Ncells
		 * must use the CVraster indices to indentify the associated record.
		 * For example LU[0][CVraster[col][row]].
		 * Note however that CVraster contains negative numbers. which correspond to the
		 * cells that have no associated info. These are outside the study area.
		 */
		std::vector<std::vector<int>> CVraster;

		std::map<std::string, NLoad> NGWLoading;

		//! LU is a 3D array of Nrow x Ncol x Years (5)
		//std::vector< std::vector< std::vector<int> > > LU;
		//! NGW is a 3D array of Nrow x Ncol x Years (8)
		//std::vector< std::vector< std::vector<float> > > NGW;
		//! UNSAT is a 2D array of Nrow x Ncol. If we decide to include more than one travel time map then this should be a 3D array
		std::map<std::string, int > UNSATscenarios;
        std::vector<std::vector<double>> UNSATData;
		//! a list of background maps
		//std::map<int, std::map<int, Polyregion> > MAPList;
		std::map<std::string, std::map<std::string, Polyregion> > MAPList;
		std::vector<std::string> backgroundMapNames;

		//! A map of well ids and wellClass
		std::map<std::string, std::map<int, wellClass> > Wellmap;

		
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

		bool readUNSAT();

		bool readCVraster();

		//void assign_point_in_sets(double x, double y, int wellid, std::string setname);
		
		bool buildLoadingFunction(Scenario &scenario, std::vector<double> &LF, int row, int col, double rch);

		float unsatTravelTime(std::map<std::string, int >::iterator& it, int gnlm_index);
	};



	Mantis::Mantis(mantisServer::options options_in)
		:
		options(options_in)
	{}

	bool Mantis::validate_msg(std::string& outmsg) {
		if ((scenario.endSimulationYear <= 1990) | (scenario.endSimulationYear > 2500))
			scenario.endSimulationYear = 2500;
		if (scenario.startReductionYear < 1945)
			scenario.startReductionYear = 2020;
		if (scenario.startReductionYear >= scenario.endReductionYear)
            scenario.endReductionYear = scenario.startReductionYear + 1;
		if (scenario.endReductionYear > scenario.endSimulationYear)
			scenario.endReductionYear = scenario.endSimulationYear;

		std::map<std::string, std::map<std::string, Polyregion> >::iterator mapit = MAPList.find(scenario.mapID);
		if (mapit == MAPList.end()) {
			outmsg += "0 ERROR: The Background map with id [";
			outmsg += scenario.mapID;
			outmsg += "] could not be found";
			return false;
		}
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

			wellscenit = regit->second.wellids.find(scenario.flowScen);
			if (wellscenit == regit->second.wellids.end()) {
				outmsg += "0 ERROR: There is no scenario with name: ";
				outmsg += scenario.flowScen;
				return false;
			}
		}

		std::map<std::string, std::map<int, wellClass> >::iterator wellscenNameit = Wellmap.find(scenario.flowScen);
		if (wellscenNameit == Wellmap.end()) {
			outmsg += "0 ERROR: There are no wells and urfs for the scenario with name: ";
			outmsg += scenario.flowScen;
			return false;
		}
		return true;
	}

	bool Mantis::parse_incoming_msg(std::string &msg, std::string &outmsg) {
		bool out = false;
		scenario.clear();
		// convert the string to a stringstream
		std::istringstream ss(msg);
		//ss << msg;
		
		int counter = 0;
		while (true) {
			std::string test;
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
			if (test == "loadScen") {
				ss >> scenario.loadScen;
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

			if (test ==  "DebugID") {
				scenario.printAdditionalInfo = true;
				int debugid;
				ss >> scenario.debugID;
				continue;
			}

			// The incoming messages express the reduction in loading
			// If the reduction is 0.2 then the Nloading 80% less compared to base case
			if (test == "Ncrops") {
				int Ncrops, cropid;
				double perc;
				ss >> Ncrops;
				for (int i = 0; i < Ncrops; i++) {
					ss >> cropid;
					ss >> perc;
					if (cropid == -9){
					    scenario.globalMod = perc;
					}
					else{
                        scenario.LoadReductionMap.insert(std::pair<int, double>(cropid, perc));
					}
				}
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
		}
		return out;


		/*
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
		std::cout << "Simulation scenario name: " << scenario.name << std::endl;
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
		*/
	}

	bool Mantis::readBackgroundMaps() {
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
					int Npoly; //number of polygons for this region

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
	}

	bool Mantis::readInputs() {

        bool tf = readCVraster();
        if (!tf) { std::cout << "Error reading Raster file" << std::endl; return false; }

		tf = readBackgroundMaps();
		if (!tf) { std::cout << "Error reading Background Maps" << std::endl; return false; }

        tf = readUNSAT();
        if (!tf) { std::cout << "Error reading UNSAT" << std::endl; return false; }

        tf = readLU_NGW();
        if (!tf) { std::cout << "Error reading LU or NGW" << std::endl; return false; }

        tf = readMultipleSets(options.WELLfile, true);
        if (!tf) { std::cout << "Error reading Wells" << std::endl; return false; }

		//if (!options.testMode) {

		//}

		//if (!options.testMode) {
		//	tf = readLU_NGW();
		//	if (!tf) { std::cout << "Error reading LU or NGW" << std::endl; return false; }
		//}
		
		//if (!options.testMode) {

		//}

		//if (!options.testMode) {
		tf = readMultipleSets(options.URFfile, false);
		if (!tf) { std::cout << "Error reading URFs" << std::endl; return false; }
		//}
		return tf;
	}

	/*
	void Mantis::assign_point_in_sets(double x, double y, int wellid, std::string setname) {
		//std::map<int, std::map<int, Polyregion> >::iterator mapit;
		std::map<std::string, std::map<std::string, Polyregion> >::iterator mapit;
		std::map<std::string, Polyregion>::iterator regit;
		//std::vector<ine_Poly_2>::iterator polyit;

		//ine_Point2 testPoint(x, y);
		for (mapit = MAPList.begin(); mapit != MAPList.end(); ++mapit) {
			bool found = false;
			for (regit = mapit->second.begin(); regit != mapit->second.end(); ++regit) {
				//for (polyit = regit->second.polys.begin(); polyit != regit->second.polys.end(); ++polyit) {
				//	switch (polyit->bounded_side(testPoint)) {
				//	case CGAL::ON_BOUNDED_SIDE:
				//		regit->second.wellids[setname].push_back(wellid);
				//		found = true;
				//		break;
				//	}
				//}
				if (found)
					break;
			}
		}
	}
	*/

	bool Mantis::readMultipleSets(std::string filename, bool isWell) {
		std::ifstream WellMasterfile;

		if (!options.bAbsolutePaths)
		    filename = options.mainPath + filename;

		WellMasterfile.open(filename);
		if (!WellMasterfile.is_open()) {
			std::cout << "Cant open file: " << filename << std::endl;
			return false;
		}
		else {
			std::string line, filename;
			while (getline(WellMasterfile, line)) {
                std::istringstream inp(line.c_str());
                inp >> filename;
			    if (filename.empty())
			        continue;
			    if (filename.front() == '#')
			        continue;

                if (!options.bAbsolutePaths)
                    filename = options.mainPath + filename;
				bool tf;
				if (isWell)
                    tf = readWellSet(filename);
				else
                    tf = readURFs(filename);

				if (!tf) {
					std::cout << "An error occured while reading " << filename << std::endl;
					return false;
				}
			}
		}
		return true;
	}

	bool Mantis::readWellSet(std::string filename) {
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
			std::string setName, line;
            getline(Welldatafile, line);
            std::istringstream inp(line.c_str());
			inp >> Nwells;
			inp >> setName;
			double xw, yw, D, SL, Q;
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
                            std::cout << "Cannot find region " <<  regionCode << " in background map " << backgroundMapNames[j] << std::endl;
                        }
                    }
                    else{
                        std::cout << " Cannot find background map with name : " << backgroundMapNames[j] << std::endl;
                    }

				    if (j == 0){

				    }
				}
				//assign_point_in_sets(xw, yw, Eid, setName);
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
        const std::string NamesNameSet("Names");
        const std::string IntsNameSet("ESIJRiv");
        const std::string FloatNameSet("MSWR");

        hf::File HDFfile(filename, hf::File::ReadOnly);

        hf::DataSet datasetNames = HDFfile.getDataSet(NamesNameSet);
        hf::DataSet datasetInts = HDFfile.getDataSet(IntsNameSet);
        hf::DataSet datasetFloat = HDFfile.getDataSet(FloatNameSet);
        std::vector<std::string> names;
        std::vector<std::vector<int>> IDS;
        std::vector<std::vector<double>> DATA;
        datasetNames.read(names);
        datasetInts.read(IDS);
        datasetFloat.read(DATA);

        URFTYPE urftype;
        if (names[1].compare("LGNRM") == 0)
            urftype = URFTYPE::LGNRM;
        else if (names[1].compare("ADE") == 0)
            urftype = URFTYPE::ADE;
        else if (names[1].compare("BOTH") == 0)
            urftype = URFTYPE::BOTH;
        else{
            std::cout << "The URF type " << names[1] << " is not valid" << std::endl;
            return false;
        }

        std::map<std::string, std::map<int, wellClass> >::iterator scenit;
        scenit = Wellmap.find(names[0]);
        if (scenit == Wellmap.end()) {
            std::cout << "The URF Scenario " << names[0] << " is not defined for the wells" << std::endl;
            return false;
        }
        else{
            std::map<int, wellClass>::iterator wellmapit;
            for (unsigned int i = 0; i < IDS[0].size(); ++i){
                wellmapit = scenit->second.find(IDS[0][i]);
                if (wellmapit != scenit->second.end()){
                    wellmapit->second.addStreamline(IDS[1][i], IDS[2][i], IDS[3][i],
                                                    DATA[2][i], DATA[3][i],urftype, IDS[4][i],
                                                    DATA[0][i],DATA[1][i]);
                }
            }
        }


        /*
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
			std::string setName, urftypestring, line;
			URFTYPE urftype;
			int Nurfs;

            getline(URFdatafile, line);
            std::istringstream inp(line.c_str());
            inp >> Nurfs;
            inp >> setName;
            inp >> urftypestring;

			// Identify urf type
            if (urftypestring == "LGNRM")
                urftype = URFTYPE::LGNRM;
            else if (urftypestring =="ADE")
                urftype = URFTYPE::ADE;
            else if (urftypestring =="BOTH")
                urftype = URFTYPE::BOTH;
            else{
                std::cout << "The URF type " << urftypestring << " is not valid" << std::endl;
                return false;
            }

            scenit = Wellmap.find(setName);
			if (scenit == Wellmap.end()) {
				std::cout << "The URF Scenario " << setName << " is not defined for the wells" << std::endl;
				return false;
			}
			else {
				int Eid, Sid, gnlm_ind, swat_ind;
				double paramA, paramB, w, rch, paramC, paramD;
				paramC = 0.0;
				paramD = 0.0;
				for (int i = 0; i < Nurfs; ++i) {
                    getline(URFdatafile, line);
                    std::istringstream inp1(line.c_str());
					inp1 >> Eid;
					inp1 >> Sid;
					inp1 >> gnlm_ind;
					inp1 >> swat_ind;
					inp1 >> paramA;
					inp1 >> paramB;
					inp1 >> w;
					inp1 >> rch;
					if (urftype == URFTYPE::BOTH){
						inp1 >> paramC;
						inp1 >> paramD;
					}

					//std::cout << Eid << ", " << Sid << ", " << ROW << ", " << COL << ", " << w << ", " << paramA << ", " << paramB << std::endl;
					//if (Eid > 699) {
					//	bool breakHere = true;
					//}
					//Eid++;// This is used since the particle tracking code writes the entities with zero based numbering
					wellmapit = scenit->second.find(Eid);
					if (wellmapit != scenit->second.end()) {
						if (swat_ind != -9)
							swat_ind = swat_ind - 1;

						wellmapit->second.addStreamline(Sid, gnlm_ind - 1, swat_ind, w, rch, urftype, paramA, paramB, paramC, paramD);
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
		*/
		return true;
	}

	/*

	void Mantis::simulate(std::string &msg, std::string &outmsg) {
		auto start = std::chrono::high_resolution_clock::now();
		outmsg.clear();
		//Scenario scenario;
		//parse_incoming_msg(msg, scenario);
		
		// Find the selected background map
		std::map<std::string, std::map<std::string, Polyregion> >::iterator mapit = MAPList.find(scenario.mapID);
		if (mapit == MAPList.end()) {
			outmsg += "ERROR: The Background map with id [";
			outmsg += scenario.mapID;
			outmsg += "] could not be found";
			return;
		}

		std::map<std::string, Polyregion>::iterator regit;
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

			wellscenit = regit->second.wellids.find(scenario.flowScen);
			if (wellscenit == regit->second.wellids.end()) {
				outmsg += "ERROR: There is no scenario with name: ";
				outmsg += scenario.flowScen;
				return;
			}

		}

		std::map<std::string, std::map<int, wellClass> >::iterator wellscenNameit = Wellmap.find(scenario.flowScen);
		if (wellscenNameit == Wellmap.end()) {
			outmsg += "ERROR: There are no wells and urfs for the scenario with name: ";
			outmsg += scenario.flowScen;
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
			wellscenit = regit->second.wellids.find(scenario.flowScen);

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
								//buildLoadingFunction(scenario, LF, strmlnit->second.row-1, strmlnit->second.col-1);
								urf.convolute(LF, BTC);
							}
							else if (urftp == URFTYPE::ADE){
								URF urf(options.nSimulationYears, strmlnit->second.sl, strmlnit->second.vel, urftp,ade_opt);
								//buildLoadingFunction(scenario, LF, strmlnit->second.row-1, strmlnit->second.col-1);
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

	*/

	bool Mantis::buildLoadingFunction(Scenario &scenario, std::vector<double> &LF, int row, int col, double rch) {
	    int lin_idx = CVraster[col][row];
		std::map<std::string, NLoad>::iterator loadit = NGWLoading.find(scenario.loadScen);
		bool out = false;
		switch (loadit->second.getLtype())
		{
		case LoadType::GNLM:
		{
			if (loadit->second.isValidIndex(lin_idx)) {
				double mult = 100.0 / rch;
				loadit->second.buildLoadingFunction(lin_idx, scenario.endSimulationYear, LF, scenario, mult);
				out = true;
				
			}
			break;
		}
		case LoadType::SWAT:
		{
			if (loadit->second.isValidIndex(lin_idx)) {
				loadit->second.buildLoadingFunction(lin_idx, scenario.endSimulationYear, LF, scenario, 1);
				out = true;
			}
			//std::map<std::string, NLoad>::iterator gnlmit = NGWLoading.find("GNLM");
			//std::vector<double> temp_LF;
			//double mult = 100.0 / rch;
			//gnlmit->second.buildLoadingFunction(gnlm_index, 1990, temp_LF, scenario, mult);
			//Then we have to join those two 
			break;
		}
		default:
			break;
		}
		return out;


		/*
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
		*/
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

	float Mantis::unsatTravelTime(std::map<std::string, int >::iterator &it, int gnlm_index) {

        return 0.0;
		//return it->second[gnlm_index];
	}

	bool Mantis::readCVraster() {
        const std::string NameSet("Raster");
        if (!options.bAbsolutePaths)
            options.CVrasterFile = options.mainPath + options.CVrasterFile;

        hf::File HDFfile(options.CVrasterFile, hf::File::ReadOnly);
        hf::DataSet dataset = HDFfile.getDataSet(NameSet);
        dataset.read(CVraster);
        return true;
	}

	bool Mantis::readUNSAT() {
        const std::string NamesNameSet("Names");
        const std::string DataNameSet("Data");

        if (!options.bAbsolutePaths)
            options.UNSATfile = options.mainPath + options.UNSATfile;

        hf::File HDFfile(options.UNSATfile, hf::File::ReadOnly);

        hf::DataSet datasetNames = HDFfile.getDataSet(NamesNameSet);
        hf::DataSet datasetData = HDFfile.getDataSet(DataNameSet);
        std::vector<std::string> names;
        datasetNames.read(names);
        datasetData.read(UNSATData);
        for (unsigned int i = 0; i < names.size(); ++i){
            UNSATscenarios.insert(std::pair<std::string, int>(names[i],i));
        }
        return true;
        /*
		auto start = std::chrono::high_resolution_clock::now();
		std::ifstream unsatfile;



		unsatfile.open(options.UNSATfile);
		if (!unsatfile.is_open()) {
			std::cout << "Cant open file: " << options.UNSATfile << std::endl;
			return false;
		}
		else {
			std::string line;
			bool readHeader = true;
			//std::vector<float> emptyvector(options.gnlmNpixels,0);
			std::vector<std::string> ScenarioNames;
			int index = 0;
			while (getline(unsatfile, line)) {
				std::istringstream inp(line.c_str());
				if (readHeader) {
					while (true) {
						std::string scenName;
						inp >> scenName;
						if (scenName.empty())
							break;
						else {
							ScenarioNames.push_back(scenName);
							UNSAT.insert(std::pair < std::string, std::vector<float> >(scenName, std::vector<float>(options.gnlmNpixels,0)));
						}
					}
					readHeader = false;
				}
				else {
					float v;
					std::map<std::string, std::vector<float> >::iterator it;

					for (unsigned int i = 0; i < ScenarioNames.size(); ++i) {
						inp >> v;
						it = UNSAT.find(ScenarioNames[i]);
						it->second.at(index) = v;
					}
					index++;
				}
			}
		}
		unsatfile.close();
		return true;
         */
	}

	bool Mantis::readLU_NGW() {
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
				std::string Ltype, Lname, Lfile;
				std::istringstream inp(line.c_str());
                inp >> Ltype;

				if (Ltype.empty())
				    continue;

				if (Ltype.front() == '#')
                    continue;

				inp >> Lname;
				inp >> Lfile;

                if (!options.bAbsolutePaths)
                    Lfile = options.mainPath + Lfile;

				if (Ltype.compare("GNLM") == 0) {
					NLoad NL;
					NL.readData(Lfile, LoadType::GNLM);
					NGWLoading.insert(std::pair<std::string, NLoad>(Lname, NL));
				}
				else if (Ltype.compare("SWAT") == 0) {
					NLoad NL;
					NL.readData(Lfile, LoadType::SWAT);
					NGWLoading.insert(std::pair<std::string, NLoad>(Lname, NL));
				}
				else {
					std::cout << "Unknown loading type " << Ltype << std::endl;
					return false;
				}
			}
		}


		/*
		// allocate memory
		//12863 x 7046 x 5 
		LU.resize(options.Nrow, std::vector< std::vector<int> >(options.Ncol, std::vector<int>(5, 0)));
		NGW.resize(options.Nrow, std::vector< std::vector<float> >(options.Ncol, std::vector<float>(8, 0)));
		UNSAT.resize(options.Nrow, std::vector<float>(options.Ncol, 0));
		std::ifstream LUdatafile, NGWdatafile, UNSATdatafile;
		
		LUdatafile.open(options.GNLM_LUfile);
		if (!LUdatafile.is_open()) {
			std::cout << "Cant open file: " << options.GNLM_LUfile << std::endl;
			return false;
		}

		NGWdatafile.open(options.GNLM_NGWfile);
		if (!NGWdatafile.is_open()) {
			std::cout << "Cant open file: " << options.GNLM_NGWfile << std::endl;
			return false;
		}

		UNSATdatafile.open(options.UNSATfile);
		if (!UNSATdatafile.is_open()) {
			std::cout << "Cant open file: " << options.UNSATfile << std::endl;
			return false;
		}
		std::cout << "Reading " << options.GNLM_LUfile << std::endl;
		std::cout << "Reading " << options.GNLM_NGWfile << std::endl;
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
		*/
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		std::cout << "Read LU and NGW in " << elapsed.count() << std::endl;
		return true;
	}


	void Mantis::simulate_with_threads(int id) {//, , std::string &outmsg

		// Get an iterator to the selected map
		std::map<std::string, std::map<std::string, Polyregion> >::iterator mapit = MAPList.find(scenario.mapID);
		// Get an iterator to the list of wells for the selected scenario
		std::map<std::string, std::map<int, wellClass> >::iterator wellscenNameit = Wellmap.find(scenario.flowScen);

		std::map<std::string, Polyregion>::iterator regit;
		std::map<std::string, std::vector<int> >::iterator wellscenit;
		std::map<std::string, int >::iterator unsatit = UNSATscenarios.find(scenario.unsatScenario);

		

		//for (int i = startWell; i < endWell; ++i) {
		//	std::cout << id << ":" << i << std::endl;
		//}

		std::map<int, streamlineClass>::iterator strmlnit;
		std::map<int, wellClass>::iterator wellit;
		int cntBTC = 0;
		for (int irg = 0; irg < static_cast<int>(scenario.regionIDs.size()); ++irg) {
			regit = mapit->second.find(scenario.regionIDs[irg]);
			wellscenit = regit->second.wellids.find(scenario.flowScen);

			// Number of wells in the selected region
			int Nwells = static_cast<int>(wellscenit->second.size());
			if (Nwells == 0)
				continue;

			int startWell, endWell;

			if (Nwells <= options.nThreads) {
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
			
			int NsimulationYears = scenario.endSimulationYear - 1945;
			
			

			std::cout << "Thread " << id << " will simulate from [" << startWell << " to " << endWell << ")" << std::endl;
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
								bool isNotZero = true;
								if (std::abs(strmlnit->second.mu - 0) > 0.00000001) {								
									//if (iw >= 7337 && iw < 7338) {
									//	std::cout << iw << ":" << strmlnit->second.mu << " " << strmlnit->second.std << " " << strmlnit->second.w << std::endl;
									//}
									URF urf(NsimulationYears, strmlnit->second.mu, strmlnit->second.std, strmlnit->second.type);
										
									isNotZero = buildLoadingFunction(scenario, LF, strmlnit->second.row, strmlnit->second.col, strmlnit->second.gwrch);
									//if (iw >= 7337 && iw < 7338)
									//	printVector<double>(LF, "LF");
									if (isNotZero) {
										urf.convolute(LF, BTC);
									}
									//if (iw >= 7337 && iw < 7338)
									//	printVector<double>(BTC, "BTC");
								}
								else if (strmlnit->second.type == URFTYPE::ADE){
									URF urf(NsimulationYears, strmlnit->second.mu, strmlnit->second.std, strmlnit->second.type, ADEoptions());
							///		isNotZero = buildLoadingFunction(scenario, LF, strmlnit->second.gnlm_index, strmlnit->second.swat_index, strmlnit->second.gwrch);
									urf.convolute(LF, BTC);
								}

								if (isNotZero) {
									// Find the travel time in the unsaturated zone
									int intTau = 0;
							///		if (unsatit != UNSATscenarios.end()) {
							///			double tau = static_cast<double>( unsatTravelTime(unsatit, strmlnit->second.gnlm_index));
							///			tau = std::floor(tau * scenario.unsatZoneMobileWaterContent);
							///			if (tau < 0)
							///				tau = 0.0;
							///			intTau = static_cast<int>(tau);
										//std::cout << tau << std::endl;
							///		}

									if (intTau < NsimulationYears) {
										int ibtc = 0;
										for (int ii = intTau; ii < NsimulationYears; ++ii) {
											//std::cout << ii << std::endl;
											weightBTC[ii] = weightBTC[ii] + BTC[ibtc] * strmlnit->second.w;
											ibtc++;
										}
									}
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
							for (int iwbtc = 0; iwbtc < NsimulationYears; ++iwbtc) {
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
		int Nyears = scenario.endSimulationYear - 1945;
		outmsg += " " + std::to_string(Nyears);

		for (int i = 0; i < static_cast<int>(replymsg.size()); ++i) {
			outmsg += " ";
			outmsg += replymsg[i];
		}
		// This is the character that indicates the end of the message
		outmsg += " ENDofMSG\n";

		//std::cout << outmsg << std::endl;
		std::cout << nBTC << "BTCs will be sent" << std::endl;
	}
}
