//#pragma once

#ifndef MANTISSERVER_MSHELPER_H
#define MANTISSERVER_MSHELPER_H

#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/lognormal_distribution.hpp>
boost::random::mt19937 gen;


//#include <boost/geometry.hpp>
//#include <boost/geometry/geometries/point_xy.hpp>
//#include <boost/geometry/geometries/polygon.hpp>
//#include <boost/geometry/algorithms/within.hpp>


#if _USEHF > 0
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#endif

// Until I find a better solution for the output log
// I will use this variable as global for redirecting the console output
//https://stackoverflow.com/questions/10150468/how-to-redirect-cin-and-cout-to-files
std::ofstream logStream;

namespace mantisServer {

    //! Set the square pi as constant
    const double sqrt2pi = std::sqrt(2*std::acos(-1));
    //! Set the pi as constant
    const double pi = std::atan(1)*4;

    double deg2rad(double d){
        return d*pi/180.0;
    }
    double cosd(double dang){
        return std::cos(deg2rad(dang));
    }
    double sind(double dang){
        return std::sin(deg2rad(dang));
    }

    /*!
     * num2Padstr converts an integer to string with padding zeros
     * @param i is the integer to convert to string
     * @param n is the number of zeros
     * @return a string. For example num2Padstr(3,3) returns 003
     */
    std::string num2Padstr(int i, int n) {
        std::stringstream ss;
        ss << std::setw(n) << std::setfill('0') << i;
        return ss.str();
    }

    /*!
	 * This prints the vectors to screen formatted for matlab.
	 * Simply copy paste the printed line to matlab workspace to create the variable
	 * @tparam T is the vector type, integer or double
	 * @param v this is the vector.
	 * @param varname is what the variable name should be in matlab
	 */
    //template<typename T>
    //void printVector(std::vector<T>& v, std::string varname) {
    //    std::cout << std::endl;
    //    std::cout << varname << " = [";
    //    for (unsigned int i = 0; i < v.size(); ++i) {
    //        std::cout << v[i] << " ";
    //    }
    //    std::cout << "];" << std::endl;
    //    std::cout << std::endl;
    //}

    //typedef boost::geometry::model::d2::point_xy<double> boost_point;
    //typedef boost::geometry::model::polygon<boost_point> boost_poly;


    /**
	 * @brief Enumeration for the type of the Unit Response function.
	 *
	 * Since it is very inefficient to hold all 200~300 double values we will keep only
	 * a number of parameters.
	 *
	 * If the URF is represented as lognormal distribution the parameteres are the mean and the
	 * standard deviation. This is what we currently used. THe other two type will be removed in future versions.
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
        ADE, /**< The URF is represented as analytical ADE with Length and velocity parameters*/
        UNKNOWN
    };

    enum class LoadUnits{
        CONC,
        MASS,
        UNKNOWN
    };

    URFTYPE string2URFTYPE(std::string str){
        if (str.compare("LGNRM") == 0){
            return URFTYPE::LGNRM;
        }
        else if (str.compare("ADE") == 0){
            return URFTYPE::ADE;
        }
        else
            return URFTYPE::UNKNOWN;
    }



    LoadUnits string2LoadUnits(std::string str){
        if (str.compare("CONC") == 0)
            return LoadUnits::CONC;
        else if (str.compare("MASS") == 0)
            return LoadUnits::MASS;
        else
            return LoadUnits::UNKNOWN;
    }

    enum class RasterOperation{
        Multiply,
        Replace,
        DONTUSE,
        UNKNOWN
    };
    RasterOperation string2RasterOperation(std::string str){
        if (str.compare("Multiply") == 0)
            return RasterOperation::Multiply;
        else if (str.compare("Replace") == 0)
            return RasterOperation::Replace;
        else
            return RasterOperation::UNKNOWN;
    }

    struct cell{
        cell(){
            row = -1;
            col = -1;
        }
        cell(int r, int c){
            row = r;
            col = c;
        }
        int row;
        int col;
        int lin_ind;
        double dist;
    };

    bool compareCellByDistance(const cell &a, const cell &b){
        return a.dist < b.dist;
    }

    std::vector<cell> SearchPattern(){
        std::vector<cell> c;
        c.push_back(cell(-1,-1));
        c.push_back(cell(-1,0));
        c.push_back(cell(-1,1));
        c.push_back(cell(0,-1));
        c.push_back(cell(0,1));
        c.push_back(cell(1,-1));
        c.push_back(cell(1,0));
        c.push_back(cell(1,1));
        return c;
    }

    void logMSG1(std::ostream & out, std::string msg){
        out << msg << std::endl;
    }

    std::string getExtension(std::string filename){
        // https://www.fluentcpp.com/2017/04/21/how-to-split-a-string-in-c/
        // Solution 1.3: stepping away from the iterators
        std::vector<std::string> tokens;
        std::string token;
        std::istringstream tokenStream(filename);
        while (std::getline(tokenStream, token, '.')){
            tokens.push_back(token);
        }
        return tokens.back();
    }

    bool isInTriangle(double ax, double ay, double bx, double by, double cx, double cy, double px, double py){
        // https://ceng2.ktu.edu.tr/~cakir/files/grafikler/Texture_Mapping.pdf page 47
        double v0x = bx - ax; double v0y = by - ay;
        double v1x = cx - ax; double v1y = cy - ay;
        double v2x = px - ax; double v2y = py - ay;

        double d00 = v0x*v0x + v0y*v0y;
        double d01 = v0x*v1x + v0y*v1y;
        double d11 = v1x*v1x + v1y*v1y;
        double d20 = v2x*v0x + v2y*v0y;
        double d21 = v2x*v1x + v2y*v1y;
        double denom = d00*d11 - d01*d01;


        if (std::abs(denom) < 0.0000000001){
            return false;
        }
        double u = (d11*d20 - d01*d21)/denom;
        double v = (d00*d21 - d01*d20)/denom;
        double w = 1.0 - u - v;

        if (u > -0.0001 && u < 1.0001 &&
            v > -0.0001 && v < 1.0001 &&
            w > -0.0001 && w < 1.0001){
            return true;
        }
        else{
            return false;
        }
    }

    class numericRange{
    public:
        numericRange(){};
        bool isInRange(double x);
        void setData(double xmn, double xmx);
        void clear();
        bool canUse(){return hasValues;}
    private:
        double xmin;
        double xmax;
        bool hasValues = false;
    };
    bool numericRange::isInRange(double x) {
        return x > xmin && x <= xmax;
    }
    void numericRange::setData(double xmn, double xmx) {
        xmin = xmn;
        xmax = xmx;
        hasValues = true;
    }
    void numericRange::clear() {
        xmin = 0.0;
        xmax = 0.0;
        hasValues = false;
    }

    class rectSelection{
    public:
        rectSelection(){};
        bool isPointIn(double x, double y);
        void setData(double xmn, double ymn, double xmx, double ymx);
        void clear();
    private:
        double xmin;
        double ymin;
        double xmax;
        double ymax;
    };
    bool rectSelection::isPointIn(double x, double y) {
        return x >= xmin && x <= xmax && y >= ymin && y <= ymax;
    }
    void rectSelection::setData(double xmn, double ymn, double xmx, double ymx){
        xmin = xmn;
        ymin = ymn;
        xmax = xmx;
        ymax = ymx;
    }

    void rectSelection::clear() {
        xmin = 0.0;
        ymin = 0.0;
        xmax = 0.0;
        ymax = 0.0;
    }

    class radialSelection{
    public:
        radialSelection(){};
        bool isPointIn(double x, double y);
        void setData(double x, double y, double rad);
        void clear();
    private:
        double cx;
        double cy;
        double r;
    };
    bool radialSelection::isPointIn(double x, double y){
        return (cx - x) * (cx - x) + (cy - y) * (cy - y) < r * r;
    }
    void radialSelection::setData(double x, double y, double rad){
        cx = x;
        cy = y;
        r = rad;
    }
    void radialSelection::clear() {
        cx = 0.0;
        cy = 0.0;
        r = 0.0;
    }

    class sourceArea{
    public:
        sourceArea(){};
        int getNpixels(int listSize);
        void setParameters(int nPix, int minPix, int maxPix, double percPix);
        void clear();
    private:
        int nPixels;
        int minPixels;
        int maxPixels;
        double percArea;
    };

    void sourceArea::setParameters(int nPix, int minPix, int maxPix, double percPix) {
        nPixels = nPix;
        minPixels = minPix;
        if (minPixels <=0)
            minPixels = 1;
        maxPixels = maxPix;
        percArea = percPix;
    }

    void sourceArea::clear() {
        nPixels = -9;
        minPixels = 1;
        maxPixels = -9;
        percArea = -9.0;
    }

    int sourceArea::getNpixels(int listSize) {
        if (listSize == 0){
            return 0;
        }
        int out = listSize;
        int tmp1 = -9;

        if (percArea > 0){
            out = std::ceil(static_cast<double>(listSize)*percArea);
        }
        if (nPixels > 0 && nPixels < out){
            out = nPixels;
        }
        if (minPixels > 0){
            if (out < minPixels && listSize >= minPixels)
                out = minPixels;
        }
        if (maxPixels < 0){
            return out;
        }

        if (out > maxPixels){
            out = std::min(maxPixels, listSize);
        }
        return out;
    }

    class LinearData{
    public:
        LinearData(){}
        bool readData(std::string filename, int Nr);
        int ScenarioIndex(std::string scenario);
        double getValue(int scenID, int lin_idx);
        void setNoDataValue(double v);
        double getNoDataValue(){return NoDataValue;}
        void multiply(double mult);
        void clear();
        bool hasScenario(std::string &unsatName);
    private:
        std::map<std::string, int > ScenarioMap;
        std::vector<std::vector<double>> Data;
        int Nrows;
        int Nscenarios;
        double NoDataValue;
        bool bDataValueSet = false;
        bool bHFread;
    };

    void LinearData::clear() {
        Nrows = 0;
        Nscenarios = 0;
        NoDataValue = 0;
        bDataValueSet = false;
        Data.clear();
        ScenarioMap.clear();
    }

    int LinearData::ScenarioIndex(std::string scenario) {
        std::map<std::string, int >::iterator it = ScenarioMap.find(scenario);
        if (it == ScenarioMap.end())
            return -1;
        else
            return it->second;
    }
    bool LinearData::hasScenario(std::string &unsatName){
        std::map<std::string, int >::iterator it = ScenarioMap.find(unsatName);
        return it != ScenarioMap.end();
    }

    void LinearData::multiply(double mult) {
        for (unsigned int i = 0; i < Data.size(); ++i){
            for (unsigned int j = 0; j < Data[i].size(); ++j){
                Data[i][j] = Data[i][j]*mult;
            }
        }
    }

    double LinearData::getValue(int scenID, int lin_idx) {
        if (scenID >=0 && scenID < Nscenarios){
            if (lin_idx >=0 && lin_idx < Nrows){
                //if (bHFread)
                //    return Data[lin_idx][scenID];
                //else
                return Data[scenID][lin_idx];
            }
            else {
                return NoDataValue;
            }
        }
        else{
            return  NoDataValue;
        }
    }

    void LinearData::setNoDataValue(double v) {
        NoDataValue = v;
        bDataValueSet = true;
    }

    bool LinearData::readData(std::string filename, int Nr) {
        if (!bDataValueSet){
            std::cout << "You must set up a no-data value before reading the file" << filename << std::endl;
            return false;
        }
        Nrows = Nr;
        auto start = std::chrono::high_resolution_clock::now();

#if _USEHF>0
        std::string ext = getExtension(filename);
        if (ext.compare("h5") == 0){
            try {
                const std::string NamesNameSet("Names");
                const std::string DataNameSet("Data");
                HighFive::File HDFfile(filename, HighFive::File::ReadOnly);
                HighFive::DataSet datasetNames = HDFfile.getDataSet(NamesNameSet);
                HighFive::DataSet datasetData = HDFfile.getDataSet(DataNameSet);

                std::vector<std::string> names;
                datasetNames.read(names);
                datasetData.read(Data);
                for (unsigned int i = 0; i < names.size(); ++i){
                    ScenarioMap.insert(std::pair<std::string, int>(names[i],i));
                }
                Nscenarios = ScenarioMap.size();
                auto finish = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = finish - start;
                std::cout << "Read data from " << filename << " in " << elapsed.count() << std::endl;
                bHFread = true;
                return true;
            }
            catch (...) {
                return false;
            }
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
            { //Get the number of Unsaturated scenarios
                getline(ifile, line);
                std::istringstream inp(line.c_str());
                inp >> Nscenarios;
                Data.clear();
                Data.resize(Nscenarios, std::vector<double>(Nr,NoDataValue));
            }
            {// Get the scenario names
                std::string name;
                for (int i = 0; i < Nscenarios; ++i){
                    getline(ifile, line);
                    std::istringstream inp(line.c_str());
                    inp >> name;
                    ScenarioMap.insert(std::pair<std::string, int>(name, i));
                }
            }
            {// Read the values
                for (int i = 0; i < Nrows; ++i){
                    getline(ifile, line);
                    std::istringstream inp(line.c_str());
                    double v;
                    for (int j = 0; j < Nscenarios; ++j){
                        inp >> v;
                        Data[j][i] = v;
                    }
                }
            }
            ifile.close();
        }
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "Read " << filename << " data in " << elapsed.count() << " sec" << std::endl;
        bHFread = false;
        return true;
    }




    class CVrasterClass{
    public:
        CVrasterClass(){};
        bool readData(std::string filename, int Nr, int Nc, int Ncells);
        int IJ(int row, int col);
        void getSurroundingPixels(int row, int col, int Ndepth, std::vector<int>& lin_inds);
    private:
        int Nrows;
        int Ncols;
        std::vector<std::vector<int>> raster;
        bool bHFread;
    };

    bool CVrasterClass::readData(std::string filename, int Nr, int Nc, int Ncells) {

        auto start = std::chrono::high_resolution_clock::now();
        Nrows = Nr;
        Ncols = Nc;

#if _USEHF>0
            std::string ext = getExtension(filename);
            if (ext.compare("h5") == 0){
                const std::string NameSet("Raster");
                HighFive::File HDFfile(filename, HighFive::File::ReadOnly);
                HighFive::DataSet dataset = HDFfile.getDataSet(NameSet);
                dataset.read(raster);
                if (raster.size() != Ncols){
                    std::cout << "The number of columns of the raster (" << raster.size() << ") is not equal to Ncols: " << Ncols << std::endl;
                }
                if (raster[0].size() != Nrows){
                    std::cout << "The number of rows of the raster (" << raster[0].size() << ") is not equal to Nrows: " << Nrows << std::endl;
                }
                bHFread = true;
                auto finish = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = finish - start;
                std::cout << "Read CV raster in " << elapsed.count() << " sec" << std::endl;
                return true;
            }
#endif

        bHFread = false;
        std::ifstream ifile;
        raster.clear();
        raster.resize(Nrows,std::vector<int>(Ncols,-1));
        ifile.open(filename);
        if (!ifile.is_open()){
            std::cout << "Cant open file: " << filename << std::endl;
            return false;
        }
        else{
            std::cout << "Reading " << filename << std::endl;
            std::string line;
            int r, c;
            for (int i = 0; i < Ncells; ++i){
                getline(ifile, line);
                std::istringstream inp(line.c_str());
                inp >> r;
                inp >> c;
                if (r < Nrows && c < Ncols)
                    raster[r][c] = i;
                else{
                    std::cout << "I can't assign pixel (" << r << "," << c << ") in raster map" << std::endl;
                }
            }
            ifile.close();
        }

        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "Read CV raster in " << elapsed.count() << " sec" << std::endl;
        return true;
    }

    int CVrasterClass::IJ(int row, int col) {

        if (row >=0 && col >=0 && row < Nrows && col < Ncols) {
            if (bHFread)
                return raster[col][row];
            else
                return raster[row][col];
        }
        else
            return -1;


    }

    void CVrasterClass::getSurroundingPixels(int row, int col, int Ndepth, std::vector<int> &lin_inds) {
        if (Ndepth == 0){
            int lid = IJ(row,col);
            if ( lid >= 0)
                lin_inds.push_back(lid);
        }
        else{
            for (int i = row - Ndepth; i <= row + Ndepth; ++i){
                int lid = IJ(row,col);
                if (lid == -1)
                    continue;
                for (int j = col - Ndepth; j < col + Ndepth; ++j){
                    lid = IJ(row,col);
                    if (lid == -1)
                        continue;
                    lin_inds.push_back(lid);
                }
            }
        }
    }

    enum class extrapMethod{
        NEAREST,
        REPEAT,
        UNKNOWN
    };

    extrapMethod string2xtrapMethod(std::string xs){
        if (xs.compare("NEAREST") == 0){
            return extrapMethod::NEAREST;
        }
        else if (xs.compare("REPEAT") == 0){
            return extrapMethod::REPEAT;
        }
        else{
            return extrapMethod::UNKNOWN;
        }
    }


    class TSgrid{
    public:
        TSgrid(){};
        void init(int sy, int N, int d, extrapMethod xm);
        void getIndices(int iyr, int &SY, int &EY, double &t);

    private:
        int StartYear;
        int EndYear;
        int Interval;
        int Nyears;
        double dStartYear;
        double dInterval;
        extrapMethod XM;

    };

    void TSgrid::init(int sy, int N, int d, extrapMethod xm) {
        StartYear = sy;
        Interval = d;
        Nyears = N;
        EndYear = StartYear + d * (Nyears - 1);
        dStartYear = static_cast<double>(StartYear);
        dInterval = static_cast<double>(Interval);
        XM = xm;
    }

    void TSgrid::getIndices(int iyr, int &SY, int &EY, double &t) {
        if (Nyears == 1){ // If there is only one time step return zero
            SY = 0;
            EY = 0;
            t = 0.0;
            return;
        }

        if (XM == extrapMethod::NEAREST){
            if (iyr < StartYear){
                SY = 0;
                EY = 0;
                t = 0.0;
                return;
            }
            else if (iyr > EndYear){
                SY = Nyears-1;
                EY = Nyears-1;
                t = 1.0;
                return;
            }
        }

        double dyr = static_cast<double>(iyr);
        int idx = static_cast<int>(std::floor((dyr - dStartYear)/dInterval)) + 1;
        int idx1;
        if (idx > 0){
            idx1 = idx % Nyears;
            if (idx1 == 0){
                SY = Nyears-1;
                EY = 0;
            }
            else{
                SY = idx1-1;
                EY = SY+1;
            }
        }
        else{
            idx1 = Nyears - std::abs(idx) % Nyears;
            if (idx1 == Nyears){
                SY = Nyears-1;
                EY = 0;
            }
            else{
                SY = idx1-1;
                EY = SY+1;
            }
        }

        int yr_bef = StartYear+ (idx - 1)*Interval;
        t = static_cast<double>(iyr - yr_bef)/dInterval;
    }

    void getStartEndIndices(int threadId, int nThreads, int N, int &startId, int &endId){
        if (N <= nThreads){
            if (threadId > 0){
                startId = 0;
                endId = 0;
            }
            else{
                startId = 0;
                endId = N;
            }
        }
        else{
            int n2calc = N / nThreads;
            startId = threadId * n2calc;
            endId = (threadId + 1) * n2calc;
            if (threadId == nThreads - 1)
                endId = N;
        }
    }

}

#endif //MANTISSERVER_MSHELPER_H
