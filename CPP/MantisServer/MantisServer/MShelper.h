#pragma once

#ifndef MANTISSERVER_MSHELPER_H
#define MANTISSERVER_MSHELPER_H

#include <fstream>
#include <iostream>
#include <map>

#if _USEHF > 0
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#endif


namespace mantisServer {

    enum class LoadUnits{
        CONC,
        MASS,
        UNKNOWN
    };

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
    };

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

    bool isInTriangle(double x1, double y1, double x2, double y2, double x3, double y3, double x, double y){
        double det = (y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3);
        if (std::abs(det) < 0.0000000001){
            return false;
        }
        double bc1 = (y2 - y3)*(x - x3) + (x3 - x2)*(y - y3);
        double bc2 = (y3 - y1)*(x - x3) + (x1 - x3)*(y - y3);
        bc1 = bc1/det;
        bc2 = bc2/det;
        double bc3 = 1 - bc1 - bc2;
        if (bc1 >= 0 && bc1 <= 1 && bc2 >= 0 && bc2 <= 1 && bc3 >= 0 && bc3 <= 1)
            return true;
        else
            return false;
    }

    class numericRange{
    public:
        numericRange(){};
        bool isInRange(double x);
        void setData(double xmn, double xmx);
        void clear();
    private:
        double xmin;
        double xmax;
    };
    bool numericRange::isInRange(double x) {
        return x > xmin && x <= xmax;
    }
    void numericRange::setData(double xmn, double xmx) {
        xmin = xmn;
        xmax = xmx;
    }
    void numericRange::clear() {
        xmin = 0.0;
        xmax = 0.0;
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
        void multiply(double mult);
        void clear();
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


    /**
	 * @brief Scenario is a struct variable that contains all the information needed for each simulation scenario.
	 *
	 * This struct is a convenient way to transfer around all inputs at once.
	 */
    struct Scenario {
        //! mapID is the id of selected background map.
        std::string mapID;

        //! flowScen is the code name of the selected flow scenario.
        std::string flowScen;

        //! This is a combination of Flow scenario and Welltype
        std::string flowWellScen;

        std::string wellType;

        bool bUseMonitoringWells;
        //bool bUseInitRFsets;
        //std::vector<std::string> RFSets;
        //std::string RFnameSet;

        //! The loading scenario from the list of initialization data
        std::string loadScen;
        int loadScenID;

        // Options for Raster loading type
        //! The name of Scenario Subset
        std::string loadSubScen;
        int loadSubScenID;

        //! The filename with the uploaded raster loading
        std::string modifierName;
        //! The operation type with respect to Scenario subscenario load option
        RasterOperation modifierType;
        //int modReplace;
        LinearData userRasterLoad;
        LoadUnits modifierUnit;

        ///! regionIDs is a list of regions to compute the breakthrough curves
        std::vector<std::string> regionIDs;

        ///! LoadReductionMap is a map that sets the nitrate loading reduction for selected land use categories.
        ///! The key value of the map is the land use category and the value is the percentage of loading reduction.
        std::map<int, double> LoadReductionMap;

        ///! If the input message contains crop id with -9 then this is applied to all crops except to the ones defined
        ///! in the LoadReductionMap
        double globalReduction;

        ///! ReductionYear is the year to start the reduction.
        int startReductionYear;

        ///! The end reduction is the year when the loading has reduced to the desired amount
        int endReductionYear;

        ///! This is a constant reduction coefficient that is applied over the simulation area on all loading types
        float constReduction;
        bool userSuppliedConstRed;

        //! This is the number of years to simulate starting from 1945
        int endSimulationYear;

        int LoadTransitionStart;
        int LoadTransitionEnd;
        std::string LoadTransitionName;
        bool buseLoadTransition;

        //! This is the name of the unsaturated flow scenario
        std::string unsatScenario;
        int unsatScenarioID;

        //! This is the unsaturated zone mobile water content (m3/m3)
        //! Typical values are 0.05,0.1 ,0.15 and 0.20
        double unsatZoneMobileWaterContent;

        /**
         * This is an id for the current scenario. e.g. test001.
         * If this variable is present then during the simulation the following files will be printed
         * - debugID_urfs.dat contains the expanded unit respondse functions
         * - debugID_lf.dat containts the loading functions
         * - debugID_btc.dat contains the breakthrough curves. The convolution between urfs and lf
         * - debugID_well_btc.dat containts the averages BTC for the wells
         */
        std::string debugID;
        //! If the #debugID is present then this is set to true;
        bool printAdditionalInfo = false;

        /**
         * In GNLM the loading is defined as kg/ha. To convert to cocnentration we divide by the groundwater
         * recharge. However sometimes this has very low values due to errors in particle tracing.
         * In those cases letting the NO3 mas divided by a very small number results in unrealistically
         * large concentration. To avoid this we set loading to zero the the recharge is less than minRecharge.
         *
         * By defaults it gets 0.000027 which corresponds to 10 mm/year
         */
        double minRecharge;
        double maxConc;

        bool bUseFlowRch;
        std::string rchName;
        int rchScenID;

        int PixelRadius;

        radialSelection RadSelect;
        rectSelection RectSelect;
        bool useRadSelect;
        bool useRectSelect;

        numericRange DepthRange;
        numericRange ScreenLengthRange;
        bool useDepthRange;
        bool useScreenLenghtRange;
        bool bNarrowSelection;
        bool printWellIds;

        sourceArea SourceArea;


        /**
         * @brief clear is making sure that the scenario has no data from a previous run.
         *
         */
        void clear() {
            mapID = "";
            flowScen = "";
            loadScen = "";
            loadSubScen = "";
            wellType = "";
            flowWellScen = "";
            modifierName = "";
            modifierType = RasterOperation::UNKNOWN;
            userRasterLoad.clear();
            regionIDs.clear();
            startReductionYear = 2020;
            endReductionYear = 2030;
            endSimulationYear = 2100;
            LoadReductionMap.clear();
            printAdditionalInfo = false;
            globalReduction = 1.0;
            constReduction = 1.0;
            userSuppliedConstRed = false;
            unsatZoneMobileWaterContent = 0.0;
            minRecharge = 0.000027; // 10 mm/year
            maxConc = 250;
            PixelRadius = 0;
            RadSelect.clear();
            RectSelect.clear();
            DepthRange.clear();
            ScreenLengthRange.clear();
            useRectSelect = false;
            useRadSelect = false;
            useDepthRange = false;
            useScreenLenghtRange = false;
            bNarrowSelection = false;
            buseLoadTransition = true;
            LoadTransitionName = "GNLM";
            LoadTransitionStart = 2005;
            LoadTransitionEnd = 2015;
            bUseMonitoringWells = false;
            bUseFlowRch = true;
            SourceArea.clear();
        }
    };

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


}

#endif //MANTISSERVER_MSHELPER_H
