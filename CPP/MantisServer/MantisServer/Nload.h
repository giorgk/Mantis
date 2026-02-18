//
// Created by giorgk on 4/28/21.
//

#ifndef MANTISSERVER_NLOAD_H
#define MANTISSERVER_NLOAD_H

#include "Scenario.h"

namespace mantisServer{
    /*!
	 * This is a list of available load types. Currently there are two types available.
	 *  - GNLM loading from GNLM
	 *  - SWAT loading from swat
	 *
	 */
    enum class LoadType {
        GTL,
        RASTER
    };

    /**
	 * Nload is a container for each nitrate loading scenario.
	 *
	 * In future versions this may need restructuring.
	 * It may be better to define a Nload base class and the each type overwrites the methods instead of overloading
	 */
    class NLoad {
    public:
        //! The constructor does nothing
        NLoad() {}
        /**
         * Reads the data from a file. The file must be an hdf5 format with the follwing datasets
         * - LU is a 2D integer array with size 2 Ncells x Nyears, where Ncells is the number of active pixels in the area
         * and Nyears is the number of yearly land use changes. For GNLM Nyears is 5 and SWAT is 1
         * - Nidx is an 1D integer type array of length Ncells. The values correspond the the row of the Nload dataset
         * - NLoad is a 2D double type array that contains the loading time series
         * @param filename is the name of the file
         * @param ltype is the loading type. Currently this must be either GNLM or SWAT
         * @return
         */
        bool readData(std::string filename, LoadType ltype, LoadUnits lunit, extrapMethod xm);

        bool readData(std::string filename, LoadType ltype, LoadUnits lunit, int Npixels);
        /**
         * This returns the N load value. This version must be used for SWAT only
         * @param index is the loading id
         * @param iyr is the year as YYYY.
         * @return
         */
        double getNload(int index,int iyr);
        /**
         * This returns the N load value. This version must be used for GNLM only
         * In GNLM the loading is defined in 15 year intervals.
         * @param index The N loading index
         * @param iyr the year in YYYY format
         * @param N1 is the N load value at the start of the period where the iyr falls in
         * @param N2 is the N load value at the end of the period where the iyr falls in
         * @param u is the parametric value that correspond the iyr. e.g iyr = StartYear *(1-u) + EndYear*u
         */
        void getNload(int index, int iyr, double &N1, double &N2, double &u);

        /**
         * This methods reads the land use information for the GNLM
         * @param index is the CV pixel
         * @param iyr is the year formatted as YYYY
         * @param LUcS The land use type at the start of the 15 year period
         * @param LUcE The land use type at the end of the 15 year period
         * @param perc the parametric value that corresponds to iyr
         */
        void getLU(int index, int iyr, int &LUcS, int &LUcE, double &perc);

        /**
         * This methods reads the land use information for the SWAT type.
         * In SWAT the land use does not change
         * @param index is the CV pixel
         * @param iyr must be always 0
         * @return
         */
        int getLU(int index, int iyr);

        /**
         * This builds the loading function.
         * @param index this is the index of #NLoad
         * @param endYear the ending year of the simulation
         * @param LF THis is the output loading function
         * @param scenario the options of the simulatied scenario
         * @param mult this is the coefficient that converts the loading to concentration for the GNLM case.
         * For SWAT this must be 1
         */
        bool buildLoadingFunction(std::vector<int>& CVindex,
                                  std::vector<double> &rch,
                                  std::vector<double> &clean_prc,
                                  std::vector<double> &unsatDepth,
                                  std::vector<double>& LF,
                                  Scenario& scenario, double& initConc);
        bool buildLoadingFromRaster(std::vector<int>& CVindex,
                                    std::vector<double> &rch,
                                    std::vector<double> &clean_prc,
                                    std::vector<double> &unsatDepth,
                                    std::vector<double>& LF,
                                    Scenario& scenario,  double& initConc);


        bool buildLoadingFromTimeSeries(std::vector<int>& cellIndex,
                                        std::vector<double> &rch,
                                        std::vector<double> &clean_prc,
                                        std::vector<double> &unsatDepth,
                                        std::vector<double>& LF,
                                        Scenario& scenario,  double& initConc);


        LoadType getLtype() {
            return loadType;
        }

        LoadUnits getLunit(){
            return loadUnits;
        }

        int getScenarioID(std::string subScen);
        //! checks if the index is not out of range
        bool isValidIndex(int index);
    private:
        //! The type of loading function GNLM, SWAT or RASTER
        LoadType loadType;
        extrapMethod extMeth;

        LoadUnits loadUnits;
        //! Container for the data
        std::vector<std::vector<double> > Ndata;
        //! Container for the Land use data
        std::vector<std::vector<int> > LU;

        //std::vector<int> LUYears;
        std::vector<int> Nidx;
        LinearData RasterLoading;
        int n_years_lu;
        int n_lu_ids;
        int nN;
        TSgrid LUtsgrid;
        TSgrid Nloadtsgrid;
        bool lu_rows_is_index = false;
        bool nload_rows_is_index = false;
    };

    bool NLoad::isValidIndex(int index) {
        if ((index >= 0) || (index < static_cast<int>(Ndata.size())))
            return true;
        else
            return false;
    }

    int NLoad::getScenarioID(std::string subScen) {
        if (loadType != LoadType::RASTER){
           return -1;
        }
        return RasterLoading.ScenarioIndex(subScen);
    }

    int NLoad::getLU(int index, int iyr) {
        // Single unsigned range-check per argument (also catches negatives).
        if ((unsigned)iyr >= (unsigned)n_years_lu || (unsigned)index >= (unsigned)n_lu_ids)
            return 0;

        // One predictable branch; keep it as the last decision before the load.
        if (lu_rows_is_index)
            return LU[index][iyr];
        else
            return LU[iyr][index];
    }
    void NLoad::getLU(int index, int iyr, int& LUcS, int& LUcE, double& perc) {

        int sy,ey;
        LUtsgrid.getIndices(iyr,sy, ey,perc);
        LUcS = getLU(index,sy);
        LUcE = getLU(index,ey);

    }

    void NLoad::getNload(int index, int iyr, double& N1, double& N2, double& u) {

        int sy, ey;
        Nloadtsgrid.getIndices(iyr,sy,ey,u);
        if (nload_rows_is_index) {
            N1 = Ndata[index][sy];
            N2 = Ndata[index][ey];
        }
        else {
            N1 = Ndata[sy][index];
            N2 = Ndata[ey][index];
        }

    }

    double NLoad::getNload(int index, int iyr) {
        double value = 0.0;
        switch (loadType)
        {
            case mantisServer::LoadType::GTL:// GNLM:
            {
                std::cout << "You can't call this method for GNLM loading type" << std::endl;
                break;
            }
            case mantisServer::LoadType::RASTER:// SWAT:
            {
                int load_index = (iyr - 1940) % 25;
                //value = Ndata[index][load_index];
                value = Ndata[load_index][index];
                break;
            }
            default:
                break;
        }
        return value;
    }

    bool NLoad::readData(std::string filename, LoadType ltype, LoadUnits lunit, int Npixels) {
        loadType = ltype;
        auto start = std::chrono::high_resolution_clock::now();
        if (ltype == LoadType::RASTER){
            RasterLoading.setNoDataValue(0.0);
            RasterLoading.readData(filename, Npixels);
            return true;
        }
        else{
            std::cout << "This method is used only for RASTER loading type" << std::endl;
            return false;
        }
    }

    bool NLoad::readData(std::string filename, LoadType ltype, LoadUnits lunit, extrapMethod xm) {
        loadType = ltype;
        loadUnits = lunit;
        auto start = std::chrono::high_resolution_clock::now();
#if _USEHF>0
        std::string ext = getExtension(filename);

        if (ext.compare("h5") == 0){
            const std::string LUNameSet("LU");
            const std::string NidxNameSet("Nidx");
            const std::string NloadNameSet("Nload");
            const std::string LUNGNameSet("LUNgrid");
            HighFive::File HDFNfile(filename, HighFive::File::ReadOnly);
            HighFive::DataSet datasetLU = HDFNfile.getDataSet(LUNameSet);
            HighFive::DataSet datasetNidx = HDFNfile.getDataSet(NidxNameSet);
            HighFive::DataSet datasetNLoad = HDFNfile.getDataSet(NloadNameSet);
            HighFive::DataSet datasetLUNG = HDFNfile.getDataSet(LUNGNameSet);
            datasetLU.read(LU);
            datasetNidx.read(Nidx);
            datasetNLoad.read(Ndata);
            std::vector<int> lung;
            datasetLUNG.read(lung);
            if (lung.size() != 6){
                std::cout << "The LUNgrid property must have 6 values" << std::endl;
                return false;
            }
            LUtsgrid.init(lung[0],lung[2],lung[1],xm);
            Nloadtsgrid.init(lung[3],lung[5],lung[4],xm);
            if (LUtsgrid.getNyears() == LU.size()) {
                lu_rows_is_index = false;
            }
            else if (LUtsgrid.getNyears() == LU[0].size()) {
                lu_rows_is_index = true;
            }
            else {
                std::cout << "The Land Use years grid property must have the same number of years as LU" << std::endl;
                return false;
            }
            if (Nloadtsgrid.getNyears() == Ndata.size()) {
                nload_rows_is_index = false;
            }
            else if (Nloadtsgrid.getNyears() == Ndata[0].size()) {
                nload_rows_is_index = true;
            }
            else {
                std::cout << "The Nload years grid property must have the same number of years as Nload" << std::endl;
                return false;
            }

            n_years_lu = LUtsgrid.getNyears();

            nN = Ndata.size();
            n_lu_ids = LU.size();
            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;
            std::cout << "Read Nload from " << filename << " in " << elapsed.count() << std::endl;
            return true;
        }
#endif


        {// First read the indices
            std::string idxlufile = filename + ".idxlu";
            std::ifstream ifile;
            ifile.open(idxlufile);
            if (!ifile.is_open()){
                std::cout << "Cant open file: " << idxlufile << std::endl;
                return false;
            }
            else{
                std::cout << "Reading " << idxlufile << std::endl;
                std::string line;
                int Nr, Sy, D, Ny; // Number of rows, Starting year, Interval, Number of years
                {// Get the dimension of the data
                    getline(ifile, line);
                    std::istringstream inp(line.c_str());
                    inp >> Nr;
                    inp >> Sy;
                    inp >> D;
                    inp >> Ny;
                    LUtsgrid.init(Sy,Ny,D,xm);
                    Nidx.clear();
                    Nidx.resize(Nr,0);
                    LU.clear();
                    std::vector<int> tmp(Nr,0);
                    for (int i = 0; i < Ny; ++i){
                        LU.push_back(tmp);
                    }
                }
                {// Read the data
                    for (int i = 0; i < Nr; ++i){
                        getline(ifile, line);
                        std::istringstream inp(line.c_str());
                        int v;
                        for (int j = 0; j < Ny+1; ++j){
                            inp >> v;
                            if (j == 0){
                                Nidx[i] = v;
                            }
                            else{
                                LU[j-1][i] = v;
                            }
                        }
                    }
                }
                ifile.close();
            }
        }

        {// Read the N load data
            std::string nloadfile = filename + ".nload";
            std::ifstream ifile;
            ifile.open(nloadfile);
            if (!ifile.is_open()){
                std::cout << "Cant open file: " << nloadfile << std::endl;
                return false;
            }
            else{
                std::cout << "Reading " << nloadfile << std::endl;
                std::string line;
                int Nr, Sy, D, Ny;
                {// Get the dimension of the data
                    getline(ifile, line);
                    std::istringstream inp(line.c_str());
                    inp >> Nr;
                    inp >> Sy;
                    inp >> D;
                    inp >> Ny;
                    Nloadtsgrid.init(Sy,Ny,D,xm);
                    Ndata.clear();
                    std::vector<double> tmp(Nr,0);
                    for (int i = 0; i < Ny; ++i){
                        Ndata.push_back(tmp);
                    }
                }
                {// Read the data
                    for (int i = 0; i < Nr; ++i){
                        getline(ifile, line);
                        std::istringstream inp(line.c_str());
                        double v;
                        for (int j = 0; j < Ny; ++j){
                            inp >> v;
                            Ndata[j][i] = v;
                        }
                    }
                }
                ifile.close();
            }
            n_years_lu = LU[0].size();
            nN = Ndata.size();
            n_lu_ids = LU.size();
        }

        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "Read Nload from " << filename << " in " << elapsed.count() << std::endl;
        return true;
    }

    bool NLoad::buildLoadingFromRaster(std::vector<int>& cellIndex,
                                       std::vector<double> &rch,
                                       std::vector<double> &clean_prc,
                                       std::vector<double> &unsatDepth,
                                       std::vector<double>& LF,
                                       Scenario& scenario,  double& initConc) {

        bool out = false;
        int Nyears = scenario.endSimulationYear - scenario.startSimulationYear;
        LF.clear();
        LF.resize(Nyears, 0.0);

        if (cellIndex.size() == 0){
            out = true;
            return out;
        }

        //double load_value = 0.0f;
        std::vector<double> load_vals;
        std::vector<double> tau_vals;
        for (unsigned int j = 0; j < cellIndex.size(); ++j){
            // Calculate the Unsaturated travel time
            double tau = scenario.unsatZoneMobileWaterContent*
                         ( (scenario.unsatDepthOffset + scenario.unsatDepthMult * unsatDepth[j]) /
                           (scenario.unsatRchOffset + scenario.unsatRchMult * rch[j]/1000.0) );
            tau_vals.push_back(static_cast<int>(std::ceil(tau)));

            // Find the loading value for each source area cell
            if (scenario.modifierType == RasterOperation::DONTUSE){// Use the scenario from the Raster Initialization Data
                double scen_value = RasterLoading.getValue(scenario.loadSubScenID, cellIndex[j]);
                if (loadUnits == LoadUnits::MASS){
                    scen_value = scen_value*100 / rch[j];
                }
                load_vals.push_back(scen_value);
                //load_value += scen_value;
            }
            else if (scenario.modifierType == RasterOperation::Replace){// Replace the value with the user data
                double user_value = scenario.userRasterLoad.getValue(0, cellIndex[j]);
                if (scenario.modifierUnit == LoadUnits::MASS){
                    user_value = user_value*100 / rch[j];
                }
                load_vals.push_back(user_value);
                //load_value += user_value;
            }
            else if (scenario.modifierType == RasterOperation::Multiply){// Multiply the user value with the initialization scenario Raster
                double user_value = scenario.userRasterLoad.getValue(0, cellIndex[j]);
                double scen_value = RasterLoading.getValue(scenario.loadSubScenID, cellIndex[j]);
                if (loadUnits == LoadUnits::MASS){
                    scen_value = scen_value*100 / rch[j];
                }
                load_vals.push_back(scen_value*user_value);
                //load_value += scen_value*user_value;
            }
        }

        int istartReduction = scenario.startReductionYear - scenario.startSimulationYear;
        int iendReduction = scenario.endReductionYear - scenario.startSimulationYear;
        double dstartReduction = static_cast<double>(istartReduction);
        double dReductionRange = static_cast<double>(iendReduction) - dstartReduction;
        double adoptionCoeff = 0.0;

        for (int iyr = 0; iyr < Nyears; ++iyr){
            if (scenario.userSuppliedConstRed){
                if ((iyr >= istartReduction) && (iyr <= iendReduction)) {
                    adoptionCoeff = (static_cast<double>(iyr) - dstartReduction) / dReductionRange;
                }
                else if (iyr > iendReduction){
                    adoptionCoeff = 1.0;
                }
                else{
                    adoptionCoeff = 0.0;
                }
            }
            for (int i = 0; i < load_vals.size(); ++i){
                if (iyr < tau_vals[i]){
                    LF[iyr] = LF[iyr] + initConc;
                }
                else{
                    LF[iyr] = LF[iyr] + load_vals[i]*(1-adoptionCoeff) + (load_vals[i]*scenario.constReduction)*adoptionCoeff;
                }
            }
        }
        double nLoadValues = static_cast<double>(cellIndex.size());
        for (int iyr = 0; iyr < Nyears; ++iyr){
            LF[iyr] = LF[iyr]/nLoadValues;
            if (scenario.maxConc > 0){
                if (LF[iyr] > scenario.maxConc){
                    LF[iyr] = scenario.maxConc;
                }
            }
        }
        return true;
    }

    bool NLoad::buildLoadingFromTimeSeries(std::vector<int>& cellIndex,
                                           std::vector<double> &rch,
                                           std::vector<double> &clean_prc,
                                           std::vector<double> &unsatDepth,
                                           std::vector<double>& LF,
                                           Scenario& scenario, double& initConc){
        bool out = false;

        int Nyears = scenario.endSimulationYear - scenario.startSimulationYear;
        LF.clear();
        LF.resize(Nyears, 0.0);
        if (cellIndex.size() == 0){
            out = true;
            return out;
        }

        int istartReduction = scenario.startReductionYear - scenario.startSimulationYear;
        int iendReduction = scenario.endReductionYear - scenario.startSimulationYear;
        double dstartReduction = static_cast<double>(istartReduction);
        double dReductionRange = static_cast<double>(iendReduction) - dstartReduction;
        double adoptionCoeff = 0;
        std::map<int, double>::iterator it;
        int lfidx;

        bool bIstauCalc = false;
        std::vector<int> tauVal;

        for (int iyr = 0; iyr < Nyears; ++iyr){
            double lf = 0;
            if ((iyr >= istartReduction) && (iyr <= iendReduction)) {
                adoptionCoeff = (static_cast<double>(iyr) - dstartReduction) / dReductionRange;
            }
            else if (iyr > iendReduction) {
                adoptionCoeff = 1.0;
            }
            else{
                adoptionCoeff = 0;
            }

            for (unsigned int j = 0; j < cellIndex.size(); ++j){
                int tauShift = 0;
                if (!bIstauCalc){
                    double tau = scenario.unsatZoneMobileWaterContent*
                                 ( (scenario.unsatDepthOffset + scenario.unsatDepthMult * unsatDepth[j]) /
                                   (scenario.unsatRchOffset + scenario.unsatRchMult * rch[j]/1000.0) );

                    tauShift = static_cast<int>(std::ceil(tau));
                    tauVal.push_back(tauShift);
                    if (j == cellIndex.size()-1){
                        bIstauCalc = true;
                    }
                    for (int itau = 0; itau < tauShift; ++itau){
                        if (itau < LF.size()){
                            LF[itau] =  LF[itau] + initConc;
                        }
                    }
                }
                else{
                    tauShift = tauVal[j];
                }

                lfidx = iyr + tauShift;
                if (lfidx >= Nyears){
                    continue;
                }
                double rs = 1.0;
                double re = 1.0;
                if (adoptionCoeff > 0){
                    int lus = 0;
                    int lue = 0;
                    double prc = 0.0;
                    getLU(cellIndex[j], iyr + scenario.startSimulationYear, lus, lue, prc);
                    rs = scenario.globalReduction;
                    it = scenario.LoadReductionMap.find(lus);
                    if (it != scenario.LoadReductionMap.end()){
                        rs = it->second;
                    }

                    if (lus == lue){
                        re = rs;
                    }
                    else{
                        re = scenario.globalReduction;
                        it = scenario.LoadReductionMap.find(lue);
                        if (it != scenario.LoadReductionMap.end()) {
                            re = it->second;
                        }
                    }
                }

                int nload_idx = Nidx[cellIndex[j]];
                if (nload_idx < 0)
                    continue;

                double N1 = 0, N2 = 0, u = 0;
                getNload(nload_idx, iyr + scenario.startSimulationYear, N1, N2, u);
                double Nbase = N1* (1 - u) + N2* u;
                if ((adoptionCoeff > 0) && ((std::abs(1 - rs) > 0.000000001) || (std::abs(1 - re) > 0.000000001))){
                    double Nred = (N1 * rs) * (1 - u) + (N2 * re) * u;
                    double tmpLoad = (Nbase * (1 - adoptionCoeff) + Nred * adoptionCoeff);
                    if (loadUnits == LoadUnits::MASS){
                        tmpLoad = tmpLoad*100 / rch[j];
                        lf = tmpLoad;
                    }
                    else{
                        lf = tmpLoad*(1 - clean_prc[j]);
                    }
                }
                else{
                    if (loadUnits == LoadUnits::MASS){
                        Nbase = Nbase*100 / rch[j];
                        lf = Nbase;
                    }
                    else{
                        lf = Nbase*(1 - clean_prc[j]);
                    }
                }
                LF[lfidx] = LF[lfidx] + lf;
            }
        }
        double nLoadValues = static_cast<double>(cellIndex.size());
        for (int iyr = 0; iyr < Nyears; ++iyr){
            LF[iyr] = LF[iyr]/nLoadValues;
            if (scenario.maxConc > 0){
                if (LF[iyr] > scenario.maxConc){
                    LF[iyr] = scenario.maxConc;
                }
            }
        }
        return true;
    }

    bool NLoad::buildLoadingFunction(std::vector<int>& CVindex,
                                     std::vector<double> &rch,
                                     std::vector<double> &clean_prc,
                                     std::vector<double> &unsatDepth,
                                     std::vector<double>& LF,
                                     Scenario& scenario,  double& initConc) {
        bool out = false;

        if (loadType == LoadType::RASTER){
            out = buildLoadingFromRaster(CVindex, rch, clean_prc, unsatDepth, LF, scenario, initConc);
            return out;
        }
        else{
            out = buildLoadingFromTimeSeries(CVindex, rch, clean_prc, unsatDepth, LF, scenario, initConc);
            return out;
        }
    }

    class NLoadList{
    public:
        NLoadList(){}

        std::map<std::string, NLoad> NLoadMaps;
        bool readData(std::string path, std::string filename, int ncells);
        bool hasLoading(std::string loadName);
    };

    bool NLoadList::hasLoading(std::string loadName) {
        std::map<std::string, NLoad>::iterator it;
        it = NLoadMaps.find(loadName);
        if (it == NLoadMaps.end()){
            return false;
        }
        return true;
    }

    bool NLoadList::readData(std::string path, std::string filename, int ncells){
        auto start = std::chrono::high_resolution_clock::now();
        filename = path + filename;
        std::ifstream no3MainFile;
        no3MainFile.open(filename);
        if (!no3MainFile.is_open()) {
            std::cout << "Cant open file: " << filename << std::endl;
            return false;
        }
        else{
            std::string line;
            while (getline(no3MainFile, line)){
                std::string Ltype, Lname, Lfile, Lunit, Lextrap;
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
                inp >> Lextrap;
                extrapMethod xm = string2xtrapMethod(Lextrap);
                Lfile = path + Lfile;

                NLoad NL;
                NLoadMaps.insert(std::pair<std::string, NLoad>(Lname,NL));
                bool tf;
                if (Ltype.compare("GTL") == 0){
                    tf = NLoadMaps[Lname].readData(Lfile, LoadType::GTL , loadunit, xm);
                }
                else if (Ltype.compare("RASTER") == 0){
                    tf = NLoadMaps[Lname].readData(Lfile, LoadType::RASTER, loadunit, ncells);
                }
                else{
                    std::cout << "Unknown loading type " << Ltype << std::endl;
                    return false;
                }
                if (!tf){
                    return false;
                }
            }
        }



        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "Read Loading maps in " << elapsed.count() << std::endl;
        return true;
    }
}

#endif //MANTISSERVER_NLOAD_H
