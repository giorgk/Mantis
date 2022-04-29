//
// Created by giorgk on 4/28/21.
//

#ifndef MANTISSERVER_NLOAD_H
#define MANTISSERVER_NLOAD_H

namespace mantisServer{
    /*!
	 * This is a list of available load types. Currently there are two types available.
	 *  - GNLM loading from GNLM
	 *  - SWAT loading from swat
	 *
	 */
    enum class LoadType {
        GNLM,
        SWAT,
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
        bool readData(std::string filename, LoadType ltype, LoadUnits lunit);

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
         * @param endYear the ending year of hte simulation
         * @param LF THis is the output loading function
         * @param scenario the options of the simulatied scenario
         * @param mult this is the coefficient that converts the loading to concentration for the GNLM case.
         * For SWAT this must be 1
         */
        bool buildLoadingFunction(std::vector<int>& CVindex, int startYear, int endYear, std::vector<double>& LF,
                                  Scenario& scenario, std::vector<double> &rch);
        bool buildLoadingFromRaster(std::vector<int>& CVindex, int startYear, int endYear, std::vector<double>& LF,
                                    Scenario& scenario, std::vector<double> &rch);
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

        LoadUnits loadUnits;
        //! Container for the data
        std::vector<std::vector<double> > Ndata;
        //! Container for the Land use data
        std::vector<std::vector<int> > LU;
        //! A map between the CV active pixels and the Ndata
        std::vector<int> Nidx;
        LinearData RasterLoading;
        int nRows;
        int nLU;
        int nN;
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
        if (iyr >=0 && iyr < nLU && index >=0 && index < nRows)
            return LU[iyr][index];
        return 0;
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
                    N1 = Ndata[0][index];
                    N2 = Ndata[0][index];
                    u = 1.0;
                }
                else if (iyr >= 2050) {
                    N1 = Ndata[7][index];
                    N2 = Ndata[7][index];
                    u = 0.0;
                }
                else {

                    int isN = static_cast<int>(std::floor((iyr - 1945) / 15)); //index of starting year
                    int ieN = isN + 1; // index of starting year
                    //std::cout << " isN=" << isN << " ieN=" << ieN;
                    u = static_cast<double>((iyr - 1945) % 15) / 15.f;
                    N1 = Ndata[isN][index];
                    N2 = Ndata[ieN][index];
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

    double NLoad::getNload(int index, int iyr) {
        double value = 0.0;
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

    bool NLoad::readData(std::string filename, LoadType ltype, LoadUnits lunit) {
        loadType = ltype;
        loadUnits = lunit;
        auto start = std::chrono::high_resolution_clock::now();
#if _USEHF>0
        std::string ext = getExtension(filename);

        if (ext.compare("h5") == 0){
            const std::string LUNameSet("LU");
            const std::string NidxNameSet("Nidx");
            const std::string NloadNameSet("NLoad");
            HighFive::File HDFNfile(filename, HighFive::File::ReadOnly);
            HighFive::DataSet datasetLU = HDFNfile.getDataSet(LUNameSet);
            HighFive::DataSet datasetNidx = HDFNfile.getDataSet(NidxNameSet);
            HighFive::DataSet datasetNLoad = HDFNfile.getDataSet(NloadNameSet);
            datasetLU.read(LU);
            datasetNidx.read(Nidx);
            datasetNLoad.read(Ndata);
            nRows = LU[0].size();
            nN = Ndata.size();
            nLU = LU.size();
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
                int Nr, Nc;
                {// Get the dimension of the data
                    getline(ifile, line);
                    std::istringstream inp(line.c_str());
                    inp >> Nr;
                    inp >> Nc;
                    Nidx.clear();
                    Nidx.resize(Nr,0);
                    LU.clear();
                    std::vector<int> tmp(Nr,0);
                    for (int i = 0; i < Nc-1; ++i){
                        LU.push_back(tmp);
                    }
                }
                {// Read the data
                    for (int i = 0; i < Nr; ++i){
                        getline(ifile, line);
                        std::istringstream inp(line.c_str());
                        int v;
                        for (int j = 0; j < Nc; ++j){
                            inp >> v;
                            if (j == 0){
                                Nidx[i] = v;
                            }
                            else{
                                LU[i][j-1] = v;
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
                int Nr, Nc;
                {// Get the dimension of the data
                    getline(ifile, line);
                    std::istringstream inp(line.c_str());
                    inp >> Nr;
                    inp >> Nc;
                    Ndata.clear();
                    Ndata.resize(Nr, std::vector<double>(Nc, 0));
                }
                {// Read the data
                    for (int i = 0; i < Nr; ++i){
                        getline(ifile, line);
                        std::istringstream inp(line.c_str());
                        double v;
                        for (int j = 0; j < Nc; ++j){
                            inp >> v;
                            Ndata[i][j] = v;
                        }
                    }
                }
                ifile.close();
            }
        }
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "Read Nload from " << filename << " in " << elapsed.count() << std::endl;
        return true;
        /*


        hf::File HDFNfile(filename, hf::File::ReadOnly);
        hf::DataSet datasetLU = HDFNfile.getDataSet(LUNameSet);
        hf::DataSet datasetNidx = HDFNfile.getDataSet(NidxNameSet);
        hf::DataSet datasetNLoad = HDFNfile.getDataSet(NloadNameSet);
        datasetLU.read(LU);
        datasetNidx.read(Nidx);
        datasetNLoad.read(Ndata);
        loadType = ltype;

        return true;
        */
    }

    bool NLoad::buildLoadingFromRaster(std::vector<int> &CVindex, int startYear, int endYear, std::vector<double> &LF,
                                       Scenario &scenario, std::vector<double> &rch) {

        double load_value = 0.0f;
        for (unsigned int j = 0; j < CVindex.size(); ++j){
            if (scenario.modifierType == RasterOperation::DONTUSE){// Use the scenario from the Raster Initialization Data
                double scen_value = RasterLoading.getValue(scenario.loadSubScenID, CVindex[j]);
                if (loadUnits == LoadUnits::MASS){
                    scen_value = scen_value*100 / rch[j];
                }
                load_value += scen_value;
            }
            else if (scenario.modifierType == RasterOperation::Replace){// Replace the value with the user data
                double user_value = scenario.userRasterLoad.getValue(0, CVindex[j]);
                if (scenario.modifierUnit == LoadUnits::MASS){
                    user_value = user_value*100 / rch[j];
                }
                load_value += user_value;
            }
            else if (scenario.modifierType == RasterOperation::Multiply){// Multiply the user value with the initialization scenario Raster
                double user_value = scenario.userRasterLoad.getValue(0, CVindex[j]);
                double scen_value = RasterLoading.getValue(scenario.loadSubScenID, CVindex[j]);
                if (loadUnits == LoadUnits::MASS){
                    scen_value = scen_value*100 / rch[j];
                }
                load_value += scen_value*user_value;
            }
        }
        load_value = load_value/static_cast<double>(CVindex.size());

        if (scenario.maxConc > 0){
            if (load_value > scenario.maxConc){
                load_value = scenario.maxConc;
            }
        }

        int Nyears = endYear - startYear;
        LF.clear();
        LF.resize(Nyears, load_value);
        if (scenario.userSuppliedConstRed){
            int istartReduction = scenario.startReductionYear - startYear;
            int iendReduction = scenario.endReductionYear - startYear;
            double dstartReduction = static_cast<double>(istartReduction);
            double dReductionRange = static_cast<double>(iendReduction) - dstartReduction;
            double adoptionCoeff = 0;
            for (int iyr = 0; iyr < Nyears; iyr++){
                if ((iyr >= istartReduction) && (iyr <= iendReduction))
                    adoptionCoeff = (static_cast<double>(iyr) - dstartReduction) / dReductionRange;
                else if (iyr > iendReduction)
                    adoptionCoeff = 1.0;
                LF[iyr] = LF[iyr] * (1 - adoptionCoeff) + LF[iyr]*scenario.constReduction * adoptionCoeff;
            }
        }
        return true;
    }

    bool NLoad::buildLoadingFunction(std::vector<int>& CVindex,
                                     int startYear,
                                     int endYear,
                                     std::vector<double>& LF,
                                     Scenario& scenario,
                                     std::vector<double> &rch) {
        bool out = false;

        if (loadType == LoadType::RASTER){
            out = buildLoadingFromRaster(CVindex, startYear, endYear, LF, scenario, rch);
            return out;
        }


        //int startYear = 1945;
        int istartReduction = scenario.startReductionYear - startYear;
        int iendReduction = scenario.endReductionYear - startYear;
        double dstartReduction = static_cast<double>(istartReduction);
        double dReductionRange = static_cast<double>(iendReduction) - dstartReduction;
        int Nyears = endYear - startYear;
        double adoptionCoeff = 0;
        LF.resize(Nyears, 0.0);

        std::vector<double> percReduction(CVindex.size(), scenario.globalReduction);
        std::map<int, double>::iterator it;
        if (loadType == LoadType::SWAT) {
            for (unsigned int i = 0; i < CVindex.size(); ++i){
                int lucode = getLU(CVindex[i], 0);
                it = scenario.LoadReductionMap.find(lucode);
                if (it != scenario.LoadReductionMap.end()) {
                    percReduction[i] = it->second;
                }
            }
        }

        for (int iyr = 0; iyr < Nyears; iyr++) {
            //std::cout << "i=" << i << " Y=" << i + startYear << std::endl;
            if ((iyr >= istartReduction) && (iyr <= iendReduction))
                adoptionCoeff = (static_cast<double>(iyr) - dstartReduction) / dReductionRange;
            else if (iyr > iendReduction)
                adoptionCoeff = 1.0;

            //std::cout << "a=" << adoptionCoeff;

            if (loadType == LoadType::GNLM) {
                double lf = 0;
                if (CVindex.size() == 0){
                    LF[iyr] = 0;
                }
                else{
                    double NvalidCells = 0;
                    for (unsigned int j = 0; j < CVindex.size(); ++j){
                        int lus = 0;
                        int lue = 0;
                        double prc = 0.0;
                        double rs = 1.0;
                        double re = 1.0;
                        getLU(CVindex[j], iyr + startYear, lus, lue, prc);
                        if (adoptionCoeff > 0) {
                            rs = scenario.globalReduction;
                            it = scenario.LoadReductionMap.find(lus);
                            if (it != scenario.LoadReductionMap.end()) {
                                rs = it->second;
                                //std::cout << " rs=" << rs;
                            }
                            re = scenario.globalReduction;
                            it = scenario.LoadReductionMap.find(lue);
                            if (it != scenario.LoadReductionMap.end()) {
                                re = it->second;
                                //std::cout << " rs=" << rs;
                            }
                        }
                        double N1 = 0, N2 = 0, u = 0;
                        int nload_idx = Nidx[CVindex[j]];
                        if (nload_idx < 0)
                            continue;
                        NvalidCells = NvalidCells + 1.0;
                        getNload(nload_idx, iyr + startYear, N1, N2, u);
                        double Nbase = N1* (1 - u) + N2* u;
                        if ((adoptionCoeff > 0) && ((std::abs(1 - rs) > 0.000000001) || (std::abs(1 - re) > 0.000000001))) {
                            double Nred = (N1 * rs) * (1 - u) + (N2 * re) * u;
                            //std::cout << " Nred=" << Nred;
                            double tmpLoad = (Nbase * (1 - adoptionCoeff) + Nred * adoptionCoeff);
                            if (loadUnits == LoadUnits::MASS){
                                tmpLoad = tmpLoad*100 / rch[j];
                            }
                            lf += tmpLoad;
                        }
                        else{
                            if (loadUnits == LoadUnits::MASS){
                                Nbase = Nbase*100 / rch[j];
                            }
                            lf += Nbase;
                        }
                    }

                    if (NvalidCells < 0.00001){
                        out = false;
                        return out;
                    }
                    else if (std::abs(NvalidCells - 1) < 0.000001){
                        LF[iyr] = lf;
                    }
                    else if (NvalidCells > 1){
                        LF[iyr] = lf/NvalidCells;
                    }
                }
            }
            else if (loadType == LoadType::SWAT) {
                double lf = 0;
                double NvalidCells = 0;
                for (unsigned int j = 0; j < CVindex.size(); ++j){
                    int nload_idx = Nidx[CVindex[j]];
                    if (nload_idx < 0)
                        continue;
                    NvalidCells = NvalidCells + 1.0;
                    double Nbase = getNload(nload_idx, iyr + startYear);
                    double Nred = percReduction[j] * Nbase;
                    double tmpLoad = (Nbase * (1 - adoptionCoeff) + Nred * adoptionCoeff);
                    if (loadUnits == LoadUnits::MASS){
                        tmpLoad = tmpLoad*100 / rch[j];
                    }
                    lf += tmpLoad;
                    //std::cout << Nbase << ", " << Nred << std::endl;
                }

                if (NvalidCells < 0.000001){
                    out = false;
                    return out;
                }
                else if (std::abs(NvalidCells - 1) < 0.000001){
                    LF[iyr] = lf;
                }
                else if (NvalidCells > 1){
                    LF[iyr] = lf/NvalidCells;
                }
            }

            if (scenario.maxConc > 0){
                if (LF[iyr] > scenario.maxConc){
                    LF[iyr] = scenario.maxConc;
                }
            }
        }
        return true;
    }
}

#endif //MANTISSERVER_NLOAD_H
