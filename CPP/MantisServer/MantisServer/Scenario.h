//
// Created by giorgk on 5/8/2023.
//

#ifndef MANTISSERVER_SCENARIO_H
#define MANTISSERVER_SCENARIO_H

#include "MShelper.h"

namespace mantisServer {

/**
	 * @brief Scenario is a struct variable that contains all the information needed for each simulation scenario.
	 *
	 * This struct is a convenient way to transfer around all inputs at once.
	 */
    class Scenario {
    public:
        Scenario() {};

        void clear();

        bool parse_incoming_msg(std::string &msg, std::string &outmsg);


        std::string modelArea;

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
        int startSimulationYear;
        int endSimulationYear;

        int LoadTransitionStart;
        int LoadTransitionEnd;
        std::string LoadTransitionName;
        bool buseLoadTransition;

        //! This is the name of the unsaturated flow scenario
        std::string unsatScenario;
        int unsatScenarioID;

        int porosity = 10;
        int porosityIndex = 0;

        //! This is the unsaturated zone mobile water content (m3/m3)
        //! Typical values are 0.05,0.1 ,0.15 and 0.20
        double unsatZoneMobileWaterContent;

        double unsatMinDepth = 1.0;
        double unsatDepthOffset = 0.0;
        double unsatDepthMult = 1.0;
        double unsatMinRch = 10.0;
        double unsatRchOffset = 0.0;
        double unsatRchMult = 1;

        /**
         * This is an id for the current scenario. e.g. test001.
         * If this variable is present then during the simulation the following files will be printed
         * - debugID_urfs.dat contains the expanded unit respondse functions
         * - debugID_lf.dat containts the loading functions
         * - debugID_btc.dat contains the breakthrough curves. The convolution between urfs and lf
         * - debugID_well_btc.dat containts the averages BTC for the wells
         */
        std::string debugID;
        std::string debugPath;
        bool printLF = false;
        bool printBTC = false;
        bool printURF = false;
        bool printWellBTC = false;
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
        int maxSourceCells;
        URFTYPE urfType = URFTYPE::LGNRM;
        double adeLambda = 0.0;
        double adeR = 1.0;
    };

/**
     * @brief clear is making sure that the scenario has no data from a previous run.
     *
     */
    void Scenario::clear() {
        modelArea = "";
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
        startSimulationYear = 1945;
        endSimulationYear = 2100;
        LoadReductionMap.clear();
        globalReduction = 1.0;
        constReduction = 1.0;
        userSuppliedConstRed = false;
        unsatZoneMobileWaterContent = 0.0;
        minRecharge = 10;// mm/year
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
        porosity = 10;
        porosityIndex = 0;
        SourceArea.clear();
        printWellIds = false;
        debugID = "";
        printLF = false;
        printURF = false;
        printBTC = false;
        printWellBTC = false;
        printAdditionalInfo = false;
        maxSourceCells = 1000;
        urfType = URFTYPE::LGNRM;
        adeLambda = 0.0;
        adeR = 1.0;
    }

    bool Scenario::parse_incoming_msg(std::string &msg, std::string &outmsg){
        bool out = false;
        clear();
        // convert the string to a stringstream
        std::istringstream ss(msg);

        int counter = 0;
        while (true){
            std::string test, tmp;
            counter = counter + 1;
            ss >> test;

            if (test == "modelArea") {
                ss >> modelArea;
                continue;
            }
            if (test == "bMap") {
                ss >> mapID;
                continue;
            }
            if (test == "Nregions") {
                int Nregions;
                ss >> Nregions;
                for (int i = 0; i < Nregions; ++i) {
                    ss >> test;
                    regionIDs.push_back(test);
                }
                continue;
            }
            if (test == "RadSelect"){
                bNarrowSelection = true;
                useRadSelect = true;
                double cx, cy, r;
                ss >> cx;
                ss >> cy;
                ss >> r;
                RadSelect.setData(cx,cy,r);
                continue;
            }

            if (test == "RectSelect"){
                bNarrowSelection = true;
                useRectSelect = true;
                double xmin, ymin, xmax, ymax;
                ss >> xmin;
                ss >> ymin;
                ss >> xmax;
                ss >> ymax;
                RectSelect.setData(xmin, ymin, xmax, ymax);
                continue;
            }

            if (test == "DepthRange"){
                bNarrowSelection = true;
                useDepthRange = true;
                double dmin, dmax;
                ss >> dmin;
                ss >> dmax;
                DepthRange.setData(dmin, dmax);
                continue;
            }

            if (test == "ScreenLenRange"){
                bNarrowSelection = true;
                useScreenLenghtRange = true;
                double slmin, slmax;
                ss >> slmin;
                ss >> slmax;
                ScreenLengthRange.setData(slmin, slmax);
                continue;
            }
            if (test == "flowScen") {
                ss >> flowScen;
                continue;
            }
            if (test == "wellType"){
                ss >> wellType;
                continue;
            }
            if (test == "unsatScen") {
                ss >> unsatScenario;
                continue;
            }
            if (test == "unsatWC") {
                ss >> unsatZoneMobileWaterContent;
                continue;
            }
            if (test == "rchMap"){
                bUseFlowRch = false;
                ss >> rchName;
                continue;
            }
            if (test == "minRch") {
                ss >> minRecharge;
                continue;
            }
            if (test == "startSimYear") {
                ss >> startSimulationYear;
                continue;
            }
            if (test == "endSimYear") {
                ss >> endSimulationYear;
                continue;
            }
            if (test == "startRed" ) {
                ss >> startReductionYear;
                continue;
            }
            if (test == "endRed") {
                ss >> endReductionYear;
                continue;
            }
            if (test == "loadScen") {
                ss >> loadScen;
                continue;
            }
            if (test == "loadSubScen") {
                ss >> loadSubScen;
                continue;
            }
            if (test == "modifierName") {
                ss >> modifierName;
                continue;
            }
            if (test == "modifierType") {
                ss >> tmp;
                modifierType = string2RasterOperation(tmp);
                continue;
            }
            if (test == "modifierUnit") {
                ss >> tmp;
                modifierUnit = string2LoadUnits(tmp);
                continue;
            }
            if (test == "Ncrops") {
                int Ncrops, cropid;
                double perc;
                ss >> Ncrops;
                for (int i = 0; i < Ncrops; ++i) {
                    ss >> cropid;
                    ss >> perc;
                    if (cropid == -9){
                        globalReduction = perc;
                    }
                    else{
                        LoadReductionMap.insert(std::pair<int, double>(cropid, perc));
                    }
                }
                continue;
            }
            if (test == "maxConc") {
                ss >> maxConc;
                continue;
            }
            if (test == "constRed"){
                userSuppliedConstRed = true;
                ss >> constReduction;
                continue;
            }
            if (test == "loadTrans"){
                ss >> LoadTransitionName;
                ss >> LoadTransitionStart;
                ss >> LoadTransitionEnd;
                if (LoadTransitionName.compare("NONE") == 0){
                    buseLoadTransition = false;
                }
                else{
                    buseLoadTransition = true;
                }

                continue;
            }

            if (test == "por"){
                ss >> porosity;
                continue;
            }

            if (test == "urfType"){
                std::string tmp;
                ss >> tmp;
                urfType = string2URFTYPE(tmp);
                continue;
            }

            if (test == "ADELR"){
                ss >> adeLambda;
                ss >> adeR;
                continue;
            }

            if (test == "maxSourceCells") {
                ss >> maxSourceCells;
                continue;
            }

            if (test == "SourceArea"){//TODO
                int nPix, minPix, maxPix;
                double percPix;
                ss >> nPix;
                ss >> minPix;
                ss >> maxPix;
                ss >> percPix;
                SourceArea.setParameters(nPix, minPix, maxPix, percPix);
                continue;
            }
            if (test ==  "PixelRadius") {//TODO
                ss >> PixelRadius;
                continue;
            }
            if (test == "getids"){
                int tmp;
                ss >> tmp;
                printWellIds = tmp != 0;
                continue;
            }
            if (test ==  "DebugID") {
                ss >> debugID;
                //if (!debugID.empty())
                //    printAdditionalInfo = true;
                continue;
            }
            if (test == "printLF"){
                int tmp;
                ss >> tmp;
                printLF = tmp != 0;
                continue;
            }
            if (test == "printURF"){
                int tmp;
                ss >> tmp;
                printURF = tmp != 0;
                continue;
            }
            if (test == "printBTC"){
                int tmp;
                ss >> tmp;
                printBTC = tmp != 0;
                continue;
            }
            if (test == "printWellBTC"){
                int tmp;
                ss >> tmp;
                printWellBTC = tmp != 0;
                continue;
            }
            /*
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
            */
            if (test == "ENDofMSG") {
                out = true;
                break;
            }
            if (test.empty()) {
                outmsg += "0 ERROR: Empty message was found\n";
                out = false;
                break;
            }
            if (counter > 500) {
                outmsg += "0 ERROR: After 200 iterations I cant find ENDofMSG flag\n";
                out = false;
                break;
            }
            outmsg += "0 ERROR: UNKNOWN option [" + test + "]\n";
            out = false;
            break;
        }
        return out;
    }

}

#endif //MANTISSERVER_SCENARIO_H
